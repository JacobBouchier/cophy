#!/usr/bin/env python3
"""
NCBI Taxonomy lookups: taxon_id -> {rank_name: scientific_name}.

We use the taxon_ids that UniProt already gave us when we fetched each
protein record. Those land in the members.tsv files (column 'taxon_id').
For each unique taxon_id we hit NCBI's `efetch -db taxonomy` once,
parse the LineageEx block, and extract phylum/class/order/family.

Batching: efetch accepts comma-separated ids, so we fetch many taxon_ids
per call. Cached per-batch on disk so re-runs are free.

Pure parsers are isolated from network for unit-testability.
"""

from __future__ import annotations

import os
import sys
import xml.etree.ElementTree as ET

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fetch import _cache_path, FetchError  # noqa: E402
from fetch_neighborhood import _entrez_get, _warn  # noqa: E402


# Ranks we extract for the output CSV.
RANKS = ("domain", "phylum", "class", "order", "family")
EMPTY_LINEAGE = {r: "unclassified" for r in RANKS}


# --------------------------------------------------------------------------
# Pure parser (testable without network)
# --------------------------------------------------------------------------

def parse_taxonomy_xml(xml_text: str) -> dict[str, dict[str, str]]:
    """
    Parse an NCBI efetch taxonomy XML response.

    Returns {taxon_id: {rank: scientific_name}} for each <Taxon> in the
    TaxaSet. Only the ranks listed in RANKS are extracted; missing ranks
    fall back to 'unclassified'. Unparseable input returns {}.

    NCBI's response shape (one <Taxon> at the top level per requested id):
        <TaxaSet>
          <Taxon>
            <TaxId>2755338</TaxId>
            <ScientificName>Aceticella autotrophica</ScientificName>
            <LineageEx>
              <Taxon>
                <TaxId>1239</TaxId>
                <ScientificName>Bacillota</ScientificName>
                <Rank>phylum</Rank>
              </Taxon>
              ...
            </LineageEx>
          </Taxon>
          ...
        </TaxaSet>
    """
    out: dict[str, dict[str, str]] = {}
    if not xml_text.strip():
        return out
    try:
        root = ET.fromstring(xml_text)
    except ET.ParseError:
        return out

    # The top-level wrapper is <TaxaSet>; iterate its direct <Taxon> children.
    # IMPORTANT: also handle the case where the root IS a <Taxon> (single-id
    # response shapes sometimes lack the TaxaSet wrapper).
    if root.tag == "TaxaSet":
        top_taxa = list(root.findall("Taxon"))
    elif root.tag == "Taxon":
        top_taxa = [root]
    else:
        return out

    for taxon in top_taxa:
        tid_el = taxon.find("TaxId")
        if tid_el is None or not tid_el.text:
            continue
        tid = tid_el.text.strip()
        ranks: dict[str, str] = dict(EMPTY_LINEAGE)

        # Walk LineageEx for the requested ranks.
        lineage_ex = taxon.find("LineageEx")
        if lineage_ex is not None:
            for entry in lineage_ex.findall("Taxon"):
                rank_el = entry.find("Rank")
                name_el = entry.find("ScientificName")
                if rank_el is None or name_el is None:
                    continue
                rank = (rank_el.text or "").strip().lower()
                name = (name_el.text or "").strip()
                if rank in ranks and name:
                    ranks[rank] = name

        out[tid] = ranks
    return out

# --------------------------------------------------------------------------
# Cached fetcher
# --------------------------------------------------------------------------

_BATCH_SIZE = 200  # NCBI accepts large batches; 200 is comfortably under limit


def fetch_lineages(taxon_ids: list[str], cache_dir: str
                   ) -> dict[str, dict[str, str]]:
    """
    Map [taxon_id, ...] -> {taxon_id: {rank: name}} using NCBI efetch.

    Per-taxon results are cached individually under
    ncbi_taxonomy_<taxon_id>.json so partial work survives interruptions.
    Returns a dict with one entry per input id (using EMPTY_LINEAGE for any
    id that couldn't be resolved).
    """
    import json

    result: dict[str, dict[str, str]] = {}
    need_fetch: list[str] = []

    for tid in taxon_ids:
        tid = str(tid).strip()
        if not tid or tid == "0":
            continue
        cache = _cache_path(cache_dir, f"ncbi_taxonomy_{tid}.json")
        if os.path.exists(cache):
            try:
                with open(cache) as fh:
                    result[tid] = json.load(fh)
                continue
            except (json.JSONDecodeError, OSError):
                pass  # corrupt cache, refetch
        need_fetch.append(tid)

    if need_fetch:
        sys.stderr.write(
            f"taxonomy: fetching lineages for {len(need_fetch)} taxon_ids "
            f"(batches of {_BATCH_SIZE})\n")

    for i in range(0, len(need_fetch), _BATCH_SIZE):
        batch = need_fetch[i:i + _BATCH_SIZE]
        try:
            text = _entrez_get("efetch.fcgi", {
                "db": "taxonomy",
                "id": ",".join(batch),
                "retmode": "xml",
            })
        except FetchError as exc:
            _warn("taxonomy_fetch_failed", f"batch starting {batch[0]}: {exc}")
            # Mark all in this batch as unclassified and cache to avoid retry.
            for tid in batch:
                result[tid] = dict(EMPTY_LINEAGE)
                with open(_cache_path(cache_dir, f"ncbi_taxonomy_{tid}.json"),
                          "w") as fh:
                    json.dump(result[tid], fh)
            continue

        parsed = parse_taxonomy_xml(text)
        for tid in batch:
            ranks = parsed.get(tid, dict(EMPTY_LINEAGE))
            result[tid] = ranks
            with open(_cache_path(cache_dir, f"ncbi_taxonomy_{tid}.json"),
                      "w") as fh:
                json.dump(ranks, fh)

        if i % (5 * _BATCH_SIZE) == 0 and i > 0:
            sys.stderr.write(
                f"taxonomy: {i + len(batch)}/{len(need_fetch)} fetched\n")

    return result

def fetch_lineages_with_domain(taxon_ids: list[str], cache_dir: str
                                ) -> dict[str, dict[str, str]]:
    """
    Like fetch_lineages, but re-fetches any cached entry that lacks the
    'domain' rank (older caches predate domain support). Use this when
    you need domain info for filtering.
    """
    import json
    # First pass: which cached files are missing 'domain'?
    needs_refetch: list[str] = []
    for tid in taxon_ids:
        tid = str(tid).strip()
        if not tid or tid == "0":
            continue
        cache = _cache_path(cache_dir, f"ncbi_taxonomy_{tid}.json")
        if os.path.exists(cache):
            try:
                with open(cache) as fh:
                    data = json.load(fh)
                if "domain" not in data:
                    needs_refetch.append(tid)
            except (json.JSONDecodeError, OSError):
                needs_refetch.append(tid)
    # Delete stale caches so fetch_lineages will re-fetch them
    for tid in needs_refetch:
        cache = _cache_path(cache_dir, f"ncbi_taxonomy_{tid}.json")
        try:
            os.unlink(cache)
        except OSError:
            pass
    if needs_refetch:
        sys.stderr.write(
            f"taxonomy: re-fetching {len(needs_refetch)} cached entries "
            f"to add 'domain' field\n")
    return fetch_lineages(taxon_ids, cache_dir)
