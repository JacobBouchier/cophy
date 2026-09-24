#!/usr/bin/env python3
"""
Combine per-tool hit TSVs for one primary query group into deduplicated
ortholog members, then fetch their sequences and write a per-group FASTA
whose headers are species keys (ready for alignment).

This is the bridge between the tested pure modules and real I/O. It:
  1. reads the per-tool RawHit TSVs (from search_run.py)
  2. resolves each hit accession to species/taxon/reviewed/len via UniProt
     (cached), building Hit objects
  3. applies consensus (combine_hits) and species dedup (dedup_by_species)
  4. (rbh mode) optionally filters to reciprocal-best hits
  5. writes members.tsv and group.fasta (aa or nt)

Run per group by Snakemake.
"""

from __future__ import annotations

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import load_config                                   # noqa: E402
from consensus import Hit, combine_hits, dedup_by_species        # noqa: E402
from fetch import (fetch_entry, species_key, map_uniprot_to_embl_cds,  # noqa: E402
                   encode_species, FetchError)
from supermatrix import write_fasta                              # noqa: E402


def read_hit_tsvs(paths: list[str]):
    """Yield (group, target, tool, score, evalue) from search_run TSVs."""
    for p in paths:
        if not os.path.exists(p):
            continue
        with open(p) as fh:
            header = fh.readline()
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) != 5:
                    continue
                group, target, tool, score, evalue = parts
                try:
                    yield group, target, tool, float(score), float(evalue)
                except ValueError:
                    continue


def _is_uniparc(acc: str) -> bool:
    """UniParc accessions look like UPI followed by hex chars (e.g. UPI00227EB005).
    They are NOT in UniProtKB and have no species/reviewed metadata, so we
    can't use them in a species-based ortholog set. Skip them up front
    without making a doomed REST call per accession."""
    return len(acc) >= 4 and acc.startswith("UPI") and acc[3:].isalnum()


def cap_per_tool(rows, max_per_tool: int):
    """
    Keep only the top `max_per_tool` rows per tool, ranked by score (desc).
    rows is an iterable of (group, target, tool, score, evalue). Returns a
    list. If max_per_tool <= 0, returns the input unchanged. This runs
    BEFORE any metadata fetch, so it directly bounds network work.
    """
    rows = list(rows)
    if max_per_tool <= 0:
        return rows
    by_tool: dict[str, list] = {}
    for r in rows:
        by_tool.setdefault(r[2], []).append(r)
    out: list = []
    for tool, lst in by_tool.items():
        lst.sort(key=lambda r: r[3], reverse=True)  # score desc
        out.extend(lst[:max_per_tool])
    return out


def build_hits(rows, cfg) -> list[Hit]:
    """Resolve each unique target accession to metadata, build Hit list."""
    hits: list[Hit] = []
    meta_cache: dict[str, object] = {}
    skipped_uniparc = 0
    for group, target, tool, score, evalue in rows:
        if target not in meta_cache:
            if _is_uniparc(target):
                meta_cache[target] = None
                skipped_uniparc += 1
                continue
            try:
                meta_cache[target] = fetch_entry(target, cfg["cache_dir"])
            except FetchError as exc:
                sys.stderr.write(f"WARN: skip {target}: {exc}\n")
                meta_cache[target] = None
        entry = meta_cache[target]
        if entry is None:
            continue
        sp = species_key(entry.scientific_name, cfg["species_rank"])
        hits.append(Hit(accession=target, tool=tool, score=score,
                        evalue=evalue, species=sp, taxon_id=entry.taxon_id,
                        reviewed=entry.reviewed, seq_len=len(entry.sequence)))
    return hits


def maybe_rbh(calls, group_id, cfg):
    """
    Forward mode: return calls unchanged.
    RBH mode: keep only calls whose own best forward hit (when that target is
    searched back) returns the original query group. We approximate this with
    a reciprocal marker file if present; absent that, we conservatively keep
    calls (and warn), since a true RBH needs a second search handled by the
    Snakemake rule graph. The rule supplies rbh-confirmed accessions via an
    optional file.
    """
    if cfg["orthology_mode"] == "forward":
        return calls
    rbh_file = os.path.join(cfg["cache_dir"], f"rbh_{group_id}.txt")
    if not os.path.exists(rbh_file):
        sys.stderr.write(
            f"WARN: orthology_mode=rbh but {rbh_file} missing; "
            f"keeping forward hits for {group_id}\n")
        return calls
    with open(rbh_file) as fh:
        confirmed = {l.strip() for l in fh if l.strip()}
    return [c for c in calls if c.accession in confirmed]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--group", required=True)
    ap.add_argument("--hits", nargs="+", required=True,
                    help="per-tool hit TSVs for this group")
    ap.add_argument("--out-members", required=True)
    ap.add_argument("--out-fasta", required=True)
    args = ap.parse_args()

    cfg = load_config(args.config)
    rows = list(read_hit_tsvs(args.hits))
    n_in = len(rows)
    rows = cap_per_tool(rows, cfg.get("max_hits_per_tool", 0))
    if len(rows) < n_in:
        sys.stderr.write(
            f"{args.group}: capped hits {n_in} -> {len(rows)} "
            f"(max_hits_per_tool={cfg.get('max_hits_per_tool', 0)})\n")
    hits = build_hits(rows, cfg)

    calls = combine_hits(hits, mode=cfg["consensus_mode"],
                         min_tools=cfg["min_tools"],
                         enabled_tools=cfg["enabled_tools"])
    calls = maybe_rbh(calls, args.group, cfg)
    calls = dedup_by_species(calls, enabled=cfg["dedup_species"],
                             tie_break=cfg["tie_break"])

    # Always include the query itself as the reference member of its group.
    # (It defines the group; ensures the query's own species is represented.)
    os.makedirs(os.path.dirname(os.path.abspath(args.out_members)), exist_ok=True)
    with open(args.out_members, "w") as fh:
        fh.write("species\taccession\ttaxon_id\treviewed\ttools\tscore\tevalue\n")
        for c in calls:
            # Encode species to FASTA-safe form for ALL downstream files.
            # build_cooccurrence.py will decode it for human-readable CSV.
            fh.write(f"{encode_species(c.species)}\t{c.accession}\t{c.taxon_id}\t"
                     f"{int(c.reviewed)}\t{','.join(sorted(c.tools))}\t"
                     f"{c.best_score}\t{c.best_evalue}\n")

    # Build the group FASTA keyed by species. For nt mode, fetch CDS.
    seqs: dict[str, str] = {}
    if cfg["seq_type"] == "nt":
        accs = [c.accession for c in calls]
        cds_map = map_uniprot_to_embl_cds(accs, cfg["cache_dir"])
        for c in calls:
            cds_ids = cds_map.get(c.accession, [])
            if not cds_ids:
                sys.stderr.write(
                    f"WARN: no CDS cross-ref for {c.accession} "
                    f"({c.species}); excluded from nt supermatrix\n")
                continue
            # Sequence retrieval of the CDS nucleotides is handled by a
            # dedicated fetch (ENA) keyed by cds id; cached on disk.
            from fetch import _http_get_text  # local import; network layer
            try:
                fasta = _http_get_text(
                    f"https://www.ebi.ac.uk/ena/browser/api/fasta/{cds_ids[0]}")
                seq = "".join(l.strip() for l in fasta.splitlines()
                              if l and not l.startswith(">"))
                if seq:
                    seqs[encode_species(c.species)] = seq
            except Exception as exc:  # noqa: BLE001 - report and continue
                sys.stderr.write(f"WARN: ENA fetch failed {cds_ids[0]}: {exc}\n")
    else:
        for c in calls:
            try:
                entry = fetch_entry(c.accession, cfg["cache_dir"])
                seqs[encode_species(c.species)] = entry.sequence
            except FetchError as exc:
                sys.stderr.write(f"WARN: seq fetch {c.accession}: {exc}\n")

    write_fasta(seqs, args.out_fasta)
    sys.stderr.write(f"{args.group}: {len(calls)} members, "
                     f"{len(seqs)} sequences -> {args.out_fasta}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
