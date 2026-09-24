#!/usr/bin/env python3
"""
Fetch per-species proteomes from UniProt for RBH-based cooccurrence detection.

For RBH we need to search each cooccur query against each tree species'
proteins. To do this:

  * `fetch_species_proteome(taxon_id, cache_dir)` returns a FASTA file path
    containing all UniProt proteins (Swiss-Prot + TrEMBL) for that taxon.
  * `fetch_species_structures(taxon_id, cache_dir)` returns a directory of
    AlphaFold PDB files for that taxon (used for foldseek searches).
    Note: not all proteins have AlphaFold models.

Caching: files are stored as `species_<taxon_id>.fasta` and
`species_<taxon_id>_structures/` to be reused across runs and across
cooccur queries.

UniProt's REST API is rate-limited; we use modest delays between fetches.
For very large proteomes (thousands of proteins) the API paginates - we
use UniProt's stream endpoint that returns all records in one call.
"""

from __future__ import annotations

import os
import sys
import time
import urllib.parse
import urllib.request

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fetch import (    # noqa: E402
    FetchError, _cache_path, _http_get, USER_AGENT,
)


UNIPROT_REST = "https://rest.uniprot.org"
AFDB_BASE = "https://alphafold.ebi.ac.uk/files"


def fetch_species_proteome(taxon_id: str, cache_dir: str,
                           include_unreviewed: bool = True) -> str:
    """
    Fetch all UniProt protein sequences for a given NCBI taxon ID.

    Returns the path to a FASTA file `species_<taxon_id>.fasta` with one
    record per protein, header format >ACCESSION SCIENTIFIC_NAME.

    On cache miss, queries UniProt's stream API which returns all records
    matching the query in a single response. For most bacterial species,
    proteomes are 1-10k proteins (few MB).

    `include_unreviewed=True` includes TrEMBL entries (the realistic case
    for non-model bacteria); set False to restrict to SwissProt only.
    """
    taxon_id = str(taxon_id).strip()
    if not taxon_id or taxon_id == "0":
        raise FetchError(f"Invalid taxon_id: {taxon_id!r}")

    cache = _cache_path(cache_dir, f"species_{taxon_id}.fasta")
    if os.path.exists(cache) and os.path.getsize(cache) > 0:
        return cache

    # Build the UniProt query. Filter on taxonomy_id (organism-only,
    # exclude lineage children — we want this exact species, not all of
    # its sub-strains, to avoid huge proteome bloat).
    if include_unreviewed:
        query = f"taxonomy_id:{taxon_id}"
    else:
        query = f"taxonomy_id:{taxon_id} AND reviewed:true"

    # /uniprotkb/stream returns all records in one shot (no pagination).
    params = {
        "query": query,
        "format": "fasta",
        "compressed": "false",
    }
    url = f"{UNIPROT_REST}/uniprotkb/stream?" + urllib.parse.urlencode(params)

    sys.stderr.write(f"proteome: fetching for taxon {taxon_id}...\n")
    # Retry on transient HTTP failures: UniProt's stream endpoint
    # occasionally drops chunked-encoding connections partway through
    # for large proteomes. Exponential backoff between attempts.
    import http.client
    last_exc: Exception | None = None
    data: bytes = b""
    for attempt in range(4):
        try:
            data = _http_get(url, accept="text/plain")
            break
        except (FetchError, http.client.IncompleteRead,
                http.client.HTTPException,
                ConnectionError, TimeoutError, OSError) as exc:
            last_exc = exc
            wait = 2 ** attempt   # 1, 2, 4, 8 seconds
            sys.stderr.write(
                f"proteome: fetch attempt {attempt+1} for taxon "
                f"{taxon_id} failed ({type(exc).__name__}: {exc}); "
                f"retrying in {wait}s\n")
            time.sleep(wait)
    else:
        # All retries exhausted.
        raise FetchError(
            f"proteome fetch failed for taxon {taxon_id} after retries: "
            f"{last_exc}")

    with open(cache, "wb") as fh:
        fh.write(data)

    # Sanity check: should look like FASTA, not an error page.
    sample = data[:200].decode("utf-8", errors="replace")
    if not sample.startswith(">"):
        # Likely got an error page. Don't keep a bad cache.
        os.remove(cache)
        raise FetchError(
            f"proteome for taxon {taxon_id} returned non-FASTA: "
            f"{sample[:100]!r}")

    n_records = data.count(b"\n>") + (1 if data.startswith(b">") else 0)
    sys.stderr.write(
        f"proteome: cached taxon {taxon_id} -> {cache} ({n_records} records)\n")
    return cache


def ensure_species_mmseqs_db(taxon_id: str, cache_dir: str) -> str:
    """Build (and cache) an mmseqs2 indexed DB for one species' proteome.

    Returns the path to the indexed DB. The DB is created from the cached
    proteome FASTA via `mmseqs createdb`. Subsequent calls return the
    cached DB path without rebuilding.

    Speeds up RBH significantly: `mmseqs easy-search` re-creates the DB
    on every call (~1-2s overhead per species). Using a pre-built DB
    with `mmseqs search` skips that.
    """
    import subprocess
    taxon_id = str(taxon_id).strip()
    db_dir = _cache_path(cache_dir, f"species_{taxon_id}_mmseqs")
    db_path = os.path.join(db_dir, "db")
    # The marker file mmseqs2 always creates on success
    marker = db_path + ".dbtype"
    if os.path.exists(marker):
        return db_path

    fasta_path = fetch_species_proteome(taxon_id, cache_dir)
    os.makedirs(db_dir, exist_ok=True)
    sys.stderr.write(f"mmseqs_db: building for taxon {taxon_id}...\n")
    cmd = ["mmseqs", "createdb", fasta_path, db_path, "--shuffle", "0"]
    try:
        subprocess.run(cmd, check=True, capture_output=True)
    except subprocess.CalledProcessError as exc:
        raise FetchError(
            f"mmseqs createdb failed for taxon {taxon_id}: "
            f"{exc.stderr.decode('utf-8', errors='replace')[:300]}")
    return db_path


def fetch_species_structures(taxon_id: str, cache_dir: str,
                             accessions: list[str] | None = None,
                             ) -> str:
    """
    Fetch AlphaFold PDB structures for proteins in a species.

    Returns the path to a directory `species_<taxon_id>_structures/`
    containing one .pdb file per protein that has an AlphaFold model.
    NOT every protein will have a model - the directory will be a subset
    of the species' full proteome.

    If `accessions` is None, this function reads the species FASTA and
    fetches structures for every accession in it. If `accessions` is
    given, only those structures are fetched.

    Uses the AlphaFold API to look up the canonical PDB URL per accession
    (same approach as fetch_structure() in fetch.py).
    """
    import json
    taxon_id = str(taxon_id).strip()
    struct_dir = _cache_path(cache_dir, f"species_{taxon_id}_structures")
    os.makedirs(struct_dir, exist_ok=True)

    if accessions is None:
        fasta_path = fetch_species_proteome(taxon_id, cache_dir)
        accessions = _parse_fasta_accessions(fasta_path)

    n_fetched = 0
    n_skip = 0
    n_no_model = 0
    for acc in accessions:
        pdb_path = os.path.join(struct_dir, f"AF-{acc}.pdb")
        if os.path.exists(pdb_path):
            n_skip += 1
            continue
        api_url = f"https://alphafold.ebi.ac.uk/api/prediction/{acc}"
        try:
            req = urllib.request.Request(
                api_url, headers={
                    "Accept": "application/json",
                    "User-Agent": USER_AGENT})
            with urllib.request.urlopen(req, timeout=30) as resp:
                payload = json.loads(resp.read())
        except Exception:    # noqa: BLE001
            n_no_model += 1
            continue
        if not payload or not isinstance(payload, list):
            n_no_model += 1
            continue
        pdb_url = payload[0].get("pdbUrl")
        if not pdb_url:
            n_no_model += 1
            continue
        try:
            data = _http_get(pdb_url, accept="text/plain")
        except FetchError:
            n_no_model += 1
            continue
        with open(pdb_path, "wb") as fh:
            fh.write(data)
        n_fetched += 1
        time.sleep(0.05)  # be polite to AFDB

    sys.stderr.write(
        f"structures: taxon {taxon_id}: "
        f"fetched={n_fetched} cached={n_skip} no_model={n_no_model}\n")
    return struct_dir


def _parse_fasta_accessions(fasta_path: str) -> list[str]:
    """Read all UniProt accessions from FASTA headers like `>sp|ACC|name`
    or `>tr|ACC|name`. Returns deduped list, preserving first-seen order.
    """
    seen: list[str] = []
    seen_set: set[str] = set()
    with open(fasta_path) as fh:
        for line in fh:
            if not line.startswith(">"):
                continue
            # UniProt FASTA headers: >sp|ACC|name OS=... OX=... etc.
            # Or >tr|ACC|name.
            f = line[1:].split("|")
            if len(f) < 2:
                continue
            acc = f[1].strip()
            if acc and acc not in seen_set:
                seen.append(acc)
                seen_set.add(acc)
    return seen
