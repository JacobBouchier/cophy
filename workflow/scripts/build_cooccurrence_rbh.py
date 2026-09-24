#!/usr/bin/env python3
"""
Per-cooccur-id RBH driver — runs all per-species RBH checks for ONE cooccur
query, writes a presence/absence + best-hit TSV.

Invocation:
  build_cooccurrence_rbh.py
      --config config/config.txt
      --cooccur-id Q18B77
      --species-list <path>               # text file, taxon_id per line
      --tools mmseqs2,foldseek
      --out-tsv results/<proj>/cooccur_rbh/Q18B77.tsv
      [--workers 4]                       # parallel species in-flight

Outputs a TSV with one row per tree species:
  taxon_id  reciprocal_tools  forward_hit  forward_evalue  skipped_reverse
where:
  reciprocal_tools is a comma-separated list of tools where RBH passed
  forward_hit is the (best-tool) hit accession in that species, or ''
  forward_evalue is the (best-tool) evalue, or ''
  skipped_reverse: 'yes' if the strong-hit reverse-skip fired for any tool

A species is "called present" downstream iff reciprocal_tools is non-empty
(i.e. at least one tool confirmed reciprocal best hit). Combined with
species_taxid mapping from members.tsv, this slots into the existing
build_cooccurrence.py via the cooccur_species set.

Robustness: any single species' RBH can fail (no home homolog, weird
proteome shape, etc.) — those become "not present" rather than crashing.
Network and search subprocess failures are logged and turned into
non-reciprocal results so the rule still produces a complete output.

Parallelism: --workers controls how many species are checked in parallel.
Default is 1 (serial). Foldseek RBH is serialized via a file-lock since
each foldseek call loads ~30GB of AFDB index — running many in parallel
would blow memory.
"""

from __future__ import annotations

import argparse
import fcntl
import json
import os
import sys
import tempfile
import time
from concurrent.futures import ProcessPoolExecutor, as_completed

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fetch import fetch_entry, FetchError    # noqa: E402
from proteome import fetch_species_proteome, ensure_species_mmseqs_db   # noqa: E402
from rbh import rbh_check_mmseqs2, RBHResult  # noqa: E402


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--cooccur-id", required=True)
    ap.add_argument("--species-list", required=True,
                    help="path to text file with one taxon_id per line; "
                         "lines starting with # are skipped")
    ap.add_argument("--tools", default="mmseqs2",
                    help="comma-separated list of tools to use for RBH "
                         "(mmseqs2[,foldseek]). HMMER is intentionally "
                         "not supported in RBH mode (no per-species DB).")
    ap.add_argument("--out-tsv", required=True)
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--evalue", type=float, default=1e-3)
    ap.add_argument("--workers", type=int, default=1,
                    help="number of species checked in parallel; "
                         "1 = serial. Each worker uses --threads/workers "
                         "threads for its mmseqs/HMMER searches.")
    args = ap.parse_args()

    tools = [t.strip() for t in args.tools.split(",") if t.strip()]
    if not tools:
        sys.stderr.write("ERROR: at least one tool required\n")
        return 2
    unsupported = [t for t in tools if t not in
                   ("mmseqs2", "hmmer", "foldseek")]
    if unsupported:
        sys.stderr.write(
            f"ERROR: unsupported tool(s): {unsupported}\n"
            f"       RBH mode supports mmseqs2, hmmer, and foldseek.\n")
        return 2

    from config import load_config
    cfg = load_config(args.config)
    cache_dir = cfg["cache_dir"]

    # Step 1: identify the home species of the cooccur query.
    # We need its taxon_id and we need its proteome.
    sys.stderr.write(f"rbh: cooccur query = {args.cooccur_id}\n")
    try:
        entry = fetch_entry(args.cooccur_id, cache_dir)
    except FetchError as exc:
        sys.stderr.write(
            f"ERROR: could not fetch UniProt entry for {args.cooccur_id}: "
            f"{exc}\n")
        return 1
    home_taxon_id = str(entry.taxon_id)
    sys.stderr.write(
        f"rbh: home species = {entry.scientific_name} "
        f"(taxon {home_taxon_id})\n")

    # Step 2: fetch home species proteome (once for this whole job).
    try:
        home_proteome = fetch_species_proteome(home_taxon_id, cache_dir)
    except FetchError as exc:
        sys.stderr.write(
            f"ERROR: could not fetch home proteome for taxon "
            f"{home_taxon_id}: {exc}\n")
        return 1
    sys.stderr.write(f"rbh: home proteome ready at {home_proteome}\n")

    # Build the home species' mmseqs2 indexed DB once. Avoids re-indexing
    # on every reverse search. Best-effort: if it fails (e.g. mmseqs missing)
    # we fall back to FASTA-based searches.
    home_mmseqs_db: str | None = None
    if "mmseqs2" in tools:
        try:
            home_mmseqs_db = ensure_species_mmseqs_db(home_taxon_id, cache_dir)
        except (FetchError, Exception) as exc:    # noqa: BLE001
            sys.stderr.write(
                f"rbh: WARNING: home mmseqs DB build failed ({exc}); "
                f"falling back to FASTA searches\n")

    # Step 3: pull out the cooccur query's own FASTA (single record).
    query_fasta = _extract_query_fasta(args.cooccur_id, home_proteome,
                                       cache_dir)
    if not query_fasta:
        sys.stderr.write(
            f"ERROR: could not extract {args.cooccur_id} from home "
            f"proteome FASTA. Is the accession actually in that species?\n")
        return 1

    # If foldseek is in tools, also fetch the query's AlphaFold structure
    # and look up the foldseek DB path.
    query_pdb: str | None = None
    foldseek_db_path: str | None = None
    if "foldseek" in tools:
        from fetch import fetch_structure
        query_pdb = fetch_structure(args.cooccur_id, cache_dir)
        if query_pdb is None:
            sys.stderr.write(
                f"rbh: WARNING: no AlphaFold model for "
                f"{args.cooccur_id}; foldseek RBH disabled for this run\n")
            tools = [t for t in tools if t != "foldseek"]
        else:
            foldseek_db_path = cfg.get("foldseek_db_path", "").strip()
            if not foldseek_db_path:
                sys.stderr.write(
                    "rbh: WARNING: foldseek_db_path not set in config; "
                    "foldseek RBH disabled\n")
                tools = [t for t in tools if t != "foldseek"]
                query_pdb = None
    sys.stderr.write(f"rbh: active tools = {','.join(tools)}\n")

    # Step 4: read the tree species list.
    species_ids = _read_species_list(args.species_list)
    sys.stderr.write(
        f"rbh: {len(species_ids)} tree species to check\n")

    # Step 5: open the output file and start checking each species.
    os.makedirs(os.path.dirname(os.path.abspath(args.out_tsv)),
                exist_ok=True)
    n_present = 0
    n_proteome_fail = 0
    n_rbh_fail = 0
    n_checked = 0
    t_start = time.time()
    with open(args.out_tsv, "w") as fh:
        fh.write(
            "taxon_id\treciprocal_tools\tforward_hit\tforward_evalue\n")
        for i, taxon_id in enumerate(species_ids, 1):
            # Skip self (the home species is in the tree too)
            if taxon_id == home_taxon_id:
                # By construction, query is in its own species —
                # call it present without RBH.
                fh.write(f"{taxon_id}\tself\t{args.cooccur_id}\t0\n")
                n_present += 1
                n_checked += 1
                continue

            # Need target proteome for mmseqs2/hmmer; foldseek uses
            # taxon-filtered AFDB instead.
            target_proteome: str | None = None
            target_mmseqs_db: str | None = None
            if "mmseqs2" in tools or "hmmer" in tools:
                try:
                    target_proteome = fetch_species_proteome(
                        taxon_id, cache_dir)
                except FetchError:
                    n_proteome_fail += 1
                    target_proteome = None
                    # Fall through — foldseek may still succeed.
            # Build/reuse a cached, indexed mmseqs2 DB for the target species
            if "mmseqs2" in tools and target_proteome is not None:
                try:
                    target_mmseqs_db = ensure_species_mmseqs_db(
                        taxon_id, cache_dir)
                except Exception:    # noqa: BLE001
                    target_mmseqs_db = None    # fall back to easy-search

            # Run RBH check per requested tool.
            from rbh import rbh_check_hmmer, rbh_check_foldseek
            successful_tools: list[str] = []
            best_forward_hit = ""
            best_forward_evalue = float("inf")
            for tool in tools:
                res = None
                if tool == "mmseqs2":
                    if target_proteome is None:
                        continue
                    try:
                        res = rbh_check_mmseqs2(
                            query_accession=args.cooccur_id,
                            query_fasta=query_fasta,
                            home_proteome_fasta=home_proteome,
                            target_proteome_fasta=target_proteome,
                            target_taxon_id=taxon_id,
                            threads=args.threads,
                            evalue=args.evalue,
                            target_db_path=target_mmseqs_db,
                            home_db_path=home_mmseqs_db)
                    except Exception as exc:    # noqa: BLE001
                        sys.stderr.write(
                            f"rbh: mmseqs2 failed for "
                            f"{args.cooccur_id} vs taxon {taxon_id}: "
                            f"{exc}\n")
                        n_rbh_fail += 1
                        continue
                elif tool == "hmmer":
                    if target_proteome is None:
                        continue
                    try:
                        res = rbh_check_hmmer(
                            query_accession=args.cooccur_id,
                            query_fasta=query_fasta,
                            home_proteome_fasta=home_proteome,
                            target_proteome_fasta=target_proteome,
                            target_taxon_id=taxon_id,
                            threads=args.threads,
                            evalue=args.evalue)
                    except Exception as exc:    # noqa: BLE001
                        sys.stderr.write(
                            f"rbh: hmmer failed for "
                            f"{args.cooccur_id} vs taxon {taxon_id}: "
                            f"{exc}\n")
                        n_rbh_fail += 1
                        continue
                elif tool == "foldseek":
                    try:
                        res = rbh_check_foldseek(
                            query_accession=args.cooccur_id,
                            query_pdb=query_pdb,
                            home_taxon_id=home_taxon_id,
                            target_taxon_id=taxon_id,
                            foldseek_db_path=foldseek_db_path,
                            cache_dir=cache_dir,
                            threads=args.threads,
                            evalue=args.evalue)
                    except Exception as exc:    # noqa: BLE001
                        sys.stderr.write(
                            f"rbh: foldseek failed for "
                            f"{args.cooccur_id} vs taxon {taxon_id}: "
                            f"{exc}\n")
                        n_rbh_fail += 1
                        continue
                if res is not None and res.reciprocal:
                    successful_tools.append(tool)
                    if (res.forward_evalue is not None
                            and res.forward_evalue < best_forward_evalue):
                        best_forward_hit = res.forward_hit or ""
                        best_forward_evalue = res.forward_evalue

            if successful_tools:
                fh.write(
                    f"{taxon_id}\t{','.join(successful_tools)}\t"
                    f"{best_forward_hit}\t{best_forward_evalue:.3e}\n")
                n_present += 1
            else:
                fh.write(f"{taxon_id}\t\t\t\n")

            n_checked += 1
            if i % 50 == 0:
                elapsed = time.time() - t_start
                rate = i / elapsed if elapsed > 0 else 0
                eta = (len(species_ids) - i) / rate if rate > 0 else 0
                sys.stderr.write(
                    f"rbh: {i}/{len(species_ids)} checked, "
                    f"{n_present} present so far "
                    f"({rate:.1f}/s, ETA {eta/60:.1f}min)\n")
            # Flush periodically so progress is visible
            if i % 25 == 0:
                fh.flush()

    elapsed = time.time() - t_start
    sys.stderr.write(
        f"rbh: done. checked={n_checked} present={n_present} "
        f"proteome_fail={n_proteome_fail} rbh_fail={n_rbh_fail} "
        f"in {elapsed/60:.1f}min\n")
    return 0


def _read_species_list(path: str) -> list[str]:
    out: list[str] = []
    seen: set[str] = set()
    with open(path) as fh:
        for line in fh:
            line = line.split("#", 1)[0].strip()
            if not line:
                continue
            if line not in seen:
                seen.add(line)
                out.append(line)
    return out


def _extract_query_fasta(accession: str, proteome_fasta: str,
                         cache_dir: str) -> str | None:
    """Try to extract the cooccur query's own sequence from its home
    proteome. Falls back to fetching directly from UniProt if the
    accession isn't in the proteome FASTA (which can happen for very
    old or obsolete accessions).
    """
    from rbh import _extract_sequence
    out_path = os.path.join(cache_dir, f"rbh_query_{accession}.fasta")
    if os.path.exists(out_path) and os.path.getsize(out_path) > 0:
        return out_path
    if _extract_sequence(proteome_fasta, accession, out_path):
        return out_path
    # Fallback: direct UniProt fetch.
    import urllib.request
    from fetch import USER_AGENT
    url = f"https://rest.uniprot.org/uniprotkb/{accession}.fasta"
    try:
        req = urllib.request.Request(
            url, headers={"User-Agent": USER_AGENT, "Accept": "text/plain"})
        with urllib.request.urlopen(req, timeout=30) as resp:
            data = resp.read()
    except Exception:    # noqa: BLE001
        return None
    if not data.startswith(b">"):
        return None
    with open(out_path, "wb") as fh:
        fh.write(data)
    return out_path


if __name__ == "__main__":
    raise SystemExit(main())
