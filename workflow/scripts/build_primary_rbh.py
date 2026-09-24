#!/usr/bin/env python3
"""
Per-primary RBH validation — runs RBH for each species already in a
primary's members.tsv, writes an RBH-annotated members file.

The existing search → consensus → dedup pipeline picks "the species
that made the cutoff." This script then runs RBH against each of those
species' proteomes to validate they truly contain a reciprocal-best
ortholog of the primary query — not just a high-scoring paralog.

Invocation:
  build_primary_rbh.py
      --config config/config.txt
      --primary-id Q18B79
      --members results/<proj>/groups/Q18B79.members.tsv
      --tools mmseqs2
      --out-tsv results/<proj>/groups/Q18B79.members_rbh.tsv

Output schema (extends members.tsv with two columns):
  species  accession  taxon_id  reviewed  tools  score  evalue
  rbh_reciprocal_tools  rbh_forward_evalue

rbh_reciprocal_tools is a comma-separated list of tools where RBH
passed for that species (e.g. 'mmseqs2,foldseek'). Empty string means
RBH failed for all tested tools.

Downstream rules will filter or annotate based on this column.
"""

from __future__ import annotations

import argparse
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fetch import fetch_entry, FetchError    # noqa: E402
from proteome import fetch_species_proteome, ensure_species_mmseqs_db   # noqa: E402
from rbh import rbh_check_mmseqs2             # noqa: E402


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--primary-id", required=True,
                    help="UniProt accession of the primary query")
    ap.add_argument("--members", required=True,
                    help="path to <primary>.members.tsv from build_group")
    ap.add_argument("--out-tsv", required=True)
    ap.add_argument("--tools", default="mmseqs2",
                    help="comma-separated list of RBH tools (mmseqs2[,foldseek])")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--evalue", type=float, default=1e-3)
    args = ap.parse_args()

    tools = [t.strip() for t in args.tools.split(",") if t.strip()]
    if not tools:
        sys.stderr.write("ERROR: at least one tool required\n")
        return 2
    unsupported = [t for t in tools if t not in
                   ("mmseqs2", "hmmer", "foldseek")]
    if unsupported:
        sys.stderr.write(
            f"ERROR: unsupported tool(s) {unsupported}; RBH supports "
            f"mmseqs2, hmmer, and foldseek\n")
        return 2

    from config import load_config
    cfg = load_config(args.config)
    cache_dir = cfg["cache_dir"]

    sys.stderr.write(f"primary_rbh: query = {args.primary_id}\n")

    # Step 1: fetch home species + proteome (once for this whole job)
    try:
        entry = fetch_entry(args.primary_id, cache_dir)
    except FetchError as exc:
        sys.stderr.write(
            f"ERROR: could not fetch UniProt entry for {args.primary_id}: "
            f"{exc}\n")
        return 1
    home_taxon = str(entry.taxon_id)
    sys.stderr.write(
        f"primary_rbh: home species = {entry.scientific_name} "
        f"(taxon {home_taxon})\n")
    try:
        home_proteome = fetch_species_proteome(home_taxon, cache_dir)
    except FetchError as exc:
        sys.stderr.write(
            f"ERROR: could not fetch home proteome: {exc}\n")
        return 1

    # Build home species' indexed mmseqs DB once (cached, indexed mmseqs DB reused across species).
    home_mmseqs_db: str | None = None
    if "mmseqs2" in tools:
        try:
            home_mmseqs_db = ensure_species_mmseqs_db(home_taxon, cache_dir)
        except Exception as exc:    # noqa: BLE001
            sys.stderr.write(
                f"primary_rbh: WARNING: home mmseqs DB build failed "
                f"({exc}); falling back to FASTA searches\n")

    # Step 2: extract the primary's own FASTA from the home proteome
    from rbh import _extract_sequence
    query_fasta = os.path.join(cache_dir, f"rbh_primary_{args.primary_id}.fasta")
    if not (os.path.exists(query_fasta) and os.path.getsize(query_fasta) > 0):
        if not _extract_sequence(home_proteome, args.primary_id, query_fasta):
            # Fall back to direct UniProt fetch
            import urllib.request
            from fetch import USER_AGENT
            url = f"https://rest.uniprot.org/uniprotkb/{args.primary_id}.fasta"
            try:
                req = urllib.request.Request(
                    url, headers={"User-Agent": USER_AGENT,
                                  "Accept": "text/plain"})
                with urllib.request.urlopen(req, timeout=30) as resp:
                    data = resp.read()
                if data.startswith(b">"):
                    with open(query_fasta, "wb") as fh:
                        fh.write(data)
            except Exception:    # noqa: BLE001
                sys.stderr.write(
                    f"ERROR: could not extract or fetch FASTA for "
                    f"{args.primary_id}\n")
                return 1

    # If foldseek is in tools, also fetch the query's AlphaFold structure.
    # If unavailable, foldseek RBH will be skipped for this primary.
    query_pdb: str | None = None
    foldseek_db_path: str | None = None
    if "foldseek" in tools:
        from fetch import fetch_structure
        query_pdb = fetch_structure(args.primary_id, cache_dir)
        if query_pdb is None:
            sys.stderr.write(
                f"primary_rbh: WARNING: no AlphaFold model for "
                f"{args.primary_id}; foldseek RBH disabled for this run\n")
            tools = [t for t in tools if t != "foldseek"]
        else:
            foldseek_db_path = cfg.get("foldseek_db_path", "").strip()
            if not foldseek_db_path:
                sys.stderr.write(
                    "primary_rbh: WARNING: foldseek_db_path not set in "
                    "config; foldseek RBH disabled\n")
                tools = [t for t in tools if t != "foldseek"]
                query_pdb = None

    # Step 3: read members.tsv, get list of species to RBH-check
    rows: list[dict[str, str]] = []
    with open(args.members) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < len(header):
                continue
            rows.append(dict(zip(header, f)))

    sys.stderr.write(
        f"primary_rbh: {len(rows)} species in members.tsv\n")
    sys.stderr.write(f"primary_rbh: active tools = {','.join(tools)}\n")

    # Step 4: open output, RBH-check each species, emit annotated row
    os.makedirs(os.path.dirname(os.path.abspath(args.out_tsv)),
                exist_ok=True)
    n_pass = 0
    n_proteome_fail = 0
    n_rbh_fail = 0
    t_start = time.time()
    from rbh import rbh_check_hmmer, rbh_check_foldseek
    with open(args.out_tsv, "w") as out:
        out_header = header + ["rbh_reciprocal_tools", "rbh_forward_evalue"]
        out.write("\t".join(out_header) + "\n")
        for i, row in enumerate(rows, 1):
            taxon = str(row.get("taxon_id", "")).strip()
            sp = row.get("species", "")
            # Species == home -> always pass (self)
            if taxon == home_taxon:
                row["rbh_reciprocal_tools"] = "self"
                row["rbh_forward_evalue"] = "0"
                n_pass += 1
                out.write("\t".join(row.get(c, "") for c in out_header) + "\n")
                continue
            if not taxon or taxon == "0":
                row["rbh_reciprocal_tools"] = ""
                row["rbh_forward_evalue"] = ""
                out.write("\t".join(row.get(c, "") for c in out_header) + "\n")
                continue
            # We need the target species' proteome FASTA for mmseqs2/hmmer RBH.
            # Foldseek doesn't need it (uses taxon-filtered AFDB instead).
            target_proteome: str | None = None
            target_mmseqs_db: str | None = None
            if "mmseqs2" in tools or "hmmer" in tools:
                try:
                    target_proteome = fetch_species_proteome(taxon, cache_dir)
                except FetchError:
                    n_proteome_fail += 1
                    # If we can't get the proteome, mmseqs2 and hmmer fail
                    # for this species — but foldseek can still try if active.
                    target_proteome = None
            if "mmseqs2" in tools and target_proteome is not None:
                try:
                    target_mmseqs_db = ensure_species_mmseqs_db(
                        taxon, cache_dir)
                except Exception:    # noqa: BLE001
                    target_mmseqs_db = None

            successful_tools: list[str] = []
            best_evalue = float("inf")
            for tool in tools:
                res = None
                if tool == "mmseqs2":
                    if target_proteome is None:
                        continue
                    try:
                        res = rbh_check_mmseqs2(
                            query_accession=args.primary_id,
                            query_fasta=query_fasta,
                            home_proteome_fasta=home_proteome,
                            target_proteome_fasta=target_proteome,
                            target_taxon_id=taxon,
                            threads=args.threads,
                            evalue=args.evalue,
                            target_db_path=target_mmseqs_db,
                            home_db_path=home_mmseqs_db)
                    except Exception as exc:    # noqa: BLE001
                        sys.stderr.write(
                            f"primary_rbh: mmseqs2 failed for "
                            f"{args.primary_id} vs {sp} (taxon {taxon}): "
                            f"{exc}\n")
                        n_rbh_fail += 1
                        continue
                elif tool == "hmmer":
                    if target_proteome is None:
                        continue
                    try:
                        res = rbh_check_hmmer(
                            query_accession=args.primary_id,
                            query_fasta=query_fasta,
                            home_proteome_fasta=home_proteome,
                            target_proteome_fasta=target_proteome,
                            target_taxon_id=taxon,
                            threads=args.threads,
                            evalue=args.evalue)
                    except Exception as exc:    # noqa: BLE001
                        sys.stderr.write(
                            f"primary_rbh: hmmer failed for "
                            f"{args.primary_id} vs {sp} (taxon {taxon}): "
                            f"{exc}\n")
                        n_rbh_fail += 1
                        continue
                elif tool == "foldseek":
                    try:
                        res = rbh_check_foldseek(
                            query_accession=args.primary_id,
                            query_pdb=query_pdb,
                            home_taxon_id=home_taxon,
                            target_taxon_id=taxon,
                            foldseek_db_path=foldseek_db_path,
                            cache_dir=cache_dir,
                            threads=args.threads,
                            evalue=args.evalue)
                    except Exception as exc:    # noqa: BLE001
                        sys.stderr.write(
                            f"primary_rbh: foldseek failed for "
                            f"{args.primary_id} vs {sp} (taxon {taxon}): "
                            f"{exc}\n")
                        n_rbh_fail += 1
                        continue
                if res is not None and res.reciprocal:
                    successful_tools.append(tool)
                    if (res.forward_evalue is not None
                            and res.forward_evalue < best_evalue):
                        best_evalue = res.forward_evalue

            # Track proteome failures separately - if a species has no proteome
            # available but foldseek (which doesn't need it) succeeded, that's
            # not really a proteome_fail in terms of outcome.
            if target_proteome is None and not successful_tools:
                # Already counted above; just record the row as fail.
                row["rbh_reciprocal_tools"] = ""
                row["rbh_forward_evalue"] = ""
                out.write("\t".join(row.get(c, "") for c in out_header) + "\n")
                continue

            if successful_tools:
                row["rbh_reciprocal_tools"] = ",".join(successful_tools)
                row["rbh_forward_evalue"] = (
                    f"{best_evalue:.3e}" if best_evalue < float("inf") else "")
                n_pass += 1
            else:
                row["rbh_reciprocal_tools"] = ""
                row["rbh_forward_evalue"] = ""

            out.write("\t".join(row.get(c, "") for c in out_header) + "\n")

            if i % 50 == 0:
                elapsed = time.time() - t_start
                rate = i / elapsed if elapsed > 0 else 0
                eta = (len(rows) - i) / rate if rate > 0 else 0
                sys.stderr.write(
                    f"primary_rbh: {i}/{len(rows)} checked, "
                    f"{n_pass} passed ({rate:.1f}/s, ETA {eta/60:.1f}min)\n")
                out.flush()

    elapsed = time.time() - t_start
    sys.stderr.write(
        f"primary_rbh: done. checked={len(rows)} passed={n_pass} "
        f"proteome_fail={n_proteome_fail} rbh_fail={n_rbh_fail} "
        f"in {elapsed/60:.1f}min\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
