#!/usr/bin/env python3
"""
Apply RBH outcomes across all primaries to filter or annotate each
primary's members file for downstream alignment/tree-building.

Algorithm (a species passes if RBH validates it under ANY primary):
  1. Read every primary's members_rbh.tsv (the per-primary RBH output).
  2. Build a set of (species, RBH-passed-anywhere?) by union across primaries:
       a species "passes" iff RBH succeeded for it in at least ONE primary.
  3. For each primary's members file, emit:
       * `filter` mode -> only rows whose species is in the passed set.
       * `annotate` mode -> all rows, with an `rbh_passed_anywhere` column.

We need step 2 because tree leaves are unioned across primaries: a species
in Q18B79.members but not Q18B78.members still appears in the tree iff
its Q18B79 RBH passed. The "ANY primary" rule says: include the species
if RBH validates ANY of its primary hits, even if a *different* primary's
RBH failed.

Invocation:
  build_primary_rbh_filter.py
      --primary-rbh A.members_rbh.tsv B.members_rbh.tsv
      --members-in A.members_rbh.tsv
      --mode filter           # or 'annotate'
      --out A.members_filtered.tsv

Each call processes ONE primary (--members-in), but the cross-primary
union is computed from ALL primary RBH files (--primary-rbh).
"""

from __future__ import annotations

import argparse
import os
import sys


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--primary-rbh", nargs="+", required=True,
                    help="all primaries' members_rbh.tsv files for the "
                         "cross-primary union calculation")
    ap.add_argument("--members-in", required=True,
                    help="the members_rbh.tsv to process (one of the files "
                         "in --primary-rbh)")
    ap.add_argument("--mode", choices=["filter", "annotate"],
                    required=True)
    ap.add_argument("--out", required=True,
                    help="output members TSV (annotated or filtered)")
    ap.add_argument("--fasta-in", default="",
                    help="optional original FASTA from build_group; when "
                         "given and mode is 'filter', also writes "
                         "--fasta-out containing only passing species")
    ap.add_argument("--fasta-out", default="",
                    help="output FASTA path (required iff --fasta-in given "
                         "and mode is 'filter')")
    ap.add_argument("--kingdom-filter", default="",
                    help="comma-separated list of NCBI domains to keep "
                         "(e.g., 'Bacteria' or 'Bacteria,Archaea'). "
                         "Empty = no taxonomy filter.")
    ap.add_argument("--cache-dir", default="resources/cache",
                    help="taxonomy cache directory (used only when "
                         "--kingdom-filter is set)")
    args = ap.parse_args()

    if args.mode == "filter" and args.fasta_in and not args.fasta_out:
        sys.stderr.write(
            "ERROR: --fasta-out required when --fasta-in given in filter mode\n")
        return 2

    # Step 1: union of species that pass RBH in ANY primary, OR were found
    # via foldseek-only in the original search (foldseek RBH is too memory-
    # expensive to do per-species, but a foldseek hit at the original search
    # stage is independently meaningful — it represents structural homology
    # that mmseqs2/HMMER can't detect for distantly-related orthologs).
    passed_anywhere: set[str] = set()
    foldseek_only_anywhere: set[str] = set()
    for path in args.primary_rbh:
        with open(path) as fh:
            header = fh.readline().rstrip("\n").split("\t")
            try:
                sp_col = header.index("species")
                rbh_col = header.index("rbh_reciprocal_tools")
                tools_col = header.index("tools")
            except ValueError:
                sys.stderr.write(
                    f"WARN: {path}: missing species/rbh/tools column\n")
                continue
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) <= max(sp_col, rbh_col, tools_col):
                    continue
                sp = f[sp_col]
                rbh_pass = bool(f[rbh_col].strip())
                tools_used = set(
                    t.strip() for t in f[tools_col].split(",")
                    if t.strip())
                if rbh_pass:
                    passed_anywhere.add(sp)
                elif tools_used == {"foldseek"}:
                    # Only foldseek found this species (no mmseqs2/HMMER).
                    # Trust the original foldseek search; RBH would
                    # by-definition fail for distant structural orthologs.
                    foldseek_only_anywhere.add(sp)
                    passed_anywhere.add(sp)

    sys.stderr.write(
        f"primary_rbh_filter: {len(passed_anywhere)} species pass "
        f"({len(passed_anywhere) - len(foldseek_only_anywhere)} via RBH, "
        f"{len(foldseek_only_anywhere)} via foldseek-only original search)\n")
    # Optional kingdom filter: look up each species' NCBI domain and drop
    # those outside the allowed list (e.g. keep only Bacteria).
    allowed_domains = set()
    if args.kingdom_filter.strip():
        allowed_domains = set(
            d.strip() for d in args.kingdom_filter.split(",") if d.strip())
        sys.stderr.write(
            f"primary_rbh_filter: applying kingdom filter -> "
            f"keep only {allowed_domains}\n")

        sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
        from taxonomy import fetch_lineages_with_domain
        cache_dir = args.cache_dir

        # Collect taxon_ids for the species we're about to filter.
        sp_to_taxon: dict[str, str] = {}
        for path in args.primary_rbh:
            with open(path) as fh:
                header = fh.readline().rstrip("\n").split("\t")
                sp_col = header.index("species")
                tid_col = header.index("taxon_id")
                for line in fh:
                    f = line.rstrip("\n").split("\t")
                    if len(f) <= max(sp_col, tid_col):
                        continue
                    sp = f[sp_col]
                    if sp in passed_anywhere and sp not in sp_to_taxon:
                        sp_to_taxon[sp] = f[tid_col].strip()

        unique_taxa = sorted(set(t for t in sp_to_taxon.values() if t))
        sys.stderr.write(
            f"primary_rbh_filter: looking up taxonomy for "
            f"{len(unique_taxa)} unique taxa\n")
        lineages = fetch_lineages_with_domain(unique_taxa, cache_dir)

        # Filter passed_anywhere by domain
        kept: set[str] = set()
        dropped_by_domain: dict[str, int] = {}
        for sp in passed_anywhere:
            tid = sp_to_taxon.get(sp, "")
            domain = lineages.get(tid, {}).get("domain", "unclassified")
            if domain in allowed_domains:
                kept.add(sp)
            else:
                dropped_by_domain[domain] = dropped_by_domain.get(domain, 0) + 1
        sys.stderr.write(
            f"primary_rbh_filter: kingdom filter kept "
            f"{len(kept)}/{len(passed_anywhere)} species\n")
        for dom, n in sorted(dropped_by_domain.items(),
                             key=lambda x: -x[1]):
            sys.stderr.write(f"    dropped {n} from domain '{dom}'\n")
        passed_anywhere = kept

    # Step 2: process the input file according to mode
    os.makedirs(os.path.dirname(os.path.abspath(args.out)),
                exist_ok=True)
    n_in = 0
    n_out = 0
    with open(args.members_in) as inp, open(args.out, "w") as out:
        header_line = inp.readline()
        header = header_line.rstrip("\n").split("\t")
        if args.mode == "annotate":
            out.write(header_line.rstrip("\n") +
                      "\trbh_passed_anywhere\n")
        else:
            out.write(header_line)
        sp_col = header.index("species")
        for line in inp:
            n_in += 1
            f = line.rstrip("\n").split("\t")
            if len(f) <= sp_col:
                continue
            sp = f[sp_col]
            passed = sp in passed_anywhere
            if args.mode == "filter":
                if passed:
                    out.write(line)
                    n_out += 1
            else:
                out.write(line.rstrip("\n") +
                          f"\t{'yes' if passed else 'no'}\n")
                n_out += 1

    sys.stderr.write(
        f"primary_rbh_filter: {args.mode}: {n_in} in, {n_out} out\n")

    # If a FASTA was provided in filter mode, also produce a filtered FASTA
    # keyed by species name (which is what the alignment step consumes).
    if args.mode == "filter" and args.fasta_in:
        os.makedirs(os.path.dirname(os.path.abspath(args.fasta_out)),
                    exist_ok=True)
        n_records_in = 0
        n_records_out = 0
        keep = False
        with open(args.fasta_in) as inp, open(args.fasta_out, "w") as out:
            for line in inp:
                if line.startswith(">"):
                    n_records_in += 1
                    # FASTA headers in our build_group output are
                    # ">SPECIES_NAME" (encoded with underscores). Match
                    # against the passed set, which also uses underscores.
                    sp = line[1:].split()[0].strip()
                    keep = sp in passed_anywhere
                    if keep:
                        n_records_out += 1
                        out.write(line)
                elif keep:
                    out.write(line)
        sys.stderr.write(
            f"primary_rbh_filter: FASTA: {n_records_in} in, "
            f"{n_records_out} out\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
