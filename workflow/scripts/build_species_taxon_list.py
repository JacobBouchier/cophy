#!/usr/bin/env python3
"""
Build the species taxon ID list used to restrict cooccur RBH searches.

Reads one or more primary members.tsv files and writes the UNION of
their taxon_ids (column 3) to a flat text file (one per line). This
list defines which species the cooccur RBH searches will check.

Invocation:
  build_species_taxon_list.py --primary-members A.tsv B.tsv --out file.txt

Output file format (consumed by build_cooccurrence_rbh.py):
  # 1234 species (union across primary members)
  # generated 2026-06-05
  562
  1280
  ...
"""

from __future__ import annotations

import argparse
import os
import sys
import time


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--primary-members", nargs="+", required=True,
                    help="one or more members.tsv files (cols include "
                         "species, accession, taxon_id, ...)")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    taxon_ids: set[str] = set()
    for path in args.primary_members:
        if not os.path.exists(path):
            sys.stderr.write(f"WARN: missing members file: {path}\n")
            continue
        with open(path) as fh:
            header = fh.readline().rstrip("\n").split("\t")
            try:
                col = header.index("taxon_id")
            except ValueError:
                sys.stderr.write(
                    f"WARN: {path}: no taxon_id column in header "
                    f"{header}\n")
                continue
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) <= col:
                    continue
                tid = f[col].strip()
                if tid and tid != "0":
                    taxon_ids.add(tid)

    sys.stderr.write(
        f"species_taxon_list: {len(taxon_ids)} unique taxon_ids across "
        f"{len(args.primary_members)} members files\n")

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    with open(args.out, "w") as fh:
        fh.write(f"# {len(taxon_ids)} species (union across primary members)\n")
        fh.write(f"# generated {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
        for tid in sorted(taxon_ids, key=lambda x: int(x) if x.isdigit()
                          else float("inf")):
            fh.write(f"{tid}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
