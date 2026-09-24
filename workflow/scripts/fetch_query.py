#!/usr/bin/env python3
"""
Write the query input file for one UniProt accession:
  --kind fasta      -> a single-sequence FASTA (for mmseqs2 / phmmer)
  --kind structure  -> the AlphaFold PDB    (for foldseek)

Used by the queries.smk rules. Network access goes through fetch.py, which
caches, so re-runs are cheap.
"""

from __future__ import annotations

import argparse
import os
import shutil
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import load_config                                  # noqa: E402
from fetch import fetch_entry, fetch_structure, FetchError      # noqa: E402


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--accession", required=True)
    ap.add_argument("--kind", required=True, choices=["fasta", "structure"])
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    cfg = load_config(args.config)
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)

    if args.kind == "fasta":
        entry = fetch_entry(args.accession, cfg["cache_dir"])
        with open(args.out, "w") as fh:
            fh.write(f">{entry.accession} {entry.scientific_name}\n")
            seq = entry.sequence
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i + 60] + "\n")
    else:
        path = fetch_structure(args.accession, cfg["cache_dir"])
        if path is None:
            sys.stderr.write(
                f"ERROR: no AlphaFold structure for {args.accession}; "
                f"foldseek cannot query it.\n")
            return 1
        shutil.copyfile(path, args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
