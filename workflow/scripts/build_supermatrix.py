#!/usr/bin/env python3
"""
Assemble the final supermatrix from per-group trimmed alignments.

Each input is a trimmed aligned FASTA (headers = species keys). We read them,
build the concatenated matrix with gap padding for missing species, and write
the supermatrix FASTA plus a RAxML-style partition file.
"""

from __future__ import annotations

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import load_config                                  # noqa: E402
from supermatrix import (make_group, build_supermatrix,          # noqa: E402
                         read_fasta, write_fasta, write_partitions)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--alignments", nargs="+", required=True,
                    help="trimmed per-group aligned FASTAs")
    ap.add_argument("--out-fasta", required=True)
    ap.add_argument("--out-partitions", required=True)
    args = ap.parse_args()

    cfg = load_config(args.config)

    groups = []
    for path in sorted(args.alignments):  # deterministic group order
        name = os.path.splitext(os.path.basename(path))[0]
        seqs = read_fasta(path)
        if not seqs:
            sys.stderr.write(f"WARN: empty alignment {path}; skipping\n")
            continue
        groups.append(make_group(name, seqs))

    if not groups:
        sys.stderr.write("ERROR: no non-empty alignments to concatenate\n")
        return 1

    matrix, partitions = build_supermatrix(groups)
    write_fasta(matrix, args.out_fasta)

    model = "GTR" if cfg["seq_type"] == "nt" else "LG"
    write_partitions(partitions, args.out_partitions, model=model)

    total = len(next(iter(matrix.values()))) if matrix else 0
    sys.stderr.write(f"supermatrix: {len(matrix)} species x {total} cols, "
                     f"{len(groups)} partitions\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
