#!/usr/bin/env python3
"""Align one ortholog group (MAFFT or PAGAN2) per the config."""

from __future__ import annotations

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import load_config                                  # noqa: E402
from align_tree import run_mafft, run_pagan2                     # noqa: E402
from supermatrix import read_fasta                               # noqa: E402


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--in", dest="inp", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--work", required=True)
    ap.add_argument("--threads", type=int, default=8)
    args = ap.parse_args()

    cfg = load_config(args.config)
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)

    seqs = read_fasta(args.inp)
    # Single-sequence (or empty) groups can't be aligned; pass through so the
    # supermatrix step can still place the column(s). MAFFT handles 1 seq, but
    # we short-circuit to avoid spurious tool invocations on empties.
    if len(seqs) == 0:
        open(args.out, "w").close()
        sys.stderr.write(f"WARN: {args.inp} empty; wrote empty alignment\n")
        return 0
    if len(seqs) == 1:
        # Nothing to align; emit as-is.
        with open(args.out, "w") as fh:
            for name, seq in seqs.items():
                fh.write(f">{name}\n{seq}\n")
        sys.stderr.write(f"NOTE: {args.inp} has 1 sequence; copied through\n")
        return 0

    if cfg["aligner"] == "mafft":
        run_mafft(args.inp, args.out, cfg, args.threads)
    else:
        run_pagan2(args.inp, args.out, cfg, args.threads, args.work)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
