#!/usr/bin/env python3
"""Trim one aligned group with ClipKIT per the config."""

from __future__ import annotations

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import load_config                                  # noqa: E402
from align_tree import run_clipkit                               # noqa: E402
from supermatrix import read_fasta                               # noqa: E402


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--in", dest="inp", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    cfg = load_config(args.config)
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)

    # ClipKIT errors on empty input; handle the degenerate cases ourselves.
    seqs = read_fasta(args.inp) if os.path.getsize(args.inp) else {}
    if len(seqs) <= 1:
        # Nothing meaningful to trim; pass through (possibly empty).
        with open(args.out, "w") as fh:
            for name, seq in seqs.items():
                fh.write(f">{name}\n{seq}\n")
        sys.stderr.write(f"NOTE: {args.inp} has <=1 seq; copied without trim\n")
        return 0

    run_clipkit(args.inp, args.out, cfg)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
