#!/usr/bin/env python3
"""Infer the species tree from the trimmed supermatrix with FastTree2."""

from __future__ import annotations

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import load_config                                  # noqa: E402
from align_tree import run_fasttree                              # noqa: E402


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--in", dest="inp", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    cfg = load_config(args.config)
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    run_fasttree(args.inp, args.out, cfg)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
