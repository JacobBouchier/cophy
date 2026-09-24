#!/usr/bin/env python3
"""
Run the enabled search tools for one query UniProt id against the configured
databases, parse + normalise the hits, and write a per-(group, tool) TSV of
RawHits. A later rule combines these across tools (consensus) and species.

The command-building functions (`cmd_*`) are pure and unit-tested; the
execution wrapper is a thin subprocess call. Nothing here needs the network.

Usage (called by Snakemake, but runnable standalone):
    search_run.py --config config.txt --group P0A7G6 --tool mmseqs2 \
        --query-fasta q.fasta --out hits.tsv --tmp /scratch/tmp --threads 8
"""

from __future__ import annotations

import argparse
import os
import shlex
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import load_config                       # noqa: E402
from parse_search import (parse_blasttab,            # noqa: E402
                          parse_phmmer_tblout, RawHit)


# --------------------------------------------------------------------------
# Command builders (pure, testable)
# --------------------------------------------------------------------------

def cmd_mmseqs2(query_fasta: str, db_path: str, out_m8: str, tmp: str,
                cfg: dict, threads: int) -> list[str]:
    p = cfg["mmseqs2"]
    cmd = [
        "mmseqs", "easy-search", query_fasta, db_path, out_m8, tmp,
        "-s", str(p["sensitivity"]),
        "-e", str(p["evalue"]),
        "--min-seq-id", str(p["min_seq_id"]),
        "-c", str(p["coverage"]),
        "--cov-mode", str(p["cov_mode"]),
        "--threads", str(threads),
        "--format-mode", "0",  # BLAST-tab .m8
    ]
    # Cap memory so a large target DB is split into chunks that fit. A value
    # of "0" (or empty) means unlimited (mmseqs2 default).
    smem = str(p.get("split_memory", "0")).strip()
    if smem and smem != "0":
        cmd += ["--split-memory-limit", smem]
    return cmd


def cmd_phmmer(query_fasta: str, db_fasta: str, tblout: str,
               cfg: dict, threads: int) -> list[str]:
    p = cfg["hmmer"]
    args = ["phmmer", "--tblout", tblout, "--cpu", str(threads),
            "-E", str(p["evalue"]), "--domE", str(p["dom_evalue"]),
            "--incE", str(p["incE"])]
    if p["bitscore"] and p["bitscore"] > 0:
        args += ["-T", str(p["bitscore"])]
    args += [query_fasta, db_fasta]
    return args


def cmd_foldseek(query_struct: str, db_path: str, out_m8: str, tmp: str,
                 cfg: dict, threads: int) -> list[str]:
    p = cfg["foldseek"]
    fmt = ("query,target,fident,alnlen,mismatch,gapopen,"
           "qstart,qend,tstart,tend,evalue,bits,alntmscore")
    args = [
        "foldseek", "easy-search", query_struct, db_path, out_m8, tmp,
        "-e", str(p["evalue"]),
        "-c", str(p["coverage"]),
        "--threads", str(threads),
        "--format-mode", "0",
        "--format-output", fmt,
    ]
    # Cap memory so the DB gets split into chunks that fit in available RAM.
    # The full AFDB/UniProt index wants ~300+ GB by default for one split;
    # capping to e.g. 100G forces foldseek into more splits (slower but
    # avoids OOM on a 120G batch node).
    if p.get("split_memory"):
        args += ["--split-memory-limit", str(p["split_memory"])]
    # Sensitivity controls: more candidates through prefilter = better recall
    # for distant homologs at the cost of slower searches. max_seqs caps the
    # prefilter output; exhaustive_search skips the prefilter entirely.
    if p.get("max_seqs"):
        args += ["--max-seqs", str(p["max_seqs"])]
    if p.get("exhaustive_search"):
        args += ["--exhaustive-search", "1"]
    if p["alntmscore"]:
        # alignment-type 1 = TMalign, needed for meaningful alntmscore
        args += ["--alignment-type", "1"]
    # Raw passthrough flags, appended verbatim (mirrors mafft_extra etc.)
    if p.get("extra"):
        import shlex
        args += shlex.split(str(p["extra"]))
    return args


# --------------------------------------------------------------------------
# Execution + parsing
# --------------------------------------------------------------------------

def _run(cmd: list[str]) -> None:
    sys.stderr.write("RUN: " + " ".join(shlex.quote(c) for c in cmd) + "\n")
    subprocess.run(cmd, check=True)


def run_tool(tool: str, group: str, query_path: str, cfg: dict,
             out_tsv: str, tmp: str, threads: int) -> list[RawHit]:
    os.makedirs(tmp, exist_ok=True)
    raw_out = out_tsv + f".{tool}.raw"

    if tool == "mmseqs2":
        _run(cmd_mmseqs2(query_path, cfg["mmseqs2_db"], raw_out, tmp,
                         cfg, threads))
        with open(raw_out) as fh:
            hits = parse_blasttab(fh.read(), group=group, tool="mmseqs2")
    elif tool == "hmmer":
        _run(cmd_phmmer(query_path, cfg["hmmer_db"], raw_out, cfg, threads))
        with open(raw_out) as fh:
            hits = parse_phmmer_tblout(fh.read(), group=group)
    elif tool == "foldseek":
        _run(cmd_foldseek(query_path, cfg["foldseek_db_path"], raw_out, tmp,
                          cfg, threads))
        with open(raw_out) as fh:
            hits = parse_blasttab(fh.read(), group=group, tool="foldseek",
                                  tmscore_col=12)
    else:
        raise ValueError(f"unknown tool: {tool}")

    return hits


def write_hits_tsv(hits: list[RawHit], path: str) -> None:
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "w") as fh:
        fh.write("group\ttarget\ttool\tscore\tevalue\n")
        for h in hits:
            fh.write(f"{h.group}\t{h.target}\t{h.tool}\t{h.score}\t{h.evalue}\n")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--group", required=True)
    ap.add_argument("--tool", required=True,
                    choices=["mmseqs2", "hmmer", "foldseek"])
    ap.add_argument("--query", required=True,
                    help="query FASTA (seq tools) or structure (foldseek)")
    ap.add_argument("--out", required=True)
    ap.add_argument("--tmp", required=True)
    ap.add_argument("--threads", type=int, default=8)
    args = ap.parse_args()

    cfg = load_config(args.config)
    hits = run_tool(args.tool, args.group, args.query, cfg,
                    args.out, args.tmp, args.threads)
    write_hits_tsv(hits, args.out)
    sys.stderr.write(f"{args.tool}/{args.group}: {len(hits)} hits -> {args.out}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
