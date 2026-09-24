#!/usr/bin/env python3
"""
Alignment / trimming / tree command builders + runners.

Stages covered:
  * MAFFT alignment (per ortholog group)
  * PAGAN2 alignment with an auto-built MAFFT+FastTree2 guide tree
  * ClipKIT trimming (per group)
  * FastTree2 tree inference (on the trimmed supermatrix)

Command builders are pure and unit-tested. Runners are thin subprocess calls.
"""

from __future__ import annotations

import os
import shlex
import subprocess
import sys


# --------------------------------------------------------------------------
# Command builders (pure)
# --------------------------------------------------------------------------

def cmd_mafft(in_fasta: str, cfg: dict, threads: int) -> list[str]:
    algo = cfg["mafft_algo"]
    args = ["mafft", "--thread", str(threads)]
    if algo == "auto":
        args.append("--auto")
    elif algo == "linsi":
        args += ["--localpair", "--maxiterate", str(cfg["mafft_maxiterate"])]
    elif algo == "ginsi":
        args += ["--globalpair", "--maxiterate", str(cfg["mafft_maxiterate"])]
    elif algo == "einsi":
        args += ["--genafpair", "--maxiterate", str(cfg["mafft_maxiterate"])]
    elif algo == "fftns":
        args += ["--retree", "2"]
    if cfg["seq_type"] == "nt":
        args.append("--nuc")
    else:
        args.append("--amino")
    if cfg["mafft_extra"].strip():
        args += shlex.split(cfg["mafft_extra"])
    args.append(in_fasta)
    return args


def cmd_fasttree(in_fasta: str, cfg: dict, *, model_override: str | None = None,
                 gamma_override: bool | None = None) -> list[str]:
    args = [cfg["fasttree_bin"]]
    model = model_override if model_override is not None else cfg["fasttree_model"]
    if cfg["seq_type"] == "nt":
        args.append("-nt")
        if model == "gtr":
            args.append("-gtr")
        # 'jc' is FastTree's default for -nt; no extra flag needed.
    else:
        if model == "lg":
            args.append("-lg")
        elif model == "wag":
            args.append("-wag")
        # 'jtt' is FastTree's default protein model; no flag needed.
    gamma = cfg["fasttree_gamma"] if gamma_override is None else gamma_override
    if gamma:
        args.append("-gamma")
    if cfg["fasttree_extra"].strip():
        args += shlex.split(cfg["fasttree_extra"])
    args.append(in_fasta)
    return args


def cmd_pagan2(in_fasta: str, guide_tree: str, out_prefix: str,
               cfg: dict) -> list[str]:
    # PAGAN2 aligns both protein and DNA from --seqfile without a type flag.
    # For codon-aware alignment of in-frame CDS, add '--codons' via
    # [align] pagan2_extra in the config.
    args = ["pagan2", "--seqfile", in_fasta, "--treefile", guide_tree,
            "--outfile", out_prefix]
    if cfg["pagan2_extra"].strip():
        args += shlex.split(cfg["pagan2_extra"])
    return args


def cmd_clipkit(in_fasta: str, out_fasta: str, cfg: dict) -> list[str]:
    args = ["clipkit", in_fasta, "-o", out_fasta, "-m", cfg["clipkit_mode"]]
    if "gappy" in cfg["clipkit_mode"]:
        args += ["-g", str(cfg["clipkit_gaps"])]
    if cfg["clipkit_extra"].strip():
        args += shlex.split(cfg["clipkit_extra"])
    return args


# --------------------------------------------------------------------------
# Runners
# --------------------------------------------------------------------------

def _run(cmd: list[str], *, stdout_path: str | None = None) -> None:
    sys.stderr.write("RUN: " + " ".join(shlex.quote(c) for c in cmd) + "\n")
    if stdout_path:
        with open(stdout_path, "w") as out:
            subprocess.run(cmd, check=True, stdout=out)
    else:
        subprocess.run(cmd, check=True)


def run_mafft(in_fasta: str, out_fasta: str, cfg: dict, threads: int) -> None:
    # MAFFT writes the alignment to stdout.
    _run(cmd_mafft(in_fasta, cfg, threads), stdout_path=out_fasta)


def run_guide_tree(in_fasta: str, out_tree: str, cfg: dict, threads: int,
                   workdir: str) -> None:
    """Build a guide tree for PAGAN2: MAFFT align -> FastTree2."""
    os.makedirs(workdir, exist_ok=True)
    guide_aln = os.path.join(workdir, "guide_mafft.fasta")
    run_mafft(in_fasta, guide_aln, cfg, threads)
    # Guide tree always uses default model/gamma off for speed.
    _run(cmd_fasttree(guide_aln, cfg, gamma_override=False),
         stdout_path=out_tree)


def run_pagan2(in_fasta: str, out_aln: str, cfg: dict, threads: int,
               workdir: str) -> None:
    os.makedirs(workdir, exist_ok=True)
    guide = os.path.join(workdir, "guide.nwk")
    run_guide_tree(in_fasta, guide, cfg, threads, workdir)
    prefix = os.path.join(workdir, "pagan_out")
    _run(cmd_pagan2(in_fasta, guide, prefix, cfg))
    # PAGAN2 writes <prefix>.fas ; normalise to out_aln.
    produced = prefix + ".fas"
    if not os.path.exists(produced):
        raise FileNotFoundError(f"PAGAN2 did not produce {produced}")
    os.replace(produced, out_aln)


def run_clipkit(in_fasta: str, out_fasta: str, cfg: dict) -> None:
    _run(cmd_clipkit(in_fasta, out_fasta, cfg))


def run_fasttree(in_fasta: str, out_tree: str, cfg: dict) -> None:
    _run(cmd_fasttree(in_fasta, cfg), stdout_path=out_tree)
