#!/usr/bin/env python3
"""
Build a concatenated supermatrix from per-ortholog-group alignments.

Input: one aligned FASTA per ortholog group, where each sequence header is
the *species key* (so rows line up across groups). Each group may be missing
some species; those cells are padded with gaps for the full alignment width
of that group. The output is a single FASTA with one row per species, plus a
partition file (RAxML-style) recording each group's column range.

Pure functions over in-memory dicts; the thin FASTA I/O is at the bottom and
is exercised by the tests too.
"""

from __future__ import annotations

import os
from dataclasses import dataclass


@dataclass
class Group:
    name: str
    # species_key -> aligned sequence (all equal length within a group)
    seqs: dict[str, str]
    width: int


def _aln_width(seqs: dict[str, str], name: str) -> int:
    widths = {len(s) for s in seqs.values()}
    if not widths:
        raise ValueError(f"group {name!r}: no sequences")
    if len(widths) != 1:
        raise ValueError(
            f"group {name!r}: aligned sequences are not all the same length "
            f"(widths seen: {sorted(widths)})")
    return widths.pop()


def make_group(name: str, seqs: dict[str, str]) -> Group:
    return Group(name=name, seqs=dict(seqs), width=_aln_width(seqs, name))


def build_supermatrix(groups: list[Group], *, gap_char: str = "-"
                      ) -> tuple[dict[str, str], list[tuple[str, int, int]]]:
    """
    Concatenate groups into a supermatrix.
    Returns (matrix, partitions) where:
      matrix: species_key -> concatenated aligned string (all equal length)
      partitions: list of (group_name, start_1based, end_1based)
    Species present in any group get a row; missing cells are gap-filled.
    Group order is preserved as given (caller should sort for determinism).
    """
    # Union of all species, in first-seen order for stable output.
    species_order: list[str] = []
    seen: set[str] = set()
    for g in groups:
        for sp in g.seqs:
            if sp not in seen:
                seen.add(sp)
                species_order.append(sp)

    matrix: dict[str, list[str]] = {sp: [] for sp in species_order}
    partitions: list[tuple[str, int, int]] = []
    cursor = 0
    for g in groups:
        start = cursor + 1
        for sp in species_order:
            seq = g.seqs.get(sp)
            if seq is None:
                seq = gap_char * g.width
            matrix[sp].append(seq)
        cursor += g.width
        partitions.append((g.name, start, cursor))

    joined = {sp: "".join(parts) for sp, parts in matrix.items()}

    # Sanity: every row must equal the total width.
    total = cursor
    for sp, s in joined.items():
        if len(s) != total:
            raise ValueError(
                f"internal error: row {sp!r} length {len(s)} != {total}")
    return joined, partitions


# --------------------------------------------------------------------------
# Minimal FASTA I/O
# --------------------------------------------------------------------------

def read_fasta(path: str) -> dict[str, str]:
    seqs: dict[str, str] = {}
    name: str | None = None
    chunks: list[str] = []
    with open(path, "r", encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(chunks)
                # Header up to first whitespace is the id (the species key).
                name = line[1:].split()[0] if line[1:].strip() else line[1:]
                chunks = []
            elif line:
                chunks.append(line.strip())
    if name is not None:
        seqs[name] = "".join(chunks)
    return seqs


def write_fasta(seqs: dict[str, str], path: str, *, width: int = 60) -> None:
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        for name, seq in seqs.items():
            fh.write(f">{name}\n")
            for i in range(0, len(seq), width):
                fh.write(seq[i:i + width] + "\n")


def write_partitions(partitions: list[tuple[str, int, int]], path: str,
                     *, model: str) -> None:
    """RAxML-style partition file: '<model>, <name> = <start>-<end>'."""
    with open(path, "w", encoding="utf-8") as fh:
        for name, start, end in partitions:
            fh.write(f"{model}, {name} = {start}-{end}\n")
