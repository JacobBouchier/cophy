#!/usr/bin/env python3
"""
Gene-neighborhood distance computation, nucleotide-based.

A *gene position* identifies one protein-coding gene's location in a
genome by the contig it sits on and its start coordinate. Distance
between two positions is the difference of start coordinates, in
nucleotides, with the sign indicating the cooccur protein's direction
relative to the primary protein along the contig.

    GenePos(assembly, contig, start, strand)

* assembly: opaque key used to constrain matches to one genome at a
  time. When coordinates come from NCBI IPG, we use the contig accession
  as the assembly key — IPG lists one row per (protein, assembly) and we
  match within a single row.
* contig: chromosome / replicon / contig id (must match between the
  primary and cooccur positions for them to be on the same molecule).
* start: 1-based start coordinate on the forward strand of the contig.
* strand: '+' or '-' (used only for sign convention; magnitude is always
  the absolute coordinate delta).

For one (species, primary_acc, cooccur_acc) call:
  * a species can have many primary-acc-mapped GenePos entries (one per
    assembly the primary hit was observed in)
  * same for the cooccur protein
  * we report the SIGNED nucleotide distance of the closest pair
    (smallest |delta|) across all available assembly/contig matches

Pure functions; no I/O. The fetch layer is in fetch_neighborhood.py.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class GenePos:
    assembly: str   # genome assembly key (we use contig accession from IPG)
    contig: str     # chromosome / replicon / contig id
    start: int      # 1-based start coordinate, forward-strand of contig
    strand: str     # '+' or '-'


def _delta(a: GenePos, b: GenePos) -> int | None:
    """
    Signed nucleotide distance from a -> b iff they share (assembly, contig).

    Sign convention: positive when b's start is downstream of a's start in
    the genomic sense relative to a's strand:
      * a on '+' strand: positive means b.start > a.start
      * a on '-' strand: positive means b.start < a.start (so 'downstream
        of a' along its reading direction)

    This is useful biologically: a positive distance means the cooccur
    protein is in the direction a's coding sequence reads, i.e. likely
    to be 'downstream' in operon terms when both are on the same strand.

    Returns None if the two positions don't share a contig.
    """
    if a.assembly != b.assembly or a.contig != b.contig:
        return None
    raw = b.start - a.start
    if a.strand == "-":
        raw = -raw
    return raw


def closest_signed_distance(primary_positions: list[GenePos],
                            cooccur_positions: list[GenePos]
                            ) -> int | None:
    """
    Find the pair (p, c) with the smallest |c.start - p.start| where p and
    c share (assembly, contig). Return the signed nucleotide delta of that
    pair, with sign convention as defined in _delta().

    Returns None if no pair shares a contig.
    """
    best: int | None = None
    for p in primary_positions:
        for c in cooccur_positions:
            d = _delta(p, c)
            if d is None:
                continue
            if best is None or abs(d) < abs(best):
                best = d
    return best


def format_distance(d: int | None, cap: int) -> str:
    """
    Render a signed nucleotide distance for the CSV.
      None        -> ''       (no shared contig / no data)
      |d| > cap   -> '>{cap}' (sign dropped at the cap, by convention)
      otherwise   -> '+N' or '-N' (or '0' if literally same coordinate)
    """
    if d is None:
        return ""
    if abs(d) > cap:
        return f">{cap}"
    if d == 0:
        return "0"
    return f"{d:+d}"


def format_distance_kb(d: int | None, cap: int) -> str:
    """
    Render an absolute kilobase distance for the CSV, for easier reading.
      None        -> ''
      |d| > cap   -> '>{cap_kb}'   (rounded down to nearest kb)
      otherwise   -> rounded to 1 decimal place if < 10 kb, else whole kb

    The kb column is always non-negative (sign lives in the bp column).
    """
    if d is None:
        return ""
    abs_d = abs(d)
    if abs_d > cap:
        return f">{cap // 1000}"
    kb = abs_d / 1000.0
    if kb < 10:
        return f"{kb:.1f}"
    return f"{int(round(kb))}"
