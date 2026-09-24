#!/usr/bin/env python3
"""
Parse the tabular outputs of mmseqs2, phmmer (HMMER), and foldseek into a
common list of (group, target_accession, tool, score, evalue) tuples.

These parsers are pure text->data functions, unit-tested against captured
sample lines, so the brittle column-counting lives in one tested place.

Identifier normalisation
-------------------------
Everything must end up keyed by UniProt accession so the three tools can be
combined:
  * mmseqs2 / phmmer search UniRef90/50 or a UniProt FASTA. Target ids look
    like 'UniRef90_P12345' or 'sp|P12345|NAME' or bare 'P12345'. We extract
    the accession.
  * foldseek searches AFDB; targets look like 'AF-P12345-F1-model_v4'. We
    extract 'P12345'.
"""

from __future__ import annotations

import re
from dataclasses import dataclass


@dataclass
class RawHit:
    group: str          # the query UniProt id whose group this hit belongs to
    target: str         # normalised UniProt accession of the hit
    tool: str
    score: float        # higher = better
    evalue: float       # lower = better


# --------------------------------------------------------------------------
# Identifier normalisation
# --------------------------------------------------------------------------

_UNIREF = re.compile(r"^UniRef(?:100|90|50)_(.+)$")
_SP_TR = re.compile(r"^(?:sp|tr)\|([^|]+)\|")
_AFDB = re.compile(r"^AF-([A-Za-z0-9]+)-F\d+-model")
# UniProt accession general shape (canonical + isoform suffix tolerated).
_ACC = re.compile(r"^[A-NR-Z0-9][A-Z0-9]{5,9}(?:-\d+)?$", re.IGNORECASE)


def normalise_target(raw: str) -> str:
    """Reduce any of the supported target id forms to a UniProt accession."""
    raw = raw.strip()
    m = _AFDB.match(raw)
    if m:
        return m.group(1)
    m = _UNIREF.match(raw)
    if m:
        raw = m.group(1)
    m = _SP_TR.match(raw)
    if m:
        return m.group(1)
    # Strip an isoform suffix's trailing junk only if clearly an accession.
    return raw


# --------------------------------------------------------------------------
# mmseqs2 / foldseek share BLAST-tab (.m8) layout
#   query target fident alnlen mismatch gapopen qstart qend tstart tend evalue bits
# foldseek (our call) appends alntmscore as a 13th column.
# --------------------------------------------------------------------------

def parse_blasttab(text: str, *, group: str, tool: str,
                   tmscore_col: int | None = None) -> list[RawHit]:
    """
    Parse BLAST-tab output. For foldseek we pass tmscore_col=12 (0-based) so
    the TM-score can be promoted to the 'score' (scaled to be comparable as
    higher=better); otherwise 'bits' (col 11) is the score.
    """
    hits: list[RawHit] = []
    for line in text.splitlines():
        line = line.rstrip("\n")
        if not line or line.startswith("#"):
            continue
        f = line.split("\t")
        if len(f) < 12:
            continue
        target = normalise_target(f[1])
        try:
            evalue = float(f[10])
            bits = float(f[11])
        except ValueError:
            continue
        if tmscore_col is not None and len(f) > tmscore_col:
            try:
                tm = float(f[tmscore_col])
                score = tm * 1000.0  # scale TM-score into a bit-like range
            except ValueError:
                score = bits
        else:
            score = bits
        hits.append(RawHit(group=group, target=target, tool=tool,
                           score=score, evalue=evalue))
    return hits


# --------------------------------------------------------------------------
# phmmer --tblout : whitespace-delimited, fixed leading columns
#   target_name accession query_name accession  full_E full_score full_bias ...
#   (col0)      (col1)    (col2)     (col3)      (col4) (col5)      (col6)
# --------------------------------------------------------------------------

def parse_phmmer_tblout(text: str, *, group: str) -> list[RawHit]:
    hits: list[RawHit] = []
    for line in text.splitlines():
        if not line or line.startswith("#"):
            continue
        f = line.split()
        if len(f) < 6:
            continue
        target = normalise_target(f[0])
        try:
            evalue = float(f[4])     # full sequence E-value
            score = float(f[5])      # full sequence bit score
        except ValueError:
            continue
        hits.append(RawHit(group=group, target=target, tool="hmmer",
                           score=score, evalue=evalue))
    return hits
