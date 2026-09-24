#!/usr/bin/env python3
"""
Combine per-tool ortholog hits into consensus calls, then de-duplicate by
species using the configured tie-break.

A "hit" is one candidate ortholog found by one tool for one query group.
All hits are normalised to a common identifier space (UniProt accession)
*before* they reach this module, so foldseek's AFDB ids must already be
translated upstream. Each hit carries the metadata needed for both the
consensus vote and the dedup tie-break.

This module is pure (no I/O, no network) so it is fully unit-testable.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Iterable


@dataclass
class Hit:
    accession: str            # UniProt accession (normalised id space)
    tool: str                 # 'mmseqs2' | 'hmmer' | 'foldseek'
    score: float              # higher = better (bitscore or TM-score*1000)
    evalue: float             # lower = better; foldseek may set 0.0
    species: str              # de-dup species key (already collapsed)
    taxon_id: int             # NCBI taxon id; -1 if unknown
    reviewed: bool            # Swiss-Prot reviewed?
    seq_len: int              # sequence length, for the 'longest' tie-break


@dataclass
class ConsensusCall:
    accession: str
    tools: set[str] = field(default_factory=set)
    best_score: float = float("-inf")
    best_evalue: float = float("inf")
    species: str = ""
    taxon_id: int = -1
    reviewed: bool = False
    seq_len: int = 0

    def absorb(self, hit: Hit) -> None:
        self.tools.add(hit.tool)
        if hit.score > self.best_score:
            self.best_score = hit.score
        if hit.evalue < self.best_evalue:
            self.best_evalue = hit.evalue
        # Metadata should be consistent across tools for the same accession;
        # take the most informative (reviewed wins, real taxon id wins).
        self.reviewed = self.reviewed or hit.reviewed
        if self.taxon_id == -1 and hit.taxon_id != -1:
            self.taxon_id = hit.taxon_id
        if not self.species and hit.species:
            self.species = hit.species
        self.seq_len = max(self.seq_len, hit.seq_len)


def combine_hits(hits: Iterable[Hit], *, mode: str, min_tools: int,
                 enabled_tools: list[str]) -> list[ConsensusCall]:
    """
    Group hits by accession and apply the consensus rule.
      mode='union'        -> keep any accession found by >=1 enabled tool
      mode='intersection' -> keep only accessions found by ALL enabled tools
      mode='consensus'    -> keep accessions found by >= min_tools tools
    Returns the surviving ConsensusCalls (one per accession).
    """
    by_acc: dict[str, ConsensusCall] = {}
    for h in hits:
        call = by_acc.setdefault(h.accession, ConsensusCall(accession=h.accession))
        call.absorb(h)

    n_enabled = len(enabled_tools)
    if mode == "union":
        threshold = 1
    elif mode == "intersection":
        threshold = n_enabled
    elif mode == "consensus":
        threshold = min_tools
    else:
        raise ValueError(f"unknown consensus mode: {mode}")

    return [c for c in by_acc.values() if len(c.tools) >= threshold]


# --------------------------------------------------------------------------
# Species de-duplication
# --------------------------------------------------------------------------

def _tie_key(call: ConsensusCall, rule: str):
    """
    Build a sort key such that the BEST representative sorts FIRST
    (ascending sort, so we negate 'higher-is-better' quantities).
    """
    # Components, smaller = better after this transform:
    not_reviewed = 0 if call.reviewed else 1     # reviewed first
    neg_score = -call.best_score                  # higher score first
    evalue = call.best_evalue                     # lower e-value first
    taxid = call.taxon_id if call.taxon_id >= 0 else float("inf")  # lowest first
    neg_len = -call.seq_len                        # longer first

    if rule == "swissprot_then_lowest_taxid":
        return (not_reviewed, taxid, evalue, neg_score)
    if rule == "swissprot_then_best_score":
        return (not_reviewed, neg_score, evalue, taxid)
    if rule == "best_score":
        return (neg_score, evalue, not_reviewed, taxid)
    if rule == "longest":
        return (neg_len, not_reviewed, neg_score, taxid)
    if rule == "lowest_taxid":
        return (taxid, not_reviewed, neg_score, evalue)
    raise ValueError(f"unknown tie_break rule: {rule}")


def dedup_by_species(calls: list[ConsensusCall], *, enabled: bool,
                     tie_break: str) -> list[ConsensusCall]:
    """
    Collapse calls so each species key appears once, keeping the best
    representative per the tie-break. If enabled is False, returns calls
    unchanged (but still sorted deterministically for reproducibility).
    Calls with an empty species key are never merged with each other.
    """
    if not enabled:
        return sorted(calls, key=lambda c: c.accession)

    best: dict[str, ConsensusCall] = {}
    passthrough: list[ConsensusCall] = []
    for call in calls:
        if not call.species:
            # No species info -> cannot dedup; keep as-is.
            passthrough.append(call)
            continue
        cur = best.get(call.species)
        if cur is None or _tie_key(call, tie_break) < _tie_key(cur, tie_break):
            best[call.species] = call

    result = list(best.values()) + passthrough
    # Deterministic output order.
    return sorted(result, key=lambda c: (c.species or "~", c.accession))
