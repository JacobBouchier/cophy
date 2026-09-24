#!/usr/bin/env python3
"""
Reciprocal-best-hit (RBH) algorithm for species-restricted cooccurrence
detection.

Given:
  * a cooccur query Q (UniProt accession, with a FASTA and optionally a PDB)
  * Q's home species' proteome (FASTA) and optionally structures
  * one tree species S's proteome (FASTA) and optionally structures

For each search tool (mmseqs2, foldseek), this module:
  1. Forward:  Q vs S's proteins  -> best hit P_S in S
  2. Reverse:  P_S vs Q's home proteins -> best hit P_H in home species
  3. Confirm:  P_H == Q  -> reciprocal best hit
                else      -> reject

Returns RBH calls per tool. The caller (build_cooccurrence_rbh.py)
combines results across tools via the configured consensus mode.

The actual search invocations are delegated to small subprocess wrappers
- this module's job is the algorithm, not the tool commands.
"""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
import tempfile
from dataclasses import dataclass


@dataclass
class RBHResult:
    """Outcome of an RBH check for one (query, target_species, tool) triple."""
    query: str                  # the cooccur query accession (Q)
    target_species: str         # taxon_id of the species being checked
    tool: str                   # 'mmseqs2' or 'foldseek'
    forward_hit: str | None     # best hit in target species, if any
    reverse_hit: str | None     # best hit in home species during reverse, if any
    reciprocal: bool            # True iff reverse_hit == query (or skipped)
    forward_evalue: float | None = None
    forward_score: float | None = None
    skipped_reverse: bool = False     # True iff forward was so strong we
                                      # skipped the reverse search


# E-value threshold below which we skip the reverse search.
# Any forward hit with e-value <= this is auto-called reciprocal.
# Set to 0 to always do the reverse search (strict RBH).
STRONG_HIT_SKIP_REVERSE = 1e-50


# --------------------------------------------------------------------------
# Pure algorithm: given parsed BLAST-tab style output, find the best hit
# --------------------------------------------------------------------------

def best_hit_from_blasttab(text: str,
                           query_target_self: str | None = None,
                           ) -> tuple[str, float, float] | None:
    """
    Parse a BLAST-tab search output and return the best non-self hit
    as (target_acc, evalue, bitscore). Returns None if no hits.

    BLAST-tab columns (foldseek/mmseqs2 default):
      query  target  fident  alnlen  mismatch  gapopen
      qstart qend  tstart  tend  evalue  bits  [alntmscore]

    `query_target_self`: if given, an exact target accession to exclude
    from consideration (used to filter self-hits when the query and DB
    overlap). For our case Q vs S has no overlap, so usually None.

    Tie-breaking: lowest evalue wins. If evalues tie, highest bitscore.
    """
    best: tuple[str, float, float] | None = None
    for line in text.splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        f = line.split("\t")
        if len(f) < 12:
            continue
        target = f[1].strip()
        if not target:
            continue
        if query_target_self and _accession_match(target, query_target_self):
            continue
        # The target field can be a full FASTA header in some tool modes.
        # Extract just the UniProt accession.
        target = _normalise_accession(target)
        try:
            evalue = float(f[10])
            bits = float(f[11])
        except ValueError:
            continue
        if best is None or (evalue, -bits) < (best[1], -best[2]):
            best = (target, evalue, bits)
    return best


def _normalise_accession(target: str) -> str:
    """Extract UniProt accession from common target id formats.

    Examples:
        'P12345'                                 -> 'P12345'
        'sp|P12345|TEST_ECOLI'                   -> 'P12345'
        'tr|A0A123|whatever'                     -> 'A0A123'
        'AF-P12345-F1-model_v4'                  -> 'P12345'
        'AF-Q18B78-F1-model_v6.pdb.gz'           -> 'Q18B78'
        'UniRef90_P12345'                        -> 'P12345'
    """
    t = target.strip()
    # Strip .pdb / .pdb.gz suffix
    if t.endswith(".pdb.gz"):
        t = t[:-7]
    elif t.endswith(".pdb"):
        t = t[:-4]
    if t.startswith("AF-"):
        # AF-<ACC>-F1-model_v<N>
        parts = t.split("-")
        if len(parts) >= 2:
            return parts[1]
    if t.startswith(("sp|", "tr|")):
        parts = t.split("|")
        if len(parts) >= 2:
            return parts[1]
    if t.startswith("UniRef"):
        # UniRef100_X, UniRef90_X, UniRef50_X
        if "_" in t:
            return t.split("_", 1)[1]
    return t


def _accession_match(a: str, b: str) -> bool:
    """True if two ids refer to the same accession after normalisation."""
    return _normalise_accession(a) == _normalise_accession(b)


# --------------------------------------------------------------------------
# Subprocess wrappers for the actual searches.
# These are kept thin so they're easy to stub in tests.
# --------------------------------------------------------------------------

def _run_mmseqs_one_vs_db(query_fasta: str, target_db_fasta: str,
                          out_tsv: str, tmpdir: str,
                          evalue: float = 1e-3,
                          threads: int = 4) -> None:
    """Run mmseqs2 easy-search with output going to out_tsv (BLAST-tab)."""
    cmd = [
        "mmseqs", "easy-search", query_fasta, target_db_fasta,
        out_tsv, tmpdir,
        "-e", str(evalue),
        "--threads", str(threads),
        "--format-mode", "0",
        # default format-output matches the BLAST-tab we parse.
        "-s", "7.5",
        "--remove-tmp-files", "1",
    ]
    subprocess.run(cmd, check=True, capture_output=True)


def _run_mmseqs_vs_indexed_db(query_fasta: str, target_db_path: str,
                              out_tsv: str, tmpdir: str,
                              evalue: float = 1e-3,
                              threads: int = 4) -> None:
    """Run mmseqs2 search using a pre-built indexed target DB.

    Skips the createdb step that easy-search does on every call,
    which is ~1-2 seconds of overhead avoided per RBH check.
    """
    # Create query DB on the fly (small, one record - fast)
    os.makedirs(tmpdir, exist_ok=True)
    qdb = os.path.join(tmpdir, "qdb")
    res_db = os.path.join(tmpdir, "result")
    search_tmp = os.path.join(tmpdir, "search_tmp")
    os.makedirs(search_tmp, exist_ok=True)
    subprocess.run(
        ["mmseqs", "createdb", query_fasta, qdb, "--shuffle", "0"],
        check=True, capture_output=True)
    subprocess.run(
        ["mmseqs", "search", qdb, target_db_path, res_db, search_tmp,
         "-e", str(evalue), "--threads", str(threads),
         "-s", "7.5", "--remove-tmp-files", "1"],
        check=True, capture_output=True)
    # Convert result DB to BLAST-tab format
    subprocess.run(
        ["mmseqs", "convertalis", qdb, target_db_path, res_db, out_tsv,
         "--format-mode", "0", "--threads", str(threads)],
        check=True, capture_output=True)


def _run_foldseek_one_vs_db(query_struct: str, target_db_path: str,
                            out_tsv: str, tmpdir: str,
                            evalue: float = 1e-3,
                            threads: int = 4,
                            taxon_list: str | None = None) -> None:
    """Run foldseek easy-search with TM-score-friendly defaults.

    If `taxon_list` is given (comma-separated taxon IDs), restrict the
    search to those taxa via --taxon-list. Used to do species-restricted
    foldseek RBH against the full AFDB without fetching per-species PDBs.
    """
    cmd = [
        "foldseek", "easy-search", query_struct, target_db_path,
        out_tsv, tmpdir,
        "-e", str(evalue),
        "--threads", str(threads),
        "--format-mode", "0",
        "--format-output",
        "query,target,fident,alnlen,mismatch,gapopen,"
        "qstart,qend,tstart,tend,evalue,bits,alntmscore",
        "--remove-tmp-files", "1",
    ]
    if taxon_list:
        cmd += ["--taxon-list", taxon_list]
    subprocess.run(cmd, check=True, capture_output=True)


def _run_phmmer_one_vs_db(query_fasta: str, target_db_fasta: str,
                          out_tblout: str,
                          evalue: float = 1e-3,
                          threads: int = 4) -> None:
    """Run HMMER's phmmer (single-sequence search), output in --tblout
    format (whitespace-delimited tabular).

    Output format columns (relevant ones):
      target_name  accession  query_name  accession  full_E  full_score ...
    """
    cmd = [
        "phmmer",
        "--tblout", out_tblout,
        "-E", str(evalue),
        "--cpu", str(threads),
        "--noali", "--notextw",
        query_fasta, target_db_fasta,
    ]
    subprocess.run(cmd, check=True, capture_output=True)


def best_hit_from_phmmer_tblout(text: str) -> tuple[str, float, float] | None:
    """Parse phmmer --tblout output and return the best hit as
    (target_acc, evalue, bitscore). Returns None if no hits.

    Format (the columns we use): target_name (0), full_E (4), full_score (5).
    Sorted by E-value ascending in the file, but we don't rely on that;
    we scan all rows for the minimum.
    """
    best: tuple[str, float, float] | None = None
    for line in text.splitlines():
        if not line or line.startswith("#"):
            continue
        f = line.split()
        if len(f) < 6:
            continue
        target = _normalise_accession(f[0])
        try:
            evalue = float(f[4])
            bits = float(f[5])
        except ValueError:
            continue
        if best is None or (evalue, -bits) < (best[1], -best[2]):
            best = (target, evalue, bits)
    return best


# --------------------------------------------------------------------------
# Top-level RBH check (per query, per target species, per tool)
# --------------------------------------------------------------------------

def rbh_check_mmseqs2(query_accession: str,
                      query_fasta: str,
                      home_proteome_fasta: str,
                      target_proteome_fasta: str,
                      target_taxon_id: str,
                      threads: int = 4,
                      evalue: float = 1e-3,
                      target_db_path: str | None = None,
                      home_db_path: str | None = None,
                      ) -> RBHResult:
    """Run mmseqs2 forward + reverse search for one (Q, S) pair.

    Args:
        query_fasta: file containing just Q's sequence (single record).
        home_proteome_fasta: Q's home species' proteome (FASTA).
        target_proteome_fasta: species S's proteome (FASTA).
        target_db_path: optional pre-built mmseqs DB for S. If given,
            skips per-call createdb overhead.
        home_db_path: same for Q's home species.
    """
    tmpdir = tempfile.mkdtemp(prefix="rbh_mmseqs2_")
    try:
        forward_tsv = os.path.join(tmpdir, "fwd.tsv")
        if target_db_path:
            _run_mmseqs_vs_indexed_db(
                query_fasta, target_db_path, forward_tsv,
                os.path.join(tmpdir, "fwd_tmp"),
                evalue=evalue, threads=threads)
        else:
            _run_mmseqs_one_vs_db(
                query_fasta, target_proteome_fasta, forward_tsv,
                os.path.join(tmpdir, "fwd_tmp"),
                evalue=evalue, threads=threads)
        with open(forward_tsv) as fh:
            fwd_text = fh.read()
        fwd_best = best_hit_from_blasttab(fwd_text)
        if fwd_best is None:
            return RBHResult(query=query_accession,
                             target_species=target_taxon_id,
                             tool="mmseqs2",
                             forward_hit=None, reverse_hit=None,
                             reciprocal=False)

        forward_hit, fwd_evalue, fwd_score = fwd_best

        # Strong-hit short-circuit: skip the reverse search when the forward hit
        # is overwhelmingly strong (e-value <= STRONG_HIT_SKIP_REVERSE).
        # Trades some false positives (paralogs that happen to score very
        # well) for ~50% time saving per species.
        if fwd_evalue <= STRONG_HIT_SKIP_REVERSE:
            return RBHResult(query=query_accession,
                             target_species=target_taxon_id,
                             tool="mmseqs2",
                             forward_hit=forward_hit,
                             reverse_hit=None,
                             reciprocal=True,
                             forward_evalue=fwd_evalue,
                             forward_score=fwd_score,
                             skipped_reverse=True)

        # Reverse: extract forward_hit's sequence from target_proteome_fasta,
        # then search it against home_proteome_fasta.
        rev_query_fasta = os.path.join(tmpdir, "rev_query.fasta")
        if not _extract_sequence(target_proteome_fasta, forward_hit,
                                 rev_query_fasta):
            # We found a hit accession but can't extract its sequence —
            # this can happen if the accession format on the target ID
            # mismatched what was in the FASTA. Mark as non-reciprocal.
            return RBHResult(query=query_accession,
                             target_species=target_taxon_id,
                             tool="mmseqs2",
                             forward_hit=forward_hit,
                             reverse_hit=None,
                             reciprocal=False,
                             forward_evalue=fwd_evalue,
                             forward_score=fwd_score)

        reverse_tsv = os.path.join(tmpdir, "rev.tsv")
        if home_db_path:
            _run_mmseqs_vs_indexed_db(
                rev_query_fasta, home_db_path, reverse_tsv,
                os.path.join(tmpdir, "rev_tmp"),
                evalue=evalue, threads=threads)
        else:
            _run_mmseqs_one_vs_db(
                rev_query_fasta, home_proteome_fasta, reverse_tsv,
                os.path.join(tmpdir, "rev_tmp"),
                evalue=evalue, threads=threads)
        with open(reverse_tsv) as fh:
            rev_text = fh.read()
        rev_best = best_hit_from_blasttab(rev_text)
        if rev_best is None:
            return RBHResult(query=query_accession,
                             target_species=target_taxon_id,
                             tool="mmseqs2",
                             forward_hit=forward_hit, reverse_hit=None,
                             reciprocal=False,
                             forward_evalue=fwd_evalue,
                             forward_score=fwd_score)

        reverse_hit, _rev_e, _rev_b = rev_best
        reciprocal = _accession_match(reverse_hit, query_accession)
        return RBHResult(query=query_accession,
                         target_species=target_taxon_id,
                         tool="mmseqs2",
                         forward_hit=forward_hit,
                         reverse_hit=reverse_hit,
                         reciprocal=reciprocal,
                         forward_evalue=fwd_evalue,
                         forward_score=fwd_score)
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)


def rbh_check_hmmer(query_accession: str,
                    query_fasta: str,
                    home_proteome_fasta: str,
                    target_proteome_fasta: str,
                    target_taxon_id: str,
                    threads: int = 4,
                    evalue: float = 1e-3,
                    ) -> RBHResult:
    """Run phmmer forward + reverse search for one (Q, S) pair.

    Same algorithm as rbh_check_mmseqs2 but with HMMER's phmmer. More
    sensitive at low identity, slower.
    """
    tmpdir = tempfile.mkdtemp(prefix="rbh_hmmer_")
    try:
        forward_tbl = os.path.join(tmpdir, "fwd.tblout")
        _run_phmmer_one_vs_db(query_fasta, target_proteome_fasta,
                              forward_tbl, evalue=evalue, threads=threads)
        with open(forward_tbl) as fh:
            fwd_text = fh.read()
        fwd_best = best_hit_from_phmmer_tblout(fwd_text)
        if fwd_best is None:
            return RBHResult(query=query_accession,
                             target_species=target_taxon_id,
                             tool="hmmer",
                             forward_hit=None, reverse_hit=None,
                             reciprocal=False)
        forward_hit, fwd_evalue, fwd_score = fwd_best

        # Strong-hit short-circuit (see STRONG_HIT_SKIP_REVERSE)
        if fwd_evalue <= STRONG_HIT_SKIP_REVERSE:
            return RBHResult(query=query_accession,
                             target_species=target_taxon_id,
                             tool="hmmer",
                             forward_hit=forward_hit,
                             reverse_hit=None,
                             reciprocal=True,
                             forward_evalue=fwd_evalue,
                             forward_score=fwd_score,
                             skipped_reverse=True)

        rev_query = os.path.join(tmpdir, "rev_query.fasta")
        if not _extract_sequence(target_proteome_fasta, forward_hit,
                                 rev_query):
            return RBHResult(query=query_accession,
                             target_species=target_taxon_id,
                             tool="hmmer",
                             forward_hit=forward_hit, reverse_hit=None,
                             reciprocal=False,
                             forward_evalue=fwd_evalue,
                             forward_score=fwd_score)

        reverse_tbl = os.path.join(tmpdir, "rev.tblout")
        _run_phmmer_one_vs_db(rev_query, home_proteome_fasta,
                              reverse_tbl, evalue=evalue, threads=threads)
        with open(reverse_tbl) as fh:
            rev_text = fh.read()
        rev_best = best_hit_from_phmmer_tblout(rev_text)
        if rev_best is None:
            return RBHResult(query=query_accession,
                             target_species=target_taxon_id,
                             tool="hmmer",
                             forward_hit=forward_hit, reverse_hit=None,
                             reciprocal=False,
                             forward_evalue=fwd_evalue,
                             forward_score=fwd_score)
        reverse_hit, _re, _rb = rev_best
        reciprocal = _accession_match(reverse_hit, query_accession)
        return RBHResult(query=query_accession,
                         target_species=target_taxon_id,
                         tool="hmmer",
                         forward_hit=forward_hit,
                         reverse_hit=reverse_hit,
                         reciprocal=reciprocal,
                         forward_evalue=fwd_evalue,
                         forward_score=fwd_score)
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)


def rbh_check_foldseek(query_accession: str,
                       query_pdb: str,
                       home_taxon_id: str,
                       target_taxon_id: str,
                       foldseek_db_path: str,
                       cache_dir: str,
                       threads: int = 4,
                       evalue: float = 1e-3,
                       ) -> RBHResult:
    """Run foldseek forward + reverse search for one (Q, S) pair using
    taxon-restricted search against the full AFDB foldseek DB.

    Args:
        query_pdb: PDB file of the query protein.
        home_taxon_id, target_taxon_id: NCBI taxon IDs used as
            --taxon-list for forward (target) and reverse (home).
        foldseek_db_path: path to the full AFDB foldseek DB.
        cache_dir: where to cache fetched PDBs.
    """
    # We need to fetch the forward_hit's AlphaFold PDB for the reverse
    # search. Import here to avoid a circular dep with fetch.py.
    import sys as _sys
    _here = os.path.dirname(os.path.abspath(__file__))
    if _here not in _sys.path:
        _sys.path.insert(0, _here)
    from fetch import fetch_structure

    tmpdir = tempfile.mkdtemp(prefix="rbh_foldseek_")
    try:
        forward_tsv = os.path.join(tmpdir, "fwd.m8")
        _run_foldseek_one_vs_db(
            query_pdb, foldseek_db_path, forward_tsv,
            os.path.join(tmpdir, "fwd_tmp"),
            evalue=evalue, threads=threads,
            taxon_list=str(target_taxon_id))
        with open(forward_tsv) as fh:
            fwd_text = fh.read()
        fwd_best = best_hit_from_blasttab(fwd_text)
        if fwd_best is None:
            return RBHResult(query=query_accession,
                             target_species=target_taxon_id,
                             tool="foldseek",
                             forward_hit=None, reverse_hit=None,
                             reciprocal=False)
        forward_hit, fwd_evalue, fwd_score = fwd_best

        # Strong-hit short-circuit (see STRONG_HIT_SKIP_REVERSE)
        if fwd_evalue <= STRONG_HIT_SKIP_REVERSE:
            return RBHResult(query=query_accession,
                             target_species=target_taxon_id,
                             tool="foldseek",
                             forward_hit=forward_hit,
                             reverse_hit=None,
                             reciprocal=True,
                             forward_evalue=fwd_evalue,
                             forward_score=fwd_score,
                             skipped_reverse=True)

        # Fetch forward_hit's PDB for the reverse search.
        rev_pdb = fetch_structure(forward_hit, cache_dir)
        if rev_pdb is None:
            return RBHResult(query=query_accession,
                             target_species=target_taxon_id,
                             tool="foldseek",
                             forward_hit=forward_hit, reverse_hit=None,
                             reciprocal=False,
                             forward_evalue=fwd_evalue,
                             forward_score=fwd_score)

        reverse_tsv = os.path.join(tmpdir, "rev.m8")
        _run_foldseek_one_vs_db(
            rev_pdb, foldseek_db_path, reverse_tsv,
            os.path.join(tmpdir, "rev_tmp"),
            evalue=evalue, threads=threads,
            taxon_list=str(home_taxon_id))
        with open(reverse_tsv) as fh:
            rev_text = fh.read()
        rev_best = best_hit_from_blasttab(rev_text)
        if rev_best is None:
            return RBHResult(query=query_accession,
                             target_species=target_taxon_id,
                             tool="foldseek",
                             forward_hit=forward_hit, reverse_hit=None,
                             reciprocal=False,
                             forward_evalue=fwd_evalue,
                             forward_score=fwd_score)
        reverse_hit, _re, _rb = rev_best
        reciprocal = _accession_match(reverse_hit, query_accession)
        return RBHResult(query=query_accession,
                         target_species=target_taxon_id,
                         tool="foldseek",
                         forward_hit=forward_hit,
                         reverse_hit=reverse_hit,
                         reciprocal=reciprocal,
                         forward_evalue=fwd_evalue,
                         forward_score=fwd_score)
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)


def _extract_sequence(fasta_path: str, accession: str,
                      out_path: str) -> bool:
    """Pull a single FASTA record matching `accession` from a multi-FASTA file
    and write it to out_path. Returns True if found.

    Matches against the normalised accession in the header (handles
    sp|ACC|name, tr|ACC|name formats)."""
    target_acc = _normalise_accession(accession)
    capture = False
    captured_lines: list[str] = []
    with open(fasta_path) as fh:
        for line in fh:
            if line.startswith(">"):
                if capture:
                    break
                # Parse header for the accession.
                hdr = line[1:].split()[0]
                if _accession_match(hdr, target_acc):
                    capture = True
                    captured_lines.append(line)
            elif capture:
                captured_lines.append(line)
    if not captured_lines:
        return False
    with open(out_path, "w") as fh:
        fh.writelines(captured_lines)
    return True
