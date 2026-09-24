#!/usr/bin/env python3
"""
Parse and validate the CoPhy config.txt.

The config is INI-style (sections in [brackets], key = value lines).
We deliberately use a hand-rolled reader rather than configparser so the
file can stay a plain `.txt` per the spec and so we control type coercion,
boolean parsing, and the inline-list / file-path conventions in one place.

Public entry point: load_config(path) -> dict with typed, validated values.
"""

from __future__ import annotations

import os
import re
import sys
from typing import Any


# --------------------------------------------------------------------------
# Low-level reading
# --------------------------------------------------------------------------

_SECTION_RE = re.compile(r"^\[(?P<name>[A-Za-z0-9_]+)\]\s*$")


def _read_raw(path: str) -> dict[str, dict[str, str]]:
    """Read the file into {section: {key: raw_string_value}}."""
    if not os.path.exists(path):
        raise FileNotFoundError(f"Config file not found: {path}")

    sections: dict[str, dict[str, str]] = {}
    current: str | None = None

    with open(path, "r", encoding="utf-8") as fh:
        for lineno, raw in enumerate(fh, start=1):
            line = raw.rstrip("\n")
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue

            m = _SECTION_RE.match(stripped)
            if m:
                current = m.group("name").lower()
                sections.setdefault(current, {})
                continue

            if "=" not in stripped:
                raise ValueError(
                    f"{path}:{lineno}: expected 'key = value' or '[section]', "
                    f"got: {stripped!r}"
                )
            if current is None:
                raise ValueError(
                    f"{path}:{lineno}: key/value before any [section] header"
                )

            key, _, value = stripped.partition("=")
            key = key.strip().lower()
            # Strip inline comments: everything from the first unquoted '#'.
            # No config value in this pipeline legitimately contains '#'.
            value = value.split("#", 1)[0].strip()
            if not key:
                raise ValueError(f"{path}:{lineno}: empty key")
            sections[current][key] = value

    return sections


# --------------------------------------------------------------------------
# Type coercion helpers
# --------------------------------------------------------------------------

_TRUE = {"true", "yes", "on", "1"}
_FALSE = {"false", "no", "off", "0"}


def _as_bool(value: str, *, ctx: str) -> bool:
    v = value.strip().lower()
    if v in _TRUE:
        return True
    if v in _FALSE:
        return False
    raise ValueError(f"{ctx}: expected a boolean (true/false), got {value!r}")


def _as_int(value: str, *, ctx: str) -> int:
    try:
        return int(value)
    except ValueError:
        raise ValueError(f"{ctx}: expected an integer, got {value!r}")


def _as_float(value: str, *, ctx: str) -> float:
    try:
        return float(value)
    except ValueError:
        raise ValueError(f"{ctx}: expected a number, got {value!r}")


def _as_id_list(inline: str, file_path: str, *, ctx: str, base_dir: str) -> list[str]:
    """
    Resolve an "inline list OR file path" pair into a clean list of IDs.
    If a file path is given and exists, it wins; one ID per line.
    Inline values may be comma- and/or whitespace-separated.
    """
    if file_path:
        resolved = file_path
        if not os.path.isabs(resolved):
            resolved = os.path.join(base_dir, resolved)
        if not os.path.exists(resolved):
            raise FileNotFoundError(f"{ctx}: id file not found: {resolved}")
        ids: list[str] = []
        with open(resolved, "r", encoding="utf-8") as fh:
            for line in fh:
                tok = line.strip()
                if tok and not tok.startswith("#"):
                    ids.append(tok)
    else:
        ids = [t for t in re.split(r"[,\s]+", inline.strip()) if t]

    # De-duplicate, preserve order.
    seen: set[str] = set()
    out: list[str] = []
    for i in ids:
        if i not in seen:
            seen.add(i)
            out.append(i)
    return out


# --------------------------------------------------------------------------
# Public API
# --------------------------------------------------------------------------

def load_config(path: str) -> dict[str, Any]:
    """Read, type, and validate the config. Returns a nested dict."""
    base_dir = os.path.dirname(os.path.abspath(path))
    raw = _read_raw(path)
    
    def get(section: str, key: str, default: str | None = None) -> str:
        sec = raw.get(section, {})
        if key in sec:
            return sec[key]
        if default is not None:
            return default
        raise KeyError(f"Missing required config key [{section}] {key}")

    cfg: dict[str, Any] = {}

    # ---- input ----
    cfg["query_ids"] = _as_id_list(
        get("input", "query_ids", ""),
        get("input", "query_ids_file", ""),
        ctx="[input] query_ids", base_dir=base_dir,
    )
    if not cfg["query_ids"]:
        raise ValueError("[input]: no query UniProt IDs provided")
    cfg["cooccurrence_ids"] = _as_id_list(
        get("input", "cooccurrence_ids", ""),
        get("input", "cooccurrence_ids_file", ""),
        ctx="[input] cooccurrence_ids", base_dir=base_dir,
    )
    cfg["kingdom_filter"] = get("search", "kingdom_filter", "").strip()
    # ---- databases ----
    cfg["uniref_version"] = _as_int(
        get("databases", "uniref_version", "90"),
        ctx="[databases] uniref_version")
    if cfg["uniref_version"] not in (50, 90):
        raise ValueError("[databases] uniref_version must be 50 or 90")
    cfg["uniref_path"] = get("databases", "uniref_path", "")
    # mmseqs2 and HMMER want different forms of the DB; each may override.
    # Fall back to uniref_path if a tool-specific key is empty.
    cfg["mmseqs2_db"] = get("databases", "mmseqs2_db", "") or cfg["uniref_path"]
    cfg["hmmer_db"] = get("databases", "hmmer_db", "") or cfg["uniref_path"]
    # Note: tool-enabled-but-no-DB validation happens after [search] is parsed.
    cfg["foldseek_db_path"] = get("databases", "foldseek_db_path", "")
    cfg["cache_dir"] = get("databases", "cache_dir", "resources/cache")

    # ---- search ----
    cfg["use_mmseqs2"] = _as_bool(get("search", "use_mmseqs2", "true"),
                                  ctx="[search] use_mmseqs2")
    cfg["use_hmmer"] = _as_bool(get("search", "use_hmmer", "true"),
                                ctx="[search] use_hmmer")
    cfg["use_foldseek"] = _as_bool(get("search", "use_foldseek", "false"),
                                   ctx="[search] use_foldseek")

    # Now validate that each enabled tool has a DB.
    if cfg["use_mmseqs2"] and not cfg["mmseqs2_db"]:
        raise ValueError(
            "[databases]: mmseqs2 is enabled but neither mmseqs2_db nor "
            "uniref_path is set")
    if cfg["use_hmmer"] and not cfg["hmmer_db"]:
        raise ValueError(
            "[databases]: HMMER is enabled but neither hmmer_db nor "
            "uniref_path is set")
    enabled = [t for t, on in (("mmseqs2", cfg["use_mmseqs2"]),
                               ("hmmer", cfg["use_hmmer"]),
                               ("foldseek", cfg["use_foldseek"])) if on]
    if not enabled:
        raise ValueError("[search]: at least one search tool must be enabled")
    cfg["enabled_tools"] = enabled

    cfg["consensus_mode"] = get("search", "consensus_mode", "consensus").lower()
    if cfg["consensus_mode"] not in ("consensus", "union", "intersection"):
        raise ValueError("[search] consensus_mode must be "
                         "consensus|union|intersection")
    cfg["min_tools"] = _as_int(get("search", "min_tools", "2"),
                               ctx="[search] min_tools")
    cfg["orthology_mode"] = get("search", "orthology_mode", "forward").lower()
    if cfg["orthology_mode"] not in ("forward", "rbh"):
        raise ValueError("[search] orthology_mode must be forward|rbh")
    cfg["max_hits_per_tool"] = _as_int(
        get("search", "max_hits_per_tool", "0"),
        ctx="[search] max_hits_per_tool")
    if cfg["max_hits_per_tool"] < 0:
        raise ValueError("[search] max_hits_per_tool must be >= 0")
    if cfg["consensus_mode"] == "consensus":
        if cfg["min_tools"] < 1 or cfg["min_tools"] > len(enabled):
            raise ValueError(
                f"[search] min_tools={cfg['min_tools']} is out of range for "
                f"{len(enabled)} enabled tool(s)")

    if cfg["use_foldseek"] and not cfg["foldseek_db_path"]:
        raise ValueError("[search] use_foldseek=true but "
                         "[databases] foldseek_db_path is empty")

    cfg["mmseqs2"] = {
        "sensitivity": _as_float(get("search", "mmseqs2_sensitivity", "7.5"),
                                 ctx="mmseqs2_sensitivity"),
        "evalue": get("search", "mmseqs2_evalue", "1e-5"),
        "min_seq_id": _as_float(get("search", "mmseqs2_min_seq_id", "0.3"),
                                ctx="mmseqs2_min_seq_id"),
        "coverage": _as_float(get("search", "mmseqs2_coverage", "0.5"),
                              ctx="mmseqs2_coverage"),
        "cov_mode": _as_int(get("search", "mmseqs2_cov_mode", "0"),
                            ctx="mmseqs2_cov_mode"),
        "split_memory": get("search", "mmseqs2_split_memory", "0"),
    }
    cfg["hmmer"] = {
        "evalue": get("search", "hmmer_evalue", "1e-5"),
        "dom_evalue": get("search", "hmmer_dom_evalue", "1e-5"),
        "incE": get("search", "hmmer_incE", "1e-5"),
        "bitscore": _as_float(get("search", "hmmer_bitscore", "0"),
                              ctx="hmmer_bitscore"),
    }
    cfg["foldseek"] = {
        "evalue": get("search", "foldseek_evalue", "1e-3"),
        "min_tmscore": _as_float(get("search", "foldseek_min_tmscore", "0.5"),
                                 ctx="foldseek_min_tmscore"),
        "alntmscore": _as_bool(get("search", "foldseek_alntmscore", "true"),
                               ctx="foldseek_alntmscore"),
        "coverage": _as_float(get("search", "foldseek_coverage", "0.5"),
                              ctx="foldseek_coverage"),
        # Memory cap that controls how many splits foldseek does; e.g. "100G".
        # On a 120 GB batch node a 100G cap gives ~20 GB working headroom.
        # Empty means use foldseek's default (which assumes ample RAM).
        "split_memory": get("search", "foldseek_split_memory", "100G").strip(),
        # Prefilter candidates kept per query. Default 1000; bump to 5000 or
        # higher to recover distant homologs missed by the prefilter's
        # cap. Empty means use foldseek's default.
        "max_seqs": get("search", "foldseek_max_seqs", "").strip(),
        # If true, foldseek skips the prefilter and aligns the query
        # against every structure in the DB. Much slower but maximum recall.
        "exhaustive_search": _as_bool(
            get("search", "foldseek_exhaustive_search", "false"),
            ctx="foldseek_exhaustive_search"),
        # Extra raw foldseek flags, appended verbatim. E.g.
        # "--prefilter-mode 1" for the ungapped prefilter, which is not
        # memory-limited and is recommended for single-query searches.
        "extra": get("search", "foldseek_extra", "").strip(),
    }

    # ---- dedup ----
    cfg["dedup_species"] = _as_bool(get("dedup", "dedup_species", "true"),
                                    ctx="[dedup] dedup_species")
    cfg["species_rank"] = get("dedup", "species_rank", "binomial").lower()
    if cfg["species_rank"] not in ("binomial", "subspecies", "none"):
        raise ValueError("[dedup] species_rank must be "
                         "binomial|subspecies|none")
    cfg["tie_break"] = get("dedup", "tie_break",
                           "swissprot_then_lowest_taxid").lower()
    valid_tb = {"swissprot_then_lowest_taxid", "best_score", "longest",
                "lowest_taxid", "swissprot_then_best_score"}
    if cfg["tie_break"] not in valid_tb:
        raise ValueError(f"[dedup] tie_break must be one of {sorted(valid_tb)}")

    # ---- sequences ----
    cfg["seq_type"] = get("sequences", "seq_type", "aa").lower()
    if cfg["seq_type"] not in ("aa", "nt"):
        raise ValueError("[sequences] seq_type must be aa|nt")

    # ---- align ----
    cfg["aligner"] = get("align", "aligner", "mafft").lower()
    if cfg["aligner"] not in ("mafft", "pagan2"):
        raise ValueError("[align] aligner must be mafft|pagan2")
    cfg["mafft_algo"] = get("align", "mafft_algo", "auto").lower()
    if cfg["mafft_algo"] not in ("auto", "linsi", "ginsi", "einsi", "fftns"):
        raise ValueError("[align] mafft_algo invalid")
    cfg["mafft_maxiterate"] = _as_int(get("align", "mafft_maxiterate", "1000"),
                                      ctx="mafft_maxiterate")
    cfg["mafft_extra"] = get("align", "mafft_extra", "")
    cfg["pagan2_extra"] = get("align", "pagan2_extra", "")

    # ---- trim ----
    cfg["clipkit_mode"] = get("trim", "clipkit_mode", "smart-gap")
    cfg["clipkit_gaps"] = _as_float(get("trim", "clipkit_gaps", "0.9"),
                                    ctx="clipkit_gaps")
    cfg["clipkit_extra"] = get("trim", "clipkit_extra", "")

    # ---- tree ----
    cfg["fasttree_model"] = get("tree", "fasttree_model", "lg").lower()
    cfg["fasttree_gamma"] = _as_bool(get("tree", "fasttree_gamma", "true"),
                                     ctx="fasttree_gamma")
    cfg["fasttree_extra"] = get("tree", "fasttree_extra", "")
    cfg["fasttree_bin"] = get("tree", "fasttree_bin", "FastTree")
    _aa_models = {"lg", "wag", "jtt"}
    _nt_models = {"gtr", "jc"}
    if cfg["seq_type"] == "aa" and cfg["fasttree_model"] not in _aa_models:
        raise ValueError(f"[tree] fasttree_model for aa must be in {_aa_models}")
    if cfg["seq_type"] == "nt" and cfg["fasttree_model"] not in _nt_models:
        raise ValueError(f"[tree] fasttree_model for nt must be in {_nt_models}")

    # ---- resources ----
    cfg["threads"] = _as_int(get("resources", "threads", "8"),
                             ctx="[resources] threads")

    # ---- neighborhood ----
    cfg["neighborhood"] = _as_bool(
        get("neighborhood", "neighborhood", "false"),
        ctx="[neighborhood] neighborhood")
    cfg["neighborhood_max_dist"] = _as_int(
        get("neighborhood", "neighborhood_max_dist", "20"),
        ctx="[neighborhood] neighborhood_max_dist")
    if cfg["neighborhood_max_dist"] < 1:
        raise ValueError("[neighborhood] neighborhood_max_dist must be >= 1")
    cfg["ncbi_api_key"] = get("neighborhood", "ncbi_api_key", "").strip()
    cfg["ncbi_email"] = get("neighborhood", "ncbi_email", "").strip()
    cfg["taxonomy"] = _as_bool(
        get("neighborhood", "taxonomy", "true"),
        ctx="[neighborhood] taxonomy")
    # Cooccurrence detection mode: 'search' (scan all of UniRef for each
    # cooccur protein, then intersect with tree species) or 'rbh'
    # (species-restricted reciprocal best hit, faster and more rigorous
    # for non-trivial proteomes).
    cfg["cooccurrence_mode"] = get(
        "neighborhood", "cooccurrence_mode", "search").strip().lower()
    if cfg["cooccurrence_mode"] not in ("search", "rbh"):
        raise ValueError(
            f"[neighborhood] cooccurrence_mode must be 'search' or 'rbh', "
            f"got {cfg['cooccurrence_mode']!r}")
    # Primary RBH validation mode: 'off' (default, no change), 'annotate'
    # (mark species with RBH pass/fail but build tree from full set), or
    # 'filter' (drop species whose RBH failed across ALL their primary
    # memberships - "ANY" rule applied across primaries).
    cfg["primary_rbh_mode"] = get(
        "search", "primary_rbh_mode", "off").strip().lower()
    if cfg["primary_rbh_mode"] not in ("off", "annotate", "filter"):
        raise ValueError(
            f"[search] primary_rbh_mode must be 'off', 'annotate', or "
            f"'filter', got {cfg['primary_rbh_mode']!r}")
    cfg["neighborhood_max_lookups_per_run"] = _as_int(
        get("neighborhood", "neighborhood_max_lookups_per_run", "0"),
        ctx="[neighborhood] neighborhood_max_lookups_per_run")
    if cfg["neighborhood_max_lookups_per_run"] < 0:
        raise ValueError("max_lookups_per_run must be >= 0")
    if cfg["neighborhood"] and not cfg["ncbi_email"]:
        raise ValueError(
            "[neighborhood] ncbi_email is required when neighborhood=true "
            "(NCBI's terms of service require an email for programmatic use)")

    return cfg


if __name__ == "__main__":
    # CLI: validate a config and pretty-print the parsed result.
    import json
    if len(sys.argv) != 2:
        print("usage: config.py <config.txt>", file=sys.stderr)
        sys.exit(2)
    try:
        parsed = load_config(sys.argv[1])
    except (ValueError, KeyError, FileNotFoundError) as exc:
        print(f"CONFIG ERROR: {exc}", file=sys.stderr)
        sys.exit(1)
    print(json.dumps(parsed, indent=2, default=str))
    print("\nOK: config is valid.", file=sys.stderr)
