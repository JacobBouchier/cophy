# Primary RBH validation — optional post-filter on the species set
# defined by build_group (i.e. on the "made the cutoff" species).
#
# Mode (config: primary_rbh_mode):
#   off       (default) - the existing members.tsv flows downstream unchanged
#   annotate  - annotate each member with rbh_passed_anywhere yes/no, no filter
#   filter    - drop species whose RBH failed for every primary they appear in
#
# Cross-primary union rule: a species passes iff RBH succeeded in ANY of its
# primary memberships (parallels consensus_mode = union).
#
# These rules are emitted only when primary_rbh_mode is enabled - otherwise
# downstream rules see members.tsv directly (no DAG change).

_PRIMARY_RBH_MODE = CFG.get("primary_rbh_mode", "off").strip().lower()
if _PRIMARY_RBH_MODE not in ("off", "annotate", "filter"):
    raise ValueError(
        f"[search] primary_rbh_mode must be 'off', 'annotate', or 'filter', "
        f"got {_PRIMARY_RBH_MODE!r}")


# Only register the rules if primary RBH is enabled. When 'off',
# downstream rules continue to use <group>.members.tsv directly.
if _PRIMARY_RBH_MODE != "off":

    # Per-primary RBH: take members.tsv, RBH-validate each species,
    # emit members_rbh.tsv (annotated with rbh_reciprocal_tools column).
    rule primary_rbh:
        input:
            members=f"{OUTDIR}/groups/{{primary}}.members.tsv",
        output:
            tsv=f"{OUTDIR}/groups/{{primary}}.members_rbh.tsv",
        params:
            cfg=CFG_PATH,
            tools="mmseqs2,hmmer",  # foldseek RBH was too memory-expensive
                                    # at scale; species found only via
                                    # foldseek at the original search stage
                                    # are passed through by the filter step.
        threads: 8
        conda:
            "../../envs/rbh.yaml"
        wildcard_constraints:
            primary="|".join(re.escape(p) for p in PRIMARY) if PRIMARY else "x",
        shell:
            r"""
            python workflow/scripts/build_primary_rbh.py \
                --config {params.cfg} \
                --primary-id {wildcards.primary} \
                --members {input.members} \
                --tools {params.tools} \
                --out-tsv {output.tsv} \
                --threads {threads}
            """

    # Cross-primary union filter / annotation. Input: all primaries' RBH
    # files. Output: per-primary members_filtered.tsv that downstream rules
    # consume in place of members.tsv, and (in filter mode) a matching
    # filtered FASTA file for the alignment step.
    rule primary_rbh_filter:
        input:
            all_rbh=expand(
                f"{OUTDIR}/groups/{{group}}.members_rbh.tsv", group=PRIMARY),
            this_rbh=f"{OUTDIR}/groups/{{primary}}.members_rbh.tsv",
            orig_fasta=f"{OUTDIR}/groups/{{primary}}.fasta",
        output:
            members=f"{OUTDIR}/groups/{{primary}}.members_filtered.tsv",
            fasta=f"{OUTDIR}/groups/{{primary}}.fasta_filtered.fasta",
        params:
            mode=_PRIMARY_RBH_MODE,
            kingdom_filter=CFG.get("kingdom_filter", "").strip(),
            cache_dir=CFG["cache_dir"],
        conda:
            "../../envs/fetch.yaml"
        wildcard_constraints:
            primary="|".join(re.escape(p) for p in PRIMARY) if PRIMARY else "x",
        shell:
            r"""
            python workflow/scripts/build_primary_rbh_filter.py \
                --primary-rbh {input.all_rbh} \
                --members-in {input.this_rbh} \
                --mode {params.mode} \
                --kingdom-filter "{params.kingdom_filter}" \
                --cache-dir {params.cache_dir} \
                --out {output.members} \
                --fasta-in {input.orig_fasta} \
                --fasta-out {output.fasta}
            """


# Helper function for other rules: given a primary, return the members
# file and FASTA downstream rules should consume. When RBH is off, that's
# the raw members.tsv + fasta. When on, it's the filtered versions.
def primary_members_for_alignment(primary: str) -> str:
    if _PRIMARY_RBH_MODE == "off":
        return f"{OUTDIR}/groups/{primary}.members.tsv"
    return f"{OUTDIR}/groups/{primary}.members_filtered.tsv"


def primary_fasta_for_alignment(primary: str) -> str:
    if _PRIMARY_RBH_MODE == "off":
        return f"{OUTDIR}/groups/{primary}.fasta"
    return f"{OUTDIR}/groups/{primary}.fasta_filtered.fasta"
