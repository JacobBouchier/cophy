# Build the co-occurrence CSV (step 9).
#
# Two input modes, switched by config key `cooccurrence_mode`:
#   * 'search' (default, legacy): consume per-cooccur members.tsv from the
#     same search → build_group pipeline as primaries.
#   * 'rbh': consume per-cooccur RBH TSVs from the species-restricted
#     reciprocal-best-hit pipeline.

_COOCCUR_MODE = CFG.get("cooccurrence_mode", "search").strip().lower()
if _COOCCUR_MODE not in ("search", "rbh"):
    raise ValueError(
        f"[cooccurrence] cooccurrence_mode must be 'search' or 'rbh', "
        f"got {_COOCCUR_MODE!r}")


def _cooccurrence_inputs(wildcards):
    # Pick the right primary members file based on whether RBH is on.
    # primary_members_for_alignment() comes from primary_rbh.smk.
    primaries = [primary_members_for_alignment(g) for g in PRIMARY]
    if _COOCCUR_MODE == "rbh":
        return {
            "primary": primaries,
            "cooccur_rbh": expand(
                f"{OUTDIR}/cooccur_rbh/{{cooccur}}.tsv", cooccur=COOCCUR),
        }
    return {
        "primary": primaries,
        "cooccur": expand(
            f"{OUTDIR}/groups/{{group}}.members.tsv", group=COOCCUR),
    }


rule cooccurrence:
    input:
        unpack(_cooccurrence_inputs),
    output:
        csv=f"{OUTDIR}/cooccurrence.csv",
        collapsed=f"{OUTDIR}/cooccurrence_collapsed.csv",
    params:
        cfg=CFG_PATH,
        cooccur_args=lambda w, input: (
            ("--cooccur-rbh " + " ".join(input.cooccur_rbh))
            if _COOCCUR_MODE == "rbh"
            else (("--cooccur-members " + " ".join(input.cooccur))
                  if COOCCUR else "")
        ),
    conda:
        "../../envs/fetch.yaml"
    shell:
        r"""
        python workflow/scripts/build_cooccurrence.py \
            --config {params.cfg} \
            --primary-members {input.primary} \
            {params.cooccur_args} \
            --out-csv {output.csv} \
            --out-collapsed-csv {output.collapsed}
        """
