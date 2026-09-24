# Build iTOL annotation files from the tree and cooccurrence CSV.
# Runs after both have been produced. Pure Python, no external tools, fast.

rule itol:
    input:
        tree=f"{OUTDIR}/tree/species_tree.nwk",
        cooccurrence=f"{OUTDIR}/cooccurrence.csv",
        collapsed=f"{OUTDIR}/cooccurrence_collapsed.csv",
    output:
        marker=f"{OUTDIR}/itol/.done",
    params:
        cooccur_ids=" ".join(COOCCUR) if COOCCUR else "",
        primary_groups=" ".join(PRIMARY) if PRIMARY else "",
        out_dir=f"{OUTDIR}/itol",
    conda:
        "../../envs/fetch.yaml"
    shell:
        r"""
        python workflow/scripts/build_itol.py \
            --tree {input.tree} \
            --cooccurrence-csv {input.cooccurrence} \
            --collapsed-csv {input.collapsed} \
            --out-dir {params.out_dir} \
            --cooccur-ids {params.cooccur_ids} \
            --primary-groups {params.primary_groups}
        touch {output.marker}
        """
