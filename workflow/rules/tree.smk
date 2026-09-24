# Infer the final species tree from the trimmed supermatrix; output newick.

rule tree:
    input:
        fasta=f"{OUTDIR}/supermatrix/supermatrix.fasta"
    output:
        nwk=f"{OUTDIR}/tree/species_tree.nwk"
    params:
        cfg=CFG_PATH,
    conda:
        "../../envs/align.yaml"
    shell:
        r"""
        python workflow/scripts/tree_cli.py \
            --config {params.cfg} --in {input.fasta} --out {output.nwk}
        """
