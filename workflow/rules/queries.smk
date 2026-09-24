# Fetch the query sequence (FASTA) and, when foldseek is enabled, the
# AlphaFold structure (PDB) for each primary + cooccurrence UniProt id.

rule fetch_query_fasta:
    output:
        f"{OUTDIR}/queries/{{group}}.fasta"
    params:
        cfg=CFG_PATH,
        cache=CFG["cache_dir"],
    conda:
        "../../envs/fetch.yaml"
    shell:
        r"""
        python workflow/scripts/fetch_query.py \
            --config {params.cfg} --accession {wildcards.group} \
            --kind fasta --out {output}
        """

rule fetch_query_structure:
    output:
        f"{OUTDIR}/queries/{{group}}.pdb"
    params:
        cfg=CFG_PATH,
    conda:
        "../../envs/fetch.yaml"
    shell:
        r"""
        python workflow/scripts/fetch_query.py \
            --config {params.cfg} --accession {wildcards.group} \
            --kind structure --out {output}
        """
