# Trim each per-group alignment with ClipKIT.

rule trim_group:
    input:
        aln=f"{OUTDIR}/align/{{group}}.aln.fasta"
    output:
        trimmed=f"{OUTDIR}/trim/{{group}}.trim.fasta"
    params:
        cfg=CFG_PATH,
    conda:
        "../../envs/align.yaml"
    shell:
        r"""
        python workflow/scripts/trim_cli.py \
            --config {params.cfg} --in {input.aln} --out {output.trimmed}
        """
