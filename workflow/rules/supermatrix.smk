# Concatenate the trimmed PRIMARY group alignments into the supermatrix
# (one row per species, gap-padded for missing groups) + partition file.

rule supermatrix:
    input:
        alns=expand(f"{OUTDIR}/trim/{{group}}.trim.fasta", group=PRIMARY)
    output:
        fasta=f"{OUTDIR}/supermatrix/supermatrix.fasta",
        parts=f"{OUTDIR}/supermatrix/partitions.txt",
    params:
        cfg=CFG_PATH,
    conda:
        "../../envs/fetch.yaml"
    shell:
        r"""
        python workflow/scripts/build_supermatrix.py \
            --config {params.cfg} --alignments {input.alns} \
            --out-fasta {output.fasta} --out-partitions {output.parts}
        """
