# Align each PRIMARY ortholog group separately (align-then-concatenate).
# Co-occurrence groups are not aligned; they only feed the CSV.
#
# When primary_rbh_mode = 'filter', the input fasta is the RBH-filtered
# version (fewer species) instead of the raw build_group output.

def _align_input_fasta(wildcards):
    # Cooccur groups (if they ever get aligned in future) always use
    # the raw fasta. Primaries use the helper which respects RBH mode.
    if wildcards.group in PRIMARY:
        return primary_fasta_for_alignment(wildcards.group)
    return f"{OUTDIR}/groups/{wildcards.group}.fasta"


rule align_group:
    input:
        fasta=_align_input_fasta
    output:
        aln=f"{OUTDIR}/align/{{group}}.aln.fasta"
    params:
        cfg=CFG_PATH,
        work=lambda w: f"{TMPDIR}/align_{w.group}",
        aligner=CFG["aligner"],
    threads: THREADS
    conda:
        "../../envs/align.yaml"
    shell:
        r"""
        python workflow/scripts/align_cli.py \
            --config {params.cfg} --in {input.fasta} --out {output.aln} \
            --work {params.work} --threads {threads}
        """
