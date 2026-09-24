# Combine the per-tool hit TSVs for a group into consensus + species-deduped
# members, and write the group FASTA (headers = species keys).

def _group_hit_tsvs(wildcards):
    return [f"{OUTDIR}/search/{wildcards.group}.{t}.tsv" for t in TOOLS]

rule build_group:
    input:
        hits=_group_hit_tsvs
    output:
        members=f"{OUTDIR}/groups/{{group}}.members.tsv",
        fasta=f"{OUTDIR}/groups/{{group}}.fasta",
    params:
        cfg=CFG_PATH,
    conda:
        "../../envs/fetch.yaml"
    shell:
        r"""
        python workflow/scripts/build_group.py \
            --config {params.cfg} --group {wildcards.group} \
            --hits {input.hits} \
            --out-members {output.members} --out-fasta {output.fasta}
        """
