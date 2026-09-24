# Run each enabled search tool for each group, producing a per-(group,tool)
# normalised hit TSV. The query input differs by tool (FASTA vs structure).

def _search_query_input(wildcards):
    return query_path(wildcards.group, wildcards.tool)

rule search:
    input:
        query=_search_query_input
    output:
        f"{OUTDIR}/search/{{group}}.{{tool}}.tsv"
    params:
        cfg=CFG_PATH,
        tmp=lambda w: f"{TMPDIR}/search_{w.group}_{w.tool}",
    threads: THREADS
    conda:
        lambda w: f"../../envs/{w.tool}.yaml"
    shell:
        r"""
        python workflow/scripts/search_run.py \
            --config {params.cfg} --group {wildcards.group} \
            --tool {wildcards.tool} --query {input.query} \
            --out {output} --tmp {params.tmp} --threads {threads}
        """
