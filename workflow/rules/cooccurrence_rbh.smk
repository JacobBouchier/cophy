# RBH-based cooccurrence detection: per-cooccur-id Slurm job that runs
# all per-species RBH checks for one cooccur query.
#
# This replaces the all-vs-all cooccur search (search.smk on cooccur ids)
# when `cooccurrence_mode = rbh` in the config.

# Rule 1: build the species taxon list (union across primary members).
# This is the input for every per-cooccur RBH job.
#
# Uses primary_members_for_alignment() so that, if primary RBH is on, the
# tree species list reflects the post-RBH species set (the species that
# will actually appear in the final tree).
rule build_species_taxon_list:
    input:
        members=lambda w: [primary_members_for_alignment(g) for g in PRIMARY]
    output:
        f"{OUTDIR}/cooccur_rbh/species_taxon_list.txt"
    conda:
        "../../envs/fetch.yaml"
    shell:
        r"""
        python workflow/scripts/build_species_taxon_list.py \
            --primary-members {input.members} \
            --out {output}
        """


# Rule 2: per-cooccur-id RBH check across all tree species.
# Snakemake parallelises across cooccur ids (7 jobs total in the
# typical case), each job loops through ~1300 species internally.
rule cooccurrence_rbh:
    input:
        taxon_list=f"{OUTDIR}/cooccur_rbh/species_taxon_list.txt",
    output:
        tsv=f"{OUTDIR}/cooccur_rbh/{{cooccur}}.tsv",
    params:
        cfg=CFG_PATH,
        tools="mmseqs2,hmmer",   # foldseek RBH dropped: --taxon-list
                                  # doesn't reduce DB-load memory, making
                                  # per-species foldseek too expensive.
    threads: 8
    conda:
        "../../envs/rbh.yaml"
    shell:
        r"""
        python workflow/scripts/build_cooccurrence_rbh.py \
            --config {params.cfg} \
            --cooccur-id {wildcards.cooccur} \
            --species-list {input.taxon_list} \
            --tools {params.tools} \
            --out-tsv {output.tsv} \
            --threads {threads}
        """
