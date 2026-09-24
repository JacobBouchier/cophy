# CoPhy

A Snakemake workflow for building concatenated (supermatrix) phylogenies of
protein families and analyzing the genomic co-occurrence of associated
proteins. Starting from one or more UniProt query accessions, CoPhy
searches for orthologs with sequence and structure tools, aligns and trims
each ortholog group, concatenates them into a supermatrix, infers a species
tree with FastTree2, and produces iTOL annotation files. It can additionally
test whether a second set of proteins co-occurs (and is genomically adjacent)
in the same species.

## Overview

For each **primary** query the pipeline:

1. **Fetches** the query sequence (and, for foldseek, its AlphaFold structure).
2. **Searches** for homologs with any combination of MMseqs2, HMMER (phmmer),
   and foldseek.
3. **Combines** the per-tool hits into a consensus ortholog set and
   **deduplicates** to one representative per species.
4. Optionally **validates** each species by reciprocal best hit (RBH) and
   **filters** by taxonomic domain (e.g. Bacteria only).
5. **Aligns** (MAFFT or PAGAN2) and **trims** (ClipKIT) each ortholog group.
6. **Concatenates** the trimmed alignments into a gap-padded supermatrix with a
   partition file.
7. **Infers** the species tree with FastTree2.
8. **Builds** iTOL annotation files (taxonomy color strips, co-occurrence
   presence, optional gene-neighborhood distance gradients).

Given **co-occurrence** query accessions, it also reports, for each species in
the tree, whether each co-occurrence protein is present (by search or RBH) and,
optionally, how far apart it is from the primary on the chromosome.

```
queries ─▶ search ─▶ build_group ─▶ [primary_rbh + filter] ─▶ align ─▶ trim ─┐
                                                                              ▼
                                                          supermatrix ─▶ tree ─▶ itol
                                    cooccurrence (search or RBH) ──────────────┘
```

## Requirements

- **Snakemake** ≥ 8 (uses the `snakemake-executor-plugin-slurm` plugin for
  cluster execution) and **conda/mamba** for per-rule environments.
- The bioinformatics tools are installed automatically by Snakemake into
  per-rule conda environments (`--use-conda`) from the specs in `envs/`:
  MMseqs2, HMMER, foldseek, MAFFT, ClipKIT, FastTree2 (and optionally PAGAN2).
- The Python scripts use only the standard library (no third-party packages).
- **Sequence/structure databases** (not included; see below):
  - MMseqs2-formatted UniRef90 (or UniRef50) database
  - the same UniRef as a plain FASTA for HMMER
  - a foldseek AlphaFold DB (only if foldseek is enabled)

### Databases

Point the `[databases]` section of your config at local copies:

```bash
# MMseqs2 UniRef90 index (mmseqs2_db)
mmseqs databases UniRef90 /path/to/databases/uniref90_db tmp

# UniRef90 FASTA for HMMER (hmmer_db) — from the UniProt release you used
# foldseek AlphaFold/UniProt DB (foldseek_db_path), only if foldseek is enabled
foldseek databases Alphafold/UniProt /path/to/databases/afdb_foldseek/afdb_uniprot_db tmp
```

## Installation

```bash
git clone https://github.com/<your-username>/cophy.git
cd cophy

# a small controller environment with Snakemake
conda create -n snakemake -c conda-forge -c bioconda snakemake mamba
conda activate snakemake
```

The heavy tools are pulled in automatically per rule via `--use-conda`.

## Configuration

Copy the template and edit it:

```bash
cp config/config.example.txt config/config.txt
```

Key settings (see `config/config.example.txt` for the fully commented list):

| Section | Key | Meaning |
|---|---|---|
| `[input]` | `project_name` | Results go to `results/<project_name>/` |
| `[input]` | `query_ids` | Primary UniProt accession(s) defining the ortholog groups |
| `[input]` | `cooccurrence_ids` | Accession(s) tested for co-occurrence |
| `[databases]` | `mmseqs2_db`, `hmmer_db`, `foldseek_db_path` | Local database paths |
| `[search]` | `use_mmseqs2`, `use_hmmer`, `use_foldseek` | Which search tools to run |
| `[search]` | `consensus_mode`, `min_tools` | How to combine tools (`union`/`consensus`/`intersection`) |
| `[search]` | `kingdom_filter` | Restrict to a domain, e.g. `Bacteria` (blank = no filter) |
| `[search]` | `primary_rbh_mode` | `off` / `annotate` / `filter` reciprocal-best-hit validation |
| `[align]` | `aligner` | `mafft` or `pagan2` |
| `[tree]` | `fasttree_bin` | `FastTree` or `FastTreeMP` |
| `[neighborhood]` | `cooccurrence_mode` | `search` or `rbh` |
| `[neighborhood]` | `neighborhood`, `ncbi_email` | Enable gene-distance analysis (requires an NCBI email) |

> **Do not commit `config/config.txt`** — it may contain your NCBI API key and
> cluster-specific paths. `.gitignore` already excludes it; only
> `config.example.txt` is tracked.

## Running

Two supported ways to run, depending on your cluster:

**A. Single-node, local cores** (simplest; recommended for small query sets).
Edit the SBATCH directives and environment-activation block in
`submit_cophy.sh`, then:

```bash
sbatch submit_cophy.sh
```

**B. Slurm executor** (one Slurm job per rule; better for large multi-query
runs). Edit `profiles/slurm/config.yaml` for your partition/account, then:

```bash
snakemake --profile profiles/slurm \
    --config cfg=config/config.txt outdir=results/<project_name>
```

To run locally without a cluster:

```bash
snakemake -s workflow/Snakefile --use-conda --conda-frontend conda -j 8 \
    --config cfg=config/config.txt outdir=results/<project_name>
```

## Outputs

Under `results/<project_name>/`:

| Path | Contents |
|---|---|
| `tree/species_tree.nwk` | The inferred species tree (Newick) |
| `supermatrix/supermatrix.fasta` | Concatenated trimmed alignment |
| `supermatrix/partitions.txt` | RAxML-style partition file |
| `cooccurrence.csv` | Per-species presence of each co-occurrence protein (+ distances) |
| `cooccurrence_collapsed.csv` | One row per species summary |
| `itol/` | iTOL annotation files (taxonomy strips, presence, gradients) |
| `groups/` | Per-group members tables and FASTAs |

## Repository layout

```
cophy/
├── workflow/
│   ├── Snakefile              # top-level workflow
│   ├── rules/                 # one .smk per stage
│   └── scripts/               # Python implementation of each step
├── envs/                      # per-rule conda environments
├── config/
│   └── config.example.txt     # annotated configuration template
├── profiles/slurm/            # Snakemake Slurm profile
├── submit_cophy.sh        # single-node driver script
├── LICENSE
└── README.md
```

## Development and AI use

CoPhy was designed and developed by the author. Anthropic's Claude (an AI
assistant) was used as a coding and debugging aid throughout development —
helping draft and refactor the Python scripts and Snakemake rules, work
through cluster and environment issues, and prepare documentation. All design
decisions, scientific choices (search tools and thresholds, orthology and
consensus logic, taxonomic filtering, tree-building parameters), and the
validation of outputs were made and reviewed by the author. The author is
responsible for the correctness of the code and results.

## Citing

If you use CoPhy, please cite this repository and the underlying tools:
MMseqs2, HMMER, foldseek, MAFFT, ClipKIT, FastTree2, and (for structure-based
trees) FoldTree. See each tool's documentation for the appropriate citation.

## License

MIT — see [LICENSE](LICENSE).
