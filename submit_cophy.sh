#!/bin/bash
# Single-node driver for CoPhy: runs Snakemake with local cores inside
# one Slurm allocation, rather than using the Slurm executor to submit each
# rule as its own job. This is simpler and robust on clusters where the
# Slurm executor is unreliable; for a single query set it is usually fast
# enough. For large multi-query runs, use `--profile profiles/slurm` instead
# (see README).
#
# Edit the SBATCH directives and the environment-activation block for your
# cluster, then:  sbatch submit_cophy.sh
#SBATCH --job-name=cophy
#SBATCH --partition=batch
#SBATCH --time=1-00:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=120000
#SBATCH --output=logs/cophy-%j.out
#SBATCH --error=logs/cophy-%j.err

set -eo pipefail

# --- pick the config for this run ---
CONFIG="config/config.txt"

# Run from the repository root (the directory containing workflow/).
cd "$(dirname "$0")"

# Put tool scratch on fast local disk if your shared filesystem is slow/full.
export TMPDIR="${TMPDIR:-/tmp/$USER}"
mkdir -p "$TMPDIR"

# --- activate the environment that has snakemake ---
# Replace this block with however your cluster provides conda + snakemake.
# Example:
#   module load miniconda
#   source "$(conda info --base)/etc/profile.d/conda.sh"
#   conda activate /path/to/snakemake_env

# --- derive results dir from project_name in the config ---
PROJECT_NAME=$(awk -F'=' '
    /^[[:space:]]*project_name[[:space:]]*=/ {
        sub(/#.*/, "", $2); gsub(/^[[:space:]]+|[[:space:]]+$/, "", $2)
        print $2; exit
    }' "$CONFIG")
if [[ -z "$PROJECT_NAME" ]]; then
    echo "ERROR: project_name missing from $CONFIG" >&2; exit 1
fi
OUTDIR="results/${PROJECT_NAME}"
mkdir -p logs "$OUTDIR"

# Validate the config before launching.
python workflow/scripts/config.py "$CONFIG" > /dev/null && echo "config: OK"

# Clear any stale lock from a previously killed run.
snakemake -s workflow/Snakefile --config cfg="$CONFIG" outdir="$OUTDIR" \
    --unlock 2>/dev/null || true

snakemake \
    -s workflow/Snakefile \
    --config cfg="$CONFIG" outdir="$OUTDIR" \
    --use-conda --conda-frontend conda \
    -j "${SLURM_CPUS_PER_TASK:-8}" \
    --rerun-triggers mtime \
    --rerun-incomplete
echo "finished: $(date)"
