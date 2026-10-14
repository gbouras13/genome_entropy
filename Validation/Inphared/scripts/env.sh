# Sourced by every Inphared ranking job: paths and the Python environment.
export PHAGE_WORK="${PHAGE_WORK:-/scratch/$PAWSEY_PROJECT/$USER/genome_entropy/work/phage_rank}"
export PHAGE_ENV="${PHAGE_ENV:-/scratch/$PAWSEY_PROJECT/$USER/software/envs/phage_rank}"
export FSENV="${FSENV:-/scratch/$PAWSEY_PROJECT/$USER/software/envs/foldseek}"
export FSDB="${FSDB:-/scratch/$PAWSEY_PROJECT/$USER/Logan/foldseek_db}"
# The checkout these scripts live in (rank_missing_orfs.py and src/ come from here).
export PHAGE_REPO="${PHAGE_REPO:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)}"
export PHAGE_SCRIPTS="$PHAGE_REPO/Validation/Inphared/scripts"
eval "$(conda shell.bash hook)"
conda activate "$PHAGE_ENV"
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
