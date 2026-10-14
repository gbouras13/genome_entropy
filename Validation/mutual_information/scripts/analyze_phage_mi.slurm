#!/bin/bash
#SBATCH --job-name=phage_mi_analyze
#SBATCH --account=pawsey1018
#SBATCH --partition=work
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=0-03:00
#SBATCH -o %x-%j.out
#SBATCH -e %x-%j.err

# Run the phage MI / length-control / AUC analysis on the per-ORF table.
#
# Usage: sbatch analyze_phage_mi.slurm PHAGE_ORFS.parquet REPO [extra args]

set -euo pipefail
TABLE=${1:?PHAGE_ORFS.parquet}
REPO=${2:?REPO}
shift 2
PY=/scratch/pawsey1018/edwa0468/software/envs/logan/bin/python
OUT="${REPO}/Validation/mutual_information"
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK}"

"$PY" "${OUT}/scripts/analyze_phage_mi.py" "${TABLE}" \
    --tables "${OUT}/tables" --figures "${OUT}/figures" "$@"
