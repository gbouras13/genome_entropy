#!/bin/bash
#SBATCH --job-name=phage_mi_relabel
#SBATCH --account=pawsey1018
#SBATCH --partition=work
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=0-03:00
#SBATCH -o %x-%j.out
#SBATCH -e %x-%j.err

# Recompute in_genbank for all phage ORFs with the current matcher.
# Usage: sbatch relabel_phage_orfs.slurm WORK_DIRECTORY REPO
# Expects WORK/chunks/phold_*.gbk.gz (fetch_phage_mi_chunks.slurm).

set -euo pipefail
WORK=${1:?WORK_DIRECTORY}
REPO=${2:?REPO}
PY=/scratch/pawsey1018/edwa0468/software/envs/logan/bin/python
export PATH="$HOME/software/bin:$PATH"   # get_orfs
"$PY" "${REPO}/Validation/mutual_information/scripts/relabel_phage_orfs.py" \
    --repo "${REPO}" --chunks "${WORK}/chunks" \
    --out "${WORK}/phage_orf_labels_v2.parquet" \
    --summary "${WORK}/relabel_chunk_summary.tsv" \
    --processes "${SLURM_CPUS_PER_TASK}"
