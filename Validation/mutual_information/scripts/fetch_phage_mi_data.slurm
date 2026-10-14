#!/bin/bash
#SBATCH --job-name=phage_mi_fetch
#SBATCH --account=pawsey1018
#SBATCH --partition=copy
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=0-06:00
#SBATCH -o %x-%j.out
#SBATCH -e %x-%j.err

# Copy the phage mutual-information inputs from the OneDrive remote to scratch.
#
# Usage: sbatch fetch_phage_mi_data.slurm WORK_DIRECTORY
#
# Fetches the aggregated TSV (784 MB) and the 145 per-chunk unified JSON
# outputs (results/phold_*.json.gz; the FASTA inputs are skipped).

set -euo pipefail

WORK=${1:?Usage: $0 WORK_DIRECTORY}
REMOTE="genome_entropy:General/Validation/mutual_information"

module load rclone/1.68.1
mkdir -p "${WORK}/results"

rclone copyto "${REMOTE}/phold_modernprost_50m_mutual_information.tsv" \
    "${WORK}/phold_modernprost_50m_mutual_information.tsv" --stats 60s -v
rclone copy "${REMOTE}/results" "${WORK}/results" --include 'phold_*.json.gz' \
    --transfers 4 --stats 60s -v
rclone copy "${REMOTE}/in_genbank_pair_analysis" "${WORK}/in_genbank_pair_analysis" -v

rclone check "${REMOTE}/results" "${WORK}/results" --include 'phold_*.json.gz' --one-way --size-only
ls "${WORK}/results" | wc -l
