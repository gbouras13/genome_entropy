#!/bin/bash
#SBATCH --job-name=phage_mi_fetch_chunks
#SBATCH --account=pawsey1018
#SBATCH --partition=copy
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=0-02:00
#SBATCH -o %x-%j.out
#SBATCH -e %x-%j.err

# Copy the 145 phold GenBank chunks (the exact inputs of the MI re-encoding).
# Usage: sbatch fetch_phage_mi_chunks.slurm WORK_DIRECTORY
set -euo pipefail
WORK=${1:?WORK_DIRECTORY}
module load rclone/1.68.1
rclone copy genome_entropy:General/Validation/mutual_information/chunks "${WORK}/chunks" \
    --transfers 4 --stats 60s -v
rclone check genome_entropy:General/Validation/mutual_information/chunks "${WORK}/chunks" --one-way --size-only
