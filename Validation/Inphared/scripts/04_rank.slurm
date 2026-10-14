#!/bin/bash
#SBATCH --job-name=phage_rank
#SBATCH --account=pawsey1018
#SBATCH --partition=work
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=115G
#SBATCH --time=04:00:00
#SBATCH -o %x-%j.out
#SBATCH -e %x-%j.err
#
# The re-ranking itself: scripts/rank_missing_orfs.py from this checkout, on
# the relabelled JSONs (02_relabel.py), with a chosen feature set.
#   sbatch 04_rank.slurm intrinsic
set -euo pipefail
source "${PHAGE_REPO:?set PHAGE_REPO to the checkout}/Validation/Inphared/scripts/env.sh"
fs="${1:-intrinsic}"
python -W ignore::FutureWarning "$PHAGE_REPO/scripts/rank_missing_orfs.py" \
    --json-dir "$PHAGE_WORK/genbank_split_json_relabelled" --output-dir "$PHAGE_WORK/rank_orfs_relabelled_$fs" \
    --folds 5 --minimum-probability 0.80 --minimum-aa-length 30 --feature-set "$fs" --device cpu
