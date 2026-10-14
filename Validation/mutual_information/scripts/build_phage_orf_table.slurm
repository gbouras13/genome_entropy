#!/bin/bash
#SBATCH --job-name=phage_mi_build
#SBATCH --account=pawsey1018
#SBATCH --partition=work
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=0-02:00
#SBATCH -o %x-%j.out
#SBATCH -e %x-%j.err

# Build the per-ORF Parquet table, check labels/ORF floor against the current
# code on three GenBank chunks, and record provenance of phold.gbk.
#
# Usage: sbatch build_phage_orf_table.slurm WORK_DIRECTORY REPO PHOLD_GBK_GZ

set -euo pipefail
WORK=${1:?WORK_DIRECTORY}
REPO=${2:?REPO}
PHOLD=${3:?PHOLD_GBK_GZ}
PY=/scratch/pawsey1018/edwa0468/software/envs/logan/bin/python
SCRIPTS="${REPO}/Validation/mutual_information/scripts"
export PATH="$HOME/software/bin:$PATH"   # get_orfs

"$PY" "${SCRIPTS}/build_phage_orf_table.py" "${WORK}/results" \
    "${WORK}/phold_modernprost_50m_mutual_information.tsv" \
    "${WORK}/phage_orfs.parquet" --processes "${SLURM_CPUS_PER_TASK}"

"$PY" "${SCRIPTS}/check_phage_labels.py" --repo "${REPO}" \
    --orf-table "${WORK}/phage_orfs.parquet" \
    --chunks "${WORK}"/provenance/phold_0000.gbk.gz "${WORK}"/provenance/phold_0072.gbk.gz \
             "${WORK}"/provenance/phold_0144.gbk.gz \
    --out "${WORK}/label_check.tsv"

# phold.gbk provenance: records, LOCUS dates, CDS count, gene caller versions.
{
    echo "sha256 $(sha256sum "${PHOLD}" | cut -d' ' -f1)"
    zcat "${PHOLD}" | awk '
        /^LOCUS/ {records++; dates[$NF]++}
        /^     CDS / {cds++}
        /\/inference=/ {inf[$0]++}
        END {
            print "records", records; print "CDS", cds;
            for (d in dates) print "locus_date", d, dates[d];
            for (i in inf) print "inference", i, inf[i];
        }'
} > "${WORK}/phold_gbk_provenance.txt"
cat "${WORK}/phold_gbk_provenance.txt"
