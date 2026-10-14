# GMSC smORFs → `genome_entropy`, on NCI Gadi

Protein, 3Di and 12-state entropy for **43,726,909 small ORFs** from the Global
Microbial smORF Catalogue — 33,647,794 positives and 10,079,115 negatives —
and how well those entropies separate real small proteins from spurious ones.

This directory is a **record of what ran**, not a portable pipeline. The
scripts are committed exactly as they were executed, so they carry Gadi paths
(`/g/data/ob80/...`, `/scratch/ob80/...`, `-P ob80`) and the conda prefixes
used. The write-up is `projects/genome_entropy_logan/smorfs/manuscript.md` in
`linsalrob/agent-handoffs` (commits bbb3f54, def54d0).

## Inputs

`positive.tar.gz` and `negative.tar.gz` from the GMSC authors, one TSV each
with columns `seqid, upstream, orf, downstream`. There are no GMSC IDs,
evidence classes or COG annotations in them, so no COG analysis is possible
from these files. They are sorted by length within each file: never test on a
file head.

| | positive | negative |
|---|---:|---:|
| sequences | 33,647,794 | 10,079,115 |
| aa length, min–max | 16–99 | 10–99 |
| aa length, median | 80 | 70 |

There is no 90 aa floor here, unlike the GTDB run (see issue #103): the set is
the curated catalogue as supplied.

## Versions

| | |
|---|---|
| `genome_entropy` | 0.3.0, installed 2026-09-01 16:34 from this repository. Source is identical to tag `v0.3.0` (45dd283); no later commit touches `src/`. |
| model | `gbouras13/modernprost-50M` (`ModernProstThreeDiEncoder`), 3Di + 12-state heads |
| encoding size | 50,000 (the 10k/50k pilot is `pilot.pbs`) |
| genetic code | table 11 |
| environment | conda `/g/data/ob80/re3494/conda/genome_entropy`: Python 3.11, torch 2.13.0+cu126, transformers 5.16.1 |
| MMseqs2 | 14.7e284 (`/g/data/ob80/re3494/conda/bioinformatics`) |
| hardware | `gpuvolta`, V100 32 GB |

## Order of execution

| step | script | job | SU | walltime |
|---|---|---|---:|---:|
| split tarballs into 1M-row chunks | `split.pbs` | 180700231 | 0.12 | 0:01:51 |
| pilot, encoding size 10k vs 50k | `pilot.pbs` → `encode_entropy.py` | 180700319 | 2.70 | 0:04:30 |
| encode all chunks, 10-way strided array | `full.pbs` → `encode_entropy.py` | 180703415[0-9] | 365.30 | 0:54–1:11 each |
| merge, summary, first plots | `merge.pbs` → `merge_plot.py` | 180703537 | 6.32 | 0:07:54 |
| amino-acid vs 3Di density | `density.pbs` → `density_aa_vs_3di.py` | 181044341 | 0.27 | 0:01:00 |
| 2.8-bit stripe analysis | `stripes.pbs` → `stripe_states.py` | 181047298 | 0.28 | 0:00:32 |
| manuscript tables, figures, models | `manuscript.pbs` → `manuscript_stats.py` | 181048769 | 3.03 | 0:05:41 |
| MMseqs2 clustering, 50% and 30% identity | `cluster.pbs` | 181057870 | 23.68 | 0:14:48 |
| cluster-split model evaluation | `cluster_models.pbs` → `cluster_models.py` | 181057871 | 9.49 | 0:08:54 |
| **total** | | | **411.19** | |

Every job exited 0. The PBS logs are in [`logs/`](logs/); the originals were
on `/scratch/ob80/re3494/smORFs/logs/`, which is purged, and a copy is on the
OneDrive remote at `genome_entropy:General/smORFs/results/logs/`.

## Outputs

On Gadi under `/g/data/ob80/re3494/Projects/smORFs/results/`, mirrored to
`genome_entropy:General/smORFs/results/`. None of it is in git.

- `smorf_entropy.tsv.gz` (913 MB): one row per smORF; `smorf_entropy_summary.tsv`
- `states/`: per-chunk 3Di and 12-state strings
- `manuscript/`: tables 2–3 (single-feature AUC, models) and figures
- `stripes/`, `density_aa_vs_3di.*`: the 2.8-bit stripe (manuscript §3.6)
- `clusters/`: MMseqs2 cluster TSVs, `cluster_split_models.tsv`

## Headline numbers

Single features, rank AUC over all 43.7 M (`manuscript/table2_single_feature_auc.tsv`):

| feature | AUC |
|---|---:|
| length | 0.640 |
| 3Di entropy | 0.605 |
| 12-state entropy | 0.602 |

All five features (length plus nt, aa, 3Di and 12-state entropy) in a
gradient-boosted model with 2 M train and 2 M test
(`manuscript/table3_models.tsv`, `clusters/cluster_split_models.tsv`):

| split | AUC | cluster-weighted |
|---|---:|---:|
| random | 0.730 | 0.700 |
| 50% identity clusters | 0.728 | 0.697 |
| 30% identity clusters | 0.721 | 0.695 |

There are 4.00 M clusters at 50% identity and 2.40 M at 30%.
