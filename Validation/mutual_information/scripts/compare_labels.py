#!/usr/bin/env python3
"""Compare the archived (legacy) in_genbank label with the recomputed one.

Writes small TSVs: label cross-tabulation by strand and contig-edge status,
and, for recomputed positives, how much of the ORF the matched CDS covers
(stop-to-stop ORFs include upstream extension; the >=90%-of-shorter-interval
rule lets a short CDS label a much longer ORF).
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("orfs", type=Path)
    parser.add_argument("labels", type=Path)
    parser.add_argument("chunk_summary", type=Path)
    parser.add_argument("--tables", type=Path, required=True)
    args = parser.parse_args()
    args.tables.mkdir(parents=True, exist_ok=True)
    keys = ["input_id", "start", "end", "strand", "frame"]
    orfs = pd.read_parquet(args.orfs, columns=[*keys, "input_dna_length",
                                               "protein_length", "in_genbank"])
    orfs["strand"] = orfs.strand.astype(str)
    labels = pd.read_parquet(args.labels)
    d = orfs.merge(labels, on=keys, how="outer", indicator=True, validate="one_to_one")
    print("ORF-set join:", d._merge.value_counts().to_dict())
    d["contig_edge"] = d.end > d.input_dna_length
    d["legacy"] = d.in_genbank.astype(bool)
    d["v2"] = d.in_genbank_v2.astype(bool)
    cross = (d.groupby(["strand", "contig_edge", "legacy", "v2"]).size()
             .rename("orfs").reset_index())
    cross.to_csv(args.tables / "label_legacy_vs_v2.tsv", sep="\t", index=False)
    print(cross.to_string())

    pos = d[d.v2]
    cover = pos.matched_cds_aa / pos.protein_length
    cover_rows = {
        "v2_positive_orfs": len(pos),
        "matched_cds_lt_90aa": int((pos.matched_cds_aa < 90).sum()),
        "cds_covers_lt_50pct_of_orf": int((cover < 0.5).sum()),
        "cds_covers_lt_90pct_of_orf": int((cover < 0.9).sum()),
        "median_cds_fraction_of_orf": float(cover.median()),
        "orfs_matching_gt1_cds": int((pos.n_cds_matched > 1).sum()),
        "matched_cds_partial": int(pos.matched_cds_partial.sum()),
    }
    summary = pd.read_csv(args.chunk_summary, sep="\t")
    for column in ("records", "cds", "cds_skipped", "cds_lt_90aa", "orfs", "orfs_labelled"):
        cover_rows[f"genbank_{column}"] = int(summary[column].sum())
    pd.Series(cover_rows).rename("value").to_csv(
        args.tables / "label_v2_summary.tsv", sep="\t", index_label="quantity")
    print(pd.Series(cover_rows).to_string())


if __name__ == "__main__":
    main()
