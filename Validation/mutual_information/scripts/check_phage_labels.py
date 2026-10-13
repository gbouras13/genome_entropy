#!/usr/bin/env python3
"""Re-derive ORFs and in_genbank labels for phold GenBank chunks with the
current genome_entropy code, and compare them with the archived run.

For each chunk this reports:
- whether the current ORF caller (``run --min-aa 30`` semantics, i.e.
  ``get_orfs -l 90``) reproduces the archived ORF set (coordinates);
- the shortest archived protein;
- agreement between the archived ``in_genbank`` label (TSV, from the legacy
  PHOLD JSON) and the label recomputed with the current matcher
  (same strand/phase, >=90% overlap of shorter interval, >=98% identity).

Usage:
  check_phage_labels.py --repo REPO --tsv-index ORFS.parquet \
      --chunks phold_0000.gbk.gz ... --out check.tsv
"""

import argparse
import sys
from pathlib import Path

import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--orf-table", type=Path, required=True,
                        help="Parquet written by build_phage_orf_table.py")
    parser.add_argument("--chunks", type=Path, nargs="+", required=True)
    parser.add_argument("--min-aa", type=int, default=30)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()

    sys.path.insert(0, str(args.repo / "src"))
    from genome_entropy.io.genbank import (extract_cds_features,
                                           match_orf_to_genbank_cds,
                                           read_genbank)
    from genome_entropy.orf.finder import find_orfs

    archived = pd.read_parquet(
        args.orf_table,
        columns=["input_id", "start", "end", "strand", "frame",
                 "protein_length", "in_genbank"],
    )

    rows = []
    for chunk in args.chunks:
        sequences = read_genbank(chunk)
        cds = extract_cds_features(chunk)
        by_parent = {}
        for feature in cds:
            by_parent.setdefault(feature.parent_id, []).append(feature)
        records = []
        for seq_id, dna in sequences.items():
            for orf in find_orfs({seq_id: dna}, table_id=11,
                                 min_nt_length=args.min_aa * 3):
                records.append({
                    "input_id": seq_id, "start": orf.start, "end": orf.end,
                    "strand": orf.strand, "frame": abs(int(orf.frame)),
                    "new_label": match_orf_to_genbank_cds(
                        orf, by_parent.get(seq_id, [])),
                    "new_protein_length": len(orf.aa_sequence.rstrip("*")),
                })
        new = pd.DataFrame(records)
        old = archived[archived.input_id.isin(set(sequences))]
        merged = old.merge(new, on=["input_id", "start", "end", "strand", "frame"],
                           how="outer", indicator=True)
        both = merged[merged._merge == "both"]
        agree = (both.in_genbank.astype(bool) == both.new_label.astype(bool))
        row = {
            "chunk": chunk.name,
            "records": len(sequences),
            "cds_features": len(cds),
            "archived_orfs": len(old),
            "recomputed_orfs": len(new),
            "orfs_in_both": len(both),
            "archived_only": int((merged._merge == "left_only").sum()),
            "recomputed_only": int((merged._merge == "right_only").sum()),
            "min_archived_protein_aa": int(old.protein_length.min()),
            "min_recomputed_protein_aa": int(new.new_protein_length.min()),
            "label_agreement": float(agree.mean()),
            "archived_true": int(both.in_genbank.astype(bool).sum()),
            "recomputed_true": int(both.new_label.astype(bool).sum()),
            "archived_true_recomputed_false": int(
                (both.in_genbank.astype(bool) & ~both.new_label.astype(bool)).sum()),
            "archived_false_recomputed_true": int(
                (~both.in_genbank.astype(bool) & both.new_label.astype(bool)).sum()),
        }
        print(row, flush=True)
        rows.append(row)
    pd.DataFrame(rows).to_csv(args.out, sep="\t", index=False)


if __name__ == "__main__":
    main()
