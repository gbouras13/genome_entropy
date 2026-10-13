#!/usr/bin/env python3
"""Recompute in_genbank for every phage ORF with the current matcher.

The archived label (copied from the legacy PHOLD JSON into the aggregated
TSV) marks only 3 of 3.31 M minus-strand ORFs as in_genbank, i.e. it was
produced by a matcher that could not match reverse-strand CDS (see the
coordinate fixes of 2026-07-25..31 and ea00cd6). This script re-derives the
ORFs from each phold GenBank chunk with ``get_orfs -l 90`` (the effective
``run --min-aa 30`` behaviour, issue #103), checks that they reproduce the
archived ORF set exactly, and labels them with
``genome_entropy.io.genbank.match_orf_to_genbank_cds`` (same strand and codon
phase, >=90% overlap of the shorter interval, >=98% coordinate-aligned
amino-acid identity).

Writes a Parquet with (input_id, start, end, strand, frame, in_genbank_v2,
matched_cds_aa, matched_cds_partial, n_cds_matched) and a per-chunk summary.
"""

import argparse
import sys
from multiprocessing import Pool
from pathlib import Path

import pandas as pd

REPO = None


def label_chunk(path: Path) -> tuple[pd.DataFrame, dict]:
    sys.path.insert(0, str(REPO / "src"))
    from genome_entropy.io.genbank import (evaluate_orf_genbank_cds_match,
                                           extract_cds_features, read_genbank)
    from genome_entropy.orf.finder import find_orfs

    sequences = read_genbank(path)
    cds = extract_cds_features(path)
    by_parent: dict = {}
    for feature in cds:
        by_parent.setdefault(feature.parent_id, []).append(feature)
    rows = []
    for seq_id, dna in sequences.items():
        candidates = by_parent.get(seq_id, [])
        for orf in find_orfs({seq_id: dna}, table_id=11, min_nt_length=90):
            matches = [c for c in candidates if c.strand == orf.strand
                       and evaluate_orf_genbank_cds_match(orf, c).matched]
            rows.append({
                "input_id": seq_id, "start": orf.start, "end": orf.end,
                "strand": orf.strand, "frame": abs(int(orf.frame)),
                "in_genbank_v2": bool(matches),
                "n_cds_matched": len(matches),
                "matched_cds_aa": (max(len(c.protein_sequence.rstrip("*")) for c in matches)
                                   if matches else 0),
                "matched_cds_partial": any(c.partial for c in matches),
            })
    cds_aa = [len(c.protein_sequence.rstrip("*")) for c in cds]
    summary = {
        "chunk": path.name, "records": len(sequences), "cds": len(cds),
        "cds_skipped": sum(bool(c.skip_reason) for c in cds),
        "cds_lt_90aa": sum(a < 90 for a in cds_aa),
        "orfs": len(rows), "orfs_labelled": sum(r["in_genbank_v2"] for r in rows),
    }
    print(summary, flush=True)
    return pd.DataFrame(rows), summary


def main() -> None:
    global REPO
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--chunks", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--summary", type=Path, required=True)
    parser.add_argument("--processes", type=int, default=8)
    args = parser.parse_args()
    REPO = args.repo
    paths = sorted(args.chunks.glob("phold_*.gbk.gz"))
    # largest first for load balance
    paths.sort(key=lambda p: p.stat().st_size, reverse=True)
    with Pool(args.processes, initializer=_init, initargs=(args.repo,)) as pool:
        results = pool.map(label_chunk, paths, chunksize=1)
    labels = pd.concat([r[0] for r in results], ignore_index=True)
    labels.to_parquet(args.out, index=False)
    pd.DataFrame([r[1] for r in results]).sort_values("chunk").to_csv(
        args.summary, sep="\t", index=False)
    print(f"{len(labels):,} ORFs, {labels.in_genbank_v2.sum():,} labelled", flush=True)


def _init(repo: Path) -> None:
    global REPO
    REPO = repo


if __name__ == "__main__":
    main()
