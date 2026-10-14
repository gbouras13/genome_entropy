#!/usr/bin/env python3
"""Emit the FULL shadow_hi pool of one chunk, with its reading-frame class.

The §6.1 clean comparator (antisense + same-strand-frameshift shadows) was
taken as a SUBSET of the shadow_hi arm after that arm had been matched 1:1 to
the candidates on length and 3Di. Same-frame and frame-undefined shadows are
long and high-3Di, so removing them left a comparator that is shorter and
lower in 3Di than the candidates (bacteria: 122 vs 134 aa, 3Di 2.68 vs 2.83,
4.6% vs 20% in the candidates' top 3Di quintile). Re-matching AFTER the frame
filter needs the whole shadow population, and 10_missed_genes.py kept only a
2,000-row random sample of it per chunk.

This script re-runs 10_missed_genes.py's classification unchanged (its
functions are imported, not copied), takes every shadow_hi ORF, and assigns
the frame class with 20_orf_context.py's context_for_contig/frame_class --
the same code that labelled the published comparator. Nothing about the
candidate/shadow split changes; the stats file written next to the pool is
checked against the published one so that a drift would fail loudly.

  39_shadow_pool.py --chunk-tsv entropy_rows/bac/bac_000.tsv.gz \
      --annotation-status genome_cds_counts_bac.tsv \
      --cds-intervals cds_intervals/bac --published missed_genes/bac \
      --out-dir missed_genes/matched_clean/bac/pool
"""
import argparse
import importlib.util
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent


def _load(name, filename):
    spec = importlib.util.spec_from_file_location(name, HERE / filename)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


m10 = _load("missed_genes", "10_missed_genes.py")
m20 = _load("orf_context", "20_orf_context.py")

POOL_COLS = ["domain", "chunk", "genome", "input_id", "orf_id", "start", "end",
             "strand", "aa_length", "protein_entropy", "three_di_entropy",
             "truncated", "g_start", "g_end", "cds_frame_class",
             "cds_frame_reference", "cds_overlap_id"]
CHECK_COLS = ["candidates", "unmatched_hi_shadow", "unmatched_hi",
              "annotated_cds", "intergenic_lo", "unannot_hi"]


def load_cds(path):
    """Deposited CDS parts, prepared exactly as 20_orf_context.py main() does."""
    cds = pd.read_csv(path, sep="\t",
                      dtype={"genome": str, "contig": str, "cds_id": str})
    n_parts = (cds.genome + "@" + cds.cds_id).value_counts()
    cds["n_parts"] = (cds.genome + "@" + cds.cds_id).map(n_parts).astype(int)
    cds["frame_reference"] = ((cds.n_parts == 1)
                              & (((cds.end - cds.start) % 3) == 0))
    return cds.sort_values(["genome", "contig", "start"], kind="stable")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--chunk-tsv", required=True)
    ap.add_argument("--annotation-status", required=True)
    ap.add_argument("--cds-intervals", required=True)
    ap.add_argument("--published", required=True,
                    help="directory holding the published <tag>.stats.tsv, "
                         "which this run must reproduce")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--threshold", type=float, default=m10.THRESH)
    args = ap.parse_args()

    tag = os.path.basename(args.chunk_tsv).split(".tsv")[0]
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    out = out_dir / f"{tag}.pool.tsv.gz"

    df, n_bad = m10.load_chunk(args.chunk_tsv)
    _, with_cds = m10.load_annotation_status(args.annotation_status)
    groups, stats = m10.classify(df, with_cds, args.cds_intervals, tag,
                                 args.threshold)
    del df

    pub = pd.read_csv(Path(args.published) / f"{tag}.stats.tsv", sep="\t")
    drift = [f"{c}: {stats[c]} vs published {int(pub[c].iloc[0])}"
             for c in CHECK_COLS if int(pub[c].iloc[0]) != stats[c]]
    if drift:
        print(f"ERROR: {tag}: classification differs from the published "
              "split -- " + "; ".join(drift), file=sys.stderr)
        return 1

    sh = groups["shadow_hi"].reset_index(drop=True)
    if len(sh) != stats["unmatched_hi_shadow"]:
        print(f"ERROR: {tag}: {len(sh)} shadow rows, stats say "
              f"{stats['unmatched_hi_shadow']}", file=sys.stderr)
        return 1

    cds = load_cds(Path(args.cds_intervals) / f"{tag}.tsv")
    by_contig = {k: v for k, v in cds.groupby(["genome", "contig"], sort=False)}
    empty = cds.iloc[:0]
    parts = []
    for (genome, contig), g in sh.groupby(["genome", "input_id"], sort=False):
        g = g.reset_index(drop=True)
        ctx = m20.context_for_contig(g, by_contig.get((genome, contig), empty))
        g["overlaps_cds"] = ctx["overlaps_cds"]
        for col in ("cds_frame_class", "cds_frame_reference", "cds_overlap_id"):
            g[col] = ctx[col]
        parts.append(g)
    pool = (pd.concat(parts, ignore_index=True) if parts
            else sh.assign(overlaps_cds=[], cds_frame_class=[],
                           cds_frame_reference=[], cds_overlap_id=[]))

    # The classifier and the context code must agree that every shadow
    # overlaps a deposited CDS (the same self-check 20_orf_context.py runs).
    n_off = int((~pool.overlaps_cds.astype(bool)).sum())
    if n_off:
        print(f"ERROR: {tag}: {n_off:,} of {len(pool):,} shadows overlap no "
              "CDS in the context code", file=sys.stderr)
        return 1

    tmp = out.with_suffix(".partial")
    pool[POOL_COLS].to_csv(tmp, sep="\t", index=False, compression="gzip")
    os.replace(tmp, out)
    vc = pool.cds_frame_class.value_counts()
    print(f"{tag}: {len(pool):,} shadows; "
          + "; ".join(f"{k}={v}" for k, v in vc.items()))
    return 0


if __name__ == "__main__":
    sys.exit(main())
