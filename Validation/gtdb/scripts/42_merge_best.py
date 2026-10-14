#!/usr/bin/env python3
"""Assemble the best-hit tables for a re-drawn query set from old + new searches.

The matched clean comparator (40_match_clean_comparator.py) reuses every ORF
already searched as shadow_hi and searches only the rest. This writes one
best/ directory covering exactly wanted_v7.tsv: rows from the published best/
tables whose query is in the new wanted list (so the dropped shadow_hi rows go),
plus every row from the new searches. The null tables are candidate-only and
are copied unchanged. A query found in both sources is an error -- the two
searches are supposed to be disjoint -- rather than silently deduplicated.

  42_merge_best.py --wanted wanted_v7.tsv --old <...>/search_pooled/best \
      --new <...>/search_shard1/best,<...>/search_shard2/best --out best_v7
"""
import argparse
import glob
import os
import shutil
import sys
from pathlib import Path

import pandas as pd


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--wanted", required=True)
    ap.add_argument("--old", required=True)
    ap.add_argument("--new", required=True, help="comma-separated best/ dirs")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    w = pd.read_csv(args.wanted, sep="\t", dtype=str)
    keep = set(w.genome + "|" + w.input_id + "|" + w.orf_id + "|" + w.group)
    new_dirs = [d for d in args.new.split(",") if d]
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    for path in sorted(glob.glob(f"{args.old}/*.m8")):
        name = os.path.basename(path)
        if name.startswith("null."):
            shutil.copyfile(path, out / name)
            print(f"{name}: copied (candidate-only null)")
            continue
        n_old = n_kept = n_new = 0
        seen = set()
        tmp = out / (name + ".partial")
        with open(tmp, "w") as fo:
            with open(path) as fh:
                for line in fh:
                    n_old += 1
                    q = line.split("\t", 1)[0]
                    if q in keep:
                        fo.write(line)
                        seen.add(q)
                        n_kept += 1
            for d in new_dirs:
                p = Path(d) / name
                if not p.exists():
                    sys.exit(f"ERROR: {p} missing")
                with open(p) as fh:
                    for line in fh:
                        q = line.split("\t", 1)[0]
                        if q in seen:
                            sys.exit(f"ERROR: {q} in both old and new {name}")
                        if q not in keep:
                            sys.exit(f"ERROR: new hit {q} is not in the wanted list")
                        seen.add(q)
                        fo.write(line)
                        n_new += 1
        os.replace(tmp, out / name)
        print(f"{name}: old {n_old:,} -> kept {n_kept:,}; new {n_new:,}; "
              f"total {n_kept + n_new:,}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
