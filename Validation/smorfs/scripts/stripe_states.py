#!/usr/bin/env python
"""Characterise the horizontal 3Di entropy stripes at log2(2), log2(3) and log2(7) bits.

For each stripe, strings within half a density-plot bin of the stripe are compared with
two flanking control bands of the same width. Reported per band: the number of distinct
3Di states used, the most common state sets, and the most repeated 3Di / amino acid strings.

Usage: stripe_states.py STATES_DIR OUT_DIR [NPROC]
"""
import csv
import glob
import gzip
import math
import os
import sys
from collections import Counter, defaultdict
from multiprocessing import Pool

STRIPES = (2, 3, 7)
HALF_BIN = 4.4 / 400 / 2  # half a bin of the density plot
FLANK = 0.03              # control band centres sit this far below and above the stripe
TOP = 15
DATASETS = ("positive", "negative")


def bands():
    for k in STRIPES:
        c = math.log2(k)
        for name, centre in (("below", c - FLANK), ("stripe", c), ("above", c + FLANK)):
            yield k, name, centre


BANDS = list(bands())


def scan(path):
    base = os.path.basename(path).split(".")[0]
    dataset = "positive" if base.startswith("pos") else "negative"
    tot = 0
    n_states, sets, tdis, aas = (defaultdict(Counter) for _ in range(4))
    rows = []
    with gzip.open(path, "rt") as fh:
        r = csv.reader(fh, delimiter="\t")
        next(r)
        for sid, aa, tdi, _ in r:
            if not tdi:
                continue
            tot += 1
            c = Counter(tdi)
            n = len(tdi)
            h = -sum(v / n * math.log2(v / n) for v in c.values())
            for k, name, centre in BANDS:
                if abs(h - centre) <= HALF_BIN:
                    key = (k, name)
                    n_states[key][len(c)] += 1
                    sets[key]["".join(sorted(c))] += 1
                    if name == "stripe":
                        tdis[key][tdi] += 1
                        aas[key][aa] += 1
                        rows.append((sid, dataset, k, f"{h:.6f}", len(c), "".join(sorted(c)), n, aa, tdi))
    # Keep only the heads of the string counters; full counts would not fit through the pool.
    for d in (tdis, aas):
        for key in d:
            d[key] = Counter(dict(d[key].most_common(200)))
    return dataset, tot, n_states, sets, tdis, aas, rows


def main():
    src, out = sys.argv[1], sys.argv[2]
    nproc = int(sys.argv[3]) if len(sys.argv) > 3 else os.cpu_count()
    os.makedirs(out, exist_ok=True)
    files = sorted(glob.glob(os.path.join(src, "*.states.tsv.gz")))
    print(f"{len(files)} states files, {nproc} processes", flush=True)

    tot = Counter()
    n_states, sets, tdis, aas = ({ds: defaultdict(Counter) for ds in DATASETS} for _ in range(4))
    members = os.path.join(out, "stripe_members.tsv.gz")
    with Pool(nproc) as pool, gzip.open(members + ".tmp", "wt", compresslevel=6) as fo:
        fo.write("seqid\tdataset\tk\tthree_di_entropy\tn_states\tstate_set\tlength\taa\tthree_di\n")
        for ds, n, ns, st, td, aa, rows in pool.imap_unordered(scan, files):
            tot[ds] += n
            for acc, part in ((n_states, ns), (sets, st), (tdis, td), (aas, aa)):
                for key, ctr in part.items():
                    acc[ds][key].update(ctr)
            for row in rows:
                fo.write("\t".join(map(str, row)) + "\n")
            print(f"done {ds}: {n} seqs, {len(rows)} stripe members", flush=True)
    os.rename(members + ".tmp", members)

    with open(os.path.join(out, "stripe_bands.tsv"), "w") as fs:
        fs.write("dataset\tk\tband\tcentre_bits\tn_dataset\tn_band\tper_million\t"
                 "pct_with_k_states\tpct_with_fewer_states\ttop_state_set\tpct_top_state_set\n")
        for ds in DATASETS:
            for k, name, centre in BANDS:
                key = (k, name)
                nb = sum(n_states[ds][key].values())
                exact = n_states[ds][key][k]
                fewer = sum(v for s, v in n_states[ds][key].items() if s < k)
                top, topn = (sets[ds][key].most_common(1) or [("", 0)])[0]
                pct = lambda x: f"{100 * x / nb:.2f}" if nb else "0"
                fs.write(f"{ds}\t{k}\t{name}\t{centre:.4f}\t{tot[ds]}\t{nb}\t{1e6 * nb / tot[ds]:.1f}\t"
                         f"{pct(exact)}\t{pct(fewer)}\t{top}\t{pct(topn)}\n")

    with open(os.path.join(out, "stripe_top_patterns.tsv"), "w") as ft:
        ft.write("dataset\tk\tfeature\trank\tvalue\tcount\tpct_of_stripe\n")
        for ds in DATASETS:
            for k in STRIPES:
                key = (k, "stripe")
                nb = sum(n_states[ds][key].values())
                for name, ctr in (("n_states", n_states[ds][key]), ("state_set", sets[ds][key]),
                                  ("three_di", tdis[ds][key]), ("aa", aas[ds][key])):
                    for i, (v, c) in enumerate(ctr.most_common(TOP), 1):
                        ft.write(f"{ds}\t{k}\t{name}\t{i}\t{v}\t{c}\t{100 * c / nb:.2f}\n")


if __name__ == "__main__":
    main()
