#!/usr/bin/env python
"""Homology-aware evaluation: cluster statistics and models trained/tested on disjoint clusters.

Usage: cluster_models.py ENTROPY_TSV_GZ OUT_DIR CLUSTER_TSV_GZ [CLUSTER_TSV_GZ ...]
Each cluster TSV is MMseqs2 createtsv output (representative<TAB>member).
"""
import os
import sys

import numpy as np
import pandas as pd
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.metrics import roc_auc_score

from manuscript_stats import ENTROPIES, FEATURES, LENGTH_BINS, MODELS, SEED

N_TRAIN = N_TEST = 2_000_000


def auc(y, p, w=None):
    return roc_auc_score(y, p, sample_weight=w) if 0 < y.sum() < len(y) else np.nan


def cluster_stats(df, name):
    g = df.groupby("cluster", sort=False)
    size = g.size()
    npos = g["y"].sum()
    mixed = (npos > 0) & (npos < size)
    rows = {"clustering": name, "n_sequences": len(df), "n_clusters": len(size),
            "n_singletons": int((size == 1).sum()), "largest_cluster": int(size.max()),
            "pct_seqs_in_clusters_ge10": 100 * size[size >= 10].sum() / len(df),
            "n_mixed_clusters": int(mixed.sum()), "pct_seqs_in_mixed_clusters": 100 * size[mixed].sum() / len(df)}
    for lab, yv in (("positive", 1), ("negative", 0)):
        sub = df[df.y == yv]
        rows[f"n_clusters_{lab}"] = sub.cluster.nunique()
        rows[f"seqs_per_cluster_{lab}"] = len(sub) / sub.cluster.nunique()
    top = size.sort_values(ascending=False).head(20)
    topdf = pd.DataFrame({"clustering": name, "representative": top.index, "size": top.values,
                          "n_positive": npos[top.index].values})
    return rows, topdf


def split(df, how, rng):
    if how == "random":
        n = min(N_TRAIN, len(df) // 2)
        idx = rng.permutation(len(df))
        return idx[:n], idx[n:2 * n]
    # Whole clusters go to train or test; sample sequences within each side afterwards.
    clusters = df[how].unique()
    test_clusters = set(rng.permutation(clusters)[: len(clusters) // 2])
    is_test = df[how].isin(test_clusters).to_numpy()
    tr, te = np.flatnonzero(~is_test), np.flatnonzero(is_test)
    return rng.permutation(tr)[:N_TRAIN], rng.permutation(te)[:N_TEST]


def main():
    src, out, *clusterings = sys.argv[1:]
    os.makedirs(out, exist_ok=True)
    df = pd.read_csv(src, sep="\t", usecols=["seqid", "dataset", *FEATURES],
                     dtype={"dataset": "category", "aa_length": "int16", **{c: "float32" for c in ENTROPIES}})
    df["y"] = (df.dataset == "positive").astype(np.int8)
    names = []
    stats, tops = [], []
    for path in clusterings:
        name = os.path.basename(path).split(".")[0]
        cl = pd.read_csv(path, sep="\t", header=None, names=["rep", "seqid"])
        codes, reps = pd.factorize(cl.rep)
        codes = pd.Series(codes, index=cl.seqid)
        df[name] = df.seqid.map(codes)
        n_missing = int(df[name].isna().sum())
        assert n_missing == 0, f"{name}: {n_missing} sequences without a cluster"
        df[name] = df[name].astype(np.int64)
        s, t = cluster_stats(df.rename(columns={name: "cluster"}), name)
        t["representative"] = reps[t["representative"].to_numpy()]
        stats.append(s)
        tops.append(t)
        names.append(name)
        print(pd.Series(s).to_string(), flush=True)
    pd.DataFrame(stats).to_csv(os.path.join(out, "cluster_stats.tsv"), sep="\t", index=False, float_format="%.3f")
    pd.concat(tops).to_csv(os.path.join(out, "cluster_largest.tsv"), sep="\t", index=False)

    rows = []
    for how in ["random", *names]:
        rng = np.random.default_rng(SEED)
        tr_i, te_i = split(df, how, rng)
        tr, te = df.iloc[tr_i], df.iloc[te_i]
        # Cluster-weighted scoring counts each test cluster once, under the finest clustering available.
        wcol = names[0] if names else None
        w = 1.0 / te.groupby(wcol)[wcol].transform("size").to_numpy() if wcol else None
        for model, cols in MODELS.items():
            clf = HistGradientBoostingClassifier(max_iter=300, learning_rate=0.1, early_stopping=True,
                                                 validation_fraction=0.1, random_state=SEED)
            clf.fit(tr[cols].to_numpy(), tr.y.to_numpy())
            p = clf.predict_proba(te[cols].to_numpy())[:, 1]
            y = te.y.to_numpy()
            row = {"split": how, "model": model, "n_train": len(tr), "n_test": len(te),
                   "test_pos_frac": float(y.mean()), "roc_auc": auc(y, p),
                   "roc_auc_cluster_weighted": auc(y, p, w)}
            for lo, hi in LENGTH_BINS:
                m = te.aa_length.between(lo, hi).to_numpy()
                row[f"auc_{lo}_{hi}"] = auc(y[m], p[m])
            rows.append(row)
            print(f"{how} | {model}: AUC {row['roc_auc']:.4f}, cluster-weighted {row['roc_auc_cluster_weighted']:.4f}",
                  flush=True)
        pd.DataFrame(rows).to_csv(os.path.join(out, "cluster_split_models.tsv"), sep="\t", index=False,
                                  float_format="%.4f")


if __name__ == "__main__":
    main()
