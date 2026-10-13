#!/usr/bin/env python
"""Length distributions, full-data entropy distributions, single-feature separation by
length, and gradient-boosted models for the smORF positive/negative comparison.

Usage: manuscript_stats.py ENTROPY_TSV_GZ STRIPE_BANDS_TSV OUT_DIR [NPROC]
"""
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.metrics import average_precision_score, roc_auc_score

SEED = 42
N_TRAIN = N_TEST = 2_000_000
ENTROPIES = ["orf_nt_entropy", "aa_entropy", "three_di_entropy", "twelve_state_entropy"]
FEATURES = ["aa_length", *ENTROPIES]
NAMES = {"aa_length": "Length (aa)", "orf_nt_entropy": "Nucleotide entropy",
         "aa_entropy": "Amino acid entropy", "three_di_entropy": "3Di entropy",
         "twelve_state_entropy": "12-state entropy"}
COLOURS = {"aa_length": "#52514e", "orf_nt_entropy": "#eda100", "aa_entropy": "#1baf7a",
           "three_di_entropy": "#2a78d6", "twelve_state_entropy": "#eb6834"}
CLASS_COLOURS = {"positive": "#2a78d6", "negative": "#eb6834"}
LENGTH_BINS = [(10, 29), (30, 49), (50, 59), (60, 69), (70, 79), (80, 89), (90, 99)]
MODELS = {
    "length": ["aa_length"],
    "nt": ["orf_nt_entropy"],
    "aa": ["aa_entropy"],
    "3Di": ["three_di_entropy"],
    "12-state": ["twelve_state_entropy"],
    "length + nt + aa": ["aa_length", "orf_nt_entropy", "aa_entropy"],
    "length + 3Di + 12-state": ["aa_length", "three_di_entropy", "twelve_state_entropy"],
    "length + aa + 3Di": ["aa_length", "aa_entropy", "three_di_entropy"],
    "all five": FEATURES,
}
SURFACE, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e6e5e1"


def style():
    plt.rcParams.update({"font.size": 9, "text.color": INK, "axes.labelcolor": INK2,
                         "xtick.color": INK2, "ytick.color": INK2, "axes.edgecolor": INK2,
                         "axes.spines.top": False, "axes.spines.right": False,
                         "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
                         "axes.axisbelow": True, "figure.facecolor": SURFACE,
                         "axes.facecolor": SURFACE, "legend.frameon": False})


def auc(y, x):
    return roc_auc_score(y, x) if 0 < y.sum() < len(y) else np.nan


def save(fig, out, name):
    for e in ("png", "pdf"):
        fig.savefig(os.path.join(out, f"{name}.{e}"), dpi=200, facecolor=SURFACE)
    plt.close(fig)


def length_figure(df, out):
    fig, ax = plt.subplots(figsize=(6.5, 3.6), layout="constrained")
    bins = np.arange(9.5, 100.5, 1)
    for lab in ("positive", "negative"):
        v = df.loc[df.dataset == lab, "aa_length"]
        ax.hist(v, bins=bins, density=True, histtype="step", linewidth=2,
                color=CLASS_COLOURS[lab], label=f"{lab} (n = {len(v):,})")
    ax.set(xlabel="Protein length (aa, stop codon excluded)", ylabel="Fraction of dataset per aa")
    ax.legend(loc="upper left")
    save(fig, out, "fig1_length_distribution")


def entropy_figure(df, out):
    fig, axes = plt.subplots(2, 2, figsize=(9, 6), layout="constrained")
    for ax, col in zip(axes.flat, ENTROPIES):
        lo, hi = df[col].min(), df[col].max()
        bins = np.linspace(lo, hi, 200)
        for lab in ("positive", "negative"):
            ax.hist(df.loc[df.dataset == lab, col], bins=bins, density=True, histtype="step",
                    linewidth=1.6, color=CLASS_COLOURS[lab], label=lab)
        ax.set(xlabel=f"{NAMES[col]} (bits)", ylabel="Density")
    axes[0, 0].set_xlim(1.5, 2.0)
    axes[0, 1].set_xlim(2.5, 4.4)
    axes[0, 0].legend(loc="upper left")
    save(fig, out, "fig2_entropy_distributions")


def separation(df, out):
    y = (df.dataset == "positive").to_numpy()
    rows = []
    for col in FEATURES:
        rows.append({"length_bin": "all", "feature": col, "n_positive": int(y.sum()),
                     "n_negative": int((~y).sum()), "roc_auc": auc(y, df[col].to_numpy())})
    for lo, hi in LENGTH_BINS:
        m = df.aa_length.between(lo, hi).to_numpy()
        ys = y[m]
        for col in FEATURES:
            rows.append({"length_bin": f"{lo}-{hi}", "feature": col, "n_positive": int(ys.sum()),
                         "n_negative": int((~ys).sum()), "roc_auc": auc(ys, df.loc[m, col].to_numpy())})
    sep = pd.DataFrame(rows)
    sep["separation"] = (sep.roc_auc - 0.5).abs() + 0.5
    sep["direction"] = np.where(sep.roc_auc >= 0.5, "higher in positive", "higher in negative")
    sep.to_csv(os.path.join(out, "table2_single_feature_auc.tsv"), sep="\t", index=False, float_format="%.4f")
    print(sep.to_string(), flush=True)

    within = sep[sep.length_bin != "all"]
    fig, ax = plt.subplots(figsize=(6.5, 4), layout="constrained")
    labels = [f"{lo}–{hi}" for lo, hi in LENGTH_BINS]
    for col in FEATURES:
        v = within[within.feature == col].set_index("length_bin").reindex([f"{lo}-{hi}" for lo, hi in LENGTH_BINS])
        ax.plot(labels, v.roc_auc, marker="o", markersize=5, linewidth=2, color=COLOURS[col], label=NAMES[col])
    ax.axhline(0.5, color=INK2, linewidth=0.8, linestyle="--")
    ax.set(xlabel="Protein length bin (aa)", ylabel="ROC AUC (positive vs negative)", ylim=(0.3, 0.9))
    ax.legend(loc="upper left", ncols=2, fontsize=8)
    save(fig, out, "fig4_auc_by_length")


def models(df, out, nproc):
    rng = np.random.default_rng(SEED)
    n = min(N_TRAIN, len(df) // 2)
    idx = rng.permutation(len(df))[:2 * n]
    tr, te = df.iloc[idx[:n]], df.iloc[idx[n:]]
    ytr, yte = (tr.dataset == "positive").to_numpy(), (te.dataset == "positive").to_numpy()
    rows, preds = [], {}
    for name, cols in MODELS.items():
        clf = HistGradientBoostingClassifier(max_iter=300, learning_rate=0.1, early_stopping=True,
                                             validation_fraction=0.1, random_state=SEED)
        clf.fit(tr[cols].to_numpy(), ytr)
        p = clf.predict_proba(te[cols].to_numpy())[:, 1]
        preds[name] = p
        row = {"model": name, "features": "+".join(cols), "roc_auc": roc_auc_score(yte, p),
               "ap_negative": average_precision_score(~yte, 1 - p),
               "negative_prevalence": float((~yte).mean())}
        for lo, hi in LENGTH_BINS:
            m = te.aa_length.between(lo, hi).to_numpy()
            row[f"auc_{lo}_{hi}"] = auc(yte[m], p[m])
        rows.append(row)
        print(f"model {name}: AUC {row['roc_auc']:.4f}", flush=True)
    res = pd.DataFrame(rows)
    res.to_csv(os.path.join(out, "table3_models.tsv"), sep="\t", index=False, float_format="%.4f")
    models_figure(res, out)


def models_figure(res, out):
    fig, ax = plt.subplots(figsize=(6.5, 4), layout="constrained")
    labels = [f"{lo}–{hi}" for lo, hi in LENGTH_BINS]
    show = {"length + nt + aa": "#1baf7a", "length + 3Di + 12-state": "#2a78d6", "all five": "#0d366b"}
    for name, colour in show.items():
        r = res.set_index("model").loc[name]
        ax.plot(labels, [r[f"auc_{lo}_{hi}"] for lo, hi in LENGTH_BINS], marker="o", markersize=5,
                linewidth=2, color=colour, label=f"{name} (overall {r.roc_auc:.3f})")
    ax.axhline(0.5, color=INK2, linewidth=0.8, linestyle="--")
    ax.set(xlabel="Protein length bin (aa)", ylabel="ROC AUC within bin (held-out)", ylim=(0.45, 1.0))
    ax.legend(loc="upper right", fontsize=8)
    save(fig, out, "fig5_models_within_length")


def stripe_figure(bands_path, out):
    b = pd.read_csv(bands_path, sep="\t")
    b = b[b.k.isin([2, 3])]
    fig, axes = plt.subplots(1, 2, figsize=(8, 3.4), sharey=True, layout="constrained")
    order = ["below", "stripe", "above"]
    width = 0.38
    for ax, k in zip(axes, (2, 3)):
        x = np.arange(len(order))
        for i, lab in enumerate(("positive", "negative")):
            v = b[(b.k == k) & (b.dataset == lab)].set_index("band").reindex(order)
            ax.bar(x + (i - 0.5) * width, v.pct_with_k_states, width - 0.04, color=CLASS_COLOURS[lab],
                   label=lab)
        centre = np.log2(k)
        ax.set_xticks(x, [f"{centre - 0.03:.3f}", f"{centre:.3f}\n(stripe)", f"{centre + 0.03:.3f}"])
        ax.set(title=f"3Di entropy band near log2({k})", xlabel="Band centre (bits)")
        ax.grid(axis="x", visible=False)
    axes[0].set_ylabel(f"% of strings using exactly k states")
    axes[0].legend(loc="upper right")
    save(fig, out, "fig6_stripe_state_counts")


def main():
    src, bands, out = sys.argv[1:4]
    nproc = int(sys.argv[4]) if len(sys.argv) > 4 else os.cpu_count()
    os.makedirs(out, exist_ok=True)
    style()
    df = pd.read_csv(src, sep="\t", usecols=["dataset", *FEATURES],
                     dtype={"dataset": "category", "aa_length": "int16",
                            **{c: "float32" for c in ENTROPIES}})
    print(df.dataset.value_counts().to_string(), flush=True)
    counts = (df.groupby(["dataset", pd.cut(df.aa_length, [b[0] - 1 for b in LENGTH_BINS] + [99])],
                         observed=False).size().unstack(0))
    counts.to_csv(os.path.join(out, "table1b_length_bins.tsv"), sep="\t")
    length_figure(df, out)
    entropy_figure(df, out)
    stripe_figure(bands, out)
    separation(df, out)
    models(df, out, nproc)


if __name__ == "__main__":
    main()
