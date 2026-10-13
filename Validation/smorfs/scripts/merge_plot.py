#!/usr/bin/env python
"""Merge per-chunk entropy tables and plot entropies coloured by dataset.

Usage: merge_plot.py RESULTS_DIR OUT_DIR [N_PER_CLASS]
"""
import glob
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

PALETTE = {"positive": "#2a78d6", "negative": "#eb6834"}
ORDER = ["positive", "negative"]
LABELS = {
    "orf_nt_entropy": "ORF nucleotide entropy (bits)",
    "aa_entropy": "Amino acid entropy (bits)",
    "three_di_entropy": "3Di entropy (bits)",
    "twelve_state_entropy": "12-state entropy (bits)",
}
SEED = 42


def main():
    res, out = sys.argv[1], sys.argv[2]
    n_per_class = int(sys.argv[3]) if len(sys.argv) > 3 else 50_000
    os.makedirs(out, exist_ok=True)

    files = sorted(glob.glob(os.path.join(res, "*.entropy.tsv")))
    print(f"{len(files)} chunk tables", flush=True)
    df = pd.concat((pd.read_csv(f, sep="\t") for f in files), ignore_index=True)
    print(df["dataset"].value_counts().to_string(), flush=True)
    assert df["seqid"].is_unique, "duplicate seqids across chunks"

    table = os.path.join(out, "smorf_entropy.tsv.gz")
    df.to_csv(table + ".tmp", sep="\t", index=False, compression={"method": "gzip", "compresslevel": 6})
    os.rename(table + ".tmp", table)

    cols = ["nt_length", "aa_length", *LABELS]
    summ = df.groupby("dataset")[cols].describe().T
    summ.to_csv(os.path.join(out, "smorf_entropy_summary.tsv"), sep="\t")
    print(summ.to_string(), flush=True)

    # Balanced random subsample; 43.7M points would only overplot into a solid block.
    sub = pd.concat(g.sample(n=min(n_per_class, len(g)), random_state=SEED)
                    for _, g in df.groupby("dataset"))
    sub = sub.sample(frac=1, random_state=SEED)  # shuffle so neither class is drawn on top
    sub.to_csv(os.path.join(out, "smorf_entropy_plot_subsample.tsv"), sep="\t", index=False)

    sns.set_theme(style="ticks", context="paper")
    fig, ax = plt.subplots(figsize=(6.5, 5.5))
    sns.scatterplot(data=sub, x="aa_entropy", y="three_di_entropy", hue="dataset",
                    hue_order=ORDER, palette=PALETTE, s=4, alpha=0.3, linewidth=0, ax=ax,
                    rasterized=True)
    ax.set(xlabel=LABELS["aa_entropy"], ylabel=LABELS["three_di_entropy"],
           title=f"smORF entropy: amino acid vs 3Di\n(random {n_per_class:,} per dataset)")
    sns.despine(ax=ax)
    ax.grid(color="0.9", linewidth=0.5)
    leg = ax.legend(title=None, markerscale=4, frameon=False)
    for h in leg.legend_handles:
        h.set_alpha(1)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(out, f"scatter_aa_vs_3di.{ext}"), dpi=200)
    plt.close(fig)

    g = sns.pairplot(sub, vars=list(LABELS), hue="dataset", hue_order=ORDER, palette=PALETTE,
                     corner=True, plot_kws=dict(s=3, alpha=0.25, linewidth=0, rasterized=True),
                     diag_kind="kde", diag_kws=dict(common_norm=False, fill=False, linewidth=1.5))
    for ax in g.axes.flat:
        if ax is None:
            continue
        ax.set_xlabel(LABELS.get(ax.get_xlabel(), ax.get_xlabel()), fontsize=8)
        ax.set_ylabel(LABELS.get(ax.get_ylabel(), ax.get_ylabel()), fontsize=8)
    g._legend.set_title(None)
    for h in g._legend.legend_handles:
        h.set_alpha(1)
        h.set_markersize(8)
    g.figure.suptitle(f"smORF entropies (random {n_per_class:,} per dataset)", y=1.01)
    for ext in ("png", "pdf"):
        g.savefig(os.path.join(out, f"pairplot_entropies.{ext}"), dpi=200)


if __name__ == "__main__":
    main()
