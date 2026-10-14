#!/usr/bin/env python
"""Binned density of amino acid vs 3Di entropy over every smORF, per dataset.

Usage: density_aa_vs_3di.py ENTROPY_TSV_GZ OUT_DIR [N_BINS]
"""
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, LogNorm, TwoSlopeNorm

ORDER = ["positive", "negative"]
BLUE = ["#cde2fb", "#86b6ef", "#3987e5", "#256abf", "#184f95", "#0d366b"]
ORANGE = ["#fde3d6", "#f6b494", "#f08a5d", "#eb6834", "#c4501f", "#8f3612"]
CMAPS = {
    "positive": LinearSegmentedColormap.from_list("blue", BLUE),
    "negative": LinearSegmentedColormap.from_list("orange", ORANGE),
}
DIVERGING = LinearSegmentedColormap.from_list("red_gray_blue", ["#b52f2e", "#e34948", "#f0efec", "#3987e5", "#184f95"])
SURFACE, INK, INK2 = "#fcfcfb", "#0b0b0b", "#52514e"
XLAB, YLAB = "Amino acid entropy (bits)", "3Di entropy (bits)"


def main():
    src, out = sys.argv[1], sys.argv[2]
    n_bins = int(sys.argv[3]) if len(sys.argv) > 3 else 400
    os.makedirs(out, exist_ok=True)

    df = pd.read_csv(src, sep="\t", usecols=["dataset", "aa_entropy", "three_di_entropy"],
                     dtype={"dataset": "category", "aa_entropy": "float32", "three_di_entropy": "float32"})
    print(df["dataset"].value_counts().to_string(), flush=True)
    n_nan = int(df[["aa_entropy", "three_di_entropy"]].isna().any(axis=1).sum())
    print(f"rows with missing entropy: {n_nan}", flush=True)

    x_edges = np.linspace(0, 4.4, n_bins + 1)
    y_edges = np.linspace(0, 4.4, n_bins + 1)
    counts, totals = {}, {}
    for lab in ORDER:
        g = df[df["dataset"] == lab]
        h, _, _ = np.histogram2d(g["aa_entropy"], g["three_di_entropy"], bins=[x_edges, y_edges])
        counts[lab], totals[lab] = h.T, len(g)  # transpose: rows = y
        assert h.sum() == len(g), f"{lab}: {len(g) - h.sum():.0f} points outside bin range"
    np.savez_compressed(os.path.join(out, "density_aa_vs_3di_counts.npz"),
                        x_edges=x_edges, y_edges=y_edges, **counts)

    frac = {lab: counts[lab] / totals[lab] for lab in ORDER}
    vmin = min(f[f > 0].min() for f in frac.values())
    vmax = max(f.max() for f in frac.values())
    diff = np.ma.masked_where((counts["positive"] + counts["negative"]) == 0,
                              frac["positive"] - frac["negative"])
    dlim = np.abs(diff).max()

    plt.rcParams.update({"font.size": 9, "text.color": INK, "axes.labelcolor": INK2,
                         "xtick.color": INK2, "ytick.color": INK2, "axes.edgecolor": INK2})
    fig, axes = plt.subplots(1, 3, figsize=(15, 5), sharex=True, sharey=True, facecolor=SURFACE,
                             layout="constrained")
    ext = [x_edges[0], x_edges[-1], y_edges[0], y_edges[-1]]
    for ax, lab in zip(axes, ORDER):
        m = np.ma.masked_equal(frac[lab], 0)
        im = ax.imshow(m, origin="lower", extent=ext, aspect="auto", cmap=CMAPS[lab],
                       norm=LogNorm(vmin, vmax), interpolation="nearest")
        ax.set_title(f"{lab.capitalize()} (n = {totals[lab]:,})", loc="left")
        fig.colorbar(im, ax=ax, shrink=0.85, label="Fraction of dataset per bin")
    im = axes[2].imshow(diff, origin="lower", extent=ext, aspect="auto", cmap=DIVERGING,
                        norm=TwoSlopeNorm(0, -dlim, dlim), interpolation="nearest")
    axes[2].set_title("Positive − negative (fraction per bin)", loc="left")
    cb = fig.colorbar(im, ax=axes[2], shrink=0.85)
    cb.set_label("← enriched in negative     enriched in positive →")
    for ax in axes:
        ax.set_facecolor(SURFACE)
        ax.set_xlabel(XLAB)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
    axes[0].set_ylabel(YLAB)
    fig.suptitle(f"smORF amino acid vs 3Di entropy, all {sum(totals.values()):,} sequences "
                 f"({n_bins}×{n_bins} bins, each dataset normalised to 1)")
    for e in ("png", "pdf"):
        fig.savefig(os.path.join(out, f"density_aa_vs_3di.{e}"), dpi=200, facecolor=SURFACE)


if __name__ == "__main__":
    main()
