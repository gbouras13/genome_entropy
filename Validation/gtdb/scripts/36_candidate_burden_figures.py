#!/usr/bin/env python3
"""Regenerate candidate_burden_{bac,arc}.png from candidate_burden_<dom>.tsv.

Replaces the figure drawn by 16_candidate_burden.py, which (i) plotted only
per-Mb burden while §6's table quotes raw-count rho, (ii) had no N50 panel
although N50 is the covariate §6 discusses, and (iii) rendered degenerately for
bacteria: the median genome carries 3 candidates, so a hexbin of a small
integer count shows horizontal stripes, and the contig_count panel collapsed to
two bands.

This version draws, per covariate, binned quantile summaries rather than a
point cloud: the covariate is cut into (up to) 20 equal-count bins and each bin
shows the median and the 25-75% and 10-90% ranges of burden. Two rows: raw
candidate count per genome (the quantity the §6 table's rho refers to) and
candidates per Mb. Both Spearman rhos are printed on every panel and written
to candidate_burden_rho.tsv so the text can quote the same numbers.

  36_candidate_burden_figures.py --burden-dir work/gtdb_extras --out-dir figures
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

INK, MUTED, PLATE = "#2b2a28", "#52514e", "#fcfcfb"
BLUE, BLUE_L, BLUE_LL = "#184f95", "#6da7ec", "#cde2fb"
COVS = [("n50_contigs", "contig N50 (bp)", True),
        ("contig_count", "contigs", True),
        ("genome_size_mb", "genome size (Mb)", False),
        ("coding_density", "coding density (%)", False),
        ("frac_orfs_in_genbank", "fraction of ORFs matched to a CDS", False),
        ("checkm2_completeness", "CheckM2 completeness (%)", False),
        ("checkm2_contamination", "CheckM2 contamination (%)", False)]
DOM = {"bac": "bacterial", "arc": "archaeal"}


def rho(x, y):
    ok = x.notna() & y.notna()
    return float(x[ok].rank().corr(y[ok].rank())), int(ok.sum())


def binned(x, y, nbins=20):
    e = np.unique(np.nanquantile(x, np.linspace(0, 1, nbins + 1)))
    b = pd.cut(x, e, include_lowest=True)
    g = pd.DataFrame({"x": x, "y": y, "b": b}).groupby("b", observed=True)
    return pd.DataFrame({"x": g.x.median(), "q10": g.y.quantile(.1),
                         "q25": g.y.quantile(.25), "q50": g.y.median(),
                         "q75": g.y.quantile(.75), "q90": g.y.quantile(.9),
                         "mean": g.y.mean(), "n": g.y.size()})


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--burden-dir", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--rho-tsv", default=None)
    a = ap.parse_args()
    sns.set_theme(style="whitegrid", context="paper")
    plt.rcParams.update({"figure.facecolor": PLATE, "savefig.facecolor": PLATE,
                         "axes.facecolor": PLATE, "axes.spines.top": False,
                         "axes.spines.right": False, "grid.color": "#e6e3de",
                         "font.size": 8.5, "text.color": INK,
                         "axes.labelcolor": INK, "savefig.dpi": 160})
    rows = []
    for dom in ("bac", "arc"):
        df = pd.read_csv(Path(a.burden_dir) / f"candidate_burden_{dom}.tsv", sep="\t")
        fig, axes = plt.subplots(2, len(COVS), figsize=(3.0 * len(COVS), 6.4))
        for c, (cov, lab, logx) in enumerate(COVS):
            r_raw, n = rho(df[cov], df.n_candidates)
            r_mb, _ = rho(df[cov], df.cand_per_mb)
            rows.append(dict(domain=dom, covariate=cov, rho_per_genome=r_raw,
                             rho_per_mb=r_mb, n_genomes=n))
            for r, (ycol, ylab, rv) in enumerate(
                    (("n_candidates", "candidates per genome", r_raw),
                     ("cand_per_mb", "candidates per Mb", r_mb))):
                ax = axes[r, c]
                bt = binned(df[cov], df[ycol])
                ax.fill_between(bt.x, bt.q10, bt.q90, color=BLUE_LL, lw=0,
                                label="10–90%", step=None)
                ax.fill_between(bt.x, bt.q25, bt.q75, color=BLUE_L, lw=0,
                                alpha=0.8, label="25–75%")
                ax.plot(bt.x, bt.q50, "-o", color=BLUE, lw=2, ms=3.5, label="median")
                ax.plot(bt.x, bt["mean"], "--", color=INK, lw=1, label="mean")
                if logx:
                    ax.set_xscale("log")
                ax.set_yscale("symlog", linthresh=1 if r == 0 else 0.1)
                ax.set_ylim(0, None)
                ax.set_xlabel(lab)
                if c == 0:
                    ax.set_ylabel(ylab)
                ax.set_title(f"{cov}\nSpearman ρ = {rv:+.3f}", loc="left",
                             fontsize=8.5, fontweight="bold" if cov == "n50_contigs" else None)
        axes[0, 0].legend(fontsize=7, loc="upper left")
        nG = len(df)
        fig.suptitle(
            f"Candidate burden against genome covariates, {nG:,} annotated {DOM[dom]} "
            f"representatives (median {df.n_candidates.median():.0f} candidates per genome, "
            f"{df.n_candidates.sum():,} in total)", x=0.01, ha="left", fontsize=11,
            fontweight="bold")
        fig.text(0.01, 0.935,
                 "Each panel: covariate cut into ≤20 equal-count bins; line = median burden "
                 "per bin, bands = 25–75% and 10–90%, dashed = mean. Top row: raw count "
                 "(the ρ quoted in §6); bottom row: per Mb. y axis symlog. "
                 "ρ computed on all genomes, not on bins.", fontsize=8, color=MUTED,
                 ha="left")
        fig.tight_layout(rect=[0, 0, 1, 0.93])
        p = Path(a.out_dir) / f"candidate_burden_{dom}.png"
        fig.savefig(p, bbox_inches="tight")
        plt.close(fig)
        print(f"  {p.name} ({p.stat().st_size/1024:.0f} kB)")
    t = pd.DataFrame(rows)
    print(t.to_string(index=False, float_format=lambda x: f"{x:+.3f}"))
    if a.rho_tsv:
        t.to_csv(a.rho_tsv, sep="\t", index=False)


if __name__ == "__main__":
    main()
