#!/usr/bin/env python3
"""Figures for report section 4.2: estimator bias against biology, and calibration.

bias_vs_biology.png
    Rise in mean entropy with length, relative to L = 100, for ORFs matched to
    a deposited CDS (stage 33 empirical reference, stratum 3) and for i.i.d.
    sequences of the same length and pooled composition (stage 31 random
    reference). The random curve is the plug-in estimator's bias and nothing
    else; the gap between the two is the part of the trend that is biology.

zscore_calibration.png
    z and decile occupancy for the figure sample's matched ORFs, scored
    against each reference by 34_calibration_scores.pbs. A calibrated
    reference gives z ~ N(0, 1) and 10% per decile.

The decomposition table behind the report's "how much is an artefact" table is
written alongside as bias_decomposition.tsv, so those numbers have a script.

  35_bias_vs_biology.py --out-dir <figures dir> [--table-dir <dir>]
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker

GD = "/g/data/ob80/re3494/gtdb_entropy"
NULLDIR = f"{GD}/entropy_null"

INK = "#2b2a28"
MUTED = "#52514e"
PLATE = "#fcfcfb"
OBSERVED = "#1f6f8b"
BIAS = "#c46b4a"
BIOLOGY = "#9fc4d2"
EMPIRICAL = "#1f6f8b"
RANDOM = "#c46b4a"

ALPHABETS = [("protein", "protein"), ("three_di", "3Di"), ("twelve_state", "12-state")]
DOMAINS = [("bac", "Bacteria"), ("arc", "Archaea")]
MATCHED = 3
REF_LENGTH = 100
REPORT_LENGTHS = (150, 200, 300, 500, 800)


def style():
    plt.rcParams.update({
        "figure.dpi": 160, "savefig.dpi": 160,
        "font.size": 8.5, "axes.labelsize": 9, "axes.titlesize": 9.5,
        "axes.edgecolor": MUTED, "axes.labelcolor": INK,
        "text.color": INK, "xtick.color": MUTED, "ytick.color": MUTED,
        "axes.spines.top": False, "axes.spines.right": False,
        "figure.facecolor": PLATE, "savefig.facecolor": PLATE,
        "axes.grid": True, "grid.color": "#e6e3de", "grid.linewidth": 0.6,
        "legend.frameon": False, "legend.fontsize": 7.6,
    })


def save(fig, out_dir, name):
    p = Path(out_dir) / name
    fig.savefig(p, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {p}", file=sys.stderr)


def mean_curves():
    """{(domain, alphabet): DataFrame(aa_length, observed, bias)} over 1 aa bins."""
    null = pd.read_csv(f"{NULLDIR}/entropy_null_quantiles.tsv", sep="\t",
                       usecols=["domain", "alphabet", "aa_length", "mean"])
    curves = {}
    for dom, _ in DOMAINS:
        emp = pd.read_csv(f"{NULLDIR}/empirical/length_entropy_{dom}.tsv", sep="\t",
                          usecols=["stratum", "alphabet", "length_lo", "length_hi",
                                   "n", "mean"])
        # Only the 1 aa bins: above 999 they are 10 and 100 aa wide and a bin
        # mean is not a point on the curve. 99.5% of ORFs are below 1000 aa.
        emp = emp[(emp.stratum == MATCHED) & (emp.length_lo == emp.length_hi)]
        for alphabet, _ in ALPHABETS:
            e = emp[emp.alphabet == alphabet].sort_values("length_lo")
            r = null[(null.domain == dom) & (null.alphabet == alphabet)] \
                .sort_values("aa_length")
            if e.empty or r.empty:
                sys.exit(f"no curve for {dom}/{alphabet}")
            lengths = e.length_lo.to_numpy()
            # The random grid is every length to 600 and every 5th to 1000; it
            # is a smooth curve, so interpolation onto the empirical lengths is
            # correct here (32_entropy_percentile.py does the same).
            bias = np.interp(lengths, r.aa_length.to_numpy(), r["mean"].to_numpy())
            curves[(dom, alphabet)] = pd.DataFrame({
                "aa_length": lengths, "n": e.n.to_numpy(),
                "observed": e["mean"].to_numpy(), "bias": bias})
    return curves


def decomposition(curves):
    rows = []
    for (dom, alphabet), c in curves.items():
        base = c[c.aa_length == REF_LENGTH].iloc[0]
        for length in REPORT_LENGTHS:
            row = c[c.aa_length == length].iloc[0]
            observed = row.observed - base.observed
            bias = row.bias - base.bias
            rows.append({
                "domain": dom, "alphabet": alphabet, "from_length": REF_LENGTH,
                "to_length": length, "observed_rise": round(observed, 4),
                "bias_rise": round(bias, 4),
                "biology_rise": round(observed - bias, 4),
                "artefact_pct": round(100 * bias / observed, 1) if observed else np.nan})
    return pd.DataFrame(rows)


def plot_bias(curves, out_dir):
    fig, axes = plt.subplots(2, 3, figsize=(10.5, 5.6), sharex=True)
    for i, (dom, dom_label) in enumerate(DOMAINS):
        for j, (alphabet, label) in enumerate(ALPHABETS):
            ax = axes[i, j]
            c = curves[(dom, alphabet)]
            base = c[c.aa_length == REF_LENGTH].iloc[0]
            obs = c.observed - base.observed
            bias = c.bias - base.bias
            ax.fill_between(c.aa_length, bias, obs, where=obs >= bias,
                            color=BIOLOGY, alpha=0.55, lw=0, label="biology")
            ax.fill_between(c.aa_length, obs, bias, where=obs < bias,
                            color=BIAS, alpha=0.18, lw=0,
                            label="bias exceeds observed")
            ax.plot(c.aa_length, obs, color=OBSERVED, lw=1.4,
                    label="observed (matched ORFs)")
            ax.plot(c.aa_length, bias, color=BIAS, lw=1.4, ls="--",
                    label="estimator bias (i.i.d. reference)")
            ax.axhline(0, color=MUTED, lw=0.6)
            ax.axvline(REF_LENGTH, color=MUTED, lw=0.6, ls=":")
            ax.set_xscale("log")
            ax.set_xlim(90, 999)
            ax.set_xticks([100, 200, 300, 500, 999])
            ax.set_xticklabels(["100", "200", "300", "500", "999"])
            ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
            at300 = c[c.aa_length == 300].iloc[0]
            o300 = at300.observed - base.observed
            b300 = at300.bias - base.bias
            ax.text(0.97, 0.05,
                    f"100→300 aa: observed {o300:+.2f},\nbias {b300:+.2f} bits",
                    transform=ax.transAxes, ha="right", va="bottom", fontsize=7.2,
                    color=MUTED, bbox=dict(facecolor=PLATE, edgecolor="#e0ddd8",
                                           boxstyle="round,pad=0.3"))
            if i == 0:
                ax.set_title(label)
            if i == 1:
                ax.set_xlabel("ORF length (aa)")
            if j == 0:
                ax.set_ylabel(f"{dom_label}\nΔ mean entropy vs L = 100 (bits)")
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=4,
               bbox_to_anchor=(0.5, -0.02))
    fig.suptitle("Entropy–length trend: plug-in estimator bias vs biology",
                 y=1.0, fontsize=10.5)
    fig.tight_layout(rect=(0, 0.05, 1, 1))
    save(fig, out_dir, "bias_vs_biology.png")


def load_scores(calib_dir, dom, alphabet):
    p = Path(calib_dir) / f"scores_{dom}_{alphabet}.tsv.gz"
    if not p.exists():
        sys.exit(f"{p} missing: run 34_calibration_scores.pbs first")
    return pd.read_csv(p, sep="\t", usecols=["z", "percentile", "emp_z",
                                              "emp_percentile"])


def plot_calibration(calib_dir, out_dir):
    fig, axes = plt.subplots(2, 3, figsize=(10.5, 5.6))
    grid = np.linspace(-6, 6, 400)
    normal = np.exp(-grid ** 2 / 2) / np.sqrt(2 * np.pi)
    summary = []
    for i, (dom, dom_label) in enumerate(DOMAINS):
        deciles = {}
        for j, (alphabet, label) in enumerate(ALPHABETS[:2]):
            s = load_scores(calib_dir, dom, alphabet)
            ax = axes[i, j]
            bins = np.linspace(-15, 8, 231)
            # Not clipped: clipping piles the i.i.d. tail into the edge bins.
            # Density stays relative to all rows, so out-of-range mass is
            # simply absent, and its share is stated in the legend.
            outside = ((s.z < bins[0]) | (s.z > bins[-1])).mean()
            weights = np.full(len(s), 1 / (len(s) * (bins[1] - bins[0])))
            ax.hist(s.z, bins=bins, weights=weights, color=RANDOM, alpha=0.55,
                    label=f"i.i.d.: mean {s.z.mean():+.2f}, sd {s.z.std():.2f} "
                          f"({100 * outside:.1f}% off axis)")
            ax.hist(s.emp_z, bins=bins, weights=weights,
                    color=EMPIRICAL, alpha=0.6,
                    label=f"empirical: mean {s.emp_z.mean():+.3f}, "
                          f"sd {s.emp_z.std():.3f}")
            ax.plot(grid, normal, color=INK, lw=1.0, ls="--", label="N(0, 1)")
            ax.set_xlim(bins[0], bins[-1])
            ax.legend(loc="upper left", fontsize=6.9)
            if i == 0:
                ax.set_title(f"{label} z")
            if i == 1:
                ax.set_xlabel("z")
            if j == 0:
                ax.set_ylabel(f"{dom_label}\ndensity  (n = {len(s):,})")
            for ref, col in (("empirical", "emp_percentile"), ("i.i.d.", "percentile")):
                occ = np.histogram(s[col].clip(0, 1), bins=np.linspace(0, 1, 11))[0]
                deciles[f"{ref} {label}"] = 100 * occ / occ.sum()
            summary.append({"domain": dom, "alphabet": alphabet, "n": len(s),
                            "iid_z_mean": s.z.mean(), "iid_z_sd": s.z.std(),
                            "emp_z_mean": s.emp_z.mean(), "emp_z_sd": s.emp_z.std()})
        ax = axes[i, 2]
        x = np.arange(10)
        colours = [EMPIRICAL, "#5aa3bd", RANDOM, "#e0a58c"]
        keys = ["empirical protein", "empirical 3Di", "i.i.d. protein", "i.i.d. 3Di"]
        w = 0.2
        for k, key in enumerate(keys):
            ax.bar(x + (k - 1.5) * w, deciles[key], width=w, color=colours[k],
                   label=key)
        ax.axhline(10, color=INK, lw=0.9, ls="--")
        ax.set_yscale("log")
        ax.set_ylim(0.3, 100)
        ax.set_yticks([0.5, 1, 2, 5, 10, 20, 50, 100])
        ax.set_yticklabels(["0.5", "1", "2", "5", "10", "20", "50", "100"])
        ax.set_xticks(x)
        ax.set_xticklabels([f"{10 * (k + 1)}" for k in x])
        if i == 0:
            ax.set_title("decile occupancy (%; 10% = calibrated)")
            ax.legend(loc="upper center", ncol=2, fontsize=6.9)
        if i == 1:
            ax.set_xlabel("percentile decile")
    fig.suptitle("Calibration of the two length-conditioned references "
                 "(matched ORFs, systematic sample)", y=1.0, fontsize=10.5)
    fig.tight_layout()
    save(fig, out_dir, "zscore_calibration.png")
    return pd.DataFrame(summary)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--table-dir", default=NULLDIR)
    ap.add_argument("--calibration-dir", default=f"{NULLDIR}/calibration")
    args = ap.parse_args()
    Path(args.out_dir).mkdir(parents=True, exist_ok=True)
    style()

    curves = mean_curves()
    table = decomposition(curves)
    table.to_csv(Path(args.table_dir) / "bias_decomposition.tsv", sep="\t", index=False)
    print(table.to_string(index=False))
    plot_bias(curves, args.out_dir)

    summary = plot_calibration(args.calibration_dir, args.out_dir)
    summary.to_csv(Path(args.table_dir) / "calibration_summary.tsv", sep="\t",
                   index=False, float_format="%.4f")
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
