#!/usr/bin/env python3
"""Phage (Inphared + Pharokka/phold) mutual-information and length-control analysis.

Input: the per-ORF Parquet table written by build_phage_orf_table.py.

Computes, for the label ``in_genbank`` (ORF matches a phold/Pharokka CDS):

1. binned plug-in mutual information I(X; label) for each single variable
   (DNA, protein, 3Di, 12-state entropy, per-ORF 3Di-12-state MI, ORF length)
   and I((A, B); label) for every pair, on (a) a class-balanced sample and
   (b) the full population, at 10/20/40 empirical-quantile bins per axis;
2. a label-permutation null (balanced sample) and bootstrap intervals;
3. length control: per-length-bin MI, AUC and AP, and the conditional MI
   I(X; label | length) for a 7-bin and a 20-quantile length stratification;
4. single-feature ROC AUC and average precision (overall and per length bin);
5. genome-grouped cross-validated gradient-boosting models (length alone vs
   length + entropies) to express the gain over length as AUC/AP;
6. figures (PNG + PDF).

All randomness is seeded (``--seed``).
"""

import argparse
import json
import math
from itertools import combinations
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.colors import LogNorm, TwoSlopeNorm
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.model_selection import GroupKFold

VARIABLES = (
    "dna_entropy",
    "protein_entropy",
    "three_di_entropy",
    "twelve_state_entropy",
    "three_di_twelve_state_mutual_information",
    "protein_length",
)
SHORT = {
    "dna_entropy": "DNA H",
    "protein_entropy": "protein H",
    "three_di_entropy": "3Di H",
    "twelve_state_entropy": "12-state H",
    "three_di_twelve_state_mutual_information": "3Di–12-state MI",
    "protein_length": "length",
}
ENTROPY_VARIABLES = VARIABLES[:5]
LENGTH_BINS = (30, 50, 75, 90, 100, 150, 200, 300, np.inf)
# 30-49, 50-74 and 75-89 are kept to document that they are empty: the
# effective ORF floor is 90 aa (issue #103).
LENGTH_LABELS = ("30-49", "50-74", "75-89", "90-99", "100-149", "150-199", "200-299",
                 "300+")
POS_COLOUR = "#0072B2"   # in_genbank (phold CDS)
NEG_COLOUR = "#D55E00"   # not annotated


# ----------------------------------------------------------------------------
# Mutual-information helpers (plug-in, bits)
# ----------------------------------------------------------------------------

def quantile_codes(values: np.ndarray, bins: int) -> np.ndarray:
    """Integer codes from empirical-quantile bins (ties merge bins)."""
    edges = np.unique(np.quantile(values, np.arange(1, bins) / bins))
    return np.searchsorted(edges, values, side="right").astype(np.int64)


def mi_from_codes(codes: np.ndarray, labels: np.ndarray) -> float:
    """Plug-in I(codes; labels) in bits for binary labels."""
    n = codes.size
    joint = np.bincount(codes * 2 + labels, minlength=2 * (codes.max() + 1))
    joint = joint.reshape(-1, 2).astype(float)
    px = joint.sum(1, keepdims=True)
    py = joint.sum(0, keepdims=True)
    nz = joint > 0
    return float((joint[nz] / n * np.log2(joint[nz] * n / (px @ py)[nz])).sum())


def label_entropy(labels: np.ndarray) -> float:
    p = labels.mean()
    return float(-(p * math.log2(p) + (1 - p) * math.log2(1 - p)))


def pair_codes(codes_a: np.ndarray, codes_b: np.ndarray) -> np.ndarray:
    return codes_a * (codes_b.max() + 1) + codes_b


def feature_sets() -> list[tuple[str, ...]]:
    return [(v,) for v in VARIABLES] + list(combinations(VARIABLES, 2))


def codes_for(frame: pd.DataFrame, features: tuple[str, ...], bins: int,
              cache: dict) -> np.ndarray:
    for name in features:
        if (name, bins) not in cache:
            cache[(name, bins)] = quantile_codes(frame[name].to_numpy(), bins)
    if len(features) == 1:
        return cache[(features[0], bins)]
    return pair_codes(cache[(features[0], bins)], cache[(features[1], bins)])


def conditional_mi(x: np.ndarray, labels: np.ndarray, strata: np.ndarray,
                   bins: int) -> tuple[float, list[float]]:
    """I(X; label | stratum) with X re-quantiled inside each stratum."""
    total = 0.0
    per_stratum = []
    n = labels.size
    for stratum in np.unique(strata):
        mask = strata == stratum
        y = labels[mask]
        if y.min() == y.max():
            per_stratum.append(0.0)
            continue
        value = mi_from_codes(quantile_codes(x[mask], bins), y)
        per_stratum.append(value)
        total += mask.sum() / n * value
    return total, per_stratum


# ----------------------------------------------------------------------------
# Analyses
# ----------------------------------------------------------------------------

def balanced_sample(frame: pd.DataFrame, per_class: int, rng) -> pd.DataFrame:
    parts = []
    for label in (False, True):
        index = np.flatnonzero(frame.in_genbank.to_numpy() == label)
        take = rng.choice(index, size=min(per_class, index.size), replace=False)
        parts.append(frame.iloc[np.sort(take)])
    return pd.concat(parts, ignore_index=True)


def mi_tables(frame: pd.DataFrame, sample: pd.DataFrame, bins_list, permutations,
              bootstraps, rng) -> pd.DataFrame:
    rows = []
    y_full = frame.in_genbank.to_numpy().astype(np.int64)
    y_samp = sample.in_genbank.to_numpy().astype(np.int64)
    h_full = label_entropy(y_full)
    h_samp = label_entropy(y_samp)
    perm_labels = [rng.permutation(y_samp) for _ in range(permutations)]
    boot_index = [rng.integers(0, y_samp.size, y_samp.size) for _ in range(bootstraps)]
    for bins in bins_list:
        cache_full, cache_samp = {}, {}
        for features in feature_sets():
            c_full = codes_for(frame, features, bins, cache_full)
            c_samp = codes_for(sample, features, bins, cache_samp)
            mi_full = mi_from_codes(c_full, y_full)
            mi_samp = mi_from_codes(c_samp, y_samp)
            row = {
                "features": "+".join(features),
                "kind": "single" if len(features) == 1 else "pair",
                "bins_per_axis": bins,
                "occupied_cells_full": int(np.unique(c_full).size),
                "mi_full_bits": mi_full,
                "mi_full_fraction_of_H": mi_full / h_full,
                "mi_balanced_bits": mi_samp,
                "mi_balanced_fraction_of_H": mi_samp / h_samp,
            }
            if permutations:
                null = np.array([mi_from_codes(c_samp, yp) for yp in perm_labels])
                row.update({
                    "null_mean_bits": null.mean(),
                    "null_p99_bits": float(np.quantile(null, 0.99)),
                    "null_max_bits": null.max(),
                    "permutation_p": (1 + (null >= mi_samp).sum()) / (1 + null.size),
                    "mi_balanced_minus_null_mean": mi_samp - null.mean(),
                })
            if bins == 20 and bootstraps:
                boot = np.array([mi_from_codes(c_samp[i], y_samp[i]) for i in boot_index])
                # basic (reverse-percentile) bootstrap interval: the plug-in
                # MI of a bootstrap resample is biased upward by duplicates.
                row["boot_lo_bits"] = 2 * mi_samp - float(np.quantile(boot, 0.975))
                row["boot_hi_bits"] = 2 * mi_samp - float(np.quantile(boot, 0.025))
                if len(features) == 2:
                    best = []
                    for i in boot_index:
                        best.append(max(
                            mi_from_codes(codes_for(sample, (f,), bins, cache_samp)[i], y_samp[i])
                            for f in features))
                    gain = boot - np.array(best)
                    point = mi_samp - max(
                        mi_from_codes(codes_for(sample, (f,), bins, cache_samp), y_samp)
                        for f in features)
                    row["gain_boot_lo_bits"] = 2 * point - float(np.quantile(gain, 0.975))
                    row["gain_boot_hi_bits"] = 2 * point - float(np.quantile(gain, 0.025))
            rows.append(row)
    table = pd.DataFrame(rows)
    # gain of each pair over its best single variable (same bins, same data)
    single = table[table.kind == "single"].set_index(["features", "bins_per_axis"])
    for column in ("mi_full_bits", "mi_balanced_bits"):
        gains, bests = [], []
        for _, row in table.iterrows():
            if row.kind == "single":
                gains.append(np.nan)
                bests.append("")
                continue
            parts = row.features.split("+")
            values = {p: single.loc[(p, row.bins_per_axis), column] for p in parts}
            best = max(values, key=values.get)
            bests.append(best)
            gains.append(row[column] - values[best])
        table[column.replace("mi_", "gain_over_best_single_")] = gains
        table["best_single_" + column.split("_")[1]] = bests
    table.attrs["h_full"] = h_full
    table.attrs["h_balanced"] = h_samp
    return table


def auc_ap(scores: np.ndarray, labels: np.ndarray) -> tuple[float, float]:
    if labels.min() == labels.max():
        return np.nan, np.nan
    return roc_auc_score(labels, scores), average_precision_score(labels, scores)


def length_tables(frame: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    y = frame.in_genbank.to_numpy().astype(np.int64)
    length = frame.protein_length.to_numpy()
    bin_index = np.digitize(length, LENGTH_BINS[1:-1], right=False)
    rows = []
    # overall
    for name in VARIABLES:
        auc, ap = auc_ap(frame[name].to_numpy(), y)
        rows.append({"length_bin": "all", "variable": name, "n": y.size,
                     "n_pos": int(y.sum()), "prevalence": y.mean(), "auc": auc,
                     "ap": ap, "mi_bits_20bins": mi_from_codes(
                         quantile_codes(frame[name].to_numpy(), 20), y)})
    for b, label in enumerate(LENGTH_LABELS):
        mask = bin_index == b
        yb = y[mask]
        for name in VARIABLES:
            if mask.sum() == 0 or yb.min() == yb.max():
                rows.append({"length_bin": label, "variable": name, "n": int(mask.sum()),
                             "n_pos": int(yb.sum()) if mask.any() else 0,
                             "prevalence": yb.mean() if mask.any() else np.nan})
                continue
            x = frame[name].to_numpy()[mask]
            auc, ap = auc_ap(x, yb)
            rows.append({"length_bin": label, "variable": name, "n": int(mask.sum()),
                         "n_pos": int(yb.sum()), "prevalence": yb.mean(), "auc": auc,
                         "ap": ap, "mi_bits_20bins": mi_from_codes(quantile_codes(x, 20), yb),
                         "h_label_bits": label_entropy(yb)})
    by_bin = pd.DataFrame(rows)

    # conditional MI given length
    strata_bins = bin_index
    strata20 = quantile_codes(length.astype(float), 20)
    h = label_entropy(y)
    i_lenbins = mi_from_codes(strata_bins, y)
    i_len20 = mi_from_codes(strata20, y)
    crow = []
    for name in ENTROPY_VARIABLES:
        x = frame[name].to_numpy()
        i_x = mi_from_codes(quantile_codes(x, 20), y)
        cb, _ = conditional_mi(x, y, strata_bins, 20)
        c20, _ = conditional_mi(x, y, strata20, 20)
        crow.append({"variable": name, "H_label_bits": h, "I_X_label": i_x,
                     "I_lenbins_label": i_lenbins, "I_X_label_given_lenbins": cb,
                     "I_len20_label": i_len20, "I_X_label_given_len20": c20,
                     "I_len20_and_X_label": i_len20 + c20,
                     "fraction_of_I_X_retained_given_len20": c20 / i_x})
    return by_bin, pd.DataFrame(crow)


def grouped_models(frame: pd.DataFrame, rng, genome_fraction: float,
                   folds: int) -> pd.DataFrame:
    genomes = frame.input_id.unique()
    chosen = set(rng.choice(genomes, size=int(len(genomes) * genome_fraction),
                            replace=False))
    sub = frame[frame.input_id.isin(chosen)].reset_index(drop=True)
    y = sub.in_genbank.to_numpy().astype(int)
    groups = sub.input_id.astype("category").cat.codes.to_numpy()
    sets = {
        "length": ["protein_length"],
        "3Di": ["three_di_entropy"],
        "length+3Di": ["protein_length", "three_di_entropy"],
        "length+protein+3Di": ["protein_length", "protein_entropy", "three_di_entropy"],
        "5 entropies": list(ENTROPY_VARIABLES),
        "length+5 entropies": ["protein_length", *ENTROPY_VARIABLES],
    }
    rows = []
    splitter = GroupKFold(n_splits=folds)
    split = list(splitter.split(np.zeros(y.size), y, groups))
    for name, columns in sets.items():
        x = sub[columns].to_numpy()
        oof = np.zeros(y.size)
        for train, test in split:
            model = HistGradientBoostingClassifier(
                max_iter=200, learning_rate=0.1, early_stopping=True,
                validation_fraction=0.1, random_state=20261014)
            model.fit(x[train], y[train])
            oof[test] = model.predict_proba(x[test])[:, 1]
        auc, ap = auc_ap(oof, y)
        rows.append({"model": name, "features": ",".join(columns), "n_orfs": y.size,
                     "n_genomes": len(chosen), "folds": folds, "auc_oof": auc,
                     "ap_oof": ap, "prevalence": y.mean()})
        print(rows[-1], flush=True)
    return pd.DataFrame(rows)


def stripe_table(frame: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for label, group in frame.groupby("in_genbank"):
        h = group.three_di_entropy.to_numpy()
        k = group.n_states_3di.to_numpy()
        rows.append({
            "in_genbank": label, "n": len(group),
            "frac_3di_H_eq_0": float((h <= 1e-12).mean()),
            "frac_3di_H_le_1": float((h <= 1 + 1e-9).mean()),
            "frac_3di_H_le_log2_3": float((h <= math.log2(3) + 1e-9).mean()),
            "frac_le_3_states": float((k <= 3).mean()),
            "median_states": float(np.median(k)),
            "mean_frac_D": float(group.frac_D_3di.mean()),
            "frac_D_majority": float((group.frac_D_3di > 0.5).mean()),
            **{f"mean_{v}": float(group[v].mean()) for v in VARIABLES},
            **{f"median_{v}": float(group[v].median()) for v in VARIABLES},
        })
    return pd.DataFrame(rows)


# ----------------------------------------------------------------------------
# Figures
# ----------------------------------------------------------------------------

def save(figure, stem: Path) -> None:
    figure.savefig(stem.with_suffix(".png"), dpi=300, bbox_inches="tight")
    figure.savefig(stem.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(figure)


def pretty(features: str) -> str:
    return " + ".join(SHORT[f] for f in features.split("+"))


def fig_mi(table: pd.DataFrame, out: Path) -> None:
    t = table[table.bins_per_axis == 20].copy()
    t["label"] = t.features.map(pretty)
    t = t.sort_values("mi_balanced_bits")
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(12.5, 8.5), width_ratios=[2.2, 1],
                                  constrained_layout=True)
    colours = [POS_COLOUR if k == "single" else "#7f7f7f" for k in t.kind]
    y = np.arange(len(t))
    ax.barh(y, t.mi_balanced_bits, color=colours, height=0.7)
    ax.errorbar(t.mi_balanced_bits, y,
                xerr=[(t.mi_balanced_bits - t.boot_lo_bits).clip(lower=0),
                      (t.boot_hi_bits - t.mi_balanced_bits).clip(lower=0)],
                fmt="none", ecolor="black", elinewidth=0.8, capsize=2)
    # best single component of each pair
    single = t[t.kind == "single"].set_index("features").mi_balanced_bits
    for yi, (_, row) in zip(y, t.iterrows()):
        if row.kind == "pair":
            ax.plot(single[row.best_single_balanced], yi, marker="|", color="black",
                    markersize=10, mew=1.5)
    null_hi = t.null_p99_bits.max()
    ax.axvspan(0, null_hi, color="#d9d9d9", zorder=0)
    ax.set_yticks(y, t.label)
    ax.set_xlabel("I(features; in_genbank) [bits], class-balanced sample (H = 1 bit)")
    ax.set_title("A  Single variables (blue) and pairs (grey), 20 quantile bins/axis",
                 loc="left", fontsize=11)
    ax.text(0.99, 0.02, f"grey band: permutation null, max 99th pct = {null_hi:.4f} bits\n"
            "| = best single component of the pair; whiskers = bootstrap 95% CI",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=8)
    ax.set_xlim(0, None)
    # bin sensitivity
    s = table.copy()
    s["label"] = s.features.map(pretty)
    keep = list(single.sort_values(ascending=False).index) + [
        "protein_entropy+three_di_entropy", "three_di_entropy+protein_length"]
    s = s[s.features.isin(keep)]
    sns.pointplot(data=s, x="bins_per_axis", y="mi_balanced_bits", hue="label",
                  ax=ax2, markersize=4, linewidth=1.2, palette="colorblind")
    ax2.set_xlabel("quantile bins per axis")
    ax2.set_ylabel("I [bits], balanced sample")
    ax2.set_title("B  Bin-count sensitivity", loc="left", fontsize=11)
    ax2.legend(fontsize=7, title=None, loc="upper left", bbox_to_anchor=(0, -0.12),
               ncol=1, frameon=False)
    save(fig, out / "fig1_mi_single_vs_pairs")


def fig_length(by_bin: pd.DataFrame, frame: pd.DataFrame, out: Path) -> None:
    b = by_bin[(by_bin.length_bin != "all") & by_bin.auc.notna()].copy()
    b["variable_label"] = b.variable.map(SHORT)
    order = [l for l in LENGTH_LABELS if l in set(b.length_bin)]
    fig, axes = plt.subplots(1, 3, figsize=(16, 4.8), constrained_layout=True)
    palette = dict(zip([SHORT[v] for v in VARIABLES], sns.color_palette("colorblind", 6)))
    for ax, metric, title in ((axes[0], "auc", "A  ROC AUC within length bin"),
                              (axes[1], "mi_bits_20bins", "B  I(X; label | length bin = l)")):
        sns.pointplot(data=b, x="length_bin", y=metric, hue="variable_label", order=order,
                      ax=ax, palette=palette, markersize=4, linewidth=1.2)
        ax.set_xlabel("ORF length (aa)")
        ax.set_title(title, loc="left", fontsize=11)
        ax.legend(fontsize=7, title=None)
    axes[0].axhline(0.5, color="black", lw=0.8, ls="--")
    axes[0].set_ylabel("AUC (in_genbank positive)")
    axes[1].set_ylabel("bits (20 within-bin quantile bins)")
    counts = frame.assign(length_bin=pd.cut(frame.protein_length, LENGTH_BINS,
                                            right=False, labels=LENGTH_LABELS))
    c = counts.groupby(["length_bin", "in_genbank"], observed=False).size().reset_index(name="n")
    c = c[c.length_bin.isin(order)]
    c["length_bin"] = c.length_bin.astype(str)
    c["label"] = c.in_genbank.map({True: "in_genbank", False: "not in_genbank"})
    sns.barplot(data=c, x="length_bin", y="n", hue="label", ax=axes[2], order=order,
                hue_order=["in_genbank", "not in_genbank"],
                palette={"in_genbank": POS_COLOUR, "not in_genbank": NEG_COLOUR})
    axes[2].set_yscale("log")
    axes[2].set_xlabel("ORF length (aa)")
    axes[2].set_ylabel("ORFs")
    axes[2].set_title("C  ORFs per length bin", loc="left", fontsize=11)
    axes[2].legend(fontsize=8, title=None)
    for ax in axes:
        ax.tick_params(axis="x", rotation=30)
    save(fig, out / "fig2_length_control")


def fig_density(frame: pd.DataFrame, out: Path) -> None:
    xbins = np.linspace(1.5, 4.4, 291)
    ybins = np.linspace(0, 4.4, 221)
    hist = {}
    for label in (True, False):
        g = frame[frame.in_genbank == label]
        h, _, _ = np.histogram2d(g.protein_entropy, g.three_di_entropy, bins=[xbins, ybins])
        hist[label] = h
    fig, axes = plt.subplots(1, 3, figsize=(17, 5.2), constrained_layout=True)
    for ax, label, title in ((axes[0], True, "A  in_genbank (phold CDS)"),
                             (axes[1], False, "B  not in_genbank")):
        h = hist[label]
        m = ax.pcolormesh(xbins, ybins, h.T, norm=LogNorm(vmin=1, vmax=h.max()),
                          cmap="viridis", rasterized=True)
        fig.colorbar(m, ax=ax, label="ORFs per cell")
        ax.set_title(f"{title}  (n = {int(h.sum()):,} in view)", loc="left", fontsize=11)
    # class-normalised densities with a 0.5 pseudocount per cell
    pos = (hist[True] + 0.5) / (hist[True].sum() + 0.5 * hist[True].size)
    neg = (hist[False] + 0.5) / (hist[False].sum() + 0.5 * hist[False].size)
    ratio = np.log2(pos / neg)
    ratio[(hist[True] + hist[False]) < 20] = np.nan
    lim = float(min(6.0, max(1.0, np.nanquantile(np.abs(ratio), 0.99))))
    m = axes[2].pcolormesh(xbins, ybins, ratio.T, cmap="RdBu_r",
                           norm=TwoSlopeNorm(0, -lim, lim), rasterized=True)
    fig.colorbar(m, ax=axes[2], label="log2 density ratio (in_genbank / not)\n"
                 "red: enriched in in_genbank; blue: in not in_genbank")
    axes[2].set_title("C  Enrichment (cells with ≥ 20 ORFs)", loc="left", fontsize=11)
    for ax in axes:
        ax.set_xlabel("protein entropy (bits)")
        ax.set_ylabel("3Di entropy (bits)")
        for k in (2, 3, 4):
            ax.axhline(math.log2(k), color="white" if ax is not axes[2] else "grey",
                       lw=0.5, ls=":")
    save(fig, out / "fig3_protein_vs_3di_density")


def fig_distributions(frame: pd.DataFrame, out: Path) -> None:
    fig, axes = plt.subplots(2, 1, figsize=(11, 7.5), constrained_layout=True)
    for ax, column, kmax, title in (
            (axes[0], "three_di_entropy", 20, "A  3Di entropy"),
            (axes[1], "twelve_state_entropy", 12, "B  12-state entropy")):
        bins = np.arange(0, math.log2(kmax) + 0.02, 0.01)
        for label, colour, name in ((True, POS_COLOUR, "in_genbank"),
                                    (False, NEG_COLOUR, "not in_genbank")):
            values = frame.loc[frame.in_genbank == label, column]
            ax.hist(values, bins=bins, density=True, histtype="step", color=colour,
                    lw=1.0, label=f"{name} (n = {len(values):,})")
        for k in range(2, kmax + 1):
            ax.axvline(math.log2(k), color="grey", lw=0.5, ls=":", zorder=0)
            if k <= 8 or k in (12, 16, 20):
                ax.text(math.log2(k), 1.0, f"log₂{k}", transform=ax.get_xaxis_transform(),
                        ha="center", va="bottom", fontsize=7, color="grey")
        ax.set_yscale("log")
        ax.set_xlabel({"three_di_entropy": "3Di entropy (bits)",
                       "twelve_state_entropy": "12-state entropy (bits)"}[column])
        ax.set_ylabel("density (log scale)")
        ax.set_title(title, loc="left", fontsize=11, pad=14)
        ax.legend(fontsize=8, loc="upper left")
    save(fig, out / "fig4_entropy_distributions_log2k")


def fig_models(models: pd.DataFrame, overall: pd.DataFrame, out: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.2), constrained_layout=True)
    o = overall.copy()
    o["label"] = o.variable.map(SHORT)
    o = o.sort_values("auc")
    axes[0].barh(o.label, o.auc, color=POS_COLOUR, height=0.6)
    axes[0].axvline(0.5, color="black", lw=0.8, ls="--")
    axes[0].set_xlim(0.3, 1)
    axes[0].set_xlabel("single-feature ROC AUC (all ORFs)")
    axes[0].set_title("A  Single features", loc="left", fontsize=11)
    for yi, v in enumerate(o.auc):
        axes[0].text(v + 0.005, yi, f"{v:.3f}", va="center", fontsize=8)
    m = models.melt(id_vars="model", value_vars=["auc_oof", "ap_oof"],
                    var_name="metric", value_name="value")
    m["metric"] = m.metric.map({"auc_oof": "ROC AUC", "ap_oof": "average precision"})
    sns.barplot(data=m, y="model", x="value", hue="metric", ax=axes[1],
                palette=["#0072B2", "#999999"])
    axes[1].axvline(models.prevalence.iloc[0], color="grey", lw=0.8, ls=":")
    axes[1].set_xlim(0, 1)
    axes[1].set_xlabel("out-of-fold score (genome-grouped 5-fold CV)")
    axes[1].set_ylabel("")
    axes[1].set_title("B  Gradient-boosting models", loc="left", fontsize=11)
    axes[1].legend(fontsize=8, title=None, loc="upper center", bbox_to_anchor=(0.5, -0.18),
                   ncol=2, frameon=False)
    save(fig, out / "fig5_auc_single_and_models")


# ----------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("table", type=Path)
    parser.add_argument("--tables", type=Path, required=True)
    parser.add_argument("--figures", type=Path, required=True)
    parser.add_argument("--per-class", type=int, default=200_000)
    parser.add_argument("--permutations", type=int, default=200)
    parser.add_argument("--bootstraps", type=int, default=100)
    parser.add_argument("--model-genome-fraction", type=float, default=0.25)
    parser.add_argument("--seed", type=int, default=20261014)
    parser.add_argument("--skip-models", action="store_true")
    parser.add_argument("--relabel", type=Path,
                        help="Parquet from relabel_phage_orfs.py (in_genbank_v2)")
    parser.add_argument("--label", choices=("legacy", "v2"), default="v2",
                        help="legacy = archived TSV label; v2 = current matcher")
    args = parser.parse_args()
    args.tables.mkdir(parents=True, exist_ok=True)
    args.figures.mkdir(parents=True, exist_ok=True)
    sns.set_theme(style="whitegrid", context="paper", font_scale=1.15)
    rng = np.random.default_rng(args.seed)

    frame = pd.read_parquet(args.table, columns=[
        "input_id", "orf_id", "start", "end", "strand", "frame", "in_genbank",
        "n_states_3di", "frac_D_3di", *VARIABLES])
    n_raw = len(frame)
    label_info = {"label": args.label}
    if args.relabel is not None:
        keys = ["input_id", "start", "end", "strand", "frame"]
        frame["strand"] = frame.strand.astype(str)
        relabel = pd.read_parquet(args.relabel, columns=[*keys, "in_genbank_v2"])
        frame = frame.merge(relabel, on=keys, how="left", validate="one_to_one")
        missing = int(frame.in_genbank_v2.isna().sum())
        legacy = frame.in_genbank.astype(bool)
        v2 = frame.in_genbank_v2.fillna(False).astype(bool)
        label_info.update({
            "relabel_missing": missing, "relabel_rows": len(relabel),
            "legacy_true": int(legacy.sum()), "v2_true": int(v2.sum()),
            "legacy_true_v2_false": int((legacy & ~v2).sum()),
            "legacy_false_v2_true": int((~legacy & v2).sum()),
            "legacy_true_minus_strand": int((legacy & (frame.strand == "-")).sum()),
            "v2_true_minus_strand": int((v2 & (frame.strand == "-")).sum()),
        })
        print(label_info, flush=True)
        if args.label == "v2":
            if missing:
                raise ValueError(f"{missing} ORFs lack a v2 label")
            frame["in_genbank"] = v2
        frame = frame.drop(columns=["in_genbank_v2"])
    elif args.label == "v2":
        raise ValueError("--label v2 requires --relabel")
    frame = frame.dropna(subset=["in_genbank", *VARIABLES]).reset_index(drop=True)
    frame["in_genbank"] = frame.in_genbank.astype(bool)
    frame["protein_length"] = frame.protein_length.astype(float)
    summary = {
        "orfs_total": n_raw, "orfs_dropped_missing": n_raw - len(frame),
        "orfs_used": len(frame), "genomes": int(frame.input_id.nunique()),
        "in_genbank_true": int(frame.in_genbank.sum()),
        "in_genbank_false": int((~frame.in_genbank).sum()),
        "prevalence": float(frame.in_genbank.mean()),
        "H_label_bits": label_entropy(frame.in_genbank.to_numpy().astype(int)),
        "min_protein_length": int(frame.protein_length.min()),
        "protein_length_quantiles": frame.protein_length.quantile(
            [0, 0.001, 0.01, 0.05, 0.5, 0.95, 0.99, 1]).to_dict(),
        "min_protein_length_in_genbank": int(frame.loc[frame.in_genbank, "protein_length"].min()),
        **label_info,
        "seed": args.seed, "per_class_sample": args.per_class,
        "permutations": args.permutations, "bootstraps": args.bootstraps,
    }
    print(json.dumps(summary, indent=1, default=str), flush=True)
    (args.tables / "dataset_summary.json").write_text(json.dumps(summary, indent=1, default=str))

    stripes = stripe_table(frame)
    stripes.to_csv(args.tables / "class_descriptives.tsv", sep="\t", index=False,
                   float_format="%.6g")
    print(stripes.T.to_string(), flush=True)

    sample = balanced_sample(frame, args.per_class, rng)
    mi = mi_tables(frame, sample, (10, 20, 40), args.permutations, args.bootstraps, rng)
    mi.to_csv(args.tables / "mi_single_and_pairs.tsv", sep="\t", index=False,
              float_format="%.6g")
    print(mi[mi.bins_per_axis == 20].sort_values("mi_balanced_bits", ascending=False)[[
        "features", "mi_balanced_bits", "mi_full_bits", "mi_full_fraction_of_H",
        "gain_over_best_single_balanced_bits", "null_p99_bits"]].to_string(), flush=True)

    by_bin, cond = length_tables(frame)
    by_bin.to_csv(args.tables / "by_length_bin_auc_ap_mi.tsv", sep="\t", index=False,
                  float_format="%.6g")
    cond.to_csv(args.tables / "conditional_mi_given_length.tsv", sep="\t", index=False,
                float_format="%.6g")
    print(by_bin.to_string(), flush=True)
    print(cond.to_string(), flush=True)

    fig_mi(mi, args.figures)
    fig_length(by_bin, frame, args.figures)
    fig_density(frame, args.figures)
    fig_distributions(frame, args.figures)

    if not args.skip_models:
        models = grouped_models(frame, rng, args.model_genome_fraction, 5)
        models.to_csv(args.tables / "grouped_cv_models.tsv", sep="\t", index=False,
                      float_format="%.6g")
        fig_models(models, by_bin[by_bin.length_bin == "all"], args.figures)
    print("done", flush=True)


if __name__ == "__main__":
    main()
