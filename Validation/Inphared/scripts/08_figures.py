#!/usr/bin/env python3
"""Figures and summary tables for the Inphared ORF-ranking analysis.

Reads the work directory (models_json/, models_relabelled/, tables/,
foldseek_<fs>/) and the committed evidence tables; writes PNG + PDF to
--figures and small TSVs to --tables.

  fig1_label_audit           JSON (0.1.15 matcher) vs recomputed labels by strand;
                             what the original 52,057 candidates are
  fig2_roc_pr                ROC and PR, JSON labels and recomputed labels
  fig3_beyond_length         pooled AUC/AP per feature set; AP within length bins
  fig4_importance_calibration gain importance and reliability curves
  fig5_score_distributions   OOF probability by ORF class (re-ranked model)
  fig6_foldseek_rates        full-length structural / sequence hit rates by arm and DB
  fig7_hit_rate_by_length    full-length structural hit rate by length bin
  fig8_replacements          replacement ORF vs the annotated ORF it would replace
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import seaborn as sns  # noqa: E402
from sklearn.metrics import precision_recall_curve, roc_curve  # noqa: E402

sns.set_theme(context="paper", style="whitegrid", font_scale=1.05)
PAL = sns.color_palette("colorblind")
SET_LABEL = {
    "orig14_stored": "original model (stored OOF, 14 features incl. start/end)",
    "orig14": "14 features incl. start/end",
    "no_coord": "no start/end (12)",
    "intrinsic": "intrinsic: entropies + lengths (8)",
    "length_only": "length only",
    "length_3di": "length + 3Di entropy",
    "length_entropies": "length + 4 entropies",
    "entropy_only": "4 entropies, no length",
}
SET_COLOR = {k: PAL[i] for i, k in enumerate(["orig14_stored", "orig14", "no_coord", "intrinsic", "length_only",
                                               "length_3di", "length_entropies", "entropy_only"])}
ARM_ORDER = ["annotated", "candidate", "frameshift_lo", "antisense_lo", "intergenic_lo", "candidate_shuffled_3di",
             "orig_candidate", "replacement", "replaced", "orig_replacement", "orig_replaced"]
ARM_LABEL = {"annotated": "annotated genes", "candidate": "candidates (re-ranked)",
             "frameshift_lo": "frame shadows, low p", "antisense_lo": "antisense, low p",
             "intergenic_lo": "intergenic, low p", "candidate_shuffled_3di": "candidates, shuffled 3Di",
             "orig_candidate": "original candidates", "replacement": "replacement ORFs",
             "replaced": "genes they replace", "orig_replacement": "original replacement ORFs",
             "orig_replaced": "genes they replace (orig.)"}
DB_LABEL = {"pdb100": "PDB100", "afdb_swissprot": "AFDB-SwissProt", "cath50": "CATH50", "bfvd": "BFVD",
            "any": "any DB"}


def bin_label(b: str) -> str:
    lo, hi = (int(float(x)) for x in b.strip("[)").split(","))
    return f">={lo}" if hi >= 10**8 else f"{lo}-{hi - 1}"


def bin_order(bins):
    return sorted(set(bins), key=lambda s: int(s.split("-")[0].lstrip(">=")))


def save(fig, out: Path, name: str):
    fig.savefig(out / f"{name}.png", dpi=200, bbox_inches="tight")
    fig.savefig(out / f"{name}.pdf", bbox_inches="tight")
    plt.close(fig)


def thin(x, y, n=2000):
    idx = np.unique(np.linspace(0, len(x) - 1, n).astype(int))
    return x[idx], y[idx]


def curves(ax_roc, ax_pr, oof, sets, title):
    y = oof.label.to_numpy()
    for s in sets:
        c = f"p_{s}"
        if c not in oof:
            continue
        p = oof[c].to_numpy()
        fpr, tpr, _ = roc_curve(y, p)
        prec, rec, _ = precision_recall_curve(y, p)
        ls = "--" if s == "orig14_stored" else "-"
        ax_roc.plot(*thin(fpr, tpr), color=SET_COLOR[s], lw=1.6, ls=ls, label=SET_LABEL[s])
        ax_pr.plot(*thin(rec, prec), color=SET_COLOR[s], lw=1.6, ls=ls, label=SET_LABEL[s])
    ax_roc.plot([0, 1], [0, 1], ls=":", color="grey", lw=0.8)
    ax_pr.axhline(y.mean(), ls=":", color="grey", lw=0.8)
    ax_roc.set(xlabel="false-positive rate", ylabel="true-positive rate", title=f"ROC: {title}")
    ax_pr.set(xlabel="recall", ylabel="precision", title=f"Precision-recall: {title}", ylim=(0, 1.02))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--work", required=True, type=Path)
    ap.add_argument("--fs", default="intrinsic")
    ap.add_argument("--figures", required=True, type=Path)
    ap.add_argument("--tables", required=True, type=Path)
    a = ap.parse_args()
    a.figures.mkdir(parents=True, exist_ok=True)
    a.tables.mkdir(parents=True, exist_ok=True)
    W = a.work

    # ---- summary tables from the model runs
    mets = []
    for regime in ("json", "relabelled"):
        m = pd.read_csv(W / f"models_{regime}" / "metrics.tsv", sep="\t").assign(labels=regime)
        mets.append(m)
    M = pd.concat(mets)
    M.to_csv(a.tables / "model_metrics.tsv", sep="\t", index=False, float_format="%.5f")
    pooled = M[M.fold == "pooled"].pivot_table(index="set", columns="labels", values=["auc", "ap", "brier"])
    pooled.columns = [f"{x}_{y}" for x, y in pooled.columns]
    pooled.to_csv(a.tables / "model_metrics_pooled.tsv", sep="\t", float_format="%.4f")
    for f in ("metrics_by_length", "thresholds", "importance", "calibration"):
        pd.concat([pd.read_csv(W / f"models_{r}" / f"{f}.tsv", sep="\t").assign(labels=r)
                   for r in ("json", "relabelled")]).to_csv(a.tables / f"model_{f}.tsv", sep="\t", index=False,
                                                            float_format="%.5f")
    checks = {r: json.load(open(W / f"models_{r}" / "fold_check.json")) for r in ("json", "relabelled")}
    json.dump(checks, open(a.tables / "model_fold_checks.json", "w"), indent=1)

    oj = pd.read_parquet(W / "models_json" / "oof.parquet")
    orl = pd.read_parquet(W / "models_relabelled" / "oof.parquet")
    o = pd.read_parquet(W / "tables" / "orfs.parquet", columns=["row", "strand", "context", "protein_length"]).sort_values(
        "row").reset_index(drop=True)
    lab = pd.read_parquet(W / "tables" / "labels.parquet").sort_values("row").reset_index(drop=True)

    # ---- fig1: label audit
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    t = pd.DataFrame({"strand": lab.strand, "JSON labels (0.1.15 matcher)": lab.in_genbank_orig,
                      "recomputed (0.2.0+ matcher)": lab.in_genbank})
    t = t.melt(id_vars="strand", var_name="labels", value_name="annotated")
    t = t[t.annotated].groupby(["labels", "strand"]).size().reset_index(name="n")
    sns.barplot(t, x="strand", y="n", hue="labels", ax=axes[0], palette=[PAL[7], PAL[0]])
    axes[0].set(ylabel="ORFs labelled in_genbank", xlabel="ORF strand",
                title="a  Annotated ORFs by strand")
    for c in axes[0].containers:
        axes[0].bar_label(c, fmt="{:,.0f}", fontsize=8)
    axes[0].legend(title=None, fontsize=8)
    oc = pd.read_csv(W / "rank_orfs_orig" / "ranked_putative_missing_orfs.tsv.gz", sep="\t", usecols=["input_id", "orf_id"])
    key = pd.read_parquet(W / "tables" / "orfs.parquet", columns=["row", "input_id", "orf_id"])
    oc = oc.merge(key, on=["input_id", "orf_id"]).merge(lab[["row", "in_genbank"]], on="row").merge(
        o[["row", "context"]], on="row")
    oc["class"] = np.where(oc.in_genbank, "annotated (recomputed)", oc.context)
    cc = oc["class"].value_counts()
    cc.to_frame("n").to_csv(a.tables / "original_candidates_recomputed.tsv", sep="\t")
    sns.barplot(x=cc.values, y=cc.index, ax=axes[1], color=PAL[1])
    axes[1].bar_label(axes[1].containers[0], fmt="{:,.0f}", fontsize=8)
    axes[1].set(xlabel="original candidates (n = 52,057)", ylabel="",
                title="b  Original candidates under recomputed labels")
    fig.tight_layout()
    save(fig, a.figures, "fig1_label_audit")

    # ---- fig2: ROC / PR
    fig, axes = plt.subplots(2, 2, figsize=(11, 9.5))
    curves(axes[0, 0], axes[0, 1], oj, ["orig14_stored", "no_coord", "intrinsic", "length_only"], "JSON labels")
    curves(axes[1, 0], axes[1, 1], orl, ["orig14", "no_coord", "intrinsic", "length_3di", "length_only"],
           "recomputed labels")
    for ax in axes.flat:
        ax.legend(fontsize=7, loc="lower right" if "ROC" in ax.get_title() else "lower left")
    fig.tight_layout()
    save(fig, a.figures, "fig2_roc_pr")

    # ---- fig3: beyond length
    fig, axes = plt.subplots(1, 3, figsize=(17, 5.2), gridspec_kw={"width_ratios": [1, 1, 1.1]})
    P = M[(M.fold == "pooled")].copy()
    P["set_label"] = P.set.map(SET_LABEL)
    order = ["length_only", "entropy_only", "length_3di", "length_entropies", "intrinsic", "no_coord", "orig14",
             "orig14_stored"]
    P = P[P.set.isin(order)]
    for ax, metric in zip(axes[:2], ("auc", "ap")):
        sns.barplot(P, y="set", x=metric, hue="labels", order=[s for s in order if s in set(P.set)], ax=ax,
                    palette=[PAL[7], PAL[0]])
        ax.set_yticks(range(len([s for s in order if s in set(P.set)])))
        ax.set_yticklabels([SET_LABEL[s] for s in order if s in set(P.set)], fontsize=8)
        ax.set(ylabel="", xlabel=metric.upper() + " (pooled out-of-fold)",
               xlim=(0.5 if metric == "auc" else 0, 1))
        for c in ax.containers:
            ax.bar_label(c, fmt="%.3f", fontsize=7)
        ax.legend(title="labels", fontsize=8, loc="lower left")
    axes[1].set_yticklabels([])
    axes[0].set_title("a  ROC AUC by feature set")
    axes[1].set_title("b  Average precision by feature set")
    L = pd.read_csv(W / "models_relabelled" / "metrics_by_length.tsv", sep="\t")
    L = L[L.set.isin(["length_only", "length_3di", "intrinsic"])].copy()
    L["length_bin"] = L.length_bin.map(bin_label)
    sns.pointplot(L, x="length_bin", y="ap", hue="set", ax=axes[2], order=bin_order(L.length_bin),
                  hue_order=["length_only", "length_3di", "intrinsic"],
                  palette={k: SET_COLOR[k] for k in ["length_only", "length_3di", "intrinsic"]})
    axes[2].set(xlabel="protein length bin (aa)", ylabel="AP within bin",
                title="c  Within length bins (recomputed labels)", ylim=(0, 1.02))
    h, l = axes[2].get_legend_handles_labels()
    axes[2].legend(h, [SET_LABEL[x] for x in l], fontsize=8, loc="lower right")
    fig.tight_layout()
    save(fig, a.figures, "fig3_beyond_length")

    # ---- fig4: importance + calibration
    imp = pd.concat([pd.read_csv(W / f"models_{r}" / "importance.tsv", sep="\t").assign(labels=r)
                     for r in ("json", "relabelled")])
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.6))
    for ax, (r, s, ttl) in zip(axes[:2], (("json", "orig14", "a  Gain importance: original 14-feature model, JSON labels"),
                                          ("relabelled", "intrinsic", "b  Gain importance: intrinsic model, recomputed labels"))):
        x = imp[(imp.labels == r) & (imp.set == s)]
        ordr = x.groupby("feature").gain.mean().sort_values(ascending=False).index
        sns.barplot(x, y="feature", x="gain", order=ordr, ax=ax, color=PAL[0], errorbar=("pi", 100))
        ax.set(title=ttl, xlabel="normalised gain (mean, range over 5 folds)", ylabel="")
        ax.title.set_fontsize(9)
    cal = pd.read_csv(W / "models_relabelled" / "calibration.tsv", sep="\t")
    for s in ["intrinsic", "length_only", "orig14"]:
        c = cal[cal.set == s]
        axes[2].plot(c.mean_pred, c.frac_pos, marker="o", ms=3, color=SET_COLOR[s], label=SET_LABEL[s])
    cj = pd.read_csv(W / "models_json" / "calibration.tsv", sep="\t")
    c = cj[cj.set == "orig14_stored"]
    axes[2].plot(c.mean_pred, c.frac_pos, marker="o", ms=3, color=SET_COLOR["orig14_stored"], ls="--",
                 label="original model vs JSON labels")
    axes[2].plot([0, 1], [0, 1], ls=":", color="grey")
    axes[2].set(xlabel="mean predicted probability (bin)", ylabel="fraction annotated",
                title="c  Calibration (recomputed labels unless stated)")
    axes[2].legend(fontsize=7)
    fig.tight_layout()
    save(fig, a.figures, "fig4_importance_calibration")

    # ---- fig5: score distributions (re-ranked model)
    pr = pd.read_parquet(W / f"rank_orfs_relabelled_{a.fs}" / "all_orf_oof_predictions.parquet",
                         columns=["probability_in_genbank"])
    d = pd.DataFrame({"p": pr.probability_in_genbank.to_numpy(), "annotated": lab.in_genbank.to_numpy(),
                      "context": o.context.to_numpy(), "length": o.protein_length.to_numpy()})
    d["class"] = np.where(d.annotated, "annotated", "unannotated: " + d.context)
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.6))
    corder = ["annotated", "unannotated: same_stop", "unannotated: same_frame", "unannotated: frameshift",
              "unannotated: antisense", "unannotated: intergenic"]
    corder = [c for c in corder if c in set(d["class"])]
    sns.violinplot(d.sample(min(len(d), 400000), random_state=1), x="p", y="class", order=corder, ax=axes[0],
                   cut=0, density_norm="width", inner="quartile", color=PAL[0], linewidth=0.7)
    axes[0].axvline(0.8, ls="--", color="k", lw=0.8)
    counts = d["class"].value_counts()
    axes[0].set_yticks(range(len(corder)))
    axes[0].set_yticklabels([f"{c} (n={counts[c]:,})" for c in corder], fontsize=8)
    axes[0].set(xlabel="out-of-fold P(annotated), intrinsic model", ylabel="",
                title="a  Score by ORF class (400k-ORF sample)")
    cand = d[(~d.annotated) & (d.p >= 0.8)]
    sns.histplot(cand, x="length", hue="context", multiple="stack", ax=axes[1], bins=np.arange(90, 610, 10),
                 palette="colorblind")
    axes[1].set(xlabel="protein length (aa)", ylabel="candidates (p >= 0.8)",
                title=f"b  Re-ranked candidates by context (n = {len(cand):,})", xlim=(90, 600))
    fig.tight_layout()
    save(fig, a.figures, "fig5_score_distributions")
    d.groupby("class").p.describe(percentiles=[0.05, 0.25, 0.5, 0.75, 0.95]).to_csv(
        a.tables / "score_by_class.tsv", sep="\t", float_format="%.4f")

    # ---- Foldseek figures (if the evidence tables exist)
    ev = a.tables / "evidence_rates.tsv"
    if ev.exists():
        R = pd.read_csv(ev, sep="\t")
        arms = [x for x in ARM_ORDER[:6] if x in set(R.arm)]
        fig, axes = plt.subplots(1, 2, figsize=(14, 5), sharey=False)
        for ax, mode, ttl in zip(axes, ("struct", "seq"), ("a  Structure (3Di+AA) search", "b  Sequence-only (AA) search")):
            x = R[(R["mode"] == mode) & R.arm.isin(arms)].copy()
            x["DB"] = x.db.map(DB_LABEL)
            sns.barplot(x, x="DB", y="fullcov_rate", hue="arm", hue_order=[h for h in arms if h in set(x.arm)],
                        order=[DB_LABEL[k] for k in ["pdb100", "afdb_swissprot", "cath50", "bfvd", "any"]],
                        ax=ax, palette="colorblind")
            ax.set(ylabel="full-length hit rate (qcov & tcov >= 0.8)", xlabel="", title=ttl, ylim=(0, 1))
            h, l = ax.get_legend_handles_labels()
            ax.legend(h, [ARM_LABEL[i] for i in l], fontsize=7, title=None)
        fig.tight_layout()
        save(fig, a.figures, "fig6_foldseek_rates")

        Lb = pd.read_csv(a.tables / "evidence_by_length.tsv", sep="\t")
        Lb = Lb[Lb.arm.isin(arms)].copy()
        Lb["length_bin"] = Lb.length_bin.map(bin_label)
        fig, ax = plt.subplots(figsize=(7.5, 4.6))
        sns.pointplot(Lb, x="length_bin", y="fullcov_struct_rate", hue="arm", order=bin_order(Lb.length_bin),
                      hue_order=[h for h in arms if h in set(Lb.arm)], ax=ax, palette="colorblind")
        ax.set(xlabel="protein length bin (aa)", ylabel="full-length structural hit rate (any DB)",
               title="Full-length structural hit rate by length", ylim=(0, 1))
        h, l = ax.get_legend_handles_labels()
        ax.legend(h, [ARM_LABEL[i] for i in l], fontsize=7)
        fig.tight_layout()
        save(fig, a.figures, "fig7_hit_rate_by_length")

        rp = a.work / f"foldseek_{a.fs}" / "replacement_pairs.tsv"
        if rp.exists():
            Pp = pd.read_csv(rp, sep="\t")
            fig, axes = plt.subplots(1, 2, figsize=(11, 4.4))
            t = Pp.groupby(["run", "outcome"]).size().reset_index(name="pairs")
            t["share"] = t.pairs / t.groupby("run").pairs.transform("sum")
            sns.barplot(t, x="outcome", y="share", hue="run", ax=axes[0], palette=[PAL[0], PAL[7]],
                        order=["replacement only", "both", "replaced only", "neither"])
            axes[0].set(ylabel="share of pairs", xlabel="full-length structural hit (any DB)",
                        title="a  Which ORF of the pair has a full-length hit")
            for c in axes[0].containers:
                axes[0].bar_label(c, fmt="%.2f", fontsize=7)
            for run, col in (("reranked", PAL[0]), ("original", PAL[7])):
                z = Pp[Pp.run == run]
                axes[1].scatter(z.replaced_bits, z.replacement_bits, s=6, alpha=0.5, color=col,
                                label=f"{run} (n={len(z):,})")
            mx = max(Pp.replaced_bits.max(), Pp.replacement_bits.max(), 10)
            axes[1].plot([0, mx], [0, mx], ls=":", color="grey")
            axes[1].set(xscale="symlog", yscale="symlog", xlabel="best structural bits: annotated gene",
                        ylabel="best structural bits: replacement ORF", title="b  Best-hit bit scores per pair")
            axes[1].legend(fontsize=8)
            fig.tight_layout()
            save(fig, a.figures, "fig8_replacements")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
