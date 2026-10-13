#!/usr/bin/env python3
"""Figures for supplementary_analyses.md, from the TSVs 34_supplementary_arms.py wrote.

Nothing is recomputed here: every panel reads a table that 34 emitted, so the
figure and the quoted number cannot drift apart.

  35_supplementary_figures.py --in-dir work/gtdb_extras/out --out-dir figures
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
# Validated categorical order (dataviz reference palette; CVD-checked).
ARM_COL = {"candidate": "#2a78d6", "clean_shadow": "#eb6834",
           "annotated_cds": "#1baf7a", "intergenic_lo": "#eda100",
           "unannot_hi": "#e87ba4", "shadow_hi": "#4a3aa7"}
ARM_LAB = {"candidate": "candidate", "clean_shadow": "clean shadow",
           "shadow_hi": "all shadows", "annotated_cds": "annotated CDS",
           "intergenic_lo": "intergenic_lo", "unannot_hi": "unannot_hi"}
ARMS = ["candidate", "clean_shadow", "shadow_hi", "annotated_cds",
        "intergenic_lo", "unannot_hi"]
DOM = {"bac": "Bacteria", "arc": "Archaea"}


def style():
    sns.set_theme(style="whitegrid", context="paper")
    plt.rcParams.update({
        "figure.dpi": 160, "savefig.dpi": 160, "font.size": 8.5,
        "axes.labelsize": 9, "axes.titlesize": 9.5, "axes.edgecolor": MUTED,
        "axes.labelcolor": INK, "text.color": INK, "xtick.color": MUTED,
        "ytick.color": MUTED, "axes.spines.top": False,
        "axes.spines.right": False, "figure.facecolor": PLATE,
        "savefig.facecolor": PLATE, "axes.facecolor": PLATE,
        "grid.color": "#e6e3de", "grid.linewidth": 0.6,
        "legend.frameon": False})


def save(fig, out, name):
    p = Path(out) / name
    fig.savefig(p, bbox_inches="tight")
    plt.close(fig)
    print(f"  {name} ({p.stat().st_size/1024:.0f} kB)")


def load(ind, name):
    return pd.concat([pd.read_csv(Path(ind) / d / name, sep="\t")
                      for d in ("bac", "arc")], ignore_index=True)


def line_arm(ax, x, y, arm, **kw):
    ls = "--" if arm == "shadow_hi" else "-"
    ax.plot(x, y, ls, color=ARM_COL[arm], lw=2, marker="o", ms=4.5,
            label=ARM_LAB[arm], **kw)


# ---------------------------------------------------------------------------
def fig_3di(ind, out):
    t = load(ind, "hit_by_3di_bin.tsv")
    t = t[t.scheme == "candidate_quintiles"]
    arms = ["candidate", "clean_shadow", "shadow_hi", "annotated_cds", "unannot_hi"]
    cols = [("any_struct_hit", "any structural hit (E < 1e-3, 4 DBs)"),
            ("est_struct_full", "full-length hit (qcov & tcov ≥ 0.8;\nSwiss-Prot/PDB/CATH)"),
            ("prodigal", "Prodigal 3′-exact agreement")]
    fig, axes = plt.subplots(2, 4, figsize=(13.5, 6.4))
    for r, dom in enumerate(("bac", "arc")):
        d = t[t.domain == dom]
        bins = d[d.arm == "candidate"].sort_values("bin")
        xl = [f"{q}\n{lo:.2f}–{hi:.2f}" for q, lo, hi in
              zip(bins.bin, bins.lo.clip(lower=2.5), bins.hi)]
        x = np.arange(len(xl))
        for c, (metric, title) in enumerate(cols):
            ax = axes[r, c]
            for a in arms:
                s = d[d.arm == a].sort_values("bin")
                if s[metric].notna().any():
                    line_arm(ax, x, s[metric] * 100, a)
            ax.set_xticks(x, xl, fontsize=7)
            ax.set_ylim(0, 102)
            ax.set_ylabel("% of arm")
            ax.set_title(f"{DOM[dom]}: {title}", loc="left")
        ax = axes[r, 3]
        sh = d[d.arm == "share_sp_clean"].sort_values("bin")
        ax.errorbar(x, sh.est_struct_full * 100,
                    yerr=[(sh.est_struct_full - sh.sp_struct_full) * 100,
                          (sh.prodigal - sh.est_struct_full) * 100],
                    fmt="o-", color=ARM_COL["candidate"], lw=2, ms=4.5, capsize=3)
        ax.axhline(0, color=MUTED, lw=0.8)
        ax.set_xticks(x, xl, fontsize=7)
        ax.set_ylabel("implied real-gene share (%)")
        ax.set_title(f"{DOM[dom]}: mixture share within bin\n(Swiss-Prot full-length vs clean shadow)",
                     loc="left")
    axes[0, 0].legend(loc="upper left", fontsize=7.5)
    for ax in axes[1]:
        ax.set_xlabel("candidate 3Di-entropy quintile (bits)")
    fig.suptitle("Any-hit rate saturates with 3Di entropy in every arm; full-length coverage "
                 "separates candidates from clean shadows only below the top quintile", x=0.01, ha="left",
                 fontsize=11, fontweight="bold")
    fig.tight_layout()
    save(fig, out, "supp_hit_rate_by_3di.png")


def fig_quality(ind, out):
    s = pd.concat([pd.read_csv(Path(ind) / d / "hit_quality_sample.tsv.gz", sep="\t")
                   for d in ("bac", "arc")])
    s = s[s.db == "afdb_swissprot"]
    arms = ["candidate", "clean_shadow", "shadow_hi", "annotated_cds",
            "intergenic_lo", "unannot_hi"]
    mets = [("fident", "sequence identity of best hit"),
            ("qcov", "query coverage"), ("tcov", "target coverage"),
            ("nlog10e", "−log10 E-value")]
    fig, axes = plt.subplots(2, 4, figsize=(14, 6.2), sharex=True)
    for r, dom in enumerate(("bac", "arc")):
        d = s[s.domain == dom]
        for c, (m, lab) in enumerate(mets):
            ax = axes[r, c]
            sns.violinplot(data=d, x="arm", y=m, hue="arm", order=arms,
                           hue_order=arms, palette=ARM_COL, cut=0, inner="quartile",
                           linewidth=0.7, density_norm="width", ax=ax, legend=False)
            ax.set_xticks(range(len(arms)), [ARM_LAB[a] for a in arms],
                          rotation=35, ha="right", fontsize=7.5)
            ax.set_xlabel("")
            ax.set_ylabel(lab)
            if m == "nlog10e":
                ax.set_ylim(0, np.nanquantile(d[m], 0.99))
            ax.set_title(f"{DOM[dom]}: {lab}", loc="left")
    fig.suptitle("Best structural hit to AFDB/Swiss-Prot, queries with a hit "
                 "(random ≤20,000 per arm; lines = quartiles)", x=0.01, ha="left",
                 fontsize=11, fontweight="bold")
    fig.tight_layout()
    save(fig, out, "supp_best_hit_quality.png")


def fig_evidence(ind, out):
    t = load(ind, "evidence_classes.tsv")
    t = t[t.db == "union"]
    arms = ["candidate", "clean_shadow", "shadow_hi", "annotated_cds",
            "intergenic_lo", "unannot_hi"]
    parts = [("both", "structure + sequence", "#2a78d6"),
             ("struct_only", "structure only", "#86b6ef"),
             ("seq_only", "sequence only", "#eda100")]
    fig, axes = plt.subplots(2, 2, figsize=(11, 6), sharey=True)
    for r, dom in enumerate(("bac", "arc")):
        for c, (lvl, title) in enumerate((("hit", "any hit (E < 1e-3)"),
                                          ("full", "full-length (qcov & tcov ≥ 0.8)"))):
            ax = axes[r, c]
            d = t[(t.domain == dom) & (t.level == lvl)].set_index("arm").loc[arms]
            left = np.zeros(len(arms))
            y = np.arange(len(arms))[::-1]
            for k, lab, col in parts:
                v = d[k].to_numpy() * 100
                ax.barh(y, v, left=left, color=col, height=0.66, label=lab,
                        edgecolor=PLATE, linewidth=1.5)
                left += v
            for yy, tot in zip(y, left):
                ax.text(tot + 1, yy, f"{tot:.1f}%", va="center", fontsize=7.5)
            ax.set_yticks(y, [ARM_LAB[a] for a in arms])
            ax.set_xlim(0, 100)
            ax.set_xlabel("% of arm (union over four databases)")
            ax.set_title(f"{DOM[dom]}: {title}", loc="left")
    axes[0, 1].legend(loc="lower right", fontsize=7.5)
    fig.suptitle("Structure-only versus sequence-supported homology by arm",
                 x=0.01, ha="left", fontsize=11, fontweight="bold")
    fig.tight_layout()
    save(fig, out, "supp_evidence_classes.png")


def fig_unannot(ind, out):
    t = load(ind, "unannot_hi.tsv")
    dbs = ["afdb_swissprot", "pdb100", "cath50", "bfvd", "union_est"]
    arms = ["intergenic_lo", "clean_shadow", "candidate", "unannot_hi", "annotated_cds"]
    fig, axes = plt.subplots(1, 2, figsize=(12, 3.8), sharey=True)
    for ax, dom in zip(axes, ("bac", "arc")):
        d = t[t.domain == dom].set_index("db").loc[dbs]
        x = np.arange(len(dbs))
        wdt = 0.16
        for i, a in enumerate(arms):
            ax.bar(x + (i - 2) * wdt, d[f"full_{a}"] * 100, wdt * 0.92,
                   color=ARM_COL[a], label=ARM_LAB[a])
        for xx, (_, r) in zip(x, d.iterrows()):
            ax.text(xx + wdt, r.full_unannot_hi * 100 + 1.2,
                    f"{r.unannot_share*100:.0f}%", ha="center", fontsize=7,
                    color=INK)
        ax.set_xticks(x, ["Swiss-Prot", "PDB100", "CATH50", "BFVD",
                          "any of SP/PDB/CATH"])
        ax.set_ylabel("full-length hit rate (%)")
        ax.set_title(f"{DOM[dom]}", loc="left", fontweight="bold")
    axes[0].legend(fontsize=7.5, loc="upper left", ncol=5, bbox_to_anchor=(0, -0.12))
    fig.suptitle("unannot_hi (high-3Di ORFs from never-annotated genomes) sits close "
                 "to annotated CDS; label above its bar = position on the "
                 "clean-shadow→annotated axis", x=0.01, ha="left", fontsize=10,
                 fontweight="bold")
    fig.tight_layout()
    save(fig, out, "supp_unannot_hi.png")


def fig_truncation(ind, out):
    t = load(ind, "truncation.tsv")
    arms = ["candidate", "clean_shadow", "shadow_hi", "annotated_cds",
            "intergenic_lo", "unannot_hi"]
    fcl = ["opposite strand", "same strand, frameshift",
           "same strand, frame undefined", "same strand, same frame"]
    fig, axes = plt.subplots(2, 2, figsize=(12, 6.4))
    for r, dom in enumerate(("bac", "arc")):
        d = t[t.domain == dom]
        ax = axes[r, 0]
        tr = d[d.truncated & d.arm.isin(arms)].set_index("arm").loc[arms]
        fc = d[d.arm.str.startswith("shadow:")].copy()
        fc["k"] = fc.arm.str.replace("shadow:", "", regex=False)
        fc = fc.set_index("k").reindex(fcl)
        labs = [ARM_LAB[a] for a in arms] + [f"shadow: {k}" for k in fcl]
        vals = list(tr.share * 100) + list(fc.share * 100)
        cols = [ARM_COL[a] for a in arms] + ["#4a3aa7"] * 4
        y = np.arange(len(labs))[::-1]
        ax.barh(y, vals, color=cols, height=0.66, alpha=1.0)
        for yy, v in zip(y, vals):
            ax.text(v + 0.8, yy, f"{v:.1f}%", va="center", fontsize=7.3)
        ax.set_yticks(y, labs, fontsize=7.5)
        ax.set_xlabel("contig-truncated (no stop codon), % of group")
        ax.set_title(f"{DOM[dom]}: truncation by arm and shadow frame class", loc="left")
        ax = axes[r, 1]
        x = np.arange(len(arms))
        comp = d[(~d.truncated) & d.arm.isin(arms)].set_index("arm").loc[arms]
        ax.bar(x - 0.2, comp.any_struct_hit * 100, 0.38, color="#86b6ef",
               label="complete ORF")
        ax.bar(x + 0.2, tr.any_struct_hit * 100, 0.38, color="#184f95",
               label="truncated ORF")
        ax.set_xticks(x, [ARM_LAB[a] for a in arms], rotation=30, ha="right")
        ax.set_ylabel("any structural hit (%)")
        ax.set_ylim(0, 100)
        ax.set_title(f"{DOM[dom]}: truncated ORFs hit far more often", loc="left")
        ax.legend(fontsize=7.5, loc="upper right")
    fig.suptitle("Truncation is concentrated in same-frame / frame-undefined shadows "
                 "— the ones the clean comparator removes", x=0.01, ha="left",
                 fontsize=11, fontweight="bold")
    fig.tight_layout()
    save(fig, out, "supp_truncation.png")


def fig_assembly(ind, out):
    t = load(ind, "assembly_strata.tsv")
    t = t[t.scheme == "pooled_query_rows"]
    fig, axes = plt.subplots(4, 3, figsize=(13, 12))
    r = 0
    for cov, clab in (("n50_contigs", "contig N50 quartile"),
                      ("checkm2_completeness", "CheckM2 completeness quartile")):
        for dom in ("bac", "arc"):
            d = t[(t.domain == dom) & (t.covariate == cov)]
            bins = d[d.arm == "candidate"].sort_values("bin")
            fmt = (lambda v: f"{v/1e3:,.0f}k") if cov == "n50_contigs" else (lambda v: f"{v:.0f}")
            xl = [f"{q}\n{fmt(lo)}–{fmt(hi)}" for q, lo, hi in zip(bins.bin, bins.lo, bins.hi)]
            x = np.arange(len(xl))
            for c, (metric, title, arms) in enumerate((
                    ("any_struct_hit", "any hit (published readout)",
                     ["candidate", "shadow_hi", "clean_shadow", "annotated_cds"]),
                    ("sp_struct_full", "Swiss-Prot full-length",
                     ["candidate", "shadow_hi", "clean_shadow", "annotated_cds"]),
                    ("prodigal", "Prodigal agreement",
                     ["candidate", "shadow_hi", "clean_shadow", "annotated_cds",
                      "intergenic_lo"]))):
                ax = axes[r, c]
                for a in arms:
                    s = d[d.arm == a].sort_values("bin")
                    line_arm(ax, x, s[metric] * 100, a)
                ax.set_xticks(x, xl, fontsize=7)
                ax.set_ylabel("% of arm")
                ax.set_ylim(0, None)
                ax.set_title(f"{DOM[dom]}, {clab}:\n{title}", loc="left", fontsize=8.5)
            sh = d[d.arm == "share_sp_clean"].sort_values("bin")
            axes[r, 1].set_title(
                f"{DOM[dom]}, {clab}: Swiss-Prot full-length\nimplied share by quartile: "
                + " / ".join(f"{v*100:.0f}%" for v in sh.est_struct_full),
                loc="left", fontsize=8.5)
            r += 1
    axes[0, 2].legend(fontsize=7.5, loc="center right")
    fig.suptitle("Assembly-quality stratification (quartiles over pooled query rows, "
                 "as in 18_pilot_analysis.py); share = mixture vs clean shadow",
                 x=0.01, ha="left", fontsize=11, fontweight="bold")
    fig.tight_layout()
    save(fig, out, "supp_assembly_strata.png")


CLASS_ORDER = ["metabolism", "other_named", "mobile_element", "transport",
               "translation", "regulation", "uncharacterized", "other classes",
               "unclassified"]
CLASS_COL = {"metabolism": "#2a78d6", "other_named": "#eb6834",
             "mobile_element": "#1baf7a", "transport": "#eda100",
             "translation": "#e87ba4", "regulation": "#008300",
             "uncharacterized": "#4a3aa7", "other classes": "#e34948",
             "unclassified": "#c9c5bf"}


def fig_function(ind, out):
    t = load(ind, "functional_composition_sets.tsv")
    t = t[t.kind == "class"].copy()
    t["key"] = t.key.where(~t.key.isin(["defence", "replication_repair",
                                        "cell_envelope"]), "other classes")
    t = t.groupby(["domain", "set", "n_set", "key"], as_index=False).share.sum()
    fl = load(ind, "functional_composition_sets.tsv")
    fl = fl[fl.key == "BFVD-only full-length"]
    fig, axes = plt.subplots(1, 2, figsize=(12, 3.6), sharex=True)
    sets = ["all candidates", "C3", "C7"]
    for ax, dom in zip(axes, ("bac", "arc")):
        d = t[t.domain == dom]
        y = np.arange(len(sets))[::-1]
        left = np.zeros(len(sets))
        for k in CLASS_ORDER:
            v = np.array([d[(d.set == s) & (d.key == k)].share.sum() * 100 for s in sets])
            ax.barh(y, v, left=left, color=CLASS_COL[k], height=0.62,
                    edgecolor=PLATE, linewidth=1.2, label=k.replace("_", " "))
            for yy, l_, vv in zip(y, left, v):
                if vv >= 6:
                    ax.text(l_ + vv / 2, yy, f"{vv:.0f}", ha="center", va="center",
                            fontsize=7, color="white" if k != "unclassified" else INK)
            left += v
        ns = [int(d[d.set == s].n_set.iloc[0]) for s in sets]
        bf = [float(fl[(fl.domain == dom) & (fl.set == s)].share.iloc[0]) * 100 for s in sets]
        ax.set_yticks(y, [f"{s}\nn = {n:,}\nBFVD-only {b:.1f}%" for s, n, b in zip(sets, ns, bf)],
                      fontsize=7.5)
        ax.set_xlim(0, 100)
        ax.set_xlabel("% of set (class from best Swiss-Prot/PDB structural hit text)")
        ax.set_title(DOM[dom], loc="left", fontweight="bold")
    axes[1].legend(fontsize=7, ncol=1, loc="center left", bbox_to_anchor=(1.01, 0.5))
    fig.suptitle("Functional composition beyond the strict set: all candidates, C3 "
                 "(any full-length hit, not truncated) and C7", x=0.01, ha="left",
                 fontsize=11, fontweight="bold")
    fig.tight_layout()
    save(fig, out, "supp_functional_sets.png")


OFF_ORDER = ["0 (exact 3' match)", "±1-2 bp (other frame)",
             "other, same frame (multiple of 3)", "other, different frame",
             "no same-strand Prodigal gene on contig"]
OFF_COL = ["#2a78d6", "#eb6834", "#1baf7a", "#c9c5bf", "#eda100"]


def fig_offsets(ind, out):
    t = load(ind, "prodigal_offsets.tsv")
    groups = ["annotated_cds", "candidate", "clean_shadow", "intergenic_lo",
              "shadow:same strand, same frame", "shadow:same strand, frame undefined",
              "shadow:same strand, frameshift", "shadow:opposite strand", "shadow_hi"]
    lab = {g: ARM_LAB.get(g, g) for g in groups}
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.2), sharey=True)
    for ax, dom in zip(axes, ("bac", "arc")):
        d = t[t.domain == dom].pivot_table(index="group", columns="offset_class",
                                           values="share").reindex(groups).fillna(0)
        y = np.arange(len(groups))[::-1]
        left = np.zeros(len(groups))
        for k, col in zip(OFF_ORDER, OFF_COL):
            v = d.get(k, pd.Series(0, index=d.index)).to_numpy() * 100
            ax.barh(y, v, left=left, color=col, height=0.66, label=k,
                    edgecolor=PLATE, linewidth=1.2)
            for yy, l_, vv in zip(y, left, v):
                if vv >= 7:
                    ax.text(l_ + vv / 2, yy, f"{vv:.0f}", ha="center", va="center",
                            fontsize=7, color="white" if col != "#c9c5bf" else INK)
            left += v
        ax.set_yticks(y, [lab[g] for g in groups], fontsize=7.6)
        ax.set_xlim(0, 100)
        ax.set_xlabel("% of group: offset to nearest same-strand Prodigal 3′ end")
        ax.set_title(DOM[dom], loc="left", fontweight="bold")
    axes[1].legend(fontsize=7, loc="center left", bbox_to_anchor=(1.01, 0.5))
    fig.suptitle("Prodigal agreement by arm, shadow frame class and 3′-offset class",
                 x=0.01, ha="left", fontsize=11, fontweight="bold")
    fig.tight_layout()
    save(fig, out, "supp_prodigal_offsets.png")


def fig_standardized(ind, out):
    t = load(ind, "standardized_mixture.tsv")
    lad = load(ind, "standardized_ladder.tsv")
    reads = ["afdb_swissprot", "pdb100", "cath50", "est_union", "prodigal"]
    rlab = ["Swiss-Prot", "PDB100", "CATH50", "any SP/PDB/CATH", "Prodigal"]
    fig, axes = plt.subplots(1, 3, figsize=(14, 4), gridspec_kw={"width_ratios": [1, 1, 1.1]})
    for ax, dom in zip(axes[:2], ("bac", "arc")):
        d = t[t.domain == dom].set_index("readout").loc[reads]
        x = np.arange(len(reads))
        for off, (k, lab, col) in zip((-0.27, 0, 0.27), (
                ("share_crude", "published (crude clean comparator)", "#86b6ef"),
                ("share_shadow_std", "clean shadow re-weighted", "#2a78d6"),
                ("share_both_std", "clean shadow + annotated re-weighted", "#184f95"))):
            v = d[k] * 100
            ax.bar(x + off, v, 0.26, color=col, label=lab)
            ax.errorbar(x + off, v, yerr=[v - d[f"{k}_lo"] * 100, d[f"{k}_hi"] * 100 - v],
                        fmt="none", ecolor=INK, lw=0.8, capsize=2)
        ax.set_xticks(x, rlab, rotation=20, ha="right")
        ax.set_ylabel("implied real-gene share (%)")
        ax.set_title(f"{DOM[dom]}: mixture share, re-weighted to candidate length × 3Di",
                     loc="left", fontsize=8.8)
    axes[0].legend(fontsize=7, loc="upper left")
    ax = axes[2]
    crit = ["C3", "C4", "C5", "C6", "C7~"]
    x = np.arange(len(crit))
    for off, dom, cols in ((-0.2, "bac", ("#86b6ef", "#2a78d6")),
                           (0.2, "arc", ("#f7b79c", "#eb6834"))):
        d = lad[lad.domain == dom].set_index("criterion").loc[crit]
        tot = 1 if dom == "bac" else 1
        ax.bar(x + off, d.excess_crude, 0.38, color=cols[0],
               label=f"{DOM[dom]} published excess")
        ax.bar(x + off, d.excess_std, 0.24, color=cols[1],
               label=f"{DOM[dom]} re-weighted excess")
    ax.set_yscale("log")
    ax.set_xticks(x, ["C3 any full", "C4 ≥2 DBs", "C5 SP/PDB/CATH", "C6 struct-only", "C7 strict"],
                  rotation=20, ha="right")
    ax.set_ylabel("candidates in excess of clean shadows (log)")
    ax.set_title("Direct-evidence excess, crude vs re-weighted", loc="left")
    ax.legend(fontsize=6.8, loc="upper right")
    fig.suptitle("The clean comparator is no longer length/3Di-matched; re-weighting it "
                 "lowers the structural estimate far more than the Prodigal one",
                 x=0.01, ha="left", fontsize=11, fontweight="bold")
    fig.tight_layout()
    save(fig, out, "supp_standardized_mixture.png")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--in-dir", required=True)
    ap.add_argument("--out-dir", required=True)
    a = ap.parse_args()
    style()
    for f in (fig_3di, fig_quality, fig_evidence, fig_unannot, fig_truncation,
              fig_assembly, fig_function, fig_offsets, fig_standardized):
        f(a.in_dir, a.out_dir)


if __name__ == "__main__":
    main()
