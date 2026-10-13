#!/usr/bin/env python3
"""Supplementary per-arm analyses of the full-population structural search.

Everything here re-reads artefacts that already exist on the project's OneDrive
mirror (`genome_entropy:General/missed_genes/`) -- nothing is re-searched and
nothing is re-encoded. It exists because `18_pilot_analysis.py` printed several
tables that the scientific report never used, and a few that it used with a
comparator (all shadows, any hit) that the rest of §6.1 had already replaced
(clean shadows, full-length mutual coverage). This stage recomputes those
tables both ways, so the old numbers are reproduced and the new ones sit beside
them on the same rows.

Inputs, per domain (all from the OneDrive mirror):

  wanted.tsv                       the query set: genome, input_id, orf_id,
                                   group (arm), aa_length, 3Di entropy,
                                   truncated
  best/query.<db>.<mode>.m8        best hit per query (max bits), Foldseek
                                   struct (3Di+AA) and seq (AA only) modes
  best/null.<db>.struct.m8         shuffled-3Di null
  prodigal_<dom>4.coincidence.tsv.gz
                                   per-ORF Prodigal 3'-end offset, coincidence
                                   and the shadow frame class (24_prodigal_overlap.py)
  func_<dom>_classified.tsv.gz     ranked candidates with the #97 criteria
                                   (22_rank_candidates.py) and functional class
                                   (25_functional_classes.py)
  candidate_burden_<dom>.tsv       per-genome GTDB r232 metadata (N50, CheckM2,
                                   genome size) as joined by 16_candidate_burden.py

Outputs: one small TSV per analysis under --out-dir (these are what the figure
script 35 reads and what supplementary_analyses.md quotes), a text report, and
a per-ORF parquet that is NOT committed (it is a convenience cache).

Arms. `clean_shadow` is the subset of `shadow_hi` whose cds_frame_class is
"opposite strand" or "same strand, frameshift" -- the comparator §6.1 quotes.
It is reported as its own arm everywhere; `shadow_hi` (all shadows) is kept
only so the tables printed by 18_pilot_analysis.py can be reproduced.

  34_supplementary_arms.py --domain bac --data-dir work/gtdb_extras/bac \
      --burden work/gtdb_extras/candidate_burden_bac.tsv --out-dir out/bac
"""
import argparse
import glob
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

DBS = ["afdb_swissprot", "pdb100", "cath50", "bfvd"]
EST_DBS = ["afdb_swissprot", "pdb100", "cath50"]
M8_COLS = ["query", "target", "fident", "alnlen", "qlen", "tlen", "qcov",
           "tcov", "evalue", "bits", "taxid", "taxname"]
ARMS = ["candidate", "clean_shadow", "shadow_hi", "annotated_cds",
        "intergenic_lo", "unannot_hi"]
CLEAN = {"opposite strand", "same strand, frameshift"}
FULL = 0.8

LINES = []


def emit(s=""):
    print(s)
    LINES.append(s)


def header(title):
    emit()
    emit("=" * 78)
    emit(title)
    emit("=" * 78)


def load_best(path):
    """Best hit per query. The bacterial best/ tables are already reduced by
    max bits (17_pilot_search.pbs); the archaeal ones were reduced the same way
    by this project's fetch job. Reduce again anyway -- it is idempotent."""
    if not os.path.getsize(path):
        return pd.DataFrame(columns=M8_COLS)
    df = pd.read_csv(path, sep="\t", names=M8_COLS, low_memory=False,
                     usecols=range(10))
    df = df.loc[df.groupby("query", sort=False)["bits"].idxmax()]
    return df


def arm_masks(df):
    """Arm -> boolean mask. clean_shadow is a subset of shadow_hi."""
    m = {a: (df.group == a) for a in ("candidate", "shadow_hi", "annotated_cds",
                                      "intergenic_lo", "unannot_hi")}
    m["clean_shadow"] = (df.group == "shadow_hi") & df.cds_frame_class.isin(CLEAN)
    return {a: m[a] for a in ARMS}


def mixture(c, s, a, n_c, n_s, n_a):
    """(c - s) / (a - s) with the first-order CI used in 18_pilot_analysis.py."""
    denom = a - s
    if not np.isfinite(denom) or denom <= 1e-9 or min(n_c, n_s, n_a) == 0:
        return np.nan, np.nan, np.nan
    f = (c - s) / denom
    var = (c * (1 - c) / n_c + s * (1 - s) / n_s * (1 - f) ** 2
           + a * (1 - a) / n_a * f ** 2) / denom ** 2
    se = var ** 0.5
    return f, f - 1.96 * se, f + 1.96 * se


def build(args):
    d = Path(args.data_dir)
    w = pd.read_csv(d / "wanted.tsv", sep="\t", dtype={"chunk": str,
                                                      "genome": str,
                                                      "input_id": str,
                                                      "orf_id": str})
    w["truncated"] = w.truncated.astype(str).str.lower().isin(("true", "1"))
    w["qid"] = w.genome + "|" + w.input_id + "|" + w.orf_id + "|" + w.group
    if w.qid.duplicated().any():
        sys.exit("ERROR: duplicate query ids in wanted.tsv")
    emit(f"query set: {len(w):,} ORFs; " + ", ".join(
        f"{k} {v:,}" for k, v in w.group.value_counts().items()))

    nulls = {}
    for db in DBS:
        for mode in ("struct", "seq"):
            p = d / "best" / f"query.{db}.{mode}.m8"
            if not p.exists():
                sys.exit(f"ERROR: missing {p}")
            h = load_best(p)
            pre = f"{db}_{mode}"
            h = h.rename(columns={c: f"{pre}_{c}" for c in M8_COLS
                                  if c != "query"}).rename(
                columns={"query": "qid"})
            keep = ["qid"] + [f"{pre}_{c}" for c in
                              ("fident", "alnlen", "qcov", "tcov", "evalue",
                               "bits")]
            n0 = len(h)
            h = h[keep]
            w = w.merge(h, on="qid", how="left")
            got = int(w[f"{pre}_bits"].notna().sum())
            if got != n0:
                emit(f"  WARNING: {pre}: {n0:,} best hits, {got:,} joined to "
                     "the query set")
            w[f"{pre}_hit"] = w[f"{pre}_bits"].notna()
            w[f"{pre}_full"] = ((w[f"{pre}_qcov"] >= FULL)
                                & (w[f"{pre}_tcov"] >= FULL)).fillna(False)
        p = d / "best" / f"null.{db}.struct.m8"
        nulls[db] = len(load_best(p)) if p.exists() else np.nan

    co = pd.read_csv(glob.glob(str(d / "prodigal_*4.coincidence.tsv.gz"))[0],
                     sep="\t", dtype={"genome": str, "input_id": str,
                                      "orf_id": str})
    co = co[["genome", "input_id", "orf_id", "group", "truncated_calc",
             "cds_frame_class", "offset3", "coincides"]]
    n_w = len(w)
    w = w.merge(co, on=["genome", "input_id", "orf_id", "group"], how="left")
    if len(w) != n_w:
        sys.exit("ERROR: coincidence join changed the row count")
    w["prodigal_known"] = w.coincides.notna()
    w["coincides"] = w.coincides.astype(str).str.lower().eq("true")

    md = pd.read_csv(args.burden, sep="\t", dtype={"genome": str},
                     usecols=["genome", "n50_contigs", "checkm2_completeness",
                              "checkm2_contamination", "genome_size_mb",
                              "contig_count"])
    w = w.merge(md, on="genome", how="left")

    fc = glob.glob(str(d / "func_*_classified.tsv.gz"))[0]
    cl = pd.read_csv(fc, sep="\t", low_memory=False,
                     dtype={"genome": str, "input_id": str, "orf_id": str},
                     usecols=["genome", "input_id", "orf_id", "strict_rule",
                              "c2", "c3", "c4", "c5", "c6", "c7",
                              "functional_class", "bfvd_only", "best_product",
                              "structure_only", "any_full_struct",
                              "n_db_full_struct", "truncated_calc"])
    cl = cl.rename(columns={"truncated_calc": "truncated_calc_ranked"})
    cl["group"] = "candidate"
    w = w.merge(cl, on=["genome", "input_id", "orf_id", "group"], how="left")
    if len(w) != n_w:
        sys.exit("ERROR: classified join changed the row count")

    # Union flags
    w["any_struct_hit"] = np.any([w[f"{db}_struct_hit"] for db in DBS], axis=0)
    w["any_seq_hit"] = np.any([w[f"{db}_seq_hit"] for db in DBS], axis=0)
    w["any_struct_full"] = np.any([w[f"{db}_struct_full"] for db in DBS], axis=0)
    w["any_seq_full"] = np.any([w[f"{db}_seq_full"] for db in DBS], axis=0)
    w["est_struct_full"] = np.any([w[f"{db}_struct_full"] for db in EST_DBS],
                                  axis=0)
    w["n_db_full"] = np.sum([w[f"{db}_struct_full"] for db in DBS], axis=0)
    return w, nulls


# ---------------------------------------------------------------------------
def t_overview(w, masks, out, dom, nulls):
    header("ARMS, AND THE CHECKS THAT TIE THIS BACK TO THE PUBLISHED TABLES")
    rows = []
    for a, m in masks.items():
        s = w[m]
        rows.append(dict(domain=dom, arm=a, n=int(m.sum()),
                         truncated=s.truncated.mean(),
                         any_struct_hit=s.any_struct_hit.mean(),
                         any_seq_hit=s.any_seq_hit.mean(),
                         any_struct_full=s.any_struct_full.mean(),
                         est_struct_full=s.est_struct_full.mean(),
                         sp_struct_full=s.afdb_swissprot_struct_full.mean(),
                         prodigal=(s.coincides[s.prodigal_known].mean()
                                   if s.prodigal_known.any() else np.nan),
                         median_aa=s.aa_length.median(),
                         median_3di=s.three_di_entropy.median()))
    t = pd.DataFrame(rows)
    t.to_csv(out / "arm_overview.tsv", sep="\t", index=False)
    emit(t.to_string(index=False, float_format=lambda x: f"{x:.4f}"))
    emit(f"  shuffled-3Di null best hits (struct): "
         + ", ".join(f"{k} {v:,}" for k, v in nulls.items()))
    # Cross-check truncation flags
    both = w.truncated_calc.notna()
    agree = (w.loc[both, "truncated"]
             == w.loc[both, "truncated_calc"].astype(str).str.lower().eq("true")).mean()
    emit(f"  wanted.truncated vs context truncated_calc agreement: {agree*100:.3f}% "
         f"over {int(both.sum()):,} ORFs")


def t_quintiles(w, masks, out, dom):
    header("1. STRUCTURAL HIT RATE BY 3Di ENTROPY BIN")
    rows = []
    # (a) the published table: quintiles over the whole query set pooled
    try:
        pooled = pd.qcut(w.three_di_entropy, 5,
                         labels=["Q1", "Q2", "Q3", "Q4", "Q5"])
    except ValueError:
        pooled = None
    # (b) quintiles of the CANDIDATE distribution, applied to every arm. The
    # pooled bins in (a) are dominated by intergenic_lo and annotated_cds, so
    # all candidates fall in Q1-Q5 by a different rule than the controls and
    # the first bin starts at 0 bits, below the 2.5-bit candidate threshold.
    cand = w.loc[masks["candidate"], "three_di_entropy"]
    edges = np.quantile(cand, [0, .2, .4, .6, .8, 1.0])
    edges[0], edges[-1] = -np.inf, np.inf
    cbin = pd.cut(w.three_di_entropy, edges,
                  labels=["Q1", "Q2", "Q3", "Q4", "Q5"])
    qe = np.quantile(cand, [0, .2, .4, .6, .8, 1.0])
    for scheme, b in (("pooled_query_set", pooled), ("candidate_quintiles", cbin)):
        if b is None:
            continue
        for q in ["Q1", "Q2", "Q3", "Q4", "Q5"]:
            inq = (b == q)
            lo, hi = (w.loc[inq, "three_di_entropy"].min(),
                      w.loc[inq, "three_di_entropy"].max())
            if scheme == "candidate_quintiles":
                i = int(q[1]) - 1
                lo, hi = qe[i], qe[i + 1]
            rec = {}
            for a, m in masks.items():
                s = w[m & inq]
                rows.append(dict(domain=dom, scheme=scheme, bin=q, lo=lo, hi=hi,
                                 arm=a, n=len(s),
                                 any_struct_hit=s.any_struct_hit.mean() if len(s) else np.nan,
                                 est_struct_full=s.est_struct_full.mean() if len(s) else np.nan,
                                 sp_struct_full=s.afdb_swissprot_struct_full.mean() if len(s) else np.nan,
                                 prodigal=(s.coincides[s.prodigal_known].mean()
                                           if s.prodigal_known.any() else np.nan)))
                rec[a] = (len(s), s.afdb_swissprot_struct_full.mean() if len(s) else np.nan)
            # mixture share within bin, Swiss-Prot full-length, clean comparator
            (nc, c), (ns, sh), (na, an) = rec["candidate"], rec["clean_shadow"], rec["annotated_cds"]
            f, flo, fhi = mixture(c, sh, an, nc, ns, na)
            rows.append(dict(domain=dom, scheme=scheme, bin=q, lo=lo, hi=hi,
                             arm="share_sp_clean", n=nc, any_struct_hit=np.nan,
                             est_struct_full=f, sp_struct_full=flo,
                             prodigal=fhi))
    t = pd.DataFrame(rows)
    t.to_csv(out / "hit_by_3di_bin.tsv", sep="\t", index=False)
    for scheme in t.scheme.unique():
        emit(f"\n  scheme: {scheme}")
        s = t[(t.scheme == scheme) & (t.arm != "share_sp_clean")]
        for metric in ("any_struct_hit", "est_struct_full", "prodigal"):
            p = s.pivot_table(index=["bin", "lo", "hi"], columns="arm",
                              values=metric, observed=True).reindex(columns=ARMS) * 100
            emit(f"  -- {metric} (%)")
            emit(p.round(1).to_string())
        p = s.pivot_table(index="bin", columns="arm", values="n").reindex(columns=ARMS)
        emit("  -- n")
        emit(p.to_string())
        sh = t[(t.scheme == scheme) & (t.arm == "share_sp_clean")]
        emit("  -- implied real-gene share in bin (Swiss-Prot full-length, clean comparator)")
        for r in sh.itertuples():
            emit(f"     {r.bin}: {r.est_struct_full*100:6.1f}% "
                 f"({r.sp_struct_full*100:.1f} - {r.prodigal*100:.1f})")


QUAL = ["fident", "qcov", "tcov", "bits", "nlog10e", "alnlen"]


def t_quality(w, masks, out, dom):
    header("2. BEST-HIT QUALITY BY ARM (struct mode; queries with a hit)")
    rows, samp = [], []
    rng = np.random.default_rng(0)
    for db in DBS:
        for mode in ("struct", "seq"):
            pre = f"{db}_{mode}"
            for a, m in masks.items():
                for subset in ("all_hits", "full_length"):
                    sel = m & w[f"{pre}_hit"]
                    if subset == "full_length":
                        sel &= w[f"{pre}_full"]
                    s = w.loc[sel]
                    if not len(s):
                        continue
                    vals = {"fident": s[f"{pre}_fident"], "qcov": s[f"{pre}_qcov"],
                            "tcov": s[f"{pre}_tcov"], "bits": s[f"{pre}_bits"],
                            "nlog10e": -np.log10(s[f"{pre}_evalue"].clip(lower=1e-300)),
                            "alnlen": s[f"{pre}_alnlen"]}
                    for k, v in vals.items():
                        q = np.nanquantile(v, [.1, .25, .5, .75, .9])
                        rows.append(dict(domain=dom, db=db, mode=mode, arm=a,
                                         subset=subset, metric=k, n=len(s),
                                         mean=float(np.nanmean(v)),
                                         q10=q[0], q25=q[1], median=q[2],
                                         q75=q[3], q90=q[4]))
                    if mode == "struct" and subset == "all_hits":
                        k = min(len(s), 20000)
                        idx = rng.choice(len(s), k, replace=False)
                        ss = pd.DataFrame({m_: np.asarray(v_)[idx]
                                           for m_, v_ in vals.items()})
                        ss["db"], ss["arm"], ss["domain"] = db, a, dom
                        ss["full"] = np.asarray(s[f"{pre}_full"])[idx]
                        samp.append(ss)
    t = pd.DataFrame(rows)
    t.to_csv(out / "hit_quality.tsv", sep="\t", index=False)
    pd.concat(samp).to_csv(out / "hit_quality_sample.tsv.gz", sep="\t",
                           index=False, compression="gzip")
    s = t[(t["mode"] == "struct") & (t.subset == "all_hits")]
    p = s.pivot_table(index=["db", "arm"], columns="metric", values="median")
    emit(p[["fident", "qcov", "tcov", "bits", "nlog10e"]].round(3).to_string())
    emit("\n  full-length hits only (median fident)")
    s = t[(t["mode"] == "struct") & (t.subset == "full_length") & (t.metric == "fident")]
    emit(s.pivot_table(index="db", columns="arm", values="median")[
        [a for a in ARMS if a in s.arm.unique()]].round(3).to_string())


def t_evidence(w, masks, out, dom):
    header("3. STRUCTURE-ONLY vs SEQUENCE EVIDENCE")
    rows = []
    for level in ("hit", "full"):
        for db in DBS + ["union"]:
            if db == "union":
                st = w[f"any_struct_{'hit' if level == 'hit' else 'full'}"]
                sq = w[f"any_seq_{'hit' if level == 'hit' else 'full'}"]
            else:
                st, sq = w[f"{db}_struct_{level}"], w[f"{db}_seq_{level}"]
            for a, m in masks.items():
                n = int(m.sum())
                rows.append(dict(domain=dom, level=level, db=db, arm=a, n=n,
                                 struct=(st & m).sum() / n,
                                 seq=(sq & m).sum() / n,
                                 struct_only=(st & ~sq & m).sum() / n,
                                 seq_only=(~st & sq & m).sum() / n,
                                 both=(st & sq & m).sum() / n,
                                 neither=(~st & ~sq & m).sum() / n))
    t = pd.DataFrame(rows)
    t.to_csv(out / "evidence_classes.tsv", sep="\t", index=False)
    for level in ("hit", "full"):
        emit(f"\n  level = {level} (union over four databases), % of arm")
        s = t[(t.level == level) & (t.db == "union")].set_index("arm")
        emit((s[["struct", "seq", "struct_only", "seq_only", "both",
                 "neither"]] * 100).round(1).to_string())
    # structure-only as share of structural evidence
    emit("\n  structure-only as a share of ORFs WITH full-length structural evidence")
    s = t[(t.level == "full") & (t.db == "union")].set_index("arm")
    emit((s.struct_only / s.struct * 100).round(1).to_string())
    # Excess over clean shadow, full-length union
    c = s.loc["candidate"]; cs = s.loc["clean_shadow"]
    for k in ("struct_only", "both"):
        emit(f"  excess {k} (candidate - clean shadow): "
             f"{(c[k]-cs[k])*100:.2f} pts = {(c[k]-cs[k])*c.n:,.0f} candidates")


def t_truncation(w, masks, out, dom):
    header("5. CONTIG-END TRUNCATION BY ARM")
    rows = []
    for a, m in masks.items():
        for tr in (False, True):
            s = w[m & (w.truncated == tr)]
            rows.append(dict(domain=dom, arm=a, truncated=tr, n=len(s),
                             share=len(s) / max(int(m.sum()), 1),
                             any_struct_hit=s.any_struct_hit.mean() if len(s) else np.nan,
                             any_struct_full=s.any_struct_full.mean() if len(s) else np.nan,
                             prodigal=(s.coincides[s.prodigal_known].mean()
                                       if s.prodigal_known.any() else np.nan)))
    # shadow truncation by frame class
    sh = w[masks["shadow_hi"]]
    for k, g in sh.groupby(sh.cds_frame_class.fillna("(no overlap)")):
        rows.append(dict(domain=dom, arm=f"shadow:{k}", truncated=True,
                         n=int(g.truncated.sum()), share=g.truncated.mean(),
                         any_struct_hit=np.nan, any_struct_full=np.nan,
                         prodigal=np.nan))
    t = pd.DataFrame(rows)
    t.to_csv(out / "truncation.tsv", sep="\t", index=False)
    emit(t.to_string(index=False, float_format=lambda x: f"{x:.4f}"))
    # The C2 row of the ladder, recomputed
    nc = int(masks["candidate"].sum()); ns = int(masks["clean_shadow"].sum())
    c2c = int((masks["candidate"] & ~w.truncated).sum())
    c2s = int((masks["clean_shadow"] & ~w.truncated).sum())
    na = int(masks["shadow_hi"].sum())
    c2a = int((masks["shadow_hi"] & ~w.truncated).sum())
    emit(f"\n  ladder C2 recomputed: candidates {c2c:,} of {nc:,} "
         f"({c2c/nc*100:.2f}% complete)")
    emit(f"    clean shadows {c2s:,} of {ns:,} ({c2s/ns*100:.2f}%), scaled "
         f"{c2s*nc/ns:,.0f} -> excess {c2c - c2s*nc/ns:,.0f}")
    emit(f"    all shadows   {c2a:,} of {na:,} ({c2a/na*100:.2f}%), scaled "
         f"{c2a*nc/na:,.0f} -> excess {c2c - c2a*nc/na:,.0f}")


def t_assembly(w, masks, out, dom):
    header("6. ASSEMBLY-QUALITY STRATIFICATION (N50, CheckM2 completeness)")
    rows = []
    has = w.n50_contigs.notna()
    for cov in ("n50_contigs", "checkm2_completeness"):
        for scheme in ("pooled_query_rows", "candidate_genomes"):
            if scheme == "pooled_query_rows":
                # exactly what 18_pilot_analysis.py did
                b = pd.Series(pd.NA, index=w.index, dtype=object)
                b[has] = pd.qcut(w.loc[has, cov], 4, labels=["Q1", "Q2", "Q3", "Q4"],
                                 duplicates="drop").astype(str)
            else:
                # quartiles of the per-genome value over genomes carrying a candidate
                g = w.loc[masks["candidate"], ["genome", cov]].drop_duplicates("genome")
                e = np.quantile(g[cov].dropna(), [0, .25, .5, .75, 1.0])
                e[0], e[-1] = -np.inf, np.inf
                e = np.unique(e)
                labs = [f"Q{i+1}" for i in range(len(e) - 1)]
                b = pd.cut(w[cov], e, labels=labs).astype(object)
            for q in sorted(pd.Series(b).dropna().unique()):
                inq = (b == q).fillna(False).astype(bool)
                lo, hi = w.loc[inq, cov].min(), w.loc[inq, cov].max()
                rec = {}
                for a, m in masks.items():
                    s = w[m & inq]
                    n = len(s)
                    rec[a] = (n, s.afdb_swissprot_struct_full.mean() if n else np.nan,
                              s.est_struct_full.mean() if n else np.nan)
                    rows.append(dict(domain=dom, covariate=cov, scheme=scheme,
                                     bin=q, lo=lo, hi=hi, arm=a, n=n,
                                     any_struct_hit=s.any_struct_hit.mean() if n else np.nan,
                                     est_struct_full=s.est_struct_full.mean() if n else np.nan,
                                     sp_struct_full=s.afdb_swissprot_struct_full.mean() if n else np.nan,
                                     prodigal=(s.coincides[s.prodigal_known].mean()
                                               if s.prodigal_known.any() else np.nan)))
                (nc, c, _), (ns, sh, _), (na, an, _) = (rec["candidate"],
                                                        rec["clean_shadow"],
                                                        rec["annotated_cds"])
                f, flo, fhi = mixture(c, sh, an, nc, ns, na)
                rows.append(dict(domain=dom, covariate=cov, scheme=scheme, bin=q,
                                 lo=lo, hi=hi, arm="share_sp_clean", n=nc,
                                 any_struct_hit=np.nan, est_struct_full=f,
                                 sp_struct_full=flo, prodigal=fhi))
    t = pd.DataFrame(rows)
    t.to_csv(out / "assembly_strata.tsv", sep="\t", index=False)
    for (cov, scheme), s in t.groupby(["covariate", "scheme"], sort=False):
        emit(f"\n  {cov} -- quartiles over {scheme}")
        s2 = s[s.arm != "share_sp_clean"]
        for metric in ("any_struct_hit", "sp_struct_full", "prodigal"):
            p = s2.pivot_table(index=["bin", "lo", "hi"], columns="arm",
                               values=metric).reindex(columns=[a for a in ARMS if a != "unannot_hi"]) * 100
            emit(f"  -- {metric} (%)")
            emit(p.round(1).to_string())
        sh = s[s.arm == "share_sp_clean"]
        emit("  -- implied real-gene share (Swiss-Prot full-length, clean comparator)")
        for r in sh.itertuples():
            emit(f"     {r.bin}: {r.est_struct_full*100:6.1f}% "
                 f"({r.sp_struct_full*100:.1f} - {r.prodigal*100:.1f})")


def offset_class(o):
    o = pd.to_numeric(o, errors="coerce")
    big = np.iinfo(np.int64).max
    out = pd.Series("no same-strand Prodigal gene on contig", index=o.index)
    k = o.notna() & (o != big)
    out[k & (o == 0)] = "0 (exact 3' match)"
    out[k & o.isin([-1, -2, 1, 2])] = "±1-2 bp (other frame)"
    out[k & (o != 0) & (o.abs() > 2) & (o % 3 == 0)] = "other, same frame (multiple of 3)"
    out[k & (o.abs() > 2) & (o % 3 != 0)] = "other, different frame"
    return out


def t_prodigal(w, masks, out, dom):
    header("8. PRODIGAL AGREEMENT BY FRAME CLASS AND 3'-OFFSET CLASS")
    k = w.prodigal_known
    w = w.copy()
    w["offset_class"] = offset_class(w.offset3)
    w.loc[~k, "offset_class"] = np.nan
    rows = []
    groups = {a: m for a, m in masks.items() if a != "unannot_hi"}
    sh = masks["shadow_hi"]
    for fc in w.loc[sh, "cds_frame_class"].dropna().unique():
        groups[f"shadow:{fc}"] = sh & (w.cds_frame_class == fc)
    for a, m in groups.items():
        s = w[m & k]
        vc = s.offset_class.value_counts()
        for cls, n in vc.items():
            rows.append(dict(domain=dom, group=a, n_group=len(s),
                             offset_class=cls, n=int(n), share=n / len(s)))
    t = pd.DataFrame(rows)
    t.to_csv(out / "prodigal_offsets.tsv", sep="\t", index=False)
    p = t.pivot_table(index="group", columns="offset_class", values="share") * 100
    emit(p.round(2).to_string())
    # exact offsets near 0
    emit("\n  exact offset histogram -3..+3, % of group")
    rows2 = []
    for a, m in groups.items():
        s = w.loc[m & k, "offset3"]
        rows2.append({"group": a, **{str(o): (s == o).mean() * 100
                                     for o in range(-3, 4)}})
    emit(pd.DataFrame(rows2).set_index("group").round(2).to_string())


def t_candidate_subsets(w, masks, out, dom):
    header("10. PRODIGAL AGREEMENT OF CANDIDATE SUBSETS (incl. C7 strict set)")
    c = w[masks["candidate"]].copy()
    rows = []
    for lab in ("c2", "c3", "c4", "c5", "c6", "c7"):
        c[lab] = c[lab].astype(str).str.lower().eq("true")
    subsets = {"all candidates": pd.Series(True, index=c.index),
               "C2 not truncated": c.c2, "C3": c.c3, "C4": c.c4, "C5": c.c5,
               "C6 structure-only": c.c6, "C7 strict": c.c7,
               "C3 and not C7": c.c3 & ~c.c7,
               "no struct hit at all": ~c.any_struct_hit,
               "BFVD-only (C3)": c.c3 & c.bfvd_only.astype(str).str.lower().eq("true")}
    cs = w[masks["clean_shadow"]].copy()
    cs_sub = {"clean shadow (all)": pd.Series(True, index=cs.index),
              "clean shadow, C3-equivalent": (~cs.truncated & cs.any_struct_full),
              "clean shadow, C7-like (>=2 db, >=1 est, >=100 aa)":
                  (~cs.truncated & (cs.n_db_full >= 2) & cs.est_struct_full
                   & (cs.aa_length >= 100))}
    for name, m in list(subsets.items()):
        s = c[m & c.prodigal_known]
        rows.append(dict(domain=dom, subset=name, n=int(m.sum()),
                         n_prodigal_known=len(s), coincide=int(s.coincides.sum()),
                         rate=s.coincides.mean()))
    for name, m in cs_sub.items():
        s = cs[m & cs.prodigal_known]
        rows.append(dict(domain=dom, subset=name, n=int(m.sum()),
                         n_prodigal_known=len(s), coincide=int(s.coincides.sum()),
                         rate=s.coincides.mean()))
    t = pd.DataFrame(rows)
    t.to_csv(out / "prodigal_candidate_subsets.tsv", sep="\t", index=False)
    emit(t.to_string(index=False, float_format=lambda x: f"{x:.4f}"))
    # consistency of our C3 with the ranked table's C3
    mine = ~c.truncated & c.any_struct_full
    emit(f"\n  C3 recomputed here {int(mine.sum()):,} vs ranked table "
         f"{int(c.c3.sum()):,}; agreement {(mine == c.c3).mean()*100:.3f}%")


def t_function(w, masks, out, dom):
    header("7. FUNCTIONAL COMPOSITION OF ALL CANDIDATES, C3 AND C7")
    c = w[masks["candidate"]].copy()
    for lab in ("c3", "c7"):
        c[lab] = c[lab].astype(str).str.lower().eq("true")
    c["bfvd_only"] = c.bfvd_only.astype(str).str.lower().eq("true")
    # Why a candidate is unclassified: no structural hit at all, or its best
    # Swiss-Prot/PDB hit is absent (only CATH/BFVD hit), or has product text.
    sp_pdb = c.afdb_swissprot_struct_hit | c.pdb100_struct_hit
    c["label_source"] = np.select(
        [~c.any_struct_hit, ~sp_pdb],
        ["no structural hit", "hit only in CATH50/BFVD"],
        default="Swiss-Prot/PDB hit (any coverage)")
    # Full-length-labelled: the class only counts if the Swiss-Prot or PDB hit
    # that supplied the product is itself full-length.
    c["labelled_full"] = c.afdb_swissprot_struct_full | c.pdb100_struct_full
    rows = []
    for name, m in (("all candidates", pd.Series(True, index=c.index)),
                    ("C3", c.c3), ("C7", c.c7)):
        s = c[m]
        for cls, n in s.functional_class.value_counts().items():
            rows.append(dict(domain=dom, set=name, n_set=len(s), kind="class",
                             key=cls, n=int(n), share=n / len(s),
                             share_full_label=((s.functional_class == cls)
                                               & s.labelled_full).mean()))
        for k2, n in s.label_source.value_counts().items():
            rows.append(dict(domain=dom, set=name, n_set=len(s), kind="label_source",
                             key=k2, n=int(n), share=n / len(s),
                             share_full_label=np.nan))
        rows.append(dict(domain=dom, set=name, n_set=len(s), kind="flag",
                         key="BFVD-only full-length", n=int(s.bfvd_only.sum()),
                         share=s.bfvd_only.mean(), share_full_label=np.nan))
        rows.append(dict(domain=dom, set=name, n_set=len(s), kind="flag",
                         key="product from full-length SP/PDB hit",
                         n=int(s.labelled_full.sum()), share=s.labelled_full.mean(),
                         share_full_label=np.nan))
        rows.append(dict(domain=dom, set=name, n_set=len(s), kind="flag",
                         key="unclassified & BFVD-only",
                         n=int((s.bfvd_only & (s.functional_class == "unclassified")).sum()),
                         share=(s.bfvd_only & (s.functional_class == "unclassified")).mean(),
                         share_full_label=np.nan))
    t = pd.DataFrame(rows)
    t.to_csv(out / "functional_composition_sets.tsv", sep="\t", index=False)
    for name in ("all candidates", "C3", "C7"):
        s = t[t.set == name]
        emit(f"\n  {name} (n = {int(s.n_set.iloc[0]):,})")
        emit(s[["kind", "key", "n", "share", "share_full_label"]].to_string(
            index=False, float_format=lambda x: f"{x*100:.1f}%"))


def t_unannot(w, masks, out, dom):
    header("4. THE unannot_hi ARM")
    rows = []
    for db in DBS + ["union_any", "union_est"]:
        if db == "union_any":
            col = "any_struct_full"
        elif db == "union_est":
            col = "est_struct_full"
        else:
            col = f"{db}_struct_full"
        r = {a: w.loc[m, col].mean() for a, m in masks.items()}
        n = {a: int(m.sum()) for a, m in masks.items()}
        # Where unannot_hi sits on the clean-shadow -> annotated axis
        f, lo, hi = mixture(r["unannot_hi"], r["clean_shadow"], r["annotated_cds"],
                            n["unannot_hi"], n["clean_shadow"], n["annotated_cds"])
        fc, lc, hc = mixture(r["candidate"], r["clean_shadow"], r["annotated_cds"],
                             n["candidate"], n["clean_shadow"], n["annotated_cds"])
        rows.append(dict(domain=dom, db=db, **{f"full_{a}": r[a] for a in ARMS},
                         unannot_share=f, unannot_lo=lo, unannot_hi_ci=hi,
                         cand_share=fc, cand_lo=lc, cand_hi=hc))
    hit = {a: w.loc[m, "any_struct_hit"].mean() for a, m in masks.items()}
    t = pd.DataFrame(rows)
    t.to_csv(out / "unannot_hi.tsv", sep="\t", index=False)
    emit("  any-hit union: " + ", ".join(f"{a} {v*100:.1f}%" for a, v in hit.items()))
    emit(t.to_string(index=False, float_format=lambda x: f"{x*100:.1f}"))
    u = w[masks["unannot_hi"]]; c = w[masks["candidate"]]
    emit(f"\n  matching check: median aa {u.aa_length.median():.0f} vs candidate "
         f"{c.aa_length.median():.0f}; median 3Di {u.three_di_entropy.median():.3f} "
         f"vs {c.three_di_entropy.median():.3f}; genomes {u.genome.nunique():,}")


def t_standardized(w, masks, out, dom, n_boot=1000):
    """Does dropping same-frame shadows break the length/3Di matching?

    shadow_hi was drawn 1:1 to candidates on length and 3Di entropy. The clean
    comparator keeps only the antisense and frameshift members of that arm, so
    it is a non-random SUBSET of a matched arm and is no longer matched. Both
    full-length coverage and Prodigal agreement rise steeply with 3Di entropy
    and with length, so a residual mismatch moves the mixture estimate.

    Direct standardisation: strata are candidate quintiles of 3Di entropy x
    candidate quintiles of aa_length (25 cells). Each comparator's rate is
    re-weighted to the candidate distribution over those cells. Two versions:
      shadow_std      : only the clean-shadow rate re-weighted
      both_std        : clean-shadow AND annotated rates re-weighted (annotated_cds
                        was matched on length only, never on 3Di)
    CIs from a parametric bootstrap of per-cell binomials.
    """
    header("STANDARDISED MIXTURE: clean comparator re-weighted to candidate "
           "length x 3Di")
    c_mask = masks["candidate"]
    cand = w[c_mask]
    e3 = np.quantile(cand.three_di_entropy, [0, .2, .4, .6, .8, 1.0])
    eL = np.unique(np.quantile(cand.aa_length, [0, .2, .4, .6, .8, 1.0]))
    e3[0], e3[-1] = -np.inf, np.inf
    eL[0], eL[-1] = -np.inf, np.inf
    cell = (pd.cut(w.three_di_entropy, e3, labels=False).astype("Int64").astype(str)
            + "_" + pd.cut(w.aa_length, eL, labels=False).astype("Int64").astype(str))
    rng = np.random.default_rng(1)
    metrics = {"afdb_swissprot": "afdb_swissprot_struct_full",
               "pdb100": "pdb100_struct_full", "cath50": "cath50_struct_full",
               "est_union": "est_struct_full", "prodigal": "coincides"}
    rows = []
    wc = cell[c_mask].value_counts(normalize=True)
    for name, col in metrics.items():
        per = {}
        for a in ("candidate", "clean_shadow", "annotated_cds"):
            m = masks[a]
            if name == "prodigal":
                m = m & w.prodigal_known
            g = w.loc[m].groupby(cell[m])[col].agg(["sum", "count"])
            per[a] = g.reindex(wc.index).fillna(0)
        crude = {a: per[a]["sum"].sum() / per[a]["count"].sum() for a in per}
        ok = (per["clean_shadow"]["count"] > 0) & (per["annotated_cds"]["count"] > 0)
        cover = wc[ok].sum()
        ww = wc[ok] / cover

        def std(k, s_, n_):
            r = (s_ / n_.where(n_ > 0))[ok]
            return float((r * ww).sum())

        c = crude["candidate"]
        s_std = std("s", per["clean_shadow"]["sum"], per["clean_shadow"]["count"])
        a_std = std("a", per["annotated_cds"]["sum"], per["annotated_cds"]["count"])
        f_crude = (c - crude["clean_shadow"]) / (crude["annotated_cds"] - crude["clean_shadow"])
        f_s = (c - s_std) / (crude["annotated_cds"] - s_std)
        f_b = (c - s_std) / (a_std - s_std)
        # parametric bootstrap
        boots = {"crude": [], "shadow_std": [], "both_std": []}
        nC = int(per["candidate"]["count"].sum())
        for _ in range(n_boot):
            cb = rng.binomial(nC, c) / nC
            sim = {}
            for a in ("clean_shadow", "annotated_cds"):
                n_ = per[a]["count"].astype(int).to_numpy()
                p_ = (per[a]["sum"] / per[a]["count"].where(per[a]["count"] > 0)).fillna(0).to_numpy()
                sim[a] = pd.Series(rng.binomial(n_, p_), index=per[a].index)
            sc = sim["clean_shadow"].sum() / per["clean_shadow"]["count"].sum()
            ac = sim["annotated_cds"].sum() / per["annotated_cds"]["count"].sum()
            ss = std("s", sim["clean_shadow"], per["clean_shadow"]["count"])
            aa = std("a", sim["annotated_cds"], per["annotated_cds"]["count"])
            boots["crude"].append((cb - sc) / (ac - sc))
            boots["shadow_std"].append((cb - ss) / (ac - ss))
            boots["both_std"].append((cb - ss) / (aa - ss))
        ci = {k: np.nanquantile(v, [.025, .975]) for k, v in boots.items()}
        rows.append(dict(domain=dom, readout=name, candidate=c,
                         clean_shadow_crude=crude["clean_shadow"],
                         clean_shadow_std=s_std,
                         annotated_crude=crude["annotated_cds"], annotated_std=a_std,
                         share_crude=f_crude, share_crude_lo=ci["crude"][0],
                         share_crude_hi=ci["crude"][1],
                         share_shadow_std=f_s, share_shadow_std_lo=ci["shadow_std"][0],
                         share_shadow_std_hi=ci["shadow_std"][1],
                         share_both_std=f_b, share_both_std_lo=ci["both_std"][0],
                         share_both_std_hi=ci["both_std"][1],
                         candidate_weight_covered=cover,
                         n_candidate=nC))
    t = pd.DataFrame(rows)
    t.to_csv(out / "standardized_mixture.tsv", sep="\t", index=False)
    emit(t.to_string(index=False, float_format=lambda x: f"{x:.4f}"))

    # The direct-evidence ladder against the same re-weighted comparator.
    # C7's "interpretable product" clause needs product text, which exists for
    # candidates only; for shadows it is approximated by "full-length Swiss-Prot
    # or PDB hit" (both databases carry text for essentially every entry) and
    # the approximation is checked against the published shadow count.
    best_e = pd.concat([w[f"{db}_struct_evalue"] for db in DBS], axis=1).min(axis=1)
    nt = ~w.truncated
    crit = {
        "C3": nt & w.any_struct_full,
        "C4": nt & (w.n_db_full >= 2),
        "C5": nt & w.est_struct_full,
        "C6": nt & w.any_struct_full & ~w.any_seq_full,
        "C7~": (nt & w.est_struct_full & (w.n_db_full >= 2)
                & (w.afdb_swissprot_struct_full | w.pdb100_struct_full)
                & (w.aa_length >= 100) & (best_e < 1e-10)),
    }
    lrows = []
    cs_m = masks["clean_shadow"]
    nC, nS = int(c_mask.sum()), int(cs_m.sum())
    for k, m in crit.items():
        cand_n = int((m & c_mask).sum())
        if k == "C7~":
            cand_n_pub = int(w.loc[c_mask, "c7"].astype(str).str.lower().eq("true").sum())
        else:
            cand_n_pub = cand_n
        sh_n = int((m & cs_m).sum())
        g = m[cs_m].groupby(cell[cs_m]).mean().reindex(wc.index)
        s_std = float((g * wc).sum() / wc[g.notna()].sum())
        lrows.append(dict(domain=dom, criterion=k, candidates=cand_n_pub,
                          candidates_recomputed=cand_n, clean_shadows=sh_n,
                          excess_crude=cand_n_pub - sh_n * nC / nS,
                          clean_shadow_rate_std=s_std,
                          excess_std=cand_n_pub - s_std * nC))
    lt = pd.DataFrame(lrows)
    lt.to_csv(out / "standardized_ladder.tsv", sep="\t", index=False)

    # Which matching variable carries the correction: 3Di alone, length alone,
    # or both. Point estimates only (Swiss-Prot full-length, est-union, Prodigal).
    q3 = pd.cut(w.three_di_entropy, e3, labels=False).astype("Int64").astype(str)
    qL = pd.cut(w.aa_length, eL, labels=False).astype("Int64").astype(str)
    q3f = pd.cut(w.three_di_entropy,
                 np.r_[-np.inf, np.quantile(cand.three_di_entropy,
                                            np.linspace(.1, .9, 9)), np.inf],
                 labels=False).astype("Int64").astype(str)
    qLf = pd.cut(w.aa_length,
                 np.unique(np.r_[-np.inf, np.quantile(cand.aa_length,
                                                      np.linspace(.1, .9, 9)), np.inf]),
                 labels=False).astype("Int64").astype(str)
    schemes = {"none": pd.Series("x", index=w.index), "3Di quintile": q3,
               "length quintile": qL, "3Di x length quintiles": q3 + "_" + qL,
               "3Di x length deciles": q3f + "_" + qLf}
    srows = []
    for sname, cl_ in schemes.items():
        wc_ = cl_[c_mask].value_counts(normalize=True)
        for name, col in (("afdb_swissprot", "afdb_swissprot_struct_full"),
                          ("est_union", "est_struct_full"),
                          ("prodigal", "coincides")):
            def rate(a):
                m = masks[a] & (w.prodigal_known if name == "prodigal" else True)
                g = w.loc[m, col].groupby(cl_[m]).mean().reindex(wc_.index)
                ok = g.notna()
                return float((g[ok] * wc_[ok]).sum() / wc_[ok].sum()), float(wc_[ok].sum())
            c = float(w.loc[c_mask & (w.prodigal_known if name == "prodigal" else True), col].mean())
            s_, cov_s = rate("clean_shadow")
            a_, cov_a = rate("annotated_cds")
            srows.append(dict(domain=dom, scheme=sname, readout=name,
                              n_cells=len(wc_), candidate=c, clean_shadow_std=s_,
                              annotated_std=a_, share=(c - s_) / (a_ - s_),
                              weight_covered=min(cov_s, cov_a)))
    st = pd.DataFrame(srows)
    st.to_csv(out / "standardized_sensitivity.tsv", sep="\t", index=False)
    emit("\n  sensitivity of the share to the standardisation variables "
         "(both comparators re-weighted)")
    emit(st.to_string(index=False, float_format=lambda x: f"{x:.4f}"))
    emit("\n  direct-evidence ladder, crude vs length x 3Di standardised clean comparator")
    emit(lt.to_string(index=False, float_format=lambda x: f"{x:,.4f}" if abs(x) < 1 else f"{x:,.0f}"))
    cs = w[masks["clean_shadow"]]
    emit(f"\n  matching after the clean filter: candidate median aa "
         f"{cand.aa_length.median():.0f}, 3Di {cand.three_di_entropy.median():.3f}; "
         f"clean shadow median aa {cs.aa_length.median():.0f}, "
         f"3Di {cs.three_di_entropy.median():.3f}")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--domain", required=True, choices=("bac", "arc"))
    ap.add_argument("--data-dir", required=True)
    ap.add_argument("--burden", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--parquet", default=None)
    args = ap.parse_args()
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    w, nulls = build(args)
    masks = arm_masks(w)
    if args.parquet:
        keep = [c for c in w.columns if c not in ("best_product",)]
        w[keep].to_parquet(args.parquet, index=False)

    t_overview(w, masks, out, args.domain, nulls)
    t_quintiles(w, masks, out, args.domain)
    t_quality(w, masks, out, args.domain)
    t_evidence(w, masks, out, args.domain)
    t_unannot(w, masks, out, args.domain)
    t_truncation(w, masks, out, args.domain)
    t_assembly(w, masks, out, args.domain)
    t_function(w, masks, out, args.domain)
    t_prodigal(w, masks, out, args.domain)
    t_candidate_subsets(w, masks, out, args.domain)
    t_standardized(w, masks, out, args.domain)
    (out / "report.txt").write_text("\n".join(LINES) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
