#!/usr/bin/env python3
"""Mixture shares, ladder floors and Prodigal share against a given comparator.

One table per domain with every §6.1 number that depends on the clean-shadow
comparator, computed the same way 18_pilot_analysis.py (coverage_table /
shadow_class_coverage), 22_rank_candidates.py (tiers, imported) and
24_prodigal_overlap.py do, but with the comparator WEIGHTED by how often each
ORF was drawn. 40_match_clean_comparator.py falls back to drawing with
replacement in the sparse high-3Di tail, and 18/22/24 count each query ID once,
so they slightly under-weight that tail; this script does not.

Run it twice per domain: once with --comparator pointing at the published
clean arm (to reproduce v6 / ranked5 / prodigal_*4 -- the check that this code
computes the published quantities) and once with the new matched arm.

  41_comparator_summary.py --domain bac --label v7 \
      --wanted <...>/wanted_v7.tsv --comparator <...>/comparator.tsv.gz \
      --context <...>/context_v7 --search-dirs <old best>,<new best> \
      --prodigal <...>/prodigal_v7_bac.coincidence.tsv.gz --out <...>.tsv
"""
import argparse
import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
_spec = importlib.util.spec_from_file_location("rank", HERE / "22_rank_candidates.py")
m22 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(m22)

KEY = ["genome", "input_id", "orf_id"]
STR = {"genome": str, "input_id": str, "orf_id": str, "chunk": str,
       "group": str}
CLEAN = {"opposite strand", "same strand, frameshift"}


def qid(d, group):
    return d.genome + "|" + d.input_id + "|" + d.orf_id + "|" + group


def share(c, n_c, s, n_s, a, n_a):
    """(c - s) / (a - s) with the first-order CI 18_pilot_analysis.py uses."""
    denom = a - s
    if denom <= 1e-9:
        return np.nan, np.nan, np.nan
    f = (c - s) / denom
    var = (c * (1 - c) / n_c + s * (1 - s) / n_s * (1 - f) ** 2
           + a * (1 - a) / n_a * f ** 2) / denom ** 2
    se = var ** 0.5
    return f, f - 1.96 * se, f + 1.96 * se


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--domain", required=True)
    ap.add_argument("--label", required=True)
    ap.add_argument("--wanted", required=True,
                    help="wanted list carrying the candidate and annotated_cds "
                         "arms (sizes and IDs)")
    ap.add_argument("--comparator", required=True,
                    help="TSV with genome/input_id/orf_id, one row per DRAW "
                         "(duplicates = weight); or 'context-clean' to take "
                         "the clean shadow_hi rows of --context, weight 1")
    ap.add_argument("--context", required=True)
    ap.add_argument("--search-dirs", required=True)
    ap.add_argument("--descriptions",
                    default=f"{m22.GD}/foldseek_db/descriptions")
    ap.add_argument("--prodigal", default=None,
                    help="24_prodigal_overlap.py per-ORF coincidence table")
    ap.add_argument("--exclude-replacement", action="store_true",
                    help="sensitivity: drop candidates matched WITH "
                         "replacement, and their draws, from every arm")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    wanted = pd.read_csv(args.wanted, sep="\t", dtype=STR)
    ctx = m22.load_context(args.context, {"candidate", "shadow_hi"})
    if args.comparator == "context-clean":
        sh = ctx[(ctx.group == "shadow_hi")
                 & ctx.cds_frame_class.isin(CLEAN)]
        draws = sh[KEY].copy()
    else:
        cols = KEY + (["match_level", "cand_genome", "cand_input_id",
                       "cand_orf_id"] if args.exclude_replacement else [])
        draws = pd.read_csv(args.comparator, sep="\t", dtype=STR,
                            usecols=cols)
    dropped = None
    if args.exclude_replacement:
        rep = draws.match_level == "with replacement"
        dropped = set(draws.loc[rep, "cand_genome"] + "|"
                      + draws.loc[rep, "cand_input_id"] + "|"
                      + draws.loc[rep, "cand_orf_id"] + "|candidate")
        draws = draws[~rep].copy()
        ctx = ctx[~ctx.qid.isin(dropped)]
    draws["qid"] = qid(draws, "shadow_hi")
    w = draws.qid.value_counts()
    n_draw = int(w.sum())
    n_eff = float(w.sum() ** 2 / (w ** 2).sum())

    cand_ids = set(qid(wanted[wanted.group == "candidate"], "candidate"))
    if dropped:
        cand_ids -= dropped
    ann_ids = set(qid(wanted[wanted.group == "annotated_cds"], "annotated_cds"))
    n_c, n_a = len(cand_ids), len(ann_ids)

    dirs = [d for d in args.search_dirs.split(",") if d]
    best = m22.load_hits(dirs, m22.DBS)

    rows = []

    def put(metric, v, lo=np.nan, hi=np.nan):
        rows.append(dict(domain=args.domain, label=args.label, metric=metric,
                         value=v, lo=lo, hi=hi))

    put("n_candidate", n_c)
    put("n_comparator_draws", n_draw)
    put("n_comparator_distinct", len(w))
    put("n_comparator_effective", n_eff)
    put("n_annotated", n_a)

    full_sets = {}
    for db in m22.DBS:
        st = best.get((db, "struct"))
        full = set(st.loc[(st.qcov >= m22.FULL) & (st.tcov >= m22.FULL),
                          "query"])
        full_sets[db] = full
    full_sets["any_sp_pdb_cath"] = set().union(
        *(full_sets[d] for d in m22.EST_DBS))
    shares = {}
    for db, full in full_sets.items():
        c = len(cand_ids & full) / n_c
        a = len(ann_ids & full) / n_a
        s = float(w[w.index.isin(full)].sum()) / n_draw
        f, lo, hi = share(c, n_c, s, n_eff, a, n_a)
        shares[db] = f
        put(f"fulllen_{db}_candidate", c)
        put(f"fulllen_{db}_comparator", s)
        put(f"fulllen_{db}_annotated", a)
        put(f"share_{db}", f, lo, hi)
    est = [shares[d] for d in m22.EST_DBS]
    put("structural_estimate", np.mean(est) * n_c,
        min(est) * n_c, max(est) * n_c)

    # Diagnostic: the comparator rate directly standardised to the candidate
    # distribution over 5 x 5 cells of candidate 3Di x length quintiles, as
    # supplementary_analyses.md §8 does. For a matched comparator this should
    # barely move; for the published clean arm it reproduces Setonix's
    # re-weighted shares (annotated left crude).
    cc = ctx[ctx.group == "candidate"]
    q3 = np.quantile(cc.three_di_entropy.astype(float), [0.2, 0.4, 0.6, 0.8])
    ql = np.quantile(cc.aa_length.astype(float), [0.2, 0.4, 0.6, 0.8])

    def cell(d):
        return (np.searchsorted(q3, d.three_di_entropy.astype(float), "right") * 5
                + np.searchsorted(ql, d.aa_length.astype(float), "right"))

    p_c = np.bincount(cell(cc), minlength=25) / len(cc)
    sc = ctx[ctx.qid.isin(w.index)].drop_duplicates("qid").set_index("qid")
    sc_cell = pd.Series(cell(sc), index=sc.index).loc[w.index].to_numpy()
    sw_all = w.to_numpy(float)
    put("comparator_cells_empty", int((np.bincount(sc_cell, minlength=25) == 0).sum()))
    for db, full in full_sets.items():
        hit = w.index.isin(full).astype(float)
        num = np.bincount(sc_cell, weights=sw_all * hit, minlength=25)
        den = np.bincount(sc_cell, weights=sw_all, minlength=25)
        rate = np.divide(num, den, out=np.zeros(25), where=den > 0)
        s_std = float((p_c * rate).sum())
        c = len(cand_ids & full) / n_c
        a = len(ann_ids & full) / n_a
        put(f"fulllen_{db}_comparator_std5x5", s_std)
        put(f"share_std5x5_{db}", (c - s_std) / (a - s_std))

    # Ladder. tiers() is imported from 22_rank_candidates.py unchanged.
    desc = m22.load_descriptions(args.descriptions, m22.DBS)
    wide = m22.widen(ctx[ctx.group.eq("candidate")
                         | ctx.qid.isin(w.index)], best, desc, m22.DBS)
    cand = wide[wide.group == "candidate"]
    shad = wide[wide.group == "shadow_hi"].drop_duplicates("qid")
    missing = set(w.index) - set(shad.qid)
    if missing:
        sys.exit(f"ERROR: {len(missing):,} comparator ORFs have no context "
                 "row; the ladder would silently drop them")
    sw = shad.qid.map(w).to_numpy(float)
    tc, ts = m22.tiers(cand), m22.tiers(shad)
    for k in tc:
        cc = int(tc[k].sum())
        ss = float((ts[k].to_numpy(bool) * sw).sum()) * n_c / n_draw
        tag = k.split()[0]
        put(f"ladder_{tag}_candidates", cc)
        put(f"ladder_{tag}_comparator_scaled", ss)
        put(f"ladder_{tag}_excess", cc - ss)

    if args.prodigal:
        p = pd.read_csv(args.prodigal, sep="\t", dtype=STR,
                        usecols=KEY + ["group", "coincides"])
        p["coincides"] = p.coincides.astype(str).map({"True": True, "False": False})
        if dropped:
            p = p[~qid(p, "candidate").isin(dropped) | (p.group != "candidate")]
        rate = {g: float(p.loc[p.group == g, "coincides"].mean())
                for g in ("candidate", "annotated_cds", "intergenic_lo")}
        n = {g: int((p.group == g).sum())
             for g in ("candidate", "annotated_cds", "intergenic_lo")}
        ps = p[p.group == "shadow_hi"].copy()
        ps["qid"] = qid(ps, "shadow_hi")
        ps = ps.drop_duplicates("qid").set_index("qid")
        if set(w.index) - set(ps.index):
            sys.exit("ERROR: comparator ORFs missing from the Prodigal table")
        s = float((ps.loc[w.index, "coincides"].to_numpy(float)
                   * w.to_numpy(float)).sum()) / n_draw
        for g in rate:
            put(f"prodigal_{g}", rate[g])
        put("prodigal_comparator", s)
        f, lo, hi = share(rate["candidate"], n["candidate"], s, n_eff,
                          rate["annotated_cds"], n["annotated_cds"])
        put("prodigal_share_clean_bg", f, lo, hi)
        f2, lo2, hi2 = share(rate["candidate"], n["candidate"],
                             rate["intergenic_lo"], n["intergenic_lo"],
                             rate["annotated_cds"], n["annotated_cds"])
        put("prodigal_share_intergenic_lo_bg", f2, lo2, hi2)
        # The report's rule: take the larger background.
        larger = f if s >= rate["intergenic_lo"] else f2
        put("prodigal_share_larger_bg", larger)
        for db in ("afdb_swissprot", "cath50"):
            put(f"ratio_prodigal_over_{db}", larger / shares[db])

    out = pd.DataFrame(rows)
    out.to_csv(args.out, sep="\t", index=False)
    with pd.option_context("display.max_rows", 200, "display.width", 160):
        print(out.to_string(index=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
