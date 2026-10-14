#!/usr/bin/env python3
"""Draw a clean-shadow comparator matched to the candidates AFTER frame filtering.

The published clean comparator is the antisense + frameshift subset of a
shadow_hi arm that had been matched to the candidates BEFORE that filter, so it
is no longer matched (supplementary_analyses.md §8). This draws a new one from
the full clean pool written by 39_shadow_pool.py, one partner per candidate,
without replacement, under three constraints:

  * same contig-end truncation status (truncated ORFs have no stop codon and
    hit full-length less often; the old clean arm was 2.4% truncated against
    7.6% of candidates);
  * 3Di entropy within --caliper-3di bits;
  * length within --caliper-len on the log scale (0.10 = about +/-10%).

Within those calipers the nearest partner is taken, preferring -- as the
original shadow_hi arm did -- a shadow from the SAME GENOME, then the same
chunk (250 genomes), then anywhere in the domain. Candidates with no partner
inside the calipers are retried at 2x and 4x the calipers; any still unmatched
get their nearest partner WITH replacement. Every fallback is counted and
reported, never silent.

Candidates are processed scarcest-first (descending 3Di), so the high-3Di
tail, which is where clean shadows are thin, gets first claim on the pool.

Outputs (in --out-dir):
  comparator.tsv.gz     the new arm, one row per candidate, with the partner
                        candidate, match level and achieved differences
  wanted_v7.tsv         the published wanted list with its shadow_hi arm
                        REPLACED by the new comparator (still labelled
                        shadow_hi, so 18/20/22/24 run unchanged)
  to_search.tsv         comparator ORFs not already searched as shadow_hi
  balance.tsv / balance.txt   candidate vs old clean vs new comparator

  40_match_clean_comparator.py --domain bac \
      --candidates-dir missed_genes/bac --pool-dir matched_clean/bac/pool \
      --old-wanted missed_genes/full_bac/wanted.tsv \
      --out-dir matched_clean/bac
"""
import argparse
import glob
import sys
from pathlib import Path

import numpy as np
import pandas as pd

CLEAN = {"opposite strand", "same strand, frameshift"}
KEY = ["genome", "input_id", "orf_id"]
WANTED_COLS = ["domain", "chunk", "genome", "input_id", "orf_id", "group",
               "aa_length", "three_di_entropy", "protein_entropy", "truncated"]
STR = {"chunk": str, "genome": str, "input_id": str, "orf_id": str,
       "domain": str}


def read_many(pattern):
    files = sorted(glob.glob(pattern))
    if not files:
        sys.exit(f"ERROR: nothing matches {pattern}")
    return pd.concat([pd.read_csv(f, sep="\t", dtype=STR) for f in files],
                     ignore_index=True)


def as_bool(s):
    return s.astype(str).map({"True": True, "False": False}).astype(bool)


def qid(d):
    return d.genome + "|" + d.input_id + "|" + d.orf_id


class Matcher:
    def __init__(self, pool, c3, cl):
        self.pool = pool
        self.c3, self.cl = c3, cl
        self.ent = pool.three_di_entropy.to_numpy(float)
        self.loglen = np.log(pool.aa_length.to_numpy(float))
        self.trunc = pool.truncated.to_numpy(bool)
        self.genome = pool.genome.to_numpy(object)
        self.used = np.zeros(len(pool), dtype=bool)

    def best(self, idx, ent, loglen, scale, genome=None):
        """Nearest unused row of idx inside scale x calipers, or -1."""
        if not len(idx):
            return -1
        idx = idx[~self.used[idx]]
        d3 = np.abs(self.ent[idx] - ent)
        dl = np.abs(self.loglen[idx] - loglen)
        ok = (d3 <= self.c3 * scale) & (dl <= self.cl * scale)
        if not ok.any():
            return -1
        idx, d = idx[ok], (d3[ok] / self.c3 + dl[ok] / self.cl)
        if genome is not None:
            same = self.genome[idx] == genome
            if same.any():
                idx, d = idx[same], d[same]
        return int(idx[int(np.argmin(d))])


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--domain", required=True)
    ap.add_argument("--candidates-dir", required=True)
    ap.add_argument("--pool-dir", required=True)
    ap.add_argument("--old-wanted", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--caliper-3di", type=float, default=0.025)
    ap.add_argument("--caliper-len", type=float, default=0.10)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    dom = args.domain

    cand = read_many(f"{args.candidates_dir}/{dom}_*.candidates.tsv.gz")
    cand["truncated"] = as_bool(cand.truncated)
    pool = read_many(f"{args.pool_dir}/{dom}_*.pool.tsv.gz")
    pool["truncated"] = as_bool(pool.truncated)
    old = pd.read_csv(args.old_wanted, sep="\t", dtype=STR)
    n_old_cand = int((old.group == "candidate").sum())
    if len(cand) != n_old_cand:
        sys.exit(f"ERROR: {len(cand):,} candidates in the tables but "
                 f"{n_old_cand:,} in {args.old_wanted}")
    for d in (cand, pool):
        d["aa_length"] = d.aa_length.astype(int)
        d["three_di_entropy"] = d.three_di_entropy.astype(float)

    lines = []

    def emit(s=""):
        print(s)
        lines.append(s)

    fc = pool.cds_frame_class.fillna("")
    emit(f"domain {dom}: {len(cand):,} candidates; full shadow_hi pool "
         f"{len(pool):,}")
    for k, v in fc.value_counts().items():
        emit(f"  {k:<32}{v:>12,}  {100 * v / len(pool):5.1f}%")
    clean = pool[fc.isin(CLEAN)].reset_index(drop=True)
    emit(f"clean pool (antisense + frameshift): {len(clean):,}")

    # Old clean arm, for the balance table: published shadow_hi rows whose
    # frame class (from the pool, same code as the context tables) is clean.
    old_sh = old[old.group == "shadow_hi"][KEY]
    old_clean = old_sh.merge(clean, on=KEY, how="inner")

    m = Matcher(clean, args.caliper_3di, args.caliper_len)
    by_chunk = {c: np.asarray(i, dtype=int)
                for c, i in clean.groupby("chunk").indices.items()}
    rng = np.random.default_rng(args.seed)
    order = rng.permutation(len(cand))
    order = order[np.argsort(-cand.three_di_entropy.to_numpy()[order],
                             kind="stable")]
    c_ent = cand.three_di_entropy.to_numpy(float)
    c_ll = np.log(cand.aa_length.to_numpy(float))
    c_tr = cand.truncated.to_numpy(bool)
    c_gen = cand.genome.to_numpy(object)
    c_chunk = cand.chunk.to_numpy(object)

    partner = np.full(len(cand), -1, dtype=np.int64)
    level = np.full(len(cand), "", dtype=object)

    # Pass 1: within chunk, same genome preferred, at 1x calipers.
    for i in order:
        idx = by_chunk.get(c_chunk[i], np.zeros(0, dtype=int))
        idx = idx[m.trunc[idx] == c_tr[i]]
        j = m.best(idx, c_ent[i], c_ll[i], 1.0, genome=c_gen[i])
        if j >= 0:
            m.used[j] = True
            partner[i] = j
            level[i] = "same genome" if m.genome[j] == c_gen[i] else "same chunk"
    n1 = int((partner >= 0).sum())
    emit(f"\npass 1 (within chunk): {n1:,} of {len(cand):,}")

    # Pass 2: domain-wide, 3Di-sorted windows, at 1x, 2x, 4x calipers.
    srt = {t: np.flatnonzero(m.trunc == t) for t in (False, True)}
    for t in srt:
        srt[t] = srt[t][np.argsort(m.ent[srt[t]], kind="stable")]
    ent_sorted = {t: m.ent[srt[t]] for t in srt}
    for scale, name in ((1.0, "domain"), (2.0, "domain 2x caliper"),
                        (4.0, "domain 4x caliper")):
        left = [i for i in order if partner[i] < 0]
        if not left:
            break
        for i in left:
            t = bool(c_tr[i])
            lo = np.searchsorted(ent_sorted[t], c_ent[i] - m.c3 * scale, "left")
            hi = np.searchsorted(ent_sorted[t], c_ent[i] + m.c3 * scale, "right")
            j = m.best(srt[t][lo:hi], c_ent[i], c_ll[i], scale)
            if j >= 0:
                m.used[j] = True
                partner[i] = j
                level[i] = name
        emit(f"pass 2 ({name}): {int((partner >= 0).sum()):,} matched")

    # Last resort: nearest clean shadow of the same truncation status, WITH
    # replacement. Counted so it cannot hide.
    left = np.flatnonzero(partner < 0)
    for i in left:
        t = bool(c_tr[i])
        idx = srt[t]
        d = (np.abs(m.ent[idx] - c_ent[i]) / m.c3
             + np.abs(m.loglen[idx] - c_ll[i]) / m.cl)
        partner[i] = int(idx[int(np.argmin(d))])
        level[i] = "with replacement"
    emit("\nmatch level:")
    for k, v in pd.Series(level).value_counts().items():
        emit(f"  {k:<24}{v:>10,}  {100 * v / len(cand):6.2f}%")

    comp = clean.iloc[partner].reset_index(drop=True).copy()
    comp["group"] = "shadow_hi"
    comp["match_level"] = level
    comp["cand_genome"] = cand.genome.to_numpy()
    comp["cand_input_id"] = cand.input_id.to_numpy()
    comp["cand_orf_id"] = cand.orf_id.to_numpy()
    comp["d_3di"] = comp.three_di_entropy.to_numpy() - c_ent
    comp["d_len"] = comp.aa_length.to_numpy() - cand.aa_length.to_numpy()
    n_dup = int(comp.duplicated(KEY).sum())
    emit(f"duplicated comparator ORFs (with-replacement draws): {n_dup:,}")
    comp.to_csv(out / "comparator.tsv.gz", sep="\t", index=False,
                compression="gzip")

    # ---- balance ---------------------------------------------------------
    q3 = np.quantile(c_ent, [0.2, 0.4, 0.6, 0.8])
    ql = np.quantile(cand.aa_length, [0.2, 0.4, 0.6, 0.8])
    rows = []
    arms = {"candidate": cand, "old clean shadow": old_clean,
            "new matched clean": comp}
    emit("\nBALANCE (top quintile = candidates' top 3Di quintile, "
         f">= {q3[-1]:.3f} bits)")
    emit(f"  {'arm':<20}{'n':>10}{'med aa':>8}{'med 3Di':>9}{'mean 3Di':>9}"
         f"{'top-Q 3Di':>10}{'trunc':>7}{'5x5 TVD':>9}")
    cell_c = None
    for name, d in arms.items():
        e = d.three_di_entropy.to_numpy(float)
        ln = d.aa_length.to_numpy(float)
        cell = (np.searchsorted(q3, e, "right") * 5
                + np.searchsorted(ql, ln, "right"))
        dist = np.bincount(cell, minlength=25) / len(cell)
        if cell_c is None:
            cell_c = dist
        tvd = 0.5 * np.abs(dist - cell_c).sum()
        r = dict(arm=name, n=len(d), median_aa=float(np.median(ln)),
                 median_3di=float(np.median(e)), mean_3di=float(e.mean()),
                 top_quintile_3di=float((e >= q3[-1]).mean()),
                 truncated=float(np.mean(d.truncated.astype(bool))),
                 tvd_5x5=float(tvd))
        rows.append(r)
        emit(f"  {name:<20}{r['n']:>10,}{r['median_aa']:>8.0f}"
             f"{r['median_3di']:>9.3f}{r['mean_3di']:>9.3f}"
             f"{100 * r['top_quintile_3di']:>9.1f}%{100 * r['truncated']:>6.1f}%"
             f"{r['tvd_5x5']:>9.4f}")
    pd.DataFrame(rows).to_csv(out / "balance.tsv", sep="\t", index=False)
    emit(f"\n  achieved |d 3Di|: median {np.median(np.abs(comp.d_3di)):.4f}, "
         f"p99 {np.quantile(np.abs(comp.d_3di), 0.99):.4f}")
    emit(f"  achieved |d len| / len: median "
         f"{np.median(np.abs(comp.d_len) / cand.aa_length):.4f}, p99 "
         f"{np.quantile(np.abs(comp.d_len) / cand.aa_length, 0.99):.4f}")
    emit(f"  frame classes of new arm: "
         + "; ".join(f"{k}={v:,}" for k, v in
                     comp.cds_frame_class.value_counts().items()))

    # ---- wanted lists ----------------------------------------------------
    new = comp.copy()
    new["domain"] = dom
    new["truncated"] = new.truncated.map({True: "True", False: "False"})
    new = new.drop_duplicates(KEY)
    keep = old[old.group != "shadow_hi"].copy()
    wanted = pd.concat([keep[WANTED_COLS], new[WANTED_COLS]],
                       ignore_index=True)
    wanted = wanted.sort_values("chunk", kind="stable")
    wanted.to_csv(out / "wanted_v7.tsv", sep="\t", index=False)
    already = set(qid(old_sh))
    todo = new[~qid(new).isin(already)]
    todo[WANTED_COLS].sort_values("chunk", kind="stable").to_csv(
        out / "to_search.tsv", sep="\t", index=False)
    emit(f"\nwanted_v7.tsv: {len(wanted):,} rows "
         f"({len(new):,} distinct comparator ORFs)")
    emit(f"already searched as shadow_hi: {len(new) - len(todo):,}; "
         f"to search: {len(todo):,}")
    (out / "balance.txt").write_text("\n".join(lines) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
