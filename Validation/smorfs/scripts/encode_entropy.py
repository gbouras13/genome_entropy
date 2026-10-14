#!/usr/bin/env python
"""Translate smORF chunks, encode 3Di/12-state with genome_entropy, write per-sequence entropies.

Usage: encode_entropy.py --outdir DIR [--encoding-size N] [--limit N] chunk.tsv [chunk.tsv ...]
Each chunk is a headerless TSV of seqid<TAB>orf_nt. Chunks with existing output are skipped.
"""
import argparse
import gzip
import os
import sys
import time

import PyGeneticCode
from Bio.Seq import Seq

from genome_entropy.encode3di import ModernProstThreeDiEncoder
from genome_entropy.entropy.shannon import shannon_entropy

TABLE = 11
SLICE = 200_000


def translate(nt):
    nt = nt.upper()
    aa = PyGeneticCode.translate(nt, TABLE) if set(nt) <= set("ACGT") else str(Seq(nt).translate(table=TABLE))
    return aa[:-1] if aa.endswith("*") else aa


def process(path, enc, outdir, encoding_size, limit):
    base = os.path.basename(path).removesuffix(".tsv")
    out = os.path.join(outdir, f"{base}.entropy.tsv")
    if os.path.exists(out):
        print(f"skip {base}: output exists", flush=True)
        return
    dataset = "positive" if base.startswith("pos") else "negative"
    ids, nts = [], []
    with open(path) as fh:
        for i, line in enumerate(fh):
            if limit and i >= limit:
                break
            sid, nt = line.rstrip("\n").split("\t")
            ids.append(sid)
            nts.append(nt)
    aas = [translate(nt) for nt in nts]
    n_internal_stop = sum("*" in a for a in aas)
    n_empty = sum(len(a) == 0 for a in aas)
    print(f"{base}: {len(ids)} seqs, {n_internal_stop} with internal stops, {n_empty} empty", flush=True)

    tmp, tmp_seq = out + ".tmp", os.path.join(outdir, f"{base}.states.tsv.gz.tmp")
    t0 = time.time()
    with open(tmp, "w") as fo, gzip.open(tmp_seq, "wt", compresslevel=3) as fs:
        fo.write("seqid\tdataset\tnt_length\taa_length\torf_nt_entropy\taa_entropy\tthree_di_entropy\ttwelve_state_entropy\n")
        fs.write("seqid\taa\tthree_di\ttwelve_state\n")
        for s in range(0, len(ids), SLICE):
            # Empty translations cannot be encoded; give them empty structural states.
            todo = [j for j in range(s, min(s + SLICE, len(ids))) if aas[j]]
            encs = dict(zip(todo, enc.encode([aas[j] for j in todo], encoding_size=encoding_size)))
            for j in range(s, min(s + SLICE, len(ids))):
                e = encs.get(j)
                tdi = e.three_di if e else ""
                tws = (e.twelve_state or "") if e else ""
                fo.write(f"{ids[j]}\t{dataset}\t{len(nts[j])}\t{len(aas[j])}\t"
                         f"{shannon_entropy(nts[j].upper()):.6f}\t{shannon_entropy(aas[j]):.6f}\t"
                         f"{shannon_entropy(tdi):.6f}\t{shannon_entropy(tws):.6f}\n")
                fs.write(f"{ids[j]}\t{aas[j]}\t{tdi}\t{tws}\n")
            done = min(s + SLICE, len(ids))
            el = time.time() - t0
            print(f"{base}: {done}/{len(ids)} in {el:.0f}s ({done / el:.0f} seq/s)", flush=True)
    os.rename(tmp_seq, tmp_seq.removesuffix(".tmp"))
    os.rename(tmp, out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--encoding-size", type=int, default=10000)
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("chunks", nargs="+")
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    enc = ModernProstThreeDiEncoder("gbouras13/modernprost-50M", device="cuda")
    for c in a.chunks:
        process(c, enc, a.outdir, a.encoding_size, a.limit)


if __name__ == "__main__":
    sys.exit(main())
