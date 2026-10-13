#!/usr/bin/env python3
"""Build one compact per-ORF table for the phage mutual-information analysis.

Reads the 145 unified JSON chunk outputs (``results/phold_*.json.gz``, the
modernprost-50M re-encoding of phold.gbk) and the aggregated TSV
(``phold_modernprost_50m_mutual_information.tsv``), which carries the
``in_genbank`` label copied from the legacy PHOLD JSON. Joins them on
(input_id, orf_id) and writes a single Parquet file with lengths,
coordinates, the five per-ORF variables, distinct-state counts and the label.

The JSON chunk outputs were produced from FASTA (no GenBank), so their own
``metadata.in_genbank`` is always False; the TSV label is authoritative.
"""

import argparse
import gzip
import json
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd

ENTROPIES = (
    "dna_entropy",
    "protein_entropy",
    "three_di_entropy",
    "twelve_state_entropy",
    "three_di_twelve_state_mutual_information",
)


def iter_records(data):
    for item in data:
        if isinstance(item, list):
            yield from iter_records(item)
        else:
            yield item


def read_chunk(path: Path) -> pd.DataFrame:
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        data = json.load(handle)
    columns = {name: [] for name in (
        "chunk", "input_id", "orf_id", "input_dna_length", "start", "end",
        "strand", "frame", "nt_length", "protein_length", "three_di_length",
        "n_states_3di", "n_states_12", "frac_D_3di", "n_X_protein",
        "starts_with_M", "model_name", "schema_version", *ENTROPIES)}
    for record in iter_records(data):
        for orf_id, feature in record.get("features", {}).items():
            protein = feature["protein"]["aa_sequence"]
            three_di = feature["three_di"]["encoding"]
            twelve = (feature.get("twelve_state") or {}).get("encoding") or ""
            entropy = feature.get("entropy", {})
            location = feature["location"]
            columns["chunk"].append(path.name)
            columns["input_id"].append(record["input_id"])
            columns["orf_id"].append(orf_id)
            columns["input_dna_length"].append(record["input_dna_length"])
            columns["start"].append(location["start"])
            columns["end"].append(location["end"])
            columns["strand"].append(location["strand"])
            columns["frame"].append(location["frame"])
            columns["nt_length"].append(feature["dna"]["length"])
            columns["protein_length"].append(len(protein))
            columns["three_di_length"].append(len(three_di))
            columns["n_states_3di"].append(len(set(three_di)))
            columns["n_states_12"].append(len(set(twelve)))
            columns["frac_D_3di"].append(
                three_di.count("D") / len(three_di) if three_di else np.nan)
            columns["n_X_protein"].append(protein.count("X"))
            columns["starts_with_M"].append(protein[:1] == "M")
            columns["model_name"].append(feature["three_di"].get("model_name"))
            columns["schema_version"].append(record.get("schema_version"))
            for name in ENTROPIES:
                value = entropy.get(name)
                columns[name].append(np.nan if value is None else float(value))
    frame = pd.DataFrame(columns)
    for name in ("chunk", "strand", "model_name", "schema_version"):
        frame[name] = frame[name].astype("category")
    return frame


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("results_directory", type=Path)
    parser.add_argument("tsv", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--processes", type=int, default=8)
    args = parser.parse_args()

    paths = sorted(args.results_directory.glob("phold_*.json.gz"))
    print(f"{len(paths)} JSON chunks", flush=True)
    with Pool(args.processes) as pool:
        frames = pool.map(read_chunk, paths, chunksize=1)
    table = pd.concat(frames, ignore_index=True)
    for name in ("chunk", "strand", "model_name", "schema_version"):
        table[name] = table[name].astype(str).astype("category")
    print(f"JSON ORFs: {len(table):,}", flush=True)

    labels = pd.read_csv(
        args.tsv, sep="\t", dtype={"input_id": str, "orf_id": str, "in_genbank": str},
        usecols=["input_id", "orf_id", *ENTROPIES, "in_genbank"],
    )
    print(f"TSV rows: {len(labels):,}", flush=True)
    print("TSV label values:", labels.in_genbank.value_counts(dropna=False).to_dict())
    duplicated = labels.duplicated(["input_id", "orf_id"]).sum()
    print(f"Duplicated (input_id, orf_id) keys in TSV: {duplicated}")
    labels = labels.rename(columns={name: f"tsv_{name}" for name in ENTROPIES})

    merged = table.merge(labels, on=["input_id", "orf_id"], how="outer",
                         indicator=True, validate="one_to_one")
    print("Join:", merged._merge.value_counts().to_dict())
    for name in ENTROPIES:
        difference = (merged[name] - merged[f"tsv_{name}"]).abs().max()
        print(f"max |JSON - TSV| {name}: {difference}")
    merged = merged[merged._merge == "both"].drop(
        columns=["_merge", *[f"tsv_{name}" for name in ENTROPIES]])
    merged["in_genbank"] = merged["in_genbank"].str.lower().map(
        {"true": True, "false": False})
    print("Label counts:", merged.in_genbank.value_counts(dropna=False).to_dict())
    print("Model names:", merged.model_name.value_counts().to_dict())
    print("Schema:", merged.schema_version.value_counts().to_dict())
    print("Protein length summary:\n", merged.protein_length.describe().to_string())
    print("Records:", merged.input_id.nunique())
    args.output.parent.mkdir(parents=True, exist_ok=True)
    merged.to_parquet(args.output, index=False)
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
