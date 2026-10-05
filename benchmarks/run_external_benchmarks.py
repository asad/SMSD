#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
"""Run checked-in external corpora with explicit local comparison settings.

The Dalke-style pair files were derived from MoleculeNet, not ChEMBL-13.
Search-only timings exclude parsing, conversion, warmup and witness validation.
Results checkpoint after each engine call; --resume skips completed pairs.
SMSD does not expose cancellation, so its timeout count is reported as unknown.
"""
import argparse
import csv
from dataclasses import asdict
import hashlib
import json
from pathlib import Path
import platform
import statistics
import time

import smsd
import smsd._smsd as native
from rdkit import Chem

from mcs_protocol import graph_from_rdkit, measure_pair, summarize

SCRIPT_DIR = Path(__file__).resolve().parent
DATA = SCRIPT_DIR / "data"
DATASETS = {
    "stress": "stress_pairs.tsv",
    "dalke-random": "dalke_random_pairs.tsv",
    "dalke-nn": "dalke_nn_pairs.tsv",
    "tautobase": "chodera_tautobase_subset.txt",
    "smarts": "ehrlich_rarey_smarts.txt",
}


def pairs(dataset):
    path = DATA / DATASETS[dataset]
    if dataset == "tautobase":
        with path.open() as handle:
            for row in csv.reader(handle, skipinitialspace=True):
                if len(row) >= 3 and row[0] != "name":
                    yield row[1].strip(), row[2].strip(), row[0].strip()
    else:
        for line in path.read_text().splitlines():
            if line.strip() and not line.startswith("#"):
                row = line.split("\t")
                if len(row) >= 2:
                    yield row[0], row[1], " / ".join(row[2:4])


def run_smarts(pattern, warmup, repeats, index):
    target = "CC(C)Cc1ccc(CC(C)C(O)=O)cc1"
    rdkit_target = Chem.MolFromSmiles(target)
    # SMARTS needs implicit-hydrogen counts as well as atom/bond flags. The
    # native SMILES parser supplies those; search timings exclude this parsing.
    graph = smsd.parse_smiles(target)
    # Compiled-query searches for both engines; compile cost is recorded separately.
    start = time.perf_counter_ns()
    query = native.compile_smarts(pattern)
    smsd_compile_us = (time.perf_counter_ns()-start)/1000
    start = time.perf_counter_ns()
    rdkit_query = Chem.MolFromSmarts(pattern)
    rdkit_compile_us = (time.perf_counter_ns()-start)/1000
    if rdkit_query is None:
        raise ValueError("RDKit rejected SMARTS")
    calls = {"smsd": lambda: query.matches(graph),
             "rdkit": lambda: rdkit_target.HasSubstructMatch(rdkit_query)}
    for _ in range(warmup):
        for call in calls.values():
            call()
    rows = []
    for trial in range(repeats):
        order = ("smsd", "rdkit") if (index+trial)%2 == 0 else ("rdkit", "smsd")
        for engine in order:
            start = time.perf_counter_ns()
            hit = bool(calls[engine]())
            rows.append({"engine": engine, "trial": trial, "hit": hit,
                         "elapsed_us": (time.perf_counter_ns()-start)/1000})
    hits = {e: sorted({r["hit"] for r in rows if r["engine"] == e}) for e in calls}
    return {"observations": rows, "smsd_hits": hits["smsd"], "rdkit_hits": hits["rdkit"],
            "hit_agreement": hits["smsd"] == hits["rdkit"],
            "smsd_compile_us": smsd_compile_us, "rdkit_compile_us": rdkit_compile_us,
            "speed_comparable": False,
            "comparison_note": "hit agreement on one target does not validate complete SMARTS dialect semantics",
            **{e+"_median_us": statistics.median(r["elapsed_us"] for r in rows if r["engine"] == e)
               for e in calls}}


def aggregate(records):
    valid = [r for r in records if not r.get("error")]
    result = {"cases": len(records), "input_or_search_errors": len(records)-len(valid),
              "speed_comparable_cases": sum(r.get("speed_comparable", False) for r in valid),
              "smsd_timeout_count": None,
              "smsd_timeout_note": "not exposed by native mapping API; wall-clock budget crossings are separate",
              "total_measured_search_seconds": sum(o["elapsed_us"] for r in valid
                  for o in r.get("observations", []))/1e6}
    for e in ("smsd", "rdkit"):
        durations = [r[e+"_median_us"] for r in valid if e+"_median_us" in r]
        result[e+"_median_pair_us"] = statistics.median(durations) if durations else None
        result[e+"_invalid_witness_cases"] = sum(r.get(e+"_valid") is False for r in valid)
        result[e+"_budget_reached_runs"] = sum(r.get(e+"_budget_reached", 0) for r in valid)
        result[e+"_near_budget_runs"] = sum(r.get(e+"_near_budget", 0) for r in valid)
        result[e+"_search_exception_runs"] = sum(bool(o.get("error")) for r in valid
            for o in r.get("observations", []) if o.get("engine") == e)
    result["rdkit_canceled_runs"] = sum(r.get("rdkit_timeouts", 0) for r in valid)
    result["hit_disagreements"] = sum(r.get("hit_agreement") is False for r in valid)
    result["equal_atom_count_cases"] = sum(r.get("smsd_atoms") == r.get("rdkit_atoms")
        for r in valid if "smsd_atoms" in r)
    comparable = [r for r in valid if r.get("speed_comparable")]
    result["paired_ratio_median"] = (statistics.median(r["rdkit_median_us"]/r["smsd_median_us"]
        for r in comparable if r["smsd_median_us"] > 0) if comparable else None)
    return result


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--datasets", nargs="+", choices=DATASETS, default=list(DATASETS))
    parser.add_argument("--policies", nargs="+", choices=("strict", "fmcs", "any"), default=["strict", "fmcs"])
    parser.add_argument("--timeout-sec", type=int, default=10)
    parser.add_argument("--warmup", type=int, default=1)
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--limit", type=int, default=0, help="0 means complete corpus")
    parser.add_argument("--output-dir", type=Path, default=Path("build/local-benchmarks"))
    parser.add_argument("--source-label", default="unspecified")
    parser.add_argument("--smsd-strategy", choices=("native", "auto", "lightweight", "coverage"), default="native")
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args(argv)
    if args.timeout_sec < 1 or args.warmup < 0 or args.repeats < 1 or args.limit < 0:
        parser.error("timeout/repeats must be positive; warmup/limit nonnegative")
    if args.smsd_strategy != "native" and "strict" in args.policies:
        parser.error("auto/lightweight/coverage cannot express the strict aromatic-atom profile; use fmcs or any")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    metadata = {"smsd_version": smsd.__version__, "smsd_package": smsd.__file__,
        "smsd_extension": native.__file__, "smsd_backend": native.gpu_device_info(),
        "rdkit_version": Chem.rdBase.rdkitVersion, "rdkit_package": Chem.__file__,
        "python": platform.python_version(), "platform": platform.platform(),
        "source_label": args.source_label, "timeout_sec": args.timeout_sec,
        "warmup": args.warmup, "repeats": args.repeats, "objective": "atoms",
        "connected_only": True, "induced": False, "parsing_timed": False,
        "protocol": "alternating engine order; independent post-timing witness validation",
        "policies": args.policies, "smsd_strategy": args.smsd_strategy,
        "dataset_sha256": {key: hashlib.sha256((DATA/DATASETS[key]).read_bytes()).hexdigest()
            for key in args.datasets}}
    meta_path = args.output_dir / "metadata.json"
    if args.resume and meta_path.exists() and json.loads(meta_path.read_text()) != metadata:
        parser.error("resume metadata differs; use a separate output directory")
    meta_path.write_text(json.dumps(metadata, indent=2)+"\n")
    summaries = {}
    wall_start = time.perf_counter()
    for dataset in args.datasets:
        dataset_policies = ["smarts"] if dataset == "smarts" else args.policies
        for policy in dataset_policies:
            path = args.output_dir / f"{dataset}-{policy}.jsonl"
            records = [json.loads(line) for line in path.read_text().splitlines()] if args.resume and path.exists() else []
            completed = {r["index"] for r in records}
            source = ([(line.split("\t")[0], "", line.split("\t")[1] if "\t" in line else "")
                       for line in (DATA/DATASETS[dataset]).read_text().splitlines()
                       if line.strip() and not line.startswith("#")]
                      if dataset == "smarts" else list(pairs(dataset)))
            if args.limit:
                source = source[:args.limit]
            with path.open("a" if args.resume else "w", buffering=1) as handle:
                for index, (smi1, smi2, name) in enumerate(source):
                    if index in completed:
                        continue
                    row = {"dataset": dataset, "policy": policy, "index": index, "name": name,
                           "smiles1": smi1, "smiles2": smi2}
                    try:
                        if dataset == "smarts":
                            row.update(run_smarts(smi1, args.warmup, args.repeats, index))
                        else:
                            observations = []
                            # A second JSONL checkpoints completed individual searches, even if a later call fails.
                            checkpoint = args.output_dir / f"{dataset}-{policy}-observations.jsonl"
                            def save_observation(observation):
                                with checkpoint.open("a") as partial:
                                    partial.write(json.dumps({"index": index, **observation})+"\n")
                            observations = measure_pair(smi1, smi2, policy, args.timeout_sec,
                                args.warmup, args.repeats, index, save_observation, args.smsd_strategy)
                            row.update(summarize(observations, policy)); row["observations"] = observations
                    except Exception as exc:
                        row["error"] = str(exc)
                    handle.write(json.dumps(row)+"\n"); records.append(row)
                    if index % 25 == 0 or index+1 == len(source):
                        print(f"{dataset}/{policy}: {index+1}/{len(source)}; elapsed {time.perf_counter()-wall_start:.1f}s", flush=True)
            summaries[f"{dataset}/{policy}"] = aggregate(records)
            (args.output_dir / "summary.json").write_text(json.dumps(summaries, indent=2)+"\n")
            print(json.dumps({"dataset": dataset, "policy": policy, **summaries[f"{dataset}/{policy}"]}), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
