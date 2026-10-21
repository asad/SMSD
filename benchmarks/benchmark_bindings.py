#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
# Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
# See the NOTICE file for attribution, trademark, and algorithm IP terms.
"""Measure Python binding setup, cache, copies and batched calls locally.

These are microbenchmarks, not whole-application speed claims. The serial and
single-thread batch MCS measurements use the same already-parsed graphs and
search options; their checksums must agree. Prepared and fresh query setup
are measured separately.
"""
import argparse
import json
from pathlib import Path
import platform
import statistics
import time

import smsd
import smsd._smsd as native
from rdkit import Chem


def measure(call, iterations, warmup, repeats):
    for _ in range(warmup):
        call()
    observations, checksum = [], None
    for _ in range(repeats):
        start = time.perf_counter_ns()
        total = 0
        for _ in range(iterations):
            total += call()
        observations.append((time.perf_counter_ns()-start)/1000/iterations)
        if checksum is not None and checksum != total:
            raise RuntimeError("unstable operation checksum")
        checksum = total
    return {"median_us_per_call": statistics.median(observations),
            "all_us_per_call": observations, "checksum": checksum,
            "iterations": iterations, "repeats": repeats}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--iterations", type=int, default=1000)
    parser.add_argument("--repeats", type=int, default=5)
    parser.add_argument("--warmup", type=int, default=5)
    parser.add_argument("--output", type=Path, default=Path("build/local-benchmarks/bindings.json"))
    parser.add_argument("--source-label", default="unspecified")
    args = parser.parse_args(argv)
    if args.iterations < 1 or args.repeats < 1 or args.warmup < 0:
        parser.error("iterations/repeats must be positive; warmup nonnegative")
    query_smiles, target_smiles = "c1ccccc1", "Cc1ccccc1"
    query, target = smsd.parse_smiles(query_smiles), smsd.parse_smiles(target_smiles)
    query.prewarm(); target.prewarm()
    rdkit_target = Chem.MolFromSmiles(target_smiles)
    targets = [target] * 32
    large_target = smsd.parse_smiles(".".join([target_smiles]*64))
    large_target.prewarm()
    large_targets = [large_target]*32
    compiled_query = smsd.Substructure.smarts("c1ccccc1")
    opts = smsd.MCS.Options()
    opts.timeout_ms = 1000
    prepared = smsd.MCS.find(query).with_options(opts)
    native_query = native.MCS.find(query).with_options(opts)
    substructure = smsd.Substructure.find(query)
    operations = {
        "parse_smiles": lambda: smsd.parse_smiles(target_smiles).n,
        "from_rdkit_uncached": lambda: smsd.from_rdkit(rdkit_target, use_cache=False).n,
        "from_rdkit_cached": lambda: smsd.from_rdkit(rdkit_target, use_cache=True).n,
        "prewarmed_scalar_n": lambda: target.n,
        "atomic_numbers_copy": lambda: len(target.atomic_num),
        "neighbors_copy": lambda: len(target.neighbors),
        "bonds_copy": lambda: len(target.bonds()),
        "native_mcs_graphs": lambda: len(native_query.match(target).matches[0].atom_mapping),
        "public_prepared_graphs": lambda: prepared.match(target).best.atom_count,
        "public_fresh_graphs": lambda: smsd.MCS.find(query).with_timeout(1000).match(target).best.atom_count,
        "serial_32_mcs_sizes": lambda: sum(prepared.match(t).best.atom_count for t in targets),
        "batch_32_mcs_sizes_one_thread": lambda: sum(result.best.atom_count for result in prepared.match_many(targets, num_threads=1)),
        "batch_32_mcs_mappings_one_thread": lambda: sum(len(result.best.atom_mapping) for result in prepared.match_many(targets, num_threads=1)),
        "serial_32_large_mcs_sizes": lambda: sum(prepared.match(t).best.atom_count for t in large_targets),
        "batch_32_large_mcs_sizes_one_thread": lambda: sum(result.best.atom_count for result in prepared.match_many(large_targets, num_threads=1)),
        "serial_32_substructure": lambda: sum(substructure.match(t).found for t in targets),
        "batch_32_substructure_one_thread": lambda: sum(result.found for result in substructure.match_many(targets, num_threads=1)),
        "serial_32_large_substructure": lambda: sum(substructure.match(t).found for t in large_targets),
        "batch_32_large_substructure_one_thread": lambda: sum(result.found for result in substructure.match_many(large_targets, num_threads=1)),
        "serial_32_compiled_smarts": lambda: sum(compiled_query.match(t).found for t in targets),
        "batch_32_compiled_smarts": lambda: sum(result.found for result in compiled_query.match_many(targets)),
        "serial_32_large_compiled_smarts": lambda: sum(compiled_query.match(t).found for t in large_targets),
        "batch_32_large_compiled_smarts": lambda: sum(result.found for result in compiled_query.match_many(large_targets)),
    }
    results = {}
    for name, call in operations.items():
        count = max(1, args.iterations//32) if "32_" in name else args.iterations
        results[name] = measure(call, count, args.warmup, args.repeats)
        print(f"{name}: {results[name]['median_us_per_call']:.3f}us", flush=True)
    checksums = [results[key]["checksum"] for key in
                 ("serial_32_mcs_sizes", "batch_32_mcs_sizes_one_thread", "batch_32_mcs_mappings_one_thread")]
    if len(set(checksums)) != 1:
        raise RuntimeError(f"serial/batch result mismatch: {checksums}")
    for serial, batch in (("serial_32_large_mcs_sizes", "batch_32_large_mcs_sizes_one_thread"),
                          ("serial_32_substructure", "batch_32_substructure_one_thread"),
                          ("serial_32_large_substructure", "batch_32_large_substructure_one_thread"),
                          ("serial_32_compiled_smarts", "batch_32_compiled_smarts"),
                          ("serial_32_large_compiled_smarts", "batch_32_large_compiled_smarts")):
        if results[serial]["checksum"] != results[batch]["checksum"]:
            raise RuntimeError(f"serial/batch result mismatch: {serial}, {batch}")
    record = {"source_label": args.source_label, "smsd_version": smsd.__version__,
              "smsd_package": smsd.__file__, "smsd_extension": native.__file__,
              "smsd_backend": native.gpu_device_info(), "rdkit_version": Chem.rdBase.rdkitVersion,
              "python": platform.python_version(), "results": results}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(record, indent=2)+"\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
