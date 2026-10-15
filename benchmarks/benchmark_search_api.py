#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
"""Measure prepared search APIs using the shared search_api_cases.tsv corpus.

Each JSONL row retains a measured call and its complete returned mappings.
Parsing and query preparation are included only in their declared phases.
Validate the output with search_api_validation.py before comparing timings.
"""
import argparse
import csv
import json
from pathlib import Path
import time

import smsd


def chemistry(profile):
    result = smsd.ChemOptions()
    result.match_atom_type = True
    result.match_bond_order = smsd.BondOrderMode.STRICT
    result.aromaticity_mode = smsd.AromaticityMode.FLEXIBLE
    result.match_formal_charge = False
    result.match_isotope = False
    result.use_chirality = False
    result.use_bond_stereo = False
    result.ring_matches_ring_only = False
    result.complete_rings_only = False
    result.ring_fusion_mode = smsd.RingFusionMode.IGNORE
    result.tautomer_aware = False
    if profile in ("strict", "rings"):
        result.aromaticity_mode = smsd.AromaticityMode.STRICT
        result.match_formal_charge = True
        result.match_isotope = True
        result.ring_matches_ring_only = True
        result.complete_rings_only = profile == "rings"
    elif profile == "any":
        result.match_bond_order = smsd.BondOrderMode.ANY
    elif profile == "tautomer":
        result.tautomer_aware = True
        result.match_bond_order = smsd.BondOrderMode.LOOSE
    elif profile == "stereo":
        result.use_chirality = result.use_bond_stereo = True
    elif profile != "default":
        raise ValueError("Unknown chemistry profile: " + profile)
    return result


def options_for(case):
    if case["operation"] == "mcs":
        options = smsd.MCS.Options()
        options.objective = getattr(smsd.MCS.Objective, case["objective"])
        options.connectivity = (smsd.MCS.Connectivity.CONNECTED
                                if case["connectivity"] == "Connected"
                                else smsd.MCS.Connectivity.FRAGMENTS)
        options.min_fragment_atoms = 1
        options.max_fragments = 2147483647
        options.atom_weights = [float(value) for value in case["weights"].split(",")] if case["weights"] else []
    else:
        options = smsd.Substructure.Options()
        options.induced = bool(case["induced"])
    options.chemistry = chemistry(case["profile"])
    options.timeout_ms = case["timeout_ms"]
    options.node_limit = 0
    return options


def prepare(case, graph, options):
    if case["operation"] == "smarts":
        return smsd.Substructure.smarts(case["query"]).with_options(options)
    factory = smsd.MCS if case["operation"] == "mcs" else smsd.Substructure
    return factory.find(graph).with_options(options)


def search(case, prepared, target):
    if case["limit"]:
        return prepared.match_all(target, case["limit"])
    return prepared.match(target)


def graph_data(graph):
    if graph is None:
        return None
    elements, charges, masses = graph.atomic_num, graph.formal_charge, graph.mass_number
    hydrogens, aromatic, rings = graph.hydrogen_count, graph.aromatic, graph.ring
    return {
        "atoms": [dict(element=elements[atom], charge=charges[atom], isotope=masses[atom],
                       hydrogens=hydrogens[atom], aromatic=bool(aromatic[atom]), ring=bool(rings[atom]))
                  for atom in range(graph.n)],
        "bonds": [[first, second, order, bool(graph.bond_aromatic(first, second)),
                   bool(graph.bond_in_ring(first, second))]
                  for first, second, order in graph.bonds()],
    }


def result_data(result):
    statuses = {smsd.SearchStatus.COMPLETED: "Completed", smsd.SearchStatus.TIMED_OUT: "TimedOut",
                smsd.SearchStatus.CANCELLED: "Cancelled", smsd.SearchStatus.LIMIT_REACHED: "LimitReached"}
    return dict(status=statuses[result.status], optimal=result.optimal, complete=result.complete,
                score=result.score, nodes_visited=result.nodes_visited,
                matches=[dict(atoms=sorted(match.atom_mapping.items()), bonds=match.bond_mapping,
                              fragments=match.fragments) for match in result.matches])


def measure_case(case, arguments, output):
    query_graph = None if case["operation"] == "smarts" else smsd.parse_smiles(case["query"])
    target_graph = smsd.parse_smiles(case["target"])
    options = options_for(case)
    targets = [target_graph] * case["batch_size"]
    for phase in case["phases"].split(","):
        if arguments.phase and arguments.phase != phase:
            continue
        prepared = prepare(case, query_graph, options) if phase in ("prepared", "batch") else None
        if phase == "prepare":
            call = lambda: prepare(case, query_graph, options)
        elif phase == "prepared":
            call = lambda: search(case, prepared, target_graph)
        elif phase == "end_to_end":
            def call():
                graph = None if case["operation"] == "smarts" else smsd.parse_smiles(case["query"])
                target = smsd.parse_smiles(case["target"])
                return search(case, prepare(case, graph, options), target)
        elif phase == "batch":
            if case["limit"]:
                raise ValueError("Batch cases require first-match searches")
            call = lambda: prepared.match_many(targets, num_threads=1)
        else:
            raise ValueError("Unknown phase: " + phase)
        for _ in range(arguments.warmup):
            call()
        for sample in range(arguments.repeats):
            result = None
            started = time.perf_counter_ns()
            result = call()
            elapsed = time.perf_counter_ns() - started
            results = [] if phase == "prepare" else result if phase == "batch" else [result]
            row = dict(case_id=case["case_id"], language="python", phase=phase, sample=sample,
                       elapsed_ns=elapsed, results=[result_data(value) for value in results])
            if sample == 0:
                row.update(query_graph=graph_data(query_graph), target_graph=graph_data(target_graph))
            output.write(json.dumps(row, separators=(",", ":")) + "\n")
        output.flush()


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cases", type=Path, default=Path("benchmarks/search_api_cases.tsv"))
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--warmup", type=int, default=2)
    parser.add_argument("--repeats", type=int, default=5)
    parser.add_argument("--phase", choices=("prepare", "prepared", "end_to_end", "batch"))
    parser.add_argument("--case")
    args = parser.parse_args(argv)
    if args.warmup < 0 or args.repeats < 1:
        parser.error("Warmups must be non-negative and repeats positive")
    with args.cases.open(newline="") as source:
        cases = list(csv.DictReader(source, delimiter="\t"))
    if args.case:
        cases = [case for case in cases if case["case_id"] == args.case]
        if not cases:
            parser.error("Unknown case")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w") as output:
        for case in cases:
            for name in ("induced", "timeout_ms", "limit", "batch_size"):
                case[name] = int(case[name])
            try:
                measure_case(case, args, output)
            except Exception as failure:
                raise RuntimeError("Benchmark failed: " + case["case_id"]) from failure


if __name__ == "__main__":
    main()
