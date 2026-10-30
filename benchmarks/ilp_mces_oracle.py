#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
# Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
# See the NOTICE file for attribution, trademark, and algorithm IP terms.
"""Bounded independent ANY-bond MCES relaxation, with optional SciPy.

Elements must match. Charges, isotope, aromaticity, rings, stereo and tautomers
are not equality predicates. Bonds are selected explicitly: mapped nonedges and
unselected edges impose no restriction. Atom mappings are injective and have no
isolated selected atoms. The objective is bonds first, then atoms.

Minimum fragment size (at least two) and a maximum fragment count describe the
original contract but are DROPPED from the model. Its optimum/dual bound is an
upper bound for that original contract. A relaxed incumbent is a lower bound
for the original contract only after its selected-edge components pass those
fragment requirements. Solver optimality alone does not prove the original
contract optimum. Numerical solver bounds are not formally verified certificates.

Graph JSON is shared with ilp_mcis_oracle. No SMSD search or production dependency
is used. Time/model/node caps preserve stopped outcomes; the time limit is
cooperative, not a hard process, memory or thread limit.
"""

import argparse
import importlib.util
import json
import math
from pathlib import Path
import time
import warnings

_SPEC = importlib.util.spec_from_file_location("_mcs_oracle_graph", Path(__file__).with_name("ilp_mcis_oracle.py"))
_SHARED = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(_SHARED)
_Graph, _Stopped = _SHARED._Graph, _SHARED._Stopped


def _fragments(mapping, edges):
    remaining = set(mapping)
    adjacent = {atom: [] for atom in mapping}
    for first, second, _, _ in edges:
        adjacent[first].append(second)
        adjacent[second].append(first)
    fragments = []
    while remaining:
        reached = {min(remaining)}
        pending = list(reached)
        while pending:
            for neighbor in adjacent[pending.pop()]:
                if neighbor not in reached:
                    reached.add(neighbor)
                    pending.append(neighbor)
        remaining -= reached
        fragments.append(sorted(reached))
    return fragments


def _check_witness(query, target, pairs, edges):
    mapping = {}
    for query_atom, target_atom in pairs:
        if (type(query_atom) is not int or type(target_atom) is not int
                or not 0 <= query_atom < query.n or not 0 <= target_atom < target.n
                or query_atom in mapping or target_atom in mapping.values()
                or query.elements[query_atom] != target.elements[target_atom]):
            raise ValueError("MCES atom mapping is not an element-compatible injection")
        mapping[query_atom] = target_atom
    source_edges, image_edges, used = set(), set(), set()
    for first, second, image_first, image_second in edges:
        source = tuple(sorted((first, second)))
        image = tuple(sorted((image_first, image_second)))
        if (any(type(value) is not int for value in (first, second, image_first, image_second))
                or mapping.get(first) != image_first or mapping.get(second) != image_second
                or source not in query.bonds or image not in target.bonds
                or source in source_edges or image in image_edges):
            raise ValueError("MCES selected edge is missing, duplicated or inconsistent with the mapping")
        source_edges.add(source)
        image_edges.add(image)
        used.update((first, second))
    if used != set(mapping):
        raise ValueError("MCES witness contains an isolated selected atom")
    return mapping, _fragments(mapping, edges)


def check_any_mces_witness(query_document, target_document, atom_mapping, bond_mapping,
                           *, min_fragment_atoms=3, max_fragments=None, max_atoms=256):
    """Validate a retained MCES witness; never infer an optimum from its presence."""
    if (type(min_fragment_atoms) is not int or min_fragment_atoms < 2
            or max_fragments is not None and (type(max_fragments) is not int or max_fragments <= 0)):
        raise ValueError("Non-isolated MCES requires min_fragment_atoms >= 2 and a positive fragment cap")
    query, target = _Graph(query_document, max_atoms), _Graph(target_document, max_atoms)
    pairs = list(atom_mapping.items()) if isinstance(atom_mapping, dict) else list(atom_mapping)
    edges = [tuple(edge) for edge in bond_mapping]
    mapping, fragments = _check_witness(query, target, pairs, edges)
    return {"mapping": mapping, "bond_mapping": edges, "atom_count": len(mapping), "bond_count": len(edges),
            "fragments": fragments, "original_contract_feasible":
                all(len(fragment) >= min_fragment_atoms for fragment in fragments)
                and (max_fragments is None or len(fragments) <= max_fragments)}


def solve_mces_any_relaxation(query_document, target_document, *, min_fragment_atoms=3,
                              max_fragments=None, time_limit=15.0, max_atoms=256,
                              max_variables=50000, max_rows=200000, max_nonzeros=1000000,
                              max_edge_checks=2000000, node_limit=100000, solver_threads=1):
    """Return raw solver limits/bounds and independently checked selected edges.

    Binary x(q,t) maps atoms. Binary z(q-edge,t-edge,orientation) selects a common
    edge and requires both endpoint x variables. Injectivity and edge uniqueness
    prevent reuse; each selected atom requires an incident selected edge. Thus
    integer solutions are precisely the unrestricted-component non-isolated MCES
    witnesses under ANY bond order. With W=min(nq,nt)+1, W*bonds+atoms is the exact
    lexicographic objective. Dropping fragment constraints gives a superset of
    every original feasible witness when the original minimum is at least two.
    """
    started = time.monotonic()

    def recorded(value):
        return value if type(value) in (str, int, bool, type(None)) or type(value) is float and math.isfinite(value) else repr(value)

    result = {"status": "input_rejected", "optimal": False, "relaxation": True,
              "original_contract_optimal": False, "original_contract_feasible": False,
              "mapping": {}, "bond_mapping": [], "atom_count": 0, "bond_count": 0,
              "objective": "MCES bonds then atoms", "chemistry": "element equality; ANY bond order; flexible aromaticity",
              "original_contract": {"connectivity": "Fragments", "min_fragment_atoms": recorded(min_fragment_atoms),
                                    "max_fragments": recorded(max_fragments)},
              "relaxed_constraints": ["minimum fragment size", "maximum fragment count"],
              "original_contract_lower_bound": {"bond_count": 0, "atom_count": 0},
              "limits": {"seconds": recorded(time_limit), "atoms_per_graph": recorded(max_atoms),
                         "variables": recorded(max_variables), "rows": recorded(max_rows),
                         "nonzeros": recorded(max_nonzeros), "edge_checks": recorded(max_edge_checks),
                         "solver_nodes": recorded(node_limit), "solver_threads": recorded(solver_threads)},
              "model": {"variables": 0, "rows": 0, "nonzeros": 0, "edge_checks": 0},
              "formal_certificate": False, "integer_tolerance": 1e-7, "bound_rounding_tolerance": 1e-6}

    def finish(status, message):
        result.update(status=status, message=str(message), elapsed_seconds=time.monotonic() - started)
        return result

    try:
        if type(time_limit) not in (int, float) or not math.isfinite(time_limit) or time_limit <= 0:
            raise ValueError("time_limit must be finite and positive")
        if (type(min_fragment_atoms) is not int or min_fragment_atoms < 2
                or max_fragments is not None and (type(max_fragments) is not int or max_fragments <= 0)):
            raise ValueError("Non-isolated MCES requires min_fragment_atoms >= 2 and a positive fragment cap")
        if any(type(value) is not int or value <= 0 for value in
               (max_atoms, max_variables, max_rows, max_nonzeros, max_edge_checks, node_limit, solver_threads)):
            raise ValueError("Model and solver limits must be positive integers")
        if solver_threads > 2147483647:
            raise ValueError("solver_threads exceeds the HiGHS integer option range")
        query, target = _Graph(query_document, max_atoms), _Graph(target_document, max_atoms)
    except (KeyError, TypeError, ValueError) as error:
        return finish("input_rejected", error)
    deadline = started + time_limit
    scale = min(query.n, target.n) + 1
    maximum_bonds = min(len(query.bonds), len(target.bonds))
    trivial_upper = scale * maximum_bonds + min(query.n, target.n)
    result.update(bond_scale=scale, upper_bound={"scalar": trivial_upper, "bond_count": maximum_bonds,
                                                "source": "graph counts", "solver_reported": False})

    def clock():
        if time.monotonic() >= deadline:
            raise _Stopped("time_limit")

    costs, row_indices, column_indices, coefficients, lower_rows, upper_rows = [], [], [], [], [], []

    def variable(cost):
        clock()
        if len(costs) >= max_variables:
            raise _Stopped("variable_limit")
        index = len(costs)
        costs.append(cost)
        result["model"]["variables"] = len(costs)
        return index

    def row(terms, upper=math.inf):
        clock()
        if len(lower_rows) >= max_rows:
            raise _Stopped("row_limit")
        if len(coefficients) + len(terms) > max_nonzeros:
            raise _Stopped("nonzero_limit")
        index = len(lower_rows)
        for column, coefficient in terms:
            row_indices.append(index)
            column_indices.append(column)
            coefficients.append(coefficient)
        lower_rows.append(-math.inf)
        upper_rows.append(upper)
        result["model"].update(rows=len(lower_rows), nonzeros=len(coefficients))

    atom_pairs, atom_index = [], {}
    query_groups, target_groups = [[] for _ in range(query.n)], [[] for _ in range(target.n)]
    edge_pairs, incident = [], {}
    query_edges = {edge: [] for edge in query.bonds}
    target_edges = {edge: [] for edge in target.bonds}
    try:
        for query_atom in range(query.n):
            for target_atom in range(target.n):
                clock()
                if query.elements[query_atom] == target.elements[target_atom]:
                    index = variable(-1)
                    atom_pairs.append((query_atom, target_atom))
                    atom_index[query_atom, target_atom] = index
                    incident[index] = []
                    query_groups[query_atom].append(index)
                    target_groups[target_atom].append(index)
        for group in query_groups + target_groups:
            if group:
                row([(index, 1) for index in group], upper=1)
        for source in sorted(query.bonds):
            first, second = source
            for image in sorted(target.bonds):
                if result["model"]["edge_checks"] >= max_edge_checks:
                    raise _Stopped("edge_check_limit")
                result["model"]["edge_checks"] += 1
                clock()
                for image_first, image_second in (image, image[::-1]):
                    endpoints = atom_index.get((first, image_first)), atom_index.get((second, image_second))
                    if None in endpoints:
                        continue
                    index = variable(-scale)
                    edge_pairs.append((index, (first, second, image_first, image_second)))
                    query_edges[source].append(index)
                    target_edges[image].append(index)
                    for endpoint in endpoints:
                        incident[endpoint].append(index)
                        row([(index, 1), (endpoint, -1)], upper=0)
        for group in list(query_edges.values()) + list(target_edges.values()):
            if group:
                row([(index, 1) for index in group], upper=1)
        for atom, edges in incident.items():
            row([(atom, 1)] + [(edge, -1) for edge in edges], upper=0)
        clock()
    except _Stopped as stopped:
        return finish("model_limit", stopped.reason)
    if not edge_pairs:
        result.update(optimal=True, original_contract_optimal=True, original_contract_feasible=True,
                      witness_checked=True, fragments=[], upper_bound={"scalar": 0, "bond_count": 0,
                                                                      "source": "no compatible edge", "solver_reported": False})
        return finish("optimal", "No compatible edge exists; non-isolated witnesses are empty")
    try:
        import numpy as np
        import scipy
        from scipy.optimize import Bounds, LinearConstraint, milp
        from scipy.sparse import coo_matrix
    except ImportError as error:
        return finish("dependency_unavailable", error)
    result.update(scipy_version=scipy.__version__, numpy_version=np.__version__)
    if deadline <= time.monotonic():
        return finish("model_limit", "time_limit")
    matrix = coo_matrix((coefficients, (np.asarray(row_indices, dtype=np.int32),
                                       np.asarray(column_indices, dtype=np.int32))),
                        shape=(len(lower_rows), len(costs))).tocsc()
    remaining = deadline - time.monotonic()
    if remaining <= 0:
        return finish("model_limit", "time_limit")
    result["solver_options"] = {"time_limit": remaining, "node_limit": node_limit,
                                "mip_rel_gap": 0.0, "presolve": True, "threads": solver_threads}
    result["solver_warnings"] = []
    try:
        with warnings.catch_warnings(record=True) as retained_warnings:
            warnings.simplefilter("always")
            try:
                solved = milp(np.asarray(costs, dtype=float), integrality=np.ones(len(costs), dtype=np.int32),
                              bounds=Bounds(np.zeros(len(costs)), np.ones(len(costs))),
                              constraints=LinearConstraint(matrix, lower_rows, upper_rows),
                              options=result["solver_options"].copy())
            finally:
                result["solver_warnings"] = [{"category": item.category.__name__, "message": str(item.message)}
                                             for item in retained_warnings]
    except (TypeError, ValueError, RuntimeError) as error:
        return finish("solver_error", error)
    result.update(solver_status=int(solved.status), solver_message=str(solved.message))
    for name in ("mip_gap", "mip_dual_bound", "mip_node_count", "fun"):
        value = getattr(solved, name, None)
        result[name] = (None if value is None or not math.isfinite(float(value))
                        else int(value) if name == "mip_node_count" else float(value))
    if result["mip_dual_bound"] is not None:
        if result["mip_dual_bound"] > 1e-6:
            return finish("solver_error", "Solver bound contradicts the feasible empty witness")
        scalar = min(trivial_upper, max(0, math.ceil(-result["mip_dual_bound"] + 1e-6)))
        result["upper_bound"] = {"scalar": scalar, "bond_count": min(maximum_bonds, scalar // scale),
                                 "source": "HiGHS numerical dual bound rounded upwards", "solver_reported": True}
    if getattr(solved, "x", None) is None:
        return finish("limit_reached" if solved.status == 1 else "solver_error", solved.message)
    try:
        values = solved.x
        if len(values) != len(costs) or any(not math.isfinite(value)
                or min(abs(value), abs(value - 1.0)) > 1e-7 for value in values):
            raise ValueError("MCES solver values are not binary within tolerance")
        pairs = [pair for index, pair in enumerate(atom_pairs) if values[index] > 0.5]
        edges = [edge for index, edge in edge_pairs if values[index] > 0.5]
        mapping, fragments = _check_witness(query, target, pairs, edges)
        scalar = scale * len(edges) + len(mapping)
        if (result["fun"] is None or abs(result["fun"] + scalar) > 1e-6
                or result["upper_bound"]["scalar"] < scalar):
            raise ValueError("MCES objective or dual bound disagrees with the integer witness")
    except ValueError as error:
        result["upper_bound"] = {"scalar": trivial_upper, "bond_count": maximum_bonds,
                                 "source": "graph counts; inconsistent solver evidence rejected", "solver_reported": False}
        return finish("invalid_witness", error)
    feasible = all(len(fragment) >= min_fragment_atoms for fragment in fragments)
    feasible = feasible and (max_fragments is None or len(fragments) <= max_fragments)
    result.update(mapping=mapping, bond_mapping=edges, fragments=fragments, atom_count=len(mapping),
                  bond_count=len(edges), integer_objective=-scalar, witness_checked=True,
                  original_contract_feasible=feasible)
    if feasible:
        result["original_contract_lower_bound"] = {"bond_count": len(edges), "atom_count": len(mapping)}
    if solved.status == 0 and solved.success:
        if result["mip_dual_bound"] is None or abs(result["mip_dual_bound"] + scalar) > 1e-6:
            return finish("solver_error", "Missing or inconsistent final solver bound")
        result.update(optimal=True, original_contract_optimal=feasible,
                      upper_bound={"scalar": scalar, "bond_count": len(edges),
                                   "source": "HiGHS optimum with checked integer witness", "solver_reported": True})
        return finish("optimal", "Relaxation optimum reported; original fragment feasibility is separate")
    return finish("limit_reached" if solved.status == 1 else "solver_error", solved.message)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("query", type=Path)
    parser.add_argument("target", type=Path)
    parser.add_argument("--seconds", type=float, default=15.0)
    parser.add_argument("--variables", type=int, default=50000)
    parser.add_argument("--rows", type=int, default=200000)
    parser.add_argument("--nonzeros", type=int, default=1000000)
    parser.add_argument("--edge-checks", type=int, default=2000000)
    parser.add_argument("--solver-nodes", type=int, default=100000)
    parser.add_argument("--solver-threads", type=int, default=1)
    parser.add_argument("--max-atoms", type=int, default=256)
    parser.add_argument("--min-fragment-atoms", type=int, default=3)
    parser.add_argument("--max-fragments", type=int)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args(argv)
    try:
        result = solve_mces_any_relaxation(json.loads(args.query.read_text()), json.loads(args.target.read_text()),
            time_limit=args.seconds, max_variables=args.variables, max_rows=args.rows,
            max_nonzeros=args.nonzeros, max_edge_checks=args.edge_checks, node_limit=args.solver_nodes,
            solver_threads=args.solver_threads,
            max_atoms=args.max_atoms, min_fragment_atoms=args.min_fragment_atoms, max_fragments=args.max_fragments)
    except (OSError, ValueError) as error:
        result = {"status": "input_rejected", "optimal": False, "relaxation": True, "message": str(error)}
    output = json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n"
    if args.output is None:
        print(output, end="")
    else:
        with args.output.open("x") as stream:
            stream.write(output)
    return 0 if result["optimal"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
