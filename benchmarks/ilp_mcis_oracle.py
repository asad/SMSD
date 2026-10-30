#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
# Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
# See the NOTICE file for attribution, trademark, and algorithm IP terms.
"""Independent bounded MILP oracle for default molecule MCIS.

Input graphs use the search benchmark JSON schema: ``atoms`` contain ``element``
and ``aromatic``; ``bonds`` contain [first, second, order, aromatic, ring]. The
default chemistry matches elements, uses strict bond order and flexible
aromaticity, and ignores charge, isotope, rings, stereo and tautomers. Connected
or unrestricted-fragment MCIS maximises atoms, then selected bonds. Other
objectives and chemistry options are outside this helper's contract.

SciPy is optional and imported only when solving. A successful result is a
HiGHS-reported optimum with an independently checked integer witness, not a
formally verified certificate. Limits and unsuccessful solver results remain
explicit. The time limit covers model construction and the solver cooperatively;
it is not a hard process or memory limit.

API references:
https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.milp.html
https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.LinearConstraint.html
"""

import argparse
import json
import math
from pathlib import Path
import time


class _Stopped(Exception):
    def __init__(self, reason):
        self.reason = reason


class _Graph:
    def __init__(self, document, max_atoms):
        atoms, bonds = document["atoms"], document["bonds"]
        if not isinstance(atoms, list) or not isinstance(bonds, list):
            raise ValueError("Graph atoms and bonds must be lists")
        if len(atoms) > max_atoms:
            raise ValueError("Graph exceeds max_atoms")
        self.elements, self.aromatic = [], []
        for atom in atoms:
            element, aromatic = atom["element"], atom.get("aromatic", False)
            if type(element) is not int or element < 0 or type(aromatic) is not bool:
                raise ValueError("Atoms require a non-negative integer element and Boolean aromatic flag")
            self.elements.append(element)
            self.aromatic.append(aromatic)
        self.n = len(atoms)
        self.bonds = {}
        self.neighbors = [[] for _ in atoms]
        if len(bonds) > max_atoms * (max_atoms - 1) // 2:
            raise ValueError("Graph has too many bonds")
        for entry in bonds:
            if not isinstance(entry, (list, tuple)) or len(entry) != 5:
                raise ValueError("Bonds require [first, second, order, aromatic, ring]")
            first, second, order, aromatic, ring = entry
            if any(type(value) is not int for value in (first, second, order)):
                raise ValueError("Bond endpoints and order must be integers")
            if not (0 <= first < self.n and 0 <= second < self.n) or first == second or not 1 <= order <= 4:
                raise ValueError("Invalid bond endpoints or order")
            if type(aromatic) is not bool or type(ring) is not bool:
                raise ValueError("Bond aromatic and ring flags must be Boolean")
            key = tuple(sorted((first, second)))
            if key in self.bonds:
                raise ValueError("Duplicate bond")
            self.bonds[key] = order, aromatic
            self.neighbors[first].append(second)
            self.neighbors[second].append(first)

    def bond(self, first, second):
        return self.bonds.get(tuple(sorted((first, second))))


def graph_from_smsd(graph):
    """Copy graph properties for an independent solve; no SMSD search is called."""
    elements, aromatic = graph.atomic_num, graph.aromatic
    return {
        "atoms": [{"element": int(elements[atom]), "aromatic": bool(aromatic[atom])}
                  for atom in range(graph.n)],
        "bonds": [[int(first), int(second), int(order), bool(graph.bond_aromatic(first, second)),
                   bool(graph.bond_in_ring(first, second))] for first, second, order in graph.bonds()],
    }


def _pair_compatible(query, target, first, second):
    query_first, target_first = first
    query_second, target_second = second
    source, image = query.bond(query_first, query_second), target.bond(target_first, target_second)
    if (source is None) != (image is None):
        return False
    if source is None:
        return True
    source_order, source_aromatic = source
    image_order, image_aromatic = image
    return (source_order == image_order
            or source_aromatic and image_aromatic
            or source_aromatic and image_order in (1, 2)
            or image_aromatic and source_order in (1, 2)
            or source_order in (1, 2) and image_order in (1, 2)
            and query.aromatic[query_first] and query.aromatic[query_second]
            and target.aromatic[target_first] and target.aromatic[target_second])


def _witness(query, target, pairs, values, connected, tolerance):
    mapping = {}
    for pair, value in zip(pairs, values):
        if not math.isfinite(value) or min(abs(value), abs(value - 1.0)) > tolerance:
            raise ValueError("Solver atom mapping is not binary within tolerance")
        if value > 0.5:
            query_atom, target_atom = pair
            if query_atom in mapping or target_atom in mapping.values():
                raise ValueError("Solver mapping is not injective")
            if query.elements[query_atom] != target.elements[target_atom]:
                raise ValueError("Solver mapped incompatible atoms")
            mapping[query_atom] = target_atom
    selected = sorted(mapping.items())
    for position, first in enumerate(selected):
        for second in selected[position + 1:]:
            if not _pair_compatible(query, target, first, second):
                raise ValueError("Solver mapping violates induced topology or bond compatibility")
    bonds = sum(first in mapping and second in mapping for first, second in query.bonds)
    if connected and mapping:
        reached = {min(mapping)}
        pending = list(reached)
        while pending:
            for neighbor in query.neighbors[pending.pop()]:
                if neighbor in mapping and neighbor not in reached:
                    reached.add(neighbor)
                    pending.append(neighbor)
        if reached != set(mapping):
            raise ValueError("Solver mapping is disconnected")
    return mapping, bonds


def solve_mcis(query_document, target_document, *, connected=True, time_limit=15.0,
               max_atoms=256, max_variables=10000, max_rows=200000,
               max_nonzeros=1000000, max_pair_checks=2000000, node_limit=100000):
    """Return a bounded solver result, preserving rejection and stopped outcomes.

    All model coefficients are integers. Atom correspondences and root variables
    are binary; selected query edges are linearised conjunctions. In connected
    mode, one selected root supplies one unit to every selected query atom over
    selected query edges. The empty mapping remains feasible.
    """
    started = time.monotonic()
    def recorded(value):
        return value if type(value) in (str, int, bool, type(None)) or type(value) is float and math.isfinite(value) else repr(value)

    result = {"status": "input_rejected", "optimal": False, "mapping": {},
              "atom_count": 0, "bond_count": 0, "connected": recorded(connected),
              "chemistry": "default", "objective": "MCIS atoms then bonds",
              "limits": {"seconds": recorded(time_limit), "atoms_per_graph": recorded(max_atoms),
                         "variables": recorded(max_variables), "rows": recorded(max_rows),
                         "nonzeros": recorded(max_nonzeros), "pair_checks": recorded(max_pair_checks),
                         "solver_nodes": recorded(node_limit)},
              "model": {"variables": 0, "rows": 0, "nonzeros": 0, "pair_checks": 0}}

    def finish(status, message):
        result.update(status=status, message=str(message), elapsed_seconds=time.monotonic() - started)
        return result

    try:
        if (type(connected) is not bool or type(time_limit) not in (int, float)
                or not math.isfinite(time_limit) or time_limit <= 0):
            raise ValueError("connected must be Boolean and time_limit finite and positive")
        if any(type(value) is not int or value <= 0 for value in
               (max_atoms, max_variables, max_rows, max_nonzeros, max_pair_checks, node_limit)):
            raise ValueError("Model and solver limits must be positive integers")
        query, target = _Graph(query_document, max_atoms), _Graph(target_document, max_atoms)
    except (KeyError, TypeError, ValueError) as error:
        return finish("input_rejected", error)
    deadline = started + time_limit

    def clock():
        if time.monotonic() >= deadline:
            raise _Stopped("time_limit")

    costs, integrality, upper = [], [], []
    row_indices, column_indices, coefficients, lower_rows, upper_rows = [], [], [], [], []

    def variable(cost=0, integer=True, ceiling=1):
        clock()
        if len(costs) >= max_variables:
            raise _Stopped("variable_limit")
        index = len(costs)
        costs.append(cost)
        integrality.append(int(integer))
        upper.append(ceiling)
        result["model"]["variables"] = len(costs)
        return index

    def row(terms, lower=-math.inf, upper=math.inf):
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
        lower_rows.append(lower)
        upper_rows.append(upper)
        result["model"].update(rows=len(lower_rows), nonzeros=len(coefficients))

    pairs, query_groups, target_groups = [], [[] for _ in range(query.n)], [[] for _ in range(target.n)]
    atom_scale = min(len(query.bonds), len(target.bonds)) + 1
    result["atom_scale"] = atom_scale
    try:
        for query_atom in range(query.n):
            clock()
            for target_atom in range(target.n):
                if query.elements[query_atom] == target.elements[target_atom]:
                    index = variable(-atom_scale)
                    pairs.append((query_atom, target_atom))
                    query_groups[query_atom].append(index)
                    target_groups[target_atom].append(index)
        if not pairs:
            clock()
            result.update(optimal=True, bound=[0, 0], witness_checked=True, formal_certificate=False)
            return finish("optimal", "No element-compatible atom pair exists")
        for group in query_groups + target_groups:
            if group:
                row([(index, 1) for index in group], upper=1)
        for position, first in enumerate(pairs):
            for second_position in range(position + 1, len(pairs)):
                if result["model"]["pair_checks"] >= max_pair_checks:
                    raise _Stopped("pair_check_limit")
                result["model"]["pair_checks"] += 1
                if (result["model"]["pair_checks"] & 1023) == 0:
                    clock()
                second = pairs[second_position]
                if first[0] != second[0] and first[1] != second[1] and not _pair_compatible(query, target, first, second):
                    row([(position, 1), (second_position, 1)], upper=1)
        edge_variables = {}
        for first, second in sorted(query.bonds):
            selected = variable(-1)
            edge_variables[first, second] = selected
            row([(selected, 1)] + [(index, -1) for index in query_groups[first]], upper=0)
            row([(selected, 1)] + [(index, -1) for index in query_groups[second]], upper=0)
            row([(selected, 1)] + [(index, -1) for index in query_groups[first] + query_groups[second]], lower=-1)
        if connected:
            capacity = min(query.n, target.n)
            roots = [variable() for _ in range(query.n)]
            supplies = [variable(integer=False, ceiling=capacity) for _ in range(query.n)]
            flows = {}
            for first, second in edge_variables:
                flows[first, second] = variable(integer=False, ceiling=capacity)
                flows[second, first] = variable(integer=False, ceiling=capacity)
            row([(root, 1) for root in roots], upper=1)
            row([(index, 1) for index in range(len(pairs))] + [(root, -capacity) for root in roots], upper=0)
            for atom in range(query.n):
                row([(roots[atom], 1)] + [(index, -1) for index in query_groups[atom]], upper=0)
                row([(supplies[atom], 1), (roots[atom], -capacity)], upper=0)
                terms = [(supplies[atom], 1)] + [(index, -1) for index in query_groups[atom]]
                for neighbor in query.neighbors[atom]:
                    terms.extend(((flows[neighbor, atom], 1), (flows[atom, neighbor], -1)))
                row(terms, lower=0, upper=0)
            for (first, second), flow in flows.items():
                selected = edge_variables[tuple(sorted((first, second)))]
                row([(flow, 1), (selected, -capacity)], upper=0)
        clock()
    except _Stopped as stopped:
        return finish("model_limit", stopped.reason)
    try:
        import numpy as np
        import scipy
        from scipy.optimize import Bounds, LinearConstraint, milp
        from scipy.sparse import coo_matrix
    except ImportError as error:
        return finish("dependency_unavailable", error)
    result["scipy_version"] = scipy.__version__
    result["numpy_version"] = np.__version__
    if deadline <= time.monotonic():
        return finish("model_limit", "time_limit")
    matrix = coo_matrix((coefficients, (np.asarray(row_indices, dtype=np.int32),
                                       np.asarray(column_indices, dtype=np.int32))),
                        shape=(len(lower_rows), len(costs))).tocsc()
    remaining = deadline - time.monotonic()
    if remaining <= 0:
        return finish("model_limit", "time_limit")
    try:
        solved = milp(np.asarray(costs, dtype=float), integrality=np.asarray(integrality),
                      bounds=Bounds(np.zeros(len(costs)), np.asarray(upper, dtype=float)),
                      constraints=LinearConstraint(matrix, lower_rows, upper_rows),
                      options={"time_limit": remaining, "node_limit": node_limit,
                               "mip_rel_gap": 0.0, "presolve": True})
    except (TypeError, ValueError, RuntimeError) as error:
        return finish("solver_error", error)
    result["solver_status"] = int(solved.status)
    result["solver_message"] = str(solved.message)
    for name in ("mip_gap", "mip_dual_bound", "mip_node_count", "fun"):
        value = getattr(solved, name, None)
        result[name] = (None if value is None or not math.isfinite(float(value))
                        else int(value) if name == "mip_node_count" else float(value))
    if getattr(solved, "x", None) is None:
        return finish("limit_reached" if solved.status == 1 else "solver_error", solved.message)
    try:
        mapping, bonds = _witness(query, target, pairs, solved.x[:len(pairs)], connected, 1e-7)
        objective = -atom_scale * len(mapping) - bonds
        if result["fun"] is None or abs(result["fun"] - objective) > 1e-6:
            raise ValueError("Solver objective differs from the checked integer witness")
    except ValueError as error:
        return finish("invalid_witness", error)
    result.update(mapping=mapping, atom_count=len(mapping), bond_count=bonds,
                  integer_objective=objective, atom_scale=atom_scale,
                  witness_checked=True, formal_certificate=False)
    if solved.status == 0 and solved.success:
        if result["mip_dual_bound"] is None or abs(result["mip_dual_bound"] - objective) > 1e-6:
            return finish("solver_error", "Missing or inconsistent final solver bound")
        result["optimal"] = True
        return finish("optimal", "HiGHS reported an optimum; integer witness independently checked")
    return finish("limit_reached" if solved.status == 1 else "solver_error", solved.message)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("query", type=Path, help="Query graph JSON")
    parser.add_argument("target", type=Path, help="Target graph JSON")
    parser.add_argument("--disconnected", action="store_true", help="Allow any number of singleton or larger fragments")
    parser.add_argument("--seconds", type=float, default=15.0)
    parser.add_argument("--variables", type=int, default=10000)
    parser.add_argument("--rows", type=int, default=200000)
    parser.add_argument("--nonzeros", type=int, default=1000000)
    parser.add_argument("--pair-checks", type=int, default=2000000)
    parser.add_argument("--solver-nodes", type=int, default=100000)
    parser.add_argument("--max-atoms", type=int, default=256)
    parser.add_argument("--output", type=Path, help="Write a new result JSON; an existing file is never replaced")
    args = parser.parse_args(argv)
    try:
        query = json.loads(args.query.read_text(encoding="utf-8"))
        target = json.loads(args.target.read_text(encoding="utf-8"))
        result = solve_mcis(query, target, connected=not args.disconnected, time_limit=args.seconds,
                            max_atoms=args.max_atoms, max_variables=args.variables, max_rows=args.rows,
                            max_nonzeros=args.nonzeros, max_pair_checks=args.pair_checks, node_limit=args.solver_nodes)
    except (OSError, ValueError) as error:
        result = {"status": "input_rejected", "optimal": False, "message": str(error)}
    output = json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n"
    if args.output is None:
        print(output, end="")
    else:
        with args.output.open("x", encoding="utf-8") as stream:
            stream.write(output)
    return 0 if result["optimal"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
