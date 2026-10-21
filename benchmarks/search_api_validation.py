#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
# Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
# See the NOTICE file for attribution, trademark, and algorithm IP terms.
"""Check search mappings and compare benchmark JSONL without toolkit dependencies.

Small MCIS/MCES proofs use exhaustive partial injections and selected-edge subsets.
SMARTS validation supports elements, atomic numbers, aromaticity, Boolean atom
expressions, charge/isotope, H/D/X/v, ring membership and recursive expressions.
Stereo and complete-ring predicates have explicit coverage warnings.
"""
import argparse
import csv
import itertools
import json
import math
import re
from collections import defaultdict


class UnsupportedValidation(ValueError):
    pass


def _edge(first, second):
    return min(first, second), max(first, second)


class Graph:
    def __init__(self, data):
        self.atoms = data["atoms"]
        self.bonds = {}
        self.neighbours = [[] for _ in self.atoms]
        for first, second, order, aromatic, ring in data["bonds"]:
            if first == second or not (0 <= first < len(self.atoms) and 0 <= second < len(self.atoms)):
                raise ValueError("Graph bond endpoint is invalid")
            key = _edge(first, second)
            if key in self.bonds or order <= 0:
                raise ValueError("Graph has a duplicate or invalid bond")
            self.bonds[key] = int(order), bool(aromatic), bool(ring)
            self.neighbours[first].append(second)
            self.neighbours[second].append(first)

    def bond(self, first, second):
        return self.bonds.get(_edge(first, second))


def _components(atoms, bonds):
    adjacent = {atom: set() for atom in atoms}
    for first, second in bonds:
        adjacent[first].add(second)
        adjacent[second].add(first)
    components = []
    while adjacent:
        pending = [min(adjacent)]
        visited = set(pending)
        while pending:
            for neighbour in adjacent[pending.pop()]:
                if neighbour not in visited:
                    visited.add(neighbour)
                    pending.append(neighbour)
        for atom in visited:
            adjacent.pop(atom)
        components.append(tuple(sorted(visited)))
    return tuple(sorted(components))


def _atom_matches(query, query_atom, target, target_atom, profile):
    first, second = query.atoms[query_atom], target.atoms[target_atom]
    if first["element"] != second["element"]:
        return False
    if profile in ("strict", "rings"):
        if any(first[key] != second[key] for key in ("charge", "aromatic", "ring")):
            return False
        if first["isotope"] and second["isotope"] and first["isotope"] != second["isotope"]:
            return False
    return True


def _bond_matches(query, first, second, target, image_first, image_second, profile):
    source, image = query.bond(first, second), target.bond(image_first, image_second)
    if source is None or image is None:
        return False
    if profile in ("strict", "rings") and source[1:] != image[1:]:
        return False
    if profile in ("any", "tautomer") or source[0] == image[0]:
        return True
    if profile not in ("strict", "rings"):
        if source[1] and image[0] in (1, 2, 4) or image[1] and source[0] in (1, 2, 4):
            return True
        if (query.atoms[first]["aromatic"] and query.atoms[second]["aromatic"]
                and target.atoms[image_first]["aromatic"] and target.atoms[image_second]["aromatic"]
                and source[0] in (1, 2) and image[0] in (1, 2)):
            return True
    return False


def _split(expression, delimiter):
    depth = 0
    start = 0
    parts = []
    for index, character in enumerate(expression):
        if character in "([":
            depth += 1
        elif character in ")]":
            depth -= 1
        elif character == delimiter and depth == 0:
            parts.append(expression[start:index])
            start = index + 1
    return parts + [expression[start:]]


_ELEMENTS = {"H": 1, "B": 5, "C": 6, "N": 7, "O": 8, "F": 9, "Na": 11,
             "Mg": 12, "Al": 13, "Si": 14, "P": 15, "S": 16, "Cl": 17,
             "Ca": 20, "Fe": 26, "Cu": 29, "Zn": 30, "As": 33, "Se": 34,
             "Br": 35, "I": 53}


def _atom_expression(expression):
    for delimiter, conjunction in ((";", True), (",", False), ("&", True)):
        parts = _split(expression, delimiter)
        if len(parts) > 1:
            predicates = [_atom_expression(part) for part in parts]
            return lambda graph, atom: (all if conjunction else any)(test(graph, atom) for test in predicates)
    predicates = []
    offset = 0
    while offset < len(expression):
        negate = expression[offset] == "!"
        if negate:
            offset += 1
        remaining = expression[offset:]
        if remaining.startswith("$("):
            depth, end = 1, 2
            while end < len(remaining) and depth:
                depth += (remaining[end] == "(") - (remaining[end] == ")")
                end += 1
            if depth:
                raise UnsupportedValidation("Unclosed recursive SMARTS")
            pattern = Pattern(remaining[2:end - 1])
            test = lambda graph, atom, pattern=pattern: bool(_embeddings(pattern, graph, False, anchor=atom, limit=1))
        elif remaining.startswith("@@") or remaining.startswith("@"):
            end = 2 if remaining.startswith("@@") else 1
            test = lambda graph, atom: True
        else:
            token = re.match(r"#[0-9]+|[0-9]+|[+-]+[0-9]*|:[0-9]+|[HDXvRr][0-9]*|Cl|Br|Si|Na|Mg|Al|Ca|Fe|Cu|Zn|As|Se|[A-Za-z*]", remaining)
            if token is None:
                raise UnsupportedValidation("Unsupported SMARTS atom primitive: " + remaining)
            text = token.group()
            end = len(text)
            if text == "H" and expression in ("H", "H+", "H-"):
                test = lambda graph, atom: graph.atoms[atom]["element"] == 1
            elif text.startswith("#"):
                number = int(text[1:])
                test = lambda graph, atom, number=number: graph.atoms[atom]["element"] == number
            elif text.isdigit():
                number = int(text)
                test = lambda graph, atom, number=number: graph.atoms[atom]["isotope"] == number
            elif text[0] in "+-":
                sign = 1 if text[0] == "+" else -1
                digits = re.search(r"[0-9]+$", text)
                number = sign * (int(digits.group()) if digits else len(text))
                test = lambda graph, atom, number=number: graph.atoms[atom]["charge"] == number
            elif text.startswith(":") or text == "*":
                test = lambda graph, atom: True
            elif text[0] in "HDXvRr" and (len(text) == 1 or text[1:].isdigit()):
                primitive = text[0]
                number = int(text[1:]) if len(text) > 1 else 1
                if primitive in "Rr" and len(text) > 1 and number != 0:
                    raise UnsupportedValidation("SMARTS ring-size/count primitives are outside the validator subset")
                def test(graph, atom, primitive=primitive, number=number):
                    value = graph.atoms[atom]
                    if primitive == "H":
                        actual = value["hydrogens"] + sum(graph.atoms[neighbour]["element"] == 1 for neighbour in graph.neighbours[atom])
                    elif primitive == "D":
                        actual = len(graph.neighbours[atom])
                    elif primitive == "X":
                        actual = len(graph.neighbours[atom]) + value["hydrogens"]
                    elif primitive == "v":
                        if value["aromatic"]:
                            raise UnsupportedValidation("Aromatic SMARTS valence is outside the validator subset")
                        actual = value["hydrogens"] + sum(graph.bond(atom, neighbour)[0] for neighbour in graph.neighbours[atom])
                    else:
                        return bool(value["ring"]) == (number != 0)
                    return actual == number
            elif text in ("a", "A"):
                test = lambda graph, atom, aromatic=text == "a": bool(graph.atoms[atom]["aromatic"]) == aromatic
            else:
                symbol = text.capitalize()
                if symbol not in _ELEMENTS:
                    raise UnsupportedValidation("Unsupported SMARTS element: " + text)
                number, aromatic = _ELEMENTS[symbol], text.islower()
                test = lambda graph, atom, number=number, aromatic=aromatic: (graph.atoms[atom]["element"] == number
                    and bool(graph.atoms[atom]["aromatic"]) == aromatic)
        offset += end
        if negate:
            predicate = test
            test = lambda graph, atom, predicate=predicate: not predicate(graph, atom)
        predicates.append(test)
    return lambda graph, atom: all(test(graph, atom) for test in predicates)


class Pattern:
    def __init__(self, expression):
        self.atoms, self.bonds = [], {}
        self.stereo = bool(re.search(r"[@/\\]", expression))
        current, branches, rings, pending, offset = None, [], {}, None, 0
        while offset < len(expression):
            character = expression[offset]
            if character == "(":
                if current is None:
                    raise UnsupportedValidation("SMARTS component grouping is outside the validator subset")
                branches.append(current)
                offset += 1
                continue
            if character == ")":
                current = branches.pop()
                offset += 1
                continue
            if character == ".":
                current, pending = None, None
                offset += 1
                continue
            if character in "-=#~:/\\":
                pending = character
                offset += 1
                continue
            if character.isdigit():
                if character in rings:
                    first, previous = rings.pop(character)
                    self.bonds[_edge(first, current)] = pending or previous or ""
                else:
                    rings[character] = current, pending
                pending = None
                offset += 1
                continue
            if character == "[":
                depth, end = 1, offset + 1
                while end < len(expression) and depth:
                    depth += (expression[end] == "[") - (expression[end] == "]")
                    end += 1
                if depth:
                    raise UnsupportedValidation("Unclosed SMARTS atom")
                atom = expression[offset + 1:end - 1]
            else:
                token = re.match(r"Cl|Br|[BCNOPSFIbcnops*]", expression[offset:])
                if token is None:
                    raise UnsupportedValidation("Unsupported SMARTS syntax: " + expression[offset:])
                atom, end = token.group(), offset + len(token.group())
            index = len(self.atoms)
            self.atoms.append(_atom_expression(atom))
            if current is not None:
                self.bonds[_edge(current, index)] = pending or ""
            current, pending, offset = index, None, end
        if branches or rings or not self.atoms:
            raise UnsupportedValidation("Incomplete SMARTS pattern")

    def bond_matches(self, first, second, target, image_first, image_second):
        bond = target.bond(image_first, image_second)
        if bond is None:
            return False
        token = self.bonds[_edge(first, second)]
        return (token == "~" or token == ":" and bond[1]
                or token == "" and (bond[0] == 1 or bond[1])
                or token in ("-", "/", "\\") and bond[0] == 1 and not bond[1]
                or token == "=" and bond[0] == 2 and not bond[1]
                or token == "#" and bond[0] == 3 and not bond[1])


def _embeddings(pattern, target, induced, anchor=None, limit=None):
    if len(pattern.atoms) > len(target.atoms):
        return set()
    candidates = [[atom for atom in range(len(target.atoms)) if predicate(target, atom)] for predicate in pattern.atoms]
    if anchor is not None:
        candidates[0] = [anchor] if anchor in candidates[0] else []
    order = sorted(range(len(candidates)), key=lambda atom: len(candidates[atom]))
    found, mapping, used = set(), {}, set()
    visits = 0
    def visit(depth):
        nonlocal visits
        visits += 1
        if visits > 200000:
            raise UnsupportedValidation("Independent embedding enumeration exceeded its validation bound")
        if depth == len(order):
            found.add(tuple(sorted(mapping.items())))
            return
        query_atom = order[depth]
        for target_atom in candidates[query_atom]:
            if target_atom in used:
                continue
            compatible = True
            for previous, image in mapping.items():
                bond = _edge(query_atom, previous) in pattern.bonds
                if ((bond and not pattern.bond_matches(query_atom, previous, target, target_atom, image))
                        or (induced and not bond and target.bond(target_atom, image))):
                    compatible = False
                    break
            if compatible:
                mapping[query_atom] = target_atom
                used.add(target_atom)
                visit(depth + 1)
                used.remove(target_atom)
                del mapping[query_atom]
                if limit is not None and len(found) >= limit:
                    return
    visit(0)
    return found


def _weights(case):
    value = case.get("weights", "")
    return [float(weight) for weight in value.split(",")] if isinstance(value, str) and value else list(value or [])


def _score(case, atoms, bonds):
    if case["operation"] != "mcs":
        return len(atoms), len(atoms), 0
    if case["objective"] == "MCES":
        return len(bonds), len(atoms), 0
    weights = _weights(case)
    primary = int(sum(weights[atom] for atom in sorted(atoms)) * 1000) if weights else len(atoms)
    return primary, len(atoms), len(bonds)


def _signature(mapping, bonds):
    return tuple(sorted(mapping.items())), tuple(sorted(_edge(*bond) for bond in bonds))


def graph_index_map(source_data, reference_data, query_weights=()):
    """Preserve emitted graph properties and any parser-indexed query weights."""
    source, reference = Graph(source_data), Graph(reference_data)
    if len(source.atoms) != len(reference.atoms):
        raise ValueError("Parsed atom counts differ")
    if query_weights and (len(query_weights) != len(source.atoms)
                          or not all(math.isfinite(weight) for weight in query_weights)):
        raise ValueError("Finite query weights must cover every parsed query atom")
    if source.atoms == reference.atoms and source.bonds == reference.bonds:
        return dict(enumerate(range(len(source.atoms))))
    fields = ("element", "charge", "isotope", "hydrogens", "aromatic", "ring")
    graphs = source, reference
    labels = [[tuple(atom[field] for field in fields) + ((query_weights[index],) if query_weights else ())
               for index, atom in enumerate(graph.atoms)] for graph in graphs]
    ranks = {value: rank for rank, value in enumerate(sorted(set(labels[0] + labels[1])))}
    labels = [[ranks[value] for value in graph_labels] for graph_labels in labels]
    for _ in range(len(source.atoms) + 1):
        signatures = [[(labels[index][atom], tuple(sorted((graph.bond(atom, neighbour), labels[index][neighbour])
                       for neighbour in graph.neighbours[atom]))) for atom in range(len(graph.atoms))]
                      for index, graph in enumerate(graphs)]
        ranks = {value: rank for rank, value in enumerate(sorted(set(signatures[0] + signatures[1])))}
        refined = [[ranks[value] for value in values] for values in signatures]
        if refined == labels:
            break
        labels = refined
    if sorted(labels[0]) != sorted(labels[1]):
        raise ValueError("Parsed atom, bond or query-weight properties differ" if query_weights
                         else "Parsed atom or bond properties differ")
    candidates = [[image for image in range(len(reference.atoms)) if labels[0][atom] == labels[1][image]]
                  for atom in range(len(source.atoms))]
    order = sorted(range(len(source.atoms)), key=lambda atom: (len(candidates[atom]), -len(source.neighbours[atom]), atom))
    mapping, used, visits = {}, set(), 0
    def assign(depth):
        nonlocal visits
        visits += 1
        if visits > 200000:
            raise UnsupportedValidation("Parser-index normalisation exceeded its validation bound")
        if depth == len(order):
            return True
        atom = order[depth]
        for image in candidates[atom]:
            if image in used or any(source.bond(atom, previous) != reference.bond(image, target)
                                    for previous, target in mapping.items()):
                continue
            mapping[atom] = image
            used.add(image)
            if assign(depth + 1):
                return True
            used.remove(image)
            del mapping[atom]
        return False
    if not assign(0):
        raise ValueError("Parsed molecular graphs are not property-preserving isomorphs")
    return mapping


def normalised_matches(matches, query_indices, target_indices):
    result = set()
    for match in matches:
        mapping = {query_indices.get(atom, atom): target_indices[image] for atom, image in match["atoms"]}
        bonds = [_edge(query_indices.get(first, first), query_indices.get(second, second))
                 for first, second, _, _ in match["bonds"]]
        result.add(_signature(mapping, bonds))
    return result


def exhaustive_mcs(case, query_data, target_data):
    """Independent unpruned oracle for graphs with at most five atoms each."""
    query, target = Graph(query_data), Graph(target_data)
    if max(len(query.atoms), len(target.atoms)) > 5 or case["profile"] in ("stereo", "rings"):
        raise UnsupportedValidation("Independent MCS oracle covers at most five atoms without stereo/ring completeness")
    optimum, matches = (0, 0, 0), set()
    for size in range(1, min(len(query.atoms), len(target.atoms)) + 1):
        for atoms in itertools.combinations(range(len(query.atoms)), size):
            for images in itertools.permutations(range(len(target.atoms)), size):
                mapping = dict(zip(atoms, images))
                if not all(_atom_matches(query, atom, target, image, case["profile"]) for atom, image in mapping.items()):
                    continue
                common, valid = [], True
                for first, second in itertools.combinations(atoms, 2):
                    source, image = query.bond(first, second), target.bond(mapping[first], mapping[second])
                    compatible = source and image and _bond_matches(query, first, second, target, mapping[first], mapping[second], case["profile"])
                    if case["objective"] == "MCIS" and (bool(source) != bool(image) or source and not compatible):
                        valid = False
                        break
                    if compatible:
                        common.append((first, second))
                if not valid:
                    continue
                selections = [common] if case["objective"] == "MCIS" else (
                    [bond for index, bond in enumerate(common) if mask & (1 << index)] for mask in range(1, 1 << len(common)))
                for bonds in selections:
                    if case["objective"] == "MCES" and {atom for bond in bonds for atom in bond} != set(atoms):
                        continue
                    components = _components(atoms, bonds)
                    if case["connectivity"] == "Connected" and len(components) != 1:
                        continue
                    score = _score(case, atoms, bonds)
                    if score > optimum:
                        optimum, matches = score, set()
                    if score == optimum:
                        matches.add(_signature(mapping, bonds))
    return optimum, matches


def validation_warnings(case, query_graph, target_graph):
    warnings = []
    if case["profile"] == "stereo" or case["operation"] == "smarts" and re.search(r"[@/\\]", case["query"]):
        warnings.append("Stereo predicates are not independently certified by the graph metadata")
    if case["profile"] == "rings":
        warnings.append("Complete-ring coverage is not independently certified")
    if case["operation"] == "smarts":
        try:
            Pattern(case["query"])
        except UnsupportedValidation as error:
            warnings.append(str(error))
        if re.search(r"\bv[0-9]", case["query"]) and any(atom["aromatic"] for atom in target_graph["atoms"]):
            warnings.append("Aromatic SMARTS valence is outside the validator subset")
    elif query_graph is not None and target_graph is not None:
        size = max(len(query_graph["atoms"]), len(target_graph["atoms"]))
        if case["operation"] == "mcs" and size > 5:
            warnings.append("Independent optimum oracle is restricted to graphs of at most five atoms")
        elif case["operation"] == "sub" and size > 8:
            warnings.append("Independent full-embedding oracle is restricted to graphs of at most eight atoms")
    return warnings


def validate_result(case, result, query_graph, target_graph):
    """Return errors only; incomplete searches are valid, separately labelled outcomes."""
    errors = []
    try:
        target = Graph(target_graph)
        pattern = None
        if case["operation"] == "smarts":
            try:
                pattern = Pattern(case["query"])
            except UnsupportedValidation:
                pass
            query = None
            query_atoms = len(pattern.atoms) if pattern else None
        else:
            query = Graph(query_graph)
            query_atoms = len(query.atoms)
        matches = result["matches"]
        status = result["status"]
        if status not in ("Completed", "TimedOut", "Cancelled", "LimitReached"):
            errors.append("Unknown stopping status")
        if result["complete"] and (status != "Completed" or not result["optimal"]):
            errors.append("Complete result lacks a completed optimal proof")
        if status == "Completed" and not result["complete"]:
            errors.append("Completed exact search lacks its completion proof")
        if str(case.get("timeout_ms", "")) == "0" and (status != "TimedOut" or matches
                or result["optimal"] or result["complete"] or result["nodes_visited"] != 0):
            errors.append("Zero timeout did not return the immediate unfinished outcome")
        if not math.isfinite(result["score"]) or result["nodes_visited"] < 0:
            errors.append("Invalid score or node count")
        limit = int(case["limit"])
        if len(matches) > (limit if limit else 1):
            errors.append("Result exceeds the requested mapping limit")
        signatures, scores = set(), []
        for match in matches:
            pairs = match["atoms"]
            mapping = dict(pairs)
            if len(mapping) != len(pairs) or len(set(mapping.values())) != len(mapping):
                errors.append("Atom mapping is not injective")
                continue
            if any(atom < 0 or query_atoms is not None and atom >= query_atoms or image < 0 or image >= len(target.atoms)
                   for atom, image in mapping.items()):
                errors.append("Atom mapping index is out of range")
                continue
            if case["operation"] != "mcs" and query_atoms is not None and set(mapping) != set(range(query_atoms)):
                errors.append("Substructure result does not cover the full query")
            for atom, image in mapping.items():
                if query and not _atom_matches(query, atom, target, image, case["profile"]):
                    errors.append("Mapped atoms violate the chemical predicates")
                if pattern:
                    try:
                        if not pattern.atoms[atom](target, image):
                            errors.append("Mapped atoms violate the SMARTS predicates")
                    except UnsupportedValidation:
                        pass
            bonds = set()
            for first, second, image_first, image_second in match["bonds"]:
                key = _edge(first, second)
                if key in bonds or mapping.get(first) != image_first or mapping.get(second) != image_second:
                    errors.append("Selected bond is duplicate or inconsistent with the atom map")
                    continue
                bonds.add(key)
                if (not target.bond(image_first, image_second)
                        or query and not _bond_matches(query, first, second, target, image_first, image_second, case["profile"])
                        or pattern and (key not in pattern.bonds or not pattern.bond_matches(first, second, target, image_first, image_second))):
                    errors.append("Selected bond violates graph adjacency or bond predicates")
            components = _components(mapping, bonds)
            supplied = tuple(sorted(tuple(sorted(component)) for component in match["fragments"]))
            if components != supplied:
                errors.append("Reported fragments differ from selected-edge components")
            if case["operation"] == "mcs":
                if case["connectivity"] == "Connected" and len(components) != 1:
                    errors.append("Connected MCS contains multiple components")
                if case["objective"] == "MCES" and {atom for bond in bonds for atom in bond} != set(mapping):
                    errors.append("MCES contains isolated mapped atoms")
            induced = case["operation"] == "mcs" and case["objective"] == "MCIS" or case["operation"] != "mcs" and int(case["induced"])
            for first, second in itertools.combinations(mapping, 2):
                source = query.bond(first, second) if query else _edge(first, second) in pattern.bonds if pattern else None
                image = target.bond(mapping[first], mapping[second])
                if induced and bool(source) != bool(image):
                    errors.append("Mapping violates induced adjacency")
                if (case["operation"] != "mcs" or case["objective"] == "MCIS") and source and _edge(first, second) not in bonds:
                    errors.append("Required query bond is absent from the selected graph")
            signature = _signature(mapping, bonds)
            if signature in signatures:
                errors.append("Duplicate retained match")
            signatures.add(signature)
            scores.append(_score(case, mapping, bonds))
        if scores and len(set(scores)) != 1:
            errors.append("Retained incumbent matches have different objective tuples")
        if scores:
            score = scores[0][0] / (1000 if _weights(case) else 1)
            if result["score"] != score:
                errors.append("Reported score differs from selected graph objective")
        elif result["score"] != 0:
            errors.append("Empty result has a non-zero score")
        if case["operation"] == "mcs":
            try:
                optimum, expected = exhaustive_mcs(case, query_graph, target_graph)
            except UnsupportedValidation:
                expected = None
        elif pattern and not pattern.stereo:
            expected_maps = _embeddings(pattern, target, bool(int(case["induced"])))
            expected = {_signature(dict(mapping), pattern.bonds) for mapping in expected_maps}
            optimum = (len(pattern.atoms), len(pattern.atoms), 0) if expected else (0, 0, 0)
        elif query and case["profile"] != "stereo" and max(len(query.atoms), len(target.atoms)) <= 8:
            class MoleculePattern:
                atoms = [lambda target_graph, image, atom=atom: _atom_matches(query, atom, target_graph, image, case["profile"])
                         for atom in range(len(query.atoms))]
                bonds = query.bonds
                @staticmethod
                def bond_matches(first, second, graph, image_first, image_second):
                    return _bond_matches(query, first, second, graph, image_first, image_second, case["profile"])
            expected_maps = _embeddings(MoleculePattern(), target, bool(int(case["induced"])))
            expected = {_signature(dict(mapping), query.bonds) for mapping in expected_maps}
            optimum = (len(query.atoms), len(query.atoms), 0) if expected else (0, 0, 0)
        else:
            expected = None
        if expected is not None:
            if result["optimal"] and (scores[0] if scores else (0, 0, 0)) != optimum:
                errors.append("Claimed optimum differs from the independent exhaustive oracle")
            if result["optimal"] and not signatures <= expected:
                errors.append("Retained match is outside the independently optimal mapping set")
            if result["complete"] and limit and signatures != expected:
                errors.append("Complete enumeration differs from the independent mapping set")
    except UnsupportedValidation:
        # Unsupported predicates are exposed separately by validation_warnings.
        pass
    except (ValueError, TypeError, KeyError, IndexError) as error:
        errors.append("Malformed search result: " + str(error))
    return sorted(set(errors))


def compare_jsonl(cases_path, result_paths):
    with open(cases_path, encoding="utf-8", newline="") as source:
        cases = {row["case_id"]: row for row in csv.DictReader(source, delimiter="\t")}
    rows, graphs, errors, warnings = [], {}, [], []
    references, index_maps, graph_differences = {}, {}, []
    # Batches of distinct targets record the graph of each target once, on their first sample.
    batch_graphs = {}
    for input_index, path in enumerate(result_paths):
        with open(path, encoding="utf-8") as source:
            for line in source:
                row = json.loads(line)
                rows.append((row, input_index))
                if row.get("target_graph") is not None:
                    key = row["case_id"], row["language"], input_index
                    graphs[key] = row.get("query_graph"), row["target_graph"]
                    references.setdefault(row["case_id"], graphs[key])
                if row.get("batch_target_graphs") is not None:
                    batch_graphs[row["case_id"], row["language"], input_index] = row["batch_target_graphs"]
    for key, (query, target) in graphs.items():
        reference_query, reference_target = references[key[0]]
        try:
            query_weights = _weights(cases[key[0]])
            index_maps[key] = (graph_index_map(query, reference_query, query_weights) if query is not None else {},
                               graph_index_map(target, reference_target))
        except UnsupportedValidation as error:
            warnings.append({"case_id": key[0], "language": key[1], "warning": str(error)})
        except ValueError as error:
            graph_differences.append({"case_id": key[0], "language": key[1], "error": str(error)})
    paired = defaultdict(list)
    for row, input_index in rows:
        case = cases[row["case_id"]]
        graph_key = row["case_id"], row["language"], input_index
        query, target = graphs.get(graph_key, (None, None))
        identity = {key: row[key] for key in ("case_id", "language", "phase", "sample")}
        if row["phase"] == "prepare":
            if row["results"]:
                errors.append({**identity, "error": "Preparation row contains a search result"})
            continue
        if target is None:
            errors.append({**identity, "error": "Missing parsed graph metadata"})
            continue
        expected_count = int(case["batch_size"]) if row["phase"] == "batch" else 1
        if len(row["results"]) != expected_count:
            errors.append({**identity, "error": "Result cardinality differs from requested calls"})
        distinct = row.get("batch_targets") == "distinct"
        targets = batch_graphs.get(graph_key) if distinct else None
        if distinct and (targets is None or len(targets) < len(row["results"])):
            errors.append({**identity, "error": "Missing batch target graph metadata"})
            continue
        for index, result in enumerate(row["results"]):
            # A distinct target is checked against its own graph; mappings are compared across
            # languages only for the case target, whose parser numbering is known.
            searched = targets[index] if distinct else target
            errors.extend({**identity, "result": index, "error": error} for error in validate_result(case, result, query, searched))
            warnings.extend({**identity, "result": index, "warning": warning} for warning in validation_warnings(case, query, searched))
            if result["optimal"]:
                matches = result["matches"]
                objective = _score(case, dict(matches[0]["atoms"]), matches[0]["bonds"]) if matches else (0, 0, 0)
                indices = None if distinct else index_maps.get(graph_key)
                normalised = normalised_matches(matches, *indices) if indices else None
                # A batch row that names no targets repeats the case target.
                paired[row["case_id"], row["phase"], row["sample"], index, row.get("batch_targets", "copies")].append(
                    (row["language"], objective, len(matches), result["complete"], normalised))
    for key, values in paired.items():
        completed = [value for value in values if value[3]]
        if (len({value[1] for value in values}) > 1
                or len({value[2] if int(cases[key[0]]["limit"]) else bool(value[2]) for value in completed}) > 1):
            errors.append({"case_id": key[0], "phase": key[1], "sample": key[2], "result": key[3],
                           "error": "Proven cross-language objectives or completed enumeration counts differ",
                           "outcomes": [value[:4] for value in values]})
        normalised = [value[4] for value in completed if value[4] is not None]
        stereo = cases[key[0]]["profile"] == "stereo" or (cases[key[0]]["operation"] == "smarts"
            and re.search(r"[@/\\]", cases[key[0]]["query"]))
        if not stereo and int(cases[key[0]]["limit"]) and normalised and any(value != normalised[0] for value in normalised[1:]):
            errors.append({"case_id": key[0], "phase": key[1], "sample": key[2], "result": key[3],
                           "error": "Complete atom/bond mapping sets differ after parser-index normalisation"})
    return {"rows": len(rows), "errors": errors, "graph_differences": graph_differences,
            "warnings": warnings, "valid": not errors and not graph_differences}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cases", required=True)
    parser.add_argument("--results", nargs="+", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    report = compare_jsonl(args.cases, args.results)
    with open(args.output, "w", encoding="utf-8") as output:
        json.dump(report, output, indent=2)
        output.write("\n")
    raise SystemExit(0 if report["valid"] else 1)


if __name__ == "__main__":
    main()
