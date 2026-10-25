# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
# Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
# See the NOTICE file for attribution, trademark, and algorithm IP terms.
"""Independent graph-oracle controls for the cross-language measurement checker."""
import importlib.util
import csv
import json
from pathlib import Path
import tempfile
import unittest


SOURCE = Path(__file__).resolve().parents[2] / "benchmarks/search_api_validation.py"
SPEC = importlib.util.spec_from_file_location("search_api_validation_tested", SOURCE)
VALIDATION = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(VALIDATION)


def graph(count, edges=()):
    return {"atoms": [{"element": 6, "charge": 0, "isotope": 0, "hydrogens": 0,
                       "aromatic": False, "ring": False} for _ in range(count)],
            "bonds": [[first, second, order, False, False] for first, second, order in edges]}


def case(**values):
    return {"operation": "mcs", "objective": "MCIS", "connectivity": "Connected",
            "profile": "default", "induced": "0", "limit": "0", "weights": "", **values}


class SearchValidationTest(unittest.TestCase):
    def test_parser_permutations_preserve_query_weight_assignments(self):
        source = graph(3, [(0, 1, 1), (1, 2, 2)])
        source["atoms"][0]["element"] = 7
        source["atoms"][2]["element"] = 8
        permuted = graph(3, [(2, 0, 1), (0, 1, 2)])
        permuted["atoms"][2]["element"] = 7
        permuted["atoms"][1]["element"] = 8
        stopped = {"status": "TimedOut", "optimal": False, "complete": False,
                   "score": 0, "nodes_visited": 0, "matches": []}
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            results = []
            for language, query in (("first", source), ("second", permuted)):
                path = directory / (language + ".jsonl")
                path.write_text(json.dumps({"case_id": "weighted", "language": language,
                    "phase": "prepared", "sample": 0, "elapsed_ns": 1,
                    "results": [stopped], "query_graph": query, "target_graph": source}) + "\n")
                results.append(path)
            for weights, comparable in (("1,1,1", True), ("1,2,3", False)):
                row = case(case_id="weighted", weights=weights, timeout_ms="500")
                cases = directory / "cases.tsv"
                with cases.open("w", newline="") as handle:
                    writer = csv.DictWriter(handle, fieldnames=list(row), delimiter="\t")
                    writer.writeheader()
                    writer.writerow(row)
                report = VALIDATION.compare_jsonl(cases, results)
                self.assertEqual(report["errors"], [])
                self.assertEqual(report["valid"], comparable)
                self.assertEqual(bool(report["graph_differences"]), not comparable)

    def test_distinct_batch_targets_are_checked_against_their_own_graphs(self):
        query, short, long = graph(2, [(0, 1, 1)]), graph(2, [(0, 1, 1)]), graph(3, [(0, 1, 1), (1, 2, 1)])
        found = {"status": "Completed", "optimal": True, "complete": True, "score": 2, "nodes_visited": 1,
                 "matches": [{"atoms": [[0, 0], [1, 1]], "bonds": [[0, 1, 0, 1]], "fragments": [[0, 1]]}]}
        outside = {**found, "matches": [{"atoms": [[0, 1], [1, 2]], "bonds": [[0, 1, 1, 2]], "fragments": [[0, 1]]}]}
        row = case(case_id="batch", batch_size="2", timeout_ms="500")
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            cases = directory / "cases.tsv"
            with cases.open("w", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=list(row), delimiter="\t")
                writer.writeheader()
                writer.writerow(row)
            for second, valid in ((outside, False), (found, True)):
                path = directory / "batch.jsonl"
                path.write_text(json.dumps({"case_id": "batch", "language": "first", "phase": "batch", "sample": 0,
                    "threads": 2, "batch_targets": "distinct", "elapsed_ns": 1, "results": [found, second],
                    "query_graph": query, "target_graph": short, "batch_target_graphs": [long, short]}) + "\n")
                report = VALIDATION.compare_jsonl(cases, [path])
                self.assertEqual(report["valid"], valid)
                if not valid:
                    self.assertEqual({error["result"] for error in report["errors"]}, {1})

    def test_batch_rows_without_their_targets_pair_with_copies(self):
        # A batch row that names no targets repeats the case target, so it pairs with rows
        # marked as copies, and differing proved objectives are reported.
        query = graph(2, [(0, 1, 1)])
        found = {"status": "Completed", "optimal": True, "complete": True, "score": 2, "nodes_visited": 1,
                 "matches": [{"atoms": [[0, 0], [1, 1]], "bonds": [[0, 1, 0, 1]], "fragments": [[0, 1]]}]}
        single = {**found, "score": 1, "matches": [{"atoms": [[0, 0]], "bonds": [], "fragments": [[0]]}]}
        row = case(case_id="batch", batch_size="1", timeout_ms="500")
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            cases = directory / "cases.tsv"
            with cases.open("w", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=list(row), delimiter="\t")
                writer.writeheader()
                writer.writerow(row)
            paths = []
            for language, result, extra in (("old", found, {}), ("new", single, {"batch_targets": "copies"})):
                path = directory / (language + ".jsonl")
                path.write_text(json.dumps({"case_id": "batch", "language": language, "phase": "batch", "sample": 0,
                    "threads": 1, **extra, "elapsed_ns": 1, "results": [result],
                    "query_graph": query, "target_graph": query}) + "\n")
                paths.append(path)
            report = VALIDATION.compare_jsonl(cases, paths)
            self.assertIn("Proven cross-language objectives or completed enumeration counts differ",
                          [error["error"] for error in report["errors"]])

    def test_parser_permutation_preserves_atom_and_selected_bond_correspondence(self):
        source = graph(3, [(0, 1, 1), (1, 2, 2)])
        source["atoms"][0]["element"] = 7
        reference = graph(3, [(2, 0, 1), (0, 1, 2)])
        reference["atoms"][2]["element"] = 7
        indices = VALIDATION.graph_index_map(source, reference)
        self.assertEqual(indices, {0: 2, 1: 0, 2: 1})
        match = {"atoms": [[0, 0], [1, 1], [2, 2]], "bonds": [[0, 1, 0, 1]]}
        transformed = VALIDATION.normalised_matches([match], indices, indices)
        self.assertEqual(transformed, {(((0, 0), (1, 1), (2, 2)), ((0, 2),))})

    def test_triangle_path_objectives_and_all_mapping_counts(self):
        triangle = graph(3, [(0, 1, 1), (1, 2, 1), (0, 2, 1)])
        path = graph(3, [(0, 1, 1), (1, 2, 1)])
        optimum, matches = VALIDATION.exhaustive_mcs(case(), triangle, path)
        self.assertEqual(optimum, (2, 2, 1))
        self.assertEqual(len(matches), 12)
        optimum, matches = VALIDATION.exhaustive_mcs(case(objective="MCES"), triangle, path)
        self.assertEqual(optimum, (2, 3, 0))
        self.assertEqual(len(matches), 6)

    def test_disconnected_no_edge_graphs_have_distinct_objectives(self):
        query, target = graph(2), graph(3, [(0, 1, 1), (1, 2, 1)])
        optimum, matches = VALIDATION.exhaustive_mcs(case(connectivity="Fragments"), query, target)
        self.assertEqual(optimum, (2, 2, 0))
        self.assertEqual(len(matches), 2)
        self.assertEqual(VALIDATION.exhaustive_mcs(case(objective="MCES", connectivity="Fragments"), query, target),
                         ((0, 0, 0), set()))

    def test_signed_and_zero_weights_use_atom_count_tie_breaks(self):
        query = graph(2, [(0, 1, 1)])
        optimum, matches = VALIDATION.exhaustive_mcs(case(weights="-0.25,0"), query, query)
        self.assertEqual(optimum, (0, 1, 0))
        self.assertEqual(len(matches), 2)
        optimum, matches = VALIDATION.exhaustive_mcs(case(weights="0,0"), query, query)
        self.assertEqual(optimum, (0, 2, 1))
        self.assertEqual(len(matches), 2)

    def test_observed_loose_and_any_orders_are_permissive_after_guards(self):
        triple = VALIDATION.Graph(graph(2, [(0, 1, 3)]))
        single = VALIDATION.Graph(graph(2, [(0, 1, 1)]))
        self.assertTrue(VALIDATION._bond_matches(triple, 0, 1, single, 0, 1, "tautomer"))
        self.assertTrue(VALIDATION._bond_matches(triple, 0, 1, single, 0, 1, "any"))
        self.assertFalse(VALIDATION._bond_matches(triple, 0, 1, single, 0, 1, "strict"))
        self.assertFalse(VALIDATION._bond_matches(triple, 0, 1, single, 0, 1, "default"))

    def test_checker_rejects_induced_nonedge_and_false_optimum(self):
        triangle = graph(3, [(0, 1, 1), (1, 2, 1), (0, 2, 1)])
        path = graph(3, [(0, 1, 1), (1, 2, 1)])
        result = {"status": "Completed", "optimal": True, "complete": True, "score": 3,
                  "nodes_visited": 1, "matches": [{"atoms": [[0, 0], [1, 1], [2, 2]],
                  "bonds": [[0, 1, 0, 1], [1, 2, 1, 2]], "fragments": [[0, 1, 2]]}]}
        errors = VALIDATION.validate_result(case(), result, triangle, path)
        self.assertIn("Mapping violates induced adjacency", errors)
        self.assertIn("Claimed optimum differs from the independent exhaustive oracle", errors)

    def test_stopped_suboptimal_witness_is_not_claimed_as_wrong_answer(self):
        query = graph(2, [(0, 1, 1)])
        result = {"status": "TimedOut", "optimal": False, "complete": False, "score": 1,
                  "nodes_visited": 1, "matches": [{"atoms": [[0, 0]], "bonds": [], "fragments": [[0]]}]}
        self.assertEqual(VALIDATION.validate_result(case(), result, query, query), [])

    def test_boolean_hydrogen_and_recursive_smarts_are_independent_predicates(self):
        target = graph(5, [(0, 1, 1), (1, 2, 2), (1, 3, 1), (3, 4, 1)])
        target["atoms"][2]["element"] = 8
        target["atoms"][3].update(element=7, hydrogens=1)
        target["atoms"][0]["hydrogens"] = 3
        methyl = VALIDATION.Pattern("[C;H3]")
        self.assertEqual(VALIDATION._embeddings(methyl, VALIDATION.Graph(target), False), {((0, 0),)})
        amide = VALIDATION.Pattern("[$([N;H1,H2]-C=O)]")
        self.assertEqual(VALIDATION._embeddings(amide, VALIDATION.Graph(target), False), {((0, 3),)})

    def test_bare_hydrogen_keeps_the_whole_boolean_atom_context(self):
        target = graph(5, [(0, 1, 1), (1, 2, 1)])
        target["atoms"][0]["hydrogens"] = 3
        target["atoms"][1].update(element=7, hydrogens=1)
        target["atoms"][2]["hydrogens"] = 3
        target["atoms"][3]["element"] = 1
        target["atoms"][4]["element"] = 17
        for expression, expected in (
                ("[N;H;X3]", {1}), ("[H,#1]", {1, 3}),
                ("[H,Cl]", {1, 4}), ("[C,H]", {0, 1, 2}),
                ("[N;!H0]", {1}), ("[N;H0;X3]", set()),
                ("[$([N;H;X3]([CX4])[CX4])]", {1})):
            with self.subTest(expression=expression):
                mappings = VALIDATION._embeddings(VALIDATION.Pattern(expression),
                                                  VALIDATION.Graph(target), False)
                self.assertEqual(mappings, {((0, atom),) for atom in expected})

    def test_atomic_hydrogen_keeps_isotope_charge_and_class_context(self):
        target = graph(6)
        for atom in target["atoms"][:5]:
            atom["element"] = 1
        target["atoms"][1]["isotope"] = 2
        target["atoms"][2]["charge"] = 1
        target["atoms"][3]["charge"] = -1
        target["atoms"][4].update(isotope=2, charge=1)
        target["atoms"][5].update(isotope=2, charge=1, hydrogens=1)
        for expression, expected in (
                ("[H]", {0, 1, 2, 3, 4}), ("[H+]", {2, 4}),
                ("[H-]", {3}), ("[2H]", {1, 4}),
                ("[H:1]", {0, 1, 2, 3, 4}), ("[2H+]", {4})):
            with self.subTest(expression=expression):
                mappings = VALIDATION._embeddings(VALIDATION.Pattern(expression),
                                                  VALIDATION.Graph(target), False)
                self.assertEqual(mappings, {((0, atom),) for atom in expected})


if __name__ == "__main__":
    unittest.main()
