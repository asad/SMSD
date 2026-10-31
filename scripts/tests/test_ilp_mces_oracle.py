# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
# Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
# See the NOTICE file for attribution, trademark, and algorithm IP terms.
"""Small selected-edge controls for the optional independent MCES relaxation."""
import importlib.util
import json
from pathlib import Path
import unittest

SOURCE = Path(__file__).resolve().parents[2] / "benchmarks/ilp_mces_oracle.py"
SPEC = importlib.util.spec_from_file_location("ilp_mces_oracle_tested", SOURCE)
ORACLE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(ORACLE)


def graph(elements, edges=()):
    return {"atoms": [{"element": element, "aromatic": False} for element in elements],
            "bonds": [[first, second, order, False, False] for first, second, order in edges]}


@unittest.skipUnless(importlib.util.find_spec("scipy"), "SciPy is required only for MILP controls")
class MCESRelaxationTest(unittest.TestCase):
    def test_selected_edges_allow_a_triangle_to_match_a_chain(self):
        query = graph([6] * 3, [(0, 1, 1), (1, 2, 1), (0, 2, 1)])
        target = graph([6] * 3, [(0, 1, 2), (1, 2, 3)])
        result = ORACLE.solve_mces_any_relaxation(query, target)
        self.assertEqual("optimal", result["status"], result)
        self.assertEqual((2, 3), (result["bond_count"], result["atom_count"]))
        self.assertTrue(result["original_contract_feasible"])
        self.assertTrue(result["original_contract_optimal"])
        self.assertEqual(11, result["upper_bound"]["scalar"])
        self.assertFalse(result["formal_certificate"])

    def test_two_atom_components_do_not_prove_the_minimum_three_contract(self):
        query = graph([6] * 4, [(0, 1, 2), (2, 3, 3)])
        target = graph([6] * 4, [(0, 1, 1), (1, 2, 1), (2, 3, 1)])
        result = ORACLE.solve_mces_any_relaxation(query, target)
        self.assertEqual("optimal", result["status"], result)
        self.assertEqual((2, 4), (result["bond_count"], result["atom_count"]))
        self.assertEqual([2, 2], sorted(map(len, result["fragments"])))
        self.assertTrue(result["optimal"])
        self.assertFalse(result["original_contract_feasible"])
        self.assertFalse(result["original_contract_optimal"])
        self.assertEqual({"bond_count": 0, "atom_count": 0}, result["original_contract_lower_bound"])

    def test_model_caps_remain_stopped_outcomes(self):
        source = graph([6] * 3, [(0, 1, 1), (1, 2, 1)])
        for limits, reason in (({"max_variables": 1}, "variable_limit"),
                               ({"max_rows": 1}, "row_limit"),
                               ({"max_nonzeros": 1}, "nonzero_limit"),
                               ({"max_edge_checks": 1}, "edge_check_limit"),
                               ({"time_limit": 1e-12}, "time_limit")):
            with self.subTest(limits=limits):
                result = ORACLE.solve_mces_any_relaxation(source, source, **limits)
                self.assertEqual(("model_limit", reason), (result["status"], result["message"]))
                self.assertFalse(result["optimal"])
                self.assertFalse(result["original_contract_optimal"])
                self.assertEqual("graph counts", result["upper_bound"]["source"])


class MCESWitnessTest(unittest.TestCase):
    def test_independently_rejects_bad_native_mappings(self):
        source = graph([6] * 3, [(0, 1, 1), (1, 2, 1)])
        valid = ORACLE.check_any_mces_witness(source, source, [(0, 0), (1, 1), (2, 2)],
                                              [(0, 1, 0, 1), (1, 2, 1, 2)])
        self.assertTrue(valid["original_contract_feasible"])
        for atoms, edges in (([(0, 0), (1, 0)], [(0, 1, 0, 0)]),
                              ([(0, 0), (1, 1)], [(0, 2, 0, 1)]),
                              ([(0, 0), (1, 1), (2, 2)], [(0, 1, 0, 1)])):
            with self.subTest(atoms=atoms, edges=edges):
                with self.assertRaises(ValueError):
                    ORACLE.check_any_mces_witness(source, source, atoms, edges)

    def test_no_compatible_edge_proves_empty_without_a_solver(self):
        result = ORACLE.solve_mces_any_relaxation(graph([6, 8], [(0, 1, 1)]),
                                                 graph([6, 7], [(0, 1, 1)]))
        self.assertTrue(result["optimal"])
        self.assertTrue(result["original_contract_optimal"])
        self.assertEqual((0, 0), (result["bond_count"], result["atom_count"]))

    def test_rejected_limits_are_serialisable(self):
        source = graph([6, 6], [(0, 1, 1)])
        for limits in ({"time_limit": float("nan")}, {"time_limit": float("inf")},
                       {"time_limit": True}, {"min_fragment_atoms": 1}, {"max_fragments": 0}):
            result = ORACLE.solve_mces_any_relaxation(source, source, **limits)
            self.assertEqual("input_rejected", result["status"])
            self.assertFalse(result["optimal"])
            json.dumps(result, allow_nan=False)


if __name__ == "__main__":
    unittest.main()
