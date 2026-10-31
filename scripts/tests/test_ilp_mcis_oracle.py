# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
# Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
# See the NOTICE file for attribution, trademark, and algorithm IP terms.
"""Small independent MILP controls; SciPy is an optional benchmark dependency."""
import importlib.util
import json
from pathlib import Path
import unittest


SOURCE = Path(__file__).resolve().parents[2] / "benchmarks/ilp_mcis_oracle.py"
SPEC = importlib.util.spec_from_file_location("ilp_mcis_oracle_tested", SOURCE)
ORACLE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(ORACLE)
HAS_SCIPY = importlib.util.find_spec("scipy") is not None


def graph(elements, edges=(), aromatic=False):
    return {"atoms": [{"element": element, "aromatic": aromatic} for element in elements],
            "bonds": [[first, second, 4 if aromatic else order, aromatic, aromatic]
                      for first, second, order in edges]}


@unittest.skipUnless(HAS_SCIPY, "SciPy is required only for MILP benchmark controls")
class MCISOracleTest(unittest.TestCase):
    def test_five_small_graph_pairs(self):
        cases = [
            (graph([6] * 3, [(0, 1, 1), (1, 2, 1)]),
             graph([6] * 3, [(0, 1, 1), (1, 2, 1), (0, 2, 1)]), (2, 1), (2, 1)),
            (graph([6, 6, 8]), graph([6, 8, 8]), (1, 0), (2, 0)),
            (graph([6, 6, 8, 7], [(0, 1, 1), (1, 2, 2), (1, 3, 1)]),
             graph([6, 6, 8, 8], [(0, 1, 1), (1, 2, 2), (1, 3, 1)]), (3, 2), (3, 2)),
            (graph([6] * 6, [(atom, (atom + 1) % 6, 1) for atom in range(6)], True),
             graph([6] * 6, [(atom, (atom + 1) % 6, 1) for atom in range(6)]), (6, 6), (6, 6)),
            (graph([6] * 4, [(0, 1, 1), (0, 2, 1), (0, 3, 1)]),
             graph([6] * 4, [(0, 1, 1), (1, 2, 1), (2, 3, 1)]), (3, 2), (3, 2)),
        ]
        for index, (query, target, connected_score, fragment_score) in enumerate(cases):
            for connected, expected in ((True, connected_score), (False, fragment_score)):
                with self.subTest(case=index, connected=connected):
                    result = ORACLE.solve_mcis(query, target, connected=connected)
                    self.assertEqual("optimal", result["status"], result)
                    self.assertTrue(result["optimal"], result)
                    self.assertTrue(result["witness_checked"])
                    self.assertFalse(result["formal_certificate"])
                    self.assertEqual(expected, (result["atom_count"], result["bond_count"]))

    def test_model_limits_remain_explicit_stopped_outcomes(self):
        source = graph([6] * 3, [(0, 1, 1), (1, 2, 1)])
        for limits, reason in (({"max_variables": 1}, "variable_limit"),
                               ({"max_rows": 1}, "row_limit"),
                               ({"max_nonzeros": 1}, "nonzero_limit"),
                               ({"max_pair_checks": 1}, "pair_check_limit"),
                               ({"time_limit": 1e-12}, "time_limit")):
            with self.subTest(limits=limits):
                result = ORACLE.solve_mcis(source, source, **limits)
                self.assertEqual("model_limit", result["status"])
                self.assertEqual(reason, result["message"])
                self.assertFalse(result["optimal"])


class MCISOracleInputTest(unittest.TestCase):
    def test_no_compatible_atom_pair_has_zero_optimum(self):
        result = ORACLE.solve_mcis(graph([7]), graph([8]))
        self.assertTrue(result["optimal"])
        self.assertEqual((0, 0), (result["atom_count"], result["bond_count"]))

    def test_rejected_inputs_do_not_become_optima(self):
        source = graph([6, 6], [(0, 1, 1)])
        duplicate = graph([6, 6], [(0, 1, 1), (1, 0, 1)])
        for target, limits in ((duplicate, {}), (source, {"max_atoms": 1}),
                               (source, {"time_limit": 0}), (source, {"connected": "yes"})):
            with self.subTest(target=target, limits=limits):
                result = ORACLE.solve_mcis(source, target, **limits)
                self.assertEqual("input_rejected", result["status"])
                self.assertFalse(result["optimal"])

    def test_nonfinite_time_limits_keep_a_serialisable_failure(self):
        source = graph([6])
        for seconds in (float("nan"), float("inf"), True):
            result = ORACLE.solve_mcis(source, source, time_limit=seconds)
            self.assertEqual("input_rejected", result["status"])
            json.dumps(result, allow_nan=False)


if __name__ == "__main__":
    unittest.main()
