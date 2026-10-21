# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
# Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
# See the NOTICE file for attribution, trademark, and algorithm IP terms.
"""Invalid native inputs must raise exceptions without terminating the process."""
import subprocess
import sys

import pytest


def _expect_error(setup, expression, error):
    script = "import smsd\nfrom smsd import _smsd as native\n" + setup + "\n"
    script += f"try:\n    {expression}\nexcept {error}:\n    pass\nelse:\n    raise AssertionError('Expected {error}')\n"
    result = subprocess.run([sys.executable, "-c", script], capture_output=True, text=True, timeout=20)
    assert result.returncode == 0, result.stderr


@pytest.mark.parametrize("operation", ["pair", "batch", "corpus"])
@pytest.mark.parametrize("weights", ["[1.0]", "[1.0, float('nan'), 1.0]", "[1.0, float('inf'), 1.0]"])
def test_parallel_mcs_rejects_invalid_weights(operation, weights):
    setup = (
        "query = smsd.parse_smiles('CCO')\n"
        "targets = [smsd.parse_smiles('CCO'), smsd.parse_smiles('CCCO')]\n"
        "options = smsd.MCS.Options()\n"
        f"options.atom_weights = {weights}\n"
    )
    if operation == "corpus":
        setup += "corpus = smsd.TargetCorpus()\ncorpus.add_targets(targets)\n"
        expression = "smsd.MCS.find(query).with_options(options).match_many(corpus, num_threads=2)"
    elif operation == "batch":
        expression = "smsd.MCS.find(query).with_options(options).match_many(targets, num_threads=2)"
    else:
        expression = "smsd.MCS.find(query).with_options(options).match(targets[0])"
    _expect_error(setup, expression, "ValueError")


@pytest.mark.parametrize("index", [-1, 3, 100000000])
@pytest.mark.parametrize("operation", ["mcs", "ring", "extract"])
def test_atom_indices_are_checked(operation, index):
    expression = {
        "mcs": f"smsd.Match({{{index}: 0}}, (), ()).to_smiles(query)",
        "ring": f"query.ring_system_of({index})",
        "extract": f"native.extract_subgraph(query, [{index}])",
    }[operation]
    _expect_error("query = smsd.parse_smiles('CCO')", expression, "IndexError")


def test_bond_normalisation_rejects_short_coordinates():
    _expect_error("query = smsd.parse_smiles('CCO')",
                  "native.normalise_bond_length(query, [(0, 0), (1, 0)])", "ValueError")


@pytest.mark.parametrize("field,value", [
    ("min_fragment_atoms", 0), ("min_fragment_atoms", -1),
    ("max_fragments", 0), ("max_fragments", -1),
    ("atom_weights", [1e300, 1, 1]),
])
def test_search_options_are_checked(field, value):
    setup = f"options = smsd.MCS.Options()\noptions.{field} = {value!r}"
    _expect_error(setup, "smsd.MCS.find('CCO').with_options(options).match('CCO')", "ValueError")


@pytest.mark.parametrize("indices", [[-1], [3], [100000000]])
def test_excluded_target_atoms_are_checked(indices):
    setup = f"options = smsd.MCS.Options()\noptions.excluded_target_atoms = {indices!r}"
    _expect_error(setup, "smsd.MCS.find('CCO').with_options(options).match('CCO')", "IndexError")


@pytest.mark.parametrize("indices", [[-1, 1], [0, 3], [0, 100000000]])
def test_selected_bond_atoms_are_checked(indices):
    first, second = indices
    expression = f"smsd.Match({{0: 0, 1: 1}}, (({first}, {second}, 0, 1),), ()).to_smiles(query)"
    _expect_error("query = smsd.parse_smiles('CCO')", expression, "IndexError")
