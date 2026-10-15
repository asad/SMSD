# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
# Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
"""Parser resource limits, integer bounds and Boolean bond expressions."""
import subprocess
import sys

import pytest
import smsd


@pytest.mark.parametrize("body", [
    "g = smsd.parse_smiles('C' + '(' * 100000 + 'O' + ')' * 100000 + 'N'); assert g.n == 3",
    "assert smsd.Substructure.smarts('[' + '!' * 100000 + 'C]').match('C').found",
    "assert smsd.Substructure.smarts('C' + '!' * 100000 + '-C').match('CC').found",
])
def test_deep_syntax_survives(body):
    result = subprocess.run([sys.executable, "-c", "import smsd; " + body], capture_output=True, text=True, timeout=30)
    assert result.returncode == 0, result.stderr


@pytest.mark.parametrize("text", [
    "[99999999999C]", "[CH99999999999]", "[C+99999999999]",
    "[C-99999999999]", "[C:99999999999]", "[R99999999999]",
    "[CH2147483647]C", "[H][CH2147483647]",
])
def test_smiles_rejects_unrepresentable_numbers(text):
    with pytest.raises(ValueError):
        smsd.parse_smiles(text)


@pytest.mark.parametrize("text", [
    "[#99999999999]", "[H99999999999]", "[X99999999999]",
    "[+99999999999]", "[-99999999999]", "[99999999999C]",
    "[:99999999999]", "[D{1-99999999999}]", "C%99999999999CC%99999999999",
])
def test_smarts_rejects_unrepresentable_numbers(text):
    with pytest.raises(ValueError):
        smsd.Substructure.smarts(text).match(smsd.parse_smiles('CC')).found


def test_representable_numeric_boundaries():
    graph = smsd.parse_smiles("[CH2147483647]")
    assert graph.hydrogen_count == [2147483647]
    assert smsd.Substructure.smarts('[H2147483647;X2147483647;v2147483647]').match(graph).found
    graph = smsd.parse_smiles("[2147483647C+2147483647:2147483647]")
    assert smsd.Substructure.smarts('[2147483647C+2147483647:2147483647]').match(graph).found
    graph = smsd.parse_smiles("[O+2147483647]")
    assert smsd.Substructure.smarts('[O+2147483647;H0;X0;v0]').match(graph).found


@pytest.mark.parametrize("query,target,expected", [
    ("C=,#@C", "C=C", True), ("C=,#@C", "C#C", False),
    ("C=,#@C", "C1#CCCCC1", True), ("C=,#@C", "CC", False),
    ("C-&!@C", "CC", True), ("C-&!@C", "C1CC1", False),
    ("C-,=;@C", "C=C", False), ("C-,=;@C", "C1CC1", True),
    ("C-,=;@C", "C1=CCCCC1", True),
])
def test_bond_boolean_precedence(query, target, expected):
    assert smsd.Substructure.smarts(query).match(smsd.parse_smiles(target)).found == expected


def test_large_chain_keeps_double_bond_stereo():
    graph = smsd.parse_smiles("C/C=C/" + "C" * 4093)
    assert graph.n == 4096
    assert smsd.Substructure.smarts('C/C=C/C').match(graph).found
    assert not smsd.Substructure.smarts('C/C=C\\C').match(graph).found
