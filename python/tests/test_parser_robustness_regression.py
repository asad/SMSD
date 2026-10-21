# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
# Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
# See the NOTICE file for attribution, trademark, and algorithm IP terms.
"""Parser resource limits, integer bounds and Boolean bond expressions."""
import subprocess
import sys

import pytest
import smsd


@pytest.mark.parametrize("body", [
    "g = smsd.parse_smiles('C' + '(' * 100000 + 'O' + ')' * 100000 + 'N'); assert g.n == 3",
    # The longest negation chains within the 100,000-character SMARTS limit.
    "assert smsd.Substructure.smarts('[' + '!' * 99996 + 'C]').match('C').found",
    "assert smsd.Substructure.smarts('C' + '!' * 99996 + '-C').match('CC').found",
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
    "[CH10]", "[CH00]", "[CH2147483647]", "[C+16]", "[C-16]", "[C+128]", "[C+016]",
    "[C++++++++++++++++]", "[1000C]", "[65536C]", "[2147483647C+2147483647:2147483647]", "[O+2147483647]",
])
def test_smiles_rejects_bracket_fields_outside_opensmiles(text):
    """A hydrogen count has one digit, a charge two within -15 to +15, and an isotope three."""
    with pytest.raises(ValueError):
        smsd.parse_smiles(text)


def test_smiles_bracket_fields_at_their_limits_round_trip():
    graph = smsd.parse_smiles("[999CH9+15:2147483647]")
    assert list(graph.mass_number) == [999] and list(graph.hydrogen_count) == [9]
    assert list(graph.formal_charge) == [15] and list(graph.atom_class) == [2147483647]
    assert smsd.to_smiles(graph) == "[999CH9+15:2147483647]"
    assert list(smsd.parse_smiles("[C-15]").formal_charge) == [-15]
    assert list(smsd.parse_smiles("[C+++++++++++++++]").formal_charge) == [15]


@pytest.mark.parametrize("text", [
    "[#99999999999]", "[H99999999999]", "[X99999999999]",
    "[+99999999999]", "[-99999999999]", "[99999999999C]",
    "[:99999999999]", "[D{1-99999999999}]", "C%99999999999CC%99999999999",
])
def test_smarts_rejects_unrepresentable_numbers(text):
    with pytest.raises(ValueError):
        smsd.Substructure.smarts(text).match(smsd.parse_smiles('CC')).found


def test_representable_numeric_boundaries():
    """SMARTS keeps its own, wider numeric fields, which a graph built directly can hold."""
    graph = smsd.MolGraphBuilder().atom_count(1).atomic_numbers([6]).hydrogen_counts([2147483647]).neighbors([[]]).build()
    assert list(graph.hydrogen_count) == [2147483647]
    assert smsd.Substructure.smarts('[H2147483647;X2147483647;v2147483647]').match(graph).found
    graph = (smsd.MolGraphBuilder().atom_count(1).atomic_numbers([8]).formal_charges([2147483647])
             .mass_numbers([2147483647]).hydrogen_counts([0]).neighbors([[]]).build())
    assert smsd.Substructure.smarts('[2147483647O+2147483647:2147483647]').match(graph).found
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


def test_smarts_input_limits():
    """SMARTS length, atoms (recursive SMARTS included) and expression nodes are bounded inclusively."""
    assert smsd.Substructure.smarts('[' + '!' * 99997 + 'C]').match('N').found
    with pytest.raises(ValueError, match="maximum length"):
        smsd.Substructure.smarts('[' + '!' * 99998 + 'C]')
    assert smsd.Substructure.smarts('C' * 10000).match('C' * 10000).found
    with pytest.raises(ValueError, match="maximum of 10000 atoms"):
        smsd.Substructure.smarts('C' * 10001)
    with pytest.raises(ValueError, match="maximum of 10000 atoms"):
        smsd.Substructure.smarts('[$(' + 'C' * 10000 + ')]')
    ring_atom = '[C' + 'R' * 4999 + ']'  # 9,999 nodes: primitives and the ANDs between them
    assert smsd.Substructure.smarts(ring_atom).match('C1CC1').found
    with pytest.raises(ValueError, match="atom or bond expression"):
        smsd.Substructure.smarts('[C' + 'R' * 5000 + ']')
    assert smsd.Substructure.smarts(ring_atom * 10)  # 99,999 nodes with the nine bonds
    with pytest.raises(ValueError, match="expression nodes"):
        smsd.Substructure.smarts(ring_atom * 10 + '.[CR]')


def test_large_macrocycle_parses_on_a_small_stack():
    """Ring perception of a macrocycle at the atom limit needs no deep call stack."""
    body = (
        "import threading, smsd\n"
        "threading.stack_size(512 * 1024)\n"
        "out = []\n"
        "t = threading.Thread(target=lambda: out.append(smsd.parse_smiles('C1' + 'C' * 9997 + 'C1').n))\n"
        "t.start(); t.join()\n"
        "assert out == [9999], out\n"
    )
    result = subprocess.run([sys.executable, "-c", body], capture_output=True, text=True, timeout=60)
    assert result.returncode == 0, result.stderr
