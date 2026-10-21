# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
# Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
# See the NOTICE file for attribution, trademark, and algorithm IP terms.
"""Fragment and atom writers preserve hydrogen counts and valid aromaticity."""
import pytest
import smsd
from .search_helpers import selected_atoms


@pytest.mark.parametrize("text", [
    "[SiH4]", "C[SiH2]C", "[n-]1cccc1", "[2H]", "[CH3]",
    "[nH]1cccc1", "c1ncc[nH]1", "[NH4+]", "[nH+]1ccccc1",
])
def test_smiles_preserves_stored_hydrogens(text):
    original = smsd.parse_smiles(text)
    output = smsd.to_smiles(original)
    restored = smsd.parse_smiles(output)
    assert sum(restored.hydrogen_count) == sum(original.hydrogen_count), output


@pytest.mark.parametrize("text", ["[SiH4]", "[n-]1cccc1", "[nH]1cccc1", "c1ncc[nH]1", "[2H]", "[CH3]"])
def test_hydrogen_qualified_smarts_matches_its_source(text):
    original = smsd.parse_smiles(text)
    output = smsd.to_smarts(original, include_h_count=True)
    assert smsd.Substructure.smarts(output).match(original).found, output


@pytest.mark.parametrize("text", ["c1ncc[nH]1", "c1ccc2[nH]ccc2c1", "O=c1cc[nH]c(=O)[nH]1", "[CH3]"])
def test_complete_mcs_preserves_hydrogens(text):
    original = smsd.parse_smiles(text)
    output = selected_atoms(original, range(original.n))
    restored = smsd.parse_smiles(output)
    assert sum(restored.hydrogen_count) == sum(original.hydrogen_count), output
    assert restored.kekulize(), output


@pytest.mark.parametrize("text,atoms", [
    ("c1ccccc1", list(range(5))),
    ("CCCO", [1, 2, 3]),
])
def test_partial_mcs_caps_cut_valences(text, atoms):
    output = selected_atoms(smsd.parse_smiles(text), atoms)
    fragment = smsd.parse_smiles(output)
    assert not any(fragment.aromatic), output
    for atom, atomic_number in enumerate(fragment.atomic_num):
        if atomic_number == 6:
            assert fragment.hydrogen_count[atom] + sum(
                order for first, second, order in fragment.bonds() if atom in (first, second)) == 4, output


def test_retained_phenyl_ring_keeps_aromaticity_and_caps_the_cut():
    output = selected_atoms(smsd.parse_smiles("Cc1ccccc1"), range(1, 7))
    fragment = smsd.parse_smiles(output)
    assert sum(fragment.aromatic) == 6
    assert sum(fragment.hydrogen_count) == 6
