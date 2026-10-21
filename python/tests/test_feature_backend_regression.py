# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
# Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
# See the NOTICE file for attribution, trademark, and algorithm IP terms.
"""Stored hydrogen and pharmacophore feature regressions."""
import pytest
import smsd


@pytest.mark.parametrize("smiles,atom", [("CC(=O)C", 2), ("c1ccncc1", 3), ("CC#N", 2)])
def test_acceptors_without_hydrogen_are_not_donors(smiles, atom):
    features = smsd.classify_pharmacophore(smsd.parse_smiles(smiles), atom)
    assert features & 2
    assert not features & 1


@pytest.mark.parametrize("smiles,atom", [("[NH4+]", 0), ("C[NH3+]", 1), ("CO", 1), ("[H]OC", 1), ("c1cc[nH]c1", 3)])
def test_attached_hydrogen_donors(smiles, atom):
    assert smsd.classify_pharmacophore(smsd.parse_smiles(smiles), atom) & 1


@pytest.mark.parametrize("smiles", ["CC#N", "C=N"])
def test_multiple_bond_nitrogen_is_not_a_saturated_basic_amine(smiles):
    molecule = smsd.parse_smiles(smiles)
    assert not smsd.classify_pharmacophore(molecule, molecule.n - 1) & 4


@pytest.mark.parametrize("first,second", [("[NH4+]", "[NH2+]"), ("C[NH3+]", "C[N+]"), ("[CH2]", "C")])
def test_ecfp_retains_stored_hydrogen_counts(first, second):
    first, second = smsd.parse_smiles(first), smsd.parse_smiles(second)
    assert smsd.circular_fingerprint(first) != smsd.circular_fingerprint(second)
    assert smsd.circular_fingerprint_counts(first) != smsd.circular_fingerprint_counts(second)


def test_ketone_and_alcohol_have_different_feature_fingerprints():
    ketone, alcohol = smsd.parse_smiles("CC(=O)C"), smsd.parse_smiles("CC(O)C")
    assert smsd.circular_fingerprint(ketone, radius=0, mode="fcfp") != smsd.circular_fingerprint(alcohol, radius=0, mode="fcfp")


def test_charged_valence_hydrogen_fallback():
    assert smsd.implicit_h(7, 0, 1) == 4
    assert smsd.implicit_h(5, 0, -1) == 4
