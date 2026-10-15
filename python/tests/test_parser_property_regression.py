# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
# Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
"""Regression checks for SMARTS atom properties and lazy canonical data."""
import pytest
import smsd


@pytest.mark.parametrize("query,target", [
    ("Cn1cccc1", "Cn1cccc1"),
    ("CSc1ccccc1", "CSc1ccccc1"),
    ("Ca", "Cc1ccccc1"),
    ("Na", "Nc1ccccc1"),
    ("[#6]Na", "CNc1ccccc1"),
    ("Ac1ccccc1", "Cc1ccccc1"),
    ("ClCBr", "ClCBr"),
    ("[Na+]", "[Na+]"),
    ("as", "c1ccsc1"),
    ("[se]", "c1cc[se]c1"),
])
def test_organic_atom_tokens(query, target):
    assert smsd.Substructure.smarts(query).match(smsd.parse_smiles(target)).found


@pytest.mark.parametrize("smiles,hydrogens", [
    ("c1ccsc1", [1, 1, 1, 0, 1]),
    ("c1ccoc1", [1, 1, 1, 0, 1]),
    ("Cn1cccc1", [3, 0, 1, 1, 1, 1]),
    ("n1ccccc1", [0, 1, 1, 1, 1, 1]),
    ("[nH]1cccc1", [1, 1, 1, 1, 1]),
])
def test_aromatic_hydrogen_counts(smiles, hydrogens):
    assert smsd.parse_smiles(smiles).hydrogen_count == hydrogens


@pytest.mark.parametrize("query,target", [
    ("[sH]", "c1ccsc1"),
    ("[nH]", "Cn1cnc2c1c(=O)n(C)c(=O)n2C"),
])
def test_aromatic_heteroatoms_without_hydrogen(query, target):
    assert [dict(match.atom_mapping) for match in smsd.Substructure.smarts(query).match_all(smsd.parse_smiles(target), 1000).matches] == []


@pytest.mark.parametrize("query,target", [
    ("[OH1]", "[H]OC"),
    ("[CH3]", "[H]C([H])([H])O[H]"),
    ("[OX2H]", "[H]OC"),
    ("[N;!H0]", "[H]N([H])C"),
    ("[N;X3;H2;v3]", "[H]N([H])C"),
    ("[C;v4]", "CC"),
    ("[N;v3]", "CN"),
    ("[c;v4]", "c1ccccc1"),
    ("[c;v4;X3]", "c1ccc2ccccc2c1"),
    ("[n;v3]", "n1ccccc1"),
    ("[nH;v3]", "[nH]1cccc1"),
    ("[n;v3;H0]", "Cn1cccc1"),
    ("[s;v2;H0]", "c1ccsc1"),
    ("[C;v3]", "[CH3]"),
    ("[N+;v4;H4]", "[NH4+]"),
    ("[C;v4;X4;H4]", "[H]C([H])([H])[H]"),
])
def test_total_hydrogen_valence_and_connectivity(query, target):
    assert smsd.Substructure.smarts(query).match(smsd.parse_smiles(target)).found


def test_explicit_hydrogen_mol_properties():
    methanol = smsd.parse_smiles("[H]C([H])([H])O[H]")
    graph = smsd.read_mol_block(smsd.write_mol_block(methanol))
    assert smsd.Substructure.smarts('[CH3;v4;X4][OH1;v2;X2]').match(graph).found


def test_canonical_hash_initialises_fresh_graphs():
    carbon = smsd.parse_smiles("C")
    oxygen = smsd.parse_smiles("O")
    assert smsd.canonical_hash(carbon) != smsd.canonical_hash(oxygen)


def test_canonical_hash_populates_lazy_properties():
    molecule = smsd.parse_smiles("CCO")
    smsd.canonical_hash(molecule)
    assert sorted(molecule.canonical_label) == [0, 1, 2]
    symmetric = smsd.parse_smiles("CC")
    smsd.canonical_hash(symmetric)
    assert len(symmetric.orbit) == 2
    assert symmetric.orbit[0] == symmetric.orbit[1]
