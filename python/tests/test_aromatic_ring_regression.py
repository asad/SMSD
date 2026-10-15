# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
# Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
"""Aromaticity and all-or-nothing Kekule conversion regressions."""
import pytest
import smsd


@pytest.mark.parametrize("text", [
    "c1cc[o+]cc1", "c1cc[s+]cc1", "c1cc[se+]cc1", "c1cc[te+]cc1",
    "c1ccc2cc[o+]cc2c1", "Cc1cc(C)[o+]c(C)c1", "c1ccoc1", "c1ccsc1", "c1cocn1",
])
def test_aromatic_chalcogen_kekule_roundtrip(text):
    graph = smsd.parse_smiles(text)
    before = smsd.canonical_hash(graph)
    smsd.kekulize(graph)
    assert not any(graph.aromatic)
    smsd.perceive_aromaticity(graph)
    assert smsd.canonical_hash(graph) == before


def test_failed_kekule_conversion_is_atomic():
    graph = smsd.parse_smiles("c1ccccc1.c1cccc1")
    before = smsd.canonical_hash(graph)
    aromatic = list(graph.aromatic)
    bonds = [(first, second, graph.bond_order(first, second), graph.bond_aromatic(first, second))
             for first in range(graph.n) for second in graph.neighbors[first] if first < second]
    with pytest.raises(ValueError):
        smsd.kekulize(graph)
    assert graph.aromatic == aromatic
    assert smsd.canonical_hash(graph) == before
    assert bonds == [(first, second, graph.bond_order(first, second), graph.bond_aromatic(first, second))
                     for first in range(graph.n) for second in graph.neighbors[first] if first < second]


@pytest.mark.parametrize("text", [
    "C1=CC=C[NH2+]C=C1", "C[N+]1(C)C=CC=CC=C1",
    "C1=CC=C[S+](C)C=C1", "C1=CC=C[O+](C)C=C1",
])
def test_saturated_cation_interrupts_conjugation(text):
    graph = smsd.parse_smiles(text)
    assert not any(graph.aromatic)


@pytest.mark.parametrize("aromatic,kekule,atoms,bonds", [
    ("c1ccc2cc3cc4cc5cc6ccccc6cc5cc4cc3cc2c1", "C1=CC=C2C=C3C=C4C=C5C=C6C=CC=CC6=CC5=CC4=CC3=CC2=C1", 26, 31),
    ("c1ccc2cc3cc4cc5cc6cc7ccccc7cc6cc5cc4cc3cc2c1", "C1=CC=C2C=C3C=C4C=C5C=C6C=C7C=CC=CC7=CC6=CC5=CC4=CC3=CC2=C1", 30, 36),
    ("c1cc2cccc3c4cccc5cccc(c(c1)c23)c54", "C1=CC2=CC=CC3=C2C(=C1)C1=CC=CC2=C1C3=CC=C2", 20, 24),
    ("c1ccc2c(c1)c1cccc3c4c5ccccc5c5cccc(c2c13)c54", "C1=CC=C2C(=C1)C1=CC=CC3=C1C2=C1C=CC=C2C4=CC=CC=C4C3=C21", 26, 32),
    ("c1ccc2cc3c(cc2c1)CCCCCCCCCCCCCCC3", "C1=CC=C2C=C3CCCCCCCCCCCCCCCC3=CC2=C1", 10, 11),
])
def test_kekule_and_aromatic_forms_have_same_identity(aromatic, kekule, atoms, bonds):
    expected = smsd.parse_smiles(aromatic)
    graph = smsd.parse_smiles(kekule)
    assert sum(graph.aromatic) == atoms
    assert sum(graph.bond_aromatic(first, second) for first in range(graph.n)
               for second in graph.neighbors[first] if first < second) == bonds
    assert smsd.canonical_hash(graph) == smsd.canonical_hash(expected)
    smsd.kekulize(graph)
    smsd.perceive_aromaticity(graph)
    assert smsd.canonical_hash(graph) == smsd.canonical_hash(expected)
