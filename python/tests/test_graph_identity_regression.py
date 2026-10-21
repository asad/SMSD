# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
# Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
# See the NOTICE file for attribution, trademark, and algorithm IP terms.
"""Chemical canonical identity and macrocycle regression checks."""
import pytest
import smsd


@pytest.mark.parametrize("first,second", [
    ("CC(=O)O", "CC(O)=O"),
    ("C=CC", "CC=C"),
    ("[13CH3]C", "C[13CH3]"),
    ("C[NH2+]CC", "CC[NH2+]C"),
    ("[C@H](F)(Cl)Br", "[C@@H](Cl)(F)Br"),
    ("C[C@H](N)C(=O)O", "N[C@@H](C)C(=O)O"),
    ("Cl/C=C/Cl", "C(=C/Cl)\\Cl"),
    ("C[C@H](F)[C@H](F)C", "C[C@@H](F)[C@@H](F)C"),
    ("C/C=C/C=C\\C", "C/C=C\\C=C\\C"),
    ("C/C=C(/F)\\C(Cl)=C/C", "FC(/C(=C\\C)Cl)=C/C"),
    ("C[C@H](F)Cl.C[C@@H](F)Cl", "C[C@@H](F)Cl.C[C@H](F)Cl"),
    ("C1[C@H](F)CCC[C@@H]1F", "F[C@H]1CCC[C@@H](F)C1"),
])
def test_equivalent_smiles_have_one_canonical_identity(first, second):
    first_graph = smsd.parse_smiles(first)
    second_graph = smsd.parse_smiles(second)
    assert smsd.canonical_hash(first_graph) == smsd.canonical_hash(second_graph)
    assert smsd.to_smiles(first_graph) == smsd.to_smiles(second_graph)
    assert smsd.same_canonical_graph(first_graph, second_graph)


@pytest.mark.parametrize("first,second", [
    ("C[C@H](N)C(=O)O", "C[C@@H](N)C(=O)O"),
    ("Cl/C=C/Cl", "Cl/C=C\\Cl"),
])
def test_stereoisomers_retain_different_hashes(first, second):
    assert smsd.canonical_hash(smsd.parse_smiles(first)) != smsd.canonical_hash(smsd.parse_smiles(second))


def test_orbits_respect_chemical_properties_and_ring_components():
    propene = smsd.parse_smiles("C=CC")
    smsd.canonical_hash(propene)
    assert propene.orbit[0] != propene.orbit[2]
    rings = smsd.parse_smiles("C1CC1.C1CCCCC1")
    smsd.canonical_hash(rings)
    assert len(set(rings.orbit[:3])) == 1
    assert len(set(rings.orbit[3:])) == 1
    assert rings.orbit[0] != rings.orbit[3]
    assert not rings.automorphism_generators_truncated()


@pytest.mark.parametrize("ring_size", [21, 22, 33])
def test_macrocycles_keep_ring_membership_and_basis(ring_size):
    graph = smsd.parse_smiles("C1" + "C" * (ring_size - 1) + "1")
    assert sum(graph.ring) == ring_size
    assert graph.ring_system_count() == 1
    assert smsd.Substructure.smarts(f'[r{ring_size}]').match(graph).found
    rings = smsd.compute_sssr(graph)
    assert len(rings) == 1 and len(rings[0]) == ring_size


def test_cyclic_peptide_keeps_its_macrocycle():
    graph = smsd.parse_smiles("O=C1NCC(=O)NCC(=O)NCC(=O)NCC(=O)NCC(=O)NCC(=O)NC1")
    assert sum(graph.ring) == 21
    assert smsd.Substructure.smarts('[r21]').match(graph).found
