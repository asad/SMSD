# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
# Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
# See the NOTICE file for attribution, trademark, and algorithm IP terms.
"""Aromaticity must not depend on a fused system's Kekule representation."""
import pytest

import smsd


@pytest.mark.parametrize("first,second", [
    ("C1=CC=C2NC=CC3=NC=CC1=C23", "C1=CC=C2NC=CC3=C2C1=CC=N3"),
    ("C1=CC=C2NC3=CC=CC=C3C3=CC=NC1=C23", "C1=CC=C2NC3=CC=CC=C3C3=C2C1=NC=C3"),
])
def test_fused_resonance_forms(first, second):
    graphs = [smsd.parse_smiles(text) for text in (first, second)]
    assert sum(graphs[0].aromatic) == sum(graphs[1].aromatic)
    counts = [sum(graph.bond_aromatic(first, second) for first, second, _ in graph.bonds())
              for graph in graphs]
    assert counts[0] == counts[1]
    assert smsd.canonical_hash(graphs[0]) == smsd.canonical_hash(graphs[1])


@pytest.mark.parametrize("text", [
    "C1=CC=C2NC=CC3=NC=CC1=C23", "C1=CC=C2NC=CC3=C2C1=CC=N3",
])
def test_fused_donor_ring_is_not_aromatic(text):
    graph = smsd.parse_smiles(text)
    assert sum(graph.aromatic) == 10


@pytest.mark.parametrize("text,count", [
    ("C1=CC=CC=C1", 6),
    ("C1=CC=C2C=CC=CC2=C1", 10),
    ("O=C1C=CC=CC1", 0),
    ("O=C1C=CC=CC=C1", 7),
    ("C1=CC=[O+]C=C1", 6),
    ("C1=CC=C[NH2+]C=C1", 0),
])
def test_external_carbonyl_and_charge_controls(text, count):
    assert sum(smsd.parse_smiles(text).aromatic) == count


@pytest.mark.parametrize("text", [
    "C=C1C=NC2=CNN=NC2=N1",
    "C1=CC(=O)C(=C2C=CN=C3NC(=N)NC=C23)C(=O)O1",
])
def test_external_carbon_double_bond_round_trip(text):
    graph = smsd.parse_smiles(text)
    assert sum(graph.aromatic) == 6
    expected = smsd.canonical_hash(graph)
    assert smsd.canonical_hash(smsd.parse_smiles(smsd.to_smiles(graph))) == expected
    assert graph.kekulize()
    assert smsd.canonical_hash(smsd.parse_smiles(smsd.to_smiles(graph))) == expected


@pytest.mark.parametrize("text,count", [
    ("C1=CC2=C3C=CC(=CC4=[O+]C(=C5C=CC(=CC1=[O+]2)O5)C=C4)O3", 20),
    ("CCC1C2=CC3=NC4=C(C(=O)OC(=O)C4=C3C)C3=NC(=CC4=NC(=CC(=N2)C1C)C(C(C)=O)=C4C)C(C)C3CCC(=O)OC", 22),
])
def test_aromatic_oxygen_donation_survives_round_trip(text, count):
    graph = smsd.parse_smiles(text)
    assert sum(graph.aromatic) == count
    assert smsd.canonical_hash(smsd.parse_smiles(smsd.to_smiles(graph))) == smsd.canonical_hash(graph)


@pytest.mark.parametrize("extra_atoms", [0, 200])
def test_aromatic_conversion_preserves_specified_stereo(extra_atoms):
    padding = "." + "C" * extra_atoms if extra_atoms else ""
    graph = smsd.parse_smiles(
        r"C1=CC2=N/C1=C\C1=N/C(=C\C3=NC4=C(C5=N/C(=C\2)CC5)CC(=O)C4=C3)C=C1.F/C=C/F" + padding)
    descriptors = smsd.assign_ez(graph)
    for focus, expected in {(4, 5): "Z", (8, 9): "Z", (16, 17): "Z", (28, 29): "E"}.items():
        assert descriptors[focus] == expected
    assert smsd.canonical_hash(smsd.parse_smiles(smsd.to_smiles(graph))) == smsd.canonical_hash(graph)


def test_mol_benzene_has_no_alkene_stereo():
    block = """Benzene
  SMSD              2D

  6  6  0  0  0  0  0  0  0  0999 V2000
    1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7500    1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7500    1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7500   -1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7500   -1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  2  0
  3  4  1  0
  4  5  2  0
  5  6  1  0
  6  1  2  0
M  END
"""
    graph = smsd.read_mol_block(block)
    assert sum(graph.aromatic) == 6
    assert smsd.canonical_hash(graph) == smsd.canonical_hash(smsd.parse_smiles("c1ccccc1"))
