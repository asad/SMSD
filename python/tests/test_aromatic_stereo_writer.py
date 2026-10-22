# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
# Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
# See the NOTICE file for attribution, trademark, and algorithm IP terms.
"""Directional aromatic carriers retain exocyclic double-bond configurations."""
import pytest
import smsd
import math


@pytest.mark.parametrize("text", [r"F/C=c1\cnc2c(n1)N=NNC=2", "F/C=c1/cnc2c(n1)N=NNC=2"])
def test_aromatic_reference_carrier_preserves_mapped_stereo(text):
    original = smsd.parse_smiles(text)
    output = smsd.to_smiles(original)
    assert output == smsd.to_smiles(original)
    restored = smsd.parse_smiles(output)
    assert smsd.canonical_hash(restored) == smsd.canonical_hash(original)
    options = smsd.Substructure.Options()
    options.chemistry.use_bond_stereo = True
    result = smsd.Substructure.find(original).with_options(options).match(restored)
    assert result.found and result.optimal and result.complete
    validation = smsd.MCS.Options()
    validation.chemistry.use_bond_stereo = True
    assert smsd.MCS.validate_mapping(original, restored, result.best, validation) == []


def test_opposite_aromatic_reference_configurations_remain_distinct():
    first = smsd.parse_smiles(r"F/C=c1\cnc2c(n1)N=NNC=2")
    second = smsd.parse_smiles("F/C=c1/cnc2c(n1)N=NNC=2")
    assert smsd.canonical_hash(first) != smsd.canonical_hash(second)
    options = smsd.Substructure.Options()
    options.chemistry.use_bond_stereo = True
    result = smsd.Substructure.find(first).with_options(options).match(second)
    assert result.optimal and result.complete and not result.found
    equivalent = smsd.parse_smiles(r"F\C=c1/cnc2c(n1)N=NNC=2")
    assert smsd.to_smiles(first) == smsd.to_smiles(equivalent)


@pytest.mark.parametrize("text,expected", [(r"F/C1=C(\Cl)C=CC=C1", "[c]1(=[c](cccc1)/F)/Cl"),
                                           (r"F/C1=C(\Cl)C=CC=CC=CC=C1", "[c]1(=[c](cccccccc1)/F)/Cl"),
                                           ("F/C1=C(/Cl)C=CC=CC=CC=C1", r"[c]1(=[c](cccccccc1)\F)/Cl")])
def test_aromatic_focus_stereo_is_written_on_single_bonds(text, expected):
    # '/' and '\' are single bonds in OpenSMILES, so the substituent bonds carry the configuration.
    # A ring double bond can be trans only in a ring of eight or more atoms, as in [10]annulene.
    graph = smsd.parse_smiles(text)
    output = smsd.to_smiles(graph)
    assert output == expected
    assert smsd.canonical_hash(smsd.parse_smiles(output)) == smsd.canonical_hash(graph)


@pytest.mark.parametrize("text", ["F/C=C/F", r"F/C=C\F", "F/C=C/C=C/F", r"F/C=C\C=C/F",
                                 "F/C=C/C(=C/F)Cl", r"F/C=C/C(=C\F)Cl"])
def test_ordinary_conjugated_stereo_still_round_trips(text):
    graph = smsd.parse_smiles(text)
    assert smsd.canonical_hash(smsd.parse_smiles(smsd.to_smiles(graph))) == smsd.canonical_hash(graph)


@pytest.mark.parametrize("size,shortcut,alternating,expected", [(6, None, True, False),
    (7, None, False, True), (8, None, False, True), (8, 5, True, False), (10, 7, False, True)])
def test_mol_geometry_uses_the_shortest_focus_cycle(size, shortcut, alternating, expected):
    lines = ["Cycloalkene", "  SMSD              2D", "",
             f"{size:3d}{size + (shortcut is not None):3d}  0  0  0  0  0  0  0  0999 V2000"]
    for atom in range(size):
        angle = 2 * math.pi * atom / size
        lines.append(f"{math.cos(angle):10.4f}{math.sin(angle):10.4f}{0.:10.4f} C   0  0  0  0  0  0  0  0  0  0  0  0")
    for atom in range(size):
        order = 2 if atom == 0 or (alternating and atom < 6 and atom % 2 == 0) else 1
        lines.append(f"{atom + 1:3d}{(atom + 1) % size + 1:3d}{order:3d}  0")
    if shortcut is not None:
        lines.append(f"  1{shortcut + 1:3d}  1  0")
    lines.append("M  END")
    molecule = smsd.read_mol_block("\n".join(lines) + "\n")
    assert ((0, 1) in smsd.assign_ez(molecule)) == expected


@pytest.mark.parametrize("text", [
    r"C1=CC2=N/C1=C\C1=N/C(=C\C3=NC4=C(C5=N/C(=C\2)CC5)CC(=O)C4=C3)C=C1",
    r"C1=CC2=N/C1=C/C1=N/C(=C\C3=NC4=C(C5=N/C(=C\2)CC5)CC(=O)C4=C3)C=C1",
])
def test_aromatic_focus_preserves_stereo_without_extra_tags(text):
    original = smsd.parse_smiles(text)
    carbon_configurations = {focus: label for focus, label in smsd.assign_ez(original).items()
                             if original.atomic_num[focus[0]] == original.atomic_num[focus[1]] == 6}
    assert len(carbon_configurations) == 3
    output = smsd.to_smiles(original)
    restored = smsd.parse_smiles(output)
    assert smsd.canonical_hash(restored) == smsd.canonical_hash(original)
    options = smsd.Substructure.Options()
    options.chemistry.use_bond_stereo = True
    result = smsd.Substructure.find(original).with_options(options).match(restored)
    assert result.found and result.optimal and result.complete
    equivalent = smsd.parse_smiles(text.translate(str.maketrans("/\\", "\\/")))
    assert smsd.to_smiles(equivalent) == output


@pytest.mark.parametrize("writer", [smsd.write_mol_block, smsd.write_mol_block_v3000])
def test_mol_coordinates_preserve_specified_aromatic_focus_stereo(writer):
    graph = smsd.parse_smiles(r"C1=CC2=N/C1=C\C1=N/C(=C\C3=NC4=C(C5=N/C(=C\2)CC5)CC(=O)C4=C3)C=C1")
    expected = {focus: label for focus, label in smsd.assign_ez(graph).items()
                if graph.atomic_num[focus[0]] == graph.atomic_num[focus[1]] == 6}
    assert len(expected) == 3
    restored = smsd.assign_ez(smsd.read_mol_block(writer(graph)))
    for focus, configuration in expected.items():
        assert restored.get(focus) == configuration


@pytest.mark.parametrize("writer", [smsd.write_mol_block, smsd.write_mol_block_v3000])
@pytest.mark.parametrize("text", [r"C/C1=C\CCCCCC1", "C/C1=C/CCCCCC1"])
def test_mol_cyclic_ligands_remain_on_opposite_sides(writer, text):
    chem = pytest.importorskip("rdkit.Chem")
    original = chem.MolFromSmiles(text)
    graph = smsd.parse_smiles(text)
    before = smsd.canonical_hash(graph)
    restored = chem.MolFromMolBlock(writer(graph), sanitize=True, removeHs=False, strictParsing=True)
    assert restored is not None
    positions = restored.GetConformer()
    first, second = positions.GetAtomPosition(1), positions.GetAtomPosition(2)
    axis_x, axis_y = second.x - first.x, second.y - first.y

    def side(ligand):
        position = positions.GetAtomPosition(ligand)
        return axis_x * (position.y - first.y) - axis_y * (position.x - first.x)

    assert side(0) * side(8) < 0
    assert restored.GetBondBetweenAtoms(1, 2).GetStereo() == original.GetBondBetweenAtoms(1, 2).GetStereo()
    assert smsd.canonical_hash(graph) == before
