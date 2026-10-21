# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
# Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
# See the NOTICE file for attribution, trademark, and algorithm IP terms.
"""Stereochemistry must survive SMILES and MOL interchange."""

import pytest
import smsd


def _rs_labels(molecule):
    return sorted(smsd.assign_rs(molecule).values())


def _ez_labels(molecule):
    return sorted(smsd.assign_ez(molecule).values())


@pytest.mark.parametrize("text", [
    "C[C@@H](N)C(=O)O", "C[C@H](N)C(=O)O", "CC[C@@H](O)C",
    "[C@@H](Br)(Cl)F", "C[C@H]1CCCCO1", "CC1.[C@H]1(N)O",
])
def test_smiles_preserves_tetrahedral_stereo(text):
    molecule = smsd.parse_smiles(text)
    assert _rs_labels(smsd.to_smiles(molecule)) == _rs_labels(molecule)


def test_tetrahedral_classes():
    assert smsd.assign_rs("N[C@TH2H](C)C(=O)O") == {1: "S"}
    assert smsd.assign_rs("N[C@TH1H](C)C(=O)O") == {1: "R"}
    with pytest.raises(ValueError):
        smsd.parse_smiles("[Pt@SP1](Cl)(F)(Br)I")


@pytest.mark.parametrize("text", [r"C(\F)=C/F", r"C(=C/Cl)\Cl"])
def test_branch_directional_stereo(text):
    assert _ez_labels(text) == ["E"]


@pytest.mark.parametrize("closing,opening,expected", [
    (r"C1CCCCC/C=C\1", "C/1CCCCC/C=C1", "Z"),
    ("C1CCCCC/C=C/1", r"C\1CCCCC/C=C1", "E"),
])
def test_ring_closure_directional_orientation(closing, opening, expected):
    molecule = smsd.parse_smiles(closing)
    equivalent = smsd.parse_smiles(opening)
    assert smsd.assign_ez(molecule) == {(6, 7): expected}
    assert molecule.get_canonical_hash() == equivalent.get_canonical_hash()
    restored = smsd.parse_smiles(smsd.to_smiles(molecule))
    assert molecule.get_canonical_hash() == restored.get_canonical_hash()


@pytest.mark.parametrize("text", [
    "C/C=C/C", "C/C=C\\C", "F/C=C(Cl)/C", "C/C=C/C=C\\C",
    "C/C=C\\C=C/C", "C/C=C\\C=C\\C", "C1/C=C\\CCCC1", "C/C=C(/F)\\C(Cl)=C/C",
])
def test_smiles_preserves_double_bond_stereo(text):
    molecule = smsd.parse_smiles(text)
    assert _ez_labels(smsd.to_smiles(molecule)) == _ez_labels(molecule)


@pytest.mark.parametrize("writer", [smsd.write_mol_block, smsd.write_mol_block_v3000])
@pytest.mark.parametrize("text", [
    "N[C@@H](C)C(=O)O", "[C@@H](Br)(Cl)F", "[H][C@](Br)(Cl)F",
    "C/C=C/C", "C/C=C\\C", "C/C=C/C=C\\C", "C/C=C/[C@H](F)Cl",
])
def test_mol_preserves_stereo(writer, text):
    molecule = smsd.parse_smiles(text)
    restored = smsd.read_mol_block(writer(molecule))
    assert _rs_labels(restored) == _rs_labels(molecule)
    assert _ez_labels(restored) == _ez_labels(molecule)


def test_known_mol_parity_and_bond_order():
    atom_lines = """M  V30 1 N 0 0 0 0
M  V30 2 C 0 0 0 0 CFG=1
M  V30 3 C 0 0 0 0
M  V30 4 C 0 0 0 0
M  V30 5 O 0 0 0 0
M  V30 6 O 0 0 0 0"""
    for first_bonds in ["M  V30 1 1 1 2\nM  V30 2 1 2 3", "M  V30 1 1 2 3\nM  V30 2 1 1 2"]:
        block = f"""Alanine


  0  0  0  0  0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 6 5 0 0 0
M  V30 BEGIN ATOM
{atom_lines}
M  V30 END ATOM
M  V30 BEGIN BOND
{first_bonds}
M  V30 3 1 2 4
M  V30 4 2 4 5
M  V30 5 1 4 6
M  V30 END BOND
M  V30 END CTAB
M  END
"""
        assert smsd.assign_rs(smsd.read_mol_block(block)) == {1: "S"}


@pytest.mark.parametrize("either", [False, True])
def test_mol_double_bond_coordinates_and_either(either):
    field = "  3" if either else "  0"
    block = """Butene


  4  3  0  0  0  0  0  0  0  0999 V2000
   -1.0000    1.0000    0.0000 C   0  0  0
    0.0000    0.0000    0.0000 C   0  0  0
    1.0000    0.0000    0.0000 C   0  0  0
    2.0000   -1.0000    0.0000 C   0  0  0
  1  2  1  0
  2  3  2""" + field + "\n  3  4  1  0\nM  END\n"
    assert smsd.assign_ez(smsd.read_mol_block(block)) == ({} if either else {(1, 2): "E"})


@pytest.mark.parametrize("text", [
    "[CH3:1]/[CH:2]=[C:3](/[F:4])\\[C:5]([Cl:6])=[CH:7]/[CH3:8]",
    "[CH3:1]/[CH:2]=[C:3](\\[F:4])/[C:5]([Cl:6])=[CH:7]/[CH3:8]",
])
def test_conjugated_branch_stereo_by_bond(text):
    def configurations(molecule):
        classes = molecule.atom_class
        return {
            tuple(sorted((classes[first], classes[second]))): label
            for (first, second), label in smsd.assign_ez(molecule).items()
        }

    molecule = smsd.parse_smiles(text)
    expected = configurations(molecule)
    assert len(expected) == 2
    restored = smsd.parse_smiles(smsd.to_smiles(molecule))
    assert configurations(restored) == expected
