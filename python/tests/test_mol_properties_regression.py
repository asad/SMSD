# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
# Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
"""MOL interchange preserves atom properties without adding query constraints."""
import pytest
import smsd


@pytest.mark.parametrize("writer", [smsd.write_mol_block, smsd.write_mol_block_v3000])
@pytest.mark.parametrize("text", ["CCO", "OCC(=O)O", "c1ccccc1", "[nH]1cccc1", "[NH4+]", "[H]OC", "[CH3]"])
def test_mol_hydrogens_round_trip(writer, text):
    original = smsd.parse_smiles(text)
    block = writer(original)
    assert "HCOUNT=" not in block
    restored = smsd.read_mol_block(block)
    assert restored.hydrogen_count == original.hydrogen_count
    assert restored.atomic_num == original.atomic_num


def test_mol_acid_hydrogens_match_smarts():
    restored = smsd.read_mol_block(smsd.write_mol_block(smsd.parse_smiles("OCC(=O)O")))
    assert smsd.Substructure.smarts('[OH1][CH2][CX3](=O)[OX2H1]').match(restored).found


@pytest.mark.parametrize("count,expected", [(2, 2), (-1, 0), (0, 3)])
def test_v3000_hydrogen_count_uses_direct_value(count, expected):
    block = f"""Ammonia


  0  0  0  0  0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 1 0 0 0 0
M  V30 BEGIN ATOM
M  V30 1 N 0 0 0 0 HCOUNT={count}
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 END BOND
M  V30 END CTAB
M  END
"""
    assert smsd.read_mol_block(block).hydrogen_count == [expected]


@pytest.mark.parametrize("text", ["c1ccc(cc1)c1ccccc1", "c1ccccc1c1ccncc1"])
def test_implicit_biaryl_bonds_are_single(text):
    molecule = smsd.parse_smiles(text)
    bridge_orders = [order for first, second, order in molecule.bonds()
                     if not molecule.bond_in_ring(first, second)]
    assert bridge_orders == [1]


@pytest.mark.parametrize("v3000", [False, True])
def test_mol_hydrogen_isotope_aliases(v3000):
    molecule = smsd.parse_smiles("[H]O[H]")
    writer = smsd.write_mol_block_v3000 if v3000 else smsd.write_mol_block
    block = writer(molecule)
    if v3000:
        block = block.replace("M  V30 1 H ", "M  V30 1 D ").replace("M  V30 3 H ", "M  V30 3 T ")
    else:
        lines = block.splitlines()
        lines[4] = lines[4][:31] + "D  " + lines[4][34:]
        lines[6] = lines[6][:31] + "T  " + lines[6][34:]
        block = "\n".join(lines) + "\n"
    assert smsd.read_mol_block(block).mass_number == [2, 0, 3]


def test_sdf_trailing_blank_lines_do_not_add_records(tmp_path):
    molecule = smsd.parse_smiles("CCO")
    path = tmp_path / "molecules.sdf"
    path.write_text((smsd.write_mol_block(molecule) + "$$$$\n") * 2 + "\n \t\n")
    assert [item.n for item in smsd.read_sdf(str(path))] == [3, 3]


def test_large_mol_output_uses_v3000():
    molecule = smsd.parse_smiles("C" * 1000)
    block = smsd.write_mol_block(molecule)
    assert "V3000" in block.splitlines()[3]
    assert smsd.read_mol_block(block).n == 1000


def test_v3000_continuation_preserves_atom_attributes():
    block = smsd.write_mol_block_v3000(smsd.parse_smiles("N"))
    block = block.replace("M  V30 1 N 0 0 0 0", "M  V30 1 N 0 0 -\nM  V30 0 0 CHG=1")
    assert smsd.read_mol_block(block).formal_charge == [1]


@pytest.mark.parametrize("bond_type,expected_bonds", [(8, 1), (9, 1), (10, 0)])
def test_v3000_extended_bond_types_are_not_truncated(bond_type, expected_bonds):
    block = smsd.write_mol_block_v3000(smsd.parse_smiles("CO"))
    block = block.replace("M  V30 1 1 1 2", f"M  V30 1 {bond_type} 1 2")
    restored = smsd.read_mol_block(block)
    assert len(restored.bonds()) == expected_bonds
    if expected_bonds:
        assert restored.bond_order(0, 1) == 1
