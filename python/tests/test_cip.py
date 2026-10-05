"""
SPDX-License-Identifier: Apache-2.0
Copyright (c) 2018-2026 BioInception PVT LTD
Algorithm Copyright (c) 2009-2026 Syed Asad Rahman

CIP stereo descriptor tests for Python bindings.
"""

import pytest
import smsd


class TestCIPRS:
    """R/S assignment on tetrahedral stereocentres."""

    def test_l_alanine_is_S(self):
        result = smsd.assign_rs("N[C@@H](C)C(=O)O")
        assert result.get(1) == 'S'

    def test_d_alanine_is_R(self):
        result = smsd.assign_rs("N[C@H](C)C(=O)O")
        assert result.get(1) == 'R'

    def test_bromochlorofluoromethane_S(self):
        result = smsd.assign_rs("[C@@H](Br)(Cl)F")
        assert result.get(0) == 'S'

    def test_bromochlorofluoromethane_R(self):
        result = smsd.assign_rs("[C@H](Br)(Cl)F")
        assert result.get(0) == 'R'

    def test_l_cysteine_is_R(self):
        result = smsd.assign_rs("N[C@@H](CS)C(=O)O")
        assert result.get(1) == 'R'

    def test_r_glyceraldehyde(self):
        result = smsd.assign_rs("OC[C@@H](O)C=O")
        assert result.get(2) == 'R'

    def test_no_stereocentre(self):
        result = smsd.assign_rs("C")
        assert len(result) == 0

    def test_citric_acid_no_stereocentre(self):
        result = smsd.assign_rs("OC(=O)CC(O)(CC(=O)O)C(=O)O")
        assert len(result) == 0


class TestCIPEZ:
    """E/Z assignment on double bonds."""

    def test_E_2_butene(self):
        result = smsd.assign_ez("C/C=C/C")
        assert (1, 2) in result
        assert result[(1, 2)] == 'E'

    def test_Z_2_butene(self):
        result = smsd.assign_ez("C/C=C\\C")
        assert (1, 2) in result
        assert result[(1, 2)] == 'Z'

    def test_no_EZ_single_bond(self):
        result = smsd.assign_ez("CC")
        assert len(result) == 0


class TestAssignCIP:
    """Batch CIP assignment using assign_cip()."""

    def test_assign_cip(self):
        desc = smsd.assign_cip("N[C@@H](C)C(=O)O")
        assert desc.rs_labels[1] == smsd.RSLabel.S
        assert len(desc.ez_bonds) == 0

    def test_repr(self):
        desc = smsd.assign_cip("N[C@@H](C)C(=O)O")
        r = repr(desc)
        assert "stereocentres=1" in r

    def test_assign_cip_mol(self):
        mol = smsd.parse_smiles("N[C@@H](C)C(=O)O")
        desc = smsd.assign_cip(mol)
        assert desc.rs_labels[1] == smsd.RSLabel.S


class TestRawBindings:
    """Test the raw C++ bindings directly."""

    def test_raw_assign_rs(self):
        from smsd._smsd import assign_rs
        mol = smsd.parse_smiles("N[C@@H](C)C(=O)O")
        label = assign_rs(mol, 1)
        assert label == smsd.RSLabel.S

    def test_raw_assign_ez(self):
        from smsd._smsd import assign_ez
        mol = smsd.parse_smiles("C/C=C/C")
        label = assign_ez(mol, 1, 2)
        assert label == smsd.EZLabel.E

    def test_raw_assign_cip_from_smiles(self):
        from smsd._smsd import assign_cip_from_smiles
        desc = assign_cip_from_smiles("N[C@@H](C)C(=O)O")
        assert desc.rs_labels[1] == smsd.RSLabel.S
