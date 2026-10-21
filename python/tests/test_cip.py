# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
# Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
# See the NOTICE file for attribution, trademark, and algorithm IP terms.
"""CIP stereo descriptor tests for the Python bindings."""

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


class TestCIPSequenceRules4And5:
    """Auxiliary descriptors, like and unlike pairs, and pseudoasymmetric centres."""

    @staticmethod
    def labels(smiles):
        return {atom: label.name for atom, label in enumerate(smsd.assign_cip(smiles).rs_labels)
                if label != smsd.RSLabel.NONE}

    def test_1_4_disubstituted_cyclohexane_is_pseudoasymmetric(self):
        assert self.labels("C[C@H]1CC[C@H](O)CC1") == {1: 'r', 4: 'r'}
        assert self.labels("C[C@H]1CC[C@@H](O)CC1") == {1: 's', 4: 's'}

    def test_acyclic_pseudoasymmetric_centre(self):
        assert self.labels("C[C@H](Cl)[C@H](O)[C@@H](C)Cl") == {1: 'S', 3: 's', 5: 'R'}
        assert self.labels("C[C@H](Cl)[C@@H](O)[C@@H](C)Cl") == {1: 'S', 3: 'r', 5: 'R'}

    def test_rule_4b_like_pairs_give_a_true_centre(self):
        labels = self.labels("CC(C)[C@H]1CC[C@H](C(=O)N[C@H](Cc2ccccc2)C(=O)[O-])CC1")
        assert labels == {3: 'R', 6: 'r', 10: 'R'}

    def test_tropane_c3_reads_the_bridgehead_descriptors(self):
        assert self.labels("C[NH+]1[C@@H]2CC[C@H]1C[C@@H](OC(=O)C(CO)c1ccccc1)C2") == {2: 'R', 5: 'S', 7: 'r'}
        assert self.labels("CN1[C@@H]2CC[C@H]1C[C@H](OC(=O)C(CO)c1ccccc1)C2") == {2: 'R', 5: 'S', 7: 's'}

    def test_assign_rs_reports_true_centres_only(self):
        assert smsd.assign_rs("C[C@H](Cl)[C@H](O)[C@@H](C)Cl") == {1: 'S', 5: 'R'}

    def test_symmetric_centre_stays_unlabelled(self):
        assert self.labels("C[C@H](C)O") == {}

    def test_unlabelled_atom_counts_its_standard_atomic_weight(self):
        # Rule 2: 18F (18) ranks below unlabelled F (18.998) in both passes
        assert self.labels("C[C@H](F)[C@H](Cl)[C@H](C)[18F]") == {1: 'S', 3: 'R', 5: 'S'}
        assert self.labels("[18F][C@H](F)[C@H](O)[C@H](F)[18F]") == {1: 'R', 3: 'r', 5: 'S'}


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
