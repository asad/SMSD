# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
# Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
# See the NOTICE file for attribution, trademark, and algorithm IP terms.
"""
Tests for force-directed layout, stress majorisation, and template matching.

Run with:  pytest tests/test_layout.py -v
"""
import math
import pytest
import smsd
from smsd import parse_smiles


# ===========================================================================
# Helpers
# ===========================================================================

def linear_layout(n):
    """Generate a linear initial layout for n atoms."""
    return [[float(i), 0.0] for i in range(n)]


# ===========================================================================
# Force-directed layout
# ===========================================================================

class TestForceDirectedLayout:
    """Test the force-directed 2D layout engine."""

    def test_benzene_stress_non_negative(self):
        mol = parse_smiles("c1ccccc1")
        coords = linear_layout(mol.n)
        stress, new_coords = smsd.force_directed_layout(mol, coords)
        assert stress >= 0.0

    def test_benzene_atoms_separated(self):
        mol = parse_smiles("c1ccccc1")
        coords = linear_layout(mol.n)
        _, new_coords = smsd.force_directed_layout(mol, coords)
        # Check no two atoms overlap
        for i in range(len(new_coords)):
            for j in range(i + 1, len(new_coords)):
                dx = new_coords[i][0] - new_coords[j][0]
                dy = new_coords[i][1] - new_coords[j][1]
                dist = math.sqrt(dx * dx + dy * dy)
                assert dist > 0.01, f"Atoms {i} and {j} overlap"

    def test_naphthalene(self):
        mol = parse_smiles("c1ccc2ccccc2c1")
        coords = linear_layout(mol.n)
        stress, new_coords = smsd.force_directed_layout(mol, coords)
        assert stress >= 0.0
        assert len(new_coords) == mol.n

    def test_ethane(self):
        mol = parse_smiles("CC")
        coords = linear_layout(mol.n)
        stress, new_coords = smsd.force_directed_layout(mol, coords)
        assert stress >= 0.0

    def test_smiles_string_input(self):
        """Accept SMILES string directly."""
        coords = linear_layout(6)
        stress, _ = smsd.force_directed_layout("c1ccccc1", coords)
        assert stress >= 0.0

    def test_custom_parameters(self):
        mol = parse_smiles("c1ccccc1")
        coords = linear_layout(mol.n)
        stress, _ = smsd.force_directed_layout(mol, coords,
                                                max_iter=100,
                                                target_bond_length=2.0)
        assert stress >= 0.0

    def test_short_coordinate_list_is_padded(self):
        mol = parse_smiles("CCC")
        stress, new_coords = smsd.force_directed_layout(mol, [[0.0, 0.0]])
        assert stress >= 0.0
        assert len(new_coords) == mol.n


# ===========================================================================
# Stress majorisation (SMACOF)
# ===========================================================================

class TestStressMajorisation:
    """Test the SMACOF stress majorisation layout engine."""

    def test_benzene_stress_non_negative(self):
        mol = parse_smiles("c1ccccc1")
        coords = linear_layout(mol.n)
        stress, new_coords = smsd.stress_majorisation(mol, coords)
        assert stress >= 0.0

    def test_naphthalene(self):
        mol = parse_smiles("c1ccc2ccccc2c1")
        coords = linear_layout(mol.n)
        stress, new_coords = smsd.stress_majorisation(mol, coords)
        assert stress >= 0.0
        assert len(new_coords) == mol.n

    def test_smiles_string_input(self):
        coords = linear_layout(6)
        stress, _ = smsd.stress_majorisation("c1ccccc1", coords)
        assert stress >= 0.0

    def test_custom_parameters(self):
        mol = parse_smiles("c1ccccc1")
        coords = linear_layout(mol.n)
        stress, _ = smsd.stress_majorisation(mol, coords,
                                              max_iter=50,
                                              target_bond_length=2.0)
        assert stress >= 0.0

    def test_atoms_not_collapsed(self):
        mol = parse_smiles("c1ccc2ccccc2c1")
        coords = linear_layout(mol.n)
        _, new_coords = smsd.stress_majorisation(mol, coords)
        for i in range(len(new_coords)):
            for j in range(i + 1, len(new_coords)):
                dx = new_coords[i][0] - new_coords[j][0]
                dy = new_coords[i][1] - new_coords[j][1]
                dist = math.sqrt(dx * dx + dy * dy)
                assert dist > 0.01, f"Atoms {i} and {j} collapsed"

    def test_short_coordinate_list_is_padded(self):
        mol = parse_smiles("c1ccccc1")
        stress, new_coords = smsd.stress_majorisation(mol, [[0.0, 0.0], [1.0, 0.0]])
        assert stress >= 0.0
        assert len(new_coords) == mol.n


# ===========================================================================
# Template matching
# ===========================================================================

class TestMatchTemplate:
    """Test the scaffold template matching engine."""

    def test_benzene_match(self):
        coords = smsd.match_template("c1ccccc1")
        assert len(coords) == 6

    def test_no_match_chain(self):
        coords = smsd.match_template("CCCCCCCC")
        assert len(coords) == 0

    def test_scaling(self):
        c1 = smsd.match_template("c1ccccc1", target_bond_length=1.0)
        c2 = smsd.match_template("c1ccccc1", target_bond_length=2.0)
        assert len(c1) == 6
        assert len(c2) == 6
        # Find a non-zero coordinate to test scaling
        for k in range(6):
            if abs(c1[k][0]) > 0.01:
                ratio = abs(c2[k][0]) / abs(c1[k][0])
                assert abs(ratio - 2.0) < 0.3
                break

    def test_template_coords_valid(self):
        """Template coordinates should have finite values."""
        coords = smsd.match_template("c1ccccc1")
        for x, y in coords:
            assert math.isfinite(x), "x coordinate must be finite"
            assert math.isfinite(y), "y coordinate must be finite"


# ===========================================================================
# Integration: pipeline tests
# ===========================================================================

class TestLayoutPipeline:
    """Test combined layout pipelines."""

    def test_force_then_crossing_reduction(self):
        mol = parse_smiles("c1ccc2ccccc2c1")  # naphthalene
        coords = linear_layout(mol.n)
        _, coords = smsd.force_directed_layout(mol, coords, max_iter=200)
        crossings, coords = smsd.reduce_crossings(mol, coords, max_iter=500)
        assert crossings >= 0

    def test_smacof_then_crossing_reduction(self):
        mol = parse_smiles("c1ccc2ccccc2c1")
        coords = linear_layout(mol.n)
        _, coords = smsd.stress_majorisation(mol, coords, max_iter=100)
        crossings, coords = smsd.reduce_crossings(mol, coords, max_iter=500)
        assert crossings >= 0

    def test_crossing_reduction_short_coordinate_list_is_padded(self):
        mol = parse_smiles("c1ccc2ccccc2c1")
        crossings, coords = smsd.reduce_crossings(mol, [[0.0, 0.0], [1.0, 0.0]])
        assert crossings >= 0
        assert len(coords) == mol.n

    def test_template_or_fallback(self):
        """Use template if available, otherwise fall back to force-directed."""
        smiles = "c1ccccc1"
        mol = parse_smiles(smiles)
        coords = smsd.match_template(mol)
        if not coords:
            init = linear_layout(mol.n)
            _, coords = smsd.force_directed_layout(mol, init)
        assert len(coords) == mol.n
