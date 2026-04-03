# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
# Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
# See the NOTICE file for attribution, trademark, and algorithm IP terms.
"""
Tests for implicit / explicit hydrogen handling in the SMSD Python bindings.

Key invariants
--------------
1. Implicit-H SMILES and bracket-H SMILES give the same heavy-atom MCS.
2. Explicit [H] graph nodes are real atoms — a query with them cannot match
   an implicit-H target of smaller size.
3. After parsing, [CH4], C, and [C@@H] all normalise to the same
   single-carbon heavy-atom representation.
4. Ring aromaticity and formal charge are unaffected by H-count notation.

Run with:  pytest python/tests/test_hydrogen_handling.py -v
"""
import pytest
from smsd import parse_smiles, find_mcs, is_substructure, ChemOptions, MCSOptions


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def mcs_size(smi1: str, smi2: str, **opts_kwargs) -> int:
    g1 = parse_smiles(smi1)
    g2 = parse_smiles(smi2)
    chem = ChemOptions()
    for k, v in opts_kwargs.items():
        setattr(chem, k, v)
    mo = MCSOptions()
    mo.timeout_ms = 10_000
    return len(find_mcs(g1, g2, chem, mo))


def is_sub(query: str, target: str) -> bool:
    q = parse_smiles(query)
    t = parse_smiles(target)
    return is_substructure(q, t, ChemOptions(), 10_000)


# ---------------------------------------------------------------------------
# 1. Implicit-H normalisation
# ---------------------------------------------------------------------------

class TestImplicitH:
    """Bracket-H notation must not change heavy-atom counts or MCS size."""

    def test_methane_implicit_vs_bracket(self):
        # [CH4] and C are the same molecule after parsing
        assert is_sub("[CH4]", "C")
        assert is_sub("C", "[CH4]")

    def test_methane_mcs_size_one(self):
        assert mcs_size("C", "[CH4]") == 1

    def test_benzene_implicit_vs_bracketed_aromatic(self):
        # [cH]1[cH][cH][cH][cH][cH]1 is benzene with explicit aromatic H counts
        assert mcs_size("c1ccccc1", "[cH]1[cH][cH][cH][cH][cH]1") == 6

    def test_ethanol_bracketed_vs_implicit(self):
        # [CH3][CH2][OH] is the same as CCO
        assert mcs_size("CCO", "[CH3][CH2][OH]") == 3

    def test_aspirin_heavy_atom_count_preserved(self):
        # One CH3 written with brackets — MCS must still cover all 13 heavy atoms
        assert mcs_size(
            "CC(=O)Oc1ccccc1C(=O)O",
            "[CH3]C(=O)Oc1ccccc1C(=O)O"
        ) == 13

    def test_water_implicit_vs_OH2(self):
        assert is_sub("O", "[OH2]")
        assert is_sub("[OH2]", "O")


# ---------------------------------------------------------------------------
# 2. Explicit [H] as graph nodes
# ---------------------------------------------------------------------------

class TestExplicitHNodes:
    """[H] atoms written as graph nodes are real atoms in the graph."""

    def test_explicit_H_methane_atom_count(self):
        # [H]C([H])([H])[H] has 5 atoms in the graph (C + 4 H)
        mol = parse_smiles("[H]C([H])([H])[H]")
        assert mol.n == 5

    def test_explicit_H_query_does_not_match_implicit_target(self):
        # 5-atom explicit-H methane cannot be substructure of 1-atom implicit methane
        assert not is_sub("[H]C([H])([H])[H]", "C")

    def test_mcs_explicit_H_vs_implicit_is_heavy_atom_only(self):
        # [H]O[H] (3 atoms) vs O (1 atom) — MCS is just the O
        assert mcs_size("[H]O[H]", "O") == 1

    def test_h2_self_match(self):
        # Hydrogen molecule matches itself
        assert is_sub("[H][H]", "[H][H]")


# ---------------------------------------------------------------------------
# 3. H-count in substructure matching
# ---------------------------------------------------------------------------

class TestHCountSubstructure:
    """Implicit H counts affect atom compatibility under strict mode."""

    def test_water_substructure_of_itself(self):
        assert is_sub("O", "O")

    def test_ammonia_matches_aniline_N_loose(self):
        # In LOOSE bond-order mode, N (NH3) should match the aniline nitrogen
        aniline = parse_smiles("Nc1ccccc1")
        nh3     = parse_smiles("N")
        chem = ChemOptions()
        from smsd import BondOrderMode
        chem.match_bond_order = BondOrderMode.LOOSE
        assert is_substructure(nh3, aniline, chem, 10_000)

    def test_caffeine_vs_theophylline_nH_notation(self):
        # [nH] vs n — MCS should be ≥ 11 atoms
        assert mcs_size(
            "Cn1cnc2c1c(=O)n(C)c(=O)n2C",
            "Cn1cnc2c1c(=O)[nH]c(=O)n2C"
        ) >= 11

    def test_cyclohexane_CH2_notation(self):
        # [CH2] bracket notation vs implicit — full 6-atom MCS
        assert mcs_size("C1CCCCC1", "[CH2]1[CH2][CH2][CH2][CH2][CH2]1") == 6

    def test_morphine_stereo_notation_mcs(self):
        # Stereo/H notation variants of morphine must give the same heavy-atom MCS
        plain  = "CN1CCC23C4C1CC5=C(C2C(C=C4)O3)C=C(C=C5)O"
        stereo = "C[N]1CC[C@@]23[C@@H]4[C@H]1C[C@@H]5=C([C@@H]2[C@H](C=C4)O3)C=C(C=C5)O"
        assert mcs_size(plain, stereo) >= 17

    def test_ibuprofen_stereo_h_notation(self):
        plain  = "CC(C)Cc1ccc(CC(C)C(=O)O)cc1"
        stereo = "[CH3][C@@H]([CH3])Cc1ccc(CC(C)C(=O)O)cc1"
        assert mcs_size(plain, stereo) >= 12


# ---------------------------------------------------------------------------
# 4. Charged species — H count vs formal charge orthogonality
# ---------------------------------------------------------------------------

class TestChargedSpecies:

    def test_ammonium_vs_ammonia_charge_ignored(self):
        # NH4+ vs NH3 — with charge matching off, MCS is the N atom
        nh3  = parse_smiles("N")
        nh4p = parse_smiles("[NH4+]")
        chem = ChemOptions()
        chem.match_formal_charge = False
        mo = MCSOptions()
        mo.timeout_ms = 5_000
        assert len(find_mcs(nh3, nh4p, chem, mo)) == 1

    def test_carboxylate_vs_carboxylic_acid(self):
        # CC(=O)[O-] vs CC(=O)O — same heavy-atom skeleton
        chem = ChemOptions()
        chem.match_formal_charge = False
        mo = MCSOptions()
        mo.timeout_ms = 5_000
        g1 = parse_smiles("CC(=O)[O-]")
        g2 = parse_smiles("CC(=O)O")
        assert len(find_mcs(g1, g2, chem, mo)) == 4

    def test_phosphate_charged_vs_neutral(self):
        # [O-]P(=O)([O-])[O-] vs OP(=O)(O)O — 4-atom P+O skeleton charge-agnostic
        chem = ChemOptions()
        chem.match_formal_charge = False
        mo = MCSOptions()
        mo.timeout_ms = 5_000
        g1 = parse_smiles("[O-]P(=O)([O-])[O-]")
        g2 = parse_smiles("OP(=O)(O)O")
        assert len(find_mcs(g1, g2, chem, mo)) >= 4
