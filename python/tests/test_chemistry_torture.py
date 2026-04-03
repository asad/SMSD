# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
"""
High-signal chemistry torture tests for representative hard classes.

These are intentionally compact and non-redundant:
- round-trip and self-MCS on representative chemistry families
- pairwise substructure/MCS on extension and charge-sensitive cases
- SMARTS motifs on biochemically relevant functionality
- one explicit expected-failure for aromatic [nH] coverage
"""

from __future__ import annotations

import pytest

import smsd


REPRESENTATIVE_MOLECULES = [
    ("glycine_zwitterion", "[NH3+]CC([O-])=O"),
    ("disconnected_salt_like", "[Na+].[O-]C(=O)C"),
    ("glutathione_like", "NCC(=O)N[C@@H](CS)C(=O)NCC(=O)O"),
    ("atp_like", "Nc1ncnc2n(cnc12)[C@@H]1O[C@H](COP(=O)(O)OP(=O)(O)O)[C@@H](O)[C@H]1O"),
    ("nad_like_fragment", "NC(=O)c1ccc[n+](c1)[C@H]1O[C@H](COP(=O)(O)O)[C@@H](O)[C@@H]1O"),
    ("glucose_like", "OC1OC(CO)C(O)C(O)C1O"),
    ("glycoside_like", "OC1OC(CO)C(O)C(O)C1OC2OC(CO)C(O)C(O)C2O"),
    ("macro_lactone_like", "CC1OC(=O)CCCCCCC(O)CO1"),
    ("steroid_like", "CC12CCC3C(C1CCC2O)CCC4=CC(=O)CCC34C"),
    ("sulfonamide_drug_like", "NS(=O)(=O)c1ccc(N)cc1"),
    ("quaternary_ammonium_like", "C[N+](C)(C)CCO"),
    ("boronate_like", "OB(O)c1ccccc1"),
    ("guanidinium_like", "NC(=[NH2+])N"),
    ("tetrazole_like", "c1nnn[nH]1"),
    ("porphyrin_like_fragment", "c1cc2cc3ccc(cc4ccc(cc1[nH]2)[nH]4)n3"),
    ("fused_natural_product_like", "COc1cc2CCC3CCCC(C3)C2cc1O"),
    ("symmetric_cage", "C1C2CC3CC1CC(C2)C3"),
]


SUBSTRUCTURE_CASES = [
    (
        "peptide_extension",
        "NCC(=O)NCC(=O)O",
        "NCC(=O)N[C@@H](CO)C(=O)NCC(=O)O",
        9,
        True,
        None,
    ),
    (
        "phosphate_chain_extension",
        "COP(=O)(O)O",
        "COP(=O)(O)OP(=O)(O)O",
        5,
        True,
        None,
    ),
    (
        "sugar_ring_extension",
        "OC1OC(CO)C(O)C(O)C1O",
        "OC1OC(CO)C(O)C(O)C1OC2OC(CO)C(O)C(O)C2O",
        12,
        True,
        None,
    ),
    (
        "sulfonamide_variant",
        "NS(=O)(=O)c1ccccc1",
        "NS(=O)(=O)c1ccc(N)cc1",
        9,
        True,
        None,
    ),
    (
        "zwitterion_vs_neutral_relaxed_charge",
        "[NH3+]CC([O-])=O",
        "NCC(=O)O",
        4,
        False,
        {"match_formal_charge": False},
    ),
    (
        "disconnected_fragment_subset",
        "c1ccccc1",
        "c1ccccc1.CC(=O)O",
        6,
        True,
        None,
    ),
]


SMARTS_CASES = [
    ("amide_in_peptide", "[NX3][CX3](=O)[#6]", "NCC(=O)N[C@@H](CO)C(=O)NCC(=O)O"),
    ("phosphate_motif", "[P](=O)(O)O", "Nc1ncnc2n(cnc12)[C@@H]1O[C@H](COP(=O)(O)OP(=O)(O)O)[C@@H](O)[C@H]1O"),
    ("sulfonamide_motif", "[NX3]S(=O)(=O)[#6]", "NS(=O)(=O)c1ccc(N)cc1"),
    ("thiol_motif", "[#16]", "NCC(=O)N[C@@H](CS)C(=O)NCC(=O)O"),
    ("carboxylate_motif", "[CX3](=O)[O-]", "[NH3+]CC([O-])=O"),
    ("quaternary_ammonium_motif", "[N+]", "C[N+](C)(C)CCO"),
    ("boronic_acid_motif", "[B]", "OB(O)c1ccccc1"),
]


class TestRepresentativeRoundTrip:
    @pytest.mark.parametrize(("name", "smiles"), REPRESENTATIVE_MOLECULES, ids=[name for name, _ in REPRESENTATIVE_MOLECULES])
    def test_roundtrip_and_self_mcs(self, name: str, smiles: str) -> None:
        mol = smsd.parse_smiles(smiles)
        assert mol.n > 0, f"{name}: parse should produce atoms"

        canonical = smsd.to_smiles(mol)
        assert canonical, f"{name}: canonical SMILES should not be empty"

        reparsed = smsd.parse_smiles(canonical)
        assert reparsed.n == mol.n, f"{name}: canonical round-trip should preserve atom count"

        mapping = smsd.find_mcs(mol, reparsed, timeout_ms=10000)
        assert len(mapping) == mol.n, f"{name}: self-equivalent round-trip should preserve full MCS"


class TestPairwiseCoverage:
    @pytest.mark.parametrize(
        ("name", "query", "target", "min_mcs", "strict_sub", "relaxed_overrides"),
        SUBSTRUCTURE_CASES,
        ids=[name for name, *_ in SUBSTRUCTURE_CASES],
    )
    def test_substructure_and_mcs(
        self,
        name: str,
        query: str,
        target: str,
        min_mcs: int,
        strict_sub: bool,
        relaxed_overrides: dict[str, object] | None,
    ) -> None:
        q = smsd.parse_smiles(query)
        t = smsd.parse_smiles(target)

        relaxed_opts = None
        if relaxed_overrides:
            relaxed_opts = smsd.ChemOptions()
            for key, value in relaxed_overrides.items():
                setattr(relaxed_opts, key, value)

        if relaxed_opts:
            mcs_opts = smsd.MCSOptions()
            mcs_opts.timeout_ms = 10000
            mcs = smsd.find_mcs(q, t, relaxed_opts, mcs_opts)
        else:
            mcs = smsd.find_mcs(q, t, timeout_ms=10000)
        assert len(mcs) >= min_mcs, f"{name}: expected MCS >= {min_mcs}, got {len(mcs)}"

        strict_hit = smsd.is_substructure(q, t)
        assert strict_hit is strict_sub, f"{name}: strict substructure expected {strict_sub}, got {strict_hit}"

        if relaxed_opts:
            relaxed_hit = smsd.is_substructure(q, t, relaxed_opts)
            assert relaxed_hit, f"{name}: relaxed substructure should succeed"


class TestSmartsCoverage:
    @pytest.mark.parametrize(("name", "smarts", "target"), SMARTS_CASES, ids=[name for name, *_ in SMARTS_CASES])
    def test_smarts_motifs(self, name: str, smarts: str, target: str) -> None:
        mol = smsd.parse_smiles(target)
        got = smsd.smarts_match(smarts, mol)
        assert got, f"{name}: SMARTS {smarts!r} should match"

    def test_aromatic_nh_smarts_matches(self) -> None:
        mol = smsd.parse_smiles("c1cc2cc3ccc(cc4ccc(cc1[nH]2)[nH]4)n3")
        assert smsd.smarts_match("[nH]", mol)
