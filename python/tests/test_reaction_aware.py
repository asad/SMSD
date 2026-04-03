# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
# Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
# See the NOTICE file for attribution, trademark, and algorithm IP terms.
"""
Tests for the reaction-aware MCS post-filter (v6.4.0).

Run with:  pytest tests/test_reaction_aware.py -v
"""
import pytest
import smsd
from smsd import (
    parse_smiles,
    MCSOptions,
    ChemOptions,
)

# Skip entire module if C++ library lacks reaction-aware support
pytestmark = pytest.mark.skipif(
    not getattr(smsd, '_HAS_REACTION_AWARE', False),
    reason="C++ library compiled without reaction-aware MCS support"
)


# SAM (S-adenosylmethionine) simplified -- includes the key S atom
SAM_SMILES = "C[S+](CCC(N)C(=O)O)CC1OC(n2cnc3c(N)ncnc32)C(O)C1O"
# SAH (S-adenosylhomocysteine) -- SAM after demethylation, S retained
# Correct structure: S bridges homocysteine to ribose via thioether
SAH_SMILES = "C(CC(N)C(=O)O)SCC1OC(n2cnc3c(N)ncnc32)C(O)C1O"


class TestReactionAwareMCS:
    """Reaction-aware MCS post-filter tests."""

    def test_sam_to_sah_includes_sulfur(self):
        """SAM -> SAH: reaction-aware MCS should include the S atom."""
        g1 = parse_smiles(SAM_SMILES)
        g2 = parse_smiles(SAH_SMILES)

        mapping = smsd.map_reaction_aware(g1, g2, timeout_ms=10000)
        assert mapping, "Reaction-aware MCS should produce a non-empty mapping"

        # Check that at least one mapped atom in g1 is sulfur (atomic number 16)
        has_sulfur = any(g1.atomic_num[qi] == 16 for qi in mapping)
        assert has_sulfur, (
            "Reaction-aware MCS for SAM->SAH should include the sulfur atom, "
            f"but mapped atoms have atomic numbers: "
            f"{[g1.atomic_num[qi] for qi in mapping]}"
        )

    def test_sam_to_sah_from_smiles(self):
        """map_reaction_aware works with SMILES strings."""
        mapping = smsd.map_reaction_aware(SAM_SMILES, SAH_SMILES, timeout_ms=10000)
        assert mapping, "Should produce a mapping from SMILES strings"

    def test_pure_hydrocarbon_falls_back(self):
        """Pure hydrocarbons: reaction-aware degrades to standard MCS size."""
        g1 = parse_smiles("C1CCCCC1")      # cyclohexane
        g2 = parse_smiles("CC1CCCC1")       # methylcyclopentane

        std_mapping = smsd.find_mcs(g1, g2)
        react_mapping = smsd.map_reaction_aware(g1, g2, timeout_ms=5000)

        # For pure hydrocarbons, scoring degrades to size -- should be same size
        assert len(react_mapping) == len(std_mapping), (
            f"Pure hydrocarbon: reaction-aware ({len(react_mapping)}) should "
            f"match standard MCS size ({len(std_mapping)})"
        )

    def test_mcs_options_reaction_aware_field(self):
        """MCSOptions has reaction_aware field, defaults to False."""
        opts = MCSOptions()
        assert opts.reaction_aware is False
        assert opts.near_mcs_delta == 2
        assert opts.near_mcs_candidates == 20

    def test_methionine_to_homocysteine(self):
        """Methionine -> homocysteine: S atom should be captured."""
        met = parse_smiles("CSCC(N)C(=O)O")   # methionine
        hcy = parse_smiles("SCC(N)C(=O)O")     # homocysteine

        mapping = smsd.map_reaction_aware(met, hcy, timeout_ms=5000)
        assert mapping, "Should find a mapping"

        has_sulfur = any(met.atomic_num[qi] == 16 for qi in mapping)
        assert has_sulfur, "Methionine->homocysteine mapping should include S"

    def test_atp_to_adp_includes_phosphorus(self):
        """ATP -> ADP: reaction-aware MCS should include a P atom."""
        atp = parse_smiles(
            "c1nc(N)c2ncn(C3OC(COP(=O)(O)OP(=O)(O)OP(=O)(O)O)C(O)C3O)c2n1")
        adp = parse_smiles(
            "c1nc(N)c2ncn(C3OC(COP(=O)(O)OP(=O)(O)O)C(O)C3O)c2n1")

        mapping = smsd.map_reaction_aware(atp, adp, timeout_ms=10000)
        assert mapping, "Reaction-aware MCS for ATP->ADP should be non-empty"

        has_phosphorus = any(atp.atomic_num[qi] == 15 for qi in mapping)
        assert has_phosphorus, (
            "Reaction-aware MCS for ATP->ADP should include a phosphorus atom, "
            f"but mapped atoms have atomic numbers: "
            f"{[atp.atomic_num[qi] for qi in mapping]}"
        )
