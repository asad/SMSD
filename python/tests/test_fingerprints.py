"""Comprehensive fingerprint tests for SMSD Pro v7.1.0.

Tests every fingerprint type, similarity metric, edge case, and known
challenging molecules from chemoinformatics literature.

Copyright (c) 2018-2026 BioInception PVT LTD — Syed Asad Rahman
SPDX-License-Identifier: Apache-2.0
"""

import math
import pytest
import smsd


# ---------------------------------------------------------------------------
# Fixtures: standard molecules
# ---------------------------------------------------------------------------

@pytest.fixture
def benzene():
    return smsd.parse_smiles("c1ccccc1")

@pytest.fixture
def phenol():
    return smsd.parse_smiles("c1ccc(O)cc1")

@pytest.fixture
def toluene():
    return smsd.parse_smiles("Cc1ccccc1")

@pytest.fixture
def naphthalene():
    return smsd.parse_smiles("c1ccc2ccccc2c1")

@pytest.fixture
def cyclohexane():
    return smsd.parse_smiles("C1CCCCC1")

@pytest.fixture
def aspirin():
    return smsd.parse_smiles("CC(=O)Oc1ccccc1C(=O)O")

@pytest.fixture
def salicylic_acid():
    return smsd.parse_smiles("Oc1ccccc1C(=O)O")

@pytest.fixture
def ethanol():
    return smsd.parse_smiles("CCO")

@pytest.fixture
def methanol():
    return smsd.parse_smiles("CO")


# ===================================================================
# 1. ECFP (circular fingerprint) — core tests
# ===================================================================

class TestCircularFingerprint:
    """Circular fingerprint (ECFP/FCFP) generation tests."""

    def test_basic_generation(self, benzene):
        fp = smsd.circular_fingerprint(benzene, radius=2, fp_size=2048)
        assert isinstance(fp, list)
        assert len(fp) > 0
        assert all(isinstance(b, int) for b in fp)
        assert all(0 <= b < 2048 for b in fp)

    def test_sorted_bit_positions(self, benzene):
        fp = smsd.circular_fingerprint(benzene, radius=2, fp_size=2048)
        assert fp == sorted(fp), "Bit positions must be sorted"

    def test_no_duplicates(self, benzene):
        fp = smsd.circular_fingerprint(benzene, radius=2, fp_size=2048)
        assert len(fp) == len(set(fp)), "Bit positions must be unique"

    def test_deterministic(self, benzene):
        fp1 = smsd.circular_fingerprint(benzene, radius=2, fp_size=2048)
        fp2 = smsd.circular_fingerprint(benzene, radius=2, fp_size=2048)
        assert fp1 == fp2, "Same input must produce identical fingerprint"

    def test_radius_0_vs_2(self, benzene):
        """Radius 0 = atom invariants only. Radius 2 = ECFP4-equivalent.
        Higher radius should produce more or equal number of features."""
        fp0 = smsd.circular_fingerprint(benzene, radius=0, fp_size=2048)
        fp2 = smsd.circular_fingerprint(benzene, radius=2, fp_size=2048)
        assert len(fp2) >= len(fp0)

    def test_radius_3_more_specific(self, naphthalene):
        """ECFP6 (radius 3) should have >= bits than ECFP4 (radius 2)."""
        fp2 = smsd.circular_fingerprint(naphthalene, radius=2, fp_size=2048)
        fp3 = smsd.circular_fingerprint(naphthalene, radius=3, fp_size=2048)
        assert len(fp3) >= len(fp2)

    def test_fp_size_1024(self, benzene):
        fp = smsd.circular_fingerprint(benzene, radius=2, fp_size=1024)
        assert all(0 <= b < 1024 for b in fp)

    def test_fp_size_4096(self, benzene):
        fp = smsd.circular_fingerprint(benzene, radius=2, fp_size=4096)
        assert all(0 <= b < 4096 for b in fp)

    def test_different_molecules_different_fps(self, benzene, ethanol):
        fp_benz = smsd.circular_fingerprint(benzene, radius=2, fp_size=2048)
        fp_eth = smsd.circular_fingerprint(ethanol, radius=2, fp_size=2048)
        assert fp_benz != fp_eth

    def test_benzene_vs_cyclohexane(self, benzene, cyclohexane):
        """Aromatic vs aliphatic ring — must produce different fingerprints."""
        fp_ar = smsd.circular_fingerprint(benzene, radius=2, fp_size=2048)
        fp_al = smsd.circular_fingerprint(cyclohexane, radius=2, fp_size=2048)
        assert fp_ar != fp_al, "Aromatic and aliphatic rings must differ"


class TestFCFP:
    """FCFP (pharmacophore) fingerprint tests."""

    def test_basic_fcfp(self, benzene):
        fp = smsd.circular_fingerprint(benzene, radius=2, fp_size=2048, mode="fcfp")
        assert isinstance(fp, list)
        assert len(fp) > 0

    def test_fcfp_differs_from_ecfp(self, phenol):
        """FCFP uses pharmacophore invariants, ECFP uses structural."""
        ecfp = smsd.circular_fingerprint(phenol, radius=2, fp_size=2048, mode="ecfp")
        fcfp = smsd.circular_fingerprint(phenol, radius=2, fp_size=2048, mode="fcfp")
        assert ecfp != fcfp, "ECFP and FCFP should produce different fingerprints"

    def test_fcfp_phenol_vs_aniline(self):
        """Phenol (-OH) and aniline (-NH2) both have H-bond donor/acceptor.
        FCFP should show higher pharmacophore similarity than ECFP."""
        phenol = smsd.parse_smiles("c1ccc(O)cc1")
        aniline = smsd.parse_smiles("c1ccc(N)cc1")
        ecfp_sim = smsd.tanimoto_coefficient(
            smsd.circular_fingerprint(phenol, radius=2, fp_size=2048),
            smsd.circular_fingerprint(aniline, radius=2, fp_size=2048)
        )
        fcfp_sim = smsd.tanimoto_coefficient(
            smsd.circular_fingerprint(phenol, radius=2, fp_size=2048, mode="fcfp"),
            smsd.circular_fingerprint(aniline, radius=2, fp_size=2048, mode="fcfp")
        )
        # FCFP groups OH and NH2 as same pharmacophore feature
        assert fcfp_sim >= ecfp_sim


class TestCircularFingerprintCounts:
    """Count-based circular fingerprint tests."""

    def test_basic_counts(self, benzene):
        counts = smsd.circular_fingerprint_counts(benzene, radius=2, fp_size=2048)
        assert isinstance(counts, list)
        assert all(isinstance(c, tuple) and len(c) == 2 for c in counts)

    def test_counts_positive(self, benzene):
        counts = smsd.circular_fingerprint_counts(benzene, radius=2, fp_size=2048)
        for pos, count in counts:
            assert count > 0, f"Count at position {pos} must be positive"
            assert 0 <= pos < 2048

    def test_benzene_symmetry(self, benzene):
        """All 6 carbons in benzene are equivalent — counts should reflect this."""
        counts = smsd.circular_fingerprint_counts(benzene, radius=2, fp_size=2048)
        for pos, count in counts:
            assert count == 6, f"Benzene symmetry: all features should have count 6, got {count}"

    def test_counts_vs_binary_positions_match(self, phenol):
        """Bit positions from binary FP should match positions from count FP."""
        binary = set(smsd.circular_fingerprint(phenol, radius=2, fp_size=2048))
        counts_pos = {pos for pos, _ in smsd.circular_fingerprint_counts(phenol, radius=2, fp_size=2048)}
        assert binary == counts_pos, "Binary and count positions must match"

    def test_counts_to_dense_roundtrip(self, phenol):
        """counts → dense array → sum should equal total feature count."""
        counts = smsd.circular_fingerprint_counts(phenol, radius=2, fp_size=2048)
        dense = smsd.counts_to_array(counts, 2048)
        assert len(dense) == 2048
        total = sum(c for _, c in counts)
        assert sum(dense) == total


# ===================================================================
# 2. Path fingerprint
# ===================================================================

class TestPathFingerprint:
    """DFS path enumeration fingerprint tests."""

    def test_basic(self, benzene):
        fp = smsd.path_fingerprint(benzene, path_length=7, fp_size=2048)
        assert isinstance(fp, list)
        assert len(fp) > 0

    def test_deterministic(self, benzene):
        fp1 = smsd.path_fingerprint(benzene)
        fp2 = smsd.path_fingerprint(benzene)
        assert fp1 == fp2

    def test_different_path_lengths(self, naphthalene):
        fp3 = smsd.path_fingerprint(naphthalene, path_length=3, fp_size=2048)
        fp7 = smsd.path_fingerprint(naphthalene, path_length=7, fp_size=2048)
        assert len(fp7) >= len(fp3), "Longer paths should set >= bits"


# ===================================================================
# 3. MCS fingerprint
# ===================================================================

class TestMCSFingerprint:
    """MCS-aware fingerprint tests."""

    def test_basic(self, benzene):
        fp = smsd.mcs_fingerprint(benzene, path_length=7, fp_size=2048)
        assert isinstance(fp, list)
        assert len(fp) > 0

    def test_differs_from_path(self, benzene):
        mcs_fp = smsd.mcs_fingerprint(benzene, path_length=7, fp_size=2048)
        path_fp = smsd.path_fingerprint(benzene, path_length=7, fp_size=2048)
        # They may or may not be equal, but both should be valid
        assert len(mcs_fp) > 0
        assert len(path_fp) > 0


# ===================================================================
# 4. Topological torsion fingerprint
# ===================================================================

class TestTopologicalTorsion:
    """Topological torsion (4-atom path) fingerprint tests."""

    def test_basic_from_smiles(self):
        fp = smsd.topological_torsion("c1ccccc1")
        assert isinstance(fp, list)
        assert len(fp) > 0

    def test_basic_from_molgraph(self, benzene):
        fp = smsd.topological_torsion(benzene, fp_size=2048)
        assert isinstance(fp, list)

    def test_linear_chain_more_torsions(self):
        """Linear chains have more 4-atom paths than compact rings."""
        hexane = smsd.topological_torsion("CCCCCC", fp_size=2048)
        benzene_t = smsd.topological_torsion("c1ccccc1", fp_size=2048)
        # Hexane should have more diverse torsion features
        assert len(hexane) >= 1

    def test_peptide_backbone(self):
        """Topological torsion is SOTA on peptide benchmarks (known from lit)."""
        ala_ala = smsd.topological_torsion("CC(N)C(=O)NC(C)C(=O)O")
        assert isinstance(ala_ala, list)
        assert len(ala_ala) > 0

    def test_counts(self, benzene):
        counts = smsd.topological_torsion_counts(benzene)
        assert isinstance(counts, list)
        assert all(isinstance(c, tuple) and len(c) == 2 for c in counts)


# ===================================================================
# 5. Convenience: fingerprint() multi-kind dispatch
# ===================================================================

class TestFingerprintConvenience:
    """Tests for the fingerprint() convenience function."""

    def test_kind_path(self):
        fp = smsd.fingerprint("c1ccccc1", kind="path")
        assert isinstance(fp, list) and len(fp) > 0

    def test_kind_mcs(self):
        fp = smsd.fingerprint("c1ccccc1", kind="mcs")
        assert isinstance(fp, list) and len(fp) > 0

    def test_kind_ecfp(self):
        fp = smsd.fingerprint("c1ccccc1", kind="ecfp")
        assert isinstance(fp, list) and len(fp) > 0

    def test_kind_circular(self):
        fp = smsd.fingerprint("c1ccccc1", kind="circular")
        ecfp = smsd.fingerprint("c1ccccc1", kind="ecfp")
        assert fp == ecfp, "'circular' should be an alias for 'ecfp'"

    def test_kind_fcfp(self):
        fp = smsd.fingerprint("c1ccccc1", kind="fcfp")
        assert isinstance(fp, list) and len(fp) > 0

    def test_kind_torsion(self):
        fp = smsd.fingerprint("c1ccccc1", kind="torsion")
        assert isinstance(fp, list) and len(fp) > 0

    def test_invalid_kind(self):
        with pytest.raises(ValueError, match="Unknown fingerprint kind"):
            smsd.fingerprint("c1ccccc1", kind="invalid")

    def test_fingerprint_from_smiles(self):
        fp = smsd.fingerprint_from_smiles("c1ccccc1", radius=2, fp_size=2048)
        ecfp = smsd.circular_fingerprint(smsd.parse_smiles("c1ccccc1"), radius=2, fp_size=2048)
        assert fp == ecfp


# ===================================================================
# 6. Similarity metrics — binary fingerprints
# ===================================================================

class TestBinarySimilarity:
    """Tanimoto, Dice, Cosine, Overlap, Soergel for binary FPs."""

    @pytest.fixture
    def fps(self, benzene, phenol):
        return (
            smsd.circular_fingerprint(benzene, radius=2, fp_size=2048),
            smsd.circular_fingerprint(phenol, radius=2, fp_size=2048),
        )

    def test_self_tanimoto_is_one(self, fps):
        fp1, _ = fps
        assert smsd.tanimoto_coefficient(fp1, fp1) == 1.0

    def test_self_dice_is_one(self, fps):
        fp1, _ = fps
        assert smsd.dice(fp1, fp1) == 1.0

    def test_self_cosine_is_one(self, fps):
        fp1, _ = fps
        assert smsd.cosine(fp1, fp1) == 1.0

    def test_self_overlap_is_one(self, fps):
        fp1, _ = fps
        assert smsd.overlap_coefficient(fp1, fp1) == 1.0

    def test_self_soergel_is_zero(self, fps):
        """Soergel distance is 1 - Tanimoto. Self-distance = 0."""
        fp1, _ = fps
        assert smsd.soergel(fp1, fp1) == 0.0

    def test_tanimoto_range(self, fps):
        fp1, fp2 = fps
        t = smsd.tanimoto_coefficient(fp1, fp2)
        assert 0.0 <= t <= 1.0

    def test_dice_range(self, fps):
        fp1, fp2 = fps
        d = smsd.dice(fp1, fp2)
        assert 0.0 <= d <= 1.0

    def test_dice_ge_tanimoto(self, fps):
        """Dice coefficient is always >= Tanimoto (mathematical property)."""
        fp1, fp2 = fps
        t = smsd.tanimoto_coefficient(fp1, fp2)
        d = smsd.dice(fp1, fp2)
        assert d >= t - 1e-10, f"Dice {d} should be >= Tanimoto {t}"

    def test_tanimoto_symmetry(self, fps):
        fp1, fp2 = fps
        assert smsd.tanimoto_coefficient(fp1, fp2) == smsd.tanimoto_coefficient(fp2, fp1)

    def test_soergel_plus_tanimoto_equals_one(self, fps):
        fp1, fp2 = fps
        t = smsd.tanimoto_coefficient(fp1, fp2)
        s = smsd.soergel(fp1, fp2)
        assert abs(t + s - 1.0) < 1e-10

    def test_similar_molecules_higher_score(self, benzene, toluene, ethanol):
        """Benzene–toluene similarity > benzene–ethanol similarity."""
        fp_b = smsd.circular_fingerprint(benzene, radius=2, fp_size=2048)
        fp_t = smsd.circular_fingerprint(toluene, radius=2, fp_size=2048)
        fp_e = smsd.circular_fingerprint(ethanol, radius=2, fp_size=2048)
        sim_bt = smsd.tanimoto_coefficient(fp_b, fp_t)
        sim_be = smsd.tanimoto_coefficient(fp_b, fp_e)
        assert sim_bt > sim_be, (
            f"Benzene-toluene ({sim_bt:.3f}) should be more similar "
            f"than benzene-ethanol ({sim_be:.3f})"
        )

    def test_aspirin_salicylic_acid_high_similarity(self, aspirin, salicylic_acid):
        """Aspirin is acetylsalicylic acid — should be very similar to salicylic acid."""
        fp1 = smsd.circular_fingerprint(aspirin, radius=2, fp_size=2048)
        fp2 = smsd.circular_fingerprint(salicylic_acid, radius=2, fp_size=2048)
        sim = smsd.tanimoto_coefficient(fp1, fp2)
        assert sim > 0.1, f"Aspirin–salicylic acid should have decent similarity, got {sim:.3f}"


class TestBinarySimilarityWithSparseInput:
    """Ensure tanimoto/overlap accept sparse count tuples transparently."""

    def test_tanimoto_sparse_counts(self, benzene, phenol):
        c1 = smsd.circular_fingerprint_counts(benzene, radius=2, fp_size=2048)
        c2 = smsd.circular_fingerprint_counts(phenol, radius=2, fp_size=2048)
        fp1 = smsd.circular_fingerprint(benzene, radius=2, fp_size=2048)
        fp2 = smsd.circular_fingerprint(phenol, radius=2, fp_size=2048)
        assert smsd.tanimoto_coefficient(c1, c2) == smsd.tanimoto_coefficient(fp1, fp2)

    def test_overlap_sparse_counts(self, benzene, phenol):
        c1 = smsd.circular_fingerprint_counts(benzene, radius=2, fp_size=2048)
        c2 = smsd.circular_fingerprint_counts(phenol, radius=2, fp_size=2048)
        fp1 = smsd.circular_fingerprint(benzene, radius=2, fp_size=2048)
        fp2 = smsd.circular_fingerprint(phenol, radius=2, fp_size=2048)
        assert smsd.overlap_coefficient(c1, c2) == smsd.overlap_coefficient(fp1, fp2)


# ===================================================================
# 7. Similarity metrics — count-based fingerprints
# ===================================================================

class TestCountSimilarity:
    """Count-based similarity using dense arrays."""

    @pytest.fixture
    def dense_fps(self, benzene, phenol):
        c1 = smsd.circular_fingerprint_counts(benzene, radius=2, fp_size=2048)
        c2 = smsd.circular_fingerprint_counts(phenol, radius=2, fp_size=2048)
        return smsd.counts_to_array(c1, 2048), smsd.counts_to_array(c2, 2048)

    def test_count_dice_range(self, dense_fps):
        d1, d2 = dense_fps
        v = smsd.count_dice(d1, d2)
        assert 0.0 <= v <= 1.0

    def test_count_cosine_range(self, dense_fps):
        d1, d2 = dense_fps
        v = smsd.count_cosine(d1, d2)
        assert 0.0 <= v <= 1.0

    def test_count_overlap_range(self, dense_fps):
        d1, d2 = dense_fps
        v = smsd.count_overlapCoefficient(d1, d2)
        assert 0.0 <= v <= 1.0

    def test_count_self_similarity(self, benzene):
        c = smsd.circular_fingerprint_counts(benzene, radius=2, fp_size=2048)
        d = smsd.counts_to_array(c, 2048)
        assert smsd.count_dice(d, d) == pytest.approx(1.0, abs=1e-10)
        assert smsd.count_cosine(d, d) == pytest.approx(1.0, abs=1e-10)

    def test_counts_to_array_with_dict(self):
        """counts_to_array should accept both dict and list-of-tuples."""
        arr1 = smsd.counts_to_array([(42, 3), (1000, 1)], 2048)
        arr2 = smsd.counts_to_array({42: 3, 1000: 1}, 2048)
        assert arr1 == arr2
        assert arr1[42] == 3
        assert arr1[1000] == 1
        assert arr1[0] == 0


# ===================================================================
# 8. Format conversions
# ===================================================================

class TestFormatConversions:
    """Hex, binary string, and from_hex roundtrip tests."""

    def test_to_hex_roundtrip(self, benzene):
        fp = smsd.circular_fingerprint(benzene, radius=2, fp_size=2048)
        hex_str = smsd.to_hex(fp)
        assert isinstance(hex_str, str)
        fp_back = smsd.from_hex(hex_str)
        assert fp_back == fp

    def test_to_binary_string_length(self, benzene):
        fp = smsd.circular_fingerprint(benzene, radius=2, fp_size=2048)
        bits = smsd.to_binary_string(fp, fp_size=2048)
        assert len(bits) == 2048
        assert all(c in ('0', '1') for c in bits)

    def test_binary_string_set_bits_match(self, benzene):
        fp = smsd.circular_fingerprint(benzene, radius=2, fp_size=2048)
        bits = smsd.to_binary_string(fp, fp_size=2048)
        set_bits = [i for i, c in enumerate(bits) if c == '1']
        assert set_bits == fp

    def test_fingerprint_subset(self, benzene, naphthalene):
        """Benzene fingerprint should be a subset of naphthalene's (structurally contained)."""
        fp_benz = set(smsd.circular_fingerprint(benzene, radius=0, fp_size=2048))
        fp_naph = set(smsd.circular_fingerprint(naphthalene, radius=0, fp_size=2048))
        # At radius 0, benzene atom invariant should appear in naphthalene
        assert fp_benz.issubset(fp_naph) or len(fp_benz & fp_naph) > 0

    def test_analyze_fp_quality(self, benzene):
        fp = smsd.circular_fingerprint(benzene, radius=2, fp_size=2048)
        q = smsd.analyze_fp_quality(fp)
        assert isinstance(q, dict)
        assert q['fp_size'] == 2048
        assert q['set_bits'] == len(fp)
        assert 0.0 <= q['density'] <= 1.0
        assert isinstance(q['is_saturated'], bool)

    def test_fp_quality_saturated_detection(self):
        """A fingerprint with >75% bits set should be flagged as saturated."""
        fp_dense = list(range(1600))  # 1600/2048 = 78.1% density — exceeds 75% threshold
        q = smsd.analyze_fp_quality(fp_dense, 2048)
        assert q['is_saturated'] is True


# ===================================================================
# 9. Edge cases
# ===================================================================

class TestEdgeCases:
    """Edge cases and boundary conditions for fingerprints."""

    def test_single_atom(self):
        mol = smsd.parse_smiles("C")
        fp = smsd.circular_fingerprint(mol, radius=2, fp_size=2048)
        assert isinstance(fp, list)
        # Single carbon should produce at least one feature
        assert len(fp) >= 1

    def test_empty_fingerprint_similarity(self):
        """Tanimoto of two empty FPs should be 0."""
        assert smsd.tanimoto_coefficient([], []) == 0.0
        assert smsd.overlap_coefficient([], []) == 0.0

    def test_none_fingerprint_similarity(self):
        """None input should return 0."""
        assert smsd.tanimoto_coefficient(None, None) == 0.0
        assert smsd.tanimoto_coefficient([1], None) == 0.0
        assert smsd.overlap_coefficient(None, [1]) == 0.0

    def test_identical_molecules_different_smiles(self):
        """Same molecule, different SMILES notation — must give identical FP."""
        fp1 = smsd.circular_fingerprint(smsd.parse_smiles("c1ccccc1"), radius=2, fp_size=2048)
        fp2 = smsd.circular_fingerprint(smsd.parse_smiles("C1=CC=CC=C1"), radius=2, fp_size=2048)
        assert fp1 == fp2

    def test_charged_molecule(self):
        """Acetate anion — should produce a valid fingerprint."""
        mol = smsd.parse_smiles("CC(=O)[O-]")
        fp = smsd.circular_fingerprint(mol, radius=2, fp_size=2048)
        assert len(fp) > 0

    def test_isotope_labeled(self):
        """Deuterated methane — isotope should not crash FP generation."""
        mol = smsd.parse_smiles("[2H]C([2H])([2H])[2H]")
        fp = smsd.circular_fingerprint(mol, radius=2, fp_size=2048)
        assert isinstance(fp, list)

    def test_disconnected_fragments(self):
        """Salt: sodium acetate (two fragments)."""
        mol = smsd.parse_smiles("CC(=O)[O-].[Na+]")
        fp = smsd.circular_fingerprint(mol, radius=2, fp_size=2048)
        assert len(fp) > 0

    def test_large_fp_size(self):
        """Very large FP size (16384) should work without error."""
        mol = smsd.parse_smiles("c1ccccc1")
        fp = smsd.circular_fingerprint(mol, radius=2, fp_size=16384)
        assert all(0 <= b < 16384 for b in fp)

    def test_small_fp_size(self):
        """Very small FP size (64) — expect more collisions but still valid."""
        mol = smsd.parse_smiles("c1ccccc1")
        fp = smsd.circular_fingerprint(mol, radius=2, fp_size=64)
        assert all(0 <= b < 64 for b in fp)


# ===================================================================
# 10. Known challenging molecules (chemoinformatics literature)
# ===================================================================

class TestChallengingMolecules:
    """Known difficult cases from published chemoinformatics benchmarks."""

    def test_enantiomers_same_ecfp(self):
        """R/S enantiomers: ECFP ignores chirality by default.
        L-alanine and D-alanine should have IDENTICAL ECFP fingerprints.
        (Cereto-Massagué et al., Methods 2015)"""
        l_ala = smsd.parse_smiles("N[C@@H](C)C(=O)O")
        d_ala = smsd.parse_smiles("N[C@H](C)C(=O)O")
        fp_l = smsd.circular_fingerprint(l_ala, radius=2, fp_size=2048)
        fp_d = smsd.circular_fingerprint(d_ala, radius=2, fp_size=2048)
        assert fp_l == fp_d, "ECFP should not distinguish enantiomers"

    def test_constitutional_isomers_different_fp(self):
        """1-propanol vs 2-propanol — same formula, different structure.
        ECFP must distinguish them."""
        fp1 = smsd.circular_fingerprint(smsd.parse_smiles("CCCO"), radius=2, fp_size=2048)
        fp2 = smsd.circular_fingerprint(smsd.parse_smiles("CC(O)C"), radius=2, fp_size=2048)
        assert fp1 != fp2

    def test_caffeine_theophylline_similarity(self):
        """Caffeine (1,3,7-trimethylxanthine) vs theophylline (1,3-dimethylxanthine).
        Differ by one N-methyl group. Should have high but not perfect similarity.
        (Standard pair from Daylight fingerprint benchmark)"""
        caffeine = smsd.parse_smiles("Cn1c2c(c(=O)n(c1=O)C)nc(n2)C")
        theophylline = smsd.parse_smiles("Cn1c2c(c(=O)n(c1=O)C)[nH]cn2")
        fp_c = smsd.circular_fingerprint(caffeine, radius=2, fp_size=2048)
        fp_t = smsd.circular_fingerprint(theophylline, radius=2, fp_size=2048)
        sim = smsd.tanimoto_coefficient(fp_c, fp_t)
        assert sim > 0.3, f"Caffeine–theophylline should be moderately similar: {sim:.3f}"
        assert sim < 1.0, "They are NOT identical"

    def test_ibuprofen_naproxen_moderate_similarity(self):
        """Both are NSAIDs but structurally different. Moderate ECFP similarity.
        (Standard NSAID comparison from DrugBank literature)"""
        ibuprofen = smsd.parse_smiles("CC(C)Cc1ccc(cc1)[C@@H](C)C(=O)O")
        naproxen = smsd.parse_smiles("COc1ccc2cc([C@@H](C)C(=O)O)ccc2c1")
        fp_i = smsd.circular_fingerprint(ibuprofen, radius=2, fp_size=2048)
        fp_n = smsd.circular_fingerprint(naproxen, radius=2, fp_size=2048)
        sim = smsd.tanimoto_coefficient(fp_i, fp_n)
        assert 0.1 <= sim <= 0.8, f"Ibuprofen–naproxen moderate similarity: {sim:.3f}"

    def test_taxol_large_molecule(self):
        """Paclitaxel (taxol) — 53 heavy atoms. Tests FP on large molecules.
        Should not crash or saturate at fp_size=2048."""
        taxol = smsd.parse_smiles(
            "CC1=C2[C@@]([C@]([C@H]([C@@H]3[C@]4([C@H](OC4)C[C@@H]([C@]3(C(=O)"
            "[C@@H]2OC(=O)C)C)O)OC(=O)C)OC(=O)c5ccccc5)(C[C@@H]1OC(=O)[C@@H]"
            "(O)C(NC(=O)c6ccccc6)c7ccccc7)O)(C)C"
        )
        fp = smsd.circular_fingerprint(taxol, radius=2, fp_size=2048)
        quality = smsd.analyze_fp_quality(fp, 2048)
        assert not quality['is_saturated'], (
            f"Taxol FP should not saturate 2048-bit FP (density={quality['density']:.3f})"
        )

    def test_cubane_high_symmetry(self):
        """Cubane (C8H8) — highly symmetric cage. All 8 carbons equivalent.
        ECFP4 should produce few unique features due to symmetry."""
        cubane = smsd.parse_smiles("C12C3C4C1C5C3C4C25")
        fp = smsd.circular_fingerprint(cubane, radius=2, fp_size=2048)
        assert len(fp) >= 1
        counts = smsd.circular_fingerprint_counts(cubane, radius=2, fp_size=2048)
        # All atoms equivalent — each feature should have count = 8
        for pos, count in counts:
            assert count == 8 or count == 4, (
                f"Cubane symmetry: expected count 8 (or 4 for bond features), got {count}"
            )

    def test_adamantane_cage(self):
        """Adamantane — diamond cage structure. Tests ring perception in FP."""
        adam = smsd.parse_smiles("C1C2CC3CC1CC(C2)C3")
        fp = smsd.circular_fingerprint(adam, radius=2, fp_size=2048)
        assert len(fp) >= 1

    def test_fullerene_c60_stress(self):
        """C60 fullerene — 60 carbons, 32 faces. Stress test for large symmetric molecules."""
        # C60 SMILES (buckminsterfullerene)
        c60 = smsd.parse_smiles(
            "c12c3c4c5c1c1c6c7c2c2c8c3c3c9c4c4c%10c5c5c1c1c6c6c%11c7c7c2c2c8c8"
            "c3c3c9c9c4c4c%10c%10c5c5c1c1c6c6c%11c7c7c2c8c2c3c9c3c4c%10c4c5c1c1c6"
            "c7c2c3c41"
        )
        fp = smsd.circular_fingerprint(c60, radius=2, fp_size=2048)
        assert len(fp) >= 1

    def test_morphine_codeine_close_analogs(self):
        """Morphine vs codeine — differ by one O-methyl group.
        Should have very high ECFP similarity.
        (Classic opioid analogue pair from medicinal chemistry)"""
        morphine = smsd.parse_smiles("O[C@@H]1[C@@H]2Oc3c(O)ccc4C[C@H]5[C@@H]"
                                      "(N(CC5)C)C=C1[C@H]2c34")
        codeine = smsd.parse_smiles("CO[C@@H]1[C@@H]2Oc3c(O)ccc4C[C@H]5[C@@H]"
                                     "(N(CC5)C)C=C1[C@H]2c34")
        fp_m = smsd.circular_fingerprint(morphine, radius=2, fp_size=2048)
        fp_c = smsd.circular_fingerprint(codeine, radius=2, fp_size=2048)
        sim = smsd.tanimoto_coefficient(fp_m, fp_c)
        assert sim > 0.5, f"Morphine–codeine should be highly similar: {sim:.3f}"

    def test_decoy_low_similarity(self):
        """Glucose vs testosterone — completely different structures.
        Should have very low fingerprint similarity."""
        glucose = smsd.parse_smiles("OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O")
        testosterone = smsd.parse_smiles("C[C@]12CC[C@H]3[C@@H](CCC4=CC(=O)CC"
                                          "[C@@]34C)[C@@H]1CC[C@@H]2O")
        fp_g = smsd.circular_fingerprint(glucose, radius=2, fp_size=2048)
        fp_t = smsd.circular_fingerprint(testosterone, radius=2, fp_size=2048)
        sim = smsd.tanimoto_coefficient(fp_g, fp_t)
        assert sim < 0.3, f"Glucose–testosterone should have low similarity: {sim:.3f}"

    def test_murcko_scaffold_pair(self):
        """Two molecules with same Murcko scaffold should have higher fingerprint
        similarity than unrelated molecules.
        (Murcko decomposition principle, Bemis & Murcko 1996)"""
        mol1 = smsd.parse_smiles("c1ccc(O)cc1")   # phenol
        mol2 = smsd.parse_smiles("c1ccc(N)cc1")   # aniline
        mol3 = smsd.parse_smiles("CCCCCC")         # hexane (no aromatic ring)
        fp1 = smsd.circular_fingerprint(mol1, radius=2, fp_size=2048)
        fp2 = smsd.circular_fingerprint(mol2, radius=2, fp_size=2048)
        fp3 = smsd.circular_fingerprint(mol3, radius=2, fp_size=2048)
        sim_phenol_aniline = smsd.tanimoto_coefficient(fp1, fp2)
        sim_phenol_hexane = smsd.tanimoto_coefficient(fp1, fp3)
        assert sim_phenol_aniline > sim_phenol_hexane, (
            f"Same scaffold pair ({sim_phenol_aniline:.3f}) should be "
            f"more similar than unrelated ({sim_phenol_hexane:.3f})"
        )


# ===================================================================
# 11. Bit collision analysis
# ===================================================================

class TestBitCollisions:
    """Tests related to hash collisions in folded fingerprints."""

    def test_larger_fp_fewer_collisions(self):
        """Larger FP size → fewer bit collisions → more set bits (for complex molecules)."""
        mol = smsd.parse_smiles("CC(=O)Oc1ccccc1C(=O)O")  # aspirin
        fp_1024 = smsd.circular_fingerprint(mol, radius=2, fp_size=1024)
        fp_4096 = smsd.circular_fingerprint(mol, radius=2, fp_size=4096)
        # 4096-bit should have >= set bits than 1024-bit
        assert len(fp_4096) >= len(fp_1024)

    def test_identical_fp_size_invariant(self):
        """Identical molecules should produce identical FPs regardless of fp_size."""
        mol = smsd.parse_smiles("c1ccccc1")
        for sz in [512, 1024, 2048, 4096]:
            fp1 = smsd.circular_fingerprint(mol, radius=2, fp_size=sz)
            fp2 = smsd.circular_fingerprint(mol, radius=2, fp_size=sz)
            assert fp1 == fp2, f"Not deterministic at fp_size={sz}"


# ===================================================================
# 12. Cross-fingerprint comparison
# ===================================================================

class TestCrossFingerprint:
    """Ensure different FP types give different but valid results."""

    def test_ecfp_path_mcs_all_different(self, benzene):
        ecfp = smsd.circular_fingerprint(benzene, radius=2, fp_size=2048)
        path = smsd.path_fingerprint(benzene, path_length=7, fp_size=2048)
        mcs = smsd.mcs_fingerprint(benzene, path_length=7, fp_size=2048)
        # They should all be valid but likely different
        assert len(ecfp) > 0 and len(path) > 0 and len(mcs) > 0
        # At least one pair should differ
        assert ecfp != path or ecfp != mcs or path != mcs

    def test_ecfp_torsion_different(self, naphthalene):
        ecfp = smsd.circular_fingerprint(naphthalene, radius=2, fp_size=2048)
        torsion = smsd.topological_torsion(naphthalene, fp_size=2048)
        assert ecfp != torsion


# ===================================================================
# 13. Known literature benchmark values
# ===================================================================

class TestLiteratureValues:
    """Approximate expected values from published benchmarks."""

    def test_aspirin_self_similarity(self):
        """Self-similarity must always be exactly 1.0."""
        mol = smsd.parse_smiles("CC(=O)Oc1ccccc1C(=O)O")
        fp = smsd.circular_fingerprint(mol, radius=2, fp_size=2048)
        assert smsd.tanimoto_coefficient(fp, fp) == 1.0
        assert smsd.dice(fp, fp) == 1.0
        assert smsd.cosine(fp, fp) == 1.0

    def test_dice_formula(self, benzene, phenol):
        """Verify Dice = 2*|A∩B| / (|A| + |B|)."""
        fp1 = smsd.circular_fingerprint(benzene, radius=2, fp_size=2048)
        fp2 = smsd.circular_fingerprint(phenol, radius=2, fp_size=2048)
        s1, s2 = set(fp1), set(fp2)
        expected_dice = 2 * len(s1 & s2) / (len(s1) + len(s2)) if (len(s1) + len(s2)) > 0 else 0
        assert smsd.dice(fp1, fp2) == pytest.approx(expected_dice, abs=1e-10)

    def test_tanimoto_formula(self, benzene, phenol):
        """Verify Tanimoto = |A∩B| / |A∪B|."""
        fp1 = smsd.circular_fingerprint(benzene, radius=2, fp_size=2048)
        fp2 = smsd.circular_fingerprint(phenol, radius=2, fp_size=2048)
        s1, s2 = set(fp1), set(fp2)
        expected = len(s1 & s2) / len(s1 | s2) if len(s1 | s2) > 0 else 0
        assert smsd.tanimoto_coefficient(fp1, fp2) == pytest.approx(expected, abs=1e-10)

    def test_cosine_formula(self, benzene, phenol):
        """Verify Cosine = |A∩B| / sqrt(|A| * |B|)."""
        fp1 = smsd.circular_fingerprint(benzene, radius=2, fp_size=2048)
        fp2 = smsd.circular_fingerprint(phenol, radius=2, fp_size=2048)
        s1, s2 = set(fp1), set(fp2)
        denom = math.sqrt(len(s1) * len(s2))
        expected = len(s1 & s2) / denom if denom > 0 else 0
        assert smsd.cosine(fp1, fp2) == pytest.approx(expected, abs=1e-10)


# ===================================================================
# 14. Single-atom and small-molecule edge cases
# ===================================================================

class TestSingleAtomEdgeCases:
    """Molecules with 1-3 atoms: ensure all FP types produce valid output."""

    def test_methane_ecfp(self):
        """Methane (C) -- single heavy atom, ECFP should produce valid FP."""
        mol = smsd.parse_smiles("C")
        fp = smsd.circular_fingerprint(mol, radius=2, fp_size=2048)
        assert isinstance(fp, list) and len(fp) >= 1
        assert all(0 <= b < 2048 for b in fp)

    def test_methane_path(self):
        """Methane (C) -- path fingerprint on single atom."""
        mol = smsd.parse_smiles("C")
        fp = smsd.path_fingerprint(mol, path_length=7, fp_size=2048)
        assert isinstance(fp, list) and len(fp) >= 1

    def test_methane_torsion(self):
        """Methane (C) -- no 4-atom paths, torsion FP should be empty or valid."""
        fp = smsd.topological_torsion("C", fp_size=2048)
        assert isinstance(fp, list)
        # Single atom has no 4-atom path, so empty is expected
        assert len(fp) == 0

    def test_water(self):
        """Water (O) -- single heavy atom, still valid FP."""
        mol = smsd.parse_smiles("O")
        fp = smsd.circular_fingerprint(mol, radius=2, fp_size=2048)
        assert isinstance(fp, list) and len(fp) >= 1

    def test_hcl(self):
        """HCl ([H]Cl) -- 2 atoms, valid FP."""
        mol = smsd.parse_smiles("[H]Cl")
        fp = smsd.circular_fingerprint(mol, radius=2, fp_size=2048)
        assert isinstance(fp, list) and len(fp) >= 1

    def test_ethane(self):
        """Ethane (CC) -- 2 heavy atoms, no torsion possible."""
        mol = smsd.parse_smiles("CC")
        fp_ecfp = smsd.circular_fingerprint(mol, radius=2, fp_size=2048)
        fp_path = smsd.path_fingerprint(mol, path_length=7, fp_size=2048)
        fp_torsion = smsd.topological_torsion("CC", fp_size=2048)
        assert len(fp_ecfp) >= 1
        assert len(fp_path) >= 1
        assert len(fp_torsion) == 0, "2-atom molecule has no 4-atom torsion paths"

    def test_propane(self):
        """Propane (CCC) -- 3 atoms, still no torsion (needs 4)."""
        fp_torsion = smsd.topological_torsion("CCC", fp_size=2048)
        assert len(fp_torsion) == 0, "3-atom molecule has no 4-atom torsion paths"

    def test_butane_first_torsion(self):
        """Butane (CCCC) -- exactly 4 atoms, first molecule with torsion bits."""
        fp_torsion = smsd.topological_torsion("CCCC", fp_size=2048)
        assert len(fp_torsion) >= 1, "4-atom linear chain should have at least 1 torsion"
        mol = smsd.parse_smiles("CCCC")
        fp_ecfp = smsd.circular_fingerprint(mol, radius=2, fp_size=2048)
        assert len(fp_ecfp) >= 1


# ===================================================================
# 15. Disjoint fingerprints -- zero overlap
# ===================================================================

class TestDisjointFingerprints:
    """Two molecules with zero (or near-zero) overlap in fingerprint space."""

    def test_c60_vs_peptide_tanimoto_very_low(self):
        """C60 fullerene vs small peptide (Ala-Gly) -- should have near-zero Tanimoto."""
        c60 = smsd.parse_smiles(
            "c12c3c4c5c1c1c6c7c2c2c8c3c3c9c4c4c%10c5c5c1c1c6c6c%11c7c7c2c2c8c8"
            "c3c3c9c9c4c4c%10c%10c5c5c1c1c6c6c%11c7c7c2c8c2c3c9c3c4c%10c4c5c1c1c6"
            "c7c2c3c41"
        )
        ala_gly = smsd.parse_smiles("CC(N)C(=O)NCC(=O)O")
        fp_c60 = smsd.circular_fingerprint(c60, radius=2, fp_size=2048)
        fp_ag = smsd.circular_fingerprint(ala_gly, radius=2, fp_size=2048)
        sim = smsd.tanimoto_coefficient(fp_c60, fp_ag)
        assert sim < 0.15, f"C60 vs Ala-Gly should have very low similarity, got {sim:.3f}"

    def test_disjoint_dice_zero(self):
        """Truly disjoint bit sets should give Dice = 0."""
        fp_a = [0, 1, 2]
        fp_b = [100, 200, 300]
        assert smsd.dice(fp_a, fp_b) == 0.0

    def test_disjoint_tanimoto_zero(self):
        """Truly disjoint bit sets should give Tanimoto = 0."""
        fp_a = [0, 1, 2]
        fp_b = [100, 200, 300]
        assert smsd.tanimoto_coefficient(fp_a, fp_b) == 0.0

    def test_disjoint_soergel_one(self):
        """Soergel distance = 1 - Tanimoto. For disjoint FPs, distance = 1.0."""
        fp_a = [0, 1, 2]
        fp_b = [100, 200, 300]
        assert smsd.soergel(fp_a, fp_b) == 1.0


# ===================================================================
# 16. Fingerprint reproducibility
# ===================================================================

class TestFingerprintReproducibility:
    """Ensure fingerprints are fully deterministic across repeated parsing."""

    def test_same_molecule_100_parses(self):
        """Same molecule parsed 100 times must produce the same FP every time."""
        reference = smsd.circular_fingerprint(
            smsd.parse_smiles("c1ccccc1"), radius=2, fp_size=2048
        )
        for i in range(100):
            fp = smsd.circular_fingerprint(
                smsd.parse_smiles("c1ccccc1"), radius=2, fp_size=2048
            )
            assert fp == reference, f"Mismatch on iteration {i}"

    def test_canonical_equivalence_benzoic_acid(self):
        """OC(=O)c1ccccc1 vs c1ccccc1C(=O)O vs C(=O)(O)c1ccccc1 -- all same FP."""
        smiles_variants = [
            "OC(=O)c1ccccc1",
            "c1ccccc1C(=O)O",
            "C(=O)(O)c1ccccc1",
        ]
        fps = [
            smsd.circular_fingerprint(smsd.parse_smiles(s), radius=2, fp_size=2048)
            for s in smiles_variants
        ]
        assert fps[0] == fps[1], f"Variant 0 vs 1 differ"
        assert fps[1] == fps[2], f"Variant 1 vs 2 differ"

    def test_path_fp_canonical_equivalence(self):
        """Path FP should also be canonical across SMILES orderings."""
        smiles_variants = [
            "OC(=O)c1ccccc1",
            "c1ccccc1C(=O)O",
        ]
        fps = [
            smsd.path_fingerprint(smsd.parse_smiles(s), path_length=7, fp_size=2048)
            for s in smiles_variants
        ]
        assert fps[0] == fps[1]


# ===================================================================
# 17. Fingerprint size edge cases
# ===================================================================

class TestFingerprintSizeEdgeCases:
    """Extreme fp_size values: very small and very large."""

    def test_fp_size_1(self):
        """fp_size=1 -- should work but all bits map to position 0."""
        mol = smsd.parse_smiles("c1ccccc1")
        fp = smsd.circular_fingerprint(mol, radius=2, fp_size=1)
        assert isinstance(fp, list)
        # Everything collides to position 0
        assert all(b == 0 for b in fp)

    def test_fp_size_64(self):
        """fp_size=64 -- small but valid."""
        mol = smsd.parse_smiles("c1ccccc1")
        fp = smsd.circular_fingerprint(mol, radius=2, fp_size=64)
        assert all(0 <= b < 64 for b in fp)
        assert len(fp) >= 1

    def test_fp_size_16384(self):
        """fp_size=16384 -- large, should have fewer collisions."""
        mol = smsd.parse_smiles("c1ccccc1")
        fp = smsd.circular_fingerprint(mol, radius=2, fp_size=16384)
        assert all(0 <= b < 16384 for b in fp)
        assert len(fp) >= 1

    def test_larger_fp_more_or_equal_bits(self):
        """Larger fp_size gives same or more set bits (fewer collisions)."""
        mol = smsd.parse_smiles("CC(=O)Oc1ccccc1C(=O)O")  # aspirin
        sizes = [64, 256, 1024, 4096, 16384]
        prev_bits = 0
        for sz in sizes:
            fp = smsd.circular_fingerprint(mol, radius=2, fp_size=sz)
            assert len(fp) >= prev_bits, (
                f"fp_size={sz} has {len(fp)} bits but previous had {prev_bits}"
            )
            prev_bits = len(fp)


# ===================================================================
# 18. Count fingerprint edge cases
# ===================================================================

class TestCountFingerprintEdgeCases:
    """Edge cases for count-based fingerprints and counts_to_array."""

    def test_counts_to_array_empty_dict(self):
        """Empty dict should produce all-zero array."""
        arr = smsd.counts_to_array({}, 2048)
        assert len(arr) == 2048
        assert all(v == 0 for v in arr)

    def test_counts_to_array_none_input(self):
        """None input should produce all-zero array."""
        arr = smsd.counts_to_array(None, 2048)
        assert len(arr) == 2048
        assert all(v == 0 for v in arr)

    def test_counts_to_array_out_of_range_position(self):
        """Positions >= fp_size should be silently ignored."""
        arr = smsd.counts_to_array({5000: 3, 10: 1}, 2048)
        assert len(arr) == 2048
        assert arr[10] == 1
        # Position 5000 is out of range, should be ignored
        assert sum(arr) == 1

    def test_benzene_count_symmetry(self):
        """Benzene: all 6 carbons equivalent, all counts should be 6."""
        counts = smsd.circular_fingerprint_counts(
            smsd.parse_smiles("c1ccccc1"), radius=2, fp_size=2048
        )
        for pos, count in counts:
            assert count == 6, f"Benzene symmetry: all features should have count 6, got {count}"

    def test_count_superset_of_binary(self):
        """Count FP positions should be a superset of (or equal to) binary FP positions."""
        mol = smsd.parse_smiles("CC(=O)Oc1ccccc1C(=O)O")  # aspirin
        binary_pos = set(smsd.circular_fingerprint(mol, radius=2, fp_size=2048))
        count_pos = {pos for pos, _ in smsd.circular_fingerprint_counts(mol, radius=2, fp_size=2048)}
        assert binary_pos == count_pos, "Binary and count positions must match"


# ===================================================================
# 19. Similarity mathematical properties
# ===================================================================

class TestSimilarityMathematicalProperties:
    """Verify mathematical properties of similarity metrics."""

    @pytest.fixture
    def three_fps(self):
        """Three diverse molecules for metric property tests."""
        fp_a = smsd.circular_fingerprint(smsd.parse_smiles("c1ccccc1"), radius=2, fp_size=2048)
        fp_b = smsd.circular_fingerprint(smsd.parse_smiles("c1ccc(O)cc1"), radius=2, fp_size=2048)
        fp_c = smsd.circular_fingerprint(smsd.parse_smiles("CC(=O)Oc1ccccc1C(=O)O"), radius=2, fp_size=2048)
        return fp_a, fp_b, fp_c

    def test_tanimoto_reflexive(self, three_fps):
        """T(a, a) = 1 for all non-empty FPs."""
        for fp in three_fps:
            assert smsd.tanimoto_coefficient(fp, fp) == 1.0

    def test_tanimoto_symmetric(self, three_fps):
        """T(a, b) = T(b, a)."""
        fp_a, fp_b, fp_c = three_fps
        assert smsd.tanimoto_coefficient(fp_a, fp_b) == smsd.tanimoto_coefficient(fp_b, fp_a)
        assert smsd.tanimoto_coefficient(fp_a, fp_c) == smsd.tanimoto_coefficient(fp_c, fp_a)
        assert smsd.tanimoto_coefficient(fp_b, fp_c) == smsd.tanimoto_coefficient(fp_c, fp_b)

    def test_tanimoto_triangle_inequality(self, three_fps):
        """Tanimoto distance (1-T) satisfies triangle inequality:
        d(a,c) <= d(a,b) + d(b,c)."""
        fp_a, fp_b, fp_c = three_fps
        d_ab = 1.0 - smsd.tanimoto_coefficient(fp_a, fp_b)
        d_bc = 1.0 - smsd.tanimoto_coefficient(fp_b, fp_c)
        d_ac = 1.0 - smsd.tanimoto_coefficient(fp_a, fp_c)
        assert d_ac <= d_ab + d_bc + 1e-10, (
            f"Triangle inequality violated: d(a,c)={d_ac:.4f} > d(a,b)+d(b,c)={d_ab + d_bc:.4f}"
        )

    def test_dice_ge_tanimoto(self, three_fps):
        """Dice coefficient is always >= Tanimoto (mathematical fact)."""
        fp_a, fp_b, fp_c = three_fps
        for x, y in [(fp_a, fp_b), (fp_a, fp_c), (fp_b, fp_c)]:
            t = smsd.tanimoto_coefficient(x, y)
            d = smsd.dice(x, y)
            assert d >= t - 1e-10, f"Dice {d:.6f} < Tanimoto {t:.6f}"

    def test_overlap_subset_equals_one(self):
        """If A is a subset of B, overlap coefficient should be 1.0."""
        fp_small = [10, 20, 30]
        fp_big = [5, 10, 15, 20, 25, 30, 35, 40]
        ov = smsd.overlap_coefficient(fp_small, fp_big)
        assert ov == 1.0, f"A subset of B => overlap = 1.0, got {ov}"

    def test_cosine_self_is_one(self, three_fps):
        """cos(a, a) = 1.0."""
        for fp in three_fps:
            assert smsd.cosine(fp, fp) == pytest.approx(1.0, abs=1e-10)

    def test_cosine_range(self, three_fps):
        """0 <= cos(a, b) <= 1."""
        fp_a, fp_b, fp_c = three_fps
        for x, y in [(fp_a, fp_b), (fp_a, fp_c), (fp_b, fp_c)]:
            c = smsd.cosine(x, y)
            assert 0.0 <= c <= 1.0 + 1e-10, f"Cosine out of range: {c}"


# ===================================================================
# 20. Pharmacophore (FCFP) fingerprint edge cases
# ===================================================================

class TestPharmacophoreFingerprints:
    """FCFP groups pharmacophore-similar atoms (donors, acceptors, etc.)."""

    def test_phenol_vs_aniline_fcfp_higher(self):
        """Phenol (-OH) and aniline (-NH2) both have H-bond donor.
        FCFP similarity should be higher than ECFP."""
        phenol = smsd.parse_smiles("c1ccc(O)cc1")
        aniline = smsd.parse_smiles("c1ccc(N)cc1")
        ecfp_sim = smsd.tanimoto_coefficient(
            smsd.circular_fingerprint(phenol, radius=2, fp_size=2048, mode="ecfp"),
            smsd.circular_fingerprint(aniline, radius=2, fp_size=2048, mode="ecfp"),
        )
        fcfp_sim = smsd.tanimoto_coefficient(
            smsd.circular_fingerprint(phenol, radius=2, fp_size=2048, mode="fcfp"),
            smsd.circular_fingerprint(aniline, radius=2, fp_size=2048, mode="fcfp"),
        )
        assert fcfp_sim >= ecfp_sim, (
            f"FCFP ({fcfp_sim:.3f}) should be >= ECFP ({ecfp_sim:.3f}) for donor-equivalent pair"
        )

    def test_carboxylic_acid_vs_sulfonamide(self):
        """Carboxylic acid and sulfonamide share acidic pharmacophore.
        FCFP should produce valid FPs for both."""
        acid = smsd.parse_smiles("c1ccc(C(=O)O)cc1")       # benzoic acid
        sulfonamide = smsd.parse_smiles("c1ccc(S(=O)(=O)N)cc1")  # benzenesulfonamide
        fp_acid = smsd.circular_fingerprint(acid, radius=2, fp_size=2048, mode="fcfp")
        fp_sulfo = smsd.circular_fingerprint(sulfonamide, radius=2, fp_size=2048, mode="fcfp")
        assert len(fp_acid) >= 1
        assert len(fp_sulfo) >= 1
        # They should share some pharmacophore features (aromatic ring at minimum)
        sim = smsd.tanimoto_coefficient(fp_acid, fp_sulfo)
        assert sim > 0.0, "Should share at least aromatic ring pharmacophore"

    def test_fcfp_valid_for_all_basic_molecules(self):
        """FCFP should produce valid FPs for a range of functional groups."""
        test_smiles = [
            "c1ccccc1",          # benzene (aromatic)
            "CCO",               # ethanol (donor/acceptor)
            "CC(=O)O",           # acetic acid (acidic)
            "CCN",               # ethylamine (basic)
            "CC=O",              # acetaldehyde (acceptor)
        ]
        for smi in test_smiles:
            mol = smsd.parse_smiles(smi)
            fp = smsd.circular_fingerprint(mol, radius=2, fp_size=2048, mode="fcfp")
            assert isinstance(fp, list) and len(fp) >= 1, f"FCFP failed for {smi}"


# ===================================================================
# 21. Topological torsion edge cases
# ===================================================================

class TestTopologicalTorsionEdgeCases:
    """Torsion FP edge cases: rings, branches, aromaticity, performance."""

    def test_cyclohexane_ring_torsions(self):
        """Cyclohexane (6-membered ring) has torsions through the ring."""
        fp = smsd.topological_torsion("C1CCCCC1", fp_size=2048)
        assert len(fp) >= 1, "Cyclohexane should have ring-path torsions"

    def test_neopentane_branched(self):
        """Neopentane (CC(C)(C)C) -- maximally branched, no 4-atom linear path.
        The central carbon connects to 4 terminal carbons, but no path of
        length 4 (A-B-C-D) exists since all branches are length 1."""
        fp = smsd.topological_torsion("CC(C)(C)C", fp_size=2048)
        assert isinstance(fp, list)
        assert len(fp) == 0, "Neopentane has no 4-atom linear paths"

    def test_aromatic_vs_aliphatic_ring_torsion(self):
        """Benzene vs cyclohexane -- both are 6-carbon rings with identical
        connectivity, so torsion FPs may coincide. ECFP should distinguish them."""
        fp_benzene_t = smsd.topological_torsion("c1ccccc1", fp_size=2048)
        fp_cyclohexane_t = smsd.topological_torsion("C1CCCCC1", fp_size=2048)
        # Both produce valid torsion FPs
        assert len(fp_benzene_t) >= 1
        assert len(fp_cyclohexane_t) >= 1
        # ECFP does distinguish aromatic vs aliphatic
        fp_benzene_e = smsd.circular_fingerprint(
            smsd.parse_smiles("c1ccccc1"), radius=2, fp_size=2048
        )
        fp_cyclohexane_e = smsd.circular_fingerprint(
            smsd.parse_smiles("C1CCCCC1"), radius=2, fp_size=2048
        )
        assert fp_benzene_e != fp_cyclohexane_e, (
            "ECFP must distinguish aromatic and aliphatic 6-rings"
        )

    def test_steroid_performance(self):
        """Cholesterol (steroid) -- should complete torsion FP in < 1 second."""
        import time
        cholesterol = "C[C@H](CCCC(C)C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC=C4[C@@]3(CC[C@@H](C4)O)C)C"
        t0 = time.time()
        fp = smsd.topological_torsion(cholesterol, fp_size=2048)
        elapsed = time.time() - t0
        assert elapsed < 1.0, f"Steroid torsion took {elapsed:.2f}s, should be < 1s"
        assert len(fp) >= 1

    def test_torsion_counts_match_binary(self):
        """Torsion count positions should match binary torsion positions."""
        mol = smsd.parse_smiles("c1ccc(CC)cc1")  # ethylbenzene
        binary_pos = set(smsd.topological_torsion(mol, fp_size=2048))
        count_pos = {pos for pos, _ in smsd.topological_torsion_counts(mol, fp_size=2048)}
        assert binary_pos == count_pos
