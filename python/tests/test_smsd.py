# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
# Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
# See the NOTICE file for attribution, trademark, and algorithm IP terms.
"""
Tests for the SMSD Python bindings.

Run with:  pytest tests/test_smsd.py -v
"""
import pytest
import smsd
from smsd import (
    parse_smiles,
    to_smiles,
    to_smarts,
    read_mol_block,
    write_mol_block,
    write_mol_block_v3000,
    write_sdf_record,
    is_substructure,
    find_substructure,
    find_mcs,
    similarity_upper_bound,
    screen_targets,
    path_fingerprint,
    mcs_fingerprint,
    fingerprint_subset,
    analyze_fp_quality,
    ChemOptions,
    McsOptions,
    MolGraph,
    MolGraphBuilder,
    BondOrderMode,
    AromaticityMode,
    RingFusionMode,
)


# ===========================================================================
# SMILES parsing
# ===========================================================================

class TestParseSMILES:
    """Test SMILES parsing for various molecules."""

    def test_benzene(self):
        mol = parse_smiles("c1ccccc1")
        assert mol.n == 6
        for i in range(6):
            assert mol.atomic_num[i] == 6
            assert mol.aromatic[i]
            assert mol.ring[i]

    def test_aspirin(self):
        mol = parse_smiles("CC(=O)Oc1ccccc1C(=O)O")
        assert mol.n == 13
        # Count carbons and oxygens
        n_c = sum(1 for z in mol.atomic_num if z == 6)
        n_o = sum(1 for z in mol.atomic_num if z == 8)
        assert n_c == 9
        assert n_o == 4

    def test_caffeine(self):
        mol = parse_smiles("Cn1cnc2c1c(=O)n(c(=O)n2C)C")
        assert mol.n == 14
        # Should have N and O atoms
        n_n = sum(1 for z in mol.atomic_num if z == 7)
        n_o = sum(1 for z in mol.atomic_num if z == 8)
        assert n_n == 4
        assert n_o == 2

    def test_ethanol(self):
        mol = parse_smiles("CCO")
        assert mol.n == 3
        assert mol.atomic_num[2] == 8

    def test_empty_returns_empty_molgraph(self):
        mol = parse_smiles("")
        assert mol.n == 0

    def test_invalid_raises(self):
        with pytest.raises(Exception):
            parse_smiles("[Zq]")

    def test_round_trip(self):
        """Parse -> write -> re-parse should give the same atom count."""
        mol1 = parse_smiles("c1ccccc1")
        smi = to_smiles(mol1)
        mol2 = parse_smiles(smi)
        assert mol2.n == mol1.n

    def test_charged_atom(self):
        mol = parse_smiles("[NH4+]")
        assert mol.n == 1
        assert mol.atomic_num[0] == 7
        assert mol.formal_charge[0] == 1

    def test_isotope(self):
        mol = parse_smiles("[13CH4]")
        assert mol.n == 1
        assert mol.mass_number[0] == 13

    def test_disconnected(self):
        mol = parse_smiles("C.C")
        assert mol.n == 2


# ===========================================================================
# SMILES writer
# ===========================================================================

class TestToSMILES:
    """Test canonical SMILES generation."""

    def test_canonical_benzene(self):
        mol = parse_smiles("c1ccccc1")
        smi = to_smiles(mol)
        assert isinstance(smi, str)
        assert len(smi) > 0

    def test_canonical_aspirin(self):
        mol = parse_smiles("CC(=O)Oc1ccccc1C(=O)O")
        smi = to_smiles(mol)
        assert isinstance(smi, str)
        # Re-parse should give same structure
        mol2 = parse_smiles(smi)
        assert mol2.n == 13


# ===========================================================================
# Substructure search
# ===========================================================================

class TestSubstructure:
    """Test substructure (subgraph isomorphism) search."""

    def test_benzene_in_phenol(self):
        benzene = parse_smiles("c1ccccc1")
        phenol = parse_smiles("c1ccc(O)cc1")
        assert is_substructure(benzene, phenol) is True

    def test_phenol_not_in_benzene(self):
        benzene = parse_smiles("c1ccccc1")
        phenol = parse_smiles("c1ccc(O)cc1")
        assert is_substructure(phenol, benzene) is False

    def test_benzene_in_naphthalene(self):
        benzene = parse_smiles("c1ccccc1")
        naphthalene = parse_smiles("c1ccc2ccccc2c1")
        assert is_substructure(benzene, naphthalene) is True

    def test_find_substructure_mapping(self):
        benzene = parse_smiles("c1ccccc1")
        phenol = parse_smiles("c1ccc(O)cc1")
        mapping = find_substructure(benzene, phenol)
        assert len(mapping) == 6  # all 6 benzene atoms mapped
        # Each pair should be (query_idx, target_idx)
        for qi, ti in mapping:
            assert 0 <= qi < 6
            assert 0 <= ti < 7

    def test_self_substructure(self):
        mol = parse_smiles("CCO")
        assert is_substructure(mol, mol) is True

    def test_empty_query_is_substructure(self):
        """An empty molecule is a substructure of anything."""
        mol = parse_smiles("C")
        empty = MolGraph()
        assert is_substructure(empty, mol) is True

    def test_small_query_fastpath_respects_formal_charge(self):
        query = parse_smiles("[NH3+]C(=O)O")
        target = parse_smiles("NC(=O)[OH2+]")
        opts = ChemOptions()
        opts.match_formal_charge = True
        assert is_substructure(query, target, opts) is False

    def test_small_query_fastpath_respects_isotopes(self):
        query = parse_smiles("[13CH3]O")
        target = parse_smiles("[12CH3]O")
        opts = ChemOptions()
        opts.match_isotope = True
        assert is_substructure(query, target, opts) is False

    def test_small_query_fastpath_respects_chirality(self):
        query = parse_smiles("N[C@@H](C)C(=O)O")
        target = parse_smiles("N[C@H](C)C(=O)O")
        opts = ChemOptions()
        opts.use_chirality = True
        assert is_substructure(query, target, opts) is False

    def test_small_query_fastpath_respects_bond_stereo(self):
        query = parse_smiles("F/C=C/Cl")
        target = parse_smiles(r"F/C=C\Cl")
        opts = ChemOptions()
        opts.use_bond_stereo = True
        assert is_substructure(query, target, opts) is False

    def test_large_sparse_target_preserves_bond_stereo(self):
        query = parse_smiles("F/C=C/Cl")
        target = parse_smiles("F/C=C\\Cl." + ("C" * 210))
        opts = ChemOptions()
        opts.use_bond_stereo = True
        assert is_substructure(query, target, opts) is False


# ===========================================================================
# MCS (Maximum Common Substructure)
# ===========================================================================

class TestMCS:
    """Test MCS search."""

    def test_benzene_toluene(self):
        benzene = parse_smiles("c1ccccc1")
        toluene = parse_smiles("Cc1ccccc1")
        mapping = find_mcs(benzene, toluene)
        assert len(mapping) == 6  # benzene ring is the MCS

    def test_aspirin_acetaminophen(self):
        aspirin = parse_smiles("CC(=O)Oc1ccccc1C(=O)O")
        acetaminophen = parse_smiles("CC(=O)Nc1ccc(O)cc1")
        mapping = find_mcs(aspirin, acetaminophen)
        assert len(mapping) >= 7  # at least the ring + some shared atoms

    def test_identical_molecules(self):
        mol = parse_smiles("c1ccccc1")
        mapping = find_mcs(mol, mol)
        assert len(mapping) == 6

    def test_mcs_returns_dict(self):
        g1 = parse_smiles("CCO")
        g2 = parse_smiles("CCCO")
        mapping = find_mcs(g1, g2)
        assert isinstance(mapping, dict)
        for k, v in mapping.items():
            assert isinstance(k, int)
            assert isinstance(v, int)

    def test_mcs_with_options(self):
        g1 = parse_smiles("c1ccccc1")
        g2 = parse_smiles("c1ccc(O)cc1")
        opts = McsOptions()
        opts.timeout_ms = 5000
        mapping = find_mcs(g1, g2, ChemOptions(), opts)
        assert len(mapping) == 6

    def test_ring_matches_ring_only_requires_ring_parity(self):
        ring = parse_smiles("C1CCCCC1")
        chain = parse_smiles("CCCCC")

        strict = ChemOptions()
        strict.ring_matches_ring_only = True
        assert len(find_mcs(ring, chain, strict)) == 0
        assert len(find_mcs(chain, ring, strict)) == 0
        assert is_substructure(chain, ring, strict) is False

        relaxed = ChemOptions()
        relaxed.ring_matches_ring_only = False
        assert len(find_mcs(ring, chain, relaxed)) == 5
        assert len(find_mcs(chain, ring, relaxed)) == 5
        assert is_substructure(chain, ring, relaxed) is True


# ===========================================================================
# RASCAL (similarity upper bound)
# ===========================================================================

class TestRASCAL:
    """Test RASCAL similarity screening."""

    def test_similar_molecules(self):
        benzene = parse_smiles("c1ccccc1")
        toluene = parse_smiles("Cc1ccccc1")
        sim = similarity_upper_bound(benzene, toluene)
        assert sim > 0.5
        assert sim <= 1.0

    def test_identical_similarity(self):
        mol = parse_smiles("c1ccccc1")
        sim = similarity_upper_bound(mol, mol)
        assert sim == pytest.approx(1.0)

    def test_dissimilar_molecules(self):
        methane = parse_smiles("C")
        cholesterol = parse_smiles(
            "CC(CCCC(C)C)C1CCC2C1(CCC3C2CC=C4C3(CCC(C4)O)C)C"
        )
        sim = similarity_upper_bound(methane, cholesterol)
        assert sim < 0.5

    def test_convenience_similarity(self):
        """Test the convenience wrapper that accepts SMILES strings."""
        sim = smsd.similarity("c1ccccc1", "Cc1ccccc1")
        assert sim > 0.5

    def test_empty_molecule(self):
        mol = parse_smiles("C")
        empty = MolGraph()
        sim = similarity_upper_bound(mol, empty)
        assert sim == 0.0


# ===========================================================================
# Fingerprints
# ===========================================================================

class TestFingerprint:
    """Test path and MCS fingerprints."""

    def test_path_fingerprint_basic(self):
        mol = parse_smiles("c1ccccc1")
        fp = path_fingerprint(mol)
        assert isinstance(fp, list)
        assert len(fp) > 0
        # All bits should be in range [0, 2048)
        assert all(0 <= b < 2048 for b in fp)

    def test_mcs_fingerprint_basic(self):
        mol = parse_smiles("c1ccccc1")
        fp = mcs_fingerprint(mol)
        assert isinstance(fp, list)
        assert len(fp) > 0
        assert all(0 <= b < 2048 for b in fp)

    def test_self_subset(self):
        """A fingerprint should be a subset of itself."""
        mol = parse_smiles("c1ccccc1")
        fp = path_fingerprint(mol)
        assert fingerprint_subset(fp, fp) is True

    def test_subset_property(self):
        """Substructure fingerprint should be subset of superstructure."""
        benzene = parse_smiles("c1ccccc1")
        phenol = parse_smiles("c1ccc(O)cc1")
        fp_benz = path_fingerprint(benzene)
        fp_phen = path_fingerprint(phenol)
        # benzene paths should be a subset of phenol paths
        assert fingerprint_subset(fp_benz, fp_phen) is True

    def test_fingerprint_custom_size(self):
        mol = parse_smiles("CCO")
        fp = path_fingerprint(mol, path_length=5, fp_size=1024)
        assert all(0 <= b < 1024 for b in fp)

    def test_mcs_fingerprint_subset(self):
        benzene = parse_smiles("c1ccccc1")
        toluene = parse_smiles("Cc1ccccc1")
        fp_benz = mcs_fingerprint(benzene)
        fp_tol = mcs_fingerprint(toluene)
        assert len(fp_benz) > 0
        assert len(fp_tol) > 0
        assert len(fp_tol) >= len(fp_benz)
        assert len(set(fp_benz) & set(fp_tol)) > 0

    def test_analyze_fp_quality(self):
        mol = parse_smiles("c1ccccc1")
        fp = path_fingerprint(mol)
        quality = analyze_fp_quality(fp)
        assert "set_bits" in quality
        assert "density" in quality
        assert "fp_size" in quality
        assert "expected_collision_rate" in quality
        assert "is_saturated" in quality
        assert quality["set_bits"] == len(fp)
        assert quality["fp_size"] == 2048
        assert 0.0 <= quality["density"] <= 1.0

    def test_convenience_fingerprint(self):
        """Test the convenience wrapper."""
        fp_path = smsd.fingerprint("CCO", kind="path")
        fp_mcs = smsd.fingerprint("CCO", kind="mcs")
        assert isinstance(fp_path, list)
        assert isinstance(fp_mcs, list)

    def test_tanimoto(self):
        fp1 = smsd.fingerprint("c1ccccc1")
        fp2 = smsd.fingerprint("Cc1ccccc1")
        sim = smsd.tanimoto(fp1, fp2)
        assert 0.0 <= sim <= 1.0
        # Same molecule should give 1.0
        assert smsd.tanimoto(fp1, fp1) == pytest.approx(1.0)


# ===========================================================================
# Tautomer-aware matching
# ===========================================================================

class TestTautomer:
    """Test tautomer-aware MCS."""

    def test_keto_enol_tautomer_profile(self):
        """Keto and enol forms should have larger MCS with tautomer-aware opts."""
        # 2-hydroxypyridine (enol) vs 2-pyridinone (keto)
        enol = parse_smiles("Oc1ccccn1")
        keto = parse_smiles("O=C1C=CC=CN1")
        # Standard MCS
        mapping_strict = find_mcs(enol, keto)
        # Tautomer-aware MCS
        taut_opts = ChemOptions.tautomer_profile()
        mapping_taut = find_mcs(enol, keto, taut_opts)
        # Tautomer-aware should find at least as many atoms
        assert len(mapping_taut) >= len(mapping_strict)

    def test_tautomer_profile_creation(self):
        opts = ChemOptions.tautomer_profile()
        assert opts.tautomer_aware is True
        assert opts.match_formal_charge is False

    def test_convenience_mcs_tautomer(self):
        """Test the convenience mcs() with tautomer_aware flag."""
        mapping = smsd.mcs("Oc1ccccn1", "O=C1C=CC=CN1", tautomer_aware=True)
        assert isinstance(mapping, dict)
        assert len(mapping) > 0


# ===========================================================================
# Batch screening
# ===========================================================================

class TestBatchScreening:
    """Test batch RASCAL screening."""

    # 20 diverse SMILES for building a screening library
    LIBRARY_SMILES = [
        "c1ccccc1",                          # benzene
        "Cc1ccccc1",                         # toluene
        "c1ccc(O)cc1",                       # phenol
        "c1ccc(N)cc1",                       # aniline
        "c1ccc(C(=O)O)cc1",                 # benzoic acid
        "CCO",                               # ethanol
        "CCCO",                              # propanol
        "CCCCO",                             # butanol
        "CC(=O)O",                           # acetic acid
        "CC(=O)Oc1ccccc1C(=O)O",            # aspirin
        "CC(=O)Nc1ccc(O)cc1",               # acetaminophen
        "c1ccc2ccccc2c1",                    # naphthalene
        "CC(C)Cc1ccc(cc1)C(C)C(=O)O",       # ibuprofen
        "Cn1cnc2c1c(=O)n(c(=O)n2C)C",       # caffeine
        "C1CCCCC1",                          # cyclohexane
        "c1ccncc1",                          # pyridine
        "c1cc[nH]c1",                        # pyrrole
        "c1ccsc1",                           # thiophene
        "C=C",                               # ethylene
        "C#C",                               # acetylene
    ]

    def test_screen_100_molecules(self):
        """Screen a library of molecules against a query."""
        # Build library by repeating and varying the 20 SMILES 5 times
        library = []
        for _ in range(5):
            for smi in self.LIBRARY_SMILES:
                library.append(parse_smiles(smi))
        assert len(library) == 100

        query = parse_smiles("c1ccccc1")  # benzene
        hits = screen_targets(query, library, 0.5)
        assert isinstance(hits, list)
        # Should find at least benzene copies (indices 0, 20, 40, 60, 80)
        assert len(hits) >= 5
        # All hits should be valid indices
        assert all(0 <= idx < 100 for idx in hits)

    def test_screen_high_threshold(self):
        """With very high threshold, only identical molecules pass."""
        library = [parse_smiles(s) for s in self.LIBRARY_SMILES]
        query = parse_smiles("c1ccccc1")
        hits = screen_targets(query, library, 0.99)
        # Only benzene (index 0) should pass
        assert 0 in hits


# ===========================================================================
# ChemOptions and McsOptions configuration
# ===========================================================================

class TestOptions:
    """Test options configuration."""

    def test_chem_options_defaults(self):
        opts = ChemOptions()
        assert opts.match_atom_type is True
        assert opts.match_formal_charge is True
        assert opts.tautomer_aware is False
        assert opts.complete_rings_only is False

    def test_chem_options_modification(self):
        opts = ChemOptions()
        opts.tautomer_aware = True
        opts.complete_rings_only = True
        assert opts.tautomer_aware is True
        assert opts.complete_rings_only is True

    def test_mcs_options_defaults(self):
        opts = McsOptions()
        assert opts.connected_only is True
        assert opts.induced is False
        assert opts.timeout_ms == -1

    def test_mcs_options_modification(self):
        opts = McsOptions()
        opts.timeout_ms = 5000
        opts.connected_only = False
        assert opts.timeout_ms == 5000
        assert opts.connected_only is False

    def test_bond_order_mode_enum(self):
        assert BondOrderMode.STRICT is not None
        assert BondOrderMode.LOOSE is not None
        assert BondOrderMode.ANY is not None

    def test_ring_fusion_mode_enum(self):
        assert RingFusionMode.IGNORE is not None
        assert RingFusionMode.PERMISSIVE is not None
        assert RingFusionMode.STRICT is not None

    def test_chem_options_repr(self):
        opts = ChemOptions()
        r = repr(opts)
        assert "ChemOptions" in r

    def test_mcs_options_repr(self):
        opts = McsOptions()
        r = repr(opts)
        assert "McsOptions" in r


# ===========================================================================
# MolGraph properties
# ===========================================================================

class TestMolGraph:
    """Test MolGraph properties and methods."""

    def test_mol_repr(self):
        mol = parse_smiles("CCO")
        assert "MolGraph" in repr(mol)
        assert "3" in repr(mol)

    def test_mol_len(self):
        mol = parse_smiles("c1ccccc1")
        assert len(mol) == 6

    def test_bond_order(self):
        mol = parse_smiles("C=C")
        assert mol.bond_order(0, 1) == 2

    def test_has_bond(self):
        mol = parse_smiles("CC")
        assert mol.has_bond(0, 1) is True

    def test_no_bond(self):
        mol = parse_smiles("C.C")
        assert mol.has_bond(0, 1) is False

    def test_builder(self):
        """Build a molecule using MolGraphBuilder."""
        builder = MolGraphBuilder()
        builder.atom_count(3)
        builder.atomic_numbers([6, 6, 8])
        builder.formal_charges([0, 0, 0])
        builder.ring_flags([0, 0, 0])
        builder.aromatic_flags([0, 0, 0])
        builder.neighbors([[1], [0, 2], [1]])
        builder.bond_orders([[1], [1, 1], [1]])
        mol = builder.build()
        assert mol.n == 3
        assert mol.atomic_num[0] == 6
        assert mol.atomic_num[2] == 8

    def test_builder_accepts_dense_bond_matrix_rows(self):
        builder = MolGraphBuilder()
        builder.atom_count(2)
        builder.atomic_numbers([6, 6])
        builder.ring_flags([0, 0])
        builder.aromatic_flags([0, 0])
        builder.neighbors([[1], [0]])
        builder.bond_orders([[0, 2], [2, 0]])

        mol = builder.build()
        assert mol.bond_order(0, 1) == 2

    def test_builder_rejects_asymmetric_adjacency(self):
        builder = MolGraphBuilder()
        builder.atom_count(2)
        builder.atomic_numbers([6, 6])
        builder.neighbors([[1], []])

        with pytest.raises(ValueError, match="symmetric"):
            builder.build()

    def test_builder_rejects_duplicate_neighbors(self):
        builder = MolGraphBuilder()
        builder.atom_count(2)
        builder.atomic_numbers([6, 6])
        builder.neighbors([[1, 1], [0]])

        with pytest.raises(ValueError, match="duplicate neighbor"):
            builder.build()

    def test_builder_rejects_inconsistent_bond_orders(self):
        builder = MolGraphBuilder()
        builder.atom_count(2)
        builder.atomic_numbers([6, 6])
        builder.neighbors([[1], [0]])
        builder.bond_orders([[1], [2]])

        with pytest.raises(ValueError, match="bondOrders must be symmetric"):
            builder.build()

    def test_builder_rejects_duplicate_atom_ids(self):
        builder = MolGraphBuilder()
        builder.atom_count(2)
        builder.atomic_numbers([6, 6])
        builder.neighbors([[1], [0]])
        builder.atom_ids([9, 9])

        with pytest.raises(ValueError, match="atomIds must be unique"):
            builder.build()

    def test_morgan_ranks(self):
        mol = parse_smiles("c1ccccc1")
        to_smiles(mol)
        assert len(mol.morgan_rank) == 6

    def test_canonical_labels(self):
        mol = parse_smiles("c1ccccc1")
        to_smiles(mol)
        assert len(mol.canonical_label) == 6


@pytest.mark.parametrize(
    ("name", "smiles"),
    [
        ("drug_like_statin", "CC(C)c1nc(C(=O)N(CCO)CCO)c(c(n1)C(C)(C)C)S(=O)(=O)N"),
        ("nucleotide_like", "Nc1ncnc2n(cnc12)[C@@H]1O[C@H](COP(=O)(O)O)[C@@H](O)[C@H]1O"),
        ("natural_product_like", "CC1CCC2C(C1)CCC3C2CCC4=CC(=O)CCC34C"),
        ("macrocycle_like", "O=C1OCC(CC(=O)NCCO)OC(=O)CCCCC1"),
        ("peptide_like", "NCC(=O)N[C@@H](CO)C(=O)NCC(=O)O"),
        ("heteroaromatic_cofactor_like", "Nc1nc2n(cc(CO)c2[nH]1)C"),
    ],
)
def test_molgraph_category_coverage(name, smiles):
    """Representative chemistry classes should parse, canonicalize, and self-match."""
    mol = parse_smiles(smiles)
    assert mol.n > 0, f"{name}: parse_smiles should produce atoms"
    can = to_smiles(mol)
    assert isinstance(can, str) and can, f"{name}: canonical SMILES should not be empty"
    reparsed = parse_smiles(can)
    assert reparsed.n == mol.n, f"{name}: canonical round-trip should preserve atom count"
    mapping = find_mcs(mol, reparsed)
    assert len(mapping) == mol.n, f"{name}: self-equivalent canonical round-trip should preserve full MCS"


# ===========================================================================
# MCS Chemical Validity Tests
# ===========================================================================
class TestMcsChemicalValidity:
    """Verify MCS results obey basic chemistry invariants."""

    @staticmethod
    def _mapping_preserves_query_bonds(query, target, mapping):
        for qi, tj in mapping.items():
            for qk in range(len(query)):
                if qk <= qi or qk not in mapping:
                    continue
                if query.bond_order(qi, qk) == 0:
                    continue
                tk = mapping[qk]
                if target.bond_order(tj, tk) == 0:
                    return False
        return True

    def test_ethanol_self_mcs(self):
        """Ethanol (CCO) self-MCS must be exactly 3 heavy atoms."""
        eth = parse_smiles("CCO")
        assert len(eth) == 3, f"Ethanol has 3 heavy atoms, got {len(eth)}"
        mcs = find_mcs(eth, eth)
        assert len(mcs) == 3, f"Ethanol self-MCS must be 3, got {len(mcs)}"

    def test_piperazine_self_mcs(self):
        """Piperazine (C1CNCCN1) self-MCS must be exactly 6 heavy atoms."""
        pip = parse_smiles("C1CNCCN1")
        assert len(pip) == 6, f"Piperazine has 6 heavy atoms, got {len(pip)}"
        mcs = find_mcs(pip, pip)
        assert len(mcs) == 6, f"Piperazine self-MCS must be 6, not 17. Got {len(mcs)}"

    def test_benzene_toluene_mcs(self):
        """Benzene / toluene MCS = 6 (the benzene ring)."""
        benz = parse_smiles("c1ccccc1")
        tol = parse_smiles("Cc1ccccc1")
        mcs = find_mcs(benz, tol)
        assert len(mcs) == 6, f"Benzene/toluene MCS must be 6, got {len(mcs)}"

    def test_strychnine_quinine_mcs_is_direction_stable(self):
        strychnine = parse_smiles("C1CN2CC3=CCOC4CC(=O)N5C6C4C3CC2C61C7=CC=CC=C75")
        quinine = parse_smiles("COC1=CC2=C(C=CN=C2C=C1)C(C3CC4CCN3CC4C=C)O")

        forward = find_mcs(strychnine, quinine, timeout_ms=10000)
        reverse = find_mcs(quinine, strychnine, timeout_ms=10000)

        assert len(forward) >= 11
        assert len(reverse) >= 11
        assert len(forward) == len(reverse)
        assert self._mapping_preserves_query_bonds(strychnine, quinine, forward)
        assert self._mapping_preserves_query_bonds(quinine, strychnine, reverse)

    def test_mcs_never_exceeds_smaller_molecule(self):
        """MCS size must never exceed min(n1, n2) heavy atoms."""
        pairs = [
            ("CCO", "CCCO"),                        # ethanol vs propanol
            ("c1ccccc1", "c1ccc2ccccc2c1"),          # benzene vs naphthalene
            ("C1CNCCN1", "C1CN(c2ccccc2)CCN1"),      # piperazine vs phenylpiperazine
            ("CC(=O)O", "CC(=O)Oc1ccccc1"),          # acetic acid vs phenyl acetate
        ]
        for smi1, smi2 in pairs:
            g1, g2 = parse_smiles(smi1), parse_smiles(smi2)
            mcs = find_mcs(g1, g2)
            min_size = min(len(g1), len(g2))
            assert len(mcs) <= min_size, (
                f"MCS({smi1}, {smi2}) = {len(mcs)} exceeds min(n1,n2) = {min_size}"
            )

    def test_strict_chirality_mcs_preserves_query_edges(self):
        query = parse_smiles("N[C@@H](C)C(=O)O")
        target = parse_smiles("N[C@H](C)C(=O)O")
        opts = ChemOptions()
        opts.use_chirality = True

        mapping = find_mcs(query, target, opts)
        assert mapping
        assert len(mapping) < len(query)
        assert self._mapping_preserves_query_bonds(query, target, mapping)

    def test_strict_bond_stereo_mcs_preserves_query_edges(self):
        query = parse_smiles("F/C=C/Cl")
        target = parse_smiles(r"F/C=C\Cl")
        opts = ChemOptions()
        opts.use_bond_stereo = True

        mapping = find_mcs(query, target, opts)
        assert mapping
        assert len(mapping) < len(query)
        assert self._mapping_preserves_query_bonds(query, target, mapping)


class TestRdkitIndexTranslation:
    """Regression tests for mcs_rdkit_native() atom index translation.

    These reproduce the v6.1.0 bug where small molecules returned
    element-mismatched mappings (C mapped to O) and out-of-bounds indices.
    """

    @staticmethod
    def _check_mapping(smi1, smi2, min_mcs=1):
        """Verify mcs_rdkit_native returns element-correct, in-bounds mapping."""
        try:
            from rdkit import Chem
        except ImportError:
            pytest.skip("RDKit not installed")

        mol1 = Chem.MolFromSmiles(smi1)
        mol2 = Chem.MolFromSmiles(smi2)
        n1, n2 = mol1.GetNumAtoms(), mol2.GetNumAtoms()

        from smsd import mcs_rdkit_native
        mapping = mcs_rdkit_native(mol1, mol2)

        # Size check
        assert len(mapping) >= min_mcs, (
            f"MCS({smi1}, {smi2}) = {len(mapping)}, expected >= {min_mcs}"
        )

        for q_idx, t_idx in mapping.items():
            # Bounds check
            assert 0 <= q_idx < n1, (
                f"Query index {q_idx} out of bounds (n1={n1}) for {smi1}"
            )
            assert 0 <= t_idx < n2, (
                f"Target index {t_idx} out of bounds (n2={n2}) for {smi2}"
            )
            # Element match check
            q_elem = mol1.GetAtomWithIdx(q_idx).GetAtomicNum()
            t_elem = mol2.GetAtomWithIdx(t_idx).GetAtomicNum()
            assert q_elem == t_elem, (
                f"Element mismatch at Q[{q_idx}]={q_elem} -> T[{t_idx}]={t_elem} "
                f"for {smi1} vs {smi2}"
            )

    def test_ester_hydrolysis_small(self):
        """v6.1.0 bug: CC(=O)OC vs CC(=O)O returned {1:3,2:1,3:0,4:2} — all wrong."""
        self._check_mapping("CC(=O)OC", "CC(=O)O", min_mcs=4)

    def test_amide_vs_acid_small(self):
        """v6.1.0 bug: CC(=O)N vs CC(=O)O returned index 4 (out of bounds)."""
        self._check_mapping("CC(=O)N", "CC(=O)O", min_mcs=3)

    def test_ethanol_self(self):
        """Ethanol self-match — simplest round-trip."""
        self._check_mapping("CCO", "CCO", min_mcs=3)

    def test_methanol_vs_ethanol(self):
        """Small molecule pair: 2 atoms vs 3 atoms."""
        self._check_mapping("CO", "CCO", min_mcs=2)

    def test_acetic_acid_vs_propionic(self):
        """4 atoms vs 5 atoms."""
        self._check_mapping("CC(=O)O", "CCC(=O)O", min_mcs=4)

    def test_benzene_phenol(self):
        """Medium molecule — should always work."""
        self._check_mapping("c1ccccc1", "c1ccc(O)cc1", min_mcs=6)

    def test_penicillin_g(self):
        """Large molecule — regression guard."""
        self._check_mapping(
            "CC1(C)SC2C(NC(=O)Cc3ccccc3)C(=O)N2C1C(=O)O",  # penicillin G
            "CC1(C)SC2C(NC(=O)c3ccccc3)C(=O)N2C1C(=O)O",   # penicillin V-like
            min_mcs=15
        )

    def test_symmetric_benzene(self):
        """Symmetric molecule — index mapping must still be element-correct."""
        self._check_mapping("c1ccccc1", "c1ccccc1", min_mcs=6)


# ===========================================================================
# Edge Cases
# ===========================================================================

class TestEdgeCases:
    """Boundary-condition tests for SMSD Python bindings."""

    # 1. parse_smiles("") -> should handle gracefully (empty MolGraph)
    def test_parse_smiles_empty_string(self):
        """Parsing an empty SMILES string should return an empty MolGraph."""
        mol = parse_smiles("")
        assert mol.n == 0

    # 2. parse_smiles("[H]") -> single H atom
    def test_parse_smiles_single_hydrogen(self):
        """Parsing explicit hydrogen should yield a single-atom molecule."""
        mol = parse_smiles("[H]")
        assert mol.n == 1


class TestRGroupDecomposition:
    """Python R-group decomposition should work without direct neighbor-array access."""

    def test_toluene_on_benzene_core_finds_methyl(self):
        results = smsd.decompose_r_groups("c1ccccc1", ["Cc1ccccc1"])
        assert len(results) == 1
        assert results[0]["core"] == [1, 2, 3, 4, 5, 6]
        assert "R1" in results[0]
        assert len(results[0]["R1"]) == 1

    def test_phenol_on_benzene_core_finds_hydroxyl(self):
        results = smsd.decompose_r_groups("c1ccccc1", ["Oc1ccccc1"])
        assert len(results) == 1
        assert "R1" in results[0]
        assert len(results[0]["R1"]) == 1

    def test_para_disubstituted_benzene_finds_two_rgroups(self):
        results = smsd.decompose_r_groups("c1ccccc1", ["Nc1ccc(O)cc1"])
        assert len(results) == 1
        assert "R1" in results[0]
        assert "R2" in results[0]
        assert len(results[0]["R1"]) == 1
        assert len(results[0]["R2"]) == 1

    def test_attachment_point_atom_map_roundtrip(self):
        mol = parse_smiles("[*:1]CC[*:2]")
        assert mol.n == 4
        assert mol.atom_class[0] == 1
        assert mol.atom_class[3] == 2
        out = to_smiles(mol)
        assert "R1" in out and "R2" in out

    def test_patent_rgroup_placeholder_reads_as_attachment_point(self):
        mol = parse_smiles("[R1]c1ccccc1")
        assert mol.n == 7
        assert mol.atomic_num[0] == 0
        assert mol.atom_class[0] == 1
        assert to_smiles(mol).startswith("[R1]")

    def test_mol_block_metadata_roundtrip(self):
        mol_block = (
            "demo-name\n"
            "  SMSD  2D\n"
            "demo-comment\n"
            "  2  1  0  0  0  0  0  0  0  0999 V2000\n"
            "    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
            "    1.5000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n"
            "  1  2  1  0  0  0  0\n"
            "M  END\n"
            "> <ID>\n"
            "cmpd-7\n"
            "\n"
        )
        mol = read_mol_block(mol_block)
        assert mol.name == "demo-name"
        assert mol.program_line == "  SMSD  2D"
        assert mol.comment == "demo-comment"
        assert mol.properties["ID"] == "cmpd-7"

        out = write_mol_block(mol)
        assert out.startswith("demo-name\n")
        assert "> <ID>\ncmpd-7\n\n" in out
        assert write_sdf_record(mol).endswith("$$$$\n")

    def test_smarts_writer_roundtrip_query_features(self):
        mol = parse_smiles("[NH3+]C1=CC=CC=C1")
        smarts = to_smarts(
            mol,
            include_aromaticity=True,
            include_charge=True,
            include_ring_member=True,
            include_h_count=True,
        )
        assert smarts
        assert smsd.smarts_match(smarts, mol)

    def test_v3000_roundtrip_metadata_and_stereo(self):
        mol = parse_smiles("N[C@@H](C)C(=O)O")
        block = write_mol_block_v3000(mol)
        assert "V3000" in block
        reparsed = read_mol_block(block)
        assert reparsed.n == mol.n
        rs = smsd.assign_rs(reparsed)
        assert rs

    def test_v3000_roundtrip_attachment_labels(self):
        mol = parse_smiles("[CH3:7][R1]")
        block = write_mol_block_v3000(mol)
        reparsed = read_mol_block(block)
        assert reparsed.atom_class[0] == 7
        assert reparsed.atom_class[1] == 1

    def test_v2000_rgroup_roundtrip(self):
        mol = parse_smiles("[R1]c1ccccc1")
        block = write_mol_block(mol)
        assert "R#" in block
        assert "M  RGP" in block
        reparsed = read_mol_block(block)
        assert reparsed.atomic_num[0] == 0
        assert reparsed.atom_class[0] == 1

    def test_v2000_alias_rgroup_reads(self):
        block = (
            "alias-rgroup\n"
            "  SMSD\n"
            "\n"
            "  2  1  0  0  0  0  0  0  0  0999 V2000\n"
            "    0.0000    0.0000    0.0000 *   0  0  0  0  0  0  0  0  0  0  0  0\n"
            "    1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
            "  1  2  1  0  0  0  0\n"
            "A    1\n"
            "R2\n"
            "M  END\n"
        )
        mol = read_mol_block(block)
        assert mol.atomic_num[0] == 0
        assert mol.atom_class[0] == 2

    def test_v2000_multi_attachment_roundtrip(self):
        mol = parse_smiles("[R1]c1ccc([R2])cc1")
        reparsed = read_mol_block(write_mol_block(mol))
        labels = [x for x in reparsed.atom_class if x in (1, 2)]
        assert len(labels) == 2

    def test_v3000_double_bond_ez_roundtrip(self):
        mol = parse_smiles("C/C=C/C")
        reparsed = read_mol_block(write_mol_block_v3000(mol))
        ez = smsd.assign_ez(reparsed)
        assert ez[(1, 2)] == "E"

    def test_v3000_multiple_stereocentres_roundtrip(self):
        mol = parse_smiles("O[C@H]1[C@@H](O)[C@H](O)[C@@H](O)[C@H](CO)O1")
        reparsed = read_mol_block(write_mol_block_v3000(mol))
        rs = smsd.assign_rs(reparsed)
        assert len(rs) >= 4

    # 3. find_mcs on two empty-ish molecules
    def test_find_mcs_empty_molecules(self):
        """MCS of two empty MolGraphs should be an empty mapping."""
        empty = MolGraph()
        mapping = find_mcs(empty, empty)
        assert isinstance(mapping, dict)
        assert len(mapping) == 0

    def test_find_mcs_empty_vs_real(self):
        """MCS of an empty MolGraph vs a real molecule should be empty."""
        empty = MolGraph()
        benzene = parse_smiles("c1ccccc1")
        mapping = find_mcs(empty, benzene)
        assert isinstance(mapping, dict)
        assert len(mapping) == 0

    # 4. fingerprint with very small molecule (radius=0 equivalent)
    def test_fingerprint_single_atom(self):
        """Fingerprint of a single atom should produce at least one set bit."""
        mol = parse_smiles("C")
        fp = path_fingerprint(mol)
        assert isinstance(fp, list)
        assert len(fp) >= 1  # at least 1 bit set

    def test_fingerprint_empty_molecule(self):
        """Fingerprint of an empty molecule should be an empty bit list."""
        empty = MolGraph()
        fp = path_fingerprint(empty)
        assert isinstance(fp, list)
        assert len(fp) == 0

    # 5. tanimoto([], []) -> 0.0
    def test_tanimoto_empty_fingerprints(self):
        """Tanimoto of two empty fingerprints should be 0.0."""
        sim = smsd.tanimoto([], [])
        assert sim == pytest.approx(0.0), (
            f"tanimoto([], []) should be 0.0, got {sim}"
        )

    def test_tanimoto_one_empty(self):
        """Tanimoto where one fingerprint is empty should be 0.0."""
        fp = smsd.fingerprint("c1ccccc1")
        sim = smsd.tanimoto(fp, [])
        assert sim == pytest.approx(0.0), (
            f"tanimoto(fp, []) should be 0.0, got {sim}"
        )

    def test_tanimoto_identical(self):
        """Tanimoto of identical non-empty fingerprints should be 1.0."""
        fp = smsd.fingerprint("c1ccccc1")
        sim = smsd.tanimoto(fp, fp)
        assert sim == pytest.approx(1.0)

    # 6. from_rdkit(None) -> should raise TypeError
    def test_from_rdkit_none_raises(self):
        """Passing None to from_rdkit should raise TypeError."""
        try:
            from smsd import from_rdkit
        except ImportError:
            pytest.skip("from_rdkit not available")
        with pytest.raises(TypeError):
            from_rdkit(None)

    # 7. mcs_rdkit_native with identical molecules -> full mapping
    def test_mcs_rdkit_native_identical(self):
        """mcs_rdkit_native on identical molecules should return a full mapping."""
        try:
            from rdkit import Chem
            from smsd import mcs_rdkit_native
        except ImportError:
            pytest.skip("RDKit not installed")

        mol = Chem.MolFromSmiles("CCO")
        n = mol.GetNumAtoms()
        mapping = mcs_rdkit_native(mol, mol)
        assert len(mapping) == n, (
            f"Identical ethanol MCS should map all {n} atoms, got {len(mapping)}"
        )
        # Every mapped atom should have matching element
        for q_idx, t_idx in mapping.items():
            assert 0 <= q_idx < n
            assert 0 <= t_idx < n
            q_elem = mol.GetAtomWithIdx(q_idx).GetAtomicNum()
            t_elem = mol.GetAtomWithIdx(t_idx).GetAtomicNum()
            assert q_elem == t_elem, (
                f"Element mismatch at {q_idx}->{t_idx}: {q_elem} vs {t_elem}"
            )

    # 8. translate_mapping with empty dict -> empty dict
    def test_translate_mapping_empty(self):
        """translate_mapping with an empty mapping should return an empty dict."""
        from smsd import translate_mapping
        result = translate_mapping({})
        assert isinstance(result, dict)
        assert len(result) == 0

    def test_translate_mapping_empty_with_graphs(self):
        """translate_mapping with empty mapping and real MolGraphs should be empty."""
        try:
            from rdkit import Chem
            from smsd import from_rdkit, translate_mapping
        except ImportError:
            pytest.skip("RDKit not installed")

        mol1 = Chem.MolFromSmiles("CCO")
        mol2 = Chem.MolFromSmiles("CCCO")
        g1 = from_rdkit(mol1)
        g2 = from_rdkit(mol2)
        result = translate_mapping({}, g1, g2)
        assert isinstance(result, dict)
        assert len(result) == 0

    # 9. Batch operations with empty target list -> empty result
    def test_batch_mcs_empty_targets(self):
        """batch_mcs with an empty target list should return an empty list."""
        from smsd import batch_mcs
        query = parse_smiles("c1ccccc1")
        results = batch_mcs(query, [])
        assert isinstance(results, list)
        assert len(results) == 0

    def test_batch_substructure_empty_targets(self):
        """batch_substructure with an empty target list should return an empty list."""
        from smsd import batch_substructure
        query = parse_smiles("c1ccccc1")
        results = batch_substructure(query, [])
        assert isinstance(results, list)
        assert len(results) == 0

    def test_screen_targets_empty_library(self):
        """screen_targets with an empty library should return an empty list."""
        query = parse_smiles("c1ccccc1")
        hits = screen_targets(query, [], 0.5)
        assert isinstance(hits, list)
        assert len(hits) == 0

    # Additional edge cases for completeness

    def test_find_mcs_single_atom_self(self):
        """MCS of a single atom with itself should be size 1."""
        mol = parse_smiles("C")
        mapping = find_mcs(mol, mol)
        assert len(mapping) == 1

    def test_find_mcs_disjoint_elements(self):
        """MCS of all-C vs all-N disconnected atoms should be empty."""
        all_c = parse_smiles("C.C.C")
        all_n = parse_smiles("[NH3].[NH3].[NH3]")
        mapping = find_mcs(all_c, all_n)
        assert len(mapping) == 0, (
            f"Disjoint elements should give empty MCS, got {len(mapping)}"
        )

    def test_mcs_never_exceeds_smaller(self):
        """MCS size must never exceed the size of the smaller molecule."""
        pairs = [
            ("C", "CCCCCC"),
            ("CCO", "c1ccccc1"),
            ("C=C", "C=CC=CC=C"),
        ]
        for smi1, smi2 in pairs:
            g1 = parse_smiles(smi1)
            g2 = parse_smiles(smi2)
            mapping = find_mcs(g1, g2)
            min_size = min(g1.n, g2.n)
            assert len(mapping) <= min_size, (
                f"MCS({smi1}, {smi2}) = {len(mapping)} exceeds min size {min_size}"
            )


# ===========================================================================
# User-Feedback Tests (v6.2.0)
# ===========================================================================

class TestUserFeedbackV620:
    """Tests from user feedback requiring coverage in v6.2.0."""

    def test_from_rdkit_round_trip_ethanol(self):
        """from_rdkit(MolFromSmiles('OCC')) and from_rdkit(MolFromSmiles('CCO'))
        should both produce MolGraphs with 3 atoms (both are ethanol)."""
        try:
            from rdkit import Chem
            from smsd import from_rdkit
        except ImportError:
            pytest.skip("RDKit not installed")

        g1 = from_rdkit(Chem.MolFromSmiles("OCC"))
        g2 = from_rdkit(Chem.MolFromSmiles("CCO"))
        assert g1.n == 3, f"OCC should have 3 heavy atoms, got {g1.n}"
        assert g2.n == 3, f"CCO should have 3 heavy atoms, got {g2.n}"
        assert g1.n == g2.n, (
            f"OCC and CCO are both ethanol; atom counts should match: "
            f"{g1.n} vs {g2.n}"
        )

    def test_empty_smiles_returns_empty_molgraph(self):
        """parse_smiles('') should return a MolGraph with len(g)==0, not raise."""
        g = parse_smiles("")
        assert len(g) == 0, f"Empty SMILES should give 0 atoms, got {len(g)}"

    def test_circular_fingerprint_dimethyl_sulfone(self):
        """circular_fingerprint on dimethyl sulfone (multi-valent S) should not crash."""
        from smsd import circular_fingerprint
        g = parse_smiles("CS(=O)(=O)C")
        fp = circular_fingerprint(g, radius=2)
        assert isinstance(fp, list), "circular_fingerprint should return a list"

    def test_circular_fingerprint_phosphoric_acid(self):
        """circular_fingerprint on H3PO4 (multi-valent P) should not crash."""
        from smsd import circular_fingerprint
        g = parse_smiles("OP(=O)(O)O")
        fp = circular_fingerprint(g, radius=2)
        assert isinstance(fp, list), "circular_fingerprint should return a list"


# ===========================================================================
# CIP (Cahn-Ingold-Prelog) stereodescriptor assignment
# ===========================================================================

class TestCIPAssignment:
    """Test CIP R/S and E/Z assignment via the Python API."""

    def test_lalanine_is_s(self):
        """L-alanine N[C@@H](C)C(=O)O => stereocentre has S configuration."""
        result = smsd.assign_rs("N[C@@H](C)C(=O)O")
        assert len(result) >= 1, "L-alanine should have at least one stereocentre"
        assert 'S' in result.values(), f"L-alanine should be S, got {result}"

    def test_dalanine_is_r(self):
        """D-alanine N[C@H](C)C(=O)O => stereocentre has R configuration."""
        result = smsd.assign_rs("N[C@H](C)C(=O)O")
        assert len(result) >= 1, "D-alanine should have at least one stereocentre"
        assert 'R' in result.values(), f"D-alanine should be R, got {result}"

    def test_lcysteine_is_r(self):
        """L-cysteine N[C@@H](CS)C(=O)O => R (sulfur changes CIP priority)."""
        result = smsd.assign_rs("N[C@@H](CS)C(=O)O")
        assert len(result) >= 1, "L-cysteine should have at least one stereocentre"
        assert 'R' in result.values(), f"L-cysteine should be R, got {result}"

    def test_no_stereocentre(self):
        """Isobutane CC(C)C has no stereocentres => empty dict."""
        result = smsd.assign_rs("CC(C)C")
        assert result == {}, f"Isobutane should have no stereocentres, got {result}"

    def test_symmetric_not_stereo(self):
        """CF2Cl2 C(F)(F)(Cl)Cl has two identical pairs => not a stereocentre."""
        result = smsd.assign_rs("C(F)(F)(Cl)Cl")
        assert result == {}, f"CF2Cl2 should have no stereocentres, got {result}"

    def test_implicit_h_lowest(self):
        """[C@@H](F)(Cl)Br => H is lowest priority, should assign R or S."""
        result = smsd.assign_rs("[C@@H](F)(Cl)Br")
        assert len(result) == 1, f"CHFClBr should have exactly 1 stereocentre, got {result}"

    def test_enantiomers_opposite(self):
        """@@ and @ on same molecule should give opposite R/S."""
        rs1 = smsd.assign_rs("[C@@H](F)(Cl)Br")
        rs2 = smsd.assign_rs("[C@H](F)(Cl)Br")
        assert len(rs1) == 1 and len(rs2) == 1
        v1 = list(rs1.values())[0]
        v2 = list(rs2.values())[0]
        assert v1 != v2, f"Enantiomers should have opposite R/S, got {v1} and {v2}"

    def test_e_butene(self):
        """(E)-2-butene C/C=C/C => E."""
        result = smsd.assign_ez("C/C=C/C")
        assert len(result) >= 1, "(E)-2-butene should have an E/Z bond"
        assert 'E' in result.values(), f"C/C=C/C should be E, got {result}"

    def test_z_butene(self):
        r"""(Z)-2-butene C/C=C\C => Z."""
        result = smsd.assign_ez(r"C/C=C\C")
        assert len(result) >= 1, "(Z)-2-butene should have an E/Z bond"
        assert 'Z' in result.values(), f"C/C=C\\C should be Z, got {result}"

    def test_ethylene_no_stereo(self):
        """Ethylene C=C without / or \\ => empty E/Z map."""
        result = smsd.assign_ez("C=C")
        assert result == {}, f"Ethylene should have no E/Z, got {result}"

    def test_assign_cip_combined(self):
        """assign_cip returns CIPDescriptors with both R/S and E/Z."""
        desc = smsd.assign_cip("N[C@@H](C)C(=O)O")
        assert hasattr(desc, 'rs_labels'), "CIPDescriptors should have rs_labels"
        assert hasattr(desc, 'ez_bonds'), "CIPDescriptors should have ez_bonds"
        # Should have at least one non-NONE rs_label
        from smsd import RSLabel
        n_stereo = sum(1 for l in desc.rs_labels if l != RSLabel.NONE)
        assert n_stereo >= 1, "L-alanine should have at least 1 stereocentre in CIPDescriptors"

    def test_assign_rs_returns_dict(self):
        """assign_rs returns a plain dict, not a CIPDescriptors object."""
        result = smsd.assign_rs("N[C@@H](C)C(=O)O")
        assert isinstance(result, dict), f"assign_rs should return dict, got {type(result)}"

    def test_assign_ez_returns_dict_with_tuples(self):
        """assign_ez returns dict with tuple keys."""
        result = smsd.assign_ez("C/C=C/C")
        assert isinstance(result, dict)
        for key in result:
            assert isinstance(key, tuple), f"E/Z key should be tuple, got {type(key)}"
            assert len(key) == 2, f"E/Z key should have 2 elements, got {len(key)}"

    def test_deep_priority_resolution(self):
        """Molecules needing depth > 1 to resolve CIP priorities."""
        # [C@@](Cl)(C)(CC)F vs [C@](Cl)(C)(CC)F
        rs1 = smsd.assign_rs("[C@@](Cl)(C)(CC)F")
        rs2 = smsd.assign_rs("[C@](Cl)(C)(CC)F")
        assert len(rs1) >= 1, "Should have a stereocentre"
        assert len(rs2) >= 1, "Should have a stereocentre"
        v1 = list(rs1.values())[0]
        v2 = list(rs2.values())[0]
        assert v1 != v2, f"@@ and @ enantiomers should differ, got {v1} and {v2}"

    def test_cholesterol_multiple_stereo(self):
        """Cholesterol should have multiple stereocentres."""
        smi = ("C([C@@H]1CC2=CC(=O)CC[C@@]2(C)[C@H]1"
               "[C@@H]1CC[C@H]([C@@H](CCCC(C)C)C)[C@@]1(C)CC)O")
        result = smsd.assign_rs(smi)
        assert len(result) >= 2, (
            f"Cholesterol should have >= 2 stereocentres, got {len(result)}"
        )

    def test_assign_rs_accepts_molgraph(self):
        """assign_rs should accept MolGraph objects, not just SMILES strings."""
        g = parse_smiles("N[C@@H](C)C(=O)O")
        result = smsd.assign_rs(g)
        assert len(result) >= 1


# ===========================================================================
# Timeout enforcement -- large molecules must never hang
# ===========================================================================

class TestTimeoutEnforcement:
    """MCS on large molecules must respect timeout_ms and return promptly."""

    # Coenzyme A (66 heavy atoms)
    COA_SMILES = (
        "CC(C)(COP(=O)(O)OP(=O)(O)OC[C@@H]1[C@H]([C@H]([C@@H](O1)"
        "n2cnc3c(N)ncnc23)O)OP(=O)(O)O)[C@@H](O)C(=O)NCCC(=O)NCCSC"
    )

    # ATP (31 heavy atoms)
    ATP_SMILES = (
        "C1=NC(=C2C(=N1)N(C=N2)[C@H]3[C@@H]([C@@H]([C@H](O3)"
        "COP(=O)(O)OP(=O)(O)OP(=O)(O)O)O)O)N"
    )

    # Paclitaxel (Taxol, 71 heavy atoms)
    TAXOL_SMILES = (
        "CC1=C2C(C(=O)C3(C(CC4C(C3C(C(C2(C)C)(CC1OC(=O)C(C(C5=CC=CC=C5)"
        "NC(=O)C6=CC=CC=C6)O)O)OC(=O)C7=CC=CC=C7)(CO4)OC(=O)C)O)C)OC(=O)C"
    )

    def test_coa_self_match_timeout(self):
        """CoA self-match with 2s timeout must return within 4s."""
        import time
        coa = parse_smiles(self.COA_SMILES)
        t0 = time.monotonic()
        mapping = find_mcs(coa, coa, ChemOptions(), McsOptions(), timeout_ms=2000)
        elapsed = time.monotonic() - t0
        assert elapsed < 4.0, (
            f"CoA self-match took {elapsed:.1f}s, expected < 4.0s (timeout_ms=2000)"
        )
        # Self-match should still produce a valid (possibly partial) mapping
        assert len(mapping) > 0, "CoA self-match should return a non-empty mapping"

    def test_atp_coa_cross_timeout(self):
        """ATP vs CoA with 2s timeout must return within 4s."""
        import time
        atp = parse_smiles(self.ATP_SMILES)
        coa = parse_smiles(self.COA_SMILES)
        t0 = time.monotonic()
        mapping = find_mcs(atp, coa, ChemOptions(), McsOptions(), timeout_ms=2000)
        elapsed = time.monotonic() - t0
        assert elapsed < 4.0, (
            f"ATP vs CoA MCS took {elapsed:.1f}s, expected < 4.0s (timeout_ms=2000)"
        )
        assert len(mapping) > 0, "ATP vs CoA should have a non-empty MCS"

    def test_taxol_self_match_timeout(self):
        """Paclitaxel (71 atoms) self-match with 2s timeout must return within 4s."""
        import time
        taxol = parse_smiles(self.TAXOL_SMILES)
        t0 = time.monotonic()
        mapping = find_mcs(taxol, taxol, ChemOptions(), McsOptions(), timeout_ms=2000)
        elapsed = time.monotonic() - t0
        assert elapsed < 4.0, (
            f"Taxol self-match took {elapsed:.1f}s, expected < 4.0s (timeout_ms=2000)"
        )
        assert len(mapping) > 0, "Taxol self-match should return a non-empty mapping"

    def test_taxol_coa_cross_timeout(self):
        """Taxol vs CoA with 2s timeout must return within 4s."""
        import time
        taxol = parse_smiles(self.TAXOL_SMILES)
        coa = parse_smiles(self.COA_SMILES)
        t0 = time.monotonic()
        mapping = find_mcs(taxol, coa, ChemOptions(), McsOptions(), timeout_ms=2000)
        elapsed = time.monotonic() - t0
        assert elapsed < 4.0, (
            f"Taxol vs CoA MCS took {elapsed:.1f}s, expected < 4.0s (timeout_ms=2000)"
        )

    def test_python_mcs_passes_timeout(self):
        """The high-level smsd.mcs() wrapper must respect timeout_ms."""
        import time
        t0 = time.monotonic()
        mapping = smsd.mcs(self.COA_SMILES, self.TAXOL_SMILES, timeout_ms=1000)
        elapsed = time.monotonic() - t0
        assert elapsed < 3.0, (
            f"smsd.mcs() CoA vs Taxol took {elapsed:.1f}s, expected < 3.0s"
        )


# ===========================================================================
# API Completeness Tests (v6.4.0 QA sweep)
# ===========================================================================

class TestAPICompleteness:
    """Tests for previously untested public APIs."""

    def test_to_hex_from_hex_roundtrip(self):
        """to_hex -> from_hex round-trip must preserve fingerprint bits."""
        fp = path_fingerprint(parse_smiles("c1ccccc1"), 7, 2048)
        hex_str = smsd.to_hex(fp, fp_size=2048)
        assert isinstance(hex_str, str)
        assert len(hex_str) > 0
        fp_back = smsd.from_hex(hex_str)
        assert sorted(fp_back) == sorted(fp), (
            "to_hex -> from_hex round-trip must preserve all set bit positions"
        )

    def test_to_binary_string(self):
        """to_binary_string must return a string of 0s and 1s of correct length."""
        fp = path_fingerprint(parse_smiles("CCO"), 7, 1024)
        bits = smsd.to_binary_string(fp, fp_size=1024)
        assert len(bits) == 1024, f"Expected 1024 chars, got {len(bits)}"
        assert all(c in ('0', '1') for c in bits), "Must contain only 0 and 1"
        # Verify set bits match
        set_positions = [i for i, c in enumerate(bits) if c == '1']
        assert sorted(set_positions) == sorted(fp), (
            "Binary string set-bit positions must match fingerprint"
        )

    def test_counts_to_array(self):
        """counts_to_array must convert sparse dict to dense list."""
        counts = {42: 3, 1000: 1, 0: 7}
        arr = smsd.counts_to_array(counts, 2048)
        assert len(arr) == 2048
        assert arr[0] == 7
        assert arr[42] == 3
        assert arr[1000] == 1
        assert arr[500] == 0
        # Out-of-range keys should be silently ignored
        counts[5000] = 99
        arr2 = smsd.counts_to_array(counts, 2048)
        assert len(arr2) == 2048

    def test_compute_sssr_benzene(self):
        """computeSSSR on benzene should return 1 ring of size 6."""
        rings = smsd.compute_sssr("c1ccccc1")
        assert len(rings) == 1, f"Benzene should have 1 SSSR ring, got {len(rings)}"
        assert len(rings[0]) == 6, f"Benzene ring should have 6 atoms, got {len(rings[0])}"

    def test_layout_sssr_naphthalene(self):
        """layoutSSSR on naphthalene should return 2 rings of size 6, sharing 2 atoms."""
        rings = smsd.layout_sssr("c1ccc2ccccc2c1")
        assert len(rings) == 2, f"Naphthalene should have 2 SSSR rings, got {len(rings)}"
        assert len(rings[0]) == 6, f"First ring should have 6 atoms, got {len(rings[0])}"
        assert len(rings[1]) == 6, f"Second ring should have 6 atoms, got {len(rings[1])}"
        # Fused rings share 2 atoms (one bond)
        shared = set(rings[0]) & set(rings[1])
        assert len(shared) == 2, (
            f"Naphthalene fused rings should share 2 atoms, got {len(shared)}"
        )

    def test_mcs_from_smiles(self):
        """mcs_from_smiles should return a MatchResult with valid fields."""
        result = smsd.mcs_from_smiles("c1ccccc1", "Cc1ccccc1")
        assert result.size >= 6, f"Benzene/toluene MCS should be >= 6, got {result.size}"
        assert result.tanimoto > 0.0, "Tanimoto should be positive"
        assert len(result.mapping) >= 6, "Mapping should have >= 6 entries"
        assert isinstance(result.mcs_smiles, str), "mcs_smiles should be a string"
        assert len(result.mcs_smiles) > 0, "mcs_smiles should not be empty"


# ===========================================================================
# Reaction-Aware MCS Post-Filter (v6.4.0)
# ===========================================================================

class TestReactionAwareMCS:
    """Tests for the reaction-aware MCS post-filter that re-ranks candidates
    by heteroatom coverage, rare-element importance, and connectivity."""

    SAM = "C[S+](CCC(N)C(=O)O)CC1OC(n2cnc3c(N)ncnc32)C(O)C1O"
    HOMOCYSTEINE = "CSCC(N)C(=O)O"
    ATP = "c1nc(N)c2ncn(C3OC(COP(=O)(O)OP(=O)(O)OP(=O)(O)O)C(O)C3O)c2n1"
    ADP = "c1nc(N)c2ncn(C3OC(COP(=O)(O)OP(=O)(O)O)C(O)C3O)c2n1"
    ACETYL_COA = (
        "CC(=O)SCC(NC(=O)CCNC(=O)C(O)C(C)(C)COP(=O)(O)OP(=O)(O)OCC1OC"
        "(n2cnc3c(N)ncnc32)C(O)C1OP(=O)(O)O)C(=O)O"
    )
    COA = (
        "OSCC(NC(=O)CCNC(=O)C(O)C(C)(C)COP(=O)(O)OP(=O)(O)OCC1OC"
        "(n2cnc3c(N)ncnc32)C(O)C1OP(=O)(O)O)C(=O)O"
    )

    @staticmethod
    def _find_mapped_element(mapping, mol, atomic_num):
        """Return index of first mapped atom in mol with given atomic number, or -1."""
        for qi in mapping:
            if mol.atomic_num[qi] == atomic_num:
                return qi
        return -1

    def test_sam_homocysteine_reaction_aware_includes_sulfur(self):
        """Reaction-aware SAM vs Homocysteine MUST include the S atom (Z=16)."""
        mapping = smsd.map_reaction_aware(self.SAM, self.HOMOCYSTEINE, timeout_ms=10000)
        assert len(mapping) > 0, "Reaction-aware SAM vs Homocysteine should produce a mapping"

        g1 = parse_smiles(self.SAM)
        s_idx = self._find_mapped_element(mapping, g1, 16)
        assert s_idx >= 0, (
            "Reaction-aware SAM->Homocysteine must include the S atom (Z=16) "
            "in the mapping"
        )
        assert g1.atomic_num[s_idx] == 16, "Mapped S atom must have atomic_num == 16"

    def test_sam_homocysteine_standard_no_assertion_on_sulfur(self):
        """Standard MCS for SAM vs Homocysteine (no assertion on S inclusion)."""
        g1 = parse_smiles(self.SAM)
        g2 = parse_smiles(self.HOMOCYSTEINE)
        chem = ChemOptions()
        chem.match_formal_charge = False
        mapping = find_mcs(g1, g2, chem, McsOptions())
        assert len(mapping) >= 3, (
            f"Standard MCS SAM vs Homocysteine should map >= 3 atoms, got {len(mapping)}"
        )

    def test_atp_adp_reaction_aware_preserves_phosphate_oxygens(self):
        """Reaction-aware ATP vs ADP should map phosphorus and most oxygens."""
        mapping = smsd.map_reaction_aware(self.ATP, self.ADP, timeout_ms=10000)
        assert len(mapping) > 0, "ATP vs ADP reaction-aware should produce a mapping"

        g1 = parse_smiles(self.ATP)

        # Count mapped oxygens (Z=8)
        o_mapped = sum(1 for qi in mapping if g1.atomic_num[qi] == 8)
        assert o_mapped >= 5, (
            f"Reaction-aware ATP->ADP should map >= 5 oxygen atoms, got {o_mapped}"
        )

        # Phosphorus (Z=15) must be mapped
        p_idx = self._find_mapped_element(mapping, g1, 15)
        assert p_idx >= 0, (
            "Reaction-aware ATP->ADP must include at least one P atom (Z=15)"
        )

    def test_acetyl_coa_coa_reaction_aware_preserves_s_c_bond(self):
        """Reaction-aware Acetyl-CoA vs CoA should preserve S atom and S-C bond."""
        mapping = smsd.map_reaction_aware(
            self.ACETYL_COA, self.COA, timeout_ms=10000
        )
        assert len(mapping) > 0, (
            "Acetyl-CoA vs CoA reaction-aware should produce a mapping"
        )

        g1 = parse_smiles(self.ACETYL_COA)

        # S atom (Z=16) must be mapped
        s_idx = self._find_mapped_element(mapping, g1, 16)
        assert s_idx >= 0, (
            "Reaction-aware Acetyl-CoA->CoA must include the S atom (Z=16)"
        )
        assert g1.atomic_num[s_idx] == 16, "Mapped S atom must have atomic_num == 16"

    def test_map_reaction_aware_accepts_smiles_strings(self):
        """map_reaction_aware should accept SMILES strings directly."""
        mapping = smsd.map_reaction_aware("CCO", "CCO", timeout_ms=5000)
        assert len(mapping) == 3, (
            f"Ethanol self-match via map_reaction_aware should map 3 atoms, got {len(mapping)}"
        )

    def test_map_reaction_aware_accepts_molgraph(self):
        """map_reaction_aware should accept MolGraph objects."""
        g1 = parse_smiles("CCO")
        g2 = parse_smiles("CCO")
        mapping = smsd.map_reaction_aware(g1, g2, timeout_ms=5000)
        assert len(mapping) == 3, (
            f"Ethanol self-match via MolGraph should map 3 atoms, got {len(mapping)}"
        )
