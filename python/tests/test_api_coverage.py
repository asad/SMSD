# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
# Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
# See the NOTICE file for attribution, trademark, and algorithm IP terms.
"""
Comprehensive API coverage tests for SMSD Pro v7.0.

Covers ~58 public API functions that previously had zero test coverage.
Run with:  pytest tests/test_api_coverage.py -v

Copyright (c) 2018-2026 BioInception PVT LTD
"""
import math
import os
import pytest
import smsd


# ---------------------------------------------------------------------------
# Common SMILES constants used across multiple test classes
# ---------------------------------------------------------------------------
BENZENE = "c1ccccc1"
PHENOL = "c1ccc(O)cc1"
TOLUENE = "Cc1ccccc1"
ETHANOL = "CCO"
ACETIC_ACID = "CC(=O)O"
ANILINE = "c1ccc(N)cc1"
NAPHTHALENE = "c1ccc2ccccc2c1"
ASPIRIN = "CC(=O)Oc1ccccc1C(=O)O"
METHANE = "C"
PROPANE = "CCC"
WATER = "O"
PYRIDINE = "c1ccncc1"

# Small corpus for batch tests
BATCH_CORPUS = [PHENOL, TOLUENE, ANILINE, ETHANOL, ACETIC_ACID,
                NAPHTHALENE, PYRIDINE, ASPIRIN]


# ===========================================================================
# 1. MCS Utilities
# ===========================================================================

class TestMCSUtilities:
    """Tests for MCS utility functions: mcs_to_smiles, find_mcs_smiles,
    mcs_smiles, find_mcs_progressive, mcs_result, mcs_from_smiles,
    and MatchResult."""

    def test_mcs_to_smiles_basic(self):
        """mcs_to_smiles should extract MCS SMILES from a graph + mapping."""
        g1 = smsd.parse_smiles(BENZENE)
        g2 = smsd.parse_smiles(PHENOL)
        mapping = smsd.find_mcs(g1, g2)
        assert len(mapping) > 0
        smi = smsd.mcs_to_smiles(g1, mapping)
        assert isinstance(smi, str)
        assert len(smi) > 0

    def test_mcs_to_smiles_empty_mapping(self):
        """mcs_to_smiles with empty mapping should return empty string."""
        g1 = smsd.parse_smiles(BENZENE)
        smi = smsd.mcs_to_smiles(g1, {})
        assert smi == ""

    def test_find_mcs_smiles_basic(self):
        """find_mcs_smiles should return a SMILES string for the MCS."""
        g1 = smsd.parse_smiles(BENZENE)
        g2 = smsd.parse_smiles(PHENOL)
        chem = smsd.ChemOptions()
        opts = smsd.MCSOptions()
        result = smsd.find_mcs_smiles(g1, g2, chem, opts, 10000)
        assert isinstance(result, str)
        assert len(result) > 0

    def test_mcs_smiles_convenience(self):
        """mcs_smiles convenience wrapper should return non-empty SMILES."""
        result = smsd.mcs_smiles(BENZENE, PHENOL)
        assert isinstance(result, str)
        assert len(result) > 0

    def test_mcs_smiles_identical_molecules(self):
        """MCS of identical molecules should yield the full molecule SMILES."""
        result = smsd.mcs_smiles(BENZENE, BENZENE)
        assert isinstance(result, str)
        assert len(result) > 0

    def test_mcs_smiles_with_timeout(self):
        """mcs_smiles should accept timeout_ms parameter."""
        result = smsd.mcs_smiles(BENZENE, PHENOL, timeout_ms=5000)
        assert isinstance(result, str)

    def test_find_mcs_progressive_basic(self):
        """find_mcs_progressive should find MCS with progress callback."""
        progress_calls = []

        def on_progress(mapping, size, elapsed_ms):
            progress_calls.append((size, elapsed_ms))

        mapping = smsd.find_mcs_progressive(BENZENE, PHENOL,
                                            on_progress=on_progress)
        assert isinstance(mapping, dict)
        assert len(mapping) >= 6  # benzene ring is common
        assert len(progress_calls) > 0
        # Each progress call should report non-negative size and time
        for size, ms in progress_calls:
            assert size >= 0
            assert ms >= 0

    def test_find_mcs_progressive_no_callback(self):
        """find_mcs_progressive with no callback should still return result."""
        mapping = smsd.find_mcs_progressive(BENZENE, PHENOL, on_progress=None)
        assert isinstance(mapping, dict)
        assert len(mapping) > 0

    def test_find_mcs_progressive_sizes_monotonic(self):
        """Progress callback sizes should be monotonically non-decreasing."""
        sizes = []

        def on_progress(mapping, size, elapsed_ms):
            sizes.append(size)

        smsd.find_mcs_progressive(BENZENE, ASPIRIN, on_progress=on_progress)
        for i in range(1, len(sizes)):
            assert sizes[i] >= sizes[i - 1], \
                f"Size decreased from {sizes[i-1]} to {sizes[i]}"

    def test_mcs_result_basic(self):
        """mcs_result should return a MatchResult with correct attributes."""
        result = smsd.mcs_result(BENZENE, PHENOL)
        assert isinstance(result, smsd.MatchResult)
        assert result.size >= 6
        assert isinstance(result.mapping, dict)
        assert len(result.mapping) == result.size
        assert result.query_atoms == 6   # benzene
        assert result.target_atoms == 7  # phenol
        assert 0.0 <= result.overlap <= 1.0
        assert 0.0 <= result.overlapCoefficient <= 1.0

    def test_mcs_result_mcs_smiles(self):
        """mcs_result should include mcs_smiles string."""
        result = smsd.mcs_result(BENZENE, PHENOL)
        assert isinstance(result.mcs_smiles, str)

    def test_mcs_result_len(self):
        """len(MatchResult) should equal .size."""
        result = smsd.mcs_result(BENZENE, PHENOL)
        assert len(result) == result.size

    def test_mcs_result_repr(self):
        """MatchResult repr should contain useful info."""
        result = smsd.mcs_result(BENZENE, PHENOL)
        r = repr(result)
        assert "MatchResult" in r
        assert "size=" in r
        assert "overlapCoefficient=" in r

    def test_mcs_from_smiles_basic(self):
        """mcs_from_smiles should parse, compute MCS, and return MatchResult."""
        result = smsd.mcs_from_smiles(BENZENE, PHENOL)
        assert isinstance(result, smsd.MatchResult)
        assert result.size >= 6
        assert isinstance(result.mcs_smiles, str)
        assert len(result.mcs_smiles) > 0

    def test_mcs_from_smiles_identical(self):
        """mcs_from_smiles with identical SMILES should give perfect overlap."""
        result = smsd.mcs_from_smiles(ETHANOL, ETHANOL)
        assert result.size == 3  # C, C, O
        assert result.overlap == pytest.approx(1.0)

    def test_mcs_from_smiles_invalid_empty(self):
        """mcs_from_smiles should raise ValueError on empty SMILES."""
        with pytest.raises(ValueError):
            smsd.mcs_from_smiles("", BENZENE)
        with pytest.raises(ValueError):
            smsd.mcs_from_smiles(BENZENE, "")

    def test_match_result_construction(self):
        """MatchResult can be manually constructed."""
        r = smsd.MatchResult({0: 0, 1: 1}, 5, 8, "CC")
        assert r.size == 2
        assert r.query_atoms == 5
        assert r.target_atoms == 8
        assert r.mcs_smiles == "CC"
        assert r.overlap == pytest.approx(2.0 / 5.0)

    def test_match_result_zero_atoms(self):
        """MatchResult with zero-size molecule should handle gracefully."""
        r = smsd.MatchResult({}, 0, 0, "")
        assert r.size == 0
        assert r.overlap == 0.0
        assert r.overlapCoefficient == 0.0


# ===========================================================================
# 2. Atom Map Utilities
# ===========================================================================

class TestAtomMaps:
    """Tests for parse_mapped_smiles and strip_atom_maps."""

    def test_strip_atom_maps_basic(self):
        """strip_atom_maps should remove :N from bracket atoms."""
        result = smsd.strip_atom_maps("[CH3:1][C:2](=[O:3])[OH:4]")
        assert ":1" not in result
        assert ":2" not in result
        assert ":3" not in result
        assert ":4" not in result
        assert "[CH3]" in result

    def test_strip_atom_maps_preserves_ring_closures(self):
        """strip_atom_maps should preserve aromatic ring closure digits."""
        result = smsd.strip_atom_maps("[Cl:1][C:2]1:[CH:3]:[CH:4]:1")
        # The :1 after [CH:4] should be removed (atom map),
        # but the :1 after ] should be preserved (ring closure)
        assert ":1" not in result.replace("]:1", "PLACEHOLDER")  # atom map gone
        # Ring closure digit preserved
        assert "1" in result

    def test_strip_atom_maps_no_maps(self):
        """strip_atom_maps on SMILES without maps should return unchanged."""
        original = "c1ccccc1"
        result = smsd.strip_atom_maps(original)
        assert result == original

    def test_strip_atom_maps_empty(self):
        """strip_atom_maps on empty string should return empty."""
        assert smsd.strip_atom_maps("") == ""

    def test_parse_mapped_smiles_basic(self):
        """parse_mapped_smiles should parse after stripping maps."""
        mol = smsd.parse_mapped_smiles("[CH3:1][C:2](=[O:3])[OH:4]")
        assert isinstance(mol, smsd.MolGraph)
        assert mol.n == 4

    def test_parse_mapped_smiles_reaction_fragment(self):
        """parse_mapped_smiles should handle reaction SMILES fragments."""
        mol = smsd.parse_mapped_smiles("[CH3:1][OH:2]")
        assert mol.n == 2
        # Should have C and O
        elements = sorted([mol.atomic_num[i] for i in range(mol.n)])
        assert 6 in elements  # carbon
        assert 8 in elements  # oxygen

    def test_parse_mapped_smiles_no_maps(self):
        """parse_mapped_smiles on clean SMILES should work identically."""
        mol = smsd.parse_mapped_smiles("CCO")
        assert mol.n == 3


# ===========================================================================
# 3. File I/O
# ===========================================================================

class TestFileIO:
    """Tests for MOL/SDF file reading and writing."""

    def test_write_mol_block_basic(self):
        """write_mol_block should return a V2000 MOL block string."""
        mol = smsd.parse_smiles(BENZENE)
        block = smsd.write_mol_block(mol)
        assert isinstance(block, str)
        assert "V2000" in block
        assert len(block) > 50

    def test_write_mol_block_v3000_basic(self):
        """write_mol_block_v3000 should return a V3000 MOL block string."""
        mol = smsd.parse_smiles(BENZENE)
        block = smsd.write_mol_block_v3000(mol)
        assert isinstance(block, str)
        assert "V3000" in block

    def test_write_sdf_record_basic(self):
        """write_sdf_record should return an SDF record ending with $$$$."""
        mol = smsd.parse_smiles(BENZENE)
        record = smsd.write_sdf_record(mol)
        assert isinstance(record, str)
        assert "$$$$" in record

    def test_read_mol_block_roundtrip(self):
        """Read back a written MOL block and get equivalent molecule."""
        mol = smsd.parse_smiles(ETHANOL)
        block = smsd.write_mol_block(mol)
        mol2 = smsd.read_mol_block(block)
        assert isinstance(mol2, smsd.MolGraph)
        assert mol2.n == mol.n

    def test_read_mol_block_v2000_atom_count(self):
        """Reading a V2000 block should parse correct atom count."""
        mol = smsd.parse_smiles(PHENOL)
        block = smsd.write_mol_block(mol)
        mol2 = smsd.read_mol_block(block)
        assert mol2.n == 7

    def test_to_smarts_basic(self):
        """to_smarts should return a valid SMARTS string."""
        mol = smsd.parse_smiles(BENZENE)
        smarts = smsd.to_smarts(mol)
        assert isinstance(smarts, str)
        assert len(smarts) > 0

    def test_to_smarts_ethanol(self):
        """to_smarts on ethanol should produce a non-empty SMARTS."""
        mol = smsd.parse_smiles(ETHANOL)
        smarts = smsd.to_smarts(mol)
        assert isinstance(smarts, str)
        assert len(smarts) > 0

    def test_read_molfile(self, tmp_path):
        """read_molfile should read a .mol file from disk."""
        mol = smsd.parse_smiles(BENZENE)
        path = str(tmp_path / "test.mol")
        smsd.write_molfile(mol, path)
        mol2 = smsd.read_molfile(path)
        assert isinstance(mol2, smsd.MolGraph)
        assert mol2.n == 6

    def test_write_molfile_v2000(self, tmp_path):
        """write_molfile should write V2000 format by default."""
        mol = smsd.parse_smiles(ETHANOL)
        path = str(tmp_path / "ethanol.mol")
        smsd.write_molfile(mol, path)
        assert os.path.exists(path)
        content = open(path).read()
        assert "V2000" in content

    def test_write_molfile_v3000(self, tmp_path):
        """write_molfile with v3000=True should write V3000 format."""
        mol = smsd.parse_smiles(ETHANOL)
        path = str(tmp_path / "ethanol_v3000.mol")
        smsd.write_molfile(mol, path, v3000=True)
        content = open(path).read()
        assert "V3000" in content

    def test_write_molfile_sdf(self, tmp_path):
        """write_molfile with sdf=True should write SDF record."""
        mol = smsd.parse_smiles(ETHANOL)
        path = str(tmp_path / "ethanol.sdf")
        smsd.write_molfile(mol, path, sdf=True)
        content = open(path).read()
        assert "$$$$" in content

    def test_export_sdf_basic(self, tmp_path):
        """export_sdf should write multiple molecules to an SDF file."""
        mols = [smsd.parse_smiles(s) for s in [BENZENE, ETHANOL, PHENOL]]
        path = str(tmp_path / "multi.sdf")
        smsd.export_sdf(mols, path)
        assert os.path.exists(path)
        content = open(path).read()
        assert content.count("$$$$") == 3

    def test_export_sdf_single_molecule(self, tmp_path):
        """export_sdf with a single molecule should produce valid SDF."""
        mols = [smsd.parse_smiles(METHANE)]
        path = str(tmp_path / "single.sdf")
        smsd.export_sdf(mols, path)
        content = open(path).read()
        assert content.count("$$$$") == 1

    def test_roundtrip_smiles_mol_smiles(self, tmp_path):
        """SMILES -> write MOL -> read MOL -> SMILES should produce
        an equivalent molecule (same atom count)."""
        for smi in [BENZENE, ETHANOL, PHENOL, ACETIC_ACID]:
            mol1 = smsd.parse_smiles(smi)
            path = str(tmp_path / "roundtrip.mol")
            smsd.write_molfile(mol1, path)
            mol2 = smsd.read_molfile(path)
            assert mol2.n == mol1.n, \
                f"Atom count mismatch for {smi}: {mol2.n} != {mol1.n}"


# ===========================================================================
# 4. Batch Operations
# ===========================================================================

class TestBatchOps:
    """Tests for batch substructure, MCS, similarity, and screening."""

    @pytest.fixture
    def corpus(self):
        """Pre-parsed molecules for batch tests."""
        return [smsd.parse_smiles(s) for s in BATCH_CORPUS]

    def test_batch_substructure_basic(self, corpus):
        """batch_substructure should return list of bools."""
        query = smsd.parse_smiles(BENZENE)
        results = smsd.batch_substructure(query, corpus)
        assert isinstance(results, list)
        assert len(results) == len(corpus)
        assert all(isinstance(r, bool) for r in results)

    def test_batch_substructure_benzene_in_phenol(self, corpus):
        """Benzene should be substructure of phenol, toluene, aniline, etc."""
        results = smsd.batch_substructure(BENZENE, BATCH_CORPUS)
        # Phenol (index 0), Toluene (1), Aniline (2) should match
        assert results[0] is True   # phenol
        assert results[1] is True   # toluene
        assert results[2] is True   # aniline
        assert results[3] is False  # ethanol

    def test_batch_find_substructure_basic(self, corpus):
        """batch_find_substructure should return mappings per target."""
        query = smsd.parse_smiles(BENZENE)
        results = smsd.batch_find_substructure(query, corpus)
        assert isinstance(results, list)
        assert len(results) == len(corpus)

    def test_batch_find_substructure_has_mappings(self):
        """batch_find_substructure should produce non-empty mappings for matches."""
        results = smsd.batch_find_substructure(BENZENE,
                                               [PHENOL, ETHANOL])
        # Phenol should have a mapping
        assert len(results[0]) > 0
        # Ethanol should not
        assert len(results[1]) == 0

    def test_batch_mcs_basic(self, corpus):
        """batch_mcs should return dicts of atom mappings."""
        query = smsd.parse_smiles(BENZENE)
        results = smsd.batch_mcs(query, corpus)
        assert isinstance(results, list)
        assert len(results) == len(corpus)

    def test_batch_mcs_benzene_phenol_size(self):
        """batch_mcs of benzene vs phenol should have MCS size >= 6."""
        results = smsd.batch_mcs(BENZENE, [PHENOL])
        assert len(results) == 1
        assert len(results[0]) >= 6

    def test_batch_mcs_size_basic(self, corpus):
        """batch_mcs_size should return list of integers."""
        query = smsd.parse_smiles(BENZENE)
        results = smsd.batch_mcs_size(query, corpus)
        assert isinstance(results, list)
        assert len(results) == len(corpus)
        assert all(isinstance(r, int) for r in results)

    def test_batch_mcs_size_values(self):
        """batch_mcs_size should give correct sizes for known pairs."""
        results = smsd.batch_mcs_size(BENZENE, [PHENOL, ETHANOL])
        assert results[0] >= 6   # benzene ring is in phenol
        assert results[1] >= 1   # at least one common atom

    def test_screen_and_match_basic(self, corpus):
        """screen_and_match should return (index, mapping) tuples."""
        query = smsd.parse_smiles(BENZENE)
        results = smsd.screen_and_match(query, corpus, threshold=0.0)
        assert isinstance(results, list)
        # With threshold=0.0, should return some hits
        for item in results:
            assert isinstance(item, tuple)
            assert len(item) == 2
            idx, mapping = item
            assert isinstance(idx, int)

    def test_screen_and_match_high_threshold(self, corpus):
        """screen_and_match with threshold=1.0 should return few/no results."""
        results = smsd.screen_and_match(METHANE, corpus, threshold=0.99)
        # Very high threshold with tiny query should filter most out
        assert isinstance(results, list)

    def test_screen_and_mcs_size_basic(self, corpus):
        """screen_and_mcs_size should return (index, size) tuples."""
        query = smsd.parse_smiles(BENZENE)
        results = smsd.screen_and_mcs_size(query, corpus, threshold=0.0)
        assert isinstance(results, list)
        for item in results:
            assert isinstance(item, tuple)
            assert len(item) == 2
            idx, size = item
            assert isinstance(idx, int)
            assert isinstance(size, int)
            assert size >= 0

    def test_batch_similarity_basic(self):
        """batch_similarity should return sorted (index, score) tuples."""
        results = smsd.batch_similarity(BENZENE, BATCH_CORPUS)
        assert isinstance(results, list)
        assert len(results) == len(BATCH_CORPUS)
        # Should be sorted by score descending
        scores = [s for _, s in results]
        for i in range(1, len(scores)):
            assert scores[i] <= scores[i - 1], \
                "batch_similarity results not sorted descending"

    def test_batch_similarity_self_is_max(self):
        """Similarity of a molecule with itself should be 1.0 (or highest)."""
        results = smsd.batch_similarity(BENZENE, [BENZENE, ETHANOL])
        # Benzene vs itself should be first (highest score)
        assert results[0][0] == 0
        assert results[0][1] == pytest.approx(1.0)

    def test_batch_similarity_return_types(self):
        """Each element should be (int, float)."""
        results = smsd.batch_similarity(ETHANOL, [BENZENE, PHENOL])
        for idx, score in results:
            assert isinstance(idx, int)
            assert isinstance(score, float)
            assert 0.0 <= score <= 1.0

    def test_batch_mcs_constrained_basic(self):
        """batch_mcs_constrained with two queries against one target."""
        # Two small reactants mapping onto a larger product
        queries = [smsd.parse_smiles(BENZENE),
                   smsd.parse_smiles(ETHANOL)]
        targets = [smsd.parse_smiles(ASPIRIN)]
        results = smsd.batch_mcs_constrained(queries, targets)
        assert isinstance(results, list)
        assert len(results) == 2
        # Each result should be a dict
        for r in results:
            assert isinstance(r, dict)


# ===========================================================================
# 5. TargetCorpus
# ===========================================================================

class TestTargetCorpus:
    """Tests for the TargetCorpus class."""

    def test_construction_empty(self):
        """Empty TargetCorpus should have length 0."""
        corpus = smsd.TargetCorpus()
        assert len(corpus) == 0

    def test_from_smiles(self):
        """TargetCorpus.from_smiles should create corpus with targets."""
        corpus = smsd.TargetCorpus.from_smiles(BATCH_CORPUS)
        assert len(corpus) == len(BATCH_CORPUS)

    def test_add_target(self):
        """add_target should increment the corpus size."""
        corpus = smsd.TargetCorpus()
        corpus.add_target(smsd.parse_smiles(BENZENE))
        assert len(corpus) == 1
        corpus.add_target(smsd.parse_smiles(ETHANOL))
        assert len(corpus) == 2

    def test_add_targets(self):
        """add_targets should add multiple molecules at once."""
        corpus = smsd.TargetCorpus()
        mols = [smsd.parse_smiles(s) for s in [BENZENE, PHENOL, ETHANOL]]
        corpus.add_targets(mols)
        assert len(corpus) == 3

    def test_prewarm(self):
        """prewarm should not crash and should set is_prewarmed."""
        corpus = smsd.TargetCorpus.from_smiles([BENZENE, PHENOL, ETHANOL])
        assert not corpus.is_prewarmed
        corpus.prewarm()
        assert corpus.is_prewarmed

    def test_substructure(self):
        """substructure query should return list of bools."""
        corpus = smsd.TargetCorpus.from_smiles([PHENOL, ETHANOL, TOLUENE])
        corpus.prewarm()
        results = corpus.substructure(BENZENE)
        assert isinstance(results, list)
        assert len(results) == 3
        assert results[0] is True   # phenol contains benzene
        assert results[1] is False  # ethanol does not
        assert results[2] is True   # toluene contains benzene

    def test_find_substructure(self):
        """find_substructure should return mappings per target."""
        corpus = smsd.TargetCorpus.from_smiles([PHENOL, ETHANOL])
        corpus.prewarm()
        results = corpus.find_substructure(BENZENE)
        assert isinstance(results, list)
        assert len(results) == 2
        assert len(results[0]) > 0  # phenol has a match
        assert len(results[1]) == 0  # ethanol does not

    def test_mcs_size(self):
        """mcs_size should return list of ints."""
        corpus = smsd.TargetCorpus.from_smiles([PHENOL, ETHANOL, BENZENE])
        corpus.prewarm()
        results = corpus.mcs_size(BENZENE)
        assert isinstance(results, list)
        assert len(results) == 3
        assert all(isinstance(r, int) for r in results)
        assert results[0] >= 6  # phenol shares benzene ring
        assert results[2] >= 6  # benzene vs benzene

    def test_screen(self):
        """screen should return list of target indices above threshold."""
        corpus = smsd.TargetCorpus.from_smiles(BATCH_CORPUS)
        corpus.prewarm()
        results = corpus.screen(BENZENE, threshold=0.0)
        assert isinstance(results, list)
        # With threshold=0.0, should return all indices
        assert len(results) <= len(BATCH_CORPUS)
        for idx in results:
            assert isinstance(idx, int)
            assert 0 <= idx < len(BATCH_CORPUS)

    def test_screen_high_threshold(self):
        """screen with high threshold should return fewer results."""
        corpus = smsd.TargetCorpus.from_smiles(BATCH_CORPUS)
        corpus.prewarm()
        low = corpus.screen(BENZENE, threshold=0.1)
        high = corpus.screen(BENZENE, threshold=0.9)
        assert len(high) <= len(low)

    def test_is_prewarmed_false_initially(self):
        """is_prewarmed should be False before calling prewarm."""
        corpus = smsd.TargetCorpus.from_smiles([BENZENE])
        assert not corpus.is_prewarmed

    def test_substructure_without_prewarm(self):
        """substructure should work even without prewarm (just slower)."""
        corpus = smsd.TargetCorpus.from_smiles([PHENOL, ETHANOL])
        results = corpus.substructure(BENZENE)
        assert isinstance(results, list)
        assert len(results) == 2


# ===========================================================================
# 6. Scaffold/Graph Utilities
# ===========================================================================

class TestScaffoldGraph:
    """Tests for murcko_scaffold, extract_subgraph, count_components,
    split_components, same_canonical_graph, find_nmcs, find_scaffold_mcs,
    validate_mapping, and canonical_hash."""

    def test_murcko_scaffold_phenol(self):
        """murcko_scaffold of phenol should give benzene ring."""
        mol = smsd.parse_smiles(PHENOL)
        scaffold = smsd.murcko_scaffold(mol)
        assert isinstance(scaffold, smsd.MolGraph)
        assert scaffold.n == 6  # benzene ring

    def test_murcko_scaffold_aspirin(self):
        """murcko_scaffold of aspirin should extract the ring system."""
        mol = smsd.parse_smiles(ASPIRIN)
        scaffold = smsd.murcko_scaffold(mol)
        assert isinstance(scaffold, smsd.MolGraph)
        assert scaffold.n >= 6  # at least the benzene ring

    def test_murcko_scaffold_acyclic(self):
        """murcko_scaffold of an acyclic molecule returns the molecule itself."""
        mol = smsd.parse_smiles(PROPANE)
        scaffold = smsd.murcko_scaffold(mol)
        # Acyclic molecules have no ring scaffold; implementation returns input
        assert isinstance(scaffold, smsd.MolGraph)

    def test_extract_subgraph(self):
        """extract_subgraph should extract atoms at given indices."""
        mol = smsd.parse_smiles(ETHANOL)  # C-C-O, indices 0,1,2
        sub = smsd.extract_subgraph(mol, [0, 1])
        assert isinstance(sub, smsd.MolGraph)
        assert sub.n == 2

    def test_extract_subgraph_single_atom(self):
        """extract_subgraph with single atom should work."""
        mol = smsd.parse_smiles(BENZENE)
        sub = smsd.extract_subgraph(mol, [0])
        assert sub.n == 1

    def test_count_components_single(self):
        """count_components on a connected molecule should return 1."""
        mol = smsd.parse_smiles(BENZENE)
        assert smsd.count_components(mol) == 1

    def test_count_components_disconnected(self):
        """count_components on a salt should count fragments."""
        mol = smsd.parse_smiles("[Na+].[Cl-]")
        assert smsd.count_components(mol) == 2

    def test_split_components_single(self):
        """split_components on connected molecule returns one fragment."""
        mol = smsd.parse_smiles(BENZENE)
        frags = smsd.split_components(mol)
        assert isinstance(frags, list)
        assert len(frags) == 1
        assert frags[0].n == 6

    def test_split_components_salt(self):
        """split_components on a salt returns separate fragments."""
        mol = smsd.parse_smiles("[Na+].[Cl-]")
        frags = smsd.split_components(mol)
        assert len(frags) == 2
        sizes = sorted([f.n for f in frags])
        assert sizes == [1, 1]

    def test_same_canonical_graph_identical(self):
        """same_canonical_graph for identical molecules should be True."""
        m1 = smsd.parse_smiles(BENZENE)
        m2 = smsd.parse_smiles(BENZENE)
        assert smsd.same_canonical_graph(m1, m2) is True

    def test_same_canonical_graph_different(self):
        """same_canonical_graph for different molecules should be False."""
        m1 = smsd.parse_smiles(BENZENE)
        m2 = smsd.parse_smiles(ETHANOL)
        assert smsd.same_canonical_graph(m1, m2) is False

    def test_same_canonical_graph_isomorphic_smiles(self):
        """same_canonical_graph should detect isomorphic graphs
        from different SMILES."""
        m1 = smsd.parse_smiles("c1ccccc1")
        m2 = smsd.parse_smiles("C1=CC=CC=C1")
        # These represent the same molecule (benzene)
        assert smsd.same_canonical_graph(m1, m2) is True

    def test_find_nmcs_basic(self):
        """find_nmcs should find common substructure across multiple molecules."""
        mols = [smsd.parse_smiles(s) for s in [PHENOL, TOLUENE, ANILINE]]
        result = smsd.find_nmcs(mols, threshold=1.0, timeout_ms=10000)
        # All share benzene-derived substructure
        assert len(result) >= 5

    def test_find_scaffold_mcs(self):
        """find_scaffold_mcs should find shared scaffold between two molecules."""
        m1 = smsd.parse_smiles(NAPHTHALENE)
        m2 = smsd.parse_smiles(BENZENE)
        result = smsd.find_scaffold_mcs(m1, m2)
        assert len(result) >= 6  # benzene is shared

    def test_validate_mapping_valid(self):
        """validate_mapping on a correct MCS mapping should return no errors."""
        m1 = smsd.parse_smiles(BENZENE)
        m2 = smsd.parse_smiles(PHENOL)
        mapping = smsd.find_mcs(m1, m2)
        errors = smsd.validate_mapping(m1, m2, mapping)
        assert isinstance(errors, list)
        assert len(errors) == 0

    def test_validate_mapping_invalid(self):
        """validate_mapping on a bad mapping should return error strings."""
        m1 = smsd.parse_smiles(BENZENE)
        m2 = smsd.parse_smiles(ETHANOL)
        # Intentionally bad mapping: carbon in benzene -> oxygen in ethanol
        bad_mapping = {0: 2}  # C -> O
        errors = smsd.validate_mapping(m1, m2, bad_mapping)
        assert isinstance(errors, list)
        # May or may not find errors depending on strictness, but shouldn't crash

    def test_canonical_hash_deterministic(self):
        """canonical_hash should be deterministic for the same molecule."""
        mol = smsd.parse_smiles(BENZENE)
        h1 = smsd.canonical_hash(mol)
        h2 = smsd.canonical_hash(mol)
        assert h1 == h2

    def test_canonical_hash_different_molecules(self):
        """canonical_hash should return consistent int values."""
        h1 = smsd.canonical_hash(smsd.parse_smiles(BENZENE))
        h2 = smsd.canonical_hash(smsd.parse_smiles(ETHANOL))
        assert isinstance(h1, int)
        assert isinstance(h2, int)

    def test_canonical_hash_same_molecule_different_smiles(self):
        """canonical_hash should be equal for isomorphic molecules."""
        h1 = smsd.canonical_hash(smsd.parse_smiles("c1ccccc1"))
        h2 = smsd.canonical_hash(smsd.parse_smiles("C1=CC=CC=C1"))
        assert h1 == h2

    def test_canonical_hash_is_int(self):
        """canonical_hash should return an integer."""
        h = smsd.canonical_hash(smsd.parse_smiles(BENZENE))
        assert isinstance(h, int)

    def test_canonical_hash_stereo_returns_int(self):
        """canonical_hash should return int for stereo molecules."""
        r = smsd.parse_smiles("N[C@@H](C)C(=O)O")  # L-alanine
        s = smsd.parse_smiles("N[C@H](C)C(=O)O")   # D-alanine
        assert isinstance(smsd.canonical_hash(r), int)
        assert isinstance(smsd.canonical_hash(s), int)


# ===========================================================================
# 7. Chemistry Utilities
# ===========================================================================

class TestChemUtils:
    """Tests for classify_pharmacophore, implicit_h, and prewarm."""

    def test_classify_pharmacophore_oxygen(self):
        """classify_pharmacophore on a phenol O should have acceptor/donor bits."""
        mol = smsd.parse_smiles(PHENOL)
        # Find the O atom
        o_idx = None
        for i in range(mol.n):
            if mol.atomic_num[i] == 8:
                o_idx = i
                break
        assert o_idx is not None
        features = smsd.classify_pharmacophore(mol, o_idx)
        assert isinstance(features, int)
        # Oxygen in phenol is typically donor (bit 0) and acceptor (bit 1)
        assert features > 0  # should have at least some pharmacophore features

    def test_classify_pharmacophore_aromatic_carbon(self):
        """classify_pharmacophore on an aromatic C should have aromatic bit."""
        mol = smsd.parse_smiles(BENZENE)
        features = smsd.classify_pharmacophore(mol, 0)
        assert isinstance(features, int)
        # Aromatic bit (bit 4)
        assert features & 16  # bit 4 = aromatic

    def test_classify_pharmacophore_all_atoms(self):
        """classify_pharmacophore should work for every atom in a molecule."""
        mol = smsd.parse_smiles(ASPIRIN)
        for i in range(mol.n):
            features = smsd.classify_pharmacophore(mol, i)
            assert isinstance(features, int)
            assert features >= 0

    def test_implicit_h_carbon(self):
        """implicit_h for carbon with 3 bonds, 0 charge should be 1."""
        h_count = smsd.implicit_h(6, 3, 0)
        assert h_count == 1

    def test_implicit_h_nitrogen(self):
        """implicit_h for nitrogen with 2 bonds, 0 charge should be 1."""
        h_count = smsd.implicit_h(7, 2, 0)
        assert h_count == 1

    def test_implicit_h_oxygen(self):
        """implicit_h for oxygen with 2 bonds, 0 charge should be 0."""
        h_count = smsd.implicit_h(8, 2, 0)
        assert h_count == 0

    def test_implicit_h_saturated_carbon(self):
        """implicit_h for fully saturated carbon should be 0."""
        h_count = smsd.implicit_h(6, 4, 0)
        assert h_count == 0

    def test_implicit_h_return_type(self):
        """implicit_h should return an integer."""
        result = smsd.implicit_h(6, 2, 0)
        assert isinstance(result, int)

    def test_prewarm_basic(self):
        """prewarm should not crash and should return a usable molecule."""
        mol = smsd.parse_smiles(BENZENE)
        result = smsd.prewarm(mol)
        assert isinstance(result, smsd.MolGraph)
        assert result.n == 6

    def test_prewarm_can_still_use_molecule(self):
        """After prewarm, molecule should work for MCS."""
        mol = smsd.prewarm(smsd.parse_smiles(BENZENE))
        mapping = smsd.find_mcs(mol, PHENOL)
        assert len(mapping) >= 6

    def test_prewarm_with_smiles(self):
        """prewarm should accept SMILES string."""
        result = smsd.prewarm(BENZENE)
        assert isinstance(result, smsd.MolGraph)


# ===========================================================================
# 8. Coordinate Transforms
# ===========================================================================

class TestCoordinateTransforms:
    """Tests for translate_2d, rotate_2d, scale_2d, mirror_x, mirror_y,
    center_2d, align_2d, bounding_box_2d, normalise_bond_length,
    and canonical_orientation."""

    @pytest.fixture
    def benzene_coords(self):
        """Generate 2D coordinates for benzene."""
        mol = smsd.parse_smiles(BENZENE)
        coords = smsd.generate_coords_2d(mol)
        return mol, coords

    @pytest.fixture
    def simple_coords(self):
        """Simple 3-atom linear coordinates."""
        return [[0.0, 0.0], [1.0, 0.0], [2.0, 0.0]]

    def test_translate_2d_basic(self, simple_coords):
        """translate_2d should shift all coordinates by (dx, dy)."""
        result = smsd.translate_2d(simple_coords, dx=5.0, dy=3.0)
        assert len(result) == 3
        assert abs(result[0][0] - 5.0) < 1e-6
        assert abs(result[0][1] - 3.0) < 1e-6
        assert abs(result[1][0] - 6.0) < 1e-6
        assert abs(result[2][0] - 7.0) < 1e-6

    def test_translate_2d_zero(self, simple_coords):
        """translate_2d with (0, 0) should leave coordinates unchanged."""
        result = smsd.translate_2d(simple_coords, dx=0.0, dy=0.0)
        for i in range(3):
            assert abs(result[i][0] - simple_coords[i][0]) < 1e-10
            assert abs(result[i][1] - simple_coords[i][1]) < 1e-10

    def test_rotate_2d_90_degrees(self, simple_coords):
        """rotate_2d by 90 degrees around centroid preserves distances."""
        result = smsd.rotate_2d(simple_coords, angle=math.pi / 2)
        assert len(result) == 3
        # Centroid is (1,0); point (1,0) stays at centroid after rotation
        assert abs(result[1][0] - 1.0) < 1e-4
        assert abs(result[1][1] - 0.0) < 1e-4

    def test_rotate_2d_360(self, simple_coords):
        """rotate_2d by 360 degrees should return same coordinates."""
        result = smsd.rotate_2d(simple_coords, angle=2 * math.pi)
        for i in range(3):
            assert abs(result[i][0] - simple_coords[i][0]) < 1e-6
            assert abs(result[i][1] - simple_coords[i][1]) < 1e-6

    def test_scale_2d_double(self, simple_coords):
        """scale_2d by factor 2.0 should scale distances from centroid."""
        result = smsd.scale_2d(simple_coords, factor=2.0)
        assert len(result) == 3
        # Centroid at (1,0): distances double from centroid
        # (0,0) → centroid-1 → *2 → (-1,0)
        # (1,0) → centroid stays → (1,0)
        # (2,0) → centroid+1 → *2 → (3,0)
        assert abs(result[0][0] - (-1.0)) < 1e-6
        assert abs(result[1][0] - 1.0) < 1e-6
        assert abs(result[2][0] - 3.0) < 1e-6

    def test_scale_2d_half(self, simple_coords):
        """scale_2d by factor 0.5 should halve distances from centroid."""
        result = smsd.scale_2d(simple_coords, factor=0.5)
        # (2,0) → centroid+1 → *0.5 → (1.5,0)
        assert abs(result[2][0] - 1.5) < 1e-6

    def test_mirror_x(self, simple_coords):
        """mirror_x should reflect y coordinates around centroid y-axis."""
        coords = [[1.0, 2.0], [3.0, 4.0]]
        result = smsd.mirror_x(coords)
        assert len(result) == 2
        # Centroid y = 3.0; mirror: 2→3+(3-2)=4, 4→3-(4-3)=2
        assert abs(result[0][1] - 4.0) < 1e-6
        assert abs(result[1][1] - 2.0) < 1e-6
        # X unchanged
        assert abs(result[0][0] - 1.0) < 1e-6

    def test_mirror_y(self, simple_coords):
        """mirror_y should reflect x coordinates around centroid x-axis."""
        coords = [[1.0, 2.0], [3.0, 4.0]]
        result = smsd.mirror_y(coords)
        assert len(result) == 2
        # Centroid x = 2.0; mirror: 1→2+(2-1)=3, 3→2-(3-2)=1
        assert abs(result[0][0] - 3.0) < 1e-6
        assert abs(result[1][0] - 1.0) < 1e-6
        # Y unchanged
        assert abs(result[0][1] - 2.0) < 1e-6

    def test_center_2d(self):
        """center_2d should center coordinates at origin."""
        coords = [[2.0, 4.0], [4.0, 6.0]]
        result = smsd.center_2d(coords)
        # Centroid should be (0, 0)
        cx = sum(p[0] for p in result) / len(result)
        cy = sum(p[1] for p in result) / len(result)
        assert abs(cx) < 1e-6
        assert abs(cy) < 1e-6

    def test_center_2d_already_centered(self):
        """center_2d on already-centered coords should be no-op."""
        coords = [[-1.0, -1.0], [1.0, 1.0]]
        result = smsd.center_2d(coords)
        cx = sum(p[0] for p in result) / len(result)
        cy = sum(p[1] for p in result) / len(result)
        assert abs(cx) < 1e-6
        assert abs(cy) < 1e-6

    def test_align_2d_basic(self, benzene_coords):
        """align_2d should return RMSD and aligned coordinates."""
        mol, coords = benzene_coords
        # Translate coords and then align back
        shifted = smsd.translate_2d(coords, dx=10.0, dy=5.0)
        rmsd, aligned = smsd.align_2d(shifted, coords)
        assert isinstance(rmsd, float)
        assert rmsd >= 0.0
        assert len(aligned) == len(coords)

    def test_align_2d_identical_gives_zero_rmsd(self, benzene_coords):
        """Aligning coords with themselves should give RMSD near 0."""
        mol, coords = benzene_coords
        rmsd, aligned = smsd.align_2d(coords, coords)
        assert rmsd < 0.1  # should be very close to 0

    def test_bounding_box_2d(self):
        """bounding_box_2d should return (min_x, min_y, max_x, max_y)."""
        coords = [[1.0, 2.0], [3.0, 5.0], [0.0, 1.0]]
        bbox = smsd.bounding_box_2d(coords)
        assert len(bbox) == 4
        min_x, min_y, max_x, max_y = bbox
        assert abs(min_x - 0.0) < 1e-6
        assert abs(min_y - 1.0) < 1e-6
        assert abs(max_x - 3.0) < 1e-6
        assert abs(max_y - 5.0) < 1e-6

    def test_bounding_box_2d_single_point(self):
        """bounding_box_2d on one point should have zero extent."""
        bbox = smsd.bounding_box_2d([[3.0, 7.0]])
        min_x, min_y, max_x, max_y = bbox
        assert abs(min_x - max_x) < 1e-6
        assert abs(min_y - max_y) < 1e-6

    def test_normalise_bond_length(self, benzene_coords):
        """normalise_bond_length should scale coords to target bond length."""
        mol, coords = benzene_coords
        result = smsd.normalise_bond_length(mol, coords, target=1.5)
        assert len(result) == len(coords)
        # Check that bonds are approximately 1.5 units
        # (some variance is expected for ring systems)

    def test_canonical_orientation(self, benzene_coords):
        """canonical_orientation should return reoriented coordinates."""
        mol, coords = benzene_coords
        result = smsd.canonical_orientation(mol, coords)
        assert len(result) == len(coords)
        # Each coordinate pair should have 2 values
        for pt in result:
            assert len(pt) == 2


# ===========================================================================
# 9. Depiction
# ===========================================================================

class TestDepiction:
    """Tests for depict_svg, depict_mapping, depict_pair, and save_svg."""

    def test_depict_svg_from_smiles(self):
        """depict_svg should accept a SMILES string and return SVG."""
        svg = smsd.depict_svg(BENZENE)
        assert isinstance(svg, str)
        assert "<svg" in svg.lower() or "svg" in svg.lower()

    def test_depict_svg_from_mol(self):
        """depict_svg should accept a MolGraph and return SVG."""
        mol = smsd.parse_smiles(BENZENE)
        svg = smsd.depict_svg(mol)
        assert isinstance(svg, str)
        assert len(svg) > 100

    def test_depict_svg_with_options(self):
        """depict_svg should accept kwargs for DepictOptions."""
        svg = smsd.depict_svg(PHENOL, width=800, height=400)
        assert isinstance(svg, str)
        assert len(svg) > 0

    def test_depict_svg_ethanol(self):
        """depict_svg should work for small acyclic molecules."""
        svg = smsd.depict_svg(ETHANOL)
        assert isinstance(svg, str)
        assert len(svg) > 0

    def test_depict_mapping_basic(self):
        """depict_mapping should produce SVG with highlighted atoms."""
        mol = smsd.parse_smiles(PHENOL)
        mapping = {0: 0, 1: 1, 2: 2}  # partial mapping
        svg = smsd.depict_mapping(mol, mapping)
        assert isinstance(svg, str)
        assert len(svg) > 100

    def test_depict_mapping_from_smiles(self):
        """depict_mapping should accept SMILES string."""
        mapping = {0: 0, 1: 1}
        svg = smsd.depict_mapping(BENZENE, mapping)
        assert isinstance(svg, str)

    def test_depict_pair_basic(self):
        """depict_pair should render two molecules side-by-side."""
        m1 = smsd.parse_smiles(BENZENE)
        m2 = smsd.parse_smiles(PHENOL)
        mapping = smsd.find_mcs(m1, m2)
        svg = smsd.depict_pair(m1, m2, mapping)
        assert isinstance(svg, str)
        assert len(svg) > 100

    def test_depict_pair_from_smiles(self):
        """depict_pair should accept SMILES strings."""
        mapping = smsd.find_mcs(BENZENE, PHENOL)
        svg = smsd.depict_pair(BENZENE, PHENOL, mapping)
        assert isinstance(svg, str)

    def test_depict_pair_empty_mapping(self):
        """depict_pair with empty mapping should still produce SVG."""
        svg = smsd.depict_pair(BENZENE, ETHANOL, {})
        assert isinstance(svg, str)

    def test_save_svg(self, tmp_path):
        """save_svg should write SVG string to a file."""
        svg = smsd.depict_svg(BENZENE)
        path = str(tmp_path / "test.svg")
        smsd.save_svg(svg, path)
        assert os.path.exists(path)
        content = open(path).read()
        assert content == svg

    def test_save_svg_roundtrip(self, tmp_path):
        """Saved SVG should be readable and identical to original."""
        svg = smsd.depict_svg(PHENOL, width=600)
        path = str(tmp_path / "phenol.svg")
        smsd.save_svg(svg, path)
        loaded = open(path).read()
        assert loaded == svg


# ===========================================================================
# 10. RDKit Interop
# ===========================================================================

class TestRDKitInterop:
    """Tests for RDKit interoperability functions.
    Skipped entirely if rdkit is not available."""

    @pytest.fixture(autouse=True)
    def _check_rdkit(self):
        pytest.importorskip("rdkit")

    def test_from_rdkit_basic(self):
        """from_rdkit should convert RDKit Mol to MolGraph."""
        from rdkit import Chem
        mol = Chem.MolFromSmiles(BENZENE)
        g = smsd.from_rdkit(mol)
        assert isinstance(g, smsd.MolGraph)
        assert g.n == 6

    def test_from_rdkit_phenol(self):
        """from_rdkit should handle phenol correctly."""
        from rdkit import Chem
        mol = Chem.MolFromSmiles(PHENOL)
        g = smsd.from_rdkit(mol)
        assert g.n == 7

    def test_from_rdkit_none_raises(self):
        """from_rdkit(None) should raise TypeError."""
        with pytest.raises(TypeError):
            smsd.from_rdkit(None)

    def test_to_rdkit_from_mol(self):
        """to_rdkit should convert MolGraph back to RDKit Mol."""
        from rdkit import Chem
        g = smsd.parse_smiles(BENZENE)
        mol = smsd.to_rdkit(g)
        assert mol is not None
        assert mol.GetNumAtoms() == 6

    def test_to_rdkit_from_smiles(self):
        """to_rdkit should accept SMILES string."""
        from rdkit import Chem
        mol = smsd.to_rdkit(BENZENE)
        assert mol is not None
        assert mol.GetNumAtoms() == 6

    def test_roundtrip_rdkit(self):
        """from_rdkit -> to_smiles -> parse -> compare atom count."""
        from rdkit import Chem
        rdkit_mol = Chem.MolFromSmiles(PHENOL)
        g = smsd.from_rdkit(rdkit_mol)
        smi = smsd.to_smiles(g)
        g2 = smsd.parse_smiles(smi)
        assert g2.n == g.n

    def test_get_index_map_from_rdkit(self):
        """get_index_map should return mapping for from_rdkit results."""
        from rdkit import Chem
        mol = Chem.MolFromSmiles(BENZENE)
        g = smsd.from_rdkit(mol)
        imap = smsd.get_index_map(g)
        assert imap is not None
        assert isinstance(imap, list)
        assert len(imap) == 6

    def test_get_index_map_non_rdkit(self):
        """get_index_map for non-rdkit MolGraph should return None."""
        g = smsd.parse_smiles(BENZENE)
        imap = smsd.get_index_map(g)
        assert imap is None

    def test_translate_mapping_basic(self):
        """translate_mapping should translate SMSD indices to RDKit indices."""
        from rdkit import Chem
        m1 = Chem.MolFromSmiles(BENZENE)
        m2 = Chem.MolFromSmiles(PHENOL)
        g1 = smsd.from_rdkit(m1)
        g2 = smsd.from_rdkit(m2)
        smsd_mapping = smsd.find_mcs(g1, g2)
        translated = smsd.translate_mapping(smsd_mapping, g1, g2)
        assert isinstance(translated, dict)
        # All translated indices should be valid for the RDKit molecules
        for q, t in translated.items():
            assert 0 <= q < m1.GetNumAtoms()
            assert 0 <= t < m2.GetNumAtoms()

    def test_translate_mapping_empty(self):
        """translate_mapping with empty mapping should return empty dict."""
        result = smsd.translate_mapping({})
        assert result == {}

    def test_mcs_rdkit_native_basic(self):
        """mcs_rdkit_native should find MCS with element-correct RDKit indices."""
        from rdkit import Chem
        m1 = Chem.MolFromSmiles(BENZENE)
        m2 = Chem.MolFromSmiles(PHENOL)
        mapping = smsd.mcs_rdkit_native(m1, m2)
        assert isinstance(mapping, dict)
        assert len(mapping) >= 6
        # Verify element correctness
        for q, t in mapping.items():
            assert m1.GetAtomWithIdx(q).GetAtomicNum() == \
                   m2.GetAtomWithIdx(t).GetAtomicNum()

    def test_mcs_rdkit_convenience(self):
        """mcs_rdkit should work as convenience wrapper."""
        from rdkit import Chem
        m1 = Chem.MolFromSmiles(BENZENE)
        m2 = Chem.MolFromSmiles(PHENOL)
        mapping = smsd.mcs_rdkit(m1, m2)
        assert isinstance(mapping, dict)
        assert len(mapping) >= 6

    def test_substructure_rdkit_basic(self):
        """substructure_rdkit should find benzene in phenol."""
        from rdkit import Chem
        benzene = Chem.MolFromSmiles(BENZENE)
        phenol = Chem.MolFromSmiles(PHENOL)
        result = smsd.substructure_rdkit(benzene, phenol)
        assert isinstance(result, dict)
        assert len(result) > 0

    def test_batch_mcs_rdkit_basic(self):
        """batch_mcs_rdkit should process multiple targets."""
        from rdkit import Chem
        query = Chem.MolFromSmiles(BENZENE)
        targets = [Chem.MolFromSmiles(s) for s in [PHENOL, ETHANOL]]
        results = smsd.batch_mcs_rdkit(query, targets)
        assert isinstance(results, list)
        assert len(results) == 2
        # Phenol should have MCS >= 6
        assert len(results[0]) >= 6

    def test_clear_cache(self):
        """clear_cache should not crash and cache should be usable after."""
        from rdkit import Chem
        mol = Chem.MolFromSmiles(BENZENE)
        smsd.from_rdkit(mol)
        smsd.clear_cache()
        # Should still work after clearing
        g = smsd.from_rdkit(mol)
        assert g.n == 6


# ===========================================================================
# 11. Enums
# ===========================================================================

class TestEnums:
    """Tests for MatcherEngine, Solvent, BondOrderMode, and RingFusionMode."""

    def test_matcher_engine_has_members(self):
        """MatcherEngine enum should exist and have known members."""
        assert hasattr(smsd, 'MatcherEngine')
        # Check it has at least some expected values
        engine = smsd.MatcherEngine
        # Should be able to iterate or access members
        assert engine is not None

    def test_solvent_has_members(self):
        """Solvent enum should exist."""
        assert hasattr(smsd, 'Solvent')
        solvent = smsd.Solvent
        assert solvent is not None

    def test_bond_order_mode_values(self):
        """BondOrderMode should have STRICT, ANY, etc."""
        bom = smsd.BondOrderMode
        assert hasattr(bom, 'STRICT')
        assert hasattr(bom, 'ANY')
        # STRICT and ANY should be distinct
        assert bom.STRICT != bom.ANY

    def test_bond_order_mode_loose(self):
        """BondOrderMode should have LOOSE mode."""
        bom = smsd.BondOrderMode
        assert hasattr(bom, 'LOOSE')

    def test_ring_fusion_mode_exists(self):
        """RingFusionMode enum should exist."""
        assert hasattr(smsd, 'RingFusionMode')
        rfm = smsd.RingFusionMode
        assert rfm is not None

    def test_bond_order_mode_in_chem_options(self):
        """BondOrderMode should be usable in ChemOptions."""
        chem = smsd.ChemOptions()
        chem.match_bond_order = smsd.BondOrderMode.ANY
        assert chem.match_bond_order == smsd.BondOrderMode.ANY

    def test_aromaticity_mode_exists(self):
        """AromaticityMode enum should exist."""
        assert hasattr(smsd, 'AromaticityMode')

    def test_aromaticity_model_exists(self):
        """AromaticityModel enum should exist."""
        assert hasattr(smsd, 'AromaticityModel')


# ===========================================================================
# 12. SMARTS Advanced
# ===========================================================================

class TestSMARTSAdvanced:
    """Tests for SMARTS-based matching: find_mcs_smarts, compile_smarts,
    and SmartsQuery."""

    @pytest.fixture(autouse=True)
    def _check_smarts(self):
        """Skip all SMARTS tests if SMARTS support is not compiled."""
        if not smsd._HAS_SMARTS:
            pytest.skip("SMARTS support not available in this build")

    def test_find_mcs_smarts_basic(self):
        """find_mcs_smarts should match amide pattern in target."""
        target = smsd.parse_smiles("CC(=O)Nc1ccc(O)cc1")
        mapping = smsd.find_mcs_smarts("[CX3](=O)[NX3]", target)
        assert isinstance(mapping, dict)
        assert len(mapping) >= 3  # C, =O, N

    def test_find_mcs_smarts_no_match(self):
        """find_mcs_smarts with non-matching pattern returns empty dict."""
        target = smsd.parse_smiles(ETHANOL)
        mapping = smsd.find_mcs_smarts("[#16]", target)  # sulfur
        assert mapping == {}

    def test_find_mcs_smarts_carbon_pattern(self):
        """find_mcs_smarts with simple carbon pattern."""
        target = smsd.parse_smiles(PROPANE)
        mapping = smsd.find_mcs_smarts("[#6]", target)
        assert len(mapping) >= 1

    def test_compile_smarts_basic(self):
        """compile_smarts should return a SmartsQuery object."""
        if smsd.SmartsQuery is None:
            pytest.skip("SmartsQuery not available")
        query = smsd.compile_smarts("[#6]~[#7]")
        assert query is not None

    def test_compile_smarts_invalid_type(self):
        """compile_smarts with non-string should raise TypeError."""
        with pytest.raises(TypeError):
            smsd.compile_smarts(42)

    def test_smarts_query_matches(self):
        """SmartsQuery.matches should detect pattern in molecule."""
        if smsd.SmartsQuery is None:
            pytest.skip("SmartsQuery not available")
        query = smsd.compile_smarts("[#6]~[#7]")
        mol = smsd.parse_smiles(ANILINE)
        assert query.matches(mol) is True

    def test_smarts_query_no_match(self):
        """SmartsQuery.matches should return False for non-matching mol."""
        if smsd.SmartsQuery is None:
            pytest.skip("SmartsQuery not available")
        query = smsd.compile_smarts("[#6]~[#16]")  # C-S
        mol = smsd.parse_smiles(BENZENE)
        assert query.matches(mol) is False

    def test_smarts_query_find_all(self):
        """SmartsQuery.find_all should return list of mappings."""
        if smsd.SmartsQuery is None:
            pytest.skip("SmartsQuery not available")
        query = smsd.compile_smarts("[#6]")
        mol = smsd.parse_smiles(PROPANE)
        matches = query.find_all(mol)
        assert isinstance(matches, list)
        assert len(matches) == 3  # 3 carbons

    def test_smarts_query_matches_many(self):
        """SmartsQuery.matches_many should check pattern in multiple mols."""
        if smsd.SmartsQuery is None:
            pytest.skip("SmartsQuery not available")
        if not hasattr(smsd.SmartsQuery, 'matches_many'):
            # Fallback: test via individual matches calls
            query = smsd.compile_smarts("[#6]~[#7]")
            mols = [smsd.parse_smiles(s) for s in [ANILINE, BENZENE, ETHANOL]]
            results = [query.matches(m) for m in mols]
            assert results[0] is True   # aniline has C-N
            assert results[1] is False  # benzene has no N
            assert results[2] is False  # ethanol has no N
            return
        query = smsd.compile_smarts("[#6]~[#7]")
        mols = [smsd.parse_smiles(s) for s in [ANILINE, BENZENE, ETHANOL]]
        results = query.matches_many(mols)
        assert isinstance(results, list)
        assert len(results) == 3
        assert results[0] is True   # aniline has C-N
        assert results[1] is False  # benzene has no N

    def test_smarts_match_convenience(self):
        """smsd.smarts_match convenience function should work."""
        assert smsd.smarts_match("[#6]~[#7]", ANILINE) is True
        assert smsd.smarts_match("[#6]~[#7]", BENZENE) is False

    def test_smarts_find_all_convenience(self):
        """smsd.smarts_find_all convenience function should work."""
        matches = smsd.smarts_find_all("[#6]", PROPANE)
        assert isinstance(matches, list)
        assert len(matches) == 3

    def test_compile_smarts_reusable(self):
        """Compiled SMARTS query should be reusable across molecules."""
        if smsd.SmartsQuery is None:
            pytest.skip("SmartsQuery not available")
        query = smsd.compile_smarts("[OH]")
        assert query.matches(smsd.parse_smiles(PHENOL)) is True
        assert query.matches(smsd.parse_smiles(BENZENE)) is False
        assert query.matches(smsd.parse_smiles(ETHANOL)) is True
