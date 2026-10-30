# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
# Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
# See the NOTICE file for attribution, trademark, and algorithm IP terms.
"""
Coverage of the public SMSD Pro Python API, one class per area.

Run with:  pytest python/tests/test_api_coverage.py -v
"""
import math
import os
import pytest
import smsd
from .search_helpers import atom_map


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
    """Prepared searches, explicit result data and selected-graph output."""

    def test_match_to_smiles_basic(self):
        query = smsd.parse_smiles(BENZENE)
        match = smsd.MCS.find(query).match(PHENOL).best
        assert match.atom_count == 6
        assert smsd.parse_smiles(match.to_smiles(query)).n == 6

    def test_empty_match_to_smiles(self):
        assert smsd.Match({}, (), ()).to_smiles(BENZENE) == ""

    def test_match_to_smiles_with_options(self):
        options = smsd.MCS.Options()
        result = smsd.MCS.find(BENZENE).with_options(options).with_timeout(10000).match(PHENOL)
        assert isinstance(result.best.to_smiles(BENZENE), str)
        assert result.best.to_smiles(BENZENE)

    def test_match_to_smiles_from_text_query(self):
        result = smsd.MCS.find(BENZENE).match(PHENOL)
        assert smsd.parse_smiles(result.best.to_smiles(BENZENE)).n == 6

    def test_match_to_smiles_identical_molecules(self):
        result = smsd.MCS.find(BENZENE).match(BENZENE)
        assert smsd.same_canonical_graph(smsd.parse_smiles(BENZENE), smsd.parse_smiles(result.best.to_smiles(BENZENE)))

    def test_match_to_smiles_with_timeout(self):
        result = smsd.MCS.find(BENZENE).with_timeout(5000).match(PHENOL)
        assert isinstance(result.best.to_smiles(BENZENE), str)
        assert result.optimal

    def test_search_reports_proof_and_statistics(self):
        result = smsd.MCS.find(BENZENE).match(PHENOL)
        assert result.best.atom_count >= 6
        assert result.optimal and result.complete
        assert result.nodes_visited >= 0 and result.elapsed_ms >= 0

    def test_search_without_observer_returns_witness(self):
        result = smsd.MCS.find(BENZENE).match(PHENOL)
        assert result.found and result.best.atom_mapping
        assert smsd.MCS.validate_mapping(BENZENE, PHENOL, result.best) == []

    def test_completed_enumeration_retains_optimum_score(self):
        query = smsd.MCS.find(BENZENE)
        first, all_results = query.match(ASPIRIN), query.match_all(ASPIRIN)
        assert first.score == all_results.score == 6
        assert all(match.atom_count == 6 for match in all_results.matches)
        assert first.optimal and all_results.complete

    def test_search_result_basic(self):
        result = smsd.MCS.find(BENZENE).match(PHENOL)
        assert isinstance(result, smsd.SearchResult)
        assert isinstance(result.best, smsd.Match)
        assert result.best.atom_count >= 6
        assert len(result.best.atom_mapping) == result.best.atom_count
        assert result.best.bond_count == 6
        assert len(result.best.fragments) == 1

    def test_search_result_smiles(self):
        result = smsd.MCS.find(BENZENE).match(PHENOL)
        assert isinstance(result.best.to_smiles(BENZENE), str)

    def test_search_result_counts(self):
        match = smsd.MCS.find(BENZENE).match(PHENOL).best
        assert len(match.atom_mapping) == match.atom_count
        assert len(match.bond_mapping) == match.bond_count

    def test_search_result_repr(self):
        rendered = repr(smsd.MCS.find(BENZENE).match(PHENOL))
        assert "SearchResult" in rendered and "status=" in rendered
        assert "optimal=" in rendered and "complete=" in rendered

    def test_smiles_search_materialises_common_graph(self):
        result = smsd.MCS.find(BENZENE).match(PHENOL)
        assert isinstance(result, smsd.SearchResult)
        assert result.best.atom_count >= 6
        assert result.best.to_smiles(BENZENE)

    def test_smiles_search_identical(self):
        result = smsd.MCS.find(ETHANOL).match(ETHANOL)
        assert result.best.atom_count == 3
        assert result.score == 3 and result.optimal

    def test_empty_smiles_is_rejected(self):
        with pytest.raises(ValueError):
            smsd.MCS.find("").match(BENZENE)
        with pytest.raises(ValueError):
            smsd.MCS.find(BENZENE).match("")

    def test_match_construction(self):
        match = smsd.Match({0: 0, 1: 1}, ((0, 1, 0, 1),), ((0, 1),))
        assert match.atom_count == 2 and match.bond_count == 1
        assert smsd.parse_smiles(match.to_smiles("CCN")).n == 2
        assert match.fragments == ((0, 1),)

    def test_match_zero_atoms(self):
        match = smsd.Match({}, (), ())
        assert match.atom_count == 0 and match.bond_count == 0
        assert match.to_smiles(smsd.parse_smiles("")) == ""


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

    def test_read_mol_file(self, tmp_path):
        """read_mol_file should read a .mol file from disk."""
        mol = smsd.parse_smiles(BENZENE)
        path = str(tmp_path / "test.mol")
        smsd.write_molfile(mol, path)
        mol2 = smsd.read_mol_file(path)
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

    def test_write_sdf_basic(self, tmp_path):
        """write_sdf should write multiple molecules to an SDF file."""
        mols = [smsd.parse_smiles(s) for s in [BENZENE, ETHANOL, PHENOL]]
        path = str(tmp_path / "multi.sdf")
        smsd.write_sdf(mols, path)
        assert os.path.exists(path)
        content = open(path).read()
        assert content.count("$$$$") == 3

    def test_write_sdf_single_molecule(self, tmp_path):
        """write_sdf with a single molecule should produce valid SDF."""
        mols = [smsd.parse_smiles(METHANE)]
        path = str(tmp_path / "single.sdf")
        smsd.write_sdf(mols, path)
        content = open(path).read()
        assert content.count("$$$$") == 1

    def test_roundtrip_smiles_mol_smiles(self, tmp_path):
        """SMILES -> write MOL -> read MOL -> SMILES should produce
        an equivalent molecule (same atom count)."""
        for smi in [BENZENE, ETHANOL, PHENOL, ACETIC_ACID]:
            mol1 = smsd.parse_smiles(smi)
            path = str(tmp_path / "roundtrip.mol")
            smsd.write_molfile(mol1, path)
            mol2 = smsd.read_mol_file(path)
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

    def test_match_many_substructure_basic(self, corpus):
        results = smsd.Substructure.find(BENZENE).match_many(corpus)
        assert isinstance(results, list) and len(results) == len(corpus)
        assert all(isinstance(result, smsd.SearchResult) for result in results)
        assert all(result.complete for result in results)

    def test_match_many_benzene_in_phenol(self, corpus):
        results = smsd.Substructure.find(BENZENE).match_many(BATCH_CORPUS)
        assert [result.found for result in results[:4]] == [True, True, True, False]

    def test_match_many_substructure_witnesses(self, corpus):
        results = smsd.Substructure.find(BENZENE).match_many(corpus)
        assert len(results) == len(corpus)
        assert all(result.best.atom_count == 6 for result in results if result.found)

    def test_match_many_substructure_presence(self):
        results = smsd.Substructure.find(BENZENE).match_many([PHENOL, ETHANOL])
        assert results[0].best.atom_mapping
        assert not results[1].found and results[1].complete

    def test_match_many_mcs_basic(self, corpus):
        results = smsd.MCS.find(BENZENE).match_many(corpus)
        assert isinstance(results, list) and len(results) == len(corpus)
        assert all(result.optimal for result in results)

    def test_match_many_mcs_benzene_phenol_size(self):
        results = smsd.MCS.find(BENZENE).match_many([PHENOL])
        assert len(results) == 1 and results[0].best.atom_count >= 6

    def test_match_many_mcs_scores(self, corpus):
        results = smsd.MCS.find(BENZENE).match_many(corpus)
        assert len(results) == len(corpus)
        assert all(result.score >= 0 for result in results)
        assert all(result.score == (result.best.atom_count if result.found else 0) for result in results)

    def test_match_many_mcs_size_values(self):
        results = smsd.MCS.find(BENZENE).match_many([PHENOL, ETHANOL])
        assert results[0].best.atom_count >= 6
        assert results[1].best.atom_count >= 1

    def test_screened_search_returns_original_target_indices(self, corpus):
        query = smsd.parse_smiles(BENZENE)
        indices = smsd.screen_targets(query, corpus, 0.0)
        results = smsd.Substructure.find(query).match_many([corpus[index] for index in indices])
        hits = [(index, result.best.atom_mapping) for index, result in zip(indices, results) if result.found]
        assert hits
        assert all(0 <= index < len(corpus) and len(mapping) == 6 for index, mapping in hits)

    def test_screened_search_high_threshold(self, corpus):
        query = smsd.parse_smiles(BENZENE)
        indices = smsd.screen_targets(query, corpus, 0.99)
        results = smsd.Substructure.find(query).match_many([corpus[index] for index in indices])
        assert isinstance(results, list)
        assert all(result.complete for result in results)

    def test_screened_mcs_reports_scores(self, corpus):
        query = smsd.parse_smiles(BENZENE)
        indices = smsd.screen_targets(query, corpus, 0.0)
        results = smsd.MCS.find(query).match_many([corpus[index] for index in indices])
        scores = [(index, result.best.atom_count if result.found else 0) for index, result in zip(indices, results)]
        assert scores
        assert all(isinstance(index, int) and isinstance(size, int) and size >= 0 for index, size in scores)

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

    def test_constrained_queries_do_not_reuse_target_atoms(self):
        queries = [smsd.parse_smiles(BENZENE), smsd.parse_smiles(ETHANOL)]
        results = smsd.MCS.match_constrained(queries, smsd.parse_smiles(ASPIRIN))
        assert len(results) == 2
        mappings = [result.best.atom_mapping if result.found else {} for result in results]
        assert mappings[0] and mappings[1]
        assert not set(mappings[0].values()) & set(mappings[1].values())


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
        corpus = smsd.TargetCorpus.from_smiles([PHENOL, ETHANOL, TOLUENE])
        corpus.prewarm()
        results = smsd.Substructure.find(BENZENE).match_many(corpus)
        assert [result.found for result in results] == [True, False, True]

    def test_substructure_mappings(self):
        corpus = smsd.TargetCorpus.from_smiles([PHENOL, ETHANOL])
        corpus.prewarm()
        results = smsd.Substructure.find(BENZENE).match_many(corpus)
        assert len(results) == 2
        assert results[0].best.atom_mapping and not results[1].found

    def test_corpus_mcs_sizes(self):
        corpus = smsd.TargetCorpus.from_smiles([PHENOL, ETHANOL, BENZENE])
        corpus.prewarm()
        results = smsd.MCS.find(BENZENE).match_many(corpus)
        assert len(results) == 3
        assert results[0].best.atom_count >= 6 and results[2].best.atom_count >= 6
        assert all(result.optimal for result in results)

    def test_screen(self):
        """screen should return list of target indices above threshold."""
        corpus = smsd.TargetCorpus.from_smiles(BATCH_CORPUS)
        corpus.prewarm()
        results = corpus.screen(BENZENE, threshold=0.0)
        assert isinstance(results, list)
        # With threshold=0.0, should return all indices
        assert sorted(results) == list(range(len(BATCH_CORPUS)))
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
        corpus = smsd.TargetCorpus.from_smiles([PHENOL, ETHANOL])
        results = smsd.Substructure.find(BENZENE).match_many(corpus)
        assert len(results) == 2 and results[0].found and not results[1].found


# ===========================================================================
# 6. Scaffold/Graph Utilities
# ===========================================================================

class TestScaffoldGraph:
    """Tests for murcko_scaffold, extract_subgraph, count_components,
    split_components, same_canonical_graph, MCS consensus, scaffold MCS,
    MCS.validate_mapping and canonical_hash."""

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

    @pytest.mark.parametrize("smiles,hydrogens,expected", [
        (TOLUENE, [1] * 6, "c1ccccc1"),
        ("CC1(C)CCCCC1", [2] * 6, "C1CCCCC1"),
        ("O=C1CCCCC1", [0, 0] + [2] * 5, "C1(CCCCC1)=O"),
        ("O=c1cccc[nH]1", [0, 0] + [1] * 5, "c1(cccc[nH]1)=O"),
        ("Cn1cccc1", [1] * 5, "c1cc[nH]c1"),
        ("C1CCC2(CC1)CCCC2", [2, 2, 2, 0, 2, 2, 2, 2, 2, 2], "C12(CCCC1)CCCCC2"),
    ])
    def test_murcko_scaffold_replaces_stripped_substituents_with_hydrogen(self, smiles, hydrogens, expected):
        scaffold = smsd.murcko_scaffold(smsd.parse_smiles(smiles))
        assert list(scaffold.hydrogen_count) == hydrogens
        assert smsd.to_smiles(scaffold) == expected

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

    def test_consensus_basic(self):
        molecules = [smsd.parse_smiles(smiles) for smiles in [PHENOL, TOLUENE, ANILINE]]
        options = smsd.MCS.Options()
        options.timeout_ms = 10000
        result = smsd.MCS.consensus(molecules, options=options, threshold=1.0)
        assert result.found and result.graph.n >= 5
        assert result.support_count == 3
        assert all(len(match.atom_mapping) == result.graph.n for match in result.mappings.values())

    def test_scaffold_mcs(self):
        query = smsd.murcko_scaffold(smsd.parse_smiles(NAPHTHALENE))
        target = smsd.murcko_scaffold(smsd.parse_smiles(BENZENE))
        result = smsd.MCS.find(query).match(target)
        assert result.best.atom_count >= 6

    def test_validate_mapping_valid(self):
        result = smsd.MCS.find(BENZENE).match(PHENOL)
        assert smsd.MCS.validate_mapping(BENZENE, PHENOL, result.best) == []

    def test_validate_mapping_invalid(self):
        errors = smsd.MCS.validate_mapping(BENZENE, ETHANOL, {0: 2})
        assert isinstance(errors, list) and errors
        # Benzene carbon 0 cannot map to the ethanol oxygen.

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
        mapping = atom_map(smsd.MCS.find(mol).match(PHENOL))
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
        """center_2d should centre coordinates on the origin."""
        coords = [[2.0, 4.0], [4.0, 6.0]]
        result = smsd.center_2d(coords)
        # Centroid should be (0, 0)
        cx = sum(p[0] for p in result) / len(result)
        cy = sum(p[1] for p in result) / len(result)
        assert abs(cx) < 1e-6
        assert abs(cy) < 1e-6

    def test_center_2d_already_centred(self):
        """center_2d on already centred coordinates should be no-op."""
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
        for first in range(mol.n):
            for second in range(first + 1, mol.n):
                if mol.has_bond(first, second):
                    assert math.dist(result[first], result[second]) == pytest.approx(1.5, abs=1e-6)

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
        mapping = atom_map(smsd.MCS.find(m1).match(m2))
        svg = smsd.depict_pair(m1, m2, mapping)
        assert isinstance(svg, str)
        assert len(svg) > 100

    def test_depict_pair_from_smiles(self):
        """depict_pair should accept SMILES strings."""
        mapping = atom_map(smsd.MCS.find(BENZENE).match(PHENOL))
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
        smsd_mapping = atom_map(smsd.MCS.find(g1).match(g2))
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

    def test_mcs_rdkit_graph_inputs(self):
        """An MCS of RDKit molecules pairs atoms of equal elements by RDKit index."""
        from rdkit import Chem
        m1 = Chem.MolFromSmiles(BENZENE)
        m2 = Chem.MolFromSmiles(PHENOL)
        mapping = atom_map(smsd.MCS.find(m1).match(m2))
        assert isinstance(mapping, dict)
        assert len(mapping) >= 6
        # Verify element correctness
        for q, t in mapping.items():
            assert m1.GetAtomWithIdx(q).GetAtomicNum() == \
                   m2.GetAtomWithIdx(t).GetAtomicNum()

    def test_mcs_rdkit_molecule_inputs(self):
        """An MCS search accepts RDKit molecules directly."""
        from rdkit import Chem
        m1 = Chem.MolFromSmiles(BENZENE)
        m2 = Chem.MolFromSmiles(PHENOL)
        mapping = atom_map(smsd.MCS.find(m1).match(m2))
        assert isinstance(mapping, dict)
        assert len(mapping) >= 6

    def test_substructure_rdkit_molecule_inputs(self):
        """Substructure.find should find benzene in phenol."""
        from rdkit import Chem
        benzene = Chem.MolFromSmiles(BENZENE)
        phenol = Chem.MolFromSmiles(PHENOL)
        result = atom_map(smsd.Substructure.find(benzene).match(phenol))
        assert isinstance(result, dict)
        assert len(result) > 0

    def test_mcs_rdkit_many_targets(self):
        """Prepared MCS should process multiple targets."""
        from rdkit import Chem
        query = Chem.MolFromSmiles(BENZENE)
        targets = [Chem.MolFromSmiles(s) for s in [PHENOL, ETHANOL]]
        results = [atom_map(result) for result in smsd.MCS.find(query).match_many(targets)]
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
    """Tests for search objectives, Solvent, BondOrderMode and RingFusionMode."""

    def test_public_objectives_replace_engine_selectors(self):
        assert not hasattr(smsd, "MatcherEngine")
        assert smsd.MCS.Objective.MCIS != smsd.MCS.Objective.MCES
        assert smsd.MCS.Connectivity.CONNECTED != smsd.MCS.Connectivity.FRAGMENTS

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

    def test_chem_options_fluent_setters_return_the_same_object(self):
        """Every fluent setter of the native and Java options sets its field and returns self."""
        options = smsd.ChemOptions()
        assert options.with_aromaticity_mode(smsd.AromaticityMode.STRICT) is options
        assert options.with_bond_order_mode(smsd.BondOrderMode.LOOSE) is options
        assert options.with_complete_rings_only(True) is options
        assert options.with_match_isotope(True) is options
        assert options.with_ring_fusion_mode(smsd.RingFusionMode.STRICT) is options
        assert options.with_tautomer_aware(True) is options
        assert options.with_solvent(smsd.Solvent.DMSO) is options
        assert options.with_ph(3.5) is options
        assert options.aromaticity_mode == smsd.AromaticityMode.STRICT
        assert options.match_bond_order == smsd.BondOrderMode.LOOSE
        assert options.complete_rings_only and options.match_isotope and options.tautomer_aware
        assert options.ring_fusion_mode == smsd.RingFusionMode.STRICT
        assert options.solvent == smsd.Solvent.DMSO and options.pH == 3.5
        text = repr(options)
        for field in ("match_atom_type", "match_formal_charge", "match_isotope", "use_chirality", "use_bond_stereo",
                      "ring_matches_ring_only", "complete_rings_only", "ring_fusion_mode=STRICT",
                      "match_bond_order=LOOSE", "aromaticity_mode=STRICT", "tautomer_aware=True", "pH=3.5",
                      "solvent=DMSO"):
            assert field in text

    def test_chem_options_copies_are_independent(self):
        """The copy constructor, copy.copy and copy.deepcopy each give an independent ChemOptions."""
        import copy
        source = smsd.ChemOptions().with_tautomer_aware(True).with_ph(2.0)
        for duplicate in (smsd.ChemOptions(source), copy.copy(source), copy.deepcopy(source)):
            assert duplicate is not source
            assert duplicate.tautomer_aware and duplicate.pH == 2.0
            duplicate.tautomer_aware = False
            duplicate.pH = 9.0
            assert source.tautomer_aware and source.pH == 2.0


# ===========================================================================
# 12. SMARTS Advanced
# ===========================================================================

class TestSMARTSAdvanced:
    """SMARTS predicate preparation, first match and complete enumeration."""

    @pytest.fixture(autouse=True)
    def _check_smarts(self):
        assert callable(smsd.Substructure.smarts)

    def test_smarts_amide_predicate(self):
        result = smsd.Substructure.smarts("[CX3](=O)[NX3]").match("CC(=O)Nc1ccc(O)cc1")
        assert result.best.atom_count == 3 and result.complete

    def test_smarts_completed_absence(self):
        result = smsd.Substructure.smarts("[#16]").match(ETHANOL)
        assert not result.found and result.complete

    def test_smarts_carbon_pattern(self):
        assert smsd.Substructure.smarts("[#6]").match(PROPANE).best.atom_count == 1

    def test_prepare_smarts(self):
        assert isinstance(smsd.Substructure.smarts("[#6]~[#7]"), smsd.Substructure)

    def test_smarts_invalid_type(self):
        with pytest.raises((TypeError, ValueError)):
            smsd.Substructure.smarts(42)

    def test_smarts_query_matches(self):
        assert smsd.Substructure.smarts("[#6]~[#7]").match(ANILINE).found

    def test_smarts_query_no_match(self):
        assert not smsd.Substructure.smarts("[#6]~[#16]").match(BENZENE).found

    def test_smarts_query_all(self):
        result = smsd.Substructure.smarts("[#6]").match_all(PROPANE)
        assert len(result.matches) == 3 and result.complete
        assert {next(iter(match.atom_mapping.values())) for match in result.matches} == {0, 1, 2}

    def test_smarts_query_many(self):
        query = smsd.Substructure.smarts("[#6]~[#7]")
        results = query.match_many([ANILINE, BENZENE, ETHANOL])
        assert isinstance(results, list)
        assert [result.found for result in results] == [True, False, False]

    def test_smarts_first(self):
        query = smsd.Substructure.smarts("[#6]~[#7]")
        assert query.match(ANILINE).found
        assert not query.match(BENZENE).found

    def test_smarts_all_limit(self):
        result = smsd.Substructure.smarts("[#6]").match_all(PROPANE, 2)
        assert len(result.matches) == 2
        assert result.status == smsd.SearchStatus.LIMIT_REACHED and not result.complete
        assert result.optimal

    def test_smarts_reusable(self):
        query = smsd.Substructure.smarts("[OH]")
        assert query.match(PHENOL).found
        assert not query.match(BENZENE).found
        assert query.match(ETHANOL).found


# ===========================================================================
# Match counts, result representation and canonical SMILES output
# ===========================================================================

class TestResultAndSmilesWriterContracts:
    """Match counts, result representation and canonical SMILES output."""

    def test_match_counts_support_distinct_similarity_metrics(self):
        size = smsd.MCS.find(PHENOL).match(BENZENE).best.atom_count
        assert size == 6
        assert size / min(7, 6) == pytest.approx(1.0, abs=1e-6)
        assert size / (7 + 6 - size) == pytest.approx(6.0 / 7.0, abs=1e-6)

    def test_search_result_repr_exposes_proof_and_completion(self):
        result = smsd.MCS.find(BENZENE).match_all(PHENOL, 1)
        rendered = repr(result)
        assert "SearchResult" in rendered
        assert "optimal=" in rendered and "complete=" in rendered
        assert result.optimal and not result.complete

    def test_to_smiles_accepts_string_and_graph(self):
        """`smsd.to_smiles()` accepts a SMILES string as well as a MolGraph,
        and gives the same text for both.
        """
        # String input should work without raising
        out = smsd.to_smiles("c1ccccc1")
        assert isinstance(out, str)
        assert out  # non-empty
        # MolGraph input gives text as well
        g = smsd.parse_smiles("c1ccccc1")
        out2 = smsd.to_smiles(g)
        assert isinstance(out2, str)
        # The two results should be identical
        assert out == out2

    def test_to_smiles_accepts_string(self):
        """`smsd.to_smiles` must accept a SMILES string on its own."""
        out = smsd.to_smiles("c1ccccc1")
        assert isinstance(out, str)
        assert out

    def test_pyrrole_type_n_emits_nh_in_canonical_smiles(self):
        """The canonical SMILES writer emits `[nH]` for pyrrole-type
        aromatic N, so that other kekulisers (RDKit, ChemAxon and so on)
        can reconstitute the structure.

        A bare `n` would describe a pyridine-type N: 4 C + 1 pyridine-N
        gives 5 pi electrons, not the 4n+2 = 6 that Hückel validation
        in RDKit requires.
        """
        # Pyrrole — the N must emit [nH]
        out_pyrrole = smsd.to_smiles("c1cc[nH]c1")
        assert "[nH]" in out_pyrrole

        # Indole — the pyrrole-side N must also emit [nH]
        out_indole = smsd.to_smiles("c1ccc2[nH]ccc2c1")
        assert "[nH]" in out_indole

        # N-methyl pyrrole — no implicit H, must NOT emit [nH]
        out_nmepyrrole = smsd.to_smiles("Cn1cccc1")
        assert "[nH]" not in out_nmepyrrole

        # Pyridine — 6-ring, never pyrrole-type, must NOT emit [nH]
        out_pyridine = smsd.to_smiles("c1ccncc1")
        assert "[nH]" not in out_pyridine

    def test_pyrrole_type_n_round_trips_via_rdkit(self):
        """SMSD's canonical SMILES output is RDKit-parseable for pyrrole,
        indole, carbazole and benzo-pyrrole systems: RDKit can kekulise a
        five-membered aromatic ring only when its pyrrole-type N is `[nH]`.
        """
        rdkit = pytest.importorskip("rdkit.Chem")
        from rdkit.Chem import MolFromSmiles

        for smi_in in [
            "c1cc[nH]c1",          # pyrrole
            "c1ccc2[nH]ccc2c1",    # indole
            "c1ccc2cc[nH]c2c1",    # indole
            "c1ccc2c(c1)[nH]c3ccccc23",  # carbazole
        ]:
            m = smsd.parse_smiles(smi_in)
            out = smsd.to_smiles(m)
            rd = MolFromSmiles(out)
            assert rd is not None, (
                f"SMSD canonical output {out!r} (from {smi_in!r}) "
                f"is not RDKit-parseable"
            )


# ===========================================================================
# Relevant-cycle families
# ===========================================================================

class TestRelevantCycleFamilies:
    """Relevant cycles grouped into families, with exact member counts."""

    def test_families_of_naphthalene(self):
        result = smsd.compute_relevant_cycle_families(NAPHTHALENE)
        assert isinstance(result, smsd.RelevantCycleFamilies)
        assert [family.size for family in result.families] == [6, 6]
        assert result.relevant_cycle_count == 2 and not result.truncated

    def test_family_record_fields(self):
        family = smsd.compute_relevant_cycle_families(BENZENE).families[0]
        assert isinstance(family, smsd.RelevantCycleFamily)
        assert (family.size, family.root, family.first, family.second, family.opposite) == (6, 5, 1, 3, 2)
        assert family.prototype[0] == family.root and family.member_count == 1

    def test_membership_and_members(self):
        mol = smsd.parse_smiles(NAPHTHALENE)
        families = smsd.compute_relevant_cycle_families(mol)
        for index, family in enumerate(families.families):
            members = smsd.relevant_cycle_family_members(mol, family, 5)
            assert members == [family.prototype]
            assert smsd.relevant_cycle_family_of(mol, families, members[0]) == index
        assert smsd.relevant_cycle_family_of(mol, families, [0, 1, 2]) == -1

    def test_acyclic_molecule_has_no_families(self):
        result = smsd.compute_relevant_cycle_families(ETHANOL)
        assert result.families == [] and result.relevant_cycle_count == 0 and not result.truncated
