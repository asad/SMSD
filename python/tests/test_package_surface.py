# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
# Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
# See the NOTICE file for attribution, trademark, and algorithm IP terms.
"""The package's export list, deprecations, accepted input forms and lightly used exports."""
import math

import pytest
import smsd

EXPORTS_8_0_0 = [
    'AromaticityMode', 'AromaticityModel', 'BondOrderMode', 'CIPDescriptors', 'CancellationToken', 'ChemOptions',
    'ConsensusResult', 'DepictOptions', 'EZLabel', 'FingerprintFamily', 'FingerprintKey', 'FingerprintOptions',
    'FragmentMode', 'MCS', 'Match', 'MolGraph', 'MolGraphBuilder', 'Point2D',
    'Point3D', 'RSLabel', 'RingFusionMode', 'ScaffoldMode', 'SearchResult', 'SearchStatus', 'Solvent', 'Substructure',
    'TargetCorpus', 'align_2d', 'analyze_fp_quality', 'assign_cip', 'assign_ez', 'assign_rs', 'batch_similarity',
    'bounding_box_2d', 'canonical_hash', 'canonical_orientation', 'canonical_smiles', 'center_2d',
    'circular_fingerprint', 'circular_fingerprint_counts', 'classify_pharmacophore', 'clear_cache', 'compute_sssr',
    'cosine', 'count_components', 'count_cosine', 'count_crossings', 'count_dice', 'count_overlapCoefficient',
    'count_overlap_coefficient', 'count_tanimoto_coefficient', 'counts_to_array', 'dearomatize', 'depict_mapping',
    'depict_mcs', 'depict_pair', 'depict_substructure', 'depict_svg', 'dice', 'export_sdf', 'extract_subgraph',
    'fingerprint', 'fingerprint_from_smiles', 'fingerprint_subset', 'force_directed_layout', 'from_hex',
    'from_openbabel', 'from_rdkit', 'generate_coords_2d', 'generate_coords_3d', 'get_index_map', 'gpu_device_info',
    'gpu_is_available', 'implicit_h', 'is_degenerate_layout', 'kekulize', 'layout_quality', 'layout_sssr',
    'match_template', 'mcs_fingerprint', 'mirror_x', 'mirror_y', 'murcko_scaffold', 'normalise_bond_length',
    'overlapCoefficient', 'overlap_coefficient', 'parse_mapped_smiles', 'parse_smiles', 'path_fingerprint',
    'perceive_aromaticity', 'prewarm', 'prewarm_graph', 'read_mol_block', 'read_mol_file', 'read_molfile',
    'read_sdf', 'reduce_crossings', 'resolve_overlaps', 'rotate_2d', 'same_canonical_graph', 'save_svg',
    'scale_2d', 'screen_targets', 'similarity', 'similarity_upper_bound', 'soergel', 'split_components',
    'stress_majorisation', 'strip_atom_maps', 'tanimoto_coefficient', 'to_binary_string', 'to_hex', 'to_openbabel',
    'to_rdkit', 'to_smarts', 'to_smiles', 'to_svg', 'topological_torsion', 'topological_torsion_counts',
    'translate_2d', 'translate_mapping', 'translate_to_atom_ids', 'write_mol_block', 'write_mol_block_v3000',
    'write_molfile', 'write_sdf', 'write_sdf_record',
]

# Names added since 8.0.0; none is removed or renamed.
ADDED_SINCE_8_0_0 = ['fingerprints', 'rings_truncated']

BENZENE, PHENOL, ETHANOL = "c1ccccc1", "c1ccc(O)cc1", "CCO"


def test_export_list_is_the_8_0_0_surface_and_its_additions():
    assert smsd.__all__ == sorted(EXPORTS_8_0_0 + ADDED_SINCE_8_0_0)
    assert all(hasattr(smsd, name) for name in smsd.__all__)
    assert not any(hasattr(smsd, name) for name in ("Path", "sys", "importlib"))


def test_native_docstrings_make_no_unsupported_claims():
    assert "guarantee" not in smsd.MolGraph.get_canonical_hash.__doc__
    for name in ("generate_coords_2d", "generate_coords_3d", "stress_majorisation", "ring_system_sig"):
        text = (getattr(smsd._smsd, name, None) or getattr(smsd.MolGraph, name)).__doc__
        assert "%%" not in text and "SMACOF" not in text and "FNV" not in text and "MDS" not in text


# --- deprecations -----------------------------------------------------------------------------

def test_prewarm_graph_is_deprecated_and_still_prewarms():
    mol = smsd.parse_smiles(BENZENE)
    with pytest.deprecated_call():
        assert smsd.prewarm_graph(mol) is None
    assert len(mol.morgan_rank) == mol.n


def test_read_molfile_and_export_sdf_are_deprecated_and_agree_with_the_native_writers(tmp_path):
    molecules = [smsd.parse_smiles(s) for s in (BENZENE, ETHANOL)]
    smsd.write_sdf(molecules, tmp_path / "native.sdf")
    with pytest.deprecated_call():
        smsd.export_sdf(molecules, tmp_path / "deprecated.sdf")
    assert (tmp_path / "native.sdf").read_bytes() == (tmp_path / "deprecated.sdf").read_bytes()
    smsd.write_molfile(molecules[0], tmp_path / "benzene.mol")
    with pytest.deprecated_call():
        deprecated = smsd.read_molfile(tmp_path / "benzene.mol")
    assert smsd.to_smiles(deprecated) == smsd.to_smiles(smsd.read_mol_file(tmp_path / "benzene.mol"))


def test_chem_keyword_is_deprecated_in_favour_of_chemistry():
    first, second = smsd.parse_smiles(BENZENE), smsd.parse_smiles(PHENOL)
    options = smsd.ChemOptions()
    expected = smsd.similarity_upper_bound(first, second, chemistry=options)
    with pytest.deprecated_call():
        assert smsd.similarity_upper_bound(first, second, chem=options) == expected
    assert smsd.similarity_upper_bound(first, second, options) == expected
    with pytest.deprecated_call():
        assert smsd.screen_targets(first, [second], 0.0, chem=options) == smsd.screen_targets(first, [second], 0.0)


def test_depict_mcs_highlight_colour_is_deprecated():
    pytest.importorskip("rdkit.Chem")
    with pytest.deprecated_call():
        smsd.depict_mcs(BENZENE, PHENOL, highlight_color=(1.0, 0.0, 0.0))


def test_import_emits_no_deprecation_warning():
    import importlib
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("error", DeprecationWarning)
        importlib.reload(smsd._fingerprints)


# --- accepted input forms ----------------------------------------------------------------------

def test_file_functions_accept_paths_and_smiles(tmp_path):
    smsd.write_sdf([ETHANOL, smsd.parse_smiles(BENZENE)], tmp_path / "mixed.sdf")
    assert [len(mol) for mol in smsd.read_sdf(tmp_path / "mixed.sdf")] == [3, 6]
    smsd.write_molfile(ETHANOL, tmp_path / "ethanol.mol")
    assert len(smsd.read_mol_file(tmp_path / "ethanol.mol")) == 3
    with pytest.raises(ValueError):
        smsd.read_mol_file(tmp_path / "missing.mol")


def test_save_svg_writes_utf8(tmp_path):
    svg = smsd.depict_svg(BENZENE)
    smsd.save_svg(svg, tmp_path / "benzene.svg")
    assert (tmp_path / "benzene.svg").read_text(encoding="utf-8") == svg


@pytest.mark.parametrize("name", ["to_smarts", "count_components", "split_components", "murcko_scaffold",
                                  "write_mol_block", "write_sdf_record", "compute_sssr", "assign_cip",
                                  "layout_sssr", "match_template", "generate_coords_2d", "generate_coords_3d",
                                  "rings_truncated"])
def test_graph_first_functions_accept_smiles(name):
    function = getattr(smsd, name)
    assert function.__name__ == name
    assert repr(function(PHENOL)) == repr(function(smsd.parse_smiles(PHENOL)))


def test_coordinate_functions_accept_tuples_and_arrays():
    numpy = pytest.importorskip("numpy")
    mol = smsd.parse_smiles(PHENOL)
    coords = smsd.generate_coords_2d(mol)
    as_tuple, as_array = tuple(coords), numpy.array(coords)
    for function in (lambda c: smsd.count_crossings(mol, c), lambda c: smsd.layout_quality(mol, c),
                     lambda c: smsd.normalise_bond_length(mol, c), lambda c: smsd.canonical_orientation(mol, c),
                     lambda c: smsd.reduce_crossings(mol, c, 10), lambda c: smsd.force_directed_layout(mol, c, 10),
                     lambda c: smsd.stress_majorisation(mol, c, 10), smsd.center_2d, smsd.mirror_x,
                     lambda c: smsd.translate_2d(c, 1.0, 2.0), smsd.bounding_box_2d, smsd.is_degenerate_layout,
                     smsd.resolve_overlaps):
        expected = function(coords)
        assert function(as_tuple) == expected
        assert function(as_array) == expected
    assert smsd.count_crossings(PHENOL, coords) == smsd.count_crossings(mol, coords)


def test_transforms_still_reject_short_pairs():
    with pytest.raises(IndexError):
        smsd.translate_2d([[0.0]], 1.0, 1.0)


def test_similarity_metrics_accept_pairs_given_as_lists_and_sets():
    first, second = [[1, 2], [5, 1], [9, 3]], [(1, 1), (9, 1)]
    assert smsd.tanimoto_coefficient(first, second) == smsd.tanimoto_coefficient([1, 5, 9], [1, 9])
    assert smsd.overlap_coefficient(first, second) == 1.0
    assert smsd.to_binary_string({1, 5, 99}, 8) == smsd.to_binary_string([1, 5], 8)


def test_count_dice_and_cosine_accept_sparse_counts():
    first = smsd.circular_fingerprint_counts(smsd.parse_smiles(BENZENE))
    second = smsd.circular_fingerprint_counts(smsd.parse_smiles(PHENOL))
    dense = smsd.counts_to_array(first, 2048), smsd.counts_to_array(second, 2048)
    assert smsd.count_dice(first, second) == pytest.approx(smsd.count_dice(*dense))
    assert smsd.count_cosine(dict(first), dict(second)) == pytest.approx(smsd.count_cosine(*dense))
    assert smsd.count_dice([], []) == 0.0


def test_circular_fingerprint_rejects_unknown_modes():
    mol = smsd.parse_smiles(BENZENE)
    with pytest.raises(ValueError):
        smsd.circular_fingerprint(mol, 2, 2048, "FCFP")
    with pytest.raises(ValueError):
        smsd.circular_fingerprint_counts(mol, mode="ecfp4")
    assert smsd.fingerprint_from_smiles(BENZENE, mode="torsion") == smsd.fingerprint(BENZENE, kind="torsion")


def test_fingerprint_sizes_reject_booleans():
    with pytest.raises(ValueError):
        smsd.fingerprint(BENZENE, fp_size=True)
    with pytest.raises(ValueError):
        smsd.MCS.find(BENZENE).match_all(PHENOL, limit=True)


def test_lenient_smiles_parsing_recovers_from_malformed_text():
    with pytest.raises(ValueError):
        smsd.parse_smiles("C(C")
    assert len(smsd.parse_smiles("C(C", lenient=True)) >= 1


def test_validate_tautomer_consistency_accepts_matches_and_mappings():
    query, target = "Oc1ncccc1", "O=c1[nH]cccc1"
    options = smsd.MCS.Options()
    options.chemistry = smsd.ChemOptions.tautomer_profile()
    result = smsd.MCS.find(query).with_options(options).match(target)
    assert result.found
    consistent = smsd.MCS.validate_tautomer_consistency(query, target, result.best)
    assert consistent == smsd.MCS.validate_tautomer_consistency(query, target, dict(result.best.atom_mapping))
    assert isinstance(consistent, bool)
    # Proton conservation counts the stored hydrogens, so the aromatic N-H of the pyridone counts one.
    identity = {atom: atom for atom in range(7)}
    assert smsd.MCS.validate_tautomer_consistency("Oc1ccccn1", "O=c1cccc[nH]1", identity) is True
    assert smsd.MCS.validate_tautomer_consistency("Oc1ccccn1", "O=C1C=CC=CN1", identity) is True
    assert smsd.MCS.validate_tautomer_consistency("CC(=O)C", "CC(O)=C", {atom: atom for atom in range(4)}) is False


# --- in-place changes of converted graphs ------------------------------------------------------

def test_kekulised_rdkit_conversion_is_converted_afresh():
    chem = pytest.importorskip("rdkit.Chem")
    molecule = chem.MolFromSmiles(BENZENE)
    graph = smsd.from_rdkit(molecule)
    assert smsd.from_rdkit(molecule) is graph
    smsd.kekulize(graph)
    assert not any(graph.aromatic)
    fresh = smsd.from_rdkit(molecule)
    assert fresh is not graph and any(fresh.aromatic)
    assert smsd.get_index_map(graph) == smsd.get_index_map(fresh)


def test_corpus_snapshot_is_reused_until_targets_are_added():
    corpus = smsd.TargetCorpus.from_smiles([BENZENE, ETHANOL])
    first = corpus._snapshot()
    assert corpus._snapshot() is first
    corpus.add_target(PHENOL)
    second = corpus._snapshot()
    assert second is not first and len(second) == 3
    results = smsd.Substructure.find(BENZENE).match_many(corpus)
    assert [result.found for result in results] == [True, False, True]


def test_corpus_snapshot_is_rebuilt_after_prewarm():
    corpus = smsd.TargetCorpus.from_smiles([BENZENE, ETHANOL])
    cold = corpus._snapshot()
    corpus.prewarm()
    assert corpus.is_prewarmed and corpus._snapshot() is not cold
    assert [result.found for result in smsd.Substructure.find(BENZENE).match_many(corpus)] == [True, False]


# --- lightly used exports ----------------------------------------------------------------------

def test_count_tanimoto_coefficient_of_sparse_counts():
    first = smsd.circular_fingerprint_counts(smsd.parse_smiles(BENZENE))
    assert smsd.count_tanimoto_coefficient(first, first) == 1.0
    assert 0.0 < smsd.count_tanimoto_coefficient(first, smsd.circular_fingerprint_counts(smsd.parse_smiles(PHENOL))) < 1.0


def test_translate_to_atom_ids_uses_builder_ids():
    builder = smsd.MolGraphBuilder()
    builder.atom_count(2).atomic_numbers([6, 8]).formal_charges([0, 0]).neighbors([[1], [0]])
    builder.bond_orders([[1], [1]]).atom_ids([10, 20])
    graph = builder.build(False)
    plain = smsd.parse_smiles("CO")
    assert smsd.translate_to_atom_ids({0: 0, 1: 1}, graph, graph) == {10: 10, 20: 20}
    assert smsd.translate_to_atom_ids({0: 0, 1: 1}, plain, plain) == {0: 0, 1: 1}


def test_layout_measures_and_overlap_resolution():
    mol = smsd.parse_smiles("c1ccc2ccccc2c1")
    coords = smsd.generate_coords_2d(mol)
    assert smsd.layout_quality(mol, coords) >= 0.0
    assert not smsd.is_degenerate_layout(coords)
    assert smsd.is_degenerate_layout([(0.0, 0.0), (1.0, 0.0), (2.0, 0.0)])
    remaining, resolved = smsd.resolve_overlaps([(0.0, 0.0), (0.0, 0.0)])
    assert remaining == 0 and math.dist(*resolved) > 0.0
    three_d = smsd.generate_coords_3d(mol)
    assert len(three_d) == mol.n and all(len(point) == 3 for point in three_d)


def test_consensus_result_fields():
    result = smsd.MCS.consensus(["N", "CCO", "CCN"], threshold=2 / 3)
    assert isinstance(result, smsd.ConsensusResult)
    assert result.found and result.support_count >= 2
    assert isinstance(result.to_smiles(), str) and result.to_smiles()
    assert set(result.mappings) <= {0, 1, 2}
