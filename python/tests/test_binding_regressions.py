# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
"""Python lifetime, option and atom-order regression checks."""
import gc
import weakref

import pytest
import smsd


def test_match_many_timeout_and_worker_count_are_separate():
    prepared = smsd.Substructure.find("C").with_timeout(37)
    results = prepared.match_many(["CC"], num_threads=2)
    assert len(results) == 1 and results[0].found
    assert results[0].status == smsd.SearchStatus.COMPLETED
    immediate = prepared.with_timeout(0).match_many(["CC"], num_threads=2)
    assert immediate[0].status == smsd.SearchStatus.TIMED_OUT
    assert not immediate[0].complete


def test_mces_preserves_flexible_aromatic_bonds():
    query = smsd.parse_smiles("O=C(Cc1ccccc1)c1ccncc1")
    target = smsd.parse_smiles("O/C(/c1ccncc1)=C\\c1ccccc1")
    options = smsd.MCS.Options()
    options.objective = smsd.MCS.Objective.MCES
    options.chemistry.match_bond_order = smsd.BondOrderMode.STRICT
    options.chemistry.aromaticity_mode = smsd.AromaticityMode.FLEXIBLE
    prepared = smsd.MCS.find(query).with_options(options).with_timeout(1000)
    for result in [prepared.match(target), *prepared.match_many([target])]:
        assert result.best.atom_count >= 10
        assert smsd.MCS.validate_mapping(query, target, result.best, options) == []


def test_match_many_options_and_signed_weights_match_single_search():
    targets = ["CCC", "CCN", "C1CC1"]
    prepared = smsd.MCS.find("CCC").with_atom_weights([10.0, -30.0, 1.0]).with_timeout(1000)
    batch = prepared.match_many(targets, num_threads=2)
    single = [prepared.match(target) for target in targets]
    assert [result.best.atom_count for result in batch] == [result.best.atom_count for result in single] == [1, 1, 1]
    assert all(set(result.best.atom_mapping) == {0} for result in batch)
    options = smsd.MCS.Options()
    options.chemistry.match_atom_type = False
    assert dict(smsd.MCS.find("C").with_options(options).match_many(["N"])[0].best.atom_mapping) == {0: 0}


@pytest.mark.parametrize("weights", [[1.0], [float("nan"), 1.0], [float("inf"), 1.0]])
def test_invalid_weights_raise_before_parallel_execution(weights):
    with pytest.raises(ValueError):
        smsd.MCS.find("CC").with_atom_weights(weights).match_many(["CCC"] * 4, num_threads=2)
    options = smsd.MCS.Options()
    options.atom_weights = weights
    with pytest.raises(ValueError):
        smsd.MCS.find("CC").with_options(options).match("CCC")


def test_repeated_graphs_remain_usable_in_parallel_searches():
    query = smsd.parse_smiles("C1CC1")
    target = smsd.parse_smiles("CC1CC1")
    prepared = smsd.MCS.find(query)
    for _ in range(3):
        assert [result.best.atom_count for result in prepared.match_many([target] * 32, num_threads=2)] == [3] * 32
    smsd.prewarm(query)
    smsd.prewarm(target)
    assert len(query.canonical_label) == len(query)
    assert len(target.canonical_label) == len(target)
    assert [result.found for result in smsd.Substructure.smarts("C1CC1").match_many([target] * 4)] == [True] * 4
    with pytest.raises((TypeError, RuntimeError)):
        prepared.match_many([target, None])


def test_search_preserves_options_and_rejects_unsupported_callbacks():
    result = smsd.MCS.find("CCC").with_atom_weights([10, -30, 1]).match("CCC")
    assert result.best.atom_count == 1 and set(result.best.atom_mapping) == {0}
    assert result.nodes_visited >= 0 and result.elapsed_ms >= 0
    token = smsd.CancellationToken()
    options = smsd.MCS.Options()
    options.cancellation = token
    prepared = smsd.MCS.find("CC").with_options(options)
    token.cancel()
    assert prepared.match("CC").status == smsd.SearchStatus.CANCELLED
    with pytest.raises(TypeError):
        prepared.match("CC", on_progress=lambda *_: None)


def _rdkit():
    return pytest.importorskip("rdkit.Chem")


def _reordered_pair():
    chem = _rdkit()
    q = chem.MolFromSmiles("NCCO")
    t = chem.MolFromSmiles("CC(O)CN")
    return chem.RenumberAtoms(q, [3, 1, 0, 2]), chem.RenumberAtoms(t, [4, 1, 3, 0, 2])


def _assert_witness(q, t, mapping):
    assert len(set(mapping.values())) == len(mapping)
    for a, b in mapping.items():
        assert q.GetAtomWithIdx(a).GetAtomicNum() == t.GetAtomWithIdx(b).GetAtomicNum()
    for bond in q.GetBonds():
        a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if a in mapping and b in mapping:
            matched = t.GetBondBetweenAtoms(mapping[a], mapping[b])
            assert matched is not None
            assert matched.GetBondType() == bond.GetBondType()


@pytest.mark.parametrize("input_kind", ["molecule", "graph", "repeated"])
def test_reordered_rdkit_indices_for_each_input_kind(input_kind):
    query, target = _reordered_pair()
    prepared = smsd.MCS.find(query)
    if input_kind == "graph":
        # The source molecule conversion supplies an original-index map.
        query_graph, target_graph = smsd.from_rdkit(query), smsd.from_rdkit(target)
        native = smsd.MCS.find(query_graph).match(target_graph)
        mapping = smsd.translate_mapping(dict(native.best.atom_mapping), query_graph, target_graph)
    else:
        result = prepared.match_many([target])[0] if input_kind == "repeated" else prepared.match(target)
        mapping = result.best.atom_mapping
    assert len(mapping) == 4
    _assert_witness(query, target, mapping)


def test_rdkit_first_all_and_many_indices():
    query, target = _reordered_pair()
    prepared_mcs = smsd.MCS.find(query)
    prepared_substructure = smsd.Substructure.find(query)
    results = [prepared_mcs.match(target), prepared_substructure.match(target),
               prepared_mcs.match_all(target, 4), prepared_substructure.match_all(target, 4),
               *prepared_mcs.match_many([target, target]), *prepared_substructure.match_many([target, target])]
    for result in results:
        assert result.found
        for match in result.matches:
            assert match.atom_count == 4
            _assert_witness(query, target, match.atom_mapping)
    chem = _rdkit()
    options = smsd.MCS.Options()
    options.chemistry.match_atom_type = False
    result = smsd.MCS.find(chem.MolFromSmiles("C")).with_options(options).match(chem.MolFromSmiles("N"))
    assert dict(result.best.atom_mapping) == {0: 0}


def test_rdkit_conversion_cache_mutation_and_lifetime():
    chem = _rdkit()
    smsd.clear_cache()
    mol = chem.RWMol(chem.MolFromSmiles("CO"))
    graph = smsd.from_rdkit(mol)
    assert smsd.from_rdkit(mol) is graph
    index_map = smsd.get_index_map(graph)
    index_map[0] = 99
    assert 99 not in smsd.get_index_map(graph)
    mol.GetAtomWithIdx(1).SetAtomicNum(7)
    mol.UpdatePropertyCache()
    changed = smsd.from_rdkit(mol)
    assert changed is not graph
    assert sorted(changed.atomic_num) == [6, 7]
    smsd.clear_cache()
    assert sorted(smsd.get_index_map(graph)) == [0, 1]
    molecule_ref = weakref.ref(mol)
    graph_ref = weakref.ref(changed)
    del mol, changed
    gc.collect()
    assert molecule_ref() is None
    assert graph_ref() is None


def test_target_selection_translates_original_atom_indices():
    chem = _rdkit()
    query = chem.MolFromSmiles("NO")
    targets = [chem.MolFromSmiles("C"), chem.RenumberAtoms(chem.MolFromSmiles("CNO"), [2, 0, 1])]
    results = smsd.MCS.find(query).match_many(targets)
    target_index, result = max(enumerate(results), key=lambda item: item[1].score)
    assert target_index == 1 and result.best.atom_count == 2
    _assert_witness(query, targets[target_index], result.best.atom_mapping)


@pytest.mark.parametrize("smiles", ["N[C@@H](C)C(=O)O", "N[C@H]1CCCCO1", "F[C@]1(Cl)CCCCO1", "CC1.[C@H]1(N)O"])
def test_stereo_survives_random_smiles_traversal(smiles):
    chem = _rdkit()
    mol = chem.MolFromSmiles(smiles)
    chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    for seed in range(8):
        random_smiles = chem.MolToRandomSmilesVect(mol, 1, randomSeed=seed + 1)[0]
        other = chem.MolFromSmiles(random_smiles)
        options = smsd.MCS.Options()
        options.chemistry.use_chirality = True
        assert smsd.MCS.find(smiles).with_options(options).match(random_smiles).best.atom_count == mol.GetNumAtoms()
        graph = smsd.from_rdkit(other, use_cache=False)
        index_map = smsd.get_index_map(graph)
        for g_index, rd_index in enumerate(index_map):
            atom = other.GetAtomWithIdx(rd_index)
            if atom.HasProp("_CIPCode"):
                assert smsd.assign_rs(graph)[g_index] == atom.GetProp("_CIPCode")


def test_rdkit_weights_follow_original_atom_order():
    chem = _rdkit()
    query = chem.RenumberAtoms(chem.MolFromSmiles("CCC"), [1, 0, 2])
    target = chem.MolFromSmiles("CCC")
    prepared = smsd.MCS.find(query).with_atom_weights([10, -30, 1]).with_timeout(1000)
    results = [prepared.match(target), prepared.match_all(target, 4), *prepared.match_many([target, target])]
    assert all(set(match.atom_mapping) == {0, 2} for result in results for match in result.matches)
    assert all(result.best.atom_count == 2 for result in results)
    assert all(result.score == 11 for result in results)


def test_rdkit_explicit_hydrogen_indices_are_retained():
    chem = _rdkit()
    molecule = chem.AddHs(chem.MolFromSmiles("CO"))
    graph = smsd.from_rdkit(molecule)
    assert len(graph) == molecule.GetNumAtoms()
    assert sorted(smsd.get_index_map(graph)) == list(range(molecule.GetNumAtoms()))
    mapping = smsd.MCS.find(molecule).match(molecule).best.atom_mapping
    assert len(mapping) == molecule.GetNumAtoms()
    _assert_witness(molecule, molecule, mapping)


@pytest.mark.parametrize("function", ["path_fingerprint", "mcs_fingerprint"])
@pytest.mark.parametrize("option", [{"fp_size": 0}, {"fp_size": -64}, {"path_length": 0}, {"path_length": -1}])
def test_raw_fingerprint_parameters_reject_invalid_storage(function, option):
    with pytest.raises(ValueError):
        getattr(smsd._smsd, function)(smsd.parse_smiles("CC"), **option)


@pytest.mark.parametrize("option", [{"fp_size": 0}, {"fp_size": -64}, {"path_length": 0}, {"path_length": -1}])
def test_batch_fingerprint_validates_before_workers(option):
    with pytest.raises(ValueError):
        smsd._smsd.batch_fingerprint([smsd.parse_smiles("CC")] * 4, num_threads=2, **option)
