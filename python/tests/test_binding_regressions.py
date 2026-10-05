# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
"""Python lifetime, option and atom-order regression checks."""
import gc
import weakref

import pytest
import smsd


def test_batch_timeout_and_worker_count_are_separate(monkeypatch):
    seen = []
    def native(query, targets, chemistry, workers, timeout):
        seen.append((workers, timeout))
        return [True] * len(targets)
    monkeypatch.setattr(smsd._smsd, "batch_substructure", native)
    assert smsd.batch_substructure("C", ["CC"], timeout_ms=37, num_threads=2) == [True]
    assert seen == [(2, 37)]


def test_coverage_preserves_flexible_aromatic_bonds():
    query = smsd.parse_smiles("O=C(Cc1ccccc1)c1ccncc1")
    target = smsd.parse_smiles("O/C(/c1ccncc1)=C\\c1ccccc1")
    chemistry = smsd.ChemOptions()
    chemistry.match_bond_order = smsd.BondOrderMode.STRICT
    chemistry.aromaticity_mode = smsd.AromaticityMode.FLEXIBLE
    mappings = [smsd._smsd.find_mcs_coverage(query, target, timeout_ms=1000),
                smsd.find_mcs(query, target, strategy="auto",
                              match_bond_order="strict", timeout_ms=1000)]
    for mapping in mappings:
        assert len(mapping) >= 10
        assert smsd._smsd.validate_mapping(query, target, mapping, chemistry) == []


def test_batch_options_and_signed_weights_match_single_search():
    targets = ["CCC", "CCN", "C1CC1"]
    options = dict(atom_weights=[10.0, -30.0, 1.0], connected_only=True,
                   induced=True, timeout_ms=1000)
    batch = smsd.batch_mcs("CCC", targets, num_threads=2, **options)
    single = [smsd.find_mcs("CCC", t, strategy="native", **options) for t in targets]
    assert [len(m) for m in batch] == [len(m) for m in single] == [1, 1, 1]
    assert all(set(m) == {0} for m in batch)
    assert smsd.batch_mcs_size("CCC", targets, num_threads=2, **options) == [1, 1, 1]
    assert smsd.batch_mcs("C", ["N"], match_atom_type=False) == [{0: 0}]


@pytest.mark.parametrize("weights", [[1.0], [float("nan"), 1.0], [float("inf"), 1.0]])
def test_invalid_batch_weights_raise_before_parallel_execution(weights):
    with pytest.raises(ValueError):
        smsd.batch_mcs("CC", ["CCC"] * 4, atom_weights=weights, num_threads=2)
    with pytest.raises(ValueError):
        smsd.batch_mcs_size("CC", ["CCC"] * 4, atom_weights=weights, num_threads=2)


def test_repeated_graphs_remain_usable_in_parallel_batches():
    q = smsd.parse_smiles("C1CC1")
    t = smsd.parse_smiles("CC1CC1")
    for _ in range(3):
        assert smsd._smsd.batch_mcs_size(q, [t] * 32, num_threads=2) == [3] * 32
    assert len(q.canonical_label) == len(q)
    assert len(t.canonical_label) == len(t)
    assert smsd.compile_smarts("C1CC1").matches_many([t] * 4) == [True] * 4
    with pytest.raises((TypeError, RuntimeError)):
        smsd._smsd.batch_mcs(q, [t, None])


def test_progressive_preserves_options_and_callback_exceptions():
    seen = []
    mapping = smsd.find_mcs_progressive("CCC", "CCC", atom_weights=[10, -30, 1],
                                       on_progress=lambda m, n, ms: seen.append((m, n, ms)))
    assert len(mapping) == 1 and set(mapping) == {0}
    assert seen[-1][0] == mapping
    assert all(n == len(m) and ms >= 0 for m, n, ms in seen)
    def failing(*args):
        raise LookupError("callback failed")
    with pytest.raises(LookupError, match="callback failed"):
        smsd.find_mcs_progressive("CC", "CC", on_progress=failing)


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


@pytest.mark.parametrize("strategy", ["native", "lightweight", "auto"])
def test_reordered_rdkit_indices_for_each_strategy(strategy):
    q, t = _reordered_pair()
    result = smsd.find_mcs(q, t, strategy=strategy)
    assert len(result) == 4
    _assert_witness(q, t, result)


def test_rdkit_single_multiple_batch_and_progressive_indices():
    q, t = _reordered_pair()
    for result in [smsd.mcs_rdkit(q, t), smsd.mcs_rdkit_native(q, t),
                   smsd.substructure_rdkit(q, t), smsd.find_mcs_progressive(q, t),
                   *smsd.find_mcs(q, t, max_results=4),
                   *smsd.find_substructure(q, t, max_results=4),
                   *smsd.batch_mcs(q, [t, t]), *smsd.batch_mcs_rdkit(q, [t, t])]:
        assert len(result) == 4
        _assert_witness(q, t, result)
    for pairs in smsd.batch_find_substructure(q, [t, t]):
        _assert_witness(q, t, dict(pairs))
    chem = _rdkit()
    assert smsd.mcs_rdkit_native(chem.MolFromSmiles("C"), chem.MolFromSmiles("N"),
                                match_atom_type=False) == {0: 0}


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


def test_constrained_batch_translates_the_selected_target():
    chem = _rdkit()
    query = chem.MolFromSmiles("NO")
    targets = [chem.MolFromSmiles("C"), chem.RenumberAtoms(chem.MolFromSmiles("CNO"), [2, 0, 1])]
    results = smsd.batch_mcs_constrained([query], targets, return_target_indices=True)
    target_index, mapping = results[0]
    assert target_index == 1 and len(mapping) == 2
    _assert_witness(query, targets[target_index], mapping)
    assert smsd.batch_mcs_constrained([query], targets) == [mapping]


@pytest.mark.parametrize("smiles", ["N[C@@H](C)C(=O)O", "N[C@H]1CCCCO1", "F[C@]1(Cl)CCCCO1", "CC1.[C@H]1(N)O"])
def test_stereo_survives_random_smiles_traversal(smiles):
    chem = _rdkit()
    mol = chem.MolFromSmiles(smiles)
    chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    for seed in range(8):
        random_smiles = chem.MolToRandomSmilesVect(mol, 1, randomSeed=seed + 1)[0]
        other = chem.MolFromSmiles(random_smiles)
        assert len(smsd.find_mcs(smiles, random_smiles, use_chirality=True, strategy="native")) == mol.GetNumAtoms()
        graph = smsd.from_rdkit(other, use_cache=False)
        index_map = smsd.get_index_map(graph)
        for g_index, rd_index in enumerate(index_map):
            atom = other.GetAtomWithIdx(rd_index)
            if atom.HasProp("_CIPCode"):
                assert smsd.assign_rs(graph)[g_index] == atom.GetProp("_CIPCode")


def test_rdkit_weights_follow_original_atom_order():
    chem = _rdkit()
    q = chem.RenumberAtoms(chem.MolFromSmiles("CCC"), [1, 0, 2])
    t = chem.MolFromSmiles("CCC")
    options = dict(atom_weights=[10, -30, 1], timeout_ms=1000)
    results = [smsd.find_mcs(q, t, strategy="native", **options),
               *smsd.find_mcs(q, t, max_results=4, **options),
               smsd.find_mcs_progressive(q, t, **options),
               *smsd.batch_mcs(q, [t], **options),
               *smsd.batch_mcs_constrained([q], [t], **options)]
    assert all(set(mapping) == {0, 2} for mapping in results)
    assert smsd.batch_mcs_size(q, [t], **options) == [2]
    with pytest.raises(ValueError, match="require one query"):
        smsd.batch_mcs_constrained([q, q], [t], **options)


def test_rdkit_explicit_hydrogen_indices_are_retained():
    chem = _rdkit()
    molecule = chem.AddHs(chem.MolFromSmiles("CO"))
    graph = smsd.from_rdkit(molecule)
    assert len(graph) == molecule.GetNumAtoms()
    assert sorted(smsd.get_index_map(graph)) == list(range(molecule.GetNumAtoms()))
    mapping = smsd.find_mcs(molecule, molecule, strategy="native")
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
