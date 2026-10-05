# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
"""Independent checks of mixed input indices and callback ownership."""

from types import SimpleNamespace

import pytest
import smsd


def _mixed_pair(graph_side):
    chem = pytest.importorskip("rdkit.Chem")
    query = chem.RenumberAtoms(chem.MolFromSmiles("NCCO"), [3, 1, 0, 2])
    target = chem.RenumberAtoms(chem.MolFromSmiles("CC(O)CN"), [4, 1, 3, 0, 2])
    if graph_side == "query":
        query = smsd.from_rdkit(query)
        assert smsd.get_index_map(query) != list(range(len(query)))
    else:
        target = smsd.from_rdkit(target)
        assert smsd.get_index_map(target) != list(range(len(target)))
    return query, target


def _atomic_numbers(molecule):
    if isinstance(molecule, smsd.MolGraph):
        return list(molecule.atomic_num)
    return [atom.GetAtomicNum() for atom in molecule.GetAtoms()]


def _bonds(molecule):
    if isinstance(molecule, smsd.MolGraph):
        return {(a, b): molecule.bond_order(a, b)
                for a in range(len(molecule)) for b in range(a + 1, len(molecule))
                if molecule.bond_order(a, b) > 0}
    return {(min(b.GetBeginAtomIdx(), b.GetEndAtomIdx()),
             max(b.GetBeginAtomIdx(), b.GetEndAtomIdx())): int(b.GetBondTypeAsDouble())
            for b in molecule.GetBonds()}


def _assert_witness(query, target, mapping):
    query_atoms, target_atoms = _atomic_numbers(query), _atomic_numbers(target)
    assert len(mapping) == len(query_atoms) == 4
    assert len(set(mapping.values())) == len(mapping)
    for a, b in mapping.items():
        assert query_atoms[a] == target_atoms[b]
    target_bonds = _bonds(target)
    for (a, b), order in _bonds(query).items():
        mapped_pair = tuple(sorted((mapping[a], mapping[b])))
        assert target_bonds.get(mapped_pair) == order


@pytest.mark.parametrize("graph_side", ["query", "target"])
@pytest.mark.parametrize("entrypoint", ["single", "multiple", "substructure",
                                         "all_substructures", "batch",
                                         "batch_substructure", "progressive",
                                         "constrained"])
def test_mixed_inputs_keep_graph_atom_indices(graph_side, entrypoint):
    query, target = _mixed_pair(graph_side)
    if entrypoint == "single":
        mappings = [smsd.find_mcs(query, target, strategy="native")]
    elif entrypoint == "multiple":
        mappings = smsd.find_mcs(query, target, max_results=4)
    elif entrypoint == "substructure":
        mappings = [smsd.find_substructure(query, target)]
    elif entrypoint == "all_substructures":
        mappings = smsd.find_substructure(query, target, max_results=4)
    elif entrypoint == "batch":
        mappings = smsd.batch_mcs(query, [target])
    elif entrypoint == "batch_substructure":
        mappings = [dict(pairs) for pairs in smsd.batch_find_substructure(query, [target])]
    elif entrypoint == "progressive":
        callbacks = []
        result = smsd.find_mcs_progressive(
            query, target, on_progress=lambda mapping, size, elapsed: callbacks.append(mapping))
        mappings = [result, *callbacks]
        assert callbacks
    else:
        selected = smsd.batch_mcs_constrained([query], [target], return_target_indices=True)
        assert selected[0][0] == 0
        mappings = [selected[0][1]]
    assert mappings
    for mapping in mappings:
        _assert_witness(query, target, mapping)


def test_explicit_hydrogen_index_order_uses_rdkit_output_metadata(monkeypatch):
    chem = pytest.importorskip("rdkit.Chem")
    molecule = chem.AddHs(chem.MolFromSmiles("N[C@@H](C)C(=O)O"))
    smiles = chem.MolToSmiles(molecule)
    expected = list(molecule.GetPropsAsDict(includePrivate=True, includeComputed=True)
                    ["_smilesAtomOutputOrder"])

    def unnecessary_round_trip(*args, **kwargs):
        raise AssertionError("Valid output-order metadata must avoid hydrogen-removing reparse")

    monkeypatch.setattr(chem, "MolFromSmiles", unnecessary_round_trip)
    assert smsd._compute_index_map(molecule, smiles) == expected
    assert sorted(expected) == list(range(molecule.GetNumAtoms()))


def test_rdkit_weight_translation_preserves_the_caller_list():
    chem = pytest.importorskip("rdkit.Chem")
    molecule = chem.RenumberAtoms(chem.MolFromSmiles("NCCO"), [3, 1, 0, 2])
    graph = smsd.from_rdkit(molecule)
    weights = [10.0, 20.0, 30.0, 40.0]
    options = SimpleNamespace(atom_weights=weights)
    smsd._remap_query_weights(options, graph, molecule)
    assert options.atom_weights == [weights[index] for index in smsd.get_index_map(graph)]
    assert weights == [10.0, 20.0, 30.0, 40.0]


def test_progress_callback_mapping_is_independent_of_returned_mapping():
    calls = []

    def mutate_callback(mapping, size, elapsed):
        calls.append((size, elapsed))
        assert size == len(mapping)
        mapping.clear()
        mapping[99] = 99

    result = smsd.find_mcs_progressive("NCCO", "NCCO", on_progress=mutate_callback)
    assert len(calls) == 1
    assert calls[0][0] == 4 and calls[0][1] >= 0
    assert result == {0: 0, 1: 1, 2: 2, 3: 3}
