# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
# Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
# See the NOTICE file for attribution, trademark, and algorithm IP terms.
"""Checks of mixed input indices and callback ownership."""

import pytest
import smsd
from smsd._interop import _compute_index_map
from .search_helpers import assert_witness, atom_map, configured_mcs


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


@pytest.mark.parametrize("graph_side", ["query", "target"])
@pytest.mark.parametrize("entrypoint", ["single", "multiple", "substructure",
                                         "all_substructures", "batch",
                                         "batch_substructure", "corpus",
                                         "corpus_substructure"])
def test_mixed_inputs_keep_graph_atom_indices(graph_side, entrypoint):
    query, target = _mixed_pair(graph_side)
    if entrypoint == "single":
        mappings = [atom_map(configured_mcs(query).match(target))]
    elif entrypoint == "multiple":
        mappings = [dict(match.atom_mapping) for match in configured_mcs(query).match_all(target, 4).matches]
    elif entrypoint == "substructure":
        mappings = [atom_map(smsd.Substructure.find(query).match(target))]
    elif entrypoint == "all_substructures":
        mappings = [dict(match.atom_mapping) for match in smsd.Substructure.find(query).match_all(target, 4).matches]
    elif entrypoint == "batch":
        mappings = [atom_map(result) for result in smsd.MCS.find(query).match_many([target])]
    elif entrypoint == "batch_substructure":
        mappings = [atom_map(result) for result in smsd.Substructure.find(query).match_many([target])]
    else:
        corpus = smsd.TargetCorpus.from_smiles([target])
        search = smsd.MCS.find(query) if entrypoint == "corpus" else smsd.Substructure.find(query)
        mappings = [atom_map(result) for result in search.match_many(corpus)]
    assert mappings
    for mapping in mappings:
        assert_witness(query, target, mapping, size=4)


def test_explicit_hydrogen_index_order_uses_rdkit_output_metadata(monkeypatch):
    chem = pytest.importorskip("rdkit.Chem")
    molecule = chem.AddHs(chem.MolFromSmiles("N[C@@H](C)C(=O)O"))
    smiles = chem.MolToSmiles(molecule)
    expected = list(molecule.GetPropsAsDict(includePrivate=True, includeComputed=True)
                    ["_smilesAtomOutputOrder"])

    def unnecessary_round_trip(*args, **kwargs):
        raise AssertionError("Valid output-order metadata must avoid hydrogen-removing reparse")

    monkeypatch.setattr(chem, "MolFromSmiles", unnecessary_round_trip)
    assert _compute_index_map(molecule, smiles) == expected
    assert sorted(expected) == list(range(molecule.GetNumAtoms()))


def test_rdkit_weight_translation_preserves_the_caller_list():
    chem = pytest.importorskip("rdkit.Chem")
    molecule = chem.RenumberAtoms(chem.MolFromSmiles("NCCO"), [3, 1, 0, 2])
    graph = smsd.from_rdkit(molecule)
    weights = [10.0, 20.0, 30.0, 40.0]
    query = smsd.MCS.find(molecule).with_atom_weights(weights)
    result = query.with_timeout(1000).match("N")
    nitrogen = next(atom.GetIdx() for atom in molecule.GetAtoms() if atom.GetAtomicNum() == 7)
    assert dict(result.best.atom_mapping) == {nitrogen: 0}
    assert result.score == weights[nitrogen]
    assert weights == [10.0, 20.0, 30.0, 40.0]


def test_returned_mapping_cannot_mutate_a_reused_query():
    query = smsd.MCS.find("NCCO")
    result = query.match("NCCO")
    with pytest.raises(TypeError):
        result.best.atom_mapping[99] = 99
    copied = dict(result.best.atom_mapping)
    copied.clear()
    assert dict(query.match("NCCO").best.atom_mapping) == {0: 0, 1: 1, 2: 2, 3: 3}
    assert result.elapsed_ms >= 0 and result.optimal
