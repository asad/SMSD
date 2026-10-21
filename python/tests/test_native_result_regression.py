# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
# Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
# See the NOTICE file for attribution, trademark, and algorithm IP terms.
"""Search results built natively equal those built in Python from the same fields."""
import dataclasses
from types import MappingProxyType

import pytest
import smsd
from smsd import _smsd

MATCH_FIELDS = ["atom_mapping", "bond_mapping", "fragments"]
RESULT_FIELDS = ["matches", "status", "optimal", "complete", "score", "nodes_visited", "elapsed_ms"]


def _mces(query):
    return smsd.MCS.find(query).with_objective(smsd.MCS.Objective.MCES)


def test_status_members_are_the_shared_enumeration_members():
    assert smsd.MCS.find("CCO").match("CCN").status is smsd.SearchStatus.COMPLETED
    assert smsd.Substructure.find("C").match_all("CC").status is smsd.SearchStatus.COMPLETED
    assert smsd.Substructure.find("C").with_timeout(0).match("CC").status is smsd.SearchStatus.TIMED_OUT
    for result in smsd.Substructure.find("C").match_many(["CC", "N"]):
        assert result.status is smsd.SearchStatus.COMPLETED


def test_native_results_equal_those_built_in_python():
    result = _mces("CC(=O)Oc1ccccc1C(=O)O").match("CC(=O)Nc1ccc(O)cc1")
    match = result.best
    assert type(result) is smsd.SearchResult and type(match) is smsd.Match
    assert isinstance(match.atom_mapping, MappingProxyType)
    assert all(type(bond) is tuple and len(bond) == 4 for bond in match.bond_mapping)
    assert all(type(fragment) is tuple for fragment in match.fragments)
    rebuilt = smsd.Match(MappingProxyType(dict(match.atom_mapping)),
                         tuple(tuple(bond) for bond in match.bond_mapping),
                         tuple(tuple(fragment) for fragment in match.fragments))
    assert match == rebuilt and rebuilt == match
    assert (match.atom_count, match.bond_count) == (rebuilt.atom_count, rebuilt.bond_count)
    assert result == smsd.SearchResult((rebuilt,), smsd.SearchStatus.COMPLETED, result.optimal,
                                       result.complete, result.score, result.nodes_visited, result.elapsed_ms)
    assert type(result.score) is float and type(result.nodes_visited) is int and type(result.elapsed_ms) is int
    assert result.optimal is True and result.complete is True


def test_native_results_are_frozen():
    result = smsd.MCS.find("CCO").match("CCN")
    with pytest.raises(dataclasses.FrozenInstanceError):
        result.status = smsd.SearchStatus.CANCELLED
    with pytest.raises(dataclasses.FrozenInstanceError):
        result.best.atom_mapping = {}
    with pytest.raises(TypeError):
        result.best.atom_mapping[0] = 0
    with pytest.raises(AttributeError):
        result.matches.append(None)


def test_native_results_set_exactly_the_dataclass_fields():
    assert [field.name for field in dataclasses.fields(smsd.Match)] == MATCH_FIELDS
    assert [field.name for field in dataclasses.fields(smsd.SearchResult)] == RESULT_FIELDS
    result = smsd.MCS.find("CCO").match("CCN")
    assert list(vars(result)) == RESULT_FIELDS
    assert list(vars(result.best)) == MATCH_FIELDS


def _dataclass(name, fields, **namespace):
    return dataclasses.make_dataclass(name, [(field, object) for field in fields], namespace=namespace, frozen=True)


@pytest.mark.parametrize("match,result", [
    (_dataclass("Renamed", ["atom_mapping", "bonds", "fragments"]), smsd.SearchResult),
    (_dataclass("Extra", MATCH_FIELDS + ["extra"]), smsd.SearchResult),
    (smsd.Match, _dataclass("Missing", RESULT_FIELDS[:-1])),
    (_dataclass("OwnNew", MATCH_FIELDS, __new__=lambda cls, *args: object.__new__(cls)), smsd.SearchResult),
    (type("Plain", (), {}), smsd.SearchResult),
    (smsd.Match, "not a class"),
])
def test_registration_rejects_classes_the_native_results_cannot_fill(match, result):
    with pytest.raises(TypeError):
        _smsd._register_search_types(match, result)
    # A rejected registration keeps the classes registered before it.
    assert type(smsd.MCS.find("CC").match("CC")) is smsd.SearchResult


def test_atom_numbering_must_hold_exact_ints():
    prepared = smsd.MCS.find("CCO")
    graph = smsd.parse_smiles("CCO")
    native = prepared._native.match(graph)
    assert dict(_smsd._search_result(native, (10, 11, 12), [20, 21, 22]).best.atom_mapping) == {
        10: 20, 11: 21, 12: 22}

    class Number(int):
        pass

    for numbering in [(10.0, 11, 12), (True, 1, 2), [Number(0), 1, 2], ("0", "1", "2")]:
        with pytest.raises(TypeError):
            _smsd._search_result(native, numbering, None)
        with pytest.raises(TypeError):
            _smsd._search_result(native, None, numbering)
        with pytest.raises(TypeError):
            prepared._native._match(graph, numbering, None)
        with pytest.raises(TypeError):
            prepared._native._match_many([(graph, numbering)], 1, None)
        with pytest.raises(TypeError):
            prepared._native._with_numbering(numbering)
        with pytest.raises(TypeError):
            prepared._search._match_many([(graph, numbering)], 1)
        with pytest.raises(TypeError):
            _smsd._python_match(native.matches[0], numbering, None)


@pytest.mark.parametrize("factory", [smsd.MCS.find, smsd.Substructure.find])
def test_numbered_queries_report_the_query_numbering(factory):
    prepared = factory("CCO")
    graph = smsd.parse_smiles("CCO")
    numbered = prepared._native._with_numbering((10, 11, 12))
    renumbered = {10: 0, 11: 1, 12: 2}
    assert dict(numbered._match(graph).best.atom_mapping) == renumbered
    assert dict(numbered._match_all(graph, 10).best.atom_mapping) == renumbered
    batch = numbered._match_many([(graph, None), (graph, (20, 21, 22))], 2)
    assert [dict(result.best.atom_mapping) for result in batch] == [renumbered, {10: 20, 11: 21, 12: 22}]
    old = prepared._native._match(graph, (10, 11, 12), None)
    new = numbered._match(graph)
    assert (new.matches, new.status, new.score, new.nodes_visited) == (old.matches, old.status, old.score,
                                                                     old.nodes_visited)
    native = prepared._native._with_numbering(None)
    assert dict(native._match(graph).best.atom_mapping) == dict(prepared.match(graph).best.atom_mapping) == {
        0: 0, 1: 1, 2: 2}


def test_changed_options_reach_graph_target_searches():
    graph = smsd.parse_smiles("CC")
    prepared = smsd.Substructure.find("C")
    stopped = prepared.with_timeout(0)
    assert stopped.match(graph).status is smsd.SearchStatus.TIMED_OUT
    assert stopped.match_all(graph).status is smsd.SearchStatus.TIMED_OUT
    assert stopped.match_many([graph])[0].status is smsd.SearchStatus.TIMED_OUT
    assert prepared.match(graph).status is smsd.SearchStatus.COMPLETED
    assert prepared.match_many([graph])[0].status is smsd.SearchStatus.COMPLETED


@pytest.mark.parametrize("factory", [smsd.MCS.find, smsd.Substructure.find])
def test_molecule_queries_keep_their_numbering_against_graphs(factory):
    chemistry = pytest.importorskip("rdkit.Chem")
    query = chemistry.RenumberAtoms(chemistry.MolFromSmiles("NCCO"), [3, 1, 0, 2])
    graph = smsd.parse_smiles("NCCO")
    prepared = factory(query)
    expected = prepared._native._match(graph, prepared._query_indices, None)
    for result in [prepared.match(graph), prepared.match_all(graph), *prepared.match_many([graph, graph])]:
        assert result.best.atom_mapping == expected.best.atom_mapping
        for query_atom, target_atom in result.best.atom_mapping.items():
            assert query.GetAtomWithIdx(query_atom).GetAtomicNum() == graph.atomic_num[target_atom]


class _FreshItems:
    """A sequence that creates each item afresh, so that only its caller holds a reference."""

    def __init__(self, smiles, pairs):
        self._smiles = smiles
        self._pairs = pairs

    def __len__(self):
        return len(self._smiles)

    def __getitem__(self, index):
        if not 0 <= index < len(self._smiles):
            raise IndexError(index)
        graph = smsd.parse_smiles(self._smiles[index])
        return (graph, tuple(range(100, 100 + graph.n))) if self._pairs else graph


def test_sequences_that_create_their_items_afresh():
    smiles = ["C", "CC", "CCO", "CCCC", "OCCO", "CCCCCC"] * 4
    prepared = _mces("CCO")
    expected = [result.best.atom_count if result.found else 0 for result in map(prepared.match, smiles)]
    assert expected[:6] == [0, 2, 3, 2, 3, 2]
    for workers in [1, 4]:
        native = prepared._native.match_many(_FreshItems(smiles, False), workers)
        assert [len(result.matches[0].atom_mapping) if result.matches else 0 for result in native] == expected
        results = prepared._native._match_many(_FreshItems(smiles, True), workers, None)
        assert [result.best.atom_count if result.found else 0 for result in results] == expected
        for result in results:
            assert all(atom >= 100 for atom in (result.best.atom_mapping.values() if result.found else ()))


@pytest.mark.parametrize("smiles,path_length,fp_size,bits", [
    ("c1ccccc1O", 7, 2048, [145, 147, 273, 275, 277, 279, 401, 403, 405, 407, 533, 535, 1815]),
    ("CC(=O)Nc1ccc(O)cc1", 7, 2048,
     [145, 146, 147, 150, 180, 273, 274, 275, 277, 278, 279, 304, 308, 401, 402, 403, 405, 406, 407,
      432, 436, 533, 534, 535, 560, 562, 690, 1110, 1236, 1238, 1240, 1362, 1364, 1368, 1490, 1492,
      1688, 1815, 1846, 1938]),
    ("C1CC2CCC1C2", 3, 64, [17, 21]),
    ("[NH3+]CC(=O)[O-]", 5, 256, [21, 22, 23, 25, 114, 145, 146, 147, 150, 151]),
    ("ClC(Br)(I)F", 7, 1024, [148, 156, 174, 192, 227, 245, 263, 533, 536, 544, 562, 580, 673, 765, 783]),
])
def test_path_fingerprint_bits_are_unchanged(smiles, path_length, fp_size, bits):
    assert _smsd.path_fingerprint(smsd.parse_smiles(smiles), path_length, fp_size) == bits
    assert smsd.path_fingerprint(smsd.parse_smiles(smiles), path_length, fp_size) == bits
