# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
"""Public search contracts are independent of solver implementation names."""
from concurrent.futures import ThreadPoolExecutor

import pytest
import smsd


@pytest.mark.parametrize("query,target", [("C1CC1", "CCC"), ("CCC", "C1CC1")])
@pytest.mark.parametrize("objective,atoms,bonds", [(smsd.MCS.Objective.MCIS, 2, 1),
                                                  (smsd.MCS.Objective.MCES, 3, 2)])
def test_induced_and_selected_edge_objectives(query, target, objective, atoms, bonds):
    prepared = smsd.MCS.find(query).with_objective(objective)
    result = prepared.match_all(target, limit=100)
    assert result.complete and result.optimal and result.status == smsd.SearchStatus.COMPLETED
    assert len(result.matches) == (12 if objective == smsd.MCS.Objective.MCIS else 6)
    options = smsd.MCS.Options()
    options.objective = objective
    for match in result.matches:
        assert (match.atom_count, match.bond_count) == (atoms, bonds)
        assert smsd.MCS.validate_mapping(query, target, match, options) == []
        output = smsd.parse_smiles(match.to_smiles(query))
        assert (output.n, len(output.bonds())) == (atoms, bonds)


@pytest.mark.parametrize("objective", [smsd.MCS.Objective.MCIS, smsd.MCS.Objective.MCES])
def test_disconnected_fragment_mapping(objective):
    search = smsd.MCS.find("CO.CN").with_objective(objective)
    connected = search.match("CO.CN")
    fragments = search.with_fragments(min_atoms=2, max_fragments=2).match("CO.CN")
    assert connected.best.atom_count == 2
    assert fragments.best.atom_count == 4 and fragments.best.bond_count == 2
    assert len(fragments.best.fragments) == 2
    assert sorted(len(fragment) for fragment in fragments.best.fragments) == [2, 2]
    assert smsd.parse_smiles(fragments.best.to_smiles("CO.CN")).n == 4


def test_mces_does_not_add_isolated_atoms():
    result = smsd.MCS.find("CC.O").with_objective(smsd.MCS.Objective.MCES).with_fragments().match("CC.O")
    assert (result.best.atom_count, result.best.bond_count) == (2, 1)
    assert not smsd.MCS.find("O").with_objective(smsd.MCS.Objective.MCES).match("O").found


@pytest.mark.parametrize("query", [smsd.MCS.find("CC"), smsd.Substructure.find("CC"),
                                    smsd.Substructure.smarts("[#6]-[#6]")])
def test_enumeration_cap_keeps_proof_separate_from_completion(query):
    result = query.match_all("CCC", limit=2)
    assert len(result.matches) == 2
    assert result.optimal and not result.complete
    assert result.status == smsd.SearchStatus.LIMIT_REACHED


@pytest.mark.parametrize("kind", ["mcs", "molecule", "smarts"])
def test_cancellation_and_zero_deadline(kind):
    query = {"mcs": smsd.MCS.find, "molecule": smsd.Substructure.find,
             "smarts": smsd.Substructure.smarts}[kind]("CC")
    options = query.Options()
    token = smsd.CancellationToken()
    options.cancellation = token
    prepared = query.with_options(options)
    token.cancel()
    stopped = prepared.match("CCC")
    assert stopped.status == smsd.SearchStatus.CANCELLED
    assert not stopped.optimal and not stopped.complete
    expired = query.with_timeout(0).match("CCC")
    assert expired.status == smsd.SearchStatus.TIMED_OUT
    assert not expired.optimal and not expired.complete


def test_configuration_is_copied_for_concurrent_reuse():
    options = smsd.MCS.Options()
    options.atom_weights = [10, -30, 1]
    prepared = smsd.MCS.find("CCC").with_options(options)
    options.atom_weights = [1, 100, 1]
    targets = ["CCC", "CCN", "C1CC1"] * 3
    with ThreadPoolExecutor(max_workers=3) as executor:
        results = list(executor.map(prepared.match, targets))
    assert all(set(result.best.atom_mapping) == {0} and result.score == 10 for result in results)
    assert all(result.optimal for result in results)


def test_mces_rejects_vertex_weights():
    with pytest.raises(ValueError):
        smsd.MCS.find("CC").with_objective(smsd.MCS.Objective.MCES).with_atom_weights([1, 2]).match("CC")


def test_excluded_atoms_use_external_input_numbering():
    chemistry = pytest.importorskip("rdkit.Chem")
    target = chemistry.RenumberAtoms(chemistry.MolFromSmiles("NCCO"), [3, 1, 0, 2])
    oxygen = next(atom.GetIdx() for atom in target.GetAtoms() if atom.GetAtomicNum() == 8)
    options = smsd.MCS.Options()
    options.excluded_target_atoms = [oxygen]
    search = smsd.MCS.find("O").with_options(options)
    assert not search.match(target).found
    assert not search.match_all(target).found
    assert not search.match_many([target], num_threads=2)[0].found


def test_corpus_captures_input_and_correspondence():
    chemistry = pytest.importorskip("rdkit.Chem")
    target = chemistry.RenumberAtoms(chemistry.MolFromSmiles("NCCO"), [3, 1, 0, 2])
    oxygen = next(atom.GetIdx() for atom in target.GetAtoms() if atom.GetAtomicNum() == 8)
    corpus = smsd.TargetCorpus.from_smiles([target])
    target.GetAtomWithIdx(oxygen).SetAtomicNum(16)
    smsd.clear_cache()
    result = smsd.Substructure.find("O").match_many(corpus)[0]
    assert dict(result.best.atom_mapping) == {0: oxygen}


def test_consensus_verifies_support_and_caps_output():
    result = smsd.MCS.consensus(["N", "CCO", "CCN"], threshold=2 / 3)
    assert result.found and result.graph.n == 2 and result.support_count == 2
    assert set(result.mappings) == {1, 2}
    assert not result.optimal
    output = smsd.parse_smiles(result.to_smiles())
    assert output.n == 2 and sum(output.hydrogen_count) == 6


def test_selected_graph_retains_context_while_smiles_caps_the_cut():
    match = smsd.MCS.find("CCO").match("CCN").best
    context = match.to_graph("CCO")
    output = smsd.parse_smiles(match.to_smiles("CCO"))
    assert context.n == output.n == 2
    assert sum(context.hydrogen_count) == 5
    assert sum(output.hydrogen_count) == 6
    assert smsd.Substructure.smarts("[CH2]").match(context).found


def test_retired_search_exports_are_absent():
    for module in (smsd, smsd._smsd):
        for name in ("find_mcs", "find_mcs_clique", "find_mcs_coverage", "MCSOptions", "MatcherEngine"):
            assert not hasattr(module, name)
