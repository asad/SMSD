# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
"""Conservative corpus filters retain exact witnesses and explicit threshold proofs."""
import ctypes
from dataclasses import FrozenInstanceError
import os
import subprocess
import sys
import threading

import pytest
import smsd
from smsd import _accelerator, _smsd

from .test_corpus_substructure_screen import Provider, _write


class _Label(ctypes.Structure):
    _fields_ = [("label", ctypes.c_uint16), ("count", ctypes.c_uint16)]


_CREATE_LABELS = ctypes.CFUNCTYPE(
    ctypes.c_int, ctypes.c_int64, ctypes.POINTER(ctypes.c_uint64), ctypes.POINTER(_Label),
    ctypes.POINTER(ctypes.c_void_p), ctypes.c_void_p, ctypes.c_size_t)
_LABEL_SCREEN = ctypes.CFUNCTYPE(
    ctypes.c_int, ctypes.c_void_p, ctypes.c_int32, ctypes.POINTER(ctypes.c_int32), ctypes.c_int32,
    ctypes.POINTER(ctypes.c_int32), ctypes.c_void_p, ctypes.c_size_t)


class CorpusProvider(Provider):
    def __init__(self):
        super().__init__()
        self.minimum = (1 << 64) - 1
        self.label_calls = 0

        def create(targets, offsets, entries, output, error, size):
            self.uploads += 1
            count = offsets[targets]
            self.records[self.uploads] = (
                "labels", targets, list(offsets[:targets + 1]),
                [(entries[index].label, entries[index].count) for index in range(count)])
            output[0] = self.uploads
            return 0

        def screen(handle, queries, counts, stride, output, error, size):
            self.label_calls += 1
            if self.before_screen is not None:
                self.before_screen()
            if self.failure:
                _write(error, size, "device lost")
                return 1
            kind, targets, offsets, entries = self.records[handle]
            assert kind == "labels"
            for query in range(queries):
                for target in range(targets):
                    output[query * targets + target] = sum(
                        min(counts[query * stride + label], count)
                        for label, count in entries[offsets[target]:offsets[target + 1]])
            return 0

        create_callback, screen_callback = _CREATE_LABELS(create), _LABEL_SCREEN(screen)
        self.callbacks += (create_callback, screen_callback)
        self.table.create_labels = ctypes.cast(create_callback, ctypes.c_void_p).value
        self.table.resident_label_overlaps = ctypes.cast(screen_callback, ctypes.c_void_p).value


@pytest.fixture
def provider(monkeypatch):
    monkeypatch.delenv("SMSD_DISABLE_GPU", raising=False)
    monkeypatch.setitem(sys.modules, "smsd_cuda", None)
    _accelerator._reset()
    value = CorpusProvider()
    value.attach()
    yield value
    for corpus in value.corpora:
        corpus.add_targets([])
    _accelerator._reset()
    assert not value.records


def search_fields(result):
    return result.matches, result.status, result.optimal, result.complete, result.score


def corpus_fields(outcome):
    return (outcome.target_index, outcome.atom_upper_bound, outcome.proven_below_threshold,
            outcome.threshold_reached, None if outcome.result is None else search_fields(outcome.result))


def compare_reference(outcomes, reference):
    assert [outcome.target_index for outcome in outcomes] == list(range(len(reference)))
    for outcome, exact in zip(outcomes, reference):
        if outcome.result is not None:
            assert search_fields(outcome.result) == search_fields(exact)
        if exact.status != smsd.SearchStatus.COMPLETED:
            assert not outcome.proven_below_threshold


@pytest.mark.parametrize("query", ["C", "CC", "CO", "C1CC1"])
def test_substructure_corpus_keeps_order_and_exact_witnesses(provider, query):
    targets = ["CCO", "NN", "CC", "C1CC1", "COC", "CCO"]
    corpus = provider.corpus(targets)
    prepared = smsd.Substructure.find(query)
    reference = prepared.match_many(targets, num_threads=1)
    cpu = prepared.match_corpus(corpus, num_threads=1)
    assert provider.calls == 0
    gpu = prepared.match_corpus(corpus, gpu_prefilter=True, num_threads=2)
    assert provider.calls > 0
    assert list(map(corpus_fields, cpu)) == list(map(corpus_fields, gpu))
    compare_reference(gpu, reference)
    for outcome, exact in zip(gpu, reference):
        assert outcome.atom_upper_bound is None
        assert outcome.threshold_reached == exact.found
        assert outcome.proven_below_threshold == (not exact.found)
        if outcome.threshold_reached:
            assert outcome.result is not None


@pytest.mark.parametrize("query,target", [("C", "CC"), ("CC", "CCC")])
def test_added_context_keeps_matches_despite_ecfp_non_subset(provider, query, target):
    query_fp = smsd.fingerprint(query, kind="ecfp")
    target_fp = smsd.fingerprint(target, kind="ecfp")
    assert not smsd.fingerprint_subset(query_fp, target_fp)
    outcome = smsd.Substructure.find(query).match_corpus(provider.corpus([target]), gpu_prefilter=True)[0]
    assert outcome.threshold_reached and not outcome.proven_below_threshold
    assert outcome.result.found


def test_mcs_atom_threshold_is_inclusive_and_does_not_invent_a_small_optimum(provider):
    targets = ["CO", "CCO", "N", "N"]
    corpus = provider.corpus(targets)
    prepared = smsd.MCS.find("CO")
    reference = prepared.match_many(targets, num_threads=1)
    cpu = prepared.match_corpus(corpus, 2, num_threads=1)
    assert provider.label_calls == 0
    gpu = prepared.match_corpus(corpus, 2, gpu_prefilter=True, num_threads=2)
    assert provider.label_calls > 0
    assert list(map(corpus_fields, cpu)) == list(map(corpus_fields, gpu))
    compare_reference(gpu, reference)
    for outcome in gpu[:2]:
        assert outcome.atom_upper_bound == 2 and outcome.threshold_reached
        assert outcome.result is not None and outcome.result.best.atom_count == 2
    for outcome in gpu[2:]:
        assert outcome.atom_upper_bound < 2 and outcome.proven_below_threshold
        assert not outcome.threshold_reached and outcome.result is None
    assert all(result.optimal for result in reference)
    with pytest.raises((FrozenInstanceError, AttributeError)):
        gpu[0].target_index = 9


def test_a_loose_atom_bound_still_requires_exact_graph_search(provider):
    prepared = smsd.MCS.find("C1CC1")
    outcome = prepared.match_corpus(provider.corpus(["CCC"]), 3, gpu_prefilter=True)[0]
    assert outcome.atom_upper_bound == 3
    assert outcome.result is not None and outcome.result.optimal
    assert outcome.result.best.atom_count == 2
    assert outcome.proven_below_threshold and not outcome.threshold_reached


def test_isotope_wildcard_and_relaxed_bonds_remain_candidates(provider):
    options = smsd.MCS.Options()
    options.chemistry.match_isotope = True
    query = smsd.MCS.find("[13C]").with_options(options)
    outcome = query.match_corpus(provider.corpus(["C"]), 1, gpu_prefilter=True)[0]
    assert outcome.threshold_reached and outcome.result.best.atom_count == 1
    options = smsd.Substructure.Options()
    options.chemistry.match_bond_order = smsd.BondOrderMode.ANY
    query = smsd.Substructure.find("C=C").with_options(options)
    outcome = query.match_corpus(provider.corpus(["CC"]), gpu_prefilter=True)[0]
    assert outcome.threshold_reached and outcome.result.found


def test_predicates_and_untyped_or_pseudo_queries_keep_exact_cpu_results(provider):
    pseudo = smsd.MolGraphBuilder().atom_count(1).atomic_numbers([0]).neighbors([[]]).build(False)
    targets = ["C", "O", pseudo]
    corpus = provider.corpus(targets)
    untyped = smsd.Substructure.Options()
    untyped.chemistry.match_atom_type = False
    queries = [smsd.Substructure.smarts("[#6,#8]"), smsd.Substructure.smarts("[$([#6]),#8]"),
               smsd.Substructure.find(pseudo), smsd.Substructure.find("[H]"),
               smsd.Substructure.find("C").with_options(untyped)]
    for query in queries:
        outcomes = query.match_corpus(corpus, gpu_prefilter=True)
        reference = query.match_many(targets)
        compare_reference(outcomes, reference)
        assert all(outcome.result is not None for outcome in outcomes)
    assert provider.calls == provider.label_calls == 0


@pytest.mark.parametrize("setting", ["mces", "weights", "excluded", "node_limit"])
def test_unsupported_mcs_bounds_keep_actual_search_results(provider, setting):
    targets = ["CCC", "NNN"]
    options = smsd.MCS.Options()
    if setting == "mces":
        options.objective = smsd.MCS.Objective.MCES
    elif setting == "weights":
        options.atom_weights = [1.0, -5.0, 1.0]
    elif setting == "excluded":
        options.excluded_target_atoms = [0]
    else:
        options.node_limit = 1
    prepared = smsd.MCS.find("CCC").with_options(options)
    reference = prepared.match_many(targets, num_threads=1)
    outcomes = prepared.match_corpus(provider.corpus(targets), 3, gpu_prefilter=True, num_threads=1)
    compare_reference(outcomes, reference)
    assert all(outcome.result is not None and outcome.atom_upper_bound is None for outcome in outcomes)
    assert provider.label_calls == 0
    if setting in ("mces", "weights"):
        assert all(not outcome.proven_below_threshold for outcome in outcomes)


@pytest.mark.parametrize("kind", ["substructure", "mcs"])
@pytest.mark.parametrize("stop", ["timeout", "cancelled"])
def test_stopped_search_is_unknown_rather_than_a_negative(provider, kind, stop):
    prepared = smsd.Substructure.find("O") if kind == "substructure" else smsd.MCS.find("O")
    if stop == "timeout":
        prepared = prepared.with_timeout(0)
        expected_status = smsd.SearchStatus.TIMED_OUT
    else:
        options = prepared.Options()
        options.cancellation = smsd.CancellationToken()
        options.cancellation.cancel()
        prepared = prepared.with_options(options)
        expected_status = smsd.SearchStatus.CANCELLED
    corpus = provider.corpus(["NN", "N"])
    outcomes = (prepared.match_corpus(corpus, gpu_prefilter=True) if kind == "substructure"
                else prepared.match_corpus(corpus, 1, gpu_prefilter=True))
    assert provider.calls == provider.label_calls == 0
    for outcome in outcomes:
        assert outcome.result.status == expected_status
        assert not outcome.result.optimal and not outcome.result.complete
        assert not outcome.proven_below_threshold and not outcome.threshold_reached


@pytest.mark.parametrize("kind", ["substructure", "mcs"])
def test_provider_failure_keeps_the_same_cpu_threshold_decisions(provider, kind):
    corpus = provider.corpus(["CCO", "NN", "CC"])
    prepared = smsd.Substructure.find("CO") if kind == "substructure" else smsd.MCS.find("CO")
    args = () if kind == "substructure" else (2,)
    expected = prepared.match_corpus(corpus, *args)
    failures = _smsd._screen_accelerator_failures()[0]
    provider.failure = True
    actual = prepared.match_corpus(corpus, *args, gpu_prefilter=True)
    assert list(map(corpus_fields, actual)) == list(map(corpus_fields, expected))
    assert _smsd._screen_accelerator_failures()[0] == failures + 1


def test_corpus_results_keep_captured_external_atom_indices(provider):
    chemistry = pytest.importorskip("rdkit.Chem")
    target = chemistry.RenumberAtoms(chemistry.MolFromSmiles("NCCO"), [3, 1, 0, 2])
    oxygen = next(atom.GetIdx() for atom in target.GetAtoms() if atom.GetAtomicNum() == 8)
    corpus = provider.corpus([target])
    target.GetAtomWithIdx(oxygen).SetAtomicNum(16)
    for prepared, args in ((smsd.Substructure.find("O"), ()), (smsd.MCS.find("O"), (1,))):
        outcome = prepared.match_corpus(corpus, *args, gpu_prefilter=True)[0]
        assert outcome.threshold_reached
        assert dict(outcome.result.best.atom_mapping) == {0: oxygen}


def _concurrent_append(kind):
    os.environ.pop("SMSD_DISABLE_GPU", None)
    sys.modules["smsd_cuda"] = None
    _accelerator._reset()
    provider = CorpusProvider()
    provider.attach()
    corpus = provider.corpus(["CCO", "NN", "CC"])
    appended = smsd.parse_smiles("O")
    entered, attempting = threading.Event(), threading.Event()
    errors = []

    def pause_screen():
        entered.set()
        assert attempting.wait(30), "append thread did not start"

    def append():
        try:
            assert entered.wait(30), "screen callback did not start"
            attempting.set()
            corpus.add_target(appended)
        except BaseException as error:
            errors.append(error)

    provider.before_screen = pause_screen
    writer = threading.Thread(target=append, daemon=True)
    writer.start()
    prepared = smsd.Substructure.find("C") if kind == "substructure" else smsd.MCS.find("C")
    args = () if kind == "substructure" else (1,)
    results = prepared.match_corpus(corpus, *args, gpu_prefilter=True)
    writer.join(30)
    assert not writer.is_alive() and not errors
    assert [outcome.target_index for outcome in results] == [0, 1, 2]
    assert [outcome.threshold_reached for outcome in results] == [True, False, True]
    assert len(corpus) == 4
    provider.before_screen = None
    assert len(prepared.match_corpus(corpus, *args, gpu_prefilter=True)) == 4
    corpus.add_targets([])
    _accelerator._reset()
    assert not provider.records


@pytest.mark.parametrize("kind", ["substructure", "mcs"])
def test_concurrent_append_keeps_one_captured_target_generation(kind):
    code = """import importlib, pathlib, sys
sys.path.insert(0, str(pathlib.Path(sys.argv[1]).resolve().parents[2]))
module = importlib.import_module('python.tests.test_corpus_search_prefilter')
module._concurrent_append(sys.argv[2])
"""
    completed = subprocess.run([sys.executable, "-c", code, __file__, kind],
                               capture_output=True, text=True, timeout=60)
    assert completed.returncode == 0, completed.stdout + completed.stderr


@pytest.mark.parametrize("minimum", [0, -1, True, 1.5])
def test_mcs_corpus_threshold_requires_a_positive_integer(minimum):
    with pytest.raises((ValueError, TypeError)):
        smsd.MCS.find("C").match_corpus(smsd.TargetCorpus.from_smiles(["C"]), minimum)
