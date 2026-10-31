# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
"""A resident integer provider preserves corpus substructure search results."""
import ctypes
import os
import subprocess
import sys
import threading

import pytest
import smsd
from smsd import _accelerator, _smsd


_DEVICE_INFO = ctypes.CFUNCTYPE(ctypes.c_int, ctypes.c_void_p, ctypes.c_size_t, ctypes.c_void_p, ctypes.c_size_t)
_LABEL_OVERLAPS = ctypes.CFUNCTYPE(ctypes.c_int, ctypes.POINTER(ctypes.c_int32), ctypes.c_int32,
                                   ctypes.c_int64, ctypes.POINTER(ctypes.c_uint64), ctypes.c_void_p,
                                   ctypes.POINTER(ctypes.c_int32), ctypes.c_void_p, ctypes.c_size_t)


def _write(buffer, size, text):
    if size:
        data = text.encode()[:size - 1] + b"\0"
        ctypes.memmove(buffer, data, len(data))


class _Hit(ctypes.Structure):
    _fields_ = [("target_index", ctypes.c_int32), ("common_bits", ctypes.c_int32)]


_MINIMUM = ctypes.CFUNCTYPE(ctypes.c_int, ctypes.c_uint32, ctypes.POINTER(ctypes.c_uint64),
                           ctypes.c_void_p, ctypes.c_size_t)
_CREATE = ctypes.CFUNCTYPE(ctypes.c_int, ctypes.c_int64, ctypes.c_int32,
                          ctypes.POINTER(ctypes.c_uint64), ctypes.POINTER(ctypes.c_int32),
                          ctypes.POINTER(ctypes.c_void_p), ctypes.c_void_p, ctypes.c_size_t)
_SCREEN = ctypes.CFUNCTYPE(ctypes.c_int, ctypes.c_void_p, ctypes.c_int32,
                          ctypes.POINTER(ctypes.c_uint64), ctypes.POINTER(ctypes.c_int32), ctypes.c_int32,
                          ctypes.POINTER(_Hit), ctypes.c_uint64, ctypes.POINTER(ctypes.c_uint64),
                          ctypes.c_void_p, ctypes.c_size_t)
_DESTROY = ctypes.CFUNCTYPE(ctypes.c_int, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_size_t)


class _ResidentTable(ctypes.Structure):
    _fields_ = [("abi_version", ctypes.c_uint32), ("struct_size", ctypes.c_uint32),
                ("smsd_version", ctypes.c_char_p), ("name", ctypes.c_char_p),
                ("device_info", _DEVICE_INFO), ("label_overlaps", _LABEL_OVERLAPS),
                ("minimum_targets", _MINIMUM), ("create_fingerprints", _CREATE),
                ("fingerprint_screen", _SCREEN), ("create_labels", ctypes.c_void_p),
                ("resident_label_overlaps", ctypes.c_void_p), ("destroy_resident", _DESTROY)]


class Provider:
    def __init__(self):
        self.uploads = self.calls = 0
        self.failure = False
        self.minimum = 0
        self.before_screen = None
        self.records = {}
        self.corpora = []

        def device_info(text, size, error, error_size):
            _write(text, size, "Test device")
            return 0

        def legacy(*args):
            return 1

        def minimum(operation, output, error, size):
            output[0] = self.minimum if operation == 3 else (1 << 64) - 1
            return 0

        def create(targets, words, bits, counts, output, error, size):
            self.uploads += 1
            record = (targets, words, list(bits[:targets * words]), list(counts[:targets]))
            self.records[self.uploads] = record
            output[0] = self.uploads
            return 0

        def screen(handle, queries, words, cutoffs, stride, hits, capacity, offsets, error, size):
            self.calls += 1
            if self.before_screen is not None:
                self.before_screen()
            if self.failure:
                _write(error, size, "device lost")
                return 1
            targets, target_words, bits, counts = self.records[handle]
            selected = 0
            offsets[0] = 0
            for query in range(queries):
                for target in range(targets):
                    common = sum((words[query * target_words + word] & bits[word * targets + target]).bit_count()
                                 for word in range(target_words))
                    if common >= cutoffs[query * stride + counts[target]]:
                        assert selected < capacity
                        hits[selected] = _Hit(target, common)
                        selected += 1
                offsets[query + 1] = selected
            return 0

        def destroy(handle, error, size):
            del self.records[handle]
            return 0

        self.callbacks = (_DEVICE_INFO(device_info), _LABEL_OVERLAPS(legacy), _MINIMUM(minimum),
                          _CREATE(create), _SCREEN(screen), _DESTROY(destroy))
        info, overlaps, minimum_callback, create_callback, screen_callback, destroy_callback = self.callbacks
        self.version = smsd.__version__.encode()
        self.table = _ResidentTable(1, ctypes.sizeof(_ResidentTable), self.version, b"Test", info, overlaps,
                                    minimum_callback, create_callback, screen_callback, None, None, destroy_callback)

    def attach(self):
        new = ctypes.pythonapi.PyCapsule_New
        new.restype = ctypes.py_object
        new.argtypes = [ctypes.c_void_p, ctypes.c_char_p, ctypes.c_void_p]
        _smsd._attach_screen_accelerator(new(ctypes.addressof(self.table), b"smsd.screen_accelerator", None))

    def corpus(self, targets):
        corpus = smsd.TargetCorpus.from_smiles(targets)
        self.corpora.append(corpus)
        return corpus


@pytest.fixture
def provider(monkeypatch):
    monkeypatch.delenv("SMSD_DISABLE_GPU", raising=False)
    monkeypatch.setitem(sys.modules, "smsd_cuda", None)
    _accelerator._reset()
    value = Provider()
    value.attach()
    yield value
    # The provider owns its callbacks until every corpus has released its handle.
    for corpus in value.corpora:
        corpus.prewarm(num_threads=1)
    _accelerator._reset()
    assert not value.records


def result_fields(result):
    return result.matches, result.status, result.optimal, result.complete, result.score


@pytest.mark.parametrize("query", ["C", "CCO", "C.O", "C1CC1"])
def test_subset_screen_preserves_original_result_slots(provider, query):
    targets = ["CCO", "CCN", "CCCO", "COC", "C1CC1", "C.O", "[H]"]
    corpus = provider.corpus(targets)
    prepared = smsd.Substructure.find(query)
    expected = prepared.match_many(targets)
    actual = prepared.match_many(corpus, num_threads=2)
    assert list(map(result_fields, actual)) == list(map(result_fields, expected))
    assert provider.calls == 1
    assert provider.uploads == 1


def test_repeated_queries_reuse_and_replace_target_paths(provider):
    prepared = smsd.Substructure.find("C")
    corpus = provider.corpus(["C", "O", "CC"])
    assert [result.found for result in prepared.match_many(corpus)] == [True, False, True]
    assert [result.found for result in prepared.match_many(corpus)] == [True, False, True]
    assert provider.uploads == 1
    corpus.add_target("C")
    assert [result.found for result in prepared.match_many(corpus)] == [True, False, True, True]
    assert provider.uploads == 2


def test_unsupported_queries_and_unmeasured_sizes_keep_cpu_route(provider):
    pseudo = smsd.MolGraphBuilder().atom_count(1).atomic_numbers([0]).neighbors([[]]).build(False)
    targets = ["C", "O", pseudo]
    corpus = provider.corpus(targets)
    options = smsd.Substructure.Options()
    options.chemistry.match_atom_type = False
    for query in (smsd.Substructure.smarts("[#6,#8]"), smsd.Substructure.find(pseudo),
                  smsd.Substructure.find("[H]"), smsd.Substructure.find("C").with_options(options)):
        assert list(map(result_fields, query.match_many(corpus))) == list(map(result_fields, query.match_many(targets)))
    assert provider.calls == 0
    assert [result.found for result in smsd.Substructure.find("C").match_many(corpus)] == [True, False, False]
    assert provider.calls == 1
    provider.minimum = (1 << 64) - 1
    assert [result.found for result in smsd.Substructure.find("C").match_many(corpus)] == [True, False, False]
    assert provider.calls == 1


def test_device_failure_counts_once_and_keeps_cpu_results(provider):
    corpus = provider.corpus(["CCO", "NN", "CC"])
    query = smsd.Substructure.find("C")
    failures = _smsd._screen_accelerator_failures()[0]
    provider.failure = True
    assert [result.found for result in query.match_many(corpus)] == [True, False, True]
    assert _smsd._screen_accelerator_failures()[0] == failures + 1


def test_immediate_timeout_bypasses_prefilter(provider):
    corpus = provider.corpus(["CCO", "NN"])
    actual = smsd.Substructure.find("O").with_timeout(0).match_many(corpus)
    assert all(result.status == smsd.SearchStatus.TIMED_OUT and not result.optimal and not result.complete for result in actual)
    assert provider.calls == 0


def test_finite_node_limit_preserves_cpu_outcomes(provider):
    targets = ["C1CC1.C", "CCCC", "CCCCC", "NNNN"]
    corpus = provider.corpus(targets)
    options = smsd.Substructure.Options()
    options.node_limit = 1
    query = smsd.Substructure.find("CCCC").with_options(options)
    expected = list(map(result_fields, query.match_many(targets, num_threads=1)))
    assert list(map(result_fields, query.match_many(corpus, num_threads=1))) == expected
    assert provider.calls == 0
    assert provider.uploads == 0


def test_cancelled_query_bypasses_prefilter(provider):
    corpus = provider.corpus(["CCO", "NN"])
    options = smsd.Substructure.Options()
    options.cancellation = smsd.CancellationToken()
    options.cancellation.cancel()
    actual = smsd.Substructure.find("O").with_options(options).match_many(corpus)
    assert all(result.status == smsd.SearchStatus.CANCELLED and not result.optimal and not result.complete for result in actual)
    assert provider.calls == 0


def _concurrent_append():
    os.environ.pop("SMSD_DISABLE_GPU", None)
    sys.modules["smsd_cuda"] = None
    _accelerator._reset()
    provider = Provider()
    provider.attach()
    corpus = provider.corpus(["CCO", "NN", "CC"])
    appended = smsd.parse_smiles("O")
    entered = threading.Event()
    attempting = threading.Event()
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
    results = smsd.Substructure.find("C").match_many(corpus)
    writer.join(30)
    assert not writer.is_alive(), "append thread did not finish"
    assert not errors, errors
    assert [result.found for result in results] == [True, False, True]
    assert len(corpus) == 4
    corpus.prewarm(num_threads=1)
    _accelerator._reset()
    assert not provider.records


def test_concurrent_append_keeps_captured_result_slots():
    code = """import importlib.util, sys
spec = importlib.util.spec_from_file_location('corpus_screen_case', sys.argv[1])
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
module._concurrent_append()
"""
    completed = subprocess.run([sys.executable, "-c", code, __file__], capture_output=True, text=True, timeout=60)
    assert completed.returncode == 0, completed.stdout + completed.stderr
