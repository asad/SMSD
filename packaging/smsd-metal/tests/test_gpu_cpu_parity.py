# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
# Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
# See the NOTICE file for attribution, trademark, and algorithm IP terms.
"""Screens that run on the Metal device give exactly the results of the same screens on the CPU.

smsd reaches the device through a forwarding copy of the installed table, whose entries call the
device and count the calls that succeed. Every device result below is compared with the CPU and
must have reached the device, so no test compares the CPU with itself. Automatic fingerprint
screens use an Apple M5 from 100,000 targets, and SMSD_METAL_MIN_TARGETS sets that minimum for any
device; label and substructure screens have no measured minimum, so the forwarding copy reports
zero where a test needs the automatic route to reach the device. The corpora are drawn at random,
with fixed seeds, from a combinatorial library, and the table entries are also checked directly
against random fingerprints and label histograms.
"""
import collections
import ctypes
import gc
import itertools
import math
import random

import pytest

smsd = pytest.importorskip("smsd")
smsd_metal = pytest.importorskip("smsd_metal")
from smsd import _accelerator, _smsd  # noqa: E402

pytestmark = pytest.mark.skipif(not smsd_metal.is_available(), reason=str(smsd_metal.unavailable_reason()))

# The measured Apple M5 minimum for automatic fingerprint screens.
LARGE = 100000

QUERIES = ["c1ccccc1", "c1ccc(O)cc1", "CCO", "CC(=O)O", "C1CCCCC1", "c1ccncc1", "O=C(O)c1ccccc1O",
           "CN1C=NC2=C1C(=O)N(C(=O)N2C)C", "[NH4+]", "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O", "c1ccc2ccccc2c1",
           "ClC(Cl)Cl", "N[C@@H](C)C(=O)O", "c1cc[nH]c1", "FC(F)(F)c1ccccc1", "OC1C(O)C(O)C(O)C(O)C1O"]

SUBSTRUCTURES = ["c1ccccc1", "CCO", "C(=O)O", "C(=O)N", "c1ccncc1", "C1CCCCC1", "FC(F)F",
                 "CS(=O)(=O)C", "Cl", "[N+](=O)[O-]", "C#N", "c1ccc2ccccc2c1"]

THRESHOLDS = (0.0, 0.2, 1 / 3, math.nextafter(1 / 3, 1.0), 0.5, math.nextafter(0.5, 0.0), 0.7, 1.0)

# A library molecule joins a prefix, a linker, a core and a suffix.
_PREFIXES = ("", "C", "N", "O", "F", "Cl", "Br", "CC", "CO", "CN", "C=C", "N#C", "OC(=O)", "FC(F)(F)",
             "CS(=O)(=O)", "[NH3+]C", "c1ccccc1", "C1CC1", "CC(C)", "OCC")
_LINKERS = ("", "C", "CC", "CCC", "CCCC", "O", "N", "C(=O)", "CO", "CN", "S")
_CORES = ("c1ccccc1", "c1ccncc1", "c1cncnc1", "c1ccc2ccccc2c1", "c1ccc2[nH]ccc2c1", "c1cnc2ccccc2c1",
          "c1ccoc1", "c1ccsc1", "c1cc[nH]c1", "C1CCCCC1", "C1CCNCC1", "C1CCOCC1", "C1CC1", "C1CCC1",
          "C1CCCC1", "C1=CCCC1", "C1CC2CCC1C2", "C(=O)N", "C(=O)O", "NC(=O)N", "S", "C=C", "C#C",
          "N=N", "OCCO")
_SUFFIXES = ("", "C", "O", "N", "F", "Cl", "Br", "I", "C(=O)O", "C#N", "OC", "N(C)C", "S", "C(F)(F)F",
             "[N+](=O)[O-]", "CC", "CCC", "C=O", "OCC", "c1ccccc1")


class LabelCount(ctypes.Structure):
    """The smsd_label_count layout."""

    _fields_ = [("label", ctypes.c_uint16), ("count", ctypes.c_uint16)]


class Hit(ctypes.Structure):
    """The smsd_fingerprint_hit layout."""

    _fields_ = [("target_index", ctypes.c_int32), ("common_bits", ctypes.c_int32)]


class Table(ctypes.Structure):
    """The smsd_screen_accelerator_v1 layout, with every pointer held as an address."""

    _fields_ = [("abi_version", ctypes.c_uint32), ("struct_size", ctypes.c_uint32),
                ("smsd_version", ctypes.c_void_p), ("name", ctypes.c_void_p),
                ("device_info", ctypes.c_void_p), ("label_overlaps", ctypes.c_void_p),
                ("minimum_targets", ctypes.c_void_p), ("create_fingerprints", ctypes.c_void_p),
                ("fingerprint_screen", ctypes.c_void_p), ("create_labels", ctypes.c_void_p),
                ("resident_label_overlaps", ctypes.c_void_p), ("destroy_resident", ctypes.c_void_p)]


_c = ctypes
# Every entry the forwarding copy wraps; each ends with its error buffer and that buffer's size.
_FORWARDED = {
    "label_overlaps": _c.CFUNCTYPE(_c.c_int, _c.c_void_p, _c.c_int32, _c.c_int64, _c.c_void_p, _c.c_void_p,
                                   _c.c_void_p, _c.c_void_p, _c.c_size_t),
    "minimum_targets": _c.CFUNCTYPE(_c.c_int, _c.c_uint32, _c.c_void_p, _c.c_void_p, _c.c_size_t),
    "create_fingerprints": _c.CFUNCTYPE(_c.c_int, _c.c_int64, _c.c_int32, _c.c_void_p, _c.c_void_p,
                                        _c.c_void_p, _c.c_void_p, _c.c_size_t),
    "fingerprint_screen": _c.CFUNCTYPE(_c.c_int, _c.c_void_p, _c.c_int32, _c.c_void_p, _c.c_void_p, _c.c_int32,
                                       _c.c_void_p, _c.c_uint64, _c.c_void_p, _c.c_void_p, _c.c_size_t),
    "create_labels": _c.CFUNCTYPE(_c.c_int, _c.c_int64, _c.c_void_p, _c.c_void_p, _c.c_void_p,
                                  _c.c_void_p, _c.c_size_t),
    "resident_label_overlaps": _c.CFUNCTYPE(_c.c_int, _c.c_void_p, _c.c_int32, _c.c_void_p, _c.c_int32,
                                            _c.c_void_p, _c.c_void_p, _c.c_size_t),
    "destroy_resident": _c.CFUNCTYPE(_c.c_int, _c.c_void_p, _c.c_void_p, _c.c_size_t),
}
_LabelOverlaps = _c.CFUNCTYPE(_c.c_int, _c.POINTER(_c.c_int32), _c.c_int32, _c.c_int64, _c.POINTER(_c.c_uint64),
                              _c.POINTER(LabelCount), _c.POINTER(_c.c_int32), _c.c_void_p, _c.c_size_t)
_CreateLabels = _c.CFUNCTYPE(_c.c_int, _c.c_int64, _c.POINTER(_c.c_uint64), _c.POINTER(LabelCount),
                             _c.POINTER(_c.c_void_p), _c.c_void_p, _c.c_size_t)
_ResidentLabelOverlaps = _c.CFUNCTYPE(_c.c_int, _c.c_void_p, _c.c_int32, _c.POINTER(_c.c_int32), _c.c_int32,
                                      _c.POINTER(_c.c_int32), _c.c_void_p, _c.c_size_t)
_Minimum = _c.CFUNCTYPE(_c.c_int, _c.c_uint32, _c.POINTER(_c.c_uint64), _c.c_void_p, _c.c_size_t)
_CreateFingerprints = _c.CFUNCTYPE(_c.c_int, _c.c_int64, _c.c_int32, _c.POINTER(_c.c_uint64),
                                   _c.POINTER(_c.c_int32), _c.POINTER(_c.c_void_p), _c.c_void_p, _c.c_size_t)
_FingerprintScreen = _c.CFUNCTYPE(_c.c_int, _c.c_void_p, _c.c_int32, _c.POINTER(_c.c_uint64),
                                  _c.POINTER(_c.c_int32), _c.c_int32, _c.POINTER(Hit), _c.c_uint64,
                                  _c.POINTER(_c.c_uint64), _c.c_void_p, _c.c_size_t)

# Forwarding copies stay valid until the end of the session, since a corpus releases its device
# handle through the table that created it, possibly after the test that attached it.
_KEPT = []


def _device_table():
    table = Table.from_address(smsd_metal._load().smsd_metal_screen_accelerator())
    assert table.abi_version == 1 and table.struct_size >= ctypes.sizeof(Table)
    return table


def _device_minimum(operation):
    minimum, error = ctypes.c_uint64(), ctypes.create_string_buffer(256)
    assert _Minimum(_device_table().minimum_targets)(operation, ctypes.byref(minimum), error, len(error)) == 0, \
        error.value
    return minimum.value


def _write(buffer, size, text):
    if buffer and size:
        data = text.encode()[:size - 1] + b"\0"
        ctypes.memmove(buffer, data, len(data))


class Forwarder:
    """A copy of the device table whose entries call the device and count the calls that succeed.

    ``zero_minimum`` names the operations whose automatic minimum is reported as zero, after the
    device has answered. ``resident=False`` publishes only the required prefix, so that smsd uses
    the single-query label entry point. ``fail`` names entries that report a failure instead of
    calling the device.
    """

    def __init__(self, zero_minimum=(), resident=True, fail=()):
        device = _device_table()
        self.calls = collections.Counter()
        self.zero_minimum = frozenset(zero_minimum)
        self.fail = frozenset(fail)
        self.failures = _smsd._screen_accelerator_failures()[0]
        self.callbacks = []
        self.table = Table(device.abi_version, ctypes.sizeof(Table), device.smsd_version, device.name,
                           device.device_info)
        for name, prototype in _FORWARDED.items():
            setattr(self.table, name, self._forward(name, prototype, prototype(getattr(device, name))))
        if not resident:
            self.table.struct_size = Table.minimum_targets.offset
        _KEPT.append(self)

    def _forward(self, name, prototype, target):
        def call(*args):
            try:
                if name in self.fail:
                    _write(args[-2], args[-1], f"the forwarding entry {name} was told to fail")
                    self.calls["failed " + name] += 1
                    return 1
                status = target(*args)
                if status == 0:
                    if name == "minimum_targets" and args[0] in self.zero_minimum:
                        ctypes.c_uint64.from_address(args[1]).value = 0
                    self.calls[name] += 1
                return status
            except BaseException as error:  # a callback that raised must not report success
                _write(args[-2], args[-1], f"the forwarding entry {name} failed: {error!r}")
                return 1

        callback = prototype(call)
        self.callbacks.append(callback)
        return ctypes.cast(callback, ctypes.c_void_p).value

    def attach(self):
        new = ctypes.pythonapi.PyCapsule_New
        new.restype = ctypes.py_object
        new.argtypes = [ctypes.c_void_p, ctypes.c_char_p, ctypes.c_void_p]
        _smsd._attach_screen_accelerator(new(ctypes.addressof(self.table), b"smsd.screen_accelerator", None))
        assert _smsd._screen_accelerator()[0] == "Metal"

    def check(self):
        """No accelerated screen ran on the CPU instead."""
        failures, reason = _smsd._screen_accelerator_failures()
        assert failures == self.failures, reason

    def release(self, *corpora):
        """Drop the corpora's device state; every handle this copy created must have been destroyed."""
        for corpus in corpora:
            corpus.add_targets([])
        gc.collect()
        assert self.calls["destroy_resident"] == self.calls["create_fingerprints"] + self.calls["create_labels"]
        self.check()


@pytest.fixture
def device(monkeypatch):
    """Attach forwarding copies of the installed device table for one test."""
    monkeypatch.delenv("SMSD_DISABLE_GPU", raising=False)
    monkeypatch.delenv("SMSD_METAL_MIN_TARGETS", raising=False)
    if "smsd_metal" not in _accelerator._PACKAGES:
        # An smsd that looks only for other accelerators is given the package to look for.
        monkeypatch.setattr(_accelerator, "_PACKAGES", _accelerator._PACKAGES + ("smsd_metal",))
    _accelerator._reset()
    # smsd attaches the installed package first; the forwarding copy then takes its place.
    _accelerator.discover()
    assert _smsd._screen_accelerator()[0] == "Metal", "smsd did not attach the installed device"

    def attach(**options):
        forwarder = Forwarder(**options)
        forwarder.attach()
        return forwarder

    yield attach
    _accelerator._reset()


def _on_device(forwarder, entry, function, *args, **kwargs):
    """function(*args, **kwargs), which must reach the device entry."""
    before = forwarder.calls[entry]
    result = function(*args, **kwargs)
    assert forwarder.calls[entry] > before, f"the call did not reach the device entry {entry}"
    return result


def _on_cpu(monkeypatch, forwarder, function, *args, **kwargs):
    """function(*args, **kwargs) with SMSD_DISABLE_GPU set, which must make no device call."""
    before = sum(forwarder.calls.values())
    monkeypatch.setenv("SMSD_DISABLE_GPU", "1")
    try:
        result = function(*args, **kwargs)
    finally:
        monkeypatch.delenv("SMSD_DISABLE_GPU")
    assert sum(forwarder.calls.values()) == before
    return result


def _random_graphs(seed, limit):
    """limit parsed library molecules drawn at random, with repeats, by a generator of this seed."""
    rng, graphs = random.Random(seed), []
    while len(graphs) < limit:
        text = rng.choice(_PREFIXES) + rng.choice(_LINKERS) + rng.choice(_CORES) + rng.choice(_SUFFIXES)
        try:
            graphs.append(smsd.parse_smiles(text))
        except Exception:  # noqa: BLE001 - a combination smsd rejects is left out of the corpus
            continue
    return graphs


def _corpus(graphs):
    corpus, graphs = smsd.TargetCorpus(), iter(graphs)
    for chunk in iter(lambda: list(itertools.islice(graphs, 10000)), []):
        corpus.add_targets(chunk)
    return corpus


@pytest.fixture(scope="module")
def random_graphs():
    return _random_graphs(20261031, 2000)


def result_fields(result):
    return result.matches, result.status, result.optimal, result.complete, result.score


def corpus_fields(outcome):
    return (outcome.target_index, outcome.atom_upper_bound, outcome.proven_below_threshold,
            outcome.threshold_reached, None if outcome.result is None else result_fields(outcome.result))


# ---------------------------------------------------------------------------------------------
# Table entries against references computed here
# ---------------------------------------------------------------------------------------------

def _random_fingerprint(rng, words):
    """Words of a random density, including empty and full fingerprints and the top bit."""
    density = rng.choice((0.0, 0.02, 0.1, 0.5, 0.9, 1.0))
    value = 0
    for bit in range(words * 64):
        if rng.random() < density:
            value |= 1 << bit
    return [(value >> (64 * word)) & ((1 << 64) - 1) for word in range(words)]


def _bits(fingerprint):
    return sum(word.bit_count() for word in fingerprint)


@pytest.mark.parametrize("seed,target_count,query_count,words", [
    (1, 1, 1, 1), (2, 255, 3, 1), (3, 256, 64, 2), (4, 257, 65, 3), (5, 3001, 4, 32), (6, 1000, 70, 16)])
def test_device_fingerprint_screen_equals_the_reference(seed, target_count, query_count, words):
    rng = random.Random(seed)
    table = _device_table()
    targets = [_random_fingerprint(rng, words) for _ in range(target_count)]
    queries = [_random_fingerprint(rng, words) for _ in range(query_count)]
    stride = words * 64 + 1
    # Cutoffs of every kind: zero, reachable, and above any intersection.
    cutoffs = [rng.choice((0, rng.randint(0, words * 64), rng.randint(0, 8), words * 64 + 1))
               for _ in range(query_count * stride)]
    popcounts = [_bits(target) for target in targets]
    expected = []
    for query_index, query in enumerate(queries):
        row = []
        for target_index, target in enumerate(targets):
            common = sum((q & t).bit_count() for q, t in zip(query, target))
            if common >= cutoffs[query_index * stride + popcounts[target_index]]:
                row.append((target_index, common))
        expected.append(row)

    word_major = (ctypes.c_uint64 * (words * target_count))(
        *(targets[target][word] for word in range(words) for target in range(target_count)))
    counts = (ctypes.c_int32 * target_count)(*popcounts)
    error = ctypes.create_string_buffer(256)
    handle = ctypes.c_void_p()
    assert _CreateFingerprints(table.create_fingerprints)(target_count, words, word_major, counts,
                                                          ctypes.byref(handle), error, len(error)) == 0, error.value
    try:
        screen = _FingerprintScreen(table.fingerprint_screen)
        query_words = (ctypes.c_uint64 * (query_count * words))(*(word for query in queries for word in query))
        minimum = (ctypes.c_int32 * len(cutoffs))(*cutoffs)
        total = sum(map(len, expected))
        # A buffer for every pair, and one that holds exactly the hits.
        for capacity in (query_count * target_count, total):
            hits = (Hit * max(1, capacity))()
            offsets = (ctypes.c_uint64 * (query_count + 1))()
            assert screen(handle, query_count, query_words, minimum, stride, hits, capacity, offsets,
                          error, len(error)) == 0, error.value
            assert offsets[0] == 0 and offsets[query_count] == total
            actual = [[(hits[position].target_index, hits[position].common_bits)
                       for position in range(offsets[query], offsets[query + 1])] for query in range(query_count)]
            assert actual == expected
        if total:
            hits = (Hit * total)()
            offsets = (ctypes.c_uint64 * (query_count + 1))()
            assert screen(handle, query_count, query_words, minimum, stride, hits, total - 1, offsets,
                          error, len(error)) != 0
            assert b"too short" in error.value
    finally:
        destroy = _FORWARDED["destroy_resident"](table.destroy_resident)
        assert destroy(handle, error, len(error)) == 0, error.value


def _histograms(rng, target_count, table_size):
    """Random label histograms with distinct labels, some at or above table_size; some targets have none."""
    histograms = []
    for _ in range(target_count):
        size = rng.choice((0, 1, 2, 5, 12, 40))
        labels = sorted(rng.sample(range(table_size + 3), min(size, table_size + 3)))
        histograms.append([(label, rng.choice((1, 2, 3, rng.randint(1, 60), 65535))) for label in labels])
    return histograms


def _label_queries(rng, query_count, table_size):
    """Query label counts, including zeros and counts above the largest histogram count."""
    return [[rng.choice((0, 0, 1, rng.randint(0, 50), 70000)) for _ in range(table_size)]
            for _ in range(query_count)]


@pytest.mark.parametrize("seed,target_count,query_count,table_size", [
    (1, 0, 3, 1024), (2, 1, 1, 1), (3, 257, 1, 7), (4, 256, 64, 1024), (5, 257, 65, 1024),
    (6, 1000, 130, 1024), (7, 4099, 9, 512)])
def test_device_label_overlaps_equal_the_reference(seed, target_count, query_count, table_size):
    rng = random.Random(seed)
    table = _device_table()
    histograms = _histograms(rng, target_count, table_size)
    queries = _label_queries(rng, query_count, table_size)
    # Labels at or above the table size count for nothing, as on the CPU.
    expected = [[sum(min(count, query[label]) for label, count in histogram if label < table_size)
                 for histogram in histograms] for query in queries]
    offsets = (ctypes.c_uint64 * (target_count + 1))(
        *itertools.accumulate((len(histogram) for histogram in histograms), initial=0))
    flat = [entry for histogram in histograms for entry in histogram]
    entries = (LabelCount * max(1, len(flat)))(*flat)
    error = ctypes.create_string_buffer(256)

    single = _LabelOverlaps(table.label_overlaps)
    for query, row in zip(queries[:3], expected):
        counts = (ctypes.c_int32 * table_size)(*query)
        overlaps = (ctypes.c_int32 * max(1, target_count))(*([-1] * max(1, target_count)))
        assert single(counts, table_size, target_count, offsets, entries, overlaps, error, len(error)) == 0, \
            error.value
        assert list(overlaps)[:target_count] == row

    handle = ctypes.c_void_p()
    assert _CreateLabels(table.create_labels)(target_count, offsets, entries, ctypes.byref(handle),
                                              error, len(error)) == 0, error.value
    assert handle.value
    try:
        resident = _ResidentLabelOverlaps(table.resident_label_overlaps)
        # The same upload answers a second batch in another order.
        for order in (list(range(query_count)), list(reversed(range(query_count)))):
            counts = (ctypes.c_int32 * (query_count * table_size))(
                *(count for query in order for count in queries[query]))
            size = max(1, query_count * target_count)
            overlaps = (ctypes.c_int32 * size)(*([-1] * size))
            assert resident(handle, query_count, counts, table_size, overlaps, error, len(error)) == 0, \
                error.value
            actual = [list(overlaps[position * target_count:(position + 1) * target_count])
                      for position in range(query_count)]
            assert actual == [expected[query] for query in order]
        # A fingerprint entry refuses a label handle.
        hits, hit_offsets = (Hit * 1)(), (ctypes.c_uint64 * 2)()
        assert _FingerprintScreen(table.fingerprint_screen)(handle, 0, None, None, 1, hits, 1, hit_offsets,
                                                            error, len(error)) != 0
    finally:
        destroy = _FORWARDED["destroy_resident"](table.destroy_resident)
        assert destroy(handle, error, len(error)) == 0, error.value


def test_invalid_label_tables_are_refused():
    table = _device_table()
    error = ctypes.create_string_buffer(256)
    handle = ctypes.c_void_p()
    create = _CreateLabels(table.create_labels)
    entries = (LabelCount * 2)((1, 1), (2, 1))
    assert create(2, (ctypes.c_uint64 * 3)(1, 1, 2), entries, ctypes.byref(handle), error, len(error)) != 0
    assert b"start at zero" in error.value and not handle.value
    assert create(2, (ctypes.c_uint64 * 3)(0, 2, 1), entries, ctypes.byref(handle), error, len(error)) != 0
    assert b"not ordered" in error.value and not handle.value
    assert create(1, (ctypes.c_uint64 * 2)(0, 2), None, ctypes.byref(handle), error, len(error)) != 0
    assert not handle.value


# ---------------------------------------------------------------------------------------------
# TargetCorpus screens through smsd
# ---------------------------------------------------------------------------------------------

@pytest.mark.timeout(600)
@pytest.mark.parametrize("seed,size", [(11, 1), (12, 63), (13, 257), (14, 2000), (15, 9973)])
def test_fingerprint_screen_on_random_corpora(device, monkeypatch, seed, size):
    forwarder = device()
    monkeypatch.setenv("SMSD_METAL_MIN_TARGETS", "1")
    graphs = _random_graphs(seed, size)
    corpus = _corpus(graphs)
    corpus.prewarm()
    rng = random.Random(seed)
    queries = rng.sample(graphs, min(len(graphs), 20)) + [smsd.parse_smiles(text) for text in QUERIES]
    for threshold in THRESHOLDS + (rng.random(), rng.random()):
        for query in queries[:3]:
            gpu = _on_device(forwarder, "fingerprint_screen", corpus.screen, query, threshold)
            assert gpu == _on_cpu(monkeypatch, forwarder, corpus.screen, query, threshold)
        gpu = _on_device(forwarder, "fingerprint_screen", corpus.screen_many, queries, threshold)
        assert gpu == _on_cpu(monkeypatch, forwarder, corpus.screen_many, queries, threshold)
    # The targets were uploaded once and reused by every screen.
    assert forwarder.calls["create_fingerprints"] == 1
    forwarder.release(corpus)


@pytest.fixture(scope="module")
def large_corpus():
    """100,000 prewarmed random library molecules."""
    corpus = _corpus(_random_graphs(100000, LARGE))
    assert len(corpus) == LARGE
    corpus.prewarm()
    yield corpus
    corpus.add_targets([])


@pytest.mark.timeout(900)
def test_automatic_fingerprint_screen_on_the_device_at_100000_targets(device, large_corpus, monkeypatch):
    forwarder = device()
    if smsd_metal.device_info().startswith("Apple M5 (Metal,"):
        assert _device_minimum(1) == LARGE  # the measured M5 minimum, without any setting
    else:
        monkeypatch.setenv("SMSD_METAL_MIN_TARGETS", str(LARGE))
    assert _device_minimum(1) == LARGE
    queries = _random_graphs(7, 70) + [smsd.parse_smiles(text) for text in QUERIES]
    for query in queries[:4] + queries[-4:]:
        for threshold in (0.0, 0.2, 1 / 3):
            gpu = _on_device(forwarder, "fingerprint_screen", large_corpus.screen, query, threshold)
            assert gpu == _on_cpu(monkeypatch, forwarder, large_corpus.screen, query, threshold)
    for threshold in (1 / 3, 0.5, math.nextafter(0.5, 1.0), 0.7, 1.0):
        gpu = _on_device(forwarder, "fingerprint_screen", large_corpus.screen_many, queries, threshold)
        assert gpu == _on_cpu(monkeypatch, forwarder, large_corpus.screen_many, queries, threshold)
    assert forwarder.calls["create_fingerprints"] == 1
    forwarder.check()


@pytest.mark.timeout(900)
def test_fingerprint_screen_below_the_device_minimum_stays_on_the_cpu(device, large_corpus, monkeypatch):
    forwarder = device()
    monkeypatch.setenv("SMSD_METAL_MIN_TARGETS", str(LARGE + 1))
    queries = [smsd.parse_smiles(text) for text in QUERIES]
    for threshold in (1 / 3, 1.0):
        attached = large_corpus.screen_many(queries, threshold)
        assert attached == _on_cpu(monkeypatch, forwarder, large_corpus.screen_many, queries, threshold)
    assert forwarder.calls["minimum_targets"] > 0
    assert forwarder.calls["create_fingerprints"] == forwarder.calls["fingerprint_screen"] == 0
    forwarder.check()


def test_failed_device_screen_gives_the_cpu_hits(device, random_graphs, monkeypatch):
    forwarder = device(fail=("fingerprint_screen",))
    monkeypatch.setenv("SMSD_METAL_MIN_TARGETS", "1")
    corpus = _corpus(random_graphs)
    corpus.prewarm()
    queries = random_graphs[::101]
    for threshold in (0.2, 0.5):
        failures = _smsd._screen_accelerator_failures()[0]
        attached = corpus.screen_many(queries, threshold)
        assert _smsd._screen_accelerator_failures()[0] == failures + 1
        assert attached == _on_cpu(monkeypatch, forwarder, corpus.screen_many, queries, threshold)
    assert forwarder.calls["failed fingerprint_screen"] == 2
    corpus.add_targets([])
    gc.collect()
    assert forwarder.calls["destroy_resident"] == forwarder.calls["create_fingerprints"]


@pytest.mark.timeout(300)
@pytest.mark.parametrize("route", ["resident", "single"])
@pytest.mark.parametrize("strict", [False, True])
def test_similarity_bound_hits_on_random_corpora(device, random_graphs, monkeypatch, route, strict):
    forwarder = device(zero_minimum=(2,)) if route == "resident" else device(resident=False)
    entry = "resident_label_overlaps" if route == "resident" else "label_overlaps"
    chemistry = smsd.ChemOptions()
    if strict:
        chemistry.aromaticity_mode = smsd.AromaticityMode.STRICT
        chemistry.ring_matches_ring_only = True
    graphs = [smsd.parse_smiles(text) for text in QUERIES]
    corpus = _corpus(graphs * 16 + random_graphs)
    for query in graphs:
        thresholds = {0.0, 1.0} | {smsd.similarity_upper_bound(query, target, chemistry) for target in graphs}
        for threshold in sorted(thresholds):
            gpu = _on_device(forwarder, entry, smsd.screen_targets, query, corpus, threshold, chemistry)
            assert gpu == _on_cpu(monkeypatch, forwarder, smsd.screen_targets, query, corpus, threshold, chemistry)
    queries = (graphs + random_graphs[::29])[:70]
    for threshold in (0.0, 0.5, 1.0):
        gpu = _on_device(forwarder, entry, corpus.screen_targets_many, queries, threshold, chemistry)
        assert gpu == _on_cpu(monkeypatch, forwarder, corpus.screen_targets_many, queries, threshold, chemistry)
    if route == "resident":
        assert forwarder.calls["create_labels"] == 1
    forwarder.release(corpus)


def test_similarity_bound_has_no_automatic_minimum(device, random_graphs, monkeypatch):
    forwarder = device()
    corpus = _corpus(random_graphs)
    query = smsd.parse_smiles("c1ccc(O)cc1")
    attached = smsd.screen_targets(query, corpus, 0.5)
    assert attached == _on_cpu(monkeypatch, forwarder, smsd.screen_targets, query, corpus, 0.5)
    assert forwarder.calls["create_labels"] == forwarder.calls["label_overlaps"] == 0
    forwarder.release(corpus)


@pytest.mark.timeout(300)
@pytest.mark.parametrize("query", ["c1ccccc1O", "CC(=O)N", "C1CCNCC1"])
def test_mcs_atom_bounds_from_device_label_overlaps(device, random_graphs, monkeypatch, query):
    forwarder = device()
    corpus = _corpus(random_graphs[:120])
    prepared = smsd.MCS.find(query)
    for minimum_atoms in (3, 7):
        gpu = _on_device(forwarder, "resident_label_overlaps", prepared.match_corpus, corpus, minimum_atoms,
                         gpu_prefilter=True, num_threads=1)
        cpu = _on_cpu(monkeypatch, forwarder, prepared.match_corpus, corpus, minimum_atoms, num_threads=1)
        assert list(map(corpus_fields, gpu)) == list(map(corpus_fields, cpu))
    forwarder.release(corpus)


@pytest.mark.timeout(300)
def test_automatic_substructure_prescreen_on_the_device(device, random_graphs, monkeypatch):
    forwarder = device(zero_minimum=(3,))
    corpus = _corpus(random_graphs)
    for query in SUBSTRUCTURES:
        prepared = smsd.Substructure.find(query)
        expected = list(map(result_fields, prepared.match_many(random_graphs, num_threads=1)))
        gpu = _on_device(forwarder, "fingerprint_screen", prepared.match_many, corpus, num_threads=1)
        assert list(map(result_fields, gpu)) == expected
        cpu = _on_cpu(monkeypatch, forwarder, prepared.match_many, corpus, num_threads=1)
        assert list(map(result_fields, cpu)) == expected
    # The element-path fingerprints of the targets were uploaded once for every query.
    assert forwarder.calls["create_fingerprints"] == 1
    forwarder.release(corpus)


@pytest.mark.timeout(300)
def test_requested_substructure_prescreen_on_the_device(device, random_graphs, monkeypatch):
    forwarder = device()
    corpus = _corpus(random_graphs)
    for query in SUBSTRUCTURES:
        prepared = smsd.Substructure.find(query)
        before = forwarder.calls["fingerprint_screen"]
        automatic = prepared.match_many(corpus, num_threads=1)
        # The device has no measured minimum for the automatic substructure pre-screen.
        assert forwarder.calls["fingerprint_screen"] == before
        gpu = _on_device(forwarder, "fingerprint_screen", prepared.match_corpus, corpus,
                         gpu_prefilter=True, num_threads=1)
        cpu = _on_cpu(monkeypatch, forwarder, prepared.match_corpus, corpus, num_threads=1)
        assert list(map(corpus_fields, gpu)) == list(map(corpus_fields, cpu))
        for outcome, result in zip(gpu, automatic):
            assert outcome.threshold_reached == result.found
            if outcome.result is not None:
                assert result_fields(outcome.result) == result_fields(result)
    forwarder.release(corpus)
