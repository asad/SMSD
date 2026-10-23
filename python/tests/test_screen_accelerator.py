# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
# Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
# See the NOTICE file for attribution, trademark, and algorithm IP terms.
"""Without smsd-cuda, smsd screens and reports its backend as before; an accelerator changes no hit.

A table of ctypes callbacks that computes the label overlaps on the CPU stands in for the device.
"""
import ctypes
import pathlib
import re
import sys
import types

import pytest
import smsd
from smsd import _accelerator, _smsd

SMILES = ["c1ccccc1", "c1ccc(O)cc1", "CCO", "CC(=O)O", "C1CCCCC1", "c1ccncc1", "O=C(O)c1ccccc1O",
          "CN1C=NC2=C1C(=O)N(C(=O)N2C)C", "[NH4+]", "[O-]C(=O)C", "ClC(Cl)Cl", "c1ccc2ccccc2c1", "C#N",
          "N[C@@H](C)C(=O)O", "c1cc[nH]c1", "FC(F)(F)c1ccccc1", "OC1C(O)C(O)C(O)C(O)C1O", "S"]


class _LabelCount(ctypes.Structure):
    _fields_ = [("label", ctypes.c_uint16), ("count", ctypes.c_uint16)]


_DEVICE_INFO = ctypes.CFUNCTYPE(ctypes.c_int, ctypes.c_void_p, ctypes.c_size_t, ctypes.c_void_p, ctypes.c_size_t)
_LABEL_OVERLAPS = ctypes.CFUNCTYPE(ctypes.c_int, ctypes.POINTER(ctypes.c_int32), ctypes.c_int32, ctypes.c_int64,
                                   ctypes.POINTER(ctypes.c_uint64), ctypes.POINTER(_LabelCount),
                                   ctypes.POINTER(ctypes.c_int32), ctypes.c_void_p, ctypes.c_size_t)


class _Table(ctypes.Structure):
    _fields_ = [("abi_version", ctypes.c_uint32), ("struct_size", ctypes.c_uint32),
                ("smsd_version", ctypes.c_char_p), ("name", ctypes.c_char_p),
                ("device_info", _DEVICE_INFO), ("label_overlaps", _LABEL_OVERLAPS)]


def _write(buffer, size, text):
    data = text.encode()[:size - 1] + b"\0"
    ctypes.memmove(buffer, data, len(data))


class FakeAccelerator:
    """The accelerator interface, computing on the CPU; fail=True reports a failure on every call."""

    def __init__(self, version=None, abi=1, fail=False):
        self.calls = 0
        self.fail = fail
        self._version = (version or smsd.__version__).encode()
        self._name = b"Test"

        def device_info(text, size, error, error_size):
            _write(text, size, "Test device (compute 0.0, 0 MB)")
            return 0

        def label_overlaps(query, table_size, targets, offsets, entries, overlaps, error, error_size):
            self.calls += 1
            if self.fail:
                _write(error, error_size, "device lost")
                return 1
            for i in range(targets):
                overlaps[i] = sum(min(entries[e].count, query[entries[e].label])
                                  for e in range(offsets[i], offsets[i + 1]) if entries[e].label < table_size)
            return 0

        self._callbacks = (_DEVICE_INFO(device_info), _LABEL_OVERLAPS(label_overlaps))
        self.table = _Table(abi, ctypes.sizeof(_Table), self._version, self._name, *self._callbacks)

    def capsule(self, name=b"smsd.screen_accelerator"):
        new = ctypes.pythonapi.PyCapsule_New
        new.restype = ctypes.py_object
        new.argtypes = [ctypes.c_void_p, ctypes.c_char_p, ctypes.c_void_p]
        self._capsule_name = name
        return new(ctypes.addressof(self.table), name, None)

    def module(self, available=True, version=None):
        """A stand-in for the smsd_cuda package."""
        module = types.ModuleType("smsd_cuda")
        module.__version__ = version or smsd.__version__
        module.__smsd_version__ = version or smsd.__version__
        module.consulted = 0

        def is_available():
            module.consulted += 1
            return available

        module.is_available = is_available
        module.capsule = self.capsule
        return module


@pytest.fixture(autouse=True)
def fresh(monkeypatch):
    """Each test starts with the GPU enabled, no smsd_cuda installed and nothing attached."""
    monkeypatch.delenv("SMSD_DISABLE_GPU", raising=False)
    monkeypatch.setitem(sys.modules, "smsd_cuda", None)
    _accelerator._reset()
    yield
    _accelerator._reset()


def corpus():
    graphs = [smsd.parse_smiles(text) for text in SMILES]
    targets = smsd.TargetCorpus()
    targets.add_targets(graphs * 20)  # several hundred, so the scores are formed on several threads
    return graphs, targets


def expected_hits(query, graphs, threshold, chemistry):
    """The CPU hits, from the pairwise bound."""
    return [i for i, target in enumerate(graphs * 20)
            if smsd.similarity_upper_bound(query, target, chemistry) >= threshold]


def chemistries():
    strict = smsd.ChemOptions()
    strict.aromaticity_mode = smsd.AromaticityMode.STRICT
    strict.ring_matches_ring_only = True
    untyped = smsd.ChemOptions()
    untyped.match_atom_type = False
    return [smsd.ChemOptions(), strict, untyped]


def screens(graphs, targets):
    """Every query, chemistry and boundary threshold, with the hits of each screen."""
    for chemistry in chemistries():
        for query in graphs[::3]:
            for threshold in sorted({0.0, 1.0} | {smsd.similarity_upper_bound(query, t, chemistry) for t in graphs}):
                yield query, chemistry, threshold, smsd.screen_targets(query, targets, threshold, chemistry)


def test_without_smsd_cuda_nothing_changes():
    assert smsd.gpu_is_available() is _smsd.gpu_is_available()
    assert smsd.gpu_device_info() == _smsd.gpu_device_info()
    graphs, targets = corpus()
    for query, chemistry, threshold, hits in screens(graphs, targets):
        assert hits == expected_hits(query, graphs, threshold, chemistry)
        assert hits == smsd.screen_targets(query, graphs * 20, threshold, chemistry)
    assert _smsd._screen_accelerator() is None


def test_attached_accelerator_gives_the_cpu_hits():
    fake = FakeAccelerator()
    _smsd._attach_screen_accelerator(fake.capsule())
    assert _smsd._screen_accelerator() == ("Test", "Test device (compute 0.0, 0 MB)")
    graphs, targets = corpus()
    failures = _smsd._screen_accelerator_failures()[0]
    for query, chemistry, threshold, hits in screens(graphs, targets):
        assert hits == expected_hits(query, graphs, threshold, chemistry)
    assert fake.calls > 0
    assert _smsd._screen_accelerator_failures()[0] == failures


def test_failing_accelerator_falls_back_to_the_cpu():
    fake = FakeAccelerator(fail=True)
    _smsd._attach_screen_accelerator(fake.capsule())
    graphs, targets = corpus()
    failures = _smsd._screen_accelerator_failures()[0]
    for query, chemistry, threshold, hits in screens(graphs, targets):
        assert hits == expected_hits(query, graphs, threshold, chemistry)
    count, reason = _smsd._screen_accelerator_failures()
    assert count == failures + fake.calls and fake.calls > 0
    assert reason == "device lost"


def test_installed_package_is_found_and_reported(monkeypatch):
    fake = FakeAccelerator()
    module = fake.module()
    monkeypatch.setitem(sys.modules, "smsd_cuda", module)
    assert smsd.gpu_is_available()
    native = _smsd.gpu_device_info()
    if native.startswith("CPU: "):
        assert smsd.gpu_device_info() == f"GPU: Test device (compute 0.0, 0 MB)  [{native[len('CPU: '):]}]"
    graphs, targets = corpus()
    assert smsd.screen_targets(graphs[1], targets, 0.5) == expected_hits(graphs[1], graphs, 0.5, smsd.ChemOptions())
    assert fake.calls > 0 and module.consulted == 1


def test_package_without_a_device_is_left_alone(monkeypatch):
    module = FakeAccelerator().module(available=False)
    monkeypatch.setitem(sys.modules, "smsd_cuda", module)
    assert smsd.gpu_is_available() is _smsd.gpu_is_available()
    assert smsd.gpu_device_info() == _smsd.gpu_device_info()
    assert _smsd._screen_accelerator() is None and module.consulted == 1


def test_package_for_another_smsd_is_refused(monkeypatch):
    module = FakeAccelerator().module(version="0.0.0")
    monkeypatch.setitem(sys.modules, "smsd_cuda", module)
    with pytest.warns(RuntimeWarning, match="smsd-cuda 0.0.0 serves smsd 0.0.0"):
        assert smsd.gpu_is_available() is _smsd.gpu_is_available()
    assert _smsd._screen_accelerator() is None and module.consulted == 0


def test_disabled_gpu_does_not_look_for_the_package(monkeypatch):
    fake = FakeAccelerator()
    module = fake.module()
    monkeypatch.setitem(sys.modules, "smsd_cuda", module)
    monkeypatch.setenv("SMSD_DISABLE_GPU", "1")
    assert smsd.gpu_is_available() is _smsd.gpu_is_available()
    assert smsd.gpu_device_info() == _smsd.gpu_device_info()
    graphs, targets = corpus()
    smsd.screen_targets(graphs[1], targets, 0.5)
    assert module.consulted == 0 and fake.calls == 0
    _smsd._attach_screen_accelerator(fake.capsule())
    assert _smsd._screen_accelerator() is None
    smsd.screen_targets(graphs[1], targets, 0.5)
    assert fake.calls == 0


@pytest.mark.parametrize("fake, name, message", [
    (FakeAccelerator(abi=2), b"smsd.screen_accelerator", "interface version 2"),
    (FakeAccelerator(version="0.0.0"), b"smsd.screen_accelerator", "built for smsd 0.0.0"),
    (FakeAccelerator(), b"another.capsule", "does not hold a screening accelerator"),
])
def test_unsuitable_tables_are_refused(fake, name, message):
    with pytest.raises(ValueError, match=message):
        _smsd._attach_screen_accelerator(fake.capsule(name))
    assert _smsd._screen_accelerator() is None


def test_smsd_cuda_is_pinned_to_this_smsd():
    package = pathlib.Path(__file__).resolve().parents[2] / "packaging" / "smsd-cuda"
    if not package.is_dir():
        pytest.skip("the smsd-cuda sources are not part of this distribution")
    project = (package / "pyproject.toml").read_text(encoding="utf-8")
    module = (package / "python" / "smsd_cuda" / "__init__.py").read_text(encoding="utf-8")
    version = smsd.__version__
    assert re.search(r'^version = "([^"]+)"', project, re.M).group(1) == version
    assert re.findall(r'"smsd==([^"]+)"', project) == [version]
    assert re.search(r'^__version__ = "([^"]+)"', module, re.M).group(1) == version
    assert re.search(r'^__smsd_version__ = "([^"]+)"', module, re.M).group(1) == version
    root = package.parents[1]
    assert (package / "LICENSE").read_bytes() == (root / "LICENSE").read_bytes()
