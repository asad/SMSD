# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
# Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
# See the NOTICE file for attribution, trademark, and algorithm IP terms.
"""smsd-metal on an Apple silicon Mac: smsd attaches it and applies its automatic minimums.

test_gpu_cpu_parity.py compares the screens that reach the device with the same screens on the CPU.
"""
import ctypes
import os

import pytest

smsd = pytest.importorskip("smsd")
smsd_metal = pytest.importorskip("smsd_metal")
from smsd import _accelerator, _smsd  # noqa: E402

pytestmark = pytest.mark.skipif(not smsd_metal.is_available(), reason=str(smsd_metal.unavailable_reason()))

NO_AUTOMATIC_USE = (1 << 64) - 1


class _Table(ctypes.Structure):
    _fields_ = [("abi_version", ctypes.c_uint32), ("struct_size", ctypes.c_uint32),
                ("smsd_version", ctypes.c_char_p), ("name", ctypes.c_char_p),
                ("device_info", ctypes.c_void_p), ("label_overlaps", ctypes.c_void_p),
                ("minimum_targets", ctypes.CFUNCTYPE(ctypes.c_int, ctypes.c_uint32,
                                                     ctypes.POINTER(ctypes.c_uint64), ctypes.c_void_p,
                                                     ctypes.c_size_t)),
                ("create_fingerprints", ctypes.c_void_p), ("fingerprint_screen", ctypes.c_void_p),
                ("create_labels", ctypes.c_void_p), ("resident_label_overlaps", ctypes.c_void_p),
                ("destroy_resident", ctypes.c_void_p)]


@pytest.fixture(autouse=True)
def gpu_enabled(monkeypatch):
    monkeypatch.delenv("SMSD_DISABLE_GPU", raising=False)
    monkeypatch.delenv("SMSD_METAL_MIN_TARGETS", raising=False)
    _accelerator._reset()
    yield
    _accelerator._reset()


def _attach():
    _smsd._attach_screen_accelerator(smsd_metal.capsule())


def test_versions_are_pinned_to_smsd():
    assert smsd_metal.__smsd_version__ == smsd.__version__ == smsd_metal.__version__


def test_smsd_finds_the_package():
    if "smsd_metal" not in _accelerator._PACKAGES:
        pytest.skip("this smsd looks for its accelerators in smsd._accelerator._PACKAGES, "
                    "which does not name smsd_metal")
    assert _accelerator.active()[0] == "Metal"


def test_smsd_attaches_the_device():
    _attach()
    assert smsd.gpu_is_available()
    info = smsd.gpu_device_info()
    assert info.startswith("GPU: ") and smsd_metal.device_info() in info
    assert _smsd._screen_accelerator() == ("Metal", smsd_metal.device_info())
    assert "(Metal, " in smsd_metal.device_info()


def test_preload_compiles_the_kernels():
    assert smsd_metal.preload() is None
    assert smsd_metal.preload() is None


def test_capsule_holds_a_table_for_this_smsd():
    capsule = smsd_metal.capsule()
    assert type(capsule).__name__ == "PyCapsule"
    table = _Table.from_address(smsd_metal._load().smsd_metal_screen_accelerator())
    # The table carries every optional member.
    assert table.abi_version == 1 and table.struct_size == ctypes.sizeof(_Table)
    assert table.smsd_version == smsd.__version__.encode() and table.name == b"Metal"


def _minimum_targets(operation):
    table = _Table.from_address(smsd_metal._load().smsd_metal_screen_accelerator())
    assert table.struct_size == ctypes.sizeof(_Table)
    minimum = ctypes.c_uint64()
    error = ctypes.create_string_buffer(256)
    assert table.minimum_targets(operation, ctypes.byref(minimum), error, len(error)) == 0, error.value
    return minimum.value


def _measured_minimum():
    """The fingerprint minimum without a setting: 100,000 targets on the measured Apple M5."""
    return 100000 if smsd_metal.device_info().startswith("Apple M5 (Metal,") else NO_AUTOMATIC_USE


@pytest.mark.parametrize("operation", [1, 2, 3])
def test_automatic_screening_minimum(operation):
    assert _minimum_targets(operation) == (_measured_minimum() if operation == 1 else NO_AUTOMATIC_USE)


@pytest.mark.parametrize("value, expected", [("1", 1), ("250000", 250000), ("0", None), ("-5", None),
                                             ("5k", None), (" 7", None), ("18446744073709551616", None)])
def test_minimum_targets_setting(monkeypatch, value, expected):
    monkeypatch.setenv("SMSD_METAL_MIN_TARGETS", value)
    assert _minimum_targets(1) == (_measured_minimum() if expected is None else expected)
    assert _minimum_targets(2) == _minimum_targets(3) == NO_AUTOMATIC_USE
    monkeypatch.setenv("SMSD_DISABLE_GPU", "1")
    assert _minimum_targets(1) == NO_AUTOMATIC_USE


def test_unknown_operation_is_refused():
    table = _Table.from_address(smsd_metal._load().smsd_metal_screen_accelerator())
    minimum, error = ctypes.c_uint64(), ctypes.create_string_buffer(256)
    assert table.minimum_targets(4, ctypes.byref(minimum), error, len(error)) != 0
    assert error.value == b"unknown screening operation"


def test_disabled_gpu_is_not_reported(monkeypatch):
    _attach()
    monkeypatch.setenv("SMSD_DISABLE_GPU", "1")
    assert _smsd._screen_accelerator() is None
    assert not smsd.gpu_is_available() or _smsd.gpu_is_available()
    assert os.environ["SMSD_DISABLE_GPU"] == "1"


def test_small_corpus_screens_on_the_cpu_by_default():
    _attach()
    failures = _smsd._screen_accelerator_failures()[0]
    corpus = smsd.TargetCorpus.from_smiles(["CCO", "CCN", "CCC"])
    corpus.prewarm(num_threads=1)
    assert corpus.screen_many(["CCO", "CCN"], 1.0, num_threads=1) == [[0], [1]]
    assert corpus.screen_targets_many(["CCO", "CCN"], 1.0, num_threads=1) == [[0], [1]]
    assert _smsd._screen_accelerator_failures()[0] == failures
