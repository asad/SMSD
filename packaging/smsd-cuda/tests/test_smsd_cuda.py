# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
# Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
# See the NOTICE file for attribution, trademark, and algorithm IP terms.
"""smsd-cuda on a machine with an NVIDIA GPU: smsd attaches it and applies its automatic minimums.

test_gpu_cpu_parity.py compares the screens that reach the device with the same screens on the CPU.
"""
import ctypes
import os

import pytest

smsd = pytest.importorskip("smsd")
smsd_cuda = pytest.importorskip("smsd_cuda")
from smsd import _accelerator, _smsd  # noqa: E402

pytestmark = pytest.mark.skipif(not smsd_cuda.is_available(), reason=str(smsd_cuda.unavailable_reason()))


@pytest.fixture(autouse=True)
def gpu_enabled(monkeypatch):
    monkeypatch.delenv("SMSD_DISABLE_GPU", raising=False)
    monkeypatch.delenv("SMSD_CUDA_MIN_TARGETS", raising=False)
    _accelerator._reset()
    yield
    _accelerator._reset()


def test_versions_are_pinned_to_smsd():
    assert smsd_cuda.__smsd_version__ == smsd.__version__ == smsd_cuda.__version__


def test_smsd_attaches_the_device():
    assert smsd.gpu_is_available()
    info = smsd.gpu_device_info()
    assert info.startswith("GPU: ") and smsd_cuda.device_info() in info
    assert _smsd._screen_accelerator()[0] == "CUDA"


def _minimum_targets(operation):
    minimum_callback = ctypes.CFUNCTYPE(ctypes.c_int, ctypes.c_uint32,
                                        ctypes.POINTER(ctypes.c_uint64), ctypes.c_void_p, ctypes.c_size_t)

    class Table(ctypes.Structure):
        _fields_ = [("abi_version", ctypes.c_uint32), ("struct_size", ctypes.c_uint32),
                    ("smsd_version", ctypes.c_char_p), ("name", ctypes.c_char_p),
                    ("device_info", ctypes.c_void_p), ("label_overlaps", ctypes.c_void_p),
                    ("minimum_targets", minimum_callback)]

    assert smsd_cuda.is_available()
    table = Table.from_address(smsd_cuda._load().smsd_cuda_screen_accelerator())
    assert table.abi_version == 1 and table.struct_size >= ctypes.sizeof(Table)
    minimum = ctypes.c_uint64()
    error = ctypes.create_string_buffer(256)
    assert table.minimum_targets(operation, ctypes.byref(minimum), error, len(error)) == 0, error.value
    return minimum.value


@pytest.mark.parametrize("operation", [1, 2, 3])
def test_automatic_screening_minimum(operation):
    assert smsd.gpu_is_available()
    expected = 100000 if operation == 1 and smsd_cuda.device_info().startswith("Tesla T4 (compute 7.5,") else (1 << 64) - 1
    assert _minimum_targets(operation) == expected


@pytest.mark.parametrize("value, expected", [("1", 1), ("250000", 250000), ("0", None), ("-5", None),
                                             ("5k", None), ("18446744073709551616", None)])
def test_minimum_targets_setting(monkeypatch, value, expected):
    default = 100000 if smsd_cuda.device_info().startswith("Tesla T4 (compute 7.5,") else (1 << 64) - 1
    monkeypatch.setenv("SMSD_CUDA_MIN_TARGETS", value)
    assert _minimum_targets(1) == (default if expected is None else expected)
    assert _minimum_targets(2) == _minimum_targets(3) == (1 << 64) - 1
    monkeypatch.setenv("SMSD_DISABLE_GPU", "1")
    assert _minimum_targets(1) == (1 << 64) - 1


def test_supported_device_is_reported():
    assert smsd_cuda._load().smsd_cuda_device_supported(None, 0) == 1


def test_disabled_gpu_is_not_reported(monkeypatch):
    monkeypatch.setenv("SMSD_DISABLE_GPU", "1")
    assert not smsd.gpu_is_available() or _smsd.gpu_is_available()
    assert _smsd._screen_accelerator() is None
    assert os.environ["SMSD_DISABLE_GPU"] == "1"
