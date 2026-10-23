# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
# Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
# See the NOTICE file for attribution, trademark, and algorithm IP terms.
"""smsd-cuda on a machine with an NVIDIA GPU: smsd attaches it and screens with the CPU hits."""
import ctypes
import os

import pytest

smsd = pytest.importorskip("smsd")
smsd_cuda = pytest.importorskip("smsd_cuda")
from smsd import _accelerator, _smsd  # noqa: E402

pytestmark = pytest.mark.skipif(not smsd_cuda.is_available(), reason=str(smsd_cuda.unavailable_reason()))

SMILES = ["c1ccccc1", "c1ccc(O)cc1", "CCO", "CC(=O)O", "C1CCCCC1", "c1ccncc1", "O=C(O)c1ccccc1O",
          "CN1C=NC2=C1C(=O)N(C(=O)N2C)C", "[NH4+]", "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O", "c1ccc2ccccc2c1",
          "ClC(Cl)Cl", "N[C@@H](C)C(=O)O", "c1cc[nH]c1", "FC(F)(F)c1ccccc1", "OC1C(O)C(O)C(O)C(O)C1O"]


@pytest.fixture(autouse=True)
def gpu_enabled(monkeypatch):
    monkeypatch.delenv("SMSD_DISABLE_GPU", raising=False)
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


@pytest.mark.parametrize("operation", [1, 2, 3])
def test_automatic_screening_minimum(operation):
    minimum_callback = ctypes.CFUNCTYPE(ctypes.c_int, ctypes.c_uint32,
                                        ctypes.POINTER(ctypes.c_uint64), ctypes.c_void_p, ctypes.c_size_t)

    class Table(ctypes.Structure):
        _fields_ = [("abi_version", ctypes.c_uint32), ("struct_size", ctypes.c_uint32),
                    ("smsd_version", ctypes.c_char_p), ("name", ctypes.c_char_p),
                    ("device_info", ctypes.c_void_p), ("label_overlaps", ctypes.c_void_p),
                    ("minimum_targets", minimum_callback)]

    assert smsd.gpu_is_available()
    table = Table.from_address(smsd_cuda._library.smsd_cuda_screen_accelerator())
    assert table.abi_version == 1 and table.struct_size >= ctypes.sizeof(Table)
    minimum = ctypes.c_uint64()
    error = ctypes.create_string_buffer(256)
    assert table.minimum_targets(operation, ctypes.byref(minimum), error, len(error)) == 0, error.value
    expected = 100000 if operation == 1 and smsd_cuda.device_info().startswith("Tesla T4 (compute 7.5,") else (1 << 64) - 1
    assert minimum.value == expected


@pytest.mark.parametrize("strict", [False, True])
def test_gpu_hits_equal_cpu_hits(monkeypatch, strict):
    chemistry = smsd.ChemOptions()
    if strict:
        chemistry.aromaticity_mode = smsd.AromaticityMode.STRICT
        chemistry.ring_matches_ring_only = True
    graphs = [smsd.parse_smiles(text) for text in SMILES]
    corpus = smsd.TargetCorpus()
    corpus.add_targets(graphs * 64)
    assert smsd.gpu_is_available()
    failures = _smsd._screen_accelerator_failures()[0]
    for query in graphs:
        thresholds = {0.0, 1.0} | {smsd.similarity_upper_bound(query, target, chemistry) for target in graphs}
        for threshold in sorted(thresholds):
            gpu = smsd.screen_targets(query, corpus, threshold, chemistry)
            monkeypatch.setenv("SMSD_DISABLE_GPU", "1")
            cpu = smsd.screen_targets(query, corpus, threshold, chemistry)
            monkeypatch.delenv("SMSD_DISABLE_GPU")
            assert gpu == cpu
    assert _smsd._screen_accelerator_failures()[0] == failures, _smsd._screen_accelerator_failures()[1]


def test_disabled_gpu_is_not_reported(monkeypatch):
    monkeypatch.setenv("SMSD_DISABLE_GPU", "1")
    assert not smsd.gpu_is_available() or _smsd.gpu_is_available()
    assert _smsd._screen_accelerator() is None
    assert os.environ["SMSD_DISABLE_GPU"] == "1"
