# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
# Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
# See the NOTICE file for attribution, trademark, and algorithm IP terms.
"""smsd-cuda availability without a GPU: a stand-in library reports the device and its capability."""
import pytest

smsd_cuda = pytest.importorskip("smsd_cuda")

OLD_DEVICE = "the CUDA device Tesla P4 has compute capability 6.1; CUDA screening needs 7.5 or later"


class _Library:
    """Answers the device count and capability calls as the screening library does."""

    def __init__(self, count, supported, message=b""):
        self.count, self.supported, self.message = count, supported, message
        self.info_calls = 0

    def smsd_cuda_device_count(self, error, size):
        if self.count < 0:
            error.value = b"cudaGetDeviceCount: no CUDA-capable device is detected"
        return self.count

    def smsd_cuda_device_supported(self, error, size):
        error.value = self.message
        return self.supported

    def smsd_cuda_device_info(self, text, text_size, error, error_size):
        self.info_calls += 1
        text.value = b"Tesla T4 (compute 7.5, 14930 MB)"
        return 0


@pytest.fixture
def library(monkeypatch):
    def install(count, supported, message=b""):
        stand_in = _Library(count, supported, message)
        monkeypatch.setattr(smsd_cuda, "_library", stand_in)
        monkeypatch.setattr(smsd_cuda, "_reason", None)
        return stand_in
    return install


def test_supported_device_is_available(library):
    stand_in = library(1, 1)
    assert smsd_cuda.unavailable_reason() is None
    assert smsd_cuda.is_available()
    assert smsd_cuda.device_info() == "Tesla T4 (compute 7.5, 14930 MB)"
    assert stand_in.info_calls == 1


def test_device_below_compute_7_5_is_unavailable_with_its_reason(library):
    stand_in = library(1, 0, OLD_DEVICE.encode())
    assert smsd_cuda.unavailable_reason() == OLD_DEVICE
    assert not smsd_cuda.is_available()
    with pytest.raises(RuntimeError, match="compute capability 6.1"):
        smsd_cuda.device_info()
    assert stand_in.info_calls == 0


def test_unsupported_device_without_a_message_still_gives_a_reason(library):
    library(1, 0)
    assert smsd_cuda.unavailable_reason() == "the CUDA device needs compute capability 7.5 or later"


def test_capability_query_failure_is_reported(library):
    library(1, -1, b"compute capability: unknown error")
    assert smsd_cuda.unavailable_reason() == "compute capability: unknown error"
    library(1, -1)
    assert smsd_cuda.unavailable_reason() == "the CUDA runtime gave no compute capability"


def test_missing_device_is_reported_before_capability(library):
    library(0, 1)
    assert smsd_cuda.unavailable_reason() == "no CUDA device was found"
    library(-1, 1)
    assert smsd_cuda.unavailable_reason().startswith("cudaGetDeviceCount:")
