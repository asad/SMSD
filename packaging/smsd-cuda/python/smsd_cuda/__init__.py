# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
# Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
# See the NOTICE file for attribution, trademark, and algorithm IP terms.
"""NVIDIA CUDA fingerprint and similarity-bound screening for SMSD Pro.

Installed beside smsd of the same version, it lets smsd screen eligible TargetCorpus libraries on
an NVIDIA GPU with the same hits as on the CPU. smsd finds it at run time;
``smsd.gpu_device_info()`` names the device, and SMSD_DISABLE_GPU=1 keeps screening on the CPU.
"""

import ctypes
import importlib.util
import os
import sys

__version__ = "8.1.0"
# The smsd version this package serves; smsd attaches it only when the two are equal.
__smsd_version__ = "8.1.0"

__all__ = ["capsule", "device_info", "is_available", "preload", "unavailable_reason"]

_HERE = os.path.dirname(os.path.abspath(__file__))
_LIBRARY = "smsd_cuda_screen.dll" if sys.platform == "win32" else "libsmsd_cuda_screen.so"
_RUNTIME = ("bin", "cudart64_12.dll") if sys.platform == "win32" else ("lib", "libcudart.so.12")
_CAPSULE_NAME = b"smsd.screen_accelerator"

_library = None
_reason = None


def _load_runtime():
    """Load the CUDA runtime of NVIDIA's nvidia-cuda-runtime-cu12 package, when it is installed."""
    try:
        spec = importlib.util.find_spec("nvidia.cuda_runtime")
    except (ImportError, ValueError):
        return
    for folder in (spec.submodule_search_locations or []) if spec is not None else []:
        path = os.path.join(folder, *_RUNTIME)
        if os.path.exists(path):
            ctypes.CDLL(path, mode=getattr(os, "RTLD_GLOBAL", 0) | getattr(os, "RTLD_NOW", 0))
            return


def _load():
    """The screening library, loaded once; None when it cannot be, with the reason kept."""
    global _library, _reason
    if _library is None and _reason is None:
        try:
            _load_runtime()
            library = ctypes.CDLL(os.path.join(_HERE, _LIBRARY))
        except OSError as error:
            _reason = f"the CUDA screening library could not be loaded: {error}"
        else:
            library.smsd_cuda_screen_accelerator.restype = ctypes.c_void_p
            library.smsd_cuda_screen_accelerator.argtypes = []
            library.smsd_cuda_device_count.restype = ctypes.c_int
            library.smsd_cuda_device_count.argtypes = [ctypes.c_char_p, ctypes.c_size_t]
            library.smsd_cuda_device_supported.restype = ctypes.c_int
            library.smsd_cuda_device_supported.argtypes = [ctypes.c_char_p, ctypes.c_size_t]
            library.smsd_cuda_device_info.restype = ctypes.c_int
            library.smsd_cuda_device_info.argtypes = [ctypes.c_char_p, ctypes.c_size_t,
                                                      ctypes.c_char_p, ctypes.c_size_t]
            library.smsd_cuda_preload.restype = ctypes.c_int
            library.smsd_cuda_preload.argtypes = [ctypes.c_char_p, ctypes.c_size_t]
            _library = library
    return _library


def unavailable_reason():
    """Why no CUDA device can be used, or None when one can.

    A device below compute capability 7.5 cannot be used; the reason names it and its capability.

    Example::

        >>> import smsd_cuda
        >>> reason = smsd_cuda.unavailable_reason()
        >>> assert reason is None or isinstance(reason, str)
    """
    library = _load()
    if library is None:
        return _reason
    error = ctypes.create_string_buffer(256)
    count = library.smsd_cuda_device_count(error, len(error))
    if count < 0:
        return error.value.decode(errors="replace") or "the CUDA runtime gave no device count"
    if count == 0:
        return "no CUDA device was found"
    error = ctypes.create_string_buffer(256)
    supported = library.smsd_cuda_device_supported(error, len(error))
    if supported > 0:
        return None
    if supported < 0:
        return error.value.decode(errors="replace") or "the CUDA runtime gave no compute capability"
    return error.value.decode(errors="replace") or "the CUDA device needs compute capability 7.5 or later"


def is_available():
    """Whether a CUDA device can be used: one is present and has compute capability 7.5 or later.

    Example::

        >>> import smsd_cuda
        >>> assert isinstance(smsd_cuda.is_available(), bool)
    """
    return unavailable_reason() is None


def device_info():
    """A description of the device, such as 'Tesla T4 (compute 7.5, 14930 MB)'.

    Example::

        >>> import smsd_cuda
        >>> if smsd_cuda.is_available():
        ...     assert "compute" in smsd_cuda.device_info()
    """
    reason = unavailable_reason()
    if reason is not None:
        raise RuntimeError(reason)
    text, error = ctypes.create_string_buffer(256), ctypes.create_string_buffer(256)
    if _library.smsd_cuda_device_info(text, len(text), error, len(error)) != 0:
        raise RuntimeError(error.value.decode(errors="replace"))
    return text.value.decode(errors="replace")


def capsule():
    """A capsule holding the accelerator's table, which smsd attaches.

    Example::

        >>> import smsd_cuda
        >>> if smsd_cuda.is_available():
        ...     handle = smsd_cuda.capsule()
        ...     assert type(handle).__name__ == "PyCapsule"
    """
    library = _load()
    if library is None:
        raise RuntimeError(_reason)
    new = ctypes.pythonapi.PyCapsule_New
    new.restype = ctypes.py_object
    new.argtypes = [ctypes.c_void_p, ctypes.c_char_p, ctypes.c_void_p]
    return new(library.smsd_cuda_screen_accelerator(), _CAPSULE_NAME, None)


def preload():
    """Initialise the current CUDA device and load screening kernels without screening data.

    This does not upload a corpus, reserve device scratch or change automatic dispatch.

    Example::

        >>> import smsd_cuda
        >>> if smsd_cuda.is_available():
        ...     assert smsd_cuda.preload() is None
    """
    library = _load()
    if library is None:
        raise RuntimeError(_reason)
    error = ctypes.create_string_buffer(256)
    if library.smsd_cuda_preload(error, len(error)) != 0:
        raise RuntimeError(error.value.decode(errors="replace"))
