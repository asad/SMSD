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

__version__ = "8.0.2"
# The smsd version this package serves; smsd attaches it only when the two are equal.
__smsd_version__ = "8.0.2"

__all__ = ["capsule", "device_info", "is_available", "unavailable_reason"]

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
            library.smsd_cuda_device_info.restype = ctypes.c_int
            library.smsd_cuda_device_info.argtypes = [ctypes.c_char_p, ctypes.c_size_t,
                                                      ctypes.c_char_p, ctypes.c_size_t]
            _library = library
    return _library


def unavailable_reason():
    """Why no CUDA device can be used, or None when one can."""
    library = _load()
    if library is None:
        return _reason
    error = ctypes.create_string_buffer(256)
    count = library.smsd_cuda_device_count(error, len(error))
    if count < 0:
        return error.value.decode(errors="replace") or "the CUDA runtime gave no device count"
    return None if count > 0 else "no CUDA device was found"


def is_available():
    """Whether a CUDA device can be used."""
    return unavailable_reason() is None


def device_info():
    """A description of the device, such as 'Tesla T4 (compute 7.5, 14930 MB)'."""
    reason = unavailable_reason()
    if reason is not None:
        raise RuntimeError(reason)
    text, error = ctypes.create_string_buffer(256), ctypes.create_string_buffer(256)
    if _library.smsd_cuda_device_info(text, len(text), error, len(error)) != 0:
        raise RuntimeError(error.value.decode(errors="replace"))
    return text.value.decode(errors="replace")


def capsule():
    """A capsule holding the accelerator's table, which smsd attaches."""
    library = _load()
    if library is None:
        raise RuntimeError(_reason)
    new = ctypes.pythonapi.PyCapsule_New
    new.restype = ctypes.py_object
    new.argtypes = [ctypes.c_void_p, ctypes.c_char_p, ctypes.c_void_p]
    return new(library.smsd_cuda_screen_accelerator(), _CAPSULE_NAME, None)
