# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
# Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
# See the NOTICE file for attribution, trademark, and algorithm IP terms.
"""Apple Metal fingerprint and similarity-bound screening for SMSD Pro.

Installed beside smsd of the same version on an Apple silicon Mac, it lets smsd screen
TargetCorpus libraries on the Apple GPU with the same hits as on the CPU. smsd finds it at run
time; ``smsd.gpu_device_info()`` names the device, and SMSD_DISABLE_GPU=1 keeps screening on the
CPU.
"""

import ctypes
import os
import sys

__version__ = "8.2.1"
# The smsd version this package serves; smsd attaches it only when the two are equal.
__smsd_version__ = "8.2.1"

__all__ = ["capsule", "device_info", "is_available", "preload", "unavailable_reason"]

_HERE = os.path.dirname(os.path.abspath(__file__))
_LIBRARY = "libsmsd_metal_screen.dylib"
_CAPSULE_NAME = b"smsd.screen_accelerator"

_library = None
_reason = None


def _load():
    """The screening library, loaded once; None when it cannot be, with the reason kept."""
    global _library, _reason
    if _library is None and _reason is None:
        if sys.platform != "darwin":
            _reason = "Metal screening needs macOS on Apple silicon"
            return None
        try:
            library = ctypes.CDLL(os.path.join(_HERE, _LIBRARY))
        except OSError as error:
            _reason = f"the Metal screening library could not be loaded: {error}"
        else:
            library.smsd_metal_screen_accelerator.restype = ctypes.c_void_p
            library.smsd_metal_screen_accelerator.argtypes = []
            library.smsd_metal_device_available.restype = ctypes.c_int
            library.smsd_metal_device_available.argtypes = [ctypes.c_char_p, ctypes.c_size_t]
            library.smsd_metal_device_info.restype = ctypes.c_int
            library.smsd_metal_device_info.argtypes = [ctypes.c_char_p, ctypes.c_size_t,
                                                       ctypes.c_char_p, ctypes.c_size_t]
            library.smsd_metal_preload.restype = ctypes.c_int
            library.smsd_metal_preload.argtypes = [ctypes.c_char_p, ctypes.c_size_t]
            _library = library
    return _library


def unavailable_reason():
    """Why no Metal device can be used, or None when one can.

    Example::

        >>> import smsd_metal
        >>> reason = smsd_metal.unavailable_reason()
        >>> assert reason is None or isinstance(reason, str)
    """
    library = _load()
    if library is None:
        return _reason
    error = ctypes.create_string_buffer(256)
    if library.smsd_metal_device_available(error, len(error)) == 1:
        return None
    return error.value.decode(errors="replace") or "no Metal device was found"


def is_available():
    """Whether a Metal device can be used: this is an Apple silicon Mac with a GPU.

    Example::

        >>> import smsd_metal
        >>> assert isinstance(smsd_metal.is_available(), bool)
    """
    return unavailable_reason() is None


def device_info():
    """A description of the device, such as 'Apple M2 Pro (Metal, 21845 MB)'.

    The figure is the memory the system recommends the GPU work within.

    Example::

        >>> import smsd_metal
        >>> if smsd_metal.is_available():
        ...     assert "Metal" in smsd_metal.device_info()
    """
    reason = unavailable_reason()
    if reason is not None:
        raise RuntimeError(reason)
    text, error = ctypes.create_string_buffer(256), ctypes.create_string_buffer(256)
    if _library.smsd_metal_device_info(text, len(text), error, len(error)) != 0:
        raise RuntimeError(error.value.decode(errors="replace"))
    return text.value.decode(errors="replace")


def capsule():
    """A capsule holding the accelerator's table, which smsd attaches.

    Example::

        >>> import smsd_metal
        >>> if smsd_metal.is_available():
        ...     handle = smsd_metal.capsule()
        ...     assert type(handle).__name__ == "PyCapsule"
    """
    library = _load()
    if library is None:
        raise RuntimeError(_reason)
    new = ctypes.pythonapi.PyCapsule_New
    new.restype = ctypes.py_object
    new.argtypes = [ctypes.c_void_p, ctypes.c_char_p, ctypes.c_void_p]
    return new(library.smsd_metal_screen_accelerator(), _CAPSULE_NAME, None)


def preload():
    """Compile the Metal screening kernels without screening data.

    The first screen otherwise compiles them. This does not upload a corpus or change
    automatic dispatch.

    Example::

        >>> import smsd_metal
        >>> if smsd_metal.is_available():
        ...     assert smsd_metal.preload() is None
    """
    library = _load()
    if library is None:
        raise RuntimeError(_reason)
    error = ctypes.create_string_buffer(256)
    if library.smsd_metal_preload(error, len(error)) != 0:
        raise RuntimeError(error.value.decode(errors="replace"))
