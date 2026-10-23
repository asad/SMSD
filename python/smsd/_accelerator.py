# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
# Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
# See the NOTICE file for attribution, trademark, and algorithm IP terms.
"""Finds an installed screening accelerator, such as the smsd-cuda package, at run time.

The smsd wheels screen on the CPU. When smsd-cuda is installed for this exact smsd version and
finds a device, the first GPU query or TargetCorpus screen attaches it, and corpus screens then
compute their label overlaps on the GPU with the same hits. Without it, or with SMSD_DISABLE_GPU
set to a value other than 0, nothing changes.
"""

import importlib
import os
import threading
import warnings

from . import _smsd

# The packages that may provide an accelerator, in order of preference.
_PACKAGES = ("smsd_cuda",)

_lock = threading.Lock()
_searched = False


def _disabled():
    """Whether SMSD_DISABLE_GPU keeps screening on the CPU, by the rule the native module applies."""
    return os.environ.get("SMSD_DISABLE_GPU", "") not in ("", "0")


def _warn(message):
    warnings.warn(f"{message}; screening runs on the CPU", RuntimeWarning, stacklevel=5)


def _attach(name):
    """Attach the accelerator of the named package; whether it was attached."""
    try:
        module = importlib.import_module(name)
    except ImportError:
        return False
    except Exception as error:  # a broken installation must not stop smsd
        _warn(f"smsd could not load {name}: {error}")
        return False
    built_for = getattr(module, "__smsd_version__", None)
    if built_for != _smsd.__version__:
        distribution = name.replace("_", "-")
        _warn(f"{distribution} {getattr(module, '__version__', '')} serves smsd {built_for}, not smsd "
              f"{_smsd.__version__}; install {distribution}=={_smsd.__version__}")
        return False
    try:
        if not module.is_available():
            return False
        _smsd._attach_screen_accelerator(module.capsule())
    except Exception as error:
        _warn(f"smsd could not attach {name}: {error}")
        return False
    return True


def discover():
    """Search once for an installed accelerator and attach the first that has a device."""
    global _searched
    if _disabled() or _searched:
        return
    with _lock:
        if not _searched:
            _searched = True
            for name in _PACKAGES:
                if _attach(name):
                    break


def active():
    """The active accelerator as (name, device), after discover(); None without one."""
    discover()
    return None if _disabled() else _smsd._screen_accelerator()


def _reset():
    """Detach any accelerator and search again on next use."""
    global _searched
    with _lock:
        _searched = False
        _smsd._detach_screen_accelerator()
