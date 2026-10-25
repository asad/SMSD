# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
# Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
# See the NOTICE file for attribution, trademark, and algorithm IP terms.
"""Chemical graph I/O, fingerprints, MCS and substructure searches.

    import smsd

    result = smsd.Substructure.find("c1ccccc1").match("c1ccc(O)cc1")
    assert result.found

    result = smsd.MCS.find("c1ccc(O)cc1").match("c1ccc(N)cc1")
    if result.found:
        print(result.best.atom_mapping, result.optimal)
"""

import importlib.machinery as _machinery
import importlib.util as _util
import os as _os
import sys as _sys

__version__ = "8.0.3"
__author__ = "Syed Asad Rahman"


def _load_native_module():
    """Load the extension beside this file, so a source tree never binds a separately installed one."""
    package = _os.path.dirname(_os.path.realpath(__file__))
    for suffix in _machinery.EXTENSION_SUFFIXES:
        candidate = _os.path.join(package, f"_smsd{suffix}")
        if not _os.path.exists(candidate):
            continue
        spec = _util.spec_from_file_location("smsd._smsd", candidate)
        if spec is None or spec.loader is None:
            continue
        module = _util.module_from_spec(spec)
        _sys.modules["smsd._smsd"] = module
        try:
            spec.loader.exec_module(module)
        except BaseException:
            _sys.modules.pop("smsd._smsd", None)
            raise
        return module

    from smsd import _smsd as module  # the installed extension
    return module


_smsd = _load_native_module()
if getattr(_smsd, "__version__", None) != __version__:
    # A native extension left over from another build would otherwise fail later and obscurely.
    _sys.modules.pop("smsd._smsd", None)
    raise ImportError(
        f"smsd {__version__} found a native extension built for "
        f"{getattr(_smsd, '__version__', 'an unknown version')} at {_smsd.__file__}; reinstall smsd")

# Native types and options.
ChemOptions = _smsd.ChemOptions
MolGraph = _smsd.MolGraph
MolGraphBuilder = _smsd.MolGraphBuilder

BondOrderMode = _smsd.BondOrderMode
AromaticityMode = _smsd.AromaticityMode
AromaticityModel = _smsd.AromaticityModel
RingFusionMode = _smsd.RingFusionMode
Solvent = _smsd.Solvent

SearchStatus = _smsd.SearchStatus
CancellationToken = _smsd.CancellationToken

FingerprintOptions = _smsd.FingerprintOptions
FingerprintFamily = _smsd.FingerprintFamily
FingerprintKey = _smsd.FingerprintKey
FragmentMode = _smsd.FragmentMode
ScaffoldMode = _smsd.ScaffoldMode

_NATIVE_TYPES = [
    "AromaticityMode", "AromaticityModel", "BondOrderMode", "CancellationToken", "ChemOptions",
    "FingerprintFamily", "FingerprintKey", "FingerprintOptions", "FragmentMode", "MolGraph",
    "MolGraphBuilder", "RingFusionMode", "ScaffoldMode", "SearchStatus", "Solvent",
]

# The public API, by purpose. The modules are private; each lists the names it exports.
from . import _chemistry, _depiction, _fingerprints, _interop, _layout, _molecule, _screening, search
from ._molecule import *
from ._interop import *
from ._chemistry import *
from ._fingerprints import *
from ._screening import *
from ._layout import *
from ._depiction import *
from .search import *

__all__ = sorted(_NATIVE_TYPES + _molecule.__all__ + _interop.__all__ + _chemistry.__all__
                 + _fingerprints.__all__ + _screening.__all__ + _layout.__all__ + _depiction.__all__
                 + search.__all__)
