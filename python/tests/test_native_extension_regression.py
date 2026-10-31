# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
# Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
# See the NOTICE file for attribution, trademark, and algorithm IP terms.
"""The native extension reports its version and lets other threads run during long calls."""
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import threading
import time

import pytest

import smsd
from smsd import _smsd


def test_extension_version_is_the_package_version():
    assert _smsd.__version__ == smsd.__version__


def test_import_rejects_an_extension_of_another_version(tmp_path):
    package = Path(smsd.__file__).resolve().parent
    copy = tmp_path / "smsd"
    shutil.copytree(package, copy, ignore=shutil.ignore_patterns("__pycache__"))
    # Wheels repaired by auditwheel or delvewheel keep their bundled libraries beside the
    # package, in smsd.libs; delocate keeps them inside it, in .dylibs, copied above.
    libraries = package.parent / "smsd.libs"
    if libraries.is_dir():
        shutil.copytree(libraries, tmp_path / "smsd.libs")
    init = copy / "__init__.py"
    text, count = re.subn(r'^__version__ = "[^"]+"$', '__version__ = "0.0.1"',
                          init.read_text(encoding="utf-8"), flags=re.M)
    assert count == 1
    init.write_text(text, encoding="utf-8")
    result = subprocess.run([sys.executable, "-c", "import smsd"], cwd=tmp_path, text=True,
                            capture_output=True, timeout=60,
                            env=dict(os.environ, PYTHONPATH=str(tmp_path)))
    assert result.returncode != 0
    assert "ImportError: smsd 0.0.1 found a native extension built for " + smsd.__version__ in result.stderr


def _main_thread_runs_during(call, attempts=10):
    """Whether this thread executes Python code in the middle half of call() on another thread.

    A call that holds the GIL blocks this thread until it returns, so no step of the loop below
    falls inside it; a call that releases the GIL lets the loop run throughout.
    """
    for _ in range(attempts):
        started, finished, span = threading.Event(), threading.Event(), []

        def worker():
            started.set()
            begin = time.perf_counter()
            try:
                call()
            finally:
                span.extend((begin, time.perf_counter()))
                finished.set()

        steps = []
        thread = threading.Thread(target=worker)
        thread.start()
        started.wait()
        while not finished.is_set():
            steps.append(time.perf_counter())
        thread.join()
        begin, end = span
        quarter = (end - begin) / 4
        if any(begin + quarter < step < end - quarter for step in steps):
            return True
    return False


@pytest.fixture(scope="module")
def workloads(tmp_path_factory):
    directory = tmp_path_factory.mktemp("released")
    records = [smsd.parse_smiles("c1ccc2ccccc2c1C(=O)NCCO")] * 4000
    sdf = str(directory / "records.sdf")
    _smsd.write_sdf(records, sdf)
    molecule = smsd.parse_smiles("c1ccccc1" + "C(C)(O)C(=O)N" * 600)
    chain, shorter = smsd.parse_smiles("C" * 1500), smsd.parse_smiles("C" * 600)
    identity = {atom: atom for atom in range(chain.n)}
    shorter_identity = {atom: atom for atom in range(shorter.n)}
    options = smsd.MCS.Options()
    return {
        "fingerprint_options": lambda: _smsd.fingerprint_options(molecule, smsd.FingerprintOptions.maccs()),
        "fingerprint_options counts": lambda: _smsd.fingerprint_options(
            molecule, smsd.FingerprintOptions(family="ecfp", counts=True, radius=6)),
        "read_sdf": lambda: _smsd.read_sdf(sdf),
        "write_sdf": lambda: _smsd.write_sdf(records, str(directory / "copy.sdf")),
        "validate_mapping": lambda: _smsd.MCS.validate_mapping(chain, chain, identity, options),
        "is_mapping_maximal": lambda: _smsd.MCS.is_mapping_maximal(shorter, shorter, shorter_identity, options),
    }


@pytest.mark.parametrize("name", ["fingerprint_options", "fingerprint_options counts", "read_sdf",
                                  "write_sdf", "validate_mapping", "is_mapping_maximal"])
def test_long_native_calls_release_the_gil(workloads, name):
    assert _main_thread_runs_during(workloads[name]), name
