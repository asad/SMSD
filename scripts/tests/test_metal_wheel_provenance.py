# SPDX-License-Identifier: Apache-2.0
"""A Metal wheel proof needs the complete native Metal suite and every Metal case passed on the device."""
import importlib.util
from pathlib import Path
import re

import pytest

SPEC = importlib.util.spec_from_file_location("metal_proof", Path(__file__).parents[1] / "metal_wheel_provenance.py")
METAL = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(METAL)

SOURCE = {"source_commit": "a" * 40, "sdist_sha256": "b" * 64, "version": "8.2.0"}
WHEEL = "smsd_metal-8.2.0-py3-none-macosx_14_0_arm64.whl"
# pytest names a suite outside its root directory by its relative path.
METAL_SUITE = "......checkout.packaging.smsd-metal.tests.test_smsd_metal"
NATIVE = [(name, name, False) for name in sorted(METAL.NATIVE_TESTS)]
PYTHON = [("test_smsd", "test_parse_smiles", False), *((METAL_SUITE, name, False) for name in sorted(METAL.DEVICE_CASES))]
IDENTITY = {"prefix": "/venv", "core_module": "/venv/smsd/__init__.py", "native_module": "/venv/smsd/_smsd.so",
            "module": "/venv/smsd_metal/__init__.py", "library": "/venv/smsd_metal/libsmsd_metal_screen.dylib",
            "versions": ["8.2.0"] * 4, "accelerator": ["Metal", "Apple M5 (Metal, 1024 MB)"]}


def junit(path, cases):
    """JUnit results of (class, name, skipped) cases."""
    rows = "".join(f'<testcase classname="{classname}" name="{name}">'
                   + ('<skipped message="no Metal device"/>' if skipped else "") + "</testcase>"
                   for classname, name, skipped in cases)
    path.write_text(f'<testsuites><testsuite name="suite" tests="{len(cases)}" '
                    f'skipped="{sum(skipped for *_, skipped in cases)}">{rows}</testsuite></testsuites>')
    return path


def proof(tmp_path, monkeypatch, native=NATIVE, python=PYTHON, **identity):
    checked = []
    monkeypatch.setattr(METAL.PROOF, "installed_files", lambda wheel, package, *paths: checked.append((wheel.name, package)))
    wheel, cpu_wheel = tmp_path / WHEEL, tmp_path / "smsd.whl"
    wheel.write_bytes(b"Metal wheel")
    result = METAL.proof_from_results(SOURCE, wheel, cpu_wheel, {**IDENTITY, **identity},
                                      junit(tmp_path / "native.xml", native), junit(tmp_path / "python.xml", python))
    assert checked == [("smsd.whl", "smsd"), (WHEEL, "smsd_metal")]
    return result


def test_complete_results_give_the_proof(tmp_path, monkeypatch):
    assert proof(tmp_path, monkeypatch) == {
        "schema_version": 1, "source_commit": "a" * 40, "sdist_sha256": "b" * 64, "smsd_version": "8.2.0",
        "wheel": WHEEL, "wheel_sha256": METAL.PROOF.sha256(tmp_path / WHEEL),
        "installed_wheel_tests": "passed", "native_gpu_tests": "passed"}


@pytest.mark.parametrize("change", [
    {"native": NATIVE[1:]},
    {"native": [*NATIVE[:-1], (*NATIVE[-1][:2], True)]},
    {"native": [*NATIVE, ("smsd_metal_unknown_tests", "smsd_metal_unknown_tests", False)]},
    {"python": [*PYTHON[:-1], (*PYTHON[-1][:2], True)]},
    {"python": PYTHON[:-1]},
    {"python": PYTHON[1:]},
    {"accelerator": None},
    {"accelerator": ["CUDA", "NVIDIA device"]},
    {"versions": ["8.2.0", "8.2.0", "8.2.0", "8.1.0"]},
], ids=["native case missing", "native case skipped", "unknown native case", "Metal case skipped",
        "device case missing", "no core suite", "no device", "another device", "another version"])
def test_incomplete_results_give_no_proof(tmp_path, monkeypatch, change):
    with pytest.raises(ValueError):
        proof(tmp_path, monkeypatch, **change)


def test_the_clean_checkout_check_covers_every_source_the_metal_build_reads():
    root = Path(__file__).parents[2]
    recipe = (root / "packaging/smsd-metal/CMakeLists.txt").read_text()
    # A path that ends a message sentence carries its full stop.
    paths = {"cpp/" + path.rstrip(".") for path in re.findall(r'\$\{SMSD_SOURCE_DIR\}/([\w./]+)', recipe)}
    assert {"cpp/metal/metal_fingerprint.mm", "cpp/include", "cpp/CMakeLists.txt"} <= paths
    assert all(any(path == source or path.startswith(source + "/") for source in METAL.METAL_SOURCES) for path in paths)
