# SPDX-License-Identifier: Apache-2.0
"""A wheel target cannot lower the minimum OS of its bundled OpenMP library."""

import importlib.util
from pathlib import Path
import subprocess

import pytest

SPEC = importlib.util.spec_from_file_location("macos_openmp", Path(__file__).parents[1] / "check_macos_openmp.py")
OPENMP = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(OPENMP)


def test_minimum_version_load_commands():
    commands = """Load command 8
          cmd LC_BUILD_VERSION
      cmdsize 32
     platform 1
        minos 14.0
          sdk 26.0
Load command 9
          cmd LC_LOAD_DYLIB
      version 99.0
"""
    assert OPENMP.minimum_versions(commands) == [(14, 0, 0)]
    assert OPENMP.minimum_versions("cmd LC_VERSION_MIN_MACOSX\nversion 11.0\nsdk 14.0\n") == [(11, 0, 0)]
    with pytest.raises(ValueError, match="no macOS"):
        OPENMP.minimum_versions("cmd LC_LOAD_DYLIB\nversion 14.0\n")


@pytest.mark.parametrize("minimum, accepted", [("11.0", True), ("14.0", True),
                                               ("14.0.1", False), ("26.0", False)])
def test_runtime_minimum_must_reach_declared_target(tmp_path, monkeypatch, minimum, accepted):
    (tmp_path / "lib").mkdir()
    (tmp_path / "include").mkdir()
    (tmp_path / "lib/libomp.dylib").write_bytes(b"fixture")
    (tmp_path / "include/omp.h").write_text("fixture")

    def output(command, **kwargs):
        if command[0] == "lipo":
            return "arm64\n"
        assert command[:4] == ["otool", "-l", "-arch", "arm64"]
        return "cmd LC_BUILD_VERSION\nminos " + minimum + "\nsdk 26.0\n"

    monkeypatch.setattr(subprocess, "check_output", output)
    if accepted:
        OPENMP.check(tmp_path, "14.0")
    else:
        with pytest.raises(ValueError, match="newer macOS"):
            OPENMP.check(tmp_path, "14.0")
