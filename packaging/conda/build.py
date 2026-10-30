# SPDX-License-Identifier: Apache-2.0
"""Build the CPU conda package using the activated native compiler."""
import os
from pathlib import Path
import re
import subprocess
import sys

try:
    import tomllib
except ModuleNotFoundError:
    import tomli as tomllib

root = Path(os.environ["SRC_DIR"])
version = os.environ["PKG_VERSION"]
project = tomllib.loads((root / "pyproject.toml").read_text())
native = re.search(r"^project\(smsd VERSION ([0-9.]+)\b", (root / "cpp/CMakeLists.txt").read_text(), re.M)
if project["project"]["version"] != version or native is None or native[1] != version:
    raise RuntimeError("Conda, Python and native source versions must match")
os.environ["CMAKE_GENERATOR"] = "Ninja"
os.environ["CMAKE_BUILD_PARALLEL_LEVEL"] = "1"
subprocess.run([sys.executable, "-m", "pip", "install", ".", "--no-deps", "--no-build-isolation",
                "--no-index", "--no-cache-dir", "-vv",
                "--config-settings=build.tool-args=-j1",
                "--config-settings=cmake.define.SMSD_BUILD_TESTS=OFF",
                "--config-settings=cmake.define.SMSD_BUILD_OPENMP=OFF",
                "--config-settings=cmake.define.SMSD_BUILD_CUDA=OFF",
                "--config-settings=cmake.define.SMSD_BUILD_METAL=OFF"], cwd=root, check=True)
