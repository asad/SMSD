# SPDX-License-Identifier: Apache-2.0
"""Build the CUDA conda package without changing the kernel or loader source."""
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
package = root / "packaging/smsd-cuda"
project_path = package / "pyproject.toml"
project = tomllib.loads(project_path.read_text())
core = tomllib.loads((root / "pyproject.toml").read_text())
if project["project"]["version"] != version or core["project"]["version"] != version:
    raise RuntimeError("Conda, core and CUDA source versions must match")
if f"smsd=={version}" not in project["project"]["dependencies"]:
    raise RuntimeError("CUDA Python metadata must require this exact core version")
native = re.search(r"^project\(smsd VERSION ([0-9.]+)\b", (root / "cpp/CMakeLists.txt").read_text(), re.M)
if native is None or native[1] != version:
    raise RuntimeError("CUDA and native source versions must match")

# Conda owns the native runtime dependency; it is not an installed PyPI distribution.
requirement = '    "nvidia-cuda-runtime-cu12>=12.8,<13",\n'
text = project_path.read_text()
if text.count(requirement) != 1:
    raise RuntimeError("Review the changed CUDA runtime metadata before adapting this recipe")
project_path.write_text(text.replace(requirement, "", 1))
os.environ["CMAKE_GENERATOR"] = "Ninja"
os.environ["CMAKE_BUILD_PARALLEL_LEVEL"] = "1"
subprocess.run([sys.executable, "-m", "pip", "install", ".", "--no-deps", "--no-build-isolation",
                "--no-index", "--no-cache-dir", "-vv",
                "--config-settings=build.tool-args=-j1",
                "--config-settings=cmake.define.SMSD_CUDA_RUNTIME_MINIMUM=12.8",
                "--config-settings=cmake.define.SMSD_CUDA_BINARY_MMA_PROTOTYPE=OFF"],
               cwd=package, check=True)
