#!/usr/bin/env bash
# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
#
# Builds the C++ GPU tests and screen_compare with CUDA, the smsd package as published (CPU and
# OpenMP only) and the smsd-cuda wheel, and installs both packages into the environment.
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"
source_manifest="$(realpath "${SMSD_CUDA_SOURCE_MANIFEST:?Set the canonical source manifest}")"
sdk_source="$(realpath "${SMSD_CUDA_SDK_SOURCE:?Set the extracted canonical SDK directory}")"
release_dir="$(realpath "${SMSD_CUDA_RELEASE_DIR:-$(dirname "$source_manifest")}")"
[[ "$source_manifest" == "$release_dir/SOURCE_MANIFEST.json" ]] || {
    printf '%s\n' 'use the canonical SOURCE_MANIFEST.json in the release directory' >&2; exit 1;
}
results="$(cd "$RESULTS" && pwd)"
case "$results/" in "$release_dir/"*) printf '%s\n' 'keep private build/test results outside the release directory' >&2; exit 1 ;; esac
use_cuda
use_python
release_commit="$(python -c 'import json,sys; print(json.load(open(sys.argv[1]))["source_commit"])' "$source_manifest")"
source_archive="$(dirname "$source_manifest")/$(python -c 'import json,sys; print(json.load(open(sys.argv[1]))["sdist"])' "$source_manifest")"
check_release_source() {
    (cd "$ROOT" && python scripts/build_python_wheels.py check-source --sdist "$source_archive" \
        --source-manifest "$source_manifest" --release-commit "$release_commit")
    python - "$ROOT" "$source_archive" "$sdk_source" <<'PY'
import importlib.util
from pathlib import Path
import sys
spec = importlib.util.spec_from_file_location("build_wheels", Path(sys.argv[1]) / "scripts/build_python_wheels.py")
build = importlib.util.module_from_spec(spec)
spec.loader.exec_module(build)
expected = build.read_source_archive(Path(sys.argv[2]))
sdk = Path(sys.argv[3])
actual = {path.relative_to(sdk).as_posix() for path in sdk.rglob("*")
          if path.is_file() and "__pycache__" not in path.parts and path.suffix != ".pyc"}
if actual != set(expected) or any((sdk / name).read_bytes() != data for name, data in expected.items()):
    raise SystemExit("Extracted SDK differs from the canonical archive")
PY
}
check_release_source > "$RESULTS/source-before-build.log" 2>&1
wheel_platform="${SMSD_CUDA_WHEEL_PLATFORM:-manylinux_2_34_x86_64}"
[[ "$wheel_platform" == manylinux_2_34_x86_64 ]] || die "release CUDA wheels require manylinux_2_34_x86_64"
python - <<'PY'
import os
import platform
import re
import subprocess
if platform.system() != "Linux" or platform.machine() != "x86_64":
    raise SystemExit("CUDA release wheels require a Linux x86_64 build host")
libc = os.confstr("CS_GNU_LIBC_VERSION").split()
if len(libc) != 2 or libc[0] != "glibc" or tuple(map(int, libc[1].split("."))) > (2, 34):
    raise SystemExit("Build CUDA release wheels on glibc 2.34 or earlier; a newer host cannot be retagged")
toolkit = subprocess.check_output([os.environ["CUDACXX"], "--version"], text=True)
if not re.search(r"release 12\.8,", toolkit):
    raise SystemExit("The declared CUDA 12.8 runtime floor requires the CUDA 12.8 toolkit")
PY

log "build: C++ GPU tests (SMSD_BUILD_CUDA=ON)"
cmake -S "$sdk_source/cpp" -B "$WORK/cpp-cuda" -G Ninja -DCMAKE_BUILD_TYPE=Release -DSMSD_BUILD_TESTS=ON \
    -DSMSD_BUILD_CUDA=ON -DSMSD_BUILD_METAL=OFF -DSMSD_BUILD_PYTHON=OFF > "$RESULTS/build-cpp.log" 2>&1
cmake --build "$WORK/cpp-cuda" -j "$JOBS" --target "${GPU_TEST_TARGETS[@]}" >> "$RESULTS/build-cpp.log" 2>&1

log "build: screen_compare"
cmake -S "$PACKAGE_DIR/tools" -B "$WORK/tools" -G Ninja -DCMAKE_BUILD_TYPE=Release > "$RESULTS/build-tools.log" 2>&1
cmake --build "$WORK/tools" -j "$JOBS" --target screen_compare >> "$RESULTS/build-tools.log" 2>&1

log "build: smsd wheel, CPU and OpenMP only"
rm -rf "$RESULTS/smsd-wheels" "$RESULTS/repaired-wheels"
python -m pip wheel --no-deps "$sdk_source" -w "$RESULTS/smsd-wheels" \
    -Ccmake.define.SMSD_BUILD_CUDA=OFF -Ccmake.define.SMSD_BUILD_METAL=OFF \
    -Cbuild-dir="$WORK/smsd-build" > "$RESULTS/build-smsd.log" 2>&1
smsd_wheels=("$RESULTS"/smsd-wheels/smsd-*.whl)
[[ ${#smsd_wheels[@]} -eq 1 && -f ${smsd_wheels[0]} ]] || die "expected one smsd wheel"
python -m pip install --no-deps --force-reinstall "${smsd_wheels[0]}" >> "$RESULTS/build-smsd.log" 2>&1

log "build: smsd-cuda wheel"
rm -rf "$RESULTS/wheels"
python -m pip wheel --no-deps "$PACKAGE_DIR" -w "$RESULTS/wheels" -Cbuild-dir="$WORK/smsd-cuda-build" \
    "-Ccmake.define.SMSD_SOURCE_DIR=$sdk_source/cpp" \
    > "$RESULTS/build-smsd-cuda.log" 2>&1
wheels=("$RESULTS"/wheels/smsd_cuda-*.whl)
[[ ${#wheels[@]} -eq 1 && -f ${wheels[0]} ]] || die "expected one smsd-cuda wheel"
wheel="${wheels[0]}"
auditwheel show "$wheel" > "$RESULTS/wheel-before-repair.txt" 2>&1
# The separately declared NVIDIA runtime remains supplied by nvidia-cuda-runtime-cu12.
# Auditwheel also checks actual symbol requirements against the declared glibc 2.34 policy.
auditwheel repair --plat "$wheel_platform" --exclude libcudart.so.12 \
    --wheel-dir "$RESULTS/repaired-wheels" "$wheel" > "$RESULTS/auditwheel-repair.log" 2>&1
repaired_wheels=("$RESULTS"/repaired-wheels/smsd_cuda-*.whl)
[[ ${#repaired_wheels[@]} -eq 1 && -f ${repaired_wheels[0]} ]] || die "expected one repaired smsd-cuda wheel"
wheel="${repaired_wheels[0]}"
# The wheel replaces any earlier build; its dependencies come after, leaving the smsd just built in place.
python -m pip install --no-deps --force-reinstall "$wheel" >> "$RESULTS/build-smsd-cuda.log" 2>&1
python -m pip install "$wheel" "nvidia-cuda-runtime-cu12>=12.8,<12.9" >> "$RESULTS/build-smsd-cuda.log" 2>&1

{
    echo "wheel: $(basename "$wheel")"
    unzip -l "$wheel" || true
    library=$(python -c 'import os, smsd_cuda; print(os.path.join(os.path.dirname(smsd_cuda.__file__), "libsmsd_cuda_screen.so"))') || true
    echo "--- dynamic section"
    readelf -d "$library" | grep -E 'NEEDED|RUNPATH|RPATH' || true
    echo "--- resolved libraries"
    ldd "$library" || true
    echo "--- device code"
    cuobjdump --list-elf --list-ptx "$library" || true
    echo "--- auditwheel"
    auditwheel show "$wheel"
    echo "--- wheel tags and declared dependencies"
    python - "$wheel" <<'PY'
import sys
import zipfile

with zipfile.ZipFile(sys.argv[1]) as archive:
    for name in archive.namelist():
        if "libcudart" in name:
            raise SystemExit("the wheel must use the separately installed NVIDIA runtime")
        if name.endswith((".dist-info/WHEEL", ".dist-info/METADATA")):
            for line in archive.read(name).decode().splitlines():
                if line.startswith(("Tag:", "Requires-Dist:")):
                    print(line)
PY
} > "$RESULTS/wheel-report.txt" 2>&1
check_release_source > "$RESULTS/source-after-build.log" 2>&1
python - "$ROOT" "$source_manifest" "$sdk_source" "${smsd_wheels[0]}" "$wheel" "$WORK/cpp-cuda" \
    "$RESULTS/build-inputs.json" <<'PY'
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
import time
spec = importlib.util.spec_from_file_location("cuda_proof", Path(sys.argv[1]) / "scripts/cuda_wheel_provenance.py")
proof = importlib.util.module_from_spec(spec)
spec.loader.exec_module(proof)
source = proof.read_json(sys.argv[2])
native = Path(sys.argv[6]).resolve()
sdk = Path(sys.argv[3]).resolve()
for directory, home, cuda in ((native, sdk / "cpp", True),
                              (Path(os.environ["WORK"]) / "smsd-cuda-build", Path(sys.argv[1]) / "packaging/smsd-cuda", True),
                              (Path(os.environ["WORK"]) / "smsd-build", sdk / "cpp", False)):
    cache = {}
    for line in (directory / "CMakeCache.txt").read_text().splitlines():
        if "=" in line and ":" in line.split("=", 1)[0] and not line.startswith(("//", "#")):
            key, value = line.split("=", 1)
            cache[key.split(":", 1)[0]] = value
    if Path(cache["CMAKE_HOME_DIRECTORY"]).resolve() != home.resolve():
        raise SystemExit("CMake cache source differs from the declared frozen build input")
    if cuda and Path(cache["CMAKE_CUDA_COMPILER"]).resolve() != Path(os.environ["CUDACXX"]).resolve():
        raise SystemExit("CMake used a different CUDA compiler")
    if home.name == "smsd-cuda" and Path(cache["SMSD_SOURCE_DIR"]).resolve() != sdk / "cpp":
        raise SystemExit("The CUDA add-on did not use the canonical SDK headers")
note = {"source_commit": source["source_commit"], "sdist_sha256": source["sdist_sha256"],
        "source_manifest": str(Path(sys.argv[2]).resolve()), "sdk_source": str(Path(sys.argv[3]).resolve()),
        "checkout": str(Path(sys.argv[1]).resolve()), "cpu_wheel": str(Path(sys.argv[4]).resolve()),
        "cpu_wheel_sha256": proof.sha256(sys.argv[4]), "cuda_wheel": str(Path(sys.argv[5]).resolve()),
        "cuda_wheel_sha256": proof.sha256(sys.argv[5]), "native_dir": str(native),
        "native_binaries": {name: proof.sha256(native / name) for name in sorted(proof.NATIVE_TESTS)},
        "cuda_toolkit": subprocess.check_output([os.environ["CUDACXX"], "--version"], text=True),
        "finished_ns": time.time_ns()}
with Path(sys.argv[7]).open("x") as stream:
    json.dump(note, stream, indent=2)
    stream.write("\n")
PY
log "build: $(basename "$wheel") installed; $(python -c 'import smsd; print(smsd.gpu_device_info())')"
