#!/usr/bin/env bash
# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
#
# Builds the C++ GPU tests and screen_compare with CUDA, the smsd package as published (CPU and
# OpenMP only) and the smsd-cuda wheel, and installs both packages into the environment.
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"
use_cuda
use_python

log "build: C++ GPU tests (SMSD_BUILD_CUDA=ON)"
cmake -S "$ROOT/cpp" -B "$WORK/cpp-cuda" -G Ninja -DCMAKE_BUILD_TYPE=Release -DSMSD_BUILD_TESTS=ON \
    -DSMSD_BUILD_CUDA=ON -DSMSD_BUILD_METAL=OFF -DSMSD_BUILD_PYTHON=OFF > "$RESULTS/build-cpp.log" 2>&1
cmake --build "$WORK/cpp-cuda" -j "$JOBS" --target "${GPU_TEST_TARGETS[@]}" >> "$RESULTS/build-cpp.log" 2>&1

log "build: screen_compare"
cmake -S "$PACKAGE_DIR/tools" -B "$WORK/tools" -G Ninja -DCMAKE_BUILD_TYPE=Release > "$RESULTS/build-tools.log" 2>&1
cmake --build "$WORK/tools" -j "$JOBS" --target screen_compare >> "$RESULTS/build-tools.log" 2>&1

log "build: smsd wheel, CPU and OpenMP only"
rm -rf "$RESULTS/smsd-wheels" "$RESULTS/repaired-wheels"
python -m pip wheel --no-deps "$ROOT" -w "$RESULTS/smsd-wheels" \
    -Ccmake.define.SMSD_BUILD_CUDA=OFF -Ccmake.define.SMSD_BUILD_METAL=OFF \
    -Cbuild-dir="$WORK/smsd-build" > "$RESULTS/build-smsd.log" 2>&1
smsd_wheels=("$RESULTS"/smsd-wheels/smsd-*.whl)
[[ ${#smsd_wheels[@]} -eq 1 && -f ${smsd_wheels[0]} ]] || die "expected one smsd wheel"
python -m pip install --no-deps --force-reinstall "${smsd_wheels[0]}" >> "$RESULTS/build-smsd.log" 2>&1

log "build: smsd-cuda wheel"
rm -rf "$RESULTS/wheels"
python -m pip wheel --no-deps "$PACKAGE_DIR" -w "$RESULTS/wheels" -Cbuild-dir="$WORK/smsd-cuda-build" \
    > "$RESULTS/build-smsd-cuda.log" 2>&1
wheels=("$RESULTS"/wheels/smsd_cuda-*.whl)
[[ ${#wheels[@]} -eq 1 && -f ${wheels[0]} ]] || die "expected one smsd-cuda wheel"
wheel="${wheels[0]}"
auditwheel show "$wheel" > "$RESULTS/wheel-before-repair.txt" 2>&1
# The separately declared NVIDIA runtime remains supplied by nvidia-cuda-runtime-cu12.
# Use this host's glibc policy unless a build environment selects an older supported policy.
wheel_platform="${SMSD_CUDA_WHEEL_PLATFORM:-$(python -c 'import os, platform; print("manylinux_" + os.confstr("CS_GNU_LIBC_VERSION").split()[-1].replace(".", "_") + "_" + platform.machine())')}"
auditwheel repair --plat "$wheel_platform" --exclude libcudart.so.12 \
    --wheel-dir "$RESULTS/repaired-wheels" "$wheel" > "$RESULTS/auditwheel-repair.log" 2>&1
repaired_wheels=("$RESULTS"/repaired-wheels/smsd_cuda-*.whl)
[[ ${#repaired_wheels[@]} -eq 1 && -f ${repaired_wheels[0]} ]] || die "expected one repaired smsd-cuda wheel"
wheel="${repaired_wheels[0]}"
# The wheel replaces any earlier build; its dependencies come after, leaving the smsd just built in place.
python -m pip install --no-deps --force-reinstall "$wheel" >> "$RESULTS/build-smsd-cuda.log" 2>&1
python -m pip install "$wheel" >> "$RESULTS/build-smsd-cuda.log" 2>&1

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
log "build: $(basename "$wheel") installed; $(python -c 'import smsd; print(smsd.gpu_device_info())')"
