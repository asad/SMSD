#!/usr/bin/env bash
# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
#
# Runs C++ and installed-wheel Python tests with GPU dispatch enabled. Individual screens still
# apply their operation-specific crossovers and CPU fallback; the resident callback tests, the
# device and CPU parity tests and the explicit packed CUDA backend tests use CUDA.
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"
unset SMSD_DISABLE_GPU
status=0
source_manifest="${SMSD_CUDA_SOURCE_MANIFEST:?Set SMSD_CUDA_SOURCE_MANIFEST to the canonical source manifest}"
source_manifest="$(realpath "$source_manifest")"
release_dir="${SMSD_CUDA_RELEASE_DIR:-$(dirname "$source_manifest")}"
release_dir="$(realpath "$release_dir")"
[[ "$source_manifest" == "$release_dir/SOURCE_MANIFEST.json" ]] || {
    printf '%s\n' 'use the canonical SOURCE_MANIFEST.json in the release directory' >&2; exit 1;
}
results="$(cd "$RESULTS" && pwd)"
case "$results/" in "$release_dir/"*) printf '%s\n' 'keep private build/test results outside the release directory' >&2; exit 1 ;; esac
use_cuda
use_python
evidence="$(mktemp -d "$results/validation.XXXXXX")"
build_note="$results/build-inputs.json"
sdk_source="$(python -c 'import json,sys; print(json.load(open(sys.argv[1]))["sdk_source"])' "$build_note")"
wheels=("$RESULTS"/repaired-wheels/smsd_cuda-*.whl)
[[ ${#wheels[@]} -eq 1 && -f ${wheels[0]} ]] || die "expected one repaired CUDA wheel"
wheel="$(realpath "${wheels[0]}")"

log "test: C++ GPU tests"
ctest_started="$(date +%s%N)"
if (cd "$WORK/cpp-cuda" && ctest --output-on-failure --no-tests=error -j 2 -R "$GPU_TEST_REGEX" \
        --output-junit "$evidence/ctest-gpu.xml") > "$evidence/ctest-gpu.log" 2>&1; then
    ctest_status=0
    log "test: $(grep 'tests passed' "$evidence/ctest-gpu.log")"
else
    ctest_status=$?
    status=1
    log "test: C++ GPU tests failed; see $evidence/ctest-gpu.log"
fi
ctest_finished="$(date +%s%N)"
log "test: Python tests with smsd-cuda"
if (cd "$WORK" && PYTHONPATH="$ROOT/scripts${PYTHONPATH:+:$PYTHONPATH}" \
        SMSD_CUDA_SOURCE_MANIFEST="$source_manifest" SMSD_CUDA_TEST_WHEEL="$wheel" \
        SMSD_CUDA_TEST_IDENTITY="$evidence/installed-cuda-identity.json" SMSD_CUDA_BUILD_NOTE="$build_note" \
        SMSD_CUDA_NATIVE_XML="$evidence/ctest-gpu.xml" SMSD_CUDA_CTEST_STATUS="$ctest_status" \
        SMSD_CUDA_CTEST_STARTED_NS="$ctest_started" SMSD_CUDA_CTEST_FINISHED_NS="$ctest_finished" \
        SMSD_TEST_PACKED_CUDA=1 \
        python -m pytest "$sdk_source/python/tests" "$PACKAGE_DIR/tests" -q -rs -p no:cacheprovider \
        -p cuda_wheel_provenance --junitxml="$evidence/pytest-gpu.xml" -c "$sdk_source/pyproject.toml" --rootdir "$sdk_source") \
        > "$evidence/pytest-gpu.log" 2>&1; then
    log "test: $(tail -1 "$evidence/pytest-gpu.log")"
else
    status=1
    log "test: Python tests failed; see $evidence/pytest-gpu.log"
fi
if [[ $status -eq 0 ]]; then
    python "$ROOT/scripts/cuda_wheel_provenance.py" --release-dir "$release_dir" --wheel "$wheel" \
        --identity "$evidence/installed-cuda-identity.json" --native-xml "$evidence/ctest-gpu.xml" \
        --python-xml "$evidence/pytest-gpu.xml" --evidence-summary "$evidence/cuda-validation-summary.json" \
        --checkout "$ROOT" || status=1
fi
exit "$status"
