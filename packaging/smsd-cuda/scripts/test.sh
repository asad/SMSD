#!/usr/bin/env bash
# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
#
# Runs the C++ GPU tests and the Python tests with the GPU enabled (SMSD_DISABLE_GPU unset), with
# smsd-cuda installed, so that every TargetCorpus screen of the Python tests runs on the GPU.
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"
use_cuda
use_python
unset SMSD_DISABLE_GPU
status=0

log "test: C++ GPU tests"
if (cd "$WORK/cpp-cuda" && ctest --output-on-failure -j 2 -R "$GPU_TEST_REGEX") > "$RESULTS/ctest-gpu.log" 2>&1; then
    log "test: $(grep 'tests passed' "$RESULTS/ctest-gpu.log")"
else
    status=1
    log "test: C++ GPU tests failed; see ctest-gpu.log"
fi
"$WORK/cpp-cuda/smsd_cuda_host_link_tests" | tee -a "$RESULTS/run.log" || status=1

log "test: Python tests with smsd-cuda"
python -c 'import smsd, smsd_cuda; print(smsd.gpu_is_available(), smsd.gpu_device_info(), smsd_cuda.unavailable_reason())' \
    | tee -a "$RESULTS/run.log" || status=1
if (cd "$WORK" && python -m pytest "$ROOT/python/tests" "$PACKAGE_DIR/tests" -q -rs -p no:cacheprovider) \
        > "$RESULTS/pytest-gpu.log" 2>&1; then
    log "test: $(tail -1 "$RESULTS/pytest-gpu.log")"
else
    status=1
    log "test: Python tests failed; see pytest-gpu.log"
fi
exit "$status"
