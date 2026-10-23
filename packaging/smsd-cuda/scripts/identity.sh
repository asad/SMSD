#!/usr/bin/env bash
# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
#
# Checks that GPU and CPU screens give the same hits over the benchmark molecules repeated to
# 10,000, 100,000 and 1,000,000 targets: smsd-cuda through smsd.screen_targets, and the CUDA
# screen of the C++ library, whose 2 KB descriptor per target limits it to 100,000 here.
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"
use_cuda
use_python
require_corpora
unset SMSD_DISABLE_GPU
status=0

log "identity: smsd-cuda against the CPU, sizes $SIZES"
if python "$PACKAGE_DIR/tools/screen_compare.py" --mode identity --sizes "$SIZES" --out "$RESULTS" "${CORPORA[@]}" \
        > "$RESULTS/identity-python.log" 2>&1; then
    log "identity: $(tail -1 "$RESULTS/identity-python.log")"
else
    status=1
    log "identity: smsd-cuda check failed; see identity-python.log"
fi

log "identity: C++ library CUDA screen against the CPU"
if "$WORK/tools/screen_compare" --sizes "${SMSD_CUDA_CPP_SIZES:-10000,100000}" --repeats 1 \
        --out "$RESULTS/identity-cpp.tsv" "${CORPORA[@]}" > "$RESULTS/identity-cpp.log" 2>&1; then
    log "identity: $(tail -1 "$RESULTS/identity-cpp.log")"
else
    status=1
    log "identity: C++ check failed; see identity-cpp.log"
fi
exit "$status"
