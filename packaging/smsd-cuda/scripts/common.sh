#!/usr/bin/env bash
# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
#
# Settings shared by the smsd-cuda build and check scripts, which source this file. Every
# setting can be overridden from the environment.
set -euo pipefail

SCRIPTS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PACKAGE_DIR="$(cd "$SCRIPTS_DIR/.." && pwd)"
ROOT="$(cd "$PACKAGE_DIR/../.." && pwd)"
WORK="${SMSD_CUDA_WORK:-$HOME/smsd-cuda-work}"
RESULTS="${SMSD_CUDA_RESULTS:-$HOME/smsd-cuda-results}"
VENV="$WORK/venv"
PYTHON_VERSION="${SMSD_CUDA_PYTHON:-3.14}"
TOOLKIT="${SMSD_CUDA_TOOLKIT:-12.8}"
SIZES="${SMSD_CUDA_SIZES:-10000,100000,1000000}"
JOBS="${SMSD_CUDA_JOBS:-$(nproc 2>/dev/null || echo 2)}"
CORPORA=("$ROOT/benchmarks/data/chembl_mcs_benchmark.smi" "$ROOT/benchmarks/diverse_molecules.txt")
GPU_TEST_TARGETS=(smsd_batch_gpu_tests smsd_cuda_host_link_tests smsd_native_gpu_domain_tests
                  smsd_gpu_domain_regression_tests smsd_screening_regression_tests smsd_screen_accelerator_tests)
GPU_TEST_REGEX='^smsd_(batch_gpu|cuda_host_link|native_gpu_domain|gpu_domain_regression|screening_regression|screen_accelerator)_tests$'
export WORK RESULTS

mkdir -p "$WORK" "$RESULTS"

log() { printf '[%s] %s\n' "$(date -u +%H:%M:%S)" "$*" | tee -a "$RESULTS/run.log"; }
die() { log "error: $*"; exit 1; }

# The nvcc of a CUDA 12 toolkit, preferring the version the package is released with.
find_nvcc() {
    local candidate
    local candidates=("/usr/local/cuda-$TOOLKIT/bin/nvcc" /usr/local/cuda/bin/nvcc)
    if command -v nvcc > /dev/null 2>&1; then candidates+=("$(command -v nvcc)"); fi
    for candidate in /usr/local/cuda-12*/bin/nvcc; do candidates+=("$candidate"); done
    for candidate in "${candidates[@]}"; do
        if [[ -x $candidate ]] && "$candidate" --version | grep -qE 'release 12\.'; then
            echo "$candidate"
            return 0
        fi
    done
    return 1
}

# Puts the CUDA toolkit on PATH for CMake and the builds.
use_cuda() {
    local nvcc
    nvcc="$(find_nvcc)" || die "no CUDA 12 toolkit was found; run setup.sh"
    CUDA_HOME="$(dirname "$(dirname "$nvcc")")"
    export CUDA_HOME CUDACXX="$nvcc"
    export PATH="$CUDA_HOME/bin:$PATH"
    export LD_LIBRARY_PATH="$CUDA_HOME/lib64${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
}

# Activates the Python environment setup.sh made.
use_python() {
    [[ -x $VENV/bin/python ]] || die "no Python environment at $VENV; run setup.sh"
    # shellcheck disable=SC1091
    source "$VENV/bin/activate"
}

require_corpora() {
    local corpus
    for corpus in "${CORPORA[@]}"; do
        [[ -s $corpus ]] || die "the benchmark molecules $corpus are missing; copy them into the checkout"
    done
}
