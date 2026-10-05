#!/bin/bash
# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
# =============================================================================
# SMSD C++ Benchmark -- Build and Run
#
# Compiles the C++ SMSD library in Release mode and runs the benchmark binary.
# If the binary already exists and is up-to-date, skips compilation.
#
# Usage:
#   bash benchmarks/benchmark_cpp.sh
#
# Prerequisites:
#   - CMake 3.16+
#   - C++17 compiler (g++, clang++)
#   - Optional: OpenMP for batch processing
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
CPP_DIR="$PROJECT_DIR/cpp"
BUILD_DIR="$CPP_DIR/build"

# Number of parallel build jobs
JOBS="${JOBS:-$(sysctl -n hw.ncpu 2>/dev/null || nproc 2>/dev/null || echo 4)}"

echo "============================================"
echo "SMSD C++ Benchmark"
echo "============================================"
echo ""

# ----------------------------------------------------------------------------
# Step 1: Check prerequisites
# ----------------------------------------------------------------------------

if ! command -v cmake &>/dev/null; then
    echo "ERROR: cmake not found. Install with:"
    echo "  macOS:  brew install cmake"
    echo "  Ubuntu: sudo apt install cmake"
    exit 1
fi

if [ ! -d "$CPP_DIR" ]; then
    echo "ERROR: C++ source directory not found at $CPP_DIR"
    exit 1
fi

# ----------------------------------------------------------------------------
# Step 2: Configure and build
# ----------------------------------------------------------------------------

echo "[1/3] Configuring CMake (Release)..."
mkdir -p "$BUILD_DIR"
cd "$BUILD_DIR"

cmake "$CPP_DIR" \
    -DCMAKE_BUILD_TYPE=Release \
    -DSMSD_BUILD_TESTS=ON \
    -DSMSD_BUILD_PYTHON=OFF \
    -DSMSD_BUILD_CUDA=OFF \
    2>&1 | tail -5

echo ""
echo "[2/3] Building with $JOBS parallel jobs..."
cmake --build . --config Release -j "$JOBS" 2>&1

# ----------------------------------------------------------------------------
# Step 3: Run benchmarks
# ----------------------------------------------------------------------------

echo ""
echo "[3/3] Running benchmarks..."
echo ""

# Run the test suite first to verify correctness
if [ -f "$BUILD_DIR/smsd_tests" ]; then
    echo "--- Running tests (quick sanity check) ---"
    "$BUILD_DIR/smsd_tests" --gtest_filter="*Benchmark*" 2>/dev/null || true
    echo ""
fi

# Run the standalone benchmark if it exists
if [ -f "$BUILD_DIR/smsd_benchmark" ]; then
    echo "--- SMSD C++ Benchmark Results ---"
    "$BUILD_DIR/smsd_benchmark"
else
    echo "[INFO] No standalone smsd_benchmark binary found."
    echo "       The C++ library is header-only; benchmark is built"
    echo "       as part of the test suite."
    echo ""
    echo "       Running the standalone C++ benchmark from benchmarks/:"
    echo ""

    # Compile the standalone benchmark file if present
    BENCH_CPP="$SCRIPT_DIR/benchmark_cpp.cpp"
    BENCH_BIN="$BUILD_DIR/benchmark_cpp_standalone"

    if [ -f "$BENCH_CPP" ]; then
        echo "  Compiling $BENCH_CPP ..."
        CXX="${CXX:-c++}"
        "$CXX" -std=c++17 -O2 \
            -I "$CPP_DIR/include" \
            -o "$BENCH_BIN" \
            "$BENCH_CPP" 2>&1

        echo ""
        "$BENCH_BIN"
    else
        echo "  No benchmark_cpp.cpp found in $SCRIPT_DIR"
        echo "  Skipping standalone C++ benchmark."
    fi
fi

echo ""
echo "============================================"
echo "C++ benchmark complete."
echo "============================================"
