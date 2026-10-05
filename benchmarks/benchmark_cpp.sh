#!/usr/bin/env bash
# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
# CPU-only native MCS benchmark; passes arguments through to benchmark_cpp.
set -euo pipefail
script_dir="$(cd "$(dirname "$0")" && pwd)"
project_dir="$(cd "$script_dir/.." && pwd)"
build_dir="${SMSD_BENCHMARK_BUILD_DIR:-$project_dir/build/local-benchmarks/cpp}"
mkdir -p "$build_dir"
"${CXX:-c++}" -std=c++17 -O3 -I "$project_dir/cpp/include" \
    "$script_dir/benchmark_cpp.cpp" -o "$build_dir/benchmark_cpp"
"$build_dir/benchmark_cpp" "$@"
