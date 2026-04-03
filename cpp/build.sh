#!/usr/bin/env bash
set -euo pipefail

ROOT="/Users/asad/tool/SMSD"
CPP_DIR="$ROOT/cpp"
BUILD_DIR="$CPP_DIR/build"

MODE="${1:-cpu}"          # cpu | gpu
BACKEND="${2:-auto}"      # auto | metal | cuda

OS="$(uname -s)"
ARCH="$(uname -m)"

METAL_FLAG="OFF"
CUDA_FLAG="OFF"
CTEST_EXCLUDE=""
JOBS="$(sysctl -n hw.logicalcpu 2>/dev/null || nproc 2>/dev/null || echo 8)"

# Default to CPU-safe runtime
export SMSD_DISABLE_GPU=1
export SMSD_FORCE_CPU=1
unset SMSD_FORCE_METAL || true

case "$MODE" in
  cpu)
    METAL_FLAG="OFF"
    CUDA_FLAG="OFF"
    CTEST_EXCLUDE="-E smsd_batch_gpu_tests"
    ;;

  gpu)
    export SMSD_DISABLE_GPU=0
    export SMSD_FORCE_CPU=0

    if [[ "$BACKEND" == "auto" ]]; then
      if [[ "$OS" == "Darwin" ]]; then
        BACKEND="metal"
      else
        BACKEND="cuda"
      fi
    fi

    case "$BACKEND" in
      metal)
        METAL_FLAG="ON"
        CUDA_FLAG="OFF"
        export SMSD_FORCE_METAL=1
        ;;
      cuda)
        METAL_FLAG="OFF"
        CUDA_FLAG="ON"
        ;;
      *)
        echo "Unknown backend: $BACKEND"
        echo "Use: auto | metal | cuda"
        exit 1
        ;;
    esac
    ;;

  *)
    echo "Unknown mode: $MODE"
    echo "Usage:"
    echo "  $0 cpu"
    echo "  $0 gpu"
    echo "  $0 gpu metal"
    echo "  $0 gpu cuda"
    exit 1
    ;;
esac

echo "== SMSD build =="
echo "Mode:     $MODE"
echo "Backend:  $BACKEND"
echo "OS/Arch:  $OS / $ARCH"
echo "Metal:    $METAL_FLAG"
echo "CUDA:     $CUDA_FLAG"
echo "Jobs:     $JOBS"
echo

rm -rf "$BUILD_DIR"

cmake -S "$CPP_DIR" -B "$BUILD_DIR" \
  -DCMAKE_BUILD_TYPE=Release \
  -DSMSD_BUILD_METAL="$METAL_FLAG" \
  -DSMSD_BUILD_CUDA="$CUDA_FLAG"

cmake --build "$BUILD_DIR" --config Release -j "$JOBS"

if [[ -n "$CTEST_EXCLUDE" ]]; then
  ctest --test-dir "$BUILD_DIR" --output-on-failure -V $CTEST_EXCLUDE
else
  ctest --test-dir "$BUILD_DIR" --output-on-failure -V
fi

cd "$ROOT"
mvn clean test

cd "$ROOT/benchmarks"
python3 benchmark_python.py

echo "ALL DONE"
