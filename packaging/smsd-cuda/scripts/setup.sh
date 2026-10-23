#!/usr/bin/env bash
# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
#
# Prepares an Ubuntu or Amazon Linux 2023 x86_64 machine with an NVIDIA GPU: checks the driver,
# installs the CUDA 12 toolkit when none is present, Python 3.14 and the build tools, adds swap
# when memory is short of the million-molecule screens, and records the environment.
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"

log "setup: $(. /etc/os-release && echo "$PRETTY_NAME"), $(nproc) CPUs, $(free -g | awk '/^Mem:/ {print $2}') GB memory"
if ! nvidia-smi > /dev/null 2>&1; then
    die "no working NVIDIA driver (nvidia-smi failed). Use an image that has it, such as a Deep Learning \
Base GPU image, or install it with: sudo apt-get install -y ubuntu-drivers-common && \
sudo ubuntu-drivers install --gpgpu; then reboot and run this again"
fi
nvidia-smi --query-gpu=name,driver_version,memory.total,compute_cap --format=csv,noheader | tee -a "$RESULTS/run.log"

if command -v apt-get > /dev/null 2>&1; then
    export DEBIAN_FRONTEND=noninteractive
    sudo apt-get update -q > "$RESULTS/apt.log" 2>&1
    sudo apt-get install -y -q build-essential cmake ninja-build git curl ca-certificates patchelf unzip \
        software-properties-common >> "$RESULTS/apt.log" 2>&1
else
    # Amazon Linux and other dnf systems; CMake and Ninja come from the Python environment below.
    sudo dnf install -y -q gcc-c++ make git patchelf unzip tar > "$RESULTS/dnf.log" 2>&1
fi

if ! find_nvcc > /dev/null; then
    log "setup: installing the CUDA $TOOLKIT toolkit from NVIDIA's repository"
    if command -v apt-get > /dev/null 2>&1; then
        distribution="ubuntu$(. /etc/os-release && echo "${VERSION_ID//./}")"
        keyring="$WORK/cuda-keyring.deb"
        curl -fsSL -o "$keyring" "https://developer.download.nvidia.com/compute/cuda/repos/$distribution/x86_64/cuda-keyring_1.1-1_all.deb"
        sudo dpkg -i "$keyring" >> "$RESULTS/apt.log" 2>&1
        sudo apt-get update -q >> "$RESULTS/apt.log" 2>&1
        sudo apt-get install -y -q "cuda-toolkit-${TOOLKIT//./-}" >> "$RESULTS/apt.log" 2>&1
    elif [[ $(. /etc/os-release && echo "$ID:$VERSION_ID") == amzn:2023 ]]; then
        repository="$WORK/cuda-amzn2023.repo"
        curl -fsSL -o "$repository" https://developer.download.nvidia.com/compute/cuda/repos/amzn2023/x86_64/cuda-amzn2023.repo
        sudo install -m 644 "$repository" /etc/yum.repos.d/cuda-amzn2023.repo
        sudo dnf install -y -q --disablerepo='cuda-rhel*' --exclude='nvidia-driver*' --exclude='cuda-drivers*' \
            "cuda-toolkit-${TOOLKIT//./-}" >> "$RESULTS/dnf.log" 2>&1
    else
        die "install the CUDA $TOOLKIT toolkit for this distribution, then run setup.sh again"
    fi
fi
use_cuda
log "setup: $("$CUDACXX" --version | grep release)"

if [[ ! -x $VENV/bin/python ]]; then
    if command -v "python$PYTHON_VERSION" > /dev/null 2>&1; then
        "python$PYTHON_VERSION" -m venv "$VENV"
    elif command -v apt-get > /dev/null 2>&1; then
        log "setup: installing Python $PYTHON_VERSION"
        sudo add-apt-repository -y ppa:deadsnakes/ppa >> "$RESULTS/apt.log" 2>&1
        sudo apt-get install -y -q "python$PYTHON_VERSION" "python$PYTHON_VERSION-venv" "python$PYTHON_VERSION-dev" \
            >> "$RESULTS/apt.log" 2>&1
        "python$PYTHON_VERSION" -m venv "$VENV"
    else
        # No packaged Python of this version: use a standalone build through uv.
        log "setup: installing Python $PYTHON_VERSION with uv"
        python3 -m pip install -q --user "uv<1" > "$RESULTS/uv.log" 2>&1
        "$HOME/.local/bin/uv" python install "$PYTHON_VERSION" >> "$RESULTS/uv.log" 2>&1
        "$HOME/.local/bin/uv" venv --seed --python "$PYTHON_VERSION" "$VENV" >> "$RESULTS/uv.log" 2>&1
    fi
fi
use_python
python -m pip install -q --upgrade pip > "$RESULTS/pip-setup.log" 2>&1
python -m pip install -q "scikit-build-core>=0.10,<2" "pybind11>=2.12,<4" "pytest>=7,<10" "pytest-timeout>=2,<3" \
    auditwheel "cmake>=3.18,<5" "ninja<2" >> "$RESULTS/pip-setup.log" 2>&1
# The tests that use NumPy or RDKit are skipped when these have no build for this Python.
python -m pip install -q numpy "rdkit<2027" >> "$RESULTS/pip-setup.log" 2>&1 || log "setup: NumPy or RDKit unavailable; their tests will be skipped"

# A million copies of the benchmark molecules take about 14 GB; swap covers a machine with 16 GB.
memory_gb=$(free -g | awk '/^Mem:/ {print $2}')
if (( memory_gb < 24 )) && [[ -z $(swapon --show --noheadings) ]]; then
    log "setup: adding 16 GB of swap"
    sudo fallocate -l 16G /swapfile-smsd && sudo chmod 600 /swapfile-smsd
    sudo mkswap /swapfile-smsd > /dev/null && sudo swapon /swapfile-smsd
fi

{
    echo "date: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
    echo "system: $(. /etc/os-release && echo "$PRETTY_NAME"), kernel $(uname -r)"
    echo "cpu: $(lscpu | awk -F: '/Model name/ {gsub(/^ +/, "", $2); print $2}'), $(nproc) CPUs"
    echo "memory: $(free -h | awk '/^Mem:/ {print $2}'), swap $(free -h | awk '/^Swap:/ {print $2}')"
    echo "gpu: $(nvidia-smi --query-gpu=name,driver_version,memory.total,compute_cap --format=csv,noheader)"
    echo "cuda: $("$CUDACXX" --version | grep release)"
    echo "compiler: $(gcc --version | head -1)"
    echo "cmake: $(cmake --version | head -1)"
    echo "python: $(python --version)"
    echo "source: $ROOT"
} > "$RESULTS/environment.txt"
cat "$RESULTS/environment.txt"
log "setup: done"
