#!/usr/bin/env bash
# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
#
# Builds the smsd-metal wheel on an Apple silicon Mac, checks that it links only system libraries,
# installs it beside an smsd wheel of the same version in a fresh environment and runs the package
# tests there. Usage: build-and-test.sh SMSD_WHEEL. Every setting can be overridden from the
# environment.
set -euo pipefail

SCRIPTS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PACKAGE_DIR="$(cd "$SCRIPTS_DIR/.." && pwd)"
WORK="${SMSD_METAL_WORK:-$HOME/smsd-metal-work}"
PYTHON="${SMSD_METAL_PYTHON:-python3.14}"
export MACOSX_DEPLOYMENT_TARGET="${MACOSX_DEPLOYMENT_TARGET:-14.0}"
smsd_wheel="$(realpath "${1:?Usage: build-and-test.sh SMSD_WHEEL}")"

[[ "$(uname -s)/$(uname -m)" == Darwin/arm64 ]] || { echo "smsd-metal is built on Apple silicon Macs" >&2; exit 1; }
rm -rf "$WORK/venv" "$WORK/dist" "$WORK/run"
mkdir -p "$WORK/dist" "$WORK/run"
"$PYTHON" -m venv "$WORK/venv"
venv="$WORK/venv/bin/python"
# The same pinned tools as the smsd-cuda build, since the wheel is built without build isolation.
"$venv" -m pip install -q --upgrade pip==26.2.2
"$venv" -m pip install -q scikit-build-core==1.1.1 pytest==9.1.1 pytest-timeout==2.4.0 \
    cmake==4.4.3 ninja==1.13.2 delocate==0.13.0

"$venv" -m pip wheel "$PACKAGE_DIR" --no-deps --no-build-isolation -w "$WORK/dist"
wheels=("$WORK"/dist/smsd_metal-*-py3-none-macosx_*_arm64.whl)
[[ ${#wheels[@]} -eq 1 && -f ${wheels[0]} ]] || { echo "expected one smsd-metal wheel" >&2; exit 1; }
wheel="${wheels[0]}"
# The library links system frameworks and libraries only, so there is nothing to bundle.
"$WORK/venv/bin/delocate-listdeps" --all "$wheel"
if [[ -n "$("$WORK/venv/bin/delocate-listdeps" "$wheel")" ]]; then
    echo "the wheel links a library outside the system" >&2
    exit 1
fi

"$venv" -m pip install -q --no-deps "$smsd_wheel" "$wheel"
(cd "$WORK/run" && "$venv" -m pytest "$PACKAGE_DIR/tests" -q -rs -p no:cacheprovider)
echo "built and tested $wheel"
