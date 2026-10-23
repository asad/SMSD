#!/usr/bin/env bash
# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
#
# Times CPU and GPU screens over the benchmark molecules repeated to each size: the median of
# repeated screens, on all CPUs and on one, and on the GPU through smsd-cuda; then the CUDA screen
# of the C++ library. Writes timing.tsv, timing-cpp.tsv and timing-summary.txt.
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"
use_cuda
use_python
require_corpora
unset SMSD_DISABLE_GPU
status=0

log "timing: smsd-cuda against the CPU, sizes $SIZES"
python "$PACKAGE_DIR/tools/screen_compare.py" --mode timing --sizes "$SIZES" --out "$RESULTS" "${CORPORA[@]}" \
    > "$RESULTS/timing-python.log" 2>&1 || { status=1; log "timing: smsd-cuda timing failed; see timing-python.log"; }

log "timing: C++ library CUDA screen against the CPU"
"$WORK/tools/screen_compare" --sizes "${SMSD_CUDA_CPP_SIZES:-10000,100000}" --repeats 3 \
    --out "$RESULTS/timing-cpp.tsv" "${CORPORA[@]}" > "$RESULTS/timing-cpp.log" 2>&1 \
    || { status=1; log "timing: C++ timing failed; see timing-cpp.log"; }

python - "$RESULTS" > "$RESULTS/timing-summary.txt" << 'EOF'
import csv, os, statistics, sys
from collections import defaultdict

folder = sys.argv[1]
for name, cpu, gpu in (("timing.tsv", "cpu_ms", "gpu_ms"), ("timing-cpp.tsv", "cpu_ms", "gpu_ms")):
    path = os.path.join(folder, name)
    if not os.path.exists(path):
        continue
    rows = defaultdict(list)
    with open(path, encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            rows[int(row["size"])].append(row)
    print(f"{name}: median over queries and thresholds, milliseconds")
    print(f"{'size':>9} {'cpu':>10} {'cpu 1 thread':>13} {'gpu':>10} {'cpu/gpu':>8}")
    for size in sorted(rows):
        cpu_ms = statistics.median(float(r[cpu]) for r in rows[size])
        gpu_ms = statistics.median(float(r[gpu]) for r in rows[size])
        one = [float(r["cpu_1_thread_ms"]) for r in rows[size] if "cpu_1_thread_ms" in r]
        one_ms = f"{statistics.median(one):13.3f}" if one else f"{'':>13}"
        print(f"{size:>9} {cpu_ms:10.3f} {one_ms} {gpu_ms:10.3f} {cpu_ms / gpu_ms:8.2f}")
    print()
EOF
tee -a "$RESULTS/run.log" < "$RESULTS/timing-summary.txt"
exit "$status"
