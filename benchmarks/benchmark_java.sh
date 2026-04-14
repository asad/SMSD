#!/bin/bash
# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
# =============================================================================
# Benchmark SMSD Java CLI on the same molecule pairs as benchmark_rdkit.py
#
# Usage:
#   cd SMSD && mvn package -DskipTests
#   bash benchmarks/benchmark_java.sh
#
# Requires: Java 11+, built SMSD jar with dependencies
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
JAR="$PROJECT_DIR/target/smsd-6.0.0-jar-with-dependencies.jar"

if [ ! -f "$JAR" ]; then
  echo "ERROR: JAR not found at $JAR"
  echo "Build with: cd $PROJECT_DIR && mvn package -DskipTests"
  exit 1
fi

NUM_RUNS=5
TIMEOUT_MS=10000

# Molecule pairs: smi1|smi2|name
PAIRS=(
  "C|CC|methane-ethane"
  "O|CO|water-methanol"
  "c1ccccc1|c1ccc(O)cc1|benzene-phenol"
  "c1ccccc1|Cc1ccccc1|benzene-toluene"
  "CC(=O)Oc1ccccc1C(=O)O|CC(=O)Nc1ccc(O)cc1|aspirin-acetaminophen"
  "Cn1cnc2c1c(=O)n(C)c(=O)n2C|Cn1cnc2c1c(=O)[nH]c(=O)n2C|caffeine-theophylline"
  "CC(C)Cc1ccc(CC(C)C(O)=O)cc1|COc1ccc2cc(CC(C)C(O)=O)ccc2c1|ibuprofen-naproxen"
  "CN1CCC23C4C1CC5=C2C(=C(C=C5)O)OC3C(C=C4)O|COC1=CC=C2C3CC4=CC(=C(C=C4C3CC5=CC1=C2O5)OC)O|morphine-codeine"
  "C1=NC(=C2C(=N1)N(C=N2)C3C(C(C(O3)COP(=O)(O)OP(=O)(O)OP(=O)(O)O)O)O)N|C1=NC2=C(N1C3C(C(C(O3)COP(=O)(O)OP(=O)(O)OP(=O)(O)O)O)O)NC(=NC2=O)N|ATP-GTP"
  "CC1=C2C(C(=O)C3(C(CC4C(C3C(C(C2(C)C)(CC1OC(=O)C(C(C5=CC=CC=C5)NC(=O)C6=CC=CC=C6)O)O)OC(=O)C7=CC=CC=C7)(CO4)OC(=O)C)O)C)OC(=O)C|CC1=C2C(C(=O)C3(C(CC4C(C3C(C(C2(C)C)(CC1OC(=O)C(C(C5=CC=CC=C5)NC(=O)OC(C)(C)C)O)O)OC(=O)C6=CC=CC=C6)(CO4)OC(=O)C)O)C)O|paclitaxel-docetaxel"
  "C1C2CC3CC1CC(C2)C3|C1C2CC3CC1CC(C2)C3|adamantane-self"
  "c1cc2ccc3ccc4ccc5ccc6ccc1c7c2c3c4c5c67|c1cc2ccc3ccc4ccc5ccc6ccc1c7c2c3c4c5c67|coronene-self"
  "CCC1C(C(C(C(=O)C(CC(C(C(C(C(C(=O)O1)C)OC2CC(C(C(O2)C)O)(C)OC)C)OC3C(C(CC(O3)C)N(C)C)O)(C)O)C)C)O)(C)O|CCC1C(C(C(N(CC(CC(C(C(C(C(C(=O)O1)C)OC2CC(C(C(O2)C)O)(C)OC)C)OC3C(C(CC(O3)C)N(C)C)O)(C)O)C)C)C)O)(C)O|erythromycin-azithromycin"
  "OCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCO|OCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCO|PEG12-PEG16"
)

printf "%-30s  %10s  %10s  %4s\n" "Pair" "Best(ms)" "Median(ms)" "MCS"
printf "%s\n" "--------------------------------------------------------------"

TOTAL=0

for entry in "${PAIRS[@]}"; do
  IFS='|' read -r smi1 smi2 name <<< "$entry"

  # Warmup run (discard)
  java -jar "$JAR" --Q SMI --q "$smi1" --T SMI --t "$smi2" \
       --mode mcs --timeout "$TIMEOUT_MS" --json - 2>/dev/null >/dev/null || true

  times=()
  mcs_size=0

  for ((i = 0; i < NUM_RUNS; i++)); do
    start_ns=$(python3 -c "import time; print(int(time.perf_counter_ns()))")

    output=$(java -jar "$JAR" --Q SMI --q "$smi1" --T SMI --t "$smi2" \
                  --mode mcs --timeout "$TIMEOUT_MS" --json - 2>/dev/null) || true

    end_ns=$(python3 -c "import time; print(int(time.perf_counter_ns()))")
    dt_ms=$(python3 -c "print(f'{($end_ns - $start_ns) / 1e6:.3f}')")
    times+=("$dt_ms")

    # Extract MCS size from JSON output
    if [ -n "$output" ]; then
      mcs_size=$(echo "$output" | python3 -c "
import sys, json
try:
    d = json.load(sys.stdin)
    print(d.get('mcsSize', d.get('mappingSize', 0)))
except:
    print(0)
" 2>/dev/null) || mcs_size=0
    fi
  done

  # Compute best and median
  best=$(printf '%s\n' "${times[@]}" | sort -g | head -1)
  median=$(printf '%s\n' "${times[@]}" | sort -g | sed -n "$((NUM_RUNS / 2 + 1))p")
  TOTAL=$(python3 -c "print(f'{$TOTAL + $best:.3f}')")

  printf "%-30s  %10s  %10s  %4s\n" "$name" "$best" "$median" "$mcs_size"
done

printf "%s\n" "--------------------------------------------------------------"
printf "%-30s  %10s ms\n" "TOTAL" "$TOTAL"
