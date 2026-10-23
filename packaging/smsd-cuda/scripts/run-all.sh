#!/usr/bin/env bash
# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
#
# Runs setup, build, test, identity and timing in turn and packs the results folder into one
# archive in the home folder. A failed setup or build stops the run; a failed check is recorded
# and the remaining checks still run. The whole run takes under an hour on a 4-CPU machine with
# an NVIDIA T4; stop the machine afterwards, as it is billed while it runs.
#
#   nohup bash packaging/smsd-cuda/scripts/run-all.sh > ~/run-all.out 2>&1 &
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"

started=$(date +%s)
failed=()
: > "$RESULTS/run.log"
log "run-all: $ROOT, results in $RESULTS"
for step in setup build test identity timing; do
    step_started=$(date +%s)
    if timeout 30m bash "$SCRIPTS_DIR/$step.sh"; then
        log "run-all: $step passed in $(( ($(date +%s) - step_started) / 60 )) min"
    else
        failed+=("$step")
        log "run-all: $step failed"
        if [[ $step == setup || $step == build ]]; then break; fi
    fi
done

minutes=$(( ($(date +%s) - started) / 60 ))
archive="$HOME/smsd-cuda-results-$(date -u +%Y%m%d-%H%M).tar.gz"
log "run-all: finished in $minutes min; failed: ${failed[*]:-none}"
tar -czf "$archive" -C "$(dirname "$RESULTS")" "$(basename "$RESULTS")"
echo "Results: $archive"
echo "Copy them back, for example: scp <user>@<host>:$archive ."
echo "Then stop or terminate the machine."
[[ ${#failed[@]} -eq 0 ]]
