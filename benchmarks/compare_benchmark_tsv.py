#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
# Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
# See the NOTICE file for attribution, trademark, and algorithm IP terms.
"""
Compare two benchmark_python.py TSV outputs and summarise before/after changes.

Usage:
    python3 benchmarks/compare_benchmark_tsv.py before.tsv after.tsv
"""

from __future__ import annotations

import argparse
import csv
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List


@dataclass
class MCSRow:
    pair: str
    category: str
    smsd_best_us: float
    smsd_median_us: float
    smsd_mcs: int
    rdkit_best_us: float
    rdkit_median_us: float
    rdkit_mcs: int
    rdkit_timed_out: bool
    timing_verified: bool = False
    contract: str = ""


@dataclass
class SubRow:
    pair: str
    category: str
    smsd_median_us: float
    smsd_hit: bool
    rdkit_median_us: float
    rdkit_hit: bool


def parse_bool(value: str) -> bool:
    return value.strip().lower() == "true"


def parse_sections(path: Path) -> tuple[Dict[str, MCSRow], Dict[str, SubRow]]:
    mcs_rows: Dict[str, MCSRow] = {}
    sub_rows: Dict[str, SubRow] = {}
    metadata_path = path.with_suffix(".json")
    metadata = json.loads(metadata_path.read_text()) if metadata_path.exists() else {}
    details = {r["pair"]: r for r in metadata.get("results", []) if "policy" in r}

    current_section = None
    header: List[str] | None = None

    for raw_line in path.read_text().splitlines():
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith("#"):
            header = None
            if "MCS Benchmark" in line:
                current_section = "mcs"
            elif "Substructure Benchmark" in line:
                current_section = "sub"
            else:
                current_section = None
            continue
        if current_section is None:
            continue
        if header is None:
            header = line.split("\t")
            continue

        row = dict(zip(header, next(csv.reader([line], delimiter="\t"))))
        if current_section == "mcs":
            item = MCSRow(
                pair=row["pair"],
                category=row["category"],
                smsd_best_us=float(row["smsd_best_us"]),
                smsd_median_us=float(row["smsd_median_us"]),
                smsd_mcs=int(row["smsd_mcs"]),
                rdkit_best_us=float(row["rdkit_best_us"]),
                rdkit_median_us=float(row["rdkit_median_us"]),
                rdkit_mcs=int(row["rdkit_mcs"]),
                rdkit_timed_out=parse_bool(row["rdkit_timeout"]),
            )
            detail = details.get(item.pair, {})
            item.timing_verified = (detail.get("smsd_valid") is True
                and detail.get("smsd_budget_reached") == 0
                and detail.get("smsd_near_budget", 0) == 0
                and len(detail.get("smsd_atoms", [])) == 1)
            item.contract = json.dumps({key: metadata.get(key) for key in
                ("timeout_sec", "warmup", "iterations")}, sort_keys=True) + str(detail.get("policy"))
            mcs_rows[item.pair] = item
        elif current_section == "sub":
            item = SubRow(
                pair=row["pair"],
                category=row["category"],
                smsd_median_us=float(row["smsd_median_us"]),
                smsd_hit=parse_bool(row["smsd_hit"]),
                rdkit_median_us=float(row["rdkit_median_us"]),
                rdkit_hit=parse_bool(row["rdkit_hit"]),
            )
            sub_rows[item.pair] = item

    return mcs_rows, sub_rows


def fmt_us(value: float) -> str:
    if value >= 1_000_000:
        return f"{value / 1_000_000:.2f}s"
    if value >= 1_000:
        return f"{value / 1_000:.2f}ms"
    return f"{value:.0f}us"


def timing_bucket(before: float, after: float, threshold: float) -> str:
    if before <= 0 or after <= 0:
        return "unknown"
    if after <= before * (1.0 - threshold):
        return "improved"
    if after >= before * (1.0 + threshold):
        return "regressed"
    return "unchanged"


def total_timing_summary(before_total: float, after_total: float, threshold: float) -> str:
    if before_total <= 0 or after_total <= 0:
        return ""
    bucket = timing_bucket(before_total, after_total, threshold)
    if bucket == "unchanged":
        return f"{fmt_us(before_total)} -> {fmt_us(after_total)} (unchanged)"
    if bucket == "improved":
        return f"{fmt_us(before_total)} -> {fmt_us(after_total)} ({before_total / after_total:.2f}x faster)"
    return f"{fmt_us(before_total)} -> {fmt_us(after_total)} ({after_total / before_total:.2f}x slower)"


def compare_mcs(before: Dict[str, MCSRow], after: Dict[str, MCSRow], threshold: float) -> None:
    common_pairs = sorted(set(before) & set(after))
    if not common_pairs:
        print("No common MCS rows found.")
        return
    all_pairs = common_pairs
    common_pairs = [p for p in common_pairs if before[p].timing_verified and after[p].timing_verified
                    and before[p].contract == after[p].contract and before[p].smsd_mcs == after[p].smsd_mcs]
    print(f"MCS atom-count changes (quality not ranked): larger={sum(after[p].smsd_mcs > before[p].smsd_mcs for p in all_pairs)} "
          f"smaller={sum(after[p].smsd_mcs < before[p].smsd_mcs for p in all_pairs)}")
    if not common_pairs:
        print("No verified, equal-size MCS timing rows with matching protocol JSON.")
        return

    improved = 0
    regressed = 0
    unchanged = 0
    size_up = 0
    size_down = 0
    size_same = 0
    deltas = []

    before_total = 0.0
    after_total = 0.0

    for pair in common_pairs:
        b = before[pair]
        a = after[pair]
        before_total += b.smsd_median_us
        after_total += a.smsd_median_us

        bucket = timing_bucket(b.smsd_median_us, a.smsd_median_us, threshold)
        if bucket == "improved":
            improved += 1
        elif bucket == "regressed":
            regressed += 1
        else:
            unchanged += 1

        if a.smsd_mcs > b.smsd_mcs:
            size_up += 1
        elif a.smsd_mcs < b.smsd_mcs:
            size_down += 1
        else:
            size_same += 1

        ratio = (b.smsd_median_us / a.smsd_median_us) if a.smsd_median_us > 0 else 0.0
        deltas.append((pair, b.category, b.smsd_median_us, a.smsd_median_us, b.smsd_mcs, a.smsd_mcs, ratio))

    print("MCS")
    print(f"  Pairs compared: {len(common_pairs)}")
    print(f"  SMSD timing: improved={improved} regressed={regressed} unchanged={unchanged} (threshold {threshold * 100:.1f}%)")
    print(f"  Timing eligibility: {len(common_pairs)}/{len(all_pairs)}; validated equal-size mappings")
    summary = total_timing_summary(before_total, after_total, threshold)
    if summary:
        print(f"  Total SMSD median time: {summary}")

    deltas.sort(key=lambda item: item[6], reverse=True)
    print("  Top timing improvements:")
    shown = 0
    for pair, category, before_us, after_us, before_mcs, after_mcs, ratio in deltas:
        if timing_bucket(before_us, after_us, threshold) != "improved":
            continue
        print(f"    {pair} [{category}]: {fmt_us(before_us)} -> {fmt_us(after_us)} ({ratio:.2f}x), MCS {before_mcs} -> {after_mcs}")
        shown += 1
        if shown == 5:
            break
    if shown == 0:
        print("    none")

    deltas.sort(key=lambda item: item[6])
    print("  Top timing regressions:")
    shown = 0
    for pair, category, before_us, after_us, before_mcs, after_mcs, ratio in deltas:
        if timing_bucket(before_us, after_us, threshold) != "regressed":
            continue
        slowdown = (after_us / before_us) if before_us > 0 else 0.0
        print(f"    {pair} [{category}]: {fmt_us(before_us)} -> {fmt_us(after_us)} ({slowdown:.2f}x slower), MCS {before_mcs} -> {after_mcs}")
        shown += 1
        if shown == 5:
            break
    if shown == 0:
        print("    none")


def compare_sub(before: Dict[str, SubRow], after: Dict[str, SubRow], threshold: float) -> None:
    common_pairs = sorted(set(before) & set(after))
    if not common_pairs:
        print("Substructure")
        print("  No common substructure rows found.")
        return

    improved = 0
    regressed = 0
    unchanged = 0
    hit_changes = 0

    before_total = 0.0
    after_total = 0.0

    for pair in common_pairs:
        b = before[pair]
        a = after[pair]
        before_total += b.smsd_median_us
        after_total += a.smsd_median_us

        bucket = timing_bucket(b.smsd_median_us, a.smsd_median_us, threshold)
        if bucket == "improved":
            improved += 1
        elif bucket == "regressed":
            regressed += 1
        else:
            unchanged += 1

        if a.smsd_hit != b.smsd_hit:
            hit_changes += 1

    print("Substructure")
    print(f"  Pairs compared: {len(common_pairs)}")
    print(f"  SMSD timing: improved={improved} regressed={regressed} unchanged={unchanged} (threshold {threshold * 100:.1f}%)")
    print(f"  SMSD hit changes: {hit_changes}")
    summary = total_timing_summary(before_total, after_total, threshold)
    if summary:
        print(f"  Total SMSD median time: {summary}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Compare two benchmark_python.py TSV outputs")
    parser.add_argument("before", type=Path, help="Baseline TSV")
    parser.add_argument("after", type=Path, help="Candidate TSV")
    parser.add_argument("--threshold", type=float, default=0.05,
                        help="Relative change threshold for timing regression/improvement (default: 0.05)")
    args = parser.parse_args()

    before_mcs, before_sub = parse_sections(args.before)
    after_mcs, after_sub = parse_sections(args.after)

    print(f"Before: {args.before}")
    print(f"After:  {args.after}")
    print()
    compare_mcs(before_mcs, after_mcs, args.threshold)
    print()
    compare_sub(before_sub, after_sub, args.threshold)


if __name__ == "__main__":
    main()
