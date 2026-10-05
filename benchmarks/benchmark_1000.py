#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
# Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
# See the NOTICE file for attribution, trademark, and algorithm IP terms.
"""
SMSD 1000-Molecule Benchmark: RDKit FindMCS vs SMSD Java
=========================================================

Reads diverse_molecules.txt (1000+ molecules), creates 1000 pairs
(500 random + 500 systematic scaffold-matched), and benchmarks both
RDKit FindMCS and SMSD Java on each pair.

Usage:
    pip install rdkit
    python benchmark_1000.py [--pairs N] [--rounds N] [--timeout N] [--smsd-jar PATH]

Output:
    benchmark_results.tsv — tab-separated results for every pair
    benchmark_summary.txt — aggregated statistics

Requirements:
    - RDKit (pip install rdkit)
    - Java 11+ on PATH (for SMSD comparison)
    - SMSD jar (built via: cd .. && mvn package -DskipTests)
"""

import argparse
import csv
import hashlib
import os
import random
import signal
import subprocess
import sys
import tempfile
import time
from collections import defaultdict
from pathlib import Path
from statistics import median, mean

# ---------------------------------------------------------------------------
# RDKit imports
# ---------------------------------------------------------------------------
try:
    from rdkit import Chem, RDLogger
    from rdkit.Chem import rdFMCS, AllChem, Descriptors
    from rdkit.Chem.Scaffolds import MurckoScaffold

    RDLogger.logger().setLevel(RDLogger.ERROR)
    HAS_RDKIT = True
except ImportError:
    HAS_RDKIT = False
    print("WARNING: RDKit not found. Only SMSD Java benchmarks will run.", file=sys.stderr)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
SCRIPT_DIR = Path(__file__).resolve().parent
MOLECULES_FILE = SCRIPT_DIR / "diverse_molecules.txt"
DEFAULT_SMSD_JAR = SCRIPT_DIR.parent / "target" / "smsd-6.4.0.jar"
OUTPUT_TSV = SCRIPT_DIR / "benchmark_results.tsv"
OUTPUT_SUMMARY = SCRIPT_DIR / "benchmark_summary.txt"

SECTION_RANGES = {
    "tiny":       (1, 50),
    "aromatic":   (51, 150),
    "druglike":   (151, 350),
    "complex":    (351, 550),
    "large":      (551, 700),
    "very_large": (701, 800),
    "tautomer":   (801, 900),
    "symmetric":  (901, 950),
    "edge":       (951, 1050),
}


# ---------------------------------------------------------------------------
# Molecule loading
# ---------------------------------------------------------------------------
def load_molecules(path: Path) -> list[tuple[str, str]]:
    """Load (smiles, name) from the benchmark file, skipping comments/blanks."""
    molecules = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t", 1)
            smi = parts[0].strip()
            name = parts[1].strip() if len(parts) > 1 else f"mol_{len(molecules)}"
            molecules.append((smi, name))
    return molecules


def validate_smiles(molecules: list[tuple[str, str]]) -> list[tuple[str, str, object]]:
    """Validate SMILES with RDKit, return (smiles, name, mol) for valid ones."""
    if not HAS_RDKIT:
        return [(s, n, None) for s, n in molecules]
    valid = []
    invalid_count = 0
    for smi, name in molecules:
        mol = Chem.MolFromSmiles(smi)
        if mol is not None:
            valid.append((smi, name, mol))
        else:
            invalid_count += 1
            print(f"  INVALID SMILES skipped: {name} -> {smi}", file=sys.stderr)
    if invalid_count:
        print(f"  {invalid_count} invalid SMILES skipped out of {len(molecules)}", file=sys.stderr)
    return valid


# ---------------------------------------------------------------------------
# Pair generation
# ---------------------------------------------------------------------------
def get_section(idx: int) -> str:
    """Return the section name for a 1-based molecule index."""
    for name, (lo, hi) in SECTION_RANGES.items():
        if lo <= idx <= hi:
            return name
    return "unknown"


def compute_scaffold(mol) -> str:
    """Return Murcko scaffold SMILES, or empty string on failure."""
    if not HAS_RDKIT or mol is None:
        return ""
    try:
        scaf = MurckoScaffold.GetScaffoldForMol(mol)
        return Chem.MolToSmiles(scaf)
    except Exception:
        return ""


def generate_pairs(
    mols: list[tuple[str, str, object]], n_random: int = 500, n_systematic: int = 500, seed: int = 42
) -> list[dict]:
    """
    Generate benchmark pairs:
      - n_random: uniformly random pairs across all molecules
      - n_systematic: pairs matched by Murcko scaffold (similar structures)
    Returns list of dicts with keys: idx_a, idx_b, smi_a, smi_b, name_a, name_b, pair_type
    """
    rng = random.Random(seed)
    n = len(mols)
    pairs = []

    # --- Random pairs ---
    seen = set()
    attempts = 0
    while len(pairs) < n_random and attempts < n_random * 20:
        i = rng.randint(0, n - 1)
        j = rng.randint(0, n - 1)
        if i == j:
            attempts += 1
            continue
        key = (min(i, j), max(i, j))
        if key in seen:
            attempts += 1
            continue
        seen.add(key)
        pairs.append({
            "idx_a": i, "idx_b": j,
            "smi_a": mols[i][0], "smi_b": mols[j][0],
            "name_a": mols[i][1], "name_b": mols[j][1],
            "pair_type": "random",
        })

    # --- Systematic scaffold-matched pairs ---
    if HAS_RDKIT:
        scaffold_groups = defaultdict(list)
        for idx, (smi, name, mol) in enumerate(mols):
            scaf = compute_scaffold(mol)
            if scaf:
                scaffold_groups[scaf].append(idx)

        # Collect groups with 2+ members
        multi_groups = [idxs for idxs in scaffold_groups.values() if len(idxs) >= 2]
        rng.shuffle(multi_groups)

        sys_count = 0
        for idxs in multi_groups:
            if sys_count >= n_systematic:
                break
            # Generate pairs within this scaffold group
            rng.shuffle(idxs)
            for k in range(0, len(idxs) - 1, 2):
                if sys_count >= n_systematic:
                    break
                i, j = idxs[k], idxs[k + 1]
                key = (min(i, j), max(i, j))
                if key in seen:
                    continue
                seen.add(key)
                pairs.append({
                    "idx_a": i, "idx_b": j,
                    "smi_a": mols[i][0], "smi_b": mols[j][0],
                    "name_a": mols[i][1], "name_b": mols[j][1],
                    "pair_type": "scaffold",
                })
                sys_count += 1

        # If not enough scaffold pairs, fill with size-matched pairs
        if sys_count < n_systematic:
            size_sorted = sorted(range(n), key=lambda x: len(mols[x][0]))
            for k in range(0, len(size_sorted) - 1, 2):
                if sys_count >= n_systematic:
                    break
                i, j = size_sorted[k], size_sorted[k + 1]
                key = (min(i, j), max(i, j))
                if key in seen:
                    continue
                seen.add(key)
                pairs.append({
                    "idx_a": i, "idx_b": j,
                    "smi_a": mols[i][0], "smi_b": mols[j][0],
                    "name_a": mols[i][1], "name_b": mols[j][1],
                    "pair_type": "size_matched",
                })
                sys_count += 1
    else:
        # Without RDKit, use size-matched pairs
        size_sorted = sorted(range(n), key=lambda x: len(mols[x][0]))
        sys_count = 0
        for k in range(0, len(size_sorted) - 1, 2):
            if sys_count >= n_systematic:
                break
            i, j = size_sorted[k], size_sorted[k + 1]
            key = (min(i, j), max(i, j))
            if key in seen:
                continue
            seen.add(key)
            pairs.append({
                "idx_a": i, "idx_b": j,
                "smi_a": mols[i][0], "smi_b": mols[j][0],
                "name_a": mols[i][1], "name_b": mols[j][1],
                "pair_type": "size_matched",
            })
            sys_count += 1

    return pairs


# ---------------------------------------------------------------------------
# RDKit MCS benchmark
# ---------------------------------------------------------------------------
def rdkit_mcs_single(smi_a: str, smi_b: str, timeout_sec: int = 10) -> tuple[float, int, bool]:
    """
    Run RDKit FindMCS on a single pair.
    Returns (elapsed_seconds, mcs_atom_count, completed_without_timeout).
    """
    if not HAS_RDKIT:
        return (-1.0, -1, False)
    mol_a = Chem.MolFromSmiles(smi_a)
    mol_b = Chem.MolFromSmiles(smi_b)
    if mol_a is None or mol_b is None:
        return (-1.0, -1, False)

    t0 = time.perf_counter()
    try:
        result = rdFMCS.FindMCS(
            [mol_a, mol_b],
            timeout=timeout_sec,
            matchValences=False,
            ringMatchesRingOnly=True,
            completeRingsOnly=False,
            bondCompare=rdFMCS.BondCompare.CompareOrder,
            atomCompare=rdFMCS.AtomCompare.CompareElements,
        )
        elapsed = time.perf_counter() - t0
        completed = not result.canceled
        mcs_size = result.numAtoms
        return (elapsed, mcs_size, completed)
    except Exception as e:
        elapsed = time.perf_counter() - t0
        return (elapsed, -1, False)


def benchmark_rdkit(pairs: list[dict], rounds: int, timeout: int) -> dict:
    """Run RDKit MCS on all pairs, 'rounds' iterations each. Returns dict of pair_key -> results."""
    if not HAS_RDKIT:
        return {}
    results = {}
    total = len(pairs)
    for pi, pair in enumerate(pairs):
        pair_key = f"{pair['name_a']}__vs__{pair['name_b']}"
        times = []
        mcs_sizes = []
        completed_count = 0
        for r in range(rounds):
            elapsed, mcs_size, completed = rdkit_mcs_single(pair["smi_a"], pair["smi_b"], timeout)
            if elapsed >= 0:
                times.append(elapsed)
            if mcs_size >= 0:
                mcs_sizes.append(mcs_size)
            if completed:
                completed_count += 1
        results[pair_key] = {
            "median_time": median(times) if times else -1,
            "mean_time": mean(times) if times else -1,
            "mcs_size": max(mcs_sizes) if mcs_sizes else -1,
            "completed": completed_count,
            "rounds": rounds,
        }
        if (pi + 1) % 50 == 0 or pi == total - 1:
            print(f"  RDKit: {pi + 1}/{total} pairs done", file=sys.stderr)
    return results


# ---------------------------------------------------------------------------
# SMSD Java MCS benchmark
# ---------------------------------------------------------------------------
def find_smsd_jar(user_path: str | None) -> Path | None:
    """Locate the SMSD jar file."""
    if user_path:
        p = Path(user_path)
        if p.exists():
            return p
    # Try default location
    if DEFAULT_SMSD_JAR.exists():
        return DEFAULT_SMSD_JAR
    # Try globbing target/
    target = SCRIPT_DIR.parent / "target"
    if target.is_dir():
        jars = sorted(target.glob("smsd-*.jar"))
        if jars:
            return jars[-1]
    return None


def smsd_mcs_single(smi_a: str, smi_b: str, smsd_jar: Path, timeout_sec: int = 10) -> tuple[float, int, bool]:
    """
    Run SMSD via CLI on a single pair.
    Returns (elapsed_seconds, mcs_atom_count, completed).
    """
    cmd = [
        "java", "-cp", str(smsd_jar),
        "com.bioinception.smsd.cli.SMSDcli",
        "-q", smi_a,
        "-t", smi_b,
        "-o", "MCS",
        "-I", "SMI",
    ]
    t0 = time.perf_counter()
    try:
        proc = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=timeout_sec + 2,
        )
        elapsed = time.perf_counter() - t0
        # Parse output for MCS size
        mcs_size = -1
        for line in proc.stdout.splitlines():
            if "MCS_SIZE" in line or "mcs_size" in line.lower():
                parts = line.split()
                for p in parts:
                    try:
                        mcs_size = int(p)
                    except ValueError:
                        pass
            elif "AtomCount" in line or "atom_count" in line.lower():
                parts = line.split()
                for p in parts:
                    try:
                        mcs_size = int(p)
                    except ValueError:
                        pass
        return (elapsed, mcs_size, proc.returncode == 0)
    except subprocess.TimeoutExpired:
        elapsed = time.perf_counter() - t0
        return (elapsed, -1, False)
    except Exception:
        elapsed = time.perf_counter() - t0
        return (elapsed, -1, False)


def benchmark_smsd(pairs: list[dict], rounds: int, timeout: int, smsd_jar: Path) -> dict:
    """Run SMSD MCS on all pairs. Returns dict of pair_key -> results."""
    results = {}
    total = len(pairs)
    for pi, pair in enumerate(pairs):
        pair_key = f"{pair['name_a']}__vs__{pair['name_b']}"
        times = []
        mcs_sizes = []
        completed_count = 0
        for r in range(rounds):
            elapsed, mcs_size, completed = smsd_mcs_single(
                pair["smi_a"], pair["smi_b"], smsd_jar, timeout
            )
            if elapsed >= 0:
                times.append(elapsed)
            if mcs_size >= 0:
                mcs_sizes.append(mcs_size)
            if completed:
                completed_count += 1
        results[pair_key] = {
            "median_time": median(times) if times else -1,
            "mean_time": mean(times) if times else -1,
            "mcs_size": max(mcs_sizes) if mcs_sizes else -1,
            "completed": completed_count,
            "rounds": rounds,
        }
        if (pi + 1) % 50 == 0 or pi == total - 1:
            print(f"  SMSD:  {pi + 1}/{total} pairs done", file=sys.stderr)
    return results


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------
def write_results(
    pairs: list[dict],
    rdkit_results: dict,
    smsd_results: dict,
    output_path: Path,
):
    """Write combined TSV with all results."""
    with open(output_path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow([
            "pair_id",
            "pair_type",
            "name_a",
            "name_b",
            "smi_a",
            "smi_b",
            "rdkit_median_s",
            "rdkit_mean_s",
            "rdkit_mcs_size",
            "rdkit_completed",
            "smsd_median_s",
            "smsd_mean_s",
            "smsd_mcs_size",
            "smsd_completed",
            "size_diff",
            "winner",
        ])
        for pi, pair in enumerate(pairs):
            pair_key = f"{pair['name_a']}__vs__{pair['name_b']}"
            rr = rdkit_results.get(pair_key, {})
            sr = smsd_results.get(pair_key, {})

            rdkit_med = rr.get("median_time", -1)
            smsd_med = sr.get("median_time", -1)

            if rdkit_med > 0 and smsd_med > 0:
                if smsd_med < rdkit_med * 0.9:
                    winner = "SMSD"
                elif rdkit_med < smsd_med * 0.9:
                    winner = "RDKit"
                else:
                    winner = "TIE"
            elif rdkit_med > 0:
                winner = "RDKit_ONLY"
            elif smsd_med > 0:
                winner = "SMSD_ONLY"
            else:
                winner = "NONE"

            mcs_diff = ""
            rdkit_mcs = rr.get("mcs_size", -1)
            smsd_mcs = sr.get("mcs_size", -1)
            if rdkit_mcs >= 0 and smsd_mcs >= 0:
                mcs_diff = smsd_mcs - rdkit_mcs

            w.writerow([
                pi + 1,
                pair["pair_type"],
                pair["name_a"],
                pair["name_b"],
                pair["smi_a"],
                pair["smi_b"],
                f"{rdkit_med:.6f}" if rdkit_med >= 0 else "N/A",
                f"{rr.get('mean_time', -1):.6f}" if rr.get('mean_time', -1) >= 0 else "N/A",
                rdkit_mcs if rdkit_mcs >= 0 else "N/A",
                f"{rr.get('completed', 0)}/{rr.get('rounds', 0)}",
                f"{smsd_med:.6f}" if smsd_med >= 0 else "N/A",
                f"{sr.get('mean_time', -1):.6f}" if sr.get('mean_time', -1) >= 0 else "N/A",
                smsd_mcs if smsd_mcs >= 0 else "N/A",
                f"{sr.get('completed', 0)}/{sr.get('rounds', 0)}",
                mcs_diff,
                winner,
            ])


def write_summary(
    pairs: list[dict],
    rdkit_results: dict,
    smsd_results: dict,
    mols: list,
    output_path: Path,
):
    """Write human-readable summary."""
    lines = []
    lines.append("=" * 70)
    lines.append("SMSD 1000-Molecule Benchmark Summary")
    lines.append("=" * 70)
    lines.append(f"Total molecules loaded: {len(mols)}")
    lines.append(f"Total pairs tested:     {len(pairs)}")
    lines.append("")

    for engine_name, results in [("RDKit", rdkit_results), ("SMSD", smsd_results)]:
        if not results:
            lines.append(f"{engine_name}: NOT TESTED (unavailable)")
            lines.append("")
            continue

        all_times = [v["median_time"] for v in results.values() if v["median_time"] >= 0]
        all_mcs = [v["mcs_size"] for v in results.values() if v["mcs_size"] >= 0]
        completed = sum(1 for v in results.values() if v["completed"] == v["rounds"])
        total = len(results)

        lines.append(f"--- {engine_name} ---")
        if all_times:
            lines.append(f"  Median time (across pairs):  {median(all_times)*1000:.2f} ms")
            lines.append(f"  Mean time (across pairs):    {mean(all_times)*1000:.2f} ms")
            lines.append(f"  Min time:                    {min(all_times)*1000:.2f} ms")
            lines.append(f"  Max time:                    {max(all_times)*1000:.2f} ms")
        if all_mcs:
            lines.append(f"  Median MCS size:             {median(all_mcs):.1f}")
            lines.append(f"  Mean MCS size:               {mean(all_mcs):.1f}")
        lines.append(f"  Completion rate:             {completed}/{total} ({100*completed/total:.1f}%)")
        lines.append("")

    # Per-section breakdown
    if rdkit_results or smsd_results:
        lines.append("--- Per-Section Breakdown (median time ms) ---")
        lines.append(f"{'Section':<15} {'RDKit':>12} {'SMSD':>12} {'Winner':>10}")
        lines.append("-" * 55)

        for section_name, (lo, hi) in SECTION_RANGES.items():
            rdkit_sec = []
            smsd_sec = []
            for pair in pairs:
                pair_key = f"{pair['name_a']}__vs__{pair['name_b']}"
                sec_a = get_section(pair["idx_a"] + 1)
                sec_b = get_section(pair["idx_b"] + 1)
                if sec_a == section_name or sec_b == section_name:
                    rr = rdkit_results.get(pair_key, {})
                    sr = smsd_results.get(pair_key, {})
                    if rr.get("median_time", -1) >= 0:
                        rdkit_sec.append(rr["median_time"])
                    if sr.get("median_time", -1) >= 0:
                        smsd_sec.append(sr["median_time"])

            r_str = f"{median(rdkit_sec)*1000:.2f}" if rdkit_sec else "N/A"
            s_str = f"{median(smsd_sec)*1000:.2f}" if smsd_sec else "N/A"
            if rdkit_sec and smsd_sec:
                r_med = median(rdkit_sec)
                s_med = median(smsd_sec)
                winner = "SMSD" if s_med < r_med * 0.9 else ("RDKit" if r_med < s_med * 0.9 else "TIE")
            else:
                winner = "-"
            lines.append(f"{section_name:<15} {r_str:>12} {s_str:>12} {winner:>10}")
        lines.append("")

    # Head-to-head
    if rdkit_results and smsd_results:
        smsd_wins = rdkit_wins = ties = 0
        smsd_bigger_mcs = rdkit_bigger_mcs = equal_mcs = 0
        for pair in pairs:
            pair_key = f"{pair['name_a']}__vs__{pair['name_b']}"
            rr = rdkit_results.get(pair_key, {})
            sr = smsd_results.get(pair_key, {})
            rt = rr.get("median_time", -1)
            st = sr.get("median_time", -1)
            if rt > 0 and st > 0:
                if st < rt * 0.9:
                    smsd_wins += 1
                elif rt < st * 0.9:
                    rdkit_wins += 1
                else:
                    ties += 1
            rm = rr.get("mcs_size", -1)
            sm = sr.get("mcs_size", -1)
            if rm >= 0 and sm >= 0:
                if sm > rm:
                    smsd_bigger_mcs += 1
                elif rm > sm:
                    rdkit_bigger_mcs += 1
                else:
                    equal_mcs += 1

        lines.append("--- Head-to-Head Comparison ---")
        lines.append(f"  Speed: SMSD faster: {smsd_wins}, RDKit faster: {rdkit_wins}, Tie: {ties}")
        lines.append(f"  Quality: SMSD larger MCS: {smsd_bigger_mcs}, RDKit larger: {rdkit_bigger_mcs}, Equal: {equal_mcs}")

    text = "\n".join(lines)
    with open(output_path, "w") as fh:
        fh.write(text + "\n")
    print(text)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description="SMSD 1000-Molecule Benchmark")
    parser.add_argument("--pairs", type=int, default=1000, help="Total pairs (split 50/50 random/systematic)")
    parser.add_argument("--rounds", type=int, default=5, help="Rounds per pair (report median)")
    parser.add_argument("--timeout", type=int, default=10, help="Timeout in seconds per pair")
    parser.add_argument("--smsd-jar", type=str, default=None, help="Path to SMSD jar")
    parser.add_argument("--rdkit-only", action="store_true", help="Skip SMSD Java benchmark")
    parser.add_argument("--smsd-only", action="store_true", help="Skip RDKit benchmark")
    parser.add_argument("--seed", type=int, default=42, help="Random seed for pair generation")
    parser.add_argument("--molecules", type=str, default=None, help="Path to molecules file")
    args = parser.parse_args()

    mol_path = Path(args.molecules) if args.molecules else MOLECULES_FILE
    print(f"Loading molecules from {mol_path} ...", file=sys.stderr)
    raw_mols = load_molecules(mol_path)
    print(f"  Loaded {len(raw_mols)} molecules", file=sys.stderr)

    print("Validating SMILES ...", file=sys.stderr)
    valid_mols = validate_smiles(raw_mols)
    print(f"  {len(valid_mols)} valid molecules", file=sys.stderr)

    if len(valid_mols) < 10:
        print("ERROR: Too few valid molecules. Aborting.", file=sys.stderr)
        sys.exit(1)

    n_random = args.pairs // 2
    n_systematic = args.pairs - n_random
    print(f"Generating {n_random} random + {n_systematic} systematic pairs ...", file=sys.stderr)
    pairs = generate_pairs(valid_mols, n_random, n_systematic, args.seed)
    print(f"  Generated {len(pairs)} pairs", file=sys.stderr)

    # --- RDKit benchmark ---
    rdkit_results = {}
    if not args.smsd_only:
        if HAS_RDKIT:
            print(f"\nRunning RDKit FindMCS benchmark ({args.rounds} rounds, {args.timeout}s timeout) ...", file=sys.stderr)
            rdkit_results = benchmark_rdkit(pairs, args.rounds, args.timeout)
        else:
            print("Skipping RDKit (not installed)", file=sys.stderr)

    # --- SMSD benchmark ---
    smsd_results = {}
    if not args.rdkit_only:
        smsd_jar = find_smsd_jar(args.smsd_jar)
        if smsd_jar:
            print(f"\nRunning SMSD Java benchmark (jar: {smsd_jar}) ...", file=sys.stderr)
            print(f"  ({args.rounds} rounds, {args.timeout}s timeout)", file=sys.stderr)
            smsd_results = benchmark_smsd(pairs, args.rounds, args.timeout, smsd_jar)
        else:
            print("Skipping SMSD Java (jar not found; build with: mvn package -DskipTests)", file=sys.stderr)

    # --- Output ---
    print(f"\nWriting results to {OUTPUT_TSV} ...", file=sys.stderr)
    write_results(pairs, rdkit_results, smsd_results, OUTPUT_TSV)

    print(f"Writing summary to {OUTPUT_SUMMARY} ...\n", file=sys.stderr)
    write_summary(pairs, rdkit_results, smsd_results, valid_mols, OUTPUT_SUMMARY)


if __name__ == "__main__":
    main()
