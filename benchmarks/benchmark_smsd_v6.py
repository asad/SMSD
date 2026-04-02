#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
# Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
# See the NOTICE file for attribution, trademark, and algorithm IP terms.
"""
SMSD v6 Performance Benchmark Suite

Times three core operations:
  1. MCS computation for 10 standard molecule pairs
  2. ECFP4 fingerprint generation for 100 molecules
  3. Batch substructure screening

Results are printed as a formatted table and saved to benchmarks/results_v6.tsv.

Usage:
    python3 benchmarks/benchmark_smsd_v6.py
"""

import statistics
import sys
import time
from pathlib import Path

# ---------------------------------------------------------------------------
# Import SMSD
# ---------------------------------------------------------------------------

try:
    import smsd
except ImportError:
    sys.exit("ERROR: smsd not installed. Run: pip install smsd")

# Circular fingerprint (ECFP) import
try:
    from smsd._smsd import circular_fingerprint as _circular_fp
except ImportError:
    _circular_fp = None

SCRIPT_DIR = Path(__file__).resolve().parent

# ---------------------------------------------------------------------------
# 10 Standard MCS Pairs
# ---------------------------------------------------------------------------

MCS_PAIRS = [
    ("c1ccccc1",            "Cc1ccccc1",           "benzene / toluene"),
    ("c1ccccc1",            "c1ccc(O)cc1",         "benzene / phenol"),
    ("c1ccncc1",            "c1ccnc(N)c1",         "pyridine / 2-aminopyridine"),
    ("C1CCCCC1",            "C1CCC(O)CC1",         "cyclohexane / cyclohexanol"),
    ("c1ccc2ccccc2c1",      "c1ccc2cc(O)ccc2c1",   "naphthalene / 2-naphthol"),
    ("CC(=O)O",             "CC(=O)OC",            "acetic acid / methyl acetate"),
    ("c1ccc(cc1)C(=O)O",    "c1ccc(cc1)C(=O)N",   "benzoic acid / benzamide"),
    # morphine / codeine (partial structures for speed)
    ("C1CC2=CC=C(C3CC4N(CC13)CCC24)O",
     "C1CC2=CC=C(C3CC4N(CC13)CCC24)OC",
     "morphine core / codeine core"),
    ("c1ccc(-c2ccccc2)cc1", "c1ccc(-c2ccc(O)cc2)cc1", "biphenyl / 4-hydroxybiphenyl"),
    ("CC(C)Cc1ccc(cc1)C(C)C(=O)O",
     "CC(C)Cc1ccc(cc1)C(C)C(=O)O",
     "ibuprofen / ibuprofen (self-match)"),
]


def load_diverse_molecules(n=100):
    """Load the first n molecules from diverse_molecules.txt."""
    mol_file = SCRIPT_DIR / "diverse_molecules.txt"
    if not mol_file.exists():
        return []
    mols = []
    with open(mol_file) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            smi = line.split("\t")[0].strip()
            if not smi:
                continue
            try:
                g = smsd.parse_smiles(smi)
                mols.append((smi, g))
            except Exception:
                continue
            if len(mols) >= n:
                break
    return mols


def time_func(func, warmup=2, iters=5):
    """Time a callable, returning (median_ms, all_times_ms)."""
    for _ in range(warmup):
        func()
    times = []
    for _ in range(iters):
        t0 = time.perf_counter()
        func()
        t1 = time.perf_counter()
        times.append((t1 - t0) * 1000.0)
    return statistics.median(times), times


# ---------------------------------------------------------------------------
# Benchmark 1: MCS for 10 standard pairs
# ---------------------------------------------------------------------------

def benchmark_mcs():
    """Benchmark MCS computation for 10 standard pairs."""
    results = []
    for smi1, smi2, name in MCS_PAIRS:
        g1 = smsd.parse_smiles(smi1)
        g2 = smsd.parse_smiles(smi2)

        def run_mcs(a=g1, b=g2):
            return smsd.find_mcs(a, b)

        median_ms, _ = time_func(run_mcs, warmup=2, iters=5)

        mapping = smsd.find_mcs(g1, g2)
        mcs_size = len(mapping)
        results.append((name, mcs_size, median_ms))
    return results


# ---------------------------------------------------------------------------
# Benchmark 2: ECFP4 fingerprint generation for 100 molecules
# ---------------------------------------------------------------------------

def benchmark_ecfp4():
    """Benchmark ECFP4 fingerprint generation for up to 100 molecules."""
    molecules = load_diverse_molecules(100)
    if not molecules:
        return None

    if _circular_fp is None:
        # Fall back to path fingerprint
        graphs = [g for _, g in molecules]

        def run_fp():
            for g in graphs:
                smsd.fingerprint(g, kind="path")

        median_ms, _ = time_func(run_fp, warmup=1, iters=3)
        return ("path_fp", len(molecules), median_ms)

    graphs = [g for _, g in molecules]

    def run_ecfp():
        for g in graphs:
            _circular_fp(g, 2, 2048, "ecfp")

    median_ms, _ = time_func(run_ecfp, warmup=1, iters=3)
    return ("ecfp4", len(molecules), median_ms)


# ---------------------------------------------------------------------------
# Benchmark 3: Batch substructure screening
# ---------------------------------------------------------------------------

def benchmark_batch_substructure():
    """Benchmark batch substructure screening: benzene ring vs 100 molecules."""
    molecules = load_diverse_molecules(100)
    if not molecules:
        return None

    query = smsd.parse_smiles("c1ccccc1")
    targets = [g for _, g in molecules]

    def run_batch():
        return smsd.batch_substructure(query, targets)

    median_ms, _ = time_func(run_batch, warmup=1, iters=3)
    hits = sum(smsd.batch_substructure(query, targets))
    return (len(molecules), hits, median_ms)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    print("=" * 72)
    print("SMSD v6 Performance Benchmark Suite")
    print(f"SMSD version: {smsd.__version__}")
    print("=" * 72)
    print()

    tsv_lines = ["benchmark\tdetail\tvalue\tunit"]

    # --- MCS Benchmark ---
    print("--- MCS: 10 standard molecule pairs ---")
    print(f"{'Pair':<45} {'MCS':>5} {'Time (ms)':>10}")
    print("-" * 62)
    mcs_results = benchmark_mcs()
    total_mcs_ms = 0.0
    for name, size, ms in mcs_results:
        print(f"{name:<45} {size:>5d} {ms:>10.2f}")
        tsv_lines.append(f"mcs\t{name}\t{ms:.3f}\tms")
        total_mcs_ms += ms
    print(f"{'TOTAL':<45} {'':>5} {total_mcs_ms:>10.2f}")
    tsv_lines.append(f"mcs\ttotal\t{total_mcs_ms:.3f}\tms")
    print()

    # --- ECFP4 Benchmark ---
    print("--- Fingerprint: ECFP4 / Path FP for 100 molecules ---")
    fp_result = benchmark_ecfp4()
    if fp_result:
        fp_kind, n_mols, median_ms = fp_result
        per_mol = median_ms / n_mols if n_mols > 0 else 0
        print(f"  Type:       {fp_kind}")
        print(f"  Molecules:  {n_mols}")
        print(f"  Total:      {median_ms:.2f} ms")
        print(f"  Per mol:    {per_mol:.3f} ms")
        tsv_lines.append(f"fingerprint\t{fp_kind}_total_{n_mols}mols\t{median_ms:.3f}\tms")
        tsv_lines.append(f"fingerprint\t{fp_kind}_per_mol\t{per_mol:.4f}\tms")
    else:
        print("  SKIPPED: diverse_molecules.txt not found")
    print()

    # --- Batch Substructure Benchmark ---
    print("--- Batch substructure: benzene ring vs 100 molecules ---")
    batch_result = benchmark_batch_substructure()
    if batch_result:
        n_targets, hits, median_ms = batch_result
        per_target = median_ms / n_targets if n_targets > 0 else 0
        print(f"  Targets:    {n_targets}")
        print(f"  Hits:       {hits}")
        print(f"  Total:      {median_ms:.2f} ms")
        print(f"  Per target: {per_target:.3f} ms")
        tsv_lines.append(f"batch_sub\ttotal_{n_targets}targets\t{median_ms:.3f}\tms")
        tsv_lines.append(f"batch_sub\thits\t{hits}\tcount")
        tsv_lines.append(f"batch_sub\tper_target\t{per_target:.4f}\tms")
    else:
        print("  SKIPPED: diverse_molecules.txt not found")
    print()

    # --- Save results ---
    out_file = SCRIPT_DIR / "results_v6.tsv"
    with open(out_file, "w") as f:
        f.write("\n".join(tsv_lines) + "\n")
    print(f"Results saved to {out_file}")
    print("=" * 72)


if __name__ == "__main__":
    main()
