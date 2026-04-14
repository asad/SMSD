#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
# Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
# See the NOTICE file for attribution, trademark, and algorithm IP terms.
"""
SMSD vs RDKit FindMCS -- Reproducible Benchmark Suite

Benchmarks RDKit FindMCS AND SMSD Java CLI on the same 20 molecule pairs,
spanning trivial (methane/ethane) to massive (paclitaxel/docetaxel) and
algorithmically hard cases (symmetric cages, self-matches, known failures).

Usage:
    pip install rdkit
    cd SMSD && mvn package -DskipTests          # build Java jar once
    python3 benchmarks/benchmark_all.py          # run everything

Output:
    - benchmarks/results_smsd_vs_rdkit_20pairs_mcs_all.tsv          (machine-readable)
    - benchmarks/results_smsd_vs_rdkit_20pairs_summary.txt  (human-readable formatted table)
    - stdout: formatted comparison table

Requirements:
    - Python 3.8+
    - rdkit  (pip install rdkit)
    - Java 11+  (for SMSD CLI)
    - Built SMSD shaded jar at target/smsd-*-jar-with-dependencies.jar

Author: Syed Asad Rahman, BioInception PVT LTD
"""

import json
import os
import platform
import subprocess
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path
from statistics import median
from typing import List, Optional, Tuple

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_DIR = SCRIPT_DIR.parent
JAR_PATH = PROJECT_DIR / "target"
# Fallback: check src/scripts/repo for older shaded jars
JAR_FALLBACK_DIR = PROJECT_DIR / "src" / "scripts" / "repo"

CPP_BUILD_DIR = PROJECT_DIR / "cpp" / "build"
CPP_BINARY = CPP_BUILD_DIR / "smsd_benchmark"

NUM_RUNS = 5
TIMEOUT_SEC = 10  # per-pair timeout for RDKit
TIMEOUT_MS = 10_000  # per-pair timeout for SMSD Java

# Vancomycin SMILES (pair #14)
VANCOMYCIN = (
    "CC1C(C(CC(O1)OC2C(C(C(OC2OC3=C4C=C5C=C3OC6=C(C=C(C=C6)C(C(C(=O)"
    "NC(C(=O)NC5C(=O)NC7C8=CC(=C(C=C8)O)C9=C(C=C(C=C9O)O)C(NC(=O)C("
    "C(C1=CC(=C(O4)C=C1)Cl)O)NC7=O)C(=O)O)CC(=O)N)NC(=O)C(CC(C)C)NC)"
    "O)Cl)CO)O)O)(C)N)O"
)

# PEG-16 SMILES (pair #17)
PEG16 = "OCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCO"

# ---------------------------------------------------------------------------
# 20 Molecule Pairs
# ---------------------------------------------------------------------------

PAIRS: List[Tuple[str, str, str, str]] = [
    # (query_smiles, target_smiles, pair_name, category)

    # 1. Trivial
    ("C", "CC", "methane-ethane", "Trivial"),

    # 2. Small aromatic
    ("c1ccccc1", "Cc1ccccc1", "benzene-toluene", "Small aromatic"),

    # 3. Heteroatom
    ("c1ccccc1", "Oc1ccccc1", "benzene-phenol", "Heteroatom"),

    # 4. Drug pair
    ("CC(=O)Oc1ccccc1C(=O)O", "CC(=O)Nc1ccc(O)cc1", "aspirin-acetaminophen", "Drug pair"),

    # 5. N-methyl diff
    ("Cn1cnc2c1c(=O)n(C)c(=O)n2C", "Cn1cnc2c1c(=O)[nH]c(=O)n2C",
     "caffeine-theophylline", "N-methyl diff"),

    # 6. Alkaloid
    ("CN1CCC23C4C1CC5=C(C2C(C=C4)O3)C=C(C=C5)O",
     "CN1CCC23C4C1CC5=C(C2C(C=C4)OC3)C=C(C=C5)O",
     "morphine-codeine", "Alkaloid"),

    # 7. NSAID
    ("CC(C)Cc1ccc(CC(C)C(=O)O)cc1",
     "COc1ccc2cc(CC(C)C(=O)O)ccc2c1",
     "ibuprofen-naproxen", "NSAID"),

    # 8. Nucleotide
    ("C1=NC(=C2C(=N1)N(C=N2)C3C(C(C(O3)COP(=O)(O)OP(=O)(O)OP(=O)(O)O)O)O)N",
     "C1=NC(=C2C(=N1)N(C=N2)C3C(C(C(O3)COP(=O)(O)OP(=O)(O)O)O)O)N",
     "ATP-ADP", "Nucleotide"),

    # 9. Cofactor
    ("C1=CC(=C[N+](=C1)C2C(C(C(O2)COP(=O)([O-])OP(=O)(O)OCC3C(C(C(O3)N4C=NC5=C(N=CN=C54)N)O)O)O)O)C(=O)N",
     "C1C=CN(C=C1C(=O)N)C2C(C(C(O2)COP(=O)(O)OP(=O)(O)OCC3C(C(C(O3)N4C=NC5=C4N=CN=C5N)O)O)O)O",
     "NAD-NADH", "Cofactor"),

    # 10. Statin
    ("CC(C)C1=C(C(=C(N1CCC(CC(CC(=O)O)O)O)C2=CC=C(C=C2)F)C3=CC=CC=C3)C(=O)NC4=CC=CC=C4",
     "CC(C)C1=NC(=NC(=C1C=CC(CC(CC(=O)O)O)O)C2=CC=C(C=C2)F)N(C)S(=O)(=O)C",
     "atorvastatin-rosuvastatin", "Statin"),

    # 11. Taxane
    ("CC1=C2C(C(=O)C3(C(CC4C(C3C(C(C2(C)C)(CC1OC(=O)C(C(C5=CC=CC=C5)NC(=O)C6=CC=CC=C6)O)O)OC(=O)C7=CC=CC=C7)(CO4)OC(=O)C)O)C)OC(=O)C",
     "CC1=C2C(C(=O)C3(C(CC4C(C3C(C(C2(C)C)(CC1OC(=O)C(C(C5=CC=CC=C5)NC(=O)OC(C)(C)C)O)O)OC(=O)C6=CC=CC=C6)(CO4)OC(=O)C)O)C)O",
     "paclitaxel-docetaxel", "Taxane"),

    # 12. Macrolide
    ("CCC1C(C(C(C(=O)C(CC(C(C(C(C(C(=O)O1)C)OC2CC(C(C(O2)C)O)(C)OC)C)OC3C(C(CC(O3)C)N(C)C)O)(C)O)C)C)O)(C)O",
     "CCC1C(C(C(N(CC(CC(C(C(C(C(C(=O)O1)C)OC2CC(C(C(O2)C)O)(C)OC)C)OC3C(C(CC(O3)C)N(C)C)O)(C)O)C)C)C)O)(C)O",
     "erythromycin-azithromycin", "Macrolide"),

    # 13. Alkaloid scaffold
    ("C1CN2CC3=CCOC4CC(=O)N5C6C4C3CC2C61C7=CC=CC=C75",
     "COC1=CC2=C(C=CN=C2C=C1)C(C3CC4CCN3CC4C=C)O",
     "strychnine-quinine", "Alkaloid scaffold"),

    # 14. Self-match large
    (VANCOMYCIN, VANCOMYCIN, "vancomycin-self", "Self-match large"),

    # 15. Symmetric
    ("C1C2CC3CC1CC(C2)C3", "C1C2CC3CC1CC(C2)C3",
     "adamantane-self", "Symmetric"),

    # 16. Cage
    ("C12C3C4C1C5C4C3C25", "C12C3C4C1C5C4C3C25",
     "cubane-self", "Cage"),

    # 17. Polymer
    ("OCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCO", PEG16,
     "PEG12-PEG16", "Polymer"),

    # 18. PAH
    ("c1cc2ccc3ccc4ccc5ccc6ccc1c7c2c3c4c5c67",
     "c1cc2ccc3ccc4ccc5ccc6ccc1c7c2c3c4c5c67",
     "coronene-self", "PAH"),

    # 19. Tautomer
    ("O=c1[nH]c(N)nc2[nH]cnc12", "Oc1nc(N)nc2[nH]cnc12",
     "guanine-keto-enol", "Tautomer"),

    # 20. Known failure (RDKit #1585)
    ("c1cc(c(c(c1)Cl)N2c3cc(cc(c3CNC2=O)c4ccc(cc4F)F)N5CCNCC5)Cl",
     "CCNc1cc(c2c(c1)N(C(=O)NC2)c3ccc(cc3)n4ccc-5ncnc5c4)c6ccnnc6",
     "rdkit-1585-pair", "Known failure"),
]


# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------

@dataclass
class BenchResult:
    pair_name: str
    category: str
    rdkit_best_ms: float = -1.0
    rdkit_median_ms: float = -1.0
    rdkit_mcs_size: int = -1
    rdkit_mcs_smarts: str = ""
    rdkit_timed_out: bool = False
    smsd_best_ms: float = -1.0
    smsd_median_ms: float = -1.0
    smsd_mcs_size: int = -1
    smsd_error: str = ""
    cpp_best_ms: float = -1.0
    cpp_median_ms: float = -1.0
    cpp_mcs_size: int = -1
    py_best_ms: float = -1.0
    py_median_ms: float = -1.0
    py_mcs_size: int = -1


# ---------------------------------------------------------------------------
# RDKit benchmark
# ---------------------------------------------------------------------------

def benchmark_rdkit(pairs) -> List[BenchResult]:
    """Run RDKit FindMCS on all pairs. Returns list of BenchResult (rdkit fields filled)."""
    try:
        from rdkit import Chem
        from rdkit.Chem import rdFMCS
    except ImportError:
        print("[WARN] RDKit not installed -- skipping RDKit benchmarks.", file=sys.stderr)
        return [BenchResult(pair_name=name, category=cat) for _, _, name, cat in pairs]

    results = []
    for smi1, smi2, name, cat in pairs:
        res = BenchResult(pair_name=name, category=cat)
        mol1 = Chem.MolFromSmiles(smi1)
        mol2 = Chem.MolFromSmiles(smi2)
        if mol1 is None or mol2 is None:
            res.rdkit_mcs_size = -1
            res.rdkit_best_ms = -1
            results.append(res)
            continue

        # Warmup
        try:
            rdFMCS.FindMCS([mol1, mol2], timeout=1)
        except Exception:
            pass

        times = []
        mcs_size = 0
        mcs_smarts = ""
        timed_out = False
        for _ in range(NUM_RUNS):
            t0 = time.perf_counter()
            try:
                result = rdFMCS.FindMCS([mol1, mol2], timeout=TIMEOUT_SEC)
                mcs_size = result.numAtoms
                mcs_smarts = result.smartsString if result.smartsString else ""
                timed_out = result.canceled
            except Exception:
                mcs_size = -1
            dt = (time.perf_counter() - t0) * 1000.0
            times.append(dt)

        res.rdkit_best_ms = min(times)
        res.rdkit_median_ms = median(times)
        res.rdkit_mcs_size = mcs_size
        res.rdkit_mcs_smarts = mcs_smarts
        res.rdkit_timed_out = timed_out
        results.append(res)

    return results


# ---------------------------------------------------------------------------
# SMSD Java benchmark
# ---------------------------------------------------------------------------

def find_smsd_jar() -> Optional[Path]:
    """Locate the SMSD fat jar."""
    if JAR_PATH.exists():
        jars = sorted(JAR_PATH.glob("smsd-*-jar-with-dependencies.jar"), reverse=True)
        if jars:
            return jars[0]
    # Try fallback jars in descending version order
    if JAR_FALLBACK_DIR.exists():
        jars = sorted(JAR_FALLBACK_DIR.glob("smsd-*-jar-with-dependencies.jar"), reverse=True)
        if jars:
            return jars[0]
    return None


def run_smsd_java(jar: Path, smi1: str, smi2: str) -> Tuple[float, int, str]:
    """
    Run SMSD Java CLI for a single pair.
    Returns (elapsed_ms, mcs_size, error_string).
    """
    # Find Java — try JAVA_HOME, then common Homebrew paths, then system PATH
    java_bin = "java"
    for java_candidate in [
        os.environ.get("JAVA_HOME", "") + "/bin/java",
        "/usr/local/opt/openjdk/bin/java",
        "/opt/homebrew/opt/openjdk/bin/java",
        "/usr/bin/java",
    ]:
        if java_candidate and os.path.isfile(java_candidate):
            java_bin = java_candidate
            break
    cmd = [
        java_bin, "-jar", str(jar),
        "--Q", "SMI", "--q", smi1,
        "--T", "SMI", "--t", smi2,
        "--mode", "mcs",
        "--timeout", str(TIMEOUT_MS),
        "--json", "-",
    ]
    t0 = time.perf_counter()
    try:
        proc = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=TIMEOUT_SEC + 30,  # generous OS-level timeout
        )
        elapsed = (time.perf_counter() - t0) * 1000.0
        if proc.returncode != 0 and not proc.stdout.strip():
            return elapsed, -1, proc.stderr.strip()[:200]
        # Parse JSON
        try:
            data = json.loads(proc.stdout)
            mcs_size = data.get("mcs_size", data.get("mappingSize", -1))
            return elapsed, mcs_size, ""
        except json.JSONDecodeError:
            return elapsed, -1, "JSON parse error"
    except subprocess.TimeoutExpired:
        elapsed = (time.perf_counter() - t0) * 1000.0
        return elapsed, -1, "OS timeout"
    except FileNotFoundError:
        return -1, -1, "java not found"


def benchmark_smsd_java(pairs, existing_results: List[BenchResult]) -> List[BenchResult]:
    """Fill in SMSD Java fields on existing results."""
    jar = find_smsd_jar()
    if jar is None:
        print(f"[WARN] SMSD jar not found at {JAR_PATH}", file=sys.stderr)
        print("       Build with: cd SMSD && mvn package -DskipTests", file=sys.stderr)
        return existing_results

    print(f"[INFO] Using SMSD jar: {jar}", file=sys.stderr)

    # Check java is available
    try:
        subprocess.run(["java", "-version"], capture_output=True, timeout=10)
    except (FileNotFoundError, subprocess.TimeoutExpired):
        print("[WARN] Java not found -- skipping SMSD Java benchmarks.", file=sys.stderr)
        return existing_results

    for i, (smi1, smi2, name, cat) in enumerate(pairs):
        res = existing_results[i]

        # Warmup run (discard)
        run_smsd_java(jar, smi1, smi2)

        times = []
        mcs_size = -1
        error = ""
        for _ in range(NUM_RUNS):
            elapsed, sz, err = run_smsd_java(jar, smi1, smi2)
            times.append(elapsed)
            if sz >= 0:
                mcs_size = sz
            if err:
                error = err

        if times and all(t >= 0 for t in times):
            res.smsd_best_ms = min(times)
            res.smsd_median_ms = median(times)
        res.smsd_mcs_size = mcs_size
        res.smsd_error = error

    return existing_results


# ---------------------------------------------------------------------------
# C++ benchmark (optional)
# ---------------------------------------------------------------------------

def benchmark_cpp_if_available(existing_results: List[BenchResult]) -> List[BenchResult]:
    """If the C++ benchmark binary exists, run it and parse results."""
    if not CPP_BINARY.exists():
        print(f"[INFO] C++ binary not found at {CPP_BINARY} -- skipping.", file=sys.stderr)
        print("       Build with: cd cpp && mkdir -p build && cd build && cmake .. -DCMAKE_BUILD_TYPE=Release && make",
              file=sys.stderr)
        return existing_results

    try:
        proc = subprocess.run(
            [str(CPP_BINARY)],
            capture_output=True, text=True, timeout=300,
        )
        # Parse the C++ benchmark output (expects TSV-like lines)
        for line in proc.stdout.splitlines():
            parts = line.split()
            if len(parts) < 4:
                continue
            pair_name = parts[0]
            for res in existing_results:
                if res.pair_name == pair_name:
                    try:
                        res.cpp_best_ms = float(parts[1])
                        res.cpp_median_ms = float(parts[2])
                        res.cpp_mcs_size = int(parts[3])
                    except (ValueError, IndexError):
                        pass
                    break
    except (subprocess.TimeoutExpired, FileNotFoundError, OSError) as e:
        print(f"[WARN] C++ benchmark failed: {e}", file=sys.stderr)

    return existing_results


# ---------------------------------------------------------------------------
# SMSD Python (pip install smsd) benchmark
# ---------------------------------------------------------------------------

def benchmark_smsd_python(pairs, existing_results: List[BenchResult]) -> List[BenchResult]:
    """Benchmark the smsd Python package (C++ bindings via pybind11)."""
    try:
        from smsd import parse_smiles, find_mcs, ChemOptions, MCSOptions
    except ImportError:
        print("[WARN] smsd Python package not found -- pip install smsd", file=sys.stderr)
        return existing_results

    for i, (smi1, smi2, name, _cat) in enumerate(pairs):
        res = existing_results[i]
        try:
            g1 = parse_smiles(smi1)
            g2 = parse_smiles(smi2)
        except Exception as e:
            res.py_mcs_size = -1
            continue

        chem = ChemOptions()
        opts = MCSOptions()
        opts.timeout_ms = TIMEOUT_MS

        # Warmup
        try:
            find_mcs(g1, g2, chem, opts)
        except Exception:
            pass

        times = []
        mcs_size = -1
        for _ in range(NUM_RUNS):
            t0 = time.perf_counter()
            try:
                mapping = find_mcs(g1, g2, chem, opts)
                mcs_size = len(mapping)
            except Exception:
                mcs_size = -1
            times.append((time.perf_counter() - t0) * 1000.0)

        if times:
            res.py_best_ms   = min(times)
            res.py_median_ms = median(times)
        res.py_mcs_size = mcs_size

    return existing_results


# ---------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------

def winner_and_speedup(rdkit_ms: float, smsd_ms: float) -> Tuple[str, str]:
    """Determine winner and speedup factor."""
    if rdkit_ms < 0 and smsd_ms < 0:
        return "N/A", ""
    if rdkit_ms < 0:
        return "SMSD", ""
    if smsd_ms < 0:
        return "RDKit", ""
    if rdkit_ms <= 0.001 and smsd_ms <= 0.001:
        return "tie", "~1x"
    if smsd_ms <= 0.001:
        return "SMSD", ">999x"
    if rdkit_ms <= 0.001:
        return "RDKit", ">999x"
    ratio = rdkit_ms / smsd_ms
    if ratio > 1.1:
        return "SMSD", f"{ratio:.1f}x"
    elif ratio < 0.9:
        return "RDKit", f"{1.0/ratio:.1f}x"
    else:
        return "tie", f"~1x"


def mcs_quality(rdkit_size: int, smsd_size: int) -> str:
    """Compare MCS sizes."""
    if rdkit_size < 0 and smsd_size < 0:
        return "both failed"
    if rdkit_size < 0:
        return "RDKit failed"
    if smsd_size < 0:
        return "SMSD failed"
    if rdkit_size == smsd_size:
        return "equal"
    diff = smsd_size - rdkit_size
    if diff > 0:
        return f"SMSD +{diff}"
    else:
        return f"RDKit +{-diff}"


def format_ms(ms: float) -> str:
    """Format milliseconds for display."""
    if ms < 0:
        return "N/A"
    if ms < 1.0:
        return f"{ms:.3f}"
    if ms < 100:
        return f"{ms:.1f}"
    return f"{ms:.0f}"


def print_results(results: List[BenchResult]):
    """Print formatted comparison table to stdout."""
    hdr = (
        f"{'#':>2s}  {'Pair':<28s}  {'Category':<18s}  "
        f"{'RDKit ms':>10s}  {'R.MCS':>5s}  "
        f"{'SMSD ms':>10s}  {'S.MCS':>5s}  "
        f"{'Py ms':>8s}  {'P.MCS':>5s}  "
        f"{'Winner':<6s}  {'Speed':>7s}  {'Quality':<14s}"
    )
    sep = "-" * len(hdr)

    print()
    print("=" * len(hdr))
    print("SMSD vs RDKit FindMCS -- Benchmark Results")
    print(f"Date: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"Platform: {platform.system()} {platform.machine()}")
    print(f"Python: {platform.python_version()}")
    try:
        from rdkit import rdBase
        print(f"RDKit: {rdBase.rdkitVersion}")
    except ImportError:
        print("RDKit: not installed")
    print(f"Runs per pair: {NUM_RUNS}  (best time reported)")
    print(f"Timeout: {TIMEOUT_SEC}s (RDKit), {TIMEOUT_MS}ms (SMSD)")
    print("=" * len(hdr))
    print()
    print(hdr)
    print(sep)

    rdkit_total = 0.0
    smsd_total = 0.0
    rdkit_wins = 0
    smsd_wins = 0
    ties = 0
    quality_equal = 0
    quality_smsd_better = 0
    quality_rdkit_better = 0

    for i, res in enumerate(results, 1):
        w, spd = winner_and_speedup(res.rdkit_best_ms, res.smsd_best_ms)
        q = mcs_quality(res.rdkit_mcs_size, res.smsd_mcs_size)

        rdkit_str = format_ms(res.rdkit_best_ms)
        smsd_str = format_ms(res.smsd_best_ms)
        rmcs = str(res.rdkit_mcs_size) if res.rdkit_mcs_size >= 0 else "ERR"
        smcs = str(res.smsd_mcs_size) if res.smsd_mcs_size >= 0 else "ERR"

        if res.rdkit_timed_out:
            rdkit_str += "*"

        py_str  = format_ms(res.py_best_ms)
        pymcs   = str(res.py_mcs_size) if res.py_mcs_size >= 0 else "ERR"
        line = (
            f"{i:>2d}  {res.pair_name:<28s}  {res.category:<18s}  "
            f"{rdkit_str:>10s}  {rmcs:>5s}  "
            f"{smsd_str:>10s}  {smcs:>5s}  "
            f"{py_str:>8s}  {pymcs:>5s}  "
            f"{w:<6s}  {spd:>7s}  {q:<14s}"
        )
        print(line)

        if res.rdkit_best_ms > 0:
            rdkit_total += res.rdkit_best_ms
        if res.smsd_best_ms > 0:
            smsd_total += res.smsd_best_ms
        if w == "SMSD":
            smsd_wins += 1
        elif w == "RDKit":
            rdkit_wins += 1
        else:
            ties += 1

        if q == "equal":
            quality_equal += 1
        elif q.startswith("SMSD"):
            quality_smsd_better += 1
        elif q.startswith("RDKit"):
            quality_rdkit_better += 1

    print(sep)
    print()
    print("SUMMARY")
    print(f"  Total RDKit time:  {rdkit_total:>10.1f} ms")
    print(f"  Total SMSD  time:  {smsd_total:>10.1f} ms")
    if smsd_total > 0:
        print(f"  Overall speedup:   {rdkit_total / smsd_total:.2f}x (SMSD vs RDKit)")
    print()
    print(f"  Speed wins:   SMSD={smsd_wins}  RDKit={rdkit_wins}  tie={ties}")
    print(f"  MCS quality:  equal={quality_equal}  SMSD better={quality_smsd_better}  RDKit better={quality_rdkit_better}")
    print()


def write_tsv(results: List[BenchResult], path: Path):
    """Write machine-readable TSV."""
    with open(path, "w") as f:
        cols = [
            "pair_name", "category",
            "rdkit_best_ms", "rdkit_median_ms", "rdkit_mcs_size", "rdkit_timed_out",
            "smsd_best_ms", "smsd_median_ms", "smsd_mcs_size", "smsd_error",
            "cpp_best_ms", "cpp_median_ms", "cpp_mcs_size",
            "winner", "speedup", "mcs_quality",
        ]
        f.write("\t".join(cols) + "\n")

        for res in results:
            w, spd = winner_and_speedup(res.rdkit_best_ms, res.smsd_best_ms)
            q = mcs_quality(res.rdkit_mcs_size, res.smsd_mcs_size)
            row = [
                res.pair_name, res.category,
                f"{res.rdkit_best_ms:.3f}", f"{res.rdkit_median_ms:.3f}",
                str(res.rdkit_mcs_size), str(res.rdkit_timed_out),
                f"{res.smsd_best_ms:.3f}", f"{res.smsd_median_ms:.3f}",
                str(res.smsd_mcs_size), res.smsd_error,
                f"{res.cpp_best_ms:.3f}", f"{res.cpp_median_ms:.3f}",
                str(res.cpp_mcs_size),
                w, spd, q,
            ]
            f.write("\t".join(row) + "\n")


def write_summary(results: List[BenchResult], path: Path):
    """Write formatted summary to file (same as stdout)."""
    import io
    old_stdout = sys.stdout
    sys.stdout = buf = io.StringIO()
    print_results(results)
    sys.stdout = old_stdout
    with open(path, "w") as f:
        f.write(buf.getvalue())


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> int:
    print("[1/5] Running RDKit FindMCS benchmarks...", file=sys.stderr)
    results = benchmark_rdkit(PAIRS)

    print("[2/5] Running SMSD Python (pip) benchmarks...", file=sys.stderr)
    results = benchmark_smsd_python(PAIRS, results)

    print("[3/5] Running SMSD Java CLI benchmarks...", file=sys.stderr)
    results = benchmark_smsd_java(PAIRS, results)

    print("[4/5] Checking for C++ benchmark binary...", file=sys.stderr)
    results = benchmark_cpp_if_available(results)

    print("[5/5] Writing results...", file=sys.stderr)
    tsv_path = SCRIPT_DIR / "results_smsd_vs_rdkit_20pairs_mcs_all.tsv"
    summary_path = SCRIPT_DIR / "results_smsd_vs_rdkit_20pairs_summary.txt"

    write_tsv(results, tsv_path)
    write_summary(results, summary_path)
    print_results(results)

    print(f"\nFiles written:", file=sys.stderr)
    print(f"  {tsv_path}", file=sys.stderr)
    print(f"  {summary_path}", file=sys.stderr)

    return 0


if __name__ == "__main__":
    sys.exit(main())
