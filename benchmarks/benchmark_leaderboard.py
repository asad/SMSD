#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
# Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
# See the NOTICE file for attribution, trademark, and algorithm IP terms.
"""
Unified SMSD benchmark driver.

This script produces a like-for-like comparison across:

* MCS: SMSD Java CLI vs RDKit FindMCS
* Substructure: SMSD Java VF2++ vs RDKit HasSubstructMatch vs CDK DfPattern

The goal is to keep the benchmark organization clean before touching any core
search heuristics. The report also highlights where the expensive paths are:
Java CLI startup, IAtomContainer -> MolGraph conversion, and the cached vs
uncached substructure path.

Outputs:
    - benchmarks/results_leaderboard_mcs.tsv
    - benchmarks/results_leaderboard_substructure.tsv
    - benchmarks/results_leaderboard_summary.txt
"""

from __future__ import annotations

import csv
import json
import os
import subprocess
import sys
import tempfile
import time
from dataclasses import dataclass
from pathlib import Path
from statistics import median
from typing import Dict, Iterable, List, Optional, Tuple


SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_DIR = SCRIPT_DIR.parent

MCS_TIMEOUT_SEC = 10
MCS_RUNS = 5
SUB_WARMUP = 20
SUB_ITERS = 100
SUB_TIMEOUT_MS = 5_000

MCS_RESULTS_TSV = SCRIPT_DIR / "results_leaderboard_mcs.tsv"
SUB_RESULTS_TSV = SCRIPT_DIR / "results_leaderboard_substructure.tsv"
SUMMARY_TXT = SCRIPT_DIR / "results_leaderboard_summary.txt"
CORE_RESULTS_TSV = SCRIPT_DIR / "results_leaderboard_core.tsv"
CORE_PY_RESULTS_TSV = SCRIPT_DIR / "results_leaderboard_core_python.tsv"
CORE_JAVA_SUB_TSV = SCRIPT_DIR / "results_leaderboard_core_java_sub.tsv"
CORE_SUMMARY_TXT = SCRIPT_DIR / "results_leaderboard_core_summary.txt"
SUBSTRUCTURE_PAIRS_FILE = SCRIPT_DIR / "substructure_pairs.tsv"
COMPARE_MODE = "defaults"
PY_WARMUP = None
PY_ITERS = None
PY_TIMEOUT_SEC = None


def find_smsd_jar() -> Optional[Path]:
    """Locate the latest SMSD jar in target/."""
    for base in [PROJECT_DIR / "target", PROJECT_DIR / "src" / "scripts" / "repo"]:
        if not base.exists():
            continue
        jars = list(base.glob("smsd-*-jar-with-dependencies.jar"))
        if jars:
            return max(jars, key=lambda p: (p.stat().st_mtime, p.name))
    return None


def run(cmd: List[str], timeout: Optional[int] = None) -> subprocess.CompletedProcess:
    return subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)


def fmt_time_ms(value: float) -> str:
    if value < 0:
        return "N/A"
    if value < 1.0:
        return f"{value:.3f}"
    if value < 100.0:
        return f"{value:.1f}"
    return f"{value:.0f}"


def fmt_time_us(value: float) -> str:
    if value < 0:
        return "N/A"
    if value < 1_000.0:
        return f"{value:.1f}"
    if value < 1_000_000.0:
        return f"{value / 1_000.0:.2f}ms"
    return f"{value / 1_000_000.0:.1f}s"


def read_substructure_pairs(path: Path) -> List[Tuple[str, str, str, str]]:
    pairs: List[Tuple[str, str, str, str]] = []
    category = ""
    for raw_line in path.read_text().splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#"):
            if line.startswith("# SECTION 1:"):
                category = "Trivial / chain"
            elif line.startswith("# SECTION 2:"):
                category = "Ring containment"
            elif line.startswith("# SECTION 3:"):
                category = "Aromatic queries"
            elif line.startswith("# SECTION 4:"):
                category = "Heterocyclic fragments"
            elif line.startswith("# SECTION 5:"):
                category = "Drug fragment searches"
            elif line.startswith("# SECTION 6:"):
                category = "Edge cases"
            elif line.startswith("# SECTION 7:"):
                category = "Negative controls"
            continue
        cols = line.split("\t", 2)
        if len(cols) != 3:
            continue
        pairs.append((cols[0].strip(), cols[1].strip(), cols[2].strip(), category))
    return pairs


@dataclass
class McsRow:
    pair: str
    category: str
    rdkit_best_ms: float
    rdkit_median_ms: float
    rdkit_mcs: int
    rdkit_timeout: bool
    smsd_best_ms: float
    smsd_median_ms: float
    smsd_mcs: int
    smsd_error: str


@dataclass
class SubRow:
    pair: str
    category: str
    query_smi: str
    target_smi: str
    rdkit_median_us: float = -1.0
    rdkit_hit: bool = False
    smsd_median_us: float = -1.0
    smsd_cached_us: float = -1.0
    cdk_median_us: float = -1.0
    smsd_hit: bool = False
    cdk_hit: bool = False


MCS_PAIRS: List[Tuple[str, str, str, str]] = [
    ("C", "CC", "methane-ethane", "Trivial"),
    ("c1ccccc1", "Cc1ccccc1", "benzene-toluene", "Small aromatic"),
    ("c1ccccc1", "Oc1ccccc1", "benzene-phenol", "Heteroatom"),
    ("CC(=O)Oc1ccccc1C(=O)O", "CC(=O)Nc1ccc(O)cc1", "aspirin-acetaminophen", "Drug pair"),
    ("Cn1cnc2c1c(=O)n(C)c(=O)n2C", "Cn1cnc2c1c(=O)[nH]c(=O)n2C", "caffeine-theophylline", "N-methyl diff"),
    ("CN1CCC23C4C1CC5=C(C2C(C=C4)O3)C=C(C=C5)O",
     "CN1CCC23C4C1CC5=C(C2C(C=C4)OC3)C=C(C=C5)O", "morphine-codeine", "Alkaloid"),
    ("CC(C)Cc1ccc(CC(C)C(=O)O)cc1", "COc1ccc2cc(CC(C)C(=O)O)ccc2c1", "ibuprofen-naproxen", "NSAID"),
    ("C1=NC(=C2C(=N1)N(C=N2)C3C(C(C(O3)COP(=O)(O)OP(=O)(O)OP(=O)(O)O)O)O)N",
     "C1=NC(=C2C(=N1)N(C=N2)C3C(C(C(O3)COP(=O)(O)OP(=O)(O)O)O)O)N", "ATP-ADP", "Nucleotide"),
    ("C1=CC(=C[N+](=C1)C2C(C(C(O2)COP(=O)([O-])OP(=O)(O)OCC3C(C(C(O3)N4C=NC5=C(N=CN=C54)N)O)O)O)O)C(=O)N",
     "C1C=CN(C=C1C(=O)N)C2C(C(C(O2)COP(=O)(O)OP(=O)(O)OCC3C(C(C(O3)N4C=NC5=C4N=CN=C5N)O)O)O)O",
     "NAD-NADH", "Cofactor"),
    ("CC(C)C1=C(C(=C(N1CCC(CC(CC(=O)O)O)O)C2=CC=C(C=C2)F)C3=CC=CC=C3)C(=O)NC4=CC=CC=C4",
     "CC(C)C1=NC(=NC(=C1C=CC(CC(CC(=O)O)O)O)C2=CC=C(C=C2)F)N(C)S(=O)(=O)C",
     "atorvastatin-rosuvastatin", "Statin"),
    ("CC1=C2C(C(=O)C3(C(CC4C(C3C(C(C2(C)C)(CC1OC(=O)C(C(C5=CC=CC=C5)NC(=O)C6=CC=CC=C6)O)O)OC(=O)C7=CC=CC=C7)(CO4)OC(=O)C)O)C)OC(=O)C",
     "CC1=C2C(C(=O)C3(C(CC4C(C3C(C(C2(C)C)(CC1OC(=O)C(C(C5=CC=CC=C5)NC(=O)OC(C)(C)C)O)O)OC(=O)C6=CC=CC=C6)(CO4)OC(=O)C)O)C)O",
     "paclitaxel-docetaxel", "Taxane"),
    ("CCC1C(C(C(C(=O)C(CC(C(C(C(C(C(=O)O1)C)OC2CC(C(C(O2)C)O)(C)OC)C)OC3C(C(CC(O3)C)N(C)C)O)(C)O)C)C)O)(C)O",
     "CCC1C(C(C(N(CC(CC(C(C(C(C(C(=O)O1)C)OC2CC(C(C(O2)C)O)(C)OC)C)OC3C(C(CC(O3)C)N(C)C)O)(C)O)C)C)C)O)(C)O",
     "erythromycin-azithromycin", "Macrolide"),
    ("C1CN2CC3=CCOC4CC(=O)N5C6C4C3CC2C61C7=CC=CC=C75",
     "COC1=CC2=C(C=CN=C2C=C1)C(C3CC4CCN3CC4C=C)O", "strychnine-quinine", "Alkaloid scaffold"),
    ("CC1C(C(CC(O1)OC2C(C(C(OC2OC3=C4C=C5C=C3OC6=C(C=C(C=C6)C(C(C(=O)"
     "NC(C(=O)NC5C(=O)NC7C8=CC(=C(C=C8)O)C9=C(C=C(C=C9O)O)C(NC(=O)C("
     "C(C1=CC(=C(O4)C=C1)Cl)O)NC7=O)C(=O)O)CC(=O)N)NC(=O)C(CC(C)C)NC)"
     "O)Cl)CO)O)O)(C)N)O", "CC1C(C(CC(O1)OC2C(C(C(OC2OC3=C4C=C5C=C3OC6=C(C=C(C=C6)C(C(C(=O)"
     "NC(C(=O)NC5C(=O)NC7C8=CC(=C(C=C8)O)C9=C(C=C(C=C9O)O)C(NC(=O)C("
     "C(C1=CC(=C(O4)C=C1)Cl)O)NC7=O)C(=O)O)CC(=O)N)NC(=O)C(CC(C)C)NC)"
     "O)Cl)CO)O)O)(C)N)O", "vancomycin-self", "Self-match large"),
    ("C1C2CC3CC1CC(C2)C3", "C1C2CC3CC1CC(C2)C3", "adamantane-self", "Symmetric"),
    ("C12C3C4C1C5C4C3C25", "C12C3C4C1C5C4C3C25", "cubane-self", "Cage"),
    ("OCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCO",
     "OCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCO", "PEG12-PEG16", "Polymer"),
    ("c1cc2ccc3ccc4ccc5ccc6ccc1c7c2c3c4c5c67", "c1cc2ccc3ccc4ccc5ccc6ccc1c7c2c3c4c5c67", "coronene-self", "PAH"),
    ("O=c1[nH]c(N)nc2[nH]cnc12", "Oc1nc(N)nc2[nH]cnc12", "guanine-keto-enol", "Tautomer"),
    ("c1cc(c(c(c1)Cl)N2c3cc(cc(c3CNC2=O)c4ccc(cc4F)F)N5CCNCC5)Cl",
     "CCNc1cc(c2c(c1)N(C(=O)NC2)c3ccc(cc3)n4ccc-5ncnc5c4)c6ccnnc6", "rdkit-1585-pair", "Known failure"),
]


def benchmark_rdkit_mcs(smi1: str, smi2: str) -> Tuple[float, float, int, bool]:
    from rdkit import Chem
    from rdkit.Chem import rdFMCS

    mol1 = Chem.MolFromSmiles(smi1)
    mol2 = Chem.MolFromSmiles(smi2)
    if mol1 is None or mol2 is None:
        return -1.0, -1.0, -1, False

    try:
        rdFMCS.FindMCS([mol1, mol2], timeout=1)
    except Exception:
        pass

    times: List[float] = []
    mcs_size = -1
    timed_out = False
    for _ in range(MCS_RUNS):
        t0 = time.perf_counter()
        try:
            result = rdFMCS.FindMCS([mol1, mol2], timeout=MCS_TIMEOUT_SEC)
            mcs_size = result.numAtoms
            timed_out = bool(result.canceled)
        except Exception:
            mcs_size = -1
        times.append((time.perf_counter() - t0) * 1000.0)

    return min(times), median(times), mcs_size, timed_out


def benchmark_smsd_java_mcs(jar: Path, smi1: str, smi2: str) -> Tuple[float, float, int, str]:
    java_bin = "java"
    for candidate in [
        os.environ.get("JAVA_HOME", "") + "/bin/java",
        "/usr/local/opt/openjdk/bin/java",
        "/opt/homebrew/opt/openjdk/bin/java",
        "/usr/bin/java",
    ]:
        if candidate and Path(candidate).exists():
            java_bin = candidate
            break

    cmd = [
        java_bin, "-jar", str(jar),
        "--Q", "SMI", "--q", smi1,
        "--T", "SMI", "--t", smi2,
        "--mode", "mcs",
        "--timeout", "10000",
        "--json", "-",
    ]

    try:
        subprocess.run(
            [java_bin, "-jar", str(jar), "--Q", "SMI", "--q", smi1, "--T", "SMI", "--t", smi2,
             "--mode", "mcs", "--timeout", "10000", "--json", "-"],
            capture_output=True, text=True, timeout=MCS_TIMEOUT_SEC + 30,
        )
    except (FileNotFoundError, subprocess.TimeoutExpired) as exc:
        return -1.0, -1.0, -1, str(exc)

    times: List[float] = []
    mcs_size = -1
    error = ""
    for _ in range(MCS_RUNS):
        t0 = time.perf_counter()
        try:
            proc = subprocess.run(cmd, capture_output=True, text=True, timeout=MCS_TIMEOUT_SEC + 30)
            elapsed = (time.perf_counter() - t0) * 1000.0
        except subprocess.TimeoutExpired:
            elapsed = (time.perf_counter() - t0) * 1000.0
            times.append(elapsed)
            error = "OS timeout"
            continue
        times.append(elapsed)
        if proc.stdout.strip():
            try:
                data = json.loads(proc.stdout)
                mcs_size = int(data.get("mcs_size", data.get("mappingSize", -1)))
            except json.JSONDecodeError:
                error = "JSON parse error"
        elif proc.returncode != 0:
            error = (proc.stderr or "").strip()[:200]

    return min(times), median(times), mcs_size, error


def run_mcs_benchmark() -> List[McsRow]:
    jar = find_smsd_jar()
    if jar is None:
        raise RuntimeError("SMSD jar not found in target/")

    rows: List[McsRow] = []
    for smi1, smi2, pair, category in MCS_PAIRS:
        rd_best, rd_med, rd_mcs, rd_timeout = benchmark_rdkit_mcs(smi1, smi2)
        sm_best, sm_med, sm_mcs, sm_err = benchmark_smsd_java_mcs(jar, smi1, smi2)
        rows.append(McsRow(
            pair=pair,
            category=category,
            rdkit_best_ms=rd_best,
            rdkit_median_ms=rd_med,
            rdkit_mcs=rd_mcs,
            rdkit_timeout=rd_timeout,
            smsd_best_ms=sm_best,
            smsd_median_ms=sm_med,
            smsd_mcs=sm_mcs,
            smsd_error=sm_err,
        ))
    return rows


def load_java_substructure_rows(path: Path) -> Dict[str, SubRow]:
    rows: Dict[str, SubRow] = {}
    with path.open() as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            pair = row["Pair"]
            rows[pair] = SubRow(
                pair=pair,
                category="",
                query_smi="",
                target_smi="",
                smsd_median_us=float(row["SMSD_us"]),
                smsd_cached_us=float(row["SMSD_cached_us"]),
                cdk_median_us=float(row["CDK_us"]),
                smsd_hit=row["SMSD_match"].strip().upper() == "Y",
                cdk_hit=row["CDK_match"].strip().upper() == "Y",
            )
    return rows


def benchmark_rdkit_substructure(pairs: List[Tuple[str, str, str, str]]) -> Dict[str, Tuple[float, bool]]:
    from rdkit import Chem

    rows: Dict[str, Tuple[float, bool]] = {}
    for smi1, smi2, label, _category in pairs:
        q = Chem.MolFromSmiles(smi1)
        t = Chem.MolFromSmiles(smi2)
        if q is None or t is None:
            rows[label] = (-1.0, False)
            continue

        for _ in range(SUB_WARMUP):
            t.HasSubstructMatch(q)

        times: List[float] = []
        hit = False
        for _ in range(SUB_ITERS):
            t0 = time.perf_counter_ns()
            hit = t.HasSubstructMatch(q)
            times.append((time.perf_counter_ns() - t0) / 1000.0)

        rows[label] = (median(times), hit)
    return rows


def run_java_substructure_benchmark(jar: Path) -> Dict[str, SubRow]:
    class_file = SCRIPT_DIR / "benchmark_substructure_java.class"
    if not class_file.exists():
        compile_cmd = [
            "javac",
            "-cp", str(jar),
            str(SCRIPT_DIR / "benchmark_substructure_java.java"),
            "-d", str(SCRIPT_DIR),
        ]
        subprocess.run(compile_cmd, check=True)

    run_cmd = [
        "java",
        "-cp", f"{jar}:{SCRIPT_DIR}",
        "benchmark_substructure_java",
        str(SUBSTRUCTURE_PAIRS_FILE),
    ]
    subprocess.run(run_cmd, check=True)
    return load_java_substructure_rows(SCRIPT_DIR / "results_substructure.tsv")


def build_substructure_rows(java_rows: Dict[str, SubRow]) -> List[SubRow]:
    pairs = read_substructure_pairs(SUBSTRUCTURE_PAIRS_FILE)
    rdkit_rows = benchmark_rdkit_substructure(pairs)

    out: List[SubRow] = []
    for smi1, smi2, label, category in pairs:
        row = java_rows[label]
        row.query_smi = smi1
        row.target_smi = smi2
        row.category = category
        rdkit_med, rdkit_hit = rdkit_rows[label]
        row.rdkit_median_us = rdkit_med
        row.rdkit_hit = rdkit_hit
        out.append(row)
    return out


def parse_python_benchmark_tsv(path: Path) -> Tuple[List[McsRow], List[SubRow]]:
    """Parse the TSV emitted by benchmarks/benchmark_python.py."""
    mcs_rows: List[McsRow] = []
    sub_rows: List[SubRow] = []
    section = None
    header: Optional[List[str]] = None

    for raw_line in path.read_text().splitlines():
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith("# === MCS Benchmark ==="):
            section = "mcs"
            header = None
            continue
        if line.startswith("# === Substructure Benchmark ==="):
            section = "sub"
            header = None
            continue
        if line.startswith("#"):
            continue
        if header is None:
            header = line.split("\t")
            continue
        cols = next(csv.reader([line], delimiter="\t"))
        row = dict(zip(header, cols))
        if section == "mcs":
            mcs_rows.append(McsRow(
                pair=row["pair"],
                category=row["category"],
                rdkit_best_ms=float(row["rdkit_best_us"]) / 1000.0,
                rdkit_median_ms=float(row["rdkit_median_us"]) / 1000.0,
                rdkit_mcs=int(row["rdkit_mcs"]),
                rdkit_timeout=row["rdkit_timeout"].strip().lower() == "true",
                smsd_best_ms=float(row["smsd_best_us"]) / 1000.0,
                smsd_median_ms=float(row["smsd_median_us"]) / 1000.0,
                smsd_mcs=int(row["smsd_mcs"]),
                smsd_error="",
            ))
        elif section == "sub":
            sub_rows.append(SubRow(
                pair=row["pair"],
                category=row["category"],
                query_smi="",
                target_smi="",
                smsd_median_us=float(row["smsd_median_us"]),
                rdkit_median_us=float(row["rdkit_median_us"]),
                smsd_hit=row["smsd_hit"].strip().lower() == "true",
                rdkit_hit=row["rdkit_hit"].strip().lower() == "true",
            ))
    return mcs_rows, sub_rows


def run_core_python_benchmark(warmup: Optional[int] = None,
                              iters: Optional[int] = None,
                              timeout_sec: Optional[int] = None) -> Tuple[List[McsRow], List[SubRow]]:
    """Run the Python benchmark for SMSD vs RDKit under the selected comparison mode."""
    py_script = SCRIPT_DIR / "benchmark_python.py"
    if not py_script.exists():
        raise RuntimeError("benchmark_python.py not found")

    with tempfile.TemporaryDirectory(prefix="smsd-bench-") as tmpdir:
        tsv_path = Path(tmpdir) / "python_benchmark.tsv"
        cmd = ["python3", str(py_script), "--output", str(tsv_path), "--compare-mode", COMPARE_MODE]
        if warmup is not None:
            cmd.extend(["--warmup", str(warmup)])
        if iters is not None:
            cmd.extend(["--iters", str(iters)])
        if timeout_sec is not None:
            cmd.extend(["--timeout-sec", str(timeout_sec)])
        subprocess.run(cmd, check=True)
        CORE_PY_RESULTS_TSV.write_text(tsv_path.read_text())
        return parse_python_benchmark_tsv(tsv_path)


def write_mcs_tsv(rows: List[McsRow], path: Path) -> None:
    with path.open("w") as f:
        f.write("pair\tcategory\trdkit_best_ms\trdkit_median_ms\trdkit_mcs\trdkit_timeout\t")
        f.write("smsd_best_ms\tsmsd_median_ms\tsmsd_mcs\tsmsd_error\twinner\tspeedup\tquality\n")
        for row in rows:
            winner = "N/A"
            speedup = ""
            if row.rdkit_median_ms > 0 and row.smsd_median_ms > 0:
                ratio = row.rdkit_median_ms / row.smsd_median_ms
                if ratio > 1.1:
                    winner, speedup = "SMSD", f"{ratio:.1f}x"
                elif ratio < 0.9:
                    winner, speedup = "RDKit", f"{1.0 / ratio:.1f}x"
                else:
                    winner, speedup = "tie", "~1x"
            quality = "both failed"
            if row.rdkit_mcs >= 0 and row.smsd_mcs >= 0:
                if row.rdkit_mcs == row.smsd_mcs:
                    quality = "equal"
                elif row.smsd_mcs > row.rdkit_mcs:
                    quality = f"SMSD +{row.smsd_mcs - row.rdkit_mcs}"
                else:
                    quality = f"RDKit +{row.rdkit_mcs - row.smsd_mcs}"
            elif row.rdkit_mcs >= 0:
                quality = "RDKit only"
            elif row.smsd_mcs >= 0:
                quality = "SMSD only"
            f.write(
                f"{row.pair}\t{row.category}\t{row.rdkit_best_ms:.3f}\t{row.rdkit_median_ms:.3f}\t"
                f"{row.rdkit_mcs}\t{row.rdkit_timeout}\t{row.smsd_best_ms:.3f}\t{row.smsd_median_ms:.3f}\t"
                f"{row.smsd_mcs}\t{row.smsd_error}\t{winner}\t{speedup}\t{quality}\n"
            )


def write_sub_tsv(rows: List[SubRow], path: Path) -> None:
    with path.open("w") as f:
        f.write("pair\tcategory\tsmsd_us\tsmsd_cached_us\trdkit_us\tcdk_us\tsmsd_hit\trdkit_hit\tcdk_hit\t")
        f.write("cached_speedup_vs_smsd\tcdk_speedup_vs_smsd\tcdk_speedup_vs_rdkit\n")
        for row in rows:
            cached_vs_smsd = row.smsd_median_us / row.smsd_cached_us if row.smsd_cached_us > 0 else 0.0
            cdk_vs_smsd = row.cdk_median_us / row.smsd_median_us if row.smsd_median_us > 0 else 0.0
            cdk_vs_rdkit = row.cdk_median_us / row.rdkit_median_us if row.rdkit_median_us > 0 else 0.0
            f.write(
                f"{row.pair}\t{row.category}\t{row.smsd_median_us:.1f}\t{row.smsd_cached_us:.1f}\t"
                f"{row.rdkit_median_us:.1f}\t{row.cdk_median_us:.1f}\t{row.smsd_hit}\t{row.rdkit_hit}\t{row.cdk_hit}\t"
                f"{cached_vs_smsd:.2f}\t{cdk_vs_smsd:.2f}\t{cdk_vs_rdkit:.2f}\n"
            )


def sub_hotspots(rows: List[SubRow]) -> List[str]:
    ranked = sorted(
        rows,
        key=lambda r: (r.smsd_median_us - r.smsd_cached_us) if r.smsd_median_us > 0 and r.smsd_cached_us > 0 else -1,
        reverse=True,
    )
    out = []
    for row in ranked[:5]:
        gap = row.smsd_median_us - row.smsd_cached_us
        out.append(
            f"{row.pair}: uncached {fmt_time_us(row.smsd_median_us)} us -> cached {fmt_time_us(row.smsd_cached_us)} us "
            f"({gap:.1f} us saved)"
        )
    return out


def mcs_hotspots(rows: List[McsRow]) -> List[str]:
    losers = [r for r in rows if r.smsd_median_ms > 0 and r.rdkit_median_ms > 0 and r.smsd_median_ms > r.rdkit_median_ms]
    losers.sort(key=lambda r: r.smsd_median_ms / r.rdkit_median_ms, reverse=True)
    out = []
    for row in losers[:5]:
        out.append(
            f"{row.pair}: SMSD {fmt_time_ms(row.smsd_median_ms)} ms vs RDKit {fmt_time_ms(row.rdkit_median_ms)} ms"
        )
    return out


def write_summary(mcs_rows: List[McsRow], sub_rows: List[SubRow], path: Path, title: str) -> None:
    mcs_smsd_wins = mcs_rdkit_wins = mcs_ties = 0
    for row in mcs_rows:
        if row.rdkit_median_ms > 0 and row.smsd_median_ms > 0:
            ratio = row.rdkit_median_ms / row.smsd_median_ms
            if ratio > 1.1:
                mcs_smsd_wins += 1
            elif ratio < 0.9:
                mcs_rdkit_wins += 1
            else:
                mcs_ties += 1

    sub_smsd_wins = sub_cdk_wins = sub_ties = 0
    agree = 0
    for row in sub_rows:
        if row.smsd_hit == row.rdkit_hit == row.cdk_hit:
            agree += 1
        if row.cdk_median_us > 0 and row.smsd_median_us > 0:
            ratio = row.cdk_median_us / row.smsd_median_us
            if ratio > 1.1:
                sub_smsd_wins += 1
            elif ratio < 0.9:
                sub_cdk_wins += 1
            else:
                sub_ties += 1

    lines = []
    lines.append(title)
    lines.append(f"Date: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    lines.append(f"Platform: {sys.platform}")
    lines.append("")
    lines.append("MCS")
    lines.append(f"  Pairs: {len(mcs_rows)}")
    lines.append(f"  Speed wins: SMSD={mcs_smsd_wins} RDKit={mcs_rdkit_wins} tie={mcs_ties}")
    lines.append("  Hotspots:")
    for item in mcs_hotspots(mcs_rows):
        lines.append(f"    {item}")
    lines.append("")
    lines.append("Substructure")
    lines.append(f"  Pairs: {len(sub_rows)}")
    lines.append(f"  Speed wins vs CDK: SMSD={sub_smsd_wins} CDK={sub_cdk_wins} tie={sub_ties}")
    lines.append(f"  Hit agreement across SMSD/RDKit/CDK: {agree}/{len(sub_rows)}")
    lines.append("  Hotspots:")
    for item in sub_hotspots(sub_rows):
        lines.append(f"    {item}")

    path.write_text("\n".join(lines) + "\n")


def write_core_summary(mcs_rows: List[McsRow], py_sub_rows: List[SubRow], java_rows: Dict[str, SubRow], path: Path) -> None:
    lines = []
    lines.append("SMSD Core Benchmark")
    lines.append(f"Date: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    lines.append(f"Platform: {sys.platform}")
    lines.append(f"Compare mode: {COMPARE_MODE}")
    if PY_WARMUP is not None or PY_ITERS is not None or PY_TIMEOUT_SEC is not None:
        warmup = "default" if PY_WARMUP is None else str(PY_WARMUP)
        iters = "default" if PY_ITERS is None else str(PY_ITERS)
        timeout = "default" if PY_TIMEOUT_SEC is None else str(PY_TIMEOUT_SEC)
        lines.append(f"Python core protocol: warmup={warmup}, iters={iters}, timeout_sec={timeout}")
    lines.append("")

    lines.append("MCS (Python core vs RDKit)")
    lines.append(f"  Pairs: {len(mcs_rows)}")
    lines.append(f"  Speed wins: SMSD={sum(1 for r in mcs_rows if r.rdkit_median_ms > 0 and r.smsd_median_ms > 0 and (r.rdkit_median_ms / r.smsd_median_ms) > 1.1)} "
                 f"RDKit={sum(1 for r in mcs_rows if r.rdkit_median_ms > 0 and r.smsd_median_ms > 0 and (r.rdkit_median_ms / r.smsd_median_ms) < 0.9)}")
    lines.append("  Hotspots:")
    for item in mcs_hotspots(mcs_rows):
        lines.append(f"    {item}")
    lines.append("")

    lines.append("Substructure (Python core vs RDKit)")
    lines.append(f"  Pairs: {len(py_sub_rows)}")
    lines.append(f"  Hotspots:")
    for row in sorted(py_sub_rows, key=lambda r: r.smsd_median_us - r.rdkit_median_us, reverse=True)[:5]:
        lines.append(
            f"    {row.pair}: SMSD {fmt_time_us(row.smsd_median_us)} us vs RDKit {fmt_time_us(row.rdkit_median_us)} us"
        )
    lines.append("")

    java_sub_rows = list(java_rows.values())
    lines.append("Substructure (Java cached path vs CDK)")
    lines.append(f"  Pairs: {len(java_sub_rows)}")
    lines.append(f"  Speed wins vs CDK: SMSD={sum(1 for r in java_sub_rows if r.cdk_median_us > 0 and r.smsd_cached_us > 0 and (r.cdk_median_us / r.smsd_cached_us) > 1.1)} "
                 f"CDK={sum(1 for r in java_sub_rows if r.cdk_median_us > 0 and r.smsd_cached_us > 0 and (r.cdk_median_us / r.smsd_cached_us) < 0.9)}")
    lines.append(f"  Hit agreement: {sum(1 for r in java_sub_rows if r.smsd_hit == r.cdk_hit)}/{len(java_sub_rows)}")
    lines.append("  Hotspots:")
    for item in sub_hotspots(java_sub_rows):
        lines.append(f"    {item}")

    path.write_text("\n".join(lines) + "\n")


def main() -> int:
    import argparse

    parser = argparse.ArgumentParser(description="Unified SMSD benchmark runner")
    parser.add_argument(
        "--mode",
        choices=["integration", "core", "profile"],
        default="integration",
        help="integration = CLI/CDK/RDKit leaderboard, core = in-process core+RDKit, profile = summary only",
    )
    parser.add_argument("--mcs-runs", type=int,
                        help="Override the number of MCS runs for integration-mode RDKit/SMSD-Java timing")
    parser.add_argument("--sub-warmup", type=int,
                        help="Override the number of substructure warmup runs")
    parser.add_argument("--sub-iters", type=int,
                        help="Override the number of measured substructure runs")
    parser.add_argument("--py-warmup", type=int,
                        help="Override the warmup count for benchmark_python.py in core mode")
    parser.add_argument("--py-iters", type=int,
                        help="Override the measured-run count for benchmark_python.py in core mode")
    parser.add_argument("--py-timeout-sec", type=int,
                        help="Override the per-pair MCS timeout passed to benchmark_python.py in core mode")
    parser.add_argument("--compare-mode", choices=["defaults", "strict", "fmcs"],
                        default="defaults",
                        help="defaults = toolkit defaults, strict = exact/ring-parity baseline, fmcs = loose FMCS-style baseline")
    args = parser.parse_args()

    global MCS_RUNS, SUB_WARMUP, SUB_ITERS, COMPARE_MODE, PY_WARMUP, PY_ITERS, PY_TIMEOUT_SEC
    if args.mcs_runs is not None:
        MCS_RUNS = args.mcs_runs
    if args.sub_warmup is not None:
        SUB_WARMUP = args.sub_warmup
    if args.sub_iters is not None:
        SUB_ITERS = args.sub_iters
    COMPARE_MODE = args.compare_mode
    PY_WARMUP = args.py_warmup
    PY_ITERS = args.py_iters
    PY_TIMEOUT_SEC = args.py_timeout_sec

    from rdkit import Chem  # noqa: F401 - fail fast if RDKit is missing

    if args.mode == "integration":
        jar = find_smsd_jar()
        if jar is None:
            print("ERROR: SMSD jar not found in target/", file=sys.stderr)
            return 1

        print(f"[1/3] Running MCS benchmark with {jar.name}...", file=sys.stderr)
        mcs_rows = run_mcs_benchmark()
        write_mcs_tsv(mcs_rows, MCS_RESULTS_TSV)

        print(f"[2/3] Running Java substructure benchmark with {jar.name}...", file=sys.stderr)
        java_rows = run_java_substructure_benchmark(jar)
        sub_rows = build_substructure_rows(java_rows)
        write_sub_tsv(sub_rows, SUB_RESULTS_TSV)

        print("[3/3] Writing summary...", file=sys.stderr)
        write_summary(mcs_rows, sub_rows, SUMMARY_TXT, "SMSD Integration Leaderboard Benchmark")

        print(f"Results written to {MCS_RESULTS_TSV}")
        print(f"Results written to {SUB_RESULTS_TSV}")
        print(f"Summary written to {SUMMARY_TXT}")
        return 0

    if args.mode == "core":
        print("[1/2] Running in-process Python core benchmark...", file=sys.stderr)
        mcs_rows, py_sub_rows = run_core_python_benchmark(
            warmup=args.py_warmup,
            iters=args.py_iters,
            timeout_sec=args.py_timeout_sec,
        )

        print("[2/2] Running Java cached substructure benchmark...", file=sys.stderr)
        jar = find_smsd_jar()
        if jar is None:
            print("ERROR: SMSD jar not found in target/", file=sys.stderr)
            return 1
        java_rows = run_java_substructure_benchmark(jar)
        sub_rows = build_substructure_rows(java_rows)
        write_mcs_tsv(mcs_rows, CORE_RESULTS_TSV.with_name("results_leaderboard_core_mcs.tsv"))
        write_sub_tsv(sub_rows, CORE_JAVA_SUB_TSV)
        write_core_summary(mcs_rows, py_sub_rows, java_rows, CORE_SUMMARY_TXT)
        print(f"Results written to {CORE_PY_RESULTS_TSV}")
        print(f"Results written to {CORE_JAVA_SUB_TSV}")
        print(f"Summary written to {CORE_SUMMARY_TXT}")
        return 0

    print("[profile] Reading existing benchmark outputs...", file=sys.stderr)
    if MCS_RESULTS_TSV.exists() and SUB_RESULTS_TSV.exists():
        print(SUMMARY_TXT.read_text())
        return 0
    print("ERROR: no benchmark outputs found to profile", file=sys.stderr)
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
