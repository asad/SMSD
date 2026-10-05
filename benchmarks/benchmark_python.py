#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
# Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
# See the NOTICE file for attribution, trademark, and algorithm IP terms.
"""
SMSD Python (C++ pybind11) and RDKit — MCS and substructure measurements

Both tools are called from the same Python process with identical SMILES
inputs. Chemistry settings, validity, quality and cancellation are reported
separately; policy differences do not receive speedup claims.

Usage:
    pip install smsd rdkit
    python3 benchmarks/benchmark_python.py

Output:
    - build/local-benchmarks/20pairs.tsv and .json (machine-readable)
    - stdout: formatted comparison table

Sections:
    Section 1: MCS benchmark  — SMSD find_mcs vs RDKit FindMCS (20 pairs)
    Section 2: Substructure   — SMSD is_substructure vs RDKit HasSubstructMatch (20 pairs)

Reference datasets (molecule pairs):
    20 curated pairs spanning trivial → glycopeptide (vancomycin, 101 atoms);
    includes self-match, tautomer, known-hard, and macrolide cases.
"""

import json
import platform
import statistics
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import List, Tuple

from mcs_protocol import chem_options, rdkit_parameters, measure_pair, summarize

# ---------------------------------------------------------------------------
# Imports
# ---------------------------------------------------------------------------

try:
    import smsd
except ImportError:
    sys.exit("ERROR: smsd not installed. Run: pip install smsd")

try:
    import smsd._smsd as smsd_native
except ImportError:
    sys.exit("ERROR: smsd native extension missing. Reinstall smsd from source.")

try:
    from rdkit import Chem
    from rdkit.Chem import rdFMCS
except ImportError:
    sys.exit("ERROR: rdkit not installed. Run: pip install rdkit")


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

# Default protocol: 1 warmup + 3 measured; report all observations and median.
WARMUP = 1
ITERS = 3
TIMEOUT_SEC = 10
COMPARE_MODE = "defaults"

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent

# ---------------------------------------------------------------------------
# 20 Molecule Pairs (identical to benchmark_all.py)
# ---------------------------------------------------------------------------

VANCOMYCIN = (
    "CC1C(C(CC(O1)OC2C(C(C(OC2OC3=C4C=C5C=C3OC6=C(C=C(C=C6)C(C(C(=O)"
    "NC(C(=O)NC5C(=O)NC7C8=CC(=C(C=C8)O)C9=C(C=C(C=C9O)O)C(NC(=O)C("
    "C(C1=CC(=C(O4)C=C1)Cl)O)NC7=O)C(=O)O)CC(=O)N)NC(=O)C(CC(C)C)NC)"
    "O)Cl)CO)O)O)(C)N)O"
)
PEG16 = "OCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCO"

PAIRS: List[Tuple[str, str, str, str]] = [
    # (smi1, smi2, name, category)
    ("C", "CC", "methane-ethane", "Trivial"),
    ("c1ccccc1", "Cc1ccccc1", "benzene-toluene", "Small aromatic"),
    ("c1ccccc1", "Oc1ccccc1", "benzene-phenol", "Heteroatom"),
    ("CC(=O)Oc1ccccc1C(=O)O", "CC(=O)Nc1ccc(O)cc1",
     "aspirin-acetaminophen", "Drug pair"),
    ("Cn1cnc2c1c(=O)n(C)c(=O)n2C", "Cn1cnc2c1c(=O)[nH]c(=O)n2C",
     "caffeine-theophylline", "N-methyl diff"),
    ("CN1CCC23C4C1CC5=C(C2C(C=C4)O3)C=C(C=C5)O",
     "CN1CCC23C4C1CC5=C(C2C(C=C4)OC3)C=C(C=C5)O",
     "morphine-codeine", "Alkaloid"),
    ("CC(C)Cc1ccc(CC(C)C(=O)O)cc1",
     "COc1ccc2cc(CC(C)C(=O)O)ccc2c1",
     "ibuprofen-naproxen", "NSAID"),
    ("C1=NC(=C2C(=N1)N(C=N2)C3C(C(C(O3)COP(=O)(O)OP(=O)(O)O)O)O)N",
     "C1=NC(=C2C(=N1)N(C=N2)C3C(C(C(O3)COP(=O)(O)OP(=O)(O)OP(=O)(O)O)O)O)N",
     "ATP-ADP", "Nucleotide"),
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
     "COC1=CC2=C(C=CN=C2C=C1)C(C3CC4CCN3CC4C=C)O",
     "strychnine-quinine", "Alkaloid scaffold"),
    (VANCOMYCIN, VANCOMYCIN, "vancomycin-self", "Self-match large"),
    ("C1C2CC3CC1CC(C2)C3", "C1C2CC3CC1CC(C2)C3",
     "adamantane-self", "Symmetric"),
    ("C12C3C4C1C5C4C3C25", "C12C3C4C1C5C4C3C25",
     "cubane-self", "Cage"),
    ("OCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCO", PEG16,
     "PEG12-PEG16", "Polymer"),
    ("c1cc2ccc3ccc4ccc5ccc6ccc1c7c2c3c4c5c67",
     "c1cc2ccc3ccc4ccc5ccc6ccc1c7c2c3c4c5c67",
     "coronene-self", "PAH"),
    ("O=c1[nH]c(N)nc2[nH]cnc12", "Oc1nc(N)nc2[nH]cnc12",
     "guanine-keto-enol", "Tautomer"),
    ("c1cc(c(c(c1)Cl)N2c3cc(cc(c3CNC2=O)c4ccc(cc4F)F)N5CCNCC5)Cl",
     "CCNc1cc(c2c(c1)N(C(=O)NC2)c3ccc(cc3)n4ccc-5ncnc5c4)c6ccnnc6",
     "rdkit-1585-pair", "Known failure"),
]


# ---------------------------------------------------------------------------
# Data class
# ---------------------------------------------------------------------------

@dataclass
class MCSResult:
    name: str
    category: str
    smsd_best_us: float  # microseconds
    smsd_median_us: float
    smsd_mcs: int
    rdkit_best_us: float
    rdkit_median_us: float
    rdkit_mcs: int
    rdkit_timed_out: bool


# Back-compat alias so existing TSV-writing code still works
Result = MCSResult


@dataclass
class SubResult:
    name: str
    category: str
    smsd_median_us: float
    smsd_hit: bool
    rdkit_median_us: float
    rdkit_hit: bool


# ---------------------------------------------------------------------------
# Benchmark functions
# ---------------------------------------------------------------------------

def make_smsd_chem_options():
    return chem_options("fmcs" if COMPARE_MODE == "defaults" else COMPARE_MODE)


def make_rdkit_mcs_call():
    params = rdkit_parameters("fmcs" if COMPARE_MODE == "defaults" else COMPARE_MODE, TIMEOUT_SEC)
    return lambda mols: rdFMCS.FindMCS(mols, params)

def bench_smsd(smi1: str, smi2: str) -> Tuple[List[float], int]:
    """Benchmark SMSD Python MCS. Returns (times_us, mcs_size)."""
    g1 = smsd.parse_smiles(smi1)
    g2 = smsd.parse_smiles(smi2)
    opts = make_smsd_chem_options()
    mcs_opts = smsd.MCSOptions()
    mcs_opts.timeout_ms = TIMEOUT_SEC * 1000

    # Warmup
    for _ in range(WARMUP):
        smsd_native.find_mcs(g1, g2, opts, mcs_opts)

    # Timed
    times = []
    mcs_size = 0
    for _ in range(ITERS):
        t0 = time.perf_counter_ns()
        mapping = smsd_native.find_mcs(g1, g2, opts, mcs_opts)
        dt = (time.perf_counter_ns() - t0) / 1000.0  # ns -> us
        times.append(dt)
        mcs_size = len(mapping)

    return times, mcs_size


def bench_rdkit(smi1: str, smi2: str) -> Tuple[List[float], int, bool]:
    """Benchmark RDKit FindMCS. Returns (times_us, mcs_size, timed_out)."""
    mol1 = Chem.MolFromSmiles(smi1)
    mol2 = Chem.MolFromSmiles(smi2)
    if mol1 is None or mol2 is None:
        return [float("inf")] * ITERS, -1, False
    find_mcs = make_rdkit_mcs_call()

    # Warmup
    for _ in range(WARMUP):
        try:
            find_mcs([mol1, mol2])
        except Exception:
            pass

    # Timed
    times = []
    mcs_size = 0
    timed_out = False
    for _ in range(ITERS):
        t0 = time.perf_counter_ns()
        try:
            result = find_mcs([mol1, mol2])
            mcs_size = result.numAtoms
            if result.canceled:
                timed_out = True
        except Exception:
            mcs_size = -1
        dt = (time.perf_counter_ns() - t0) / 1000.0
        times.append(dt)

    return times, mcs_size, timed_out


def bench_smsd_sub(query_smi: str, target_smi: str) -> Tuple[List[float], bool]:
    """Benchmark SMSD Python substructure search.

    Returns (times_us, is_match) where is_match is True if query is found
    as a substructure in target.
    """
    g_query = smsd.parse_smiles(query_smi)
    g_target = smsd.parse_smiles(target_smi)
    opts = make_smsd_chem_options()

    for _ in range(WARMUP):
        smsd_native.is_substructure(g_query, g_target, opts, TIMEOUT_SEC * 1000)

    times = []
    match = False
    for _ in range(ITERS):
        t0 = time.perf_counter_ns()
        match = smsd_native.is_substructure(g_query, g_target, opts, TIMEOUT_SEC * 1000)
        dt = (time.perf_counter_ns() - t0) / 1000.0
        times.append(dt)

    return times, match


def bench_rdkit_sub(query_smi: str, target_smi: str) -> Tuple[List[float], bool]:
    """Benchmark RDKit HasSubstructMatch.

    Returns (times_us, is_match) where is_match is True if query is found
    as a substructure in target.
    """
    mol_query = Chem.MolFromSmiles(query_smi)
    mol_target = Chem.MolFromSmiles(target_smi)
    if mol_query is None or mol_target is None:
        return [float("inf")] * ITERS, False

    for _ in range(WARMUP):
        mol_target.HasSubstructMatch(mol_query)

    times = []
    match = False
    for _ in range(ITERS):
        t0 = time.perf_counter_ns()
        match = mol_target.HasSubstructMatch(mol_query)
        dt = (time.perf_counter_ns() - t0) / 1000.0
        times.append(dt)

    return times, match


def get_smsd_paths() -> Tuple[Path, Path]:
    """Return resolved paths for the Python package and native extension."""
    return Path(smsd.__file__).resolve(), Path(smsd_native.__file__).resolve()


def is_local_smsd_import() -> bool:
    """True if both the Python package and native extension resolve inside this repo."""
    pkg_path, native_path = get_smsd_paths()
    return pkg_path.is_relative_to(REPO_ROOT) and native_path.is_relative_to(REPO_ROOT)


def write_results_tsv(tsv_path: Path, results: List[Result], sub_results: List[SubResult]) -> None:
    """Write whichever benchmark sections were executed to a TSV file."""
    tsv_path.parent.mkdir(parents=True, exist_ok=True)
    with open(tsv_path, "w") as f:
        f.write(f"# compare_mode={COMPARE_MODE}\n")
        if results:
            f.write("# === MCS Benchmark ===\n")
            f.write("pair\tcategory\tsmsd_best_us\tsmsd_median_us\tsmsd_mcs\t"
                    "rdkit_best_us\trdkit_median_us\trdkit_mcs\trdkit_timeout\n")
            for r in results:
                f.write(f"{r.name}\t{r.category}\t"
                        f"{r.smsd_best_us:.1f}\t{r.smsd_median_us:.1f}\t{r.smsd_mcs}\t"
                        f"{r.rdkit_best_us:.1f}\t{r.rdkit_median_us:.1f}\t{r.rdkit_mcs}\t"
                        f"{r.rdkit_timed_out}\n")
            f.write("\n")

        if sub_results:
            f.write("# === Substructure Benchmark ===\n")
            f.write("pair\tcategory\tsmsd_median_us\tsmsd_hit\trdkit_median_us\trdkit_hit\n")
            for r in sub_results:
                f.write(f"{r.name}\t{r.category}\t"
                        f"{r.smsd_median_us:.1f}\t{r.smsd_hit}\t"
                        f"{r.rdkit_median_us:.1f}\t{r.rdkit_hit}\n")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main(sub_only=False, mcs_only=False, output_path=None,
         print_smsd_path=False, require_local_smsd=False,
         warmup=None, iters=None, timeout_sec=None, compare_mode=None):
    global WARMUP, ITERS, TIMEOUT_SEC, COMPARE_MODE
    if warmup is not None: WARMUP = warmup
    if iters is not None: ITERS = iters
    if timeout_sec is not None: TIMEOUT_SEC = timeout_sec
    if compare_mode is not None: COMPARE_MODE = compare_mode
    if WARMUP < 0 or ITERS < 1 or TIMEOUT_SEC < 1:
        raise ValueError("warmup must be nonnegative; iterations/timeout positive")
    pkg_path, native_path = get_smsd_paths()
    if require_local_smsd and not is_local_smsd_import():
        raise RuntimeError(f"SMSD package/extension outside repository: {pkg_path}, {native_path}")
    print(f"SMSD {smsd.__version__}, RDKit {Chem.rdBase.rdkitVersion}, Python {platform.python_version()}")
    print(f"Protocol: {WARMUP} warmup + {ITERS} measured, timeout {TIMEOUT_SEC}s, atom objective, connected")
    print(f"SMSD package: {pkg_path}; extension: {native_path}; {smsd_native.gpu_device_info()}")
    policy = "fmcs" if COMPARE_MODE == "defaults" else COMPARE_MODE
    print(f"Policy: {policy}. Strict RDKit atom comparison uses Python callbacks; FMCS chemistry differs.")
    print("Search-only medians. Invalid, canceled, unequal-quality or different-policy rows have no speed ratio.")
    results, sub_results, detailed = [], [], []
    if not sub_only:
        for index, (smi1, smi2, name, category) in enumerate(PAIRS):
            detail = {"pair": name, "category": category, "policy": policy}
            try:
                rows = measure_pair(smi1, smi2, policy, TIMEOUT_SEC, WARMUP, ITERS, index)
                detail.update(summarize(rows, policy)); detail["observations"] = rows
                values = {engine: [r for r in rows if r["engine"] == engine] for engine in ("smsd", "rdkit")}
                results.append(Result(name, category,
                    min(r["elapsed_us"] for r in values["smsd"]), detail["smsd_median_us"],
                    min(detail["smsd_atoms"]), min(r["elapsed_us"] for r in values["rdkit"]),
                    detail["rdkit_median_us"], min(detail["rdkit_atoms"]), detail["rdkit_timeouts"] > 0))
                print(f"{name}: SMSD {detail['smsd_median_us']:.1f}us {detail['smsd_atoms']} atoms; "
                      f"RDKit {detail['rdkit_median_us']:.1f}us {detail['rdkit_atoms']} atoms; "
                      f"valid {detail['smsd_valid']}/{detail['rdkit_valid']}; comparable={detail['speed_comparable']}", flush=True)
            except Exception as exc:
                detail["error"] = str(exc)
                print(f"{name}: input/search error: {exc}", flush=True)
                results.append(Result(name, category, float("inf"), float("inf"), -1,
                                      float("inf"), float("inf"), -1, False))
            detailed.append(detail)
    if not mcs_only:
        for smi1, smi2, name, category in PAIRS:
            try:
                st, sh = bench_smsd_sub(smi1, smi2)
                rt, rh = bench_rdkit_sub(smi1, smi2)
                sub_results.append(SubResult(name, category, statistics.median(st), sh,
                                             statistics.median(rt), rh))
                print(f"substructure {name}: hit agreement={sh == rh}", flush=True)
            except Exception as exc:
                detailed.append({"pair": name, "operation": "substructure", "error": str(exc)})
    path = Path(output_path) if output_path else REPO_ROOT / "build/local-benchmarks/20pairs.tsv"
    write_results_tsv(path, results, sub_results)
    path.with_suffix(".json").write_text(json.dumps({"smsd_version": smsd.__version__,
        "smsd_backend": smsd_native.gpu_device_info(), "smsd_package": str(pkg_path),
        "smsd_extension": str(native_path), "rdkit_version": Chem.rdBase.rdkitVersion,
        "timeout_sec": TIMEOUT_SEC, "warmup": WARMUP, "iterations": ITERS,
        "smsd_timeout_status": "unknown: mapping API does not expose cancellation",
        "results": detailed}, indent=2)+"\n")
    print(f"Results: {path}")


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="SMSD vs RDKit benchmark")
    parser.add_argument("--sub-only", action="store_true",
                        help="Run only the substructure benchmark (skip MCS)")
    parser.add_argument("--mcs-only", action="store_true",
                        help="Run only the MCS benchmark (skip substructure)")
    parser.add_argument("--output", type=Path,
                        help="Write TSV output to this path instead of overwriting the default file")
    parser.add_argument("--print-smsd-path", action="store_true",
                        help="Print the resolved smsd package and native extension paths")
    parser.add_argument("--require-local-smsd", action="store_true",
                        help="Fail unless both smsd and smsd._smsd resolve inside this repo")
    parser.add_argument("--warmup", type=int,
                        help="Override the number of warmup runs")
    parser.add_argument("--iters", type=int,
                        help="Override the number of measured runs")
    parser.add_argument("--timeout-sec", type=int,
                        help="Override the per-pair MCS timeout in seconds")
    parser.add_argument("--compare-mode", choices=["defaults", "strict", "fmcs", "any"],
                        help="defaults = fmcs; strict = exact/aromatic/ring/charge, fmcs = toolkit chemistry comparison, any = common bond-any policy")
    args = parser.parse_args()
    main(sub_only=args.sub_only,
         mcs_only=args.mcs_only,
         output_path=args.output,
         print_smsd_path=args.print_smsd_path,
         require_local_smsd=args.require_local_smsd,
         warmup=args.warmup,
         iters=args.iters,
         timeout_sec=args.timeout_sec,
         compare_mode=args.compare_mode)
