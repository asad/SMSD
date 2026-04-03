#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
# Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
# See the NOTICE file for attribution, trademark, and algorithm IP terms.
"""
SMSD Python (C++ pybind11) vs RDKit — Fair Python-to-Python MCS + Substructure Benchmark

Both tools are called from the same Python process with identical SMILES
inputs and matching parameters. This eliminates cross-language overhead
(no subprocess, no JVM startup) and gives a like-for-like comparison.

Usage:
    pip install smsd rdkit
    python3 benchmarks/benchmark_python.py

Output:
    - benchmarks/results_smsd_vs_rdkit_20pairs_mcs_sub.tsv   (machine-readable)
    - stdout: formatted comparison table

Sections:
    Section 1: MCS benchmark  — SMSD find_mcs vs RDKit FindMCS (20 pairs)
    Section 2: Substructure   — SMSD is_substructure vs RDKit HasSubstructMatch (20 pairs)

Reference datasets (molecule pairs):
    20 curated pairs spanning trivial → glycopeptide (vancomycin, 101 atoms);
    includes self-match, tautomer, known-hard, and macrolide cases.
"""

import platform
import statistics
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import List, Tuple

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

# Protocol: 3 warmup + 10 measured, report median.
# Rationale: 10 runs exceeds the McCreesh/Glasgow (2017) standard of 5 while
# remaining practical for large-molecule pairs (paclitaxel, vancomycin).
# Median across 10 runs gives <0.5% coefficient of variation for timings >=10us.
WARMUP = 3
ITERS = 10
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
class McsResult:
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
Result = McsResult


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
    """Return ChemOptions for the selected comparison mode."""
    if COMPARE_MODE == "strict":
        return smsd.ChemOptions.profile("strict")
    if COMPARE_MODE == "fmcs":
        opts = smsd.ChemOptions()
        opts.match_bond_order = smsd.BondOrderMode.LOOSE
        opts.aromaticity_mode = smsd.AromaticityMode.FLEXIBLE
        opts.ring_matches_ring_only = False
        opts.complete_rings_only = False
        opts.match_formal_charge = False
        return opts
    return smsd.ChemOptions()


def make_rdkit_mcs_call():
    """Return an RDKit FindMCS callable for the selected comparison mode."""
    if COMPARE_MODE == "defaults":
        return lambda mols: rdFMCS.FindMCS(mols, timeout=TIMEOUT_SEC)

    params = rdFMCS.MCSParameters()
    params.Timeout = TIMEOUT_SEC
    params.AtomTyper = rdFMCS.AtomCompare.CompareElements
    atom = params.AtomCompareParameters
    bond = params.BondCompareParameters
    atom.MatchChiralTag = False
    atom.MatchIsotope = False
    atom.MatchValences = False
    bond.MatchStereo = False
    bond.MatchFusedRings = False
    bond.MatchFusedRingsStrict = False

    if COMPARE_MODE == "strict":
        params.BondTyper = rdFMCS.BondCompare.CompareOrderExact
        atom.MatchFormalCharge = True
        atom.RingMatchesRingOnly = True
        bond.RingMatchesRingOnly = True
        atom.CompleteRingsOnly = False
        bond.CompleteRingsOnly = False
    elif COMPARE_MODE == "fmcs":
        params.BondTyper = rdFMCS.BondCompare.CompareOrder
        atom.MatchFormalCharge = False
        atom.RingMatchesRingOnly = False
        bond.RingMatchesRingOnly = False
        atom.CompleteRingsOnly = False
        bond.CompleteRingsOnly = False
    else:
        raise ValueError(f"Unsupported compare mode: {COMPARE_MODE}")

    return lambda mols: rdFMCS.FindMCS(mols, params)

def bench_smsd(smi1: str, smi2: str) -> Tuple[List[float], int]:
    """Benchmark SMSD Python MCS. Returns (times_us, mcs_size)."""
    g1 = smsd.parse_smiles(smi1)
    g2 = smsd.parse_smiles(smi2)
    opts = make_smsd_chem_options()
    mcs_opts = smsd.McsOptions()
    mcs_opts.timeout_ms = TIMEOUT_SEC * 1000

    # Warmup
    for _ in range(WARMUP):
        smsd.find_mcs(g1, g2, opts, mcs_opts)

    # Timed
    times = []
    mcs_size = 0
    for _ in range(ITERS):
        t0 = time.perf_counter_ns()
        mapping = smsd.find_mcs(g1, g2, opts, mcs_opts)
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
        smsd.is_substructure(g_query, g_target, opts)

    times = []
    match = False
    for _ in range(ITERS):
        t0 = time.perf_counter_ns()
        match = smsd.is_substructure(g_query, g_target, opts)
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
    if warmup is not None:
        WARMUP = warmup
    if iters is not None:
        ITERS = iters
    if timeout_sec is not None:
        TIMEOUT_SEC = timeout_sec
    if compare_mode is not None:
        COMPARE_MODE = compare_mode

    pkg_path, native_path = get_smsd_paths()
    local_smsd = is_local_smsd_import()
    if require_local_smsd and not local_smsd:
        sys.exit(
            "ERROR: smsd did not resolve to this repo.\n"
            f"  package: {pkg_path}\n"
            f"  native:  {native_path}\n"
            "Reinstall the local tree first, for example: python3 -m pip install -e ."
        )

    print("=" * 120)
    print("SMSD Python (C++ pybind11) vs RDKit FindMCS")
    print(f"Platform: {platform.system()} {platform.machine()}")
    print(f"Python:   {platform.python_version()}")
    print(f"SMSD:     {smsd.__version__}")
    print(f"RDKit:    {Chem.rdBase.rdkitVersion}")
    print(f"Protocol: {WARMUP} warmup + {ITERS} measured, median reported")
    print(f"Timeout:  {TIMEOUT_SEC}s per pair")
    print(f"Mode:     {COMPARE_MODE}")
    if print_smsd_path or require_local_smsd:
        print(f"SMSD pkg: {pkg_path}")
        print(f"SMSD ext: {native_path}")
        print(f"Local:    {'yes' if local_smsd else 'no'}")
    print("=" * 120)
    print()

    results: List[Result] = []
    sub_results: List[SubResult] = []
    smsd_wins = rdkit_wins = ties = 0
    mcs_equal = mcs_smsd_better = mcs_rdkit_better = 0

    if not sub_only:
        hdr = (f" # {'Pair':<30s} {'Category':<20s} "
               f"{'SMSD(us)':>10s} {'RDKit(us)':>10s} "
               f"{'S.MCS':>5s} {'R.MCS':>5s} "
               f"{'Speedup':>10s} {'MCS Quality':<15s}")
        print(hdr)
        print("-" * 120)

        for idx, (smi1, smi2, name, category) in enumerate(PAIRS, 1):
            # --- SMSD ---
            try:
                s_times, s_mcs = bench_smsd(smi1, smi2)
                s_best = min(s_times)
                s_median = statistics.median(s_times)
            except Exception as e:
                s_best = s_median = float("inf")
                s_mcs = -1

            # --- RDKit ---
            try:
                r_times, r_mcs, r_timeout = bench_rdkit(smi1, smi2)
                r_best = min(r_times)
                r_median = statistics.median(r_times)
            except Exception as e:
                r_best = r_median = float("inf")
                r_mcs = -1
                r_timeout = False

            # --- Speedup ---
            if r_timeout:
                speedup_str = "timeout"
                smsd_wins += 1
            elif s_median == 0 or r_median == 0:
                speedup_str = "N/A"
                ties += 1
            else:
                ratio = r_median / s_median
                if ratio > 1.1:
                    speedup_str = f"SMSD {ratio:.1f}x"
                    smsd_wins += 1
                elif ratio < 0.9:
                    speedup_str = f"RDKit {1/ratio:.1f}x"
                    rdkit_wins += 1
                else:
                    speedup_str = "~tie"
                    ties += 1

            # --- MCS quality ---
            if s_mcs > 0 and r_mcs > 0:
                if s_mcs > r_mcs:
                    quality = f"SMSD +{s_mcs - r_mcs}"
                    mcs_smsd_better += 1
                elif r_mcs > s_mcs:
                    quality = f"RDKit +{r_mcs - s_mcs}"
                    mcs_rdkit_better += 1
                else:
                    quality = "equal"
                    mcs_equal += 1
            elif s_mcs > 0 and (r_mcs <= 0 or r_timeout):
                quality = f"SMSD only ({s_mcs})"
                mcs_smsd_better += 1
            elif r_mcs > 0 and s_mcs <= 0:
                quality = f"RDKit only ({r_mcs})"
                mcs_rdkit_better += 1
            else:
                quality = "both failed"

            # --- Format times ---
            def fmt_us(val):
                if val == float("inf"):
                    return "ERR"
                if val >= 1_000_000:
                    return f"{val/1_000_000:.1f}s"
                if val >= 1_000:
                    return f"{val/1_000:.2f}ms"
                return f"{val:.0f}us"

            r_display = "timeout" if r_timeout else fmt_us(r_median)

            print(f"{idx:2d} {name:<30s} {category:<20s} "
                  f"{fmt_us(s_median):>10s} {r_display:>10s} "
                  f"{s_mcs:>5d} {r_mcs:>5d} "
                  f"{speedup_str:>10s} {quality:<15s}")

            results.append(Result(
                name=name, category=category,
                smsd_best_us=s_best, smsd_median_us=s_median, smsd_mcs=s_mcs,
                rdkit_best_us=r_best, rdkit_median_us=r_median, rdkit_mcs=r_mcs,
                rdkit_timed_out=r_timeout,
            ))

        print("-" * 120)
        print()
        print("MCS SUMMARY")
        print(f"  Speed wins:  SMSD={smsd_wins}  RDKit={rdkit_wins}  tie={ties}")
        print(f"  MCS quality: equal={mcs_equal}  SMSD better={mcs_smsd_better}  RDKit better={mcs_rdkit_better}")

    if mcs_only:
        tsv_path = Path(output_path) if output_path else SCRIPT_DIR / "results_smsd_vs_rdkit_20pairs_mcs_sub.tsv"
        write_results_tsv(tsv_path, results, sub_results)
        print(f"\n  Results written to {tsv_path}")
        return

    # =========================================================================
    # Substructure benchmark: SMSD is_substructure vs RDKit HasSubstructMatch
    # =========================================================================
    print()
    print("=" * 120)
    print("SUBSTRUCTURE SEARCH — SMSD Python (C++ pybind11) vs RDKit HasSubstructMatch")
    print(f"Protocol: {WARMUP} warmup + {ITERS} measured, median reported")
    print(f"Query = smi1 searched within Target = smi2  (all {len(PAIRS)} pairs)")
    print("=" * 120)
    print()

    sub_hdr = (f" # {'Pair':<30s} {'Category':<20s} "
               f"{'SMSD(us)':>10s} {'RDKit(us)':>10s} "
               f"{'SMSD hit':>8s} {'RDKit hit':>9s} "
               f"{'Speedup':>10s}")
    print(sub_hdr)
    print("-" * 100)

    sub_smsd_wins = sub_rdkit_wins = sub_ties = 0
    agree = disagree = 0

    for idx, (smi1, smi2, name, category) in enumerate(PAIRS, 1):
        try:
            s_times, s_hit = bench_smsd_sub(smi1, smi2)
            s_med = statistics.median(s_times)
        except Exception:
            s_med = float("inf")
            s_hit = False

        try:
            r_times, r_hit = bench_rdkit_sub(smi1, smi2)
            r_med = statistics.median(r_times)
        except Exception:
            r_med = float("inf")
            r_hit = False

        # Agreement check
        if s_hit == r_hit:
            agree += 1
        else:
            disagree += 1

        # Speedup
        if s_med == float("inf") or r_med == float("inf") or s_med == 0 or r_med == 0:
            speedup_str = "N/A"
            sub_ties += 1
        else:
            ratio = r_med / s_med
            if ratio > 1.1:
                speedup_str = f"SMSD {ratio:.1f}x"
                sub_smsd_wins += 1
            elif ratio < 0.9:
                speedup_str = f"RDKit {1/ratio:.1f}x"
                sub_rdkit_wins += 1
            else:
                speedup_str = "~tie"
                sub_ties += 1

        def fmt_us(val):
            if val == float("inf"):
                return "ERR"
            if val >= 1_000_000:
                return f"{val/1_000_000:.1f}s"
            if val >= 1_000:
                return f"{val/1_000:.2f}ms"
            return f"{val:.1f}us"

        print(f"{idx:2d} {name:<30s} {category:<20s} "
              f"{fmt_us(s_med):>10s} {fmt_us(r_med):>10s} "
              f"{'yes' if s_hit else 'no':>8s} {'yes' if r_hit else 'no':>9s} "
              f"{speedup_str:>10s}")

        sub_results.append(SubResult(
            name=name, category=category,
            smsd_median_us=s_med, smsd_hit=s_hit,
            rdkit_median_us=r_med, rdkit_hit=r_hit,
        ))

    print("-" * 100)
    print()
    print("SUBSTRUCTURE SUMMARY")
    print(f"  Speed wins:  SMSD={sub_smsd_wins}  RDKit={sub_rdkit_wins}  tie={sub_ties}")
    print(f"  Hit agreement: {agree}/{len(PAIRS)} pairs agree")
    if disagree > 0:
        print(f"  WARNING: {disagree} pair(s) disagree — verify SMILES semantics")

    tsv_path = Path(output_path) if output_path else SCRIPT_DIR / "results_smsd_vs_rdkit_20pairs_mcs_sub.tsv"
    write_results_tsv(tsv_path, results, sub_results)
    print(f"\n  Results written to {tsv_path}")


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
    parser.add_argument("--compare-mode", choices=["defaults", "strict", "fmcs"],
                        help="defaults = toolkit defaults, strict = exact/ring-parity baseline, fmcs = loose FMCS-style baseline")
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
