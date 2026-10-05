#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
# Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
# See the NOTICE file for attribution, trademark, and algorithm IP terms.
"""
Historical ten-pair compatibility entrypoint; uses the maintained protocol.

Both libraries called purely from Python on the same molecule pairs,
same hardware, same process — no JVM, no subprocess overhead.

Usage:
    pip install smsd rdkit
    python benchmarks/benchmark_python_vs_rdkit.py

Requirements:
    smsd (the intended baseline or candidate wheel)
    rdkit >= 2022.03 (pip install rdkit)
"""

import sys
import time
import platform
from statistics import median as stat_median

# ---------------------------------------------------------------------------
# Molecule pairs — identical to benchmark_all.py
# ---------------------------------------------------------------------------

VANCOMYCIN = (
    "CC1C(C(CC(O1)OC2C(C(C(OC2OC3=C4C=C5C=C3OC6=C(C=C(C=C6)C(C(C(=O)"
    "NC(C(=O)NC5C(=O)NC7C8=CC(=C(C=C8)O)C9=C(C=C(C=C9O)O)C(NC(=O)C("
    "C(C1=CC(=C(O4)C=C1)Cl)O)NC7=O)C(=O)O)CC(=O)N)NC(=O)C(CC(C)C)NC)"
    "O)Cl)CO)O)O)(C)N)O"
)

PAIRS = [
    # (smiles_1, smiles_2, name, category)
    ("c1ccccc1",  "Oc1ccccc1",                                 "benzene-phenol",           "Small aromatic"),
    ("Cn1cnc2c1c(=O)n(C)c(=O)n2C",
     "Cn1cnc2c1c(=O)[nH]c(=O)n2C",                            "caffeine-theophylline",    "N-methyl diff"),
    ("CC(=O)Oc1ccccc1C(=O)O",
     "CC(=O)Nc1ccc(O)cc1",                                     "aspirin-acetaminophen",    "Drug pair"),
    ("CC(C)Cc1ccc(CC(C)C(=O)O)cc1",
     "COc1ccc2cc(CC(C)C(=O)O)ccc2c1",                          "ibuprofen-naproxen",       "NSAID"),
    ("CN1CCC23C4C1CC5=C(C2C(C=C4)O3)C=C(C=C5)O",
     "CN1CCC23C4C1CC5=C(C2C(C=C4)OC3)C=C(C=C5)O",             "morphine-codeine",         "Alkaloid"),
    ("C1=NC(=C2C(=N1)N(C=N2)C3C(C(C(O3)COP(=O)(O)OP(=O)(O)OP(=O)(O)O)O)O)N",
     "C1=NC(=C2C(=N1)N(C=N2)C3C(C(C(O3)COP(=O)(O)OP(=O)(O)O)O)O)N",
                                                                "ATP-ADP",                  "Nucleotide"),
    ("CC(C)C1=C(C(=C(N1CCC(CC(CC(=O)O)O)O)C2=CC=C(C=C2)F)C3=CC=CC=C3)C(=O)NC4=CC=CC=C4",
     "CC(C)C1=NC(=NC(=C1C=CC(CC(CC(=O)O)O)O)C2=CC=C(C=C2)F)N(C)S(=O)(=O)C",
                                                                "atorvastatin-rosuvastatin","Statin"),
    ("CC1=C2C(C(=O)C3(C(CC4C(C3C(C(C2(C)C)(CC1OC(=O)C(C(C5=CC=CC=C5)NC(=O)C6=CC=CC=C6)O)O)"
     "OC(=O)C7=CC=CC=C7)(CO4)OC(=O)C)O)C)OC(=O)C",
     "CC1=C2C(C(=O)C3(C(CC4C(C3C(C(C2(C)C)(CC1OC(=O)C(C(C5=CC=CC=C5)NC(=O)OC(C)(C)C)O)O)"
     "OC(=O)C6=CC=CC=C6)(CO4)OC(=O)C)O)C)O",                  "paclitaxel-docetaxel",     "Taxane"),
    ("CCC1C(C(C(C(=O)C(CC(C(C(C(C(C(=O)O1)C)OC2CC(C(C(O2)C)O)(C)OC)C)OC3C(C(CC(O3)C)N(C)C)O)(C)O)C)C)O)(C)O",
     "CCC1C(C(C(N(CC(CC(C(C(C(C(C(=O)O1)C)OC2CC(C(C(O2)C)O)(C)OC)C)OC3C(C(CC(O3)C)N(C)C)O)(C)O)C)C)C)O)(C)O",
                                                                "erythromycin-azithromycin","Macrolide"),
    (VANCOMYCIN, VANCOMYCIN,                                    "vancomycin-self",          "Self-match"),
]

NUM_RUNS   = 5
TIMEOUT_S  = 10  # RDKit timeout per pair


# ---------------------------------------------------------------------------
# RDKit runner
# ---------------------------------------------------------------------------

def run_rdkit(mol1, mol2):
    from rdkit.Chem import rdFMCS
    result = rdFMCS.FindMCS([mol1, mol2], timeout=TIMEOUT_S)
    return result.numAtoms, result.canceled


def bench_rdkit(pairs):
    from rdkit import Chem
    rows = []
    for smi1, smi2, name, cat in pairs:
        m1 = Chem.MolFromSmiles(smi1)
        m2 = Chem.MolFromSmiles(smi2)
        if m1 is None or m2 is None:
            rows.append((name, cat, -1.0, -1, False))
            continue
        # warmup
        run_rdkit(m1, m2)
        times = []
        mcs_n, timed_out = -1, False
        for _ in range(NUM_RUNS):
            t0 = time.perf_counter()
            mcs_n, timed_out = run_rdkit(m1, m2)
            times.append((time.perf_counter() - t0) * 1000.0)
        rows.append((name, cat, min(times), mcs_n, timed_out))
    return rows


# ---------------------------------------------------------------------------
# SMSD Python runner
# ---------------------------------------------------------------------------

def run_smsd(g1, g2):
    from smsd import ChemOptions, MCSOptions
    from smsd._smsd import find_mcs
    opts = MCSOptions()
    opts.timeout_ms = int(TIMEOUT_S * 1000)
    mapping = find_mcs(g1, g2, ChemOptions(), opts)
    return len(mapping)


def bench_smsd(pairs):
    from smsd import parse_smiles
    rows = []
    for smi1, smi2, name, cat in pairs:
        try:
            g1 = parse_smiles(smi1)
            g2 = parse_smiles(smi2)
        except Exception:
            rows.append((name, cat, -1.0, -1))
            continue
        # warmup
        run_smsd(g1, g2)
        times = []
        mcs_n = -1
        for _ in range(NUM_RUNS):
            t0 = time.perf_counter()
            mcs_n = run_smsd(g1, g2)
            times.append((time.perf_counter() - t0) * 1000.0)
        rows.append((name, cat, min(times), mcs_n))
    return rows


# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------

def fmt(ms):
    if ms < 0:    return "  N/A"
    if ms < 1.0:  return f"{ms:6.3f}"
    if ms < 100:  return f"{ms:6.1f}"
    return f"{ms:6.0f}"


def speedup(rdkit_ms, smsd_ms):
    if rdkit_ms <= 0 or smsd_ms <= 0:
        return "  N/A"
    r = rdkit_ms / smsd_ms
    if r > 1.05:   return f"{r:5.1f}×"
    if r < 0.95:   return f"0.{int(100/r):02d}×"
    return " ~1×"


def main():
    # Compatibility entrypoint retaining its historical ten-pair corpus.
    import runpy
    import benchmark_python
    benchmark_python.PAIRS = PAIRS
    import argparse
    parser = argparse.ArgumentParser(description="Historical ten-pair corpus using the maintained validity-aware protocol")
    parser.add_argument("--timeout-sec", type=int, default=10)
    parser.add_argument("--iters", type=int, default=3)
    parser.add_argument("--warmup", type=int, default=1)
    parser.add_argument("--compare-mode", choices=("strict", "fmcs", "any"), default="fmcs")
    parser.add_argument("--output")
    args = parser.parse_args()
    return benchmark_python.main(mcs_only=True, timeout_sec=args.timeout_sec,
        iters=args.iters, warmup=args.warmup, compare_mode=args.compare_mode, output_path=args.output)


if __name__ == "__main__":
    sys.exit(main())
