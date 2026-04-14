#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
# Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
# See the NOTICE file for attribution, trademark, and algorithm IP terms.
"""
Benchmark RDKit FindMCS on diverse molecule pairs.

Usage:
    pip install rdkit
    python benchmark_rdkit.py

Compares against SMSD Java benchmark (benchmark_java.sh) using
identical molecule pairs ranging from tiny (2 atoms) to complex
macrocycles and symmetric graphs.
"""

from rdkit import Chem
from rdkit.Chem import rdFMCS
import time
import sys

# 15 molecule pairs from tiny to most complex:
PAIRS = [
    # Tiny (2-6 atoms)
    ("C", "CC", "methane-ethane"),
    ("O", "CO", "water-methanol"),

    # Small aromatic (6-12 atoms)
    ("c1ccccc1", "c1ccc(O)cc1", "benzene-phenol"),
    ("c1ccccc1", "Cc1ccccc1", "benzene-toluene"),

    # Medium drugs (13-25 atoms)
    ("CC(=O)Oc1ccccc1C(=O)O", "CC(=O)Nc1ccc(O)cc1", "aspirin-acetaminophen"),
    ("Cn1cnc2c1c(=O)n(C)c(=O)n2C", "Cn1cnc2c1c(=O)[nH]c(=O)n2C", "caffeine-theophylline"),
    ("CC(C)Cc1ccc(CC(C)C(O)=O)cc1", "COc1ccc2cc(CC(C)C(O)=O)ccc2c1", "ibuprofen-naproxen"),

    # Large drugs (20-30 atoms)
    ("CN1CCC23C4C1CC5=C2C(=C(C=C5)O)OC3C(C=C4)O",
     "COC1=CC=C2C3CC4=CC(=C(C=C4C3CC5=CC1=C2O5)OC)O", "morphine-codeine"),

    # Nucleotides (25-40 atoms)
    ("C1=NC(=C2C(=N1)N(C=N2)C3C(C(C(O3)COP(=O)(O)OP(=O)(O)OP(=O)(O)O)O)O)N",
     "C1=NC2=C(N1C3C(C(C(O3)COP(=O)(O)OP(=O)(O)OP(=O)(O)O)O)O)NC(=NC2=O)N", "ATP-GTP"),

    # Very large (40-60 atoms)
    ("CC1=C2C(C(=O)C3(C(CC4C(C3C(C(C2(C)C)(CC1OC(=O)C(C(C5=CC=CC=C5)NC(=O)C6=CC=CC=C6)O)O)OC(=O)C7=CC=CC=C7)(CO4)OC(=O)C)O)C)OC(=O)C",
     "CC1=C2C(C(=O)C3(C(CC4C(C3C(C(C2(C)C)(CC1OC(=O)C(C(C5=CC=CC=C5)NC(=O)OC(C)(C)C)O)O)OC(=O)C6=CC=CC=C6)(CO4)OC(=O)C)O)C)O", "paclitaxel-docetaxel"),

    # Symmetric nightmares
    ("C1C2CC3CC1CC(C2)C3", "C1C2CC3CC1CC(C2)C3", "adamantane-self"),
    ("c1cc2ccc3ccc4ccc5ccc6ccc1c7c2c3c4c5c67", "c1cc2ccc3ccc4ccc5ccc6ccc1c7c2c3c4c5c67", "coronene-self"),

    # Macrocycles
    ("CCC1C(C(C(C(=O)C(CC(C(C(C(C(C(=O)O1)C)OC2CC(C(C(O2)C)O)(C)OC)C)OC3C(C(CC(O3)C)N(C)C)O)(C)O)C)C)O)(C)O",
     "CCC1C(C(C(N(CC(CC(C(C(C(C(C(=O)O1)C)OC2CC(C(C(O2)C)O)(C)OC)C)OC3C(C(CC(O3)C)N(C)C)O)(C)O)C)C)C)O)(C)O", "erythromycin-azithromycin"),

    # Polymer
    ("OCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCO", "OCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCO", "PEG12-PEG16"),
]

TIMEOUT_SEC = 10
NUM_RUNS = 5


def main():
    print(f"{'Pair':30s}  {'Best(ms)':>10s}  {'Median(ms)':>10s}  {'MCS':>4s}")
    print("-" * 62)

    total_time = 0.0
    errors = 0

    for smi1, smi2, name in PAIRS:
        mol1 = Chem.MolFromSmiles(smi1)
        mol2 = Chem.MolFromSmiles(smi2)
        if mol1 is None or mol2 is None:
            print(f"{name:30s}  PARSE_ERROR")
            errors += 1
            continue

        # Warmup
        try:
            rdFMCS.FindMCS([mol1, mol2], timeout=1)
        except Exception:
            pass

        # Timed runs
        times = []
        mcs_size = 0
        for _ in range(NUM_RUNS):
            t0 = time.perf_counter()
            try:
                result = rdFMCS.FindMCS([mol1, mol2], timeout=TIMEOUT_SEC)
                mcs_size = result.numAtoms
            except Exception:
                mcs_size = -1
            dt = (time.perf_counter() - t0) * 1000
            times.append(dt)

        best = min(times)
        median = sorted(times)[len(times) // 2]
        total_time += best
        print(f"{name:30s}  {best:10.3f}  {median:10.3f}  {mcs_size:>4d}")

    print("-" * 62)
    print(f"{'TOTAL':30s}  {total_time:10.3f} ms")
    if errors:
        print(f"  ({errors} parse errors)")

    return 0


if __name__ == "__main__":
    sys.exit(main())
