#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
"""
Generate Dalke-style MCS benchmark pairs from a molecule collection.

Produces two benchmark sets:
  1. Random pairs (low similarity) — 1000 pairs
  2. Nearest-neighbor pairs (high similarity, k=2) — 1000 pairs

Usage:
    python benchmarks/generate_dalke_pairs.py

Input:  benchmarks/data/chembl_mcs_benchmark.smi
Output: benchmarks/data/dalke_random_pairs.tsv
        benchmarks/data/dalke_nn_pairs.tsv
"""
import random
import os
import sys

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(SCRIPT_DIR, "data")

def main():
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem, DataStructs
    except ImportError:
        print("ERROR: RDKit required. Install with: pip install rdkit")
        sys.exit(1)

    smi_file = os.path.join(DATA_DIR, "chembl_mcs_benchmark.smi")
    if not os.path.exists(smi_file):
        print(f"ERROR: {smi_file} not found")
        sys.exit(1)

    # Load molecules
    print("Loading molecules...")
    mols = []
    smiles_list = []
    names = []
    with open(smi_file) as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) < 1:
                continue
            smi = parts[0]
            name = parts[1] if len(parts) > 1 else f"mol_{len(mols)}"
            mol = Chem.MolFromSmiles(smi)
            if mol is not None and mol.GetNumHeavyAtoms() >= 5:
                mols.append(mol)
                smiles_list.append(Chem.MolToSmiles(mol))
                names.append(name)

    print(f"Loaded {len(mols)} valid molecules (>= 5 heavy atoms)")

    # Cap at 5000 for manageable fingerprint computation
    if len(mols) > 5000:
        idx = random.sample(range(len(mols)), 5000)
        mols = [mols[i] for i in idx]
        smiles_list = [smiles_list[i] for i in idx]
        names = [names[i] for i in idx]
        print(f"Subsampled to {len(mols)} molecules")

    # === 1. Random pairs (1000) ===
    print("Generating 1000 random pairs...")
    random.seed(42)
    random_out = os.path.join(DATA_DIR, "dalke_random_pairs.tsv")
    with open(random_out, "w") as f:
        f.write("# Dalke-style random MCS benchmark pairs\n")
        f.write("# Source: chembl_mcs_benchmark.smi (MoleculeNet drug collections)\n")
        f.write("# SMILES1\\tSMILES2\\tName1\\tName2\n")
        seen = set()
        count = 0
        while count < 1000:
            i, j = random.sample(range(len(mols)), 2)
            key = (min(i, j), max(i, j))
            if key in seen:
                continue
            seen.add(key)
            f.write(f"{smiles_list[i]}\t{smiles_list[j]}\t{names[i]}\t{names[j]}\n")
            count += 1
    print(f"  Written to {random_out}")

    # === 2. Nearest-neighbor pairs (1000) ===
    print("Computing Morgan fingerprints for nearest-neighbor search...")
    fps = [AllChem.GetMorganFingerprintAsBitVect(m, 2, 1024) for m in mols]

    nn_out = os.path.join(DATA_DIR, "dalke_nn_pairs.tsv")
    # Pick 1000 random query molecules, find their nearest neighbor
    query_indices = random.sample(range(len(mols)), min(1000, len(mols)))

    print("Generating nearest-neighbor pairs...")
    with open(nn_out, "w") as f:
        f.write("# Dalke-style nearest-neighbor MCS benchmark pairs (k=2)\n")
        f.write("# Source: chembl_mcs_benchmark.smi (MoleculeNet drug collections)\n")
        f.write("# SMILES1\\tSMILES2\\tName1\\tName2\\tTanimoto\n")
        for qi in query_indices:
            sims = DataStructs.BulkTanimotoSimilarity(fps[qi], fps)
            # Sort by descending similarity, skip self (index 0)
            ranked = sorted(range(len(sims)), key=lambda x: sims[x], reverse=True)
            nn = ranked[1]  # nearest neighbor (not self)
            tc = sims[nn]
            f.write(f"{smiles_list[qi]}\t{smiles_list[nn]}\t{names[qi]}\t{names[nn]}\t{tc:.4f}\n")
    print(f"  Written to {nn_out}")

    print("Done!")

if __name__ == "__main__":
    main()
