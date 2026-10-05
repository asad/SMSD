#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
"""Generate reproducible Dalke-style pairs from a specified molecule collection.

These are derived corpora, not the original ChEMBL-13 FMCS benchmark. A nearest
neighbor is the highest Morgan-radius-2/1024-bit Tanimoto match after explicitly
excluding the query index; no similarity cutoff is applied. Existing checked-in
corpora are not overwritten by default.
"""
import argparse
import hashlib
from pathlib import Path
import random

DATA = Path(__file__).resolve().parent / "data"


def sample_pool(items, maximum, rng):
    return [items[i] for i in rng.sample(range(len(items)), maximum)] if len(items) > maximum else list(items)


def random_pair_indices(count, pairs, rng):
    if count < 2 or pairs > count*(count-1)//2:
        raise ValueError("not enough distinct unordered index pairs")
    seen, result = set(), []
    while len(result) < pairs:
        i, j = rng.sample(range(count), 2)
        key = tuple(sorted((i, j)))
        if key not in seen:
            seen.add(key); result.append((i, j))
    return result


def nearest_neighbor_index(query, similarities):
    candidates = (i for i in range(len(similarities)) if i != query)
    # Break ties by index, rather than assuming the query sorts before its neighbors.
    return max(candidates, key=lambda i: (similarities[i], -i))


def main(argv=None):
    from rdkit import Chem, DataStructs
    from rdkit.Chem import rdFingerprintGenerator
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DATA / "chembl_mcs_benchmark.smi")
    parser.add_argument("--output-dir", type=Path, default=Path("build/local-benchmarks/generated-pairs"))
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--pairs", type=int, default=1000)
    parser.add_argument("--max-molecules", type=int, default=5000)
    args = parser.parse_args(argv)
    if args.pairs < 1 or args.max_molecules < 2:
        parser.error("pairs must be positive and max-molecules at least 2")
    rng = random.Random(args.seed)
    items = []
    for line in args.input.read_text().splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        row = line.split()
        mol = Chem.MolFromSmiles(row[0])
        if mol is not None and mol.GetNumHeavyAtoms() >= 5:
            items.append((mol, Chem.MolToSmiles(mol), row[1] if len(row)>1 else f"mol_{len(items)}"))
    items = sample_pool(items, args.max_molecules, rng)
    pairs = random_pair_indices(len(items), args.pairs, rng)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    metadata = (f"# Derived Dalke-style corpus; not the original ChEMBL-13 data\n"
                f"# source_file={args.input.name}; source_sha256={hashlib.sha256(args.input.read_bytes()).hexdigest()}\n"
                f"# seed={args.seed}; pool_size={len(items)}; rdkit={Chem.rdBase.rdkitVersion}\n")
    with (args.output_dir / "dalke_random_pairs.tsv").open("w") as handle:
        handle.write(metadata+"# SMILES1\tSMILES2\tName1\tName2\n")
        for i, j in pairs:
            handle.write(f"{items[i][1]}\t{items[j][1]}\t{items[i][2]}\t{items[j][2]}\n")
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=1024)
    fingerprints = [generator.GetFingerprint(item[0]) for item in items]
    queries = rng.sample(range(len(items)), min(args.pairs, len(items)))
    with (args.output_dir / "dalke_nn_pairs.tsv").open("w") as handle:
        handle.write(metadata+"# One nearest neighbor per query; query index excluded; no similarity cutoff\n")
        handle.write("# SMILES1\tSMILES2\tName1\tName2\tTanimoto\n")
        for query in queries:
            similarities = DataStructs.BulkTanimotoSimilarity(fingerprints[query], fingerprints)
            neighbor = nearest_neighbor_index(query, similarities)
            handle.write(f"{items[query][1]}\t{items[neighbor][1]}\t{items[query][2]}\t{items[neighbor][2]}\t{similarities[neighbor]:.4f}\n")
    print(f"Wrote {len(pairs)} random and {len(queries)} nearest-neighbor pairs to {args.output_dir}")


if __name__ == "__main__":
    main()
