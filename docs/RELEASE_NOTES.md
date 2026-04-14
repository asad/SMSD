# SMSD Pro 7.1.0

**Syed Asad Rahman — BioInception PVT LTD**

Major release: unified Python and Java API, clean break from legacy aliases,
full Java parity for convenience methods, all tests green.


## What's New

### Unified Python API — clean break

Two entry points replace all previous aliases:

```python
import smsd

# MCS — returns dict (single) or list[dict] (multiple)
mcs = smsd.find_mcs("c1ccccc1", "c1ccc(O)cc1")
mcs = smsd.find_mcs(mol1, mol2, max_results=5)

# Substructure — returns dict (single) or list[dict] (multiple)
hit = smsd.find_substructure("c1ccccc1", "c1ccc(O)cc1")
hit = smsd.find_substructure(query, target, max_results=3)

# Boolean convenience check
if smsd.is_substructure(query, target):
    print("Query is a substructure of target")
```

Both functions accept SMILES strings, MolGraph objects, or RDKit Mol objects.

### Removed legacy aliases

The following names are no longer available in v7.0.0:

| Removed | Replacement |
|---------|-------------|
| `smsd.mcs()` | `smsd.find_mcs()` |
| `smsd.substructure_search()` | `smsd.find_substructure()` |
| `smsd.all_mcs()` | `smsd.find_mcs(mol1, mol2, max_results=N)` |
| `smsd.overlapCoefficient()` | `smsd.overlap_coefficient()` |
| `smsd.tanimoto()` | `smsd.tanimoto_coefficient()` |
| `smsd.count_overlap_coefficient()` | (use C++ binding directly) |
| `smsd.count_tanimoto()` | (use C++ binding directly) |

### Java parity — unified convenience methods

```java
import com.bioinception.smsd.core.SearchEngine;

// MCS with default options
Map<Integer, Integer> mcs = SearchEngine.findMCS(g1, g2);

// Substructure with default options
Map<Integer, Integer> hit = SearchEngine.findSubstructure(query, target);

// With custom options and timeout
Map<Integer, Integer> hit = SearchEngine.findSubstructure(query, target, chemOpts, 10_000L);
```

All overloads work with both `MolGraph` and CDK `IAtomContainer` inputs.

### Internal improvements

- Raw C++ bindings renamed to `_native_find_mcs`, `_native_find_all_mcs`,
  `_native_is_substructure`, `_native_find_substructure` — clearly internal.
- All internal calls use the unified API (`find_mcs`, `find_substructure`).
- `mcs_from_smiles()`, `mcs_rdkit()`, `substructure_rdkit()`, and
  `depict_mcs()` all updated to the new names.


## Migration Guide

### Python

```diff
- mapping = smsd.mcs(mol1, mol2)
+ mapping = smsd.find_mcs(mol1, mol2)

- mapping = smsd.substructure_search(query, target)
+ mapping = smsd.find_substructure(query, target)

- results = smsd.all_mcs(mol1, mol2, max_results=5)
+ results = smsd.find_mcs(mol1, mol2, max_results=5)
```

### Java

No breaking changes. New convenience methods added; existing methods unchanged.


## Benchmark

Dalke NN dataset (1,000 high-similarity ChEMBL pairs):

| Metric | SMSD Pro 7.0.0 | RDKit FindMCS 2026.03 |
|--------|:---------------:|:---------------------:|
| Total time | **40 s** | 213 s |
| Median time | 0.6 ms | 0.4 ms |
| Mean MCS size | **25.8 atoms** | 25.0 atoms |
| Timeouts | **0** | 8 |
| Larger-MCS wins | **211 (21 %)** | 29 (3 %) |

5x faster overall, finds larger MCS 7x more often, zero timeouts.


## Test Status

- Python: 310 passed, 0 failed
- Java: BUILD SUCCESS (581+ tests)


## Compatibility

- Java 25+, C++17, Python 3.10-3.13
- GPU: Metal (Apple Silicon), CUDA (Volta+)
- Platforms: macOS (arm64, x86_64), Linux (x86_64, aarch64), Windows (AMD64)


## Copyright

Copyright (c) 2018-2026 BioInception PVT LTD
Algorithm copyright (c) 2009-2026 Syed Asad Rahman
Licensed under Apache License 2.0. See NOTICE for details.
