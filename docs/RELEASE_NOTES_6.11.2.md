# SMSD Pro 6.11.2

**Syed Asad Rahman — BioInception PVT LTD**

This release focuses on correctness fixes, performance improvements, and
cleaning up some long-standing naming issues in the API.


## What changed

### Bug fixes

Fixed a memory leak in `SearchEngine` where the `IdentityHashMap` cache
kept growing indefinitely during batch workloads. Switched to `WeakHashMap`.

The Bron-Kerbosch color bound was using a single 64-bit mask, which broke
silently when the product graph needed more than 64 colors. Replaced with
a proper multi-word array.

CIP assignment was missing Rule 3 (Z > E) from the IUPAC 2013 spec. Added
it in both Java and C++.

Added `volatile` to the lazy-init fields in `MolGraph` (`adj`, NLF caches,
pharmacophore features) to prevent torn reads in multi-threaded MCS searches.

### Chemistry

Updated tautomer equilibrium weights based on literature review:
nitroso-oxime to 0.95 (was 0.70), nitro-aci to 0.95 (was 0.25),
pyridone-hydroxypyridine to 0.95 (was 0.85). Added selenium to the
tautomer compatibility table and iodine to reaction-aware rarity scoring.
Simplified a redundant `Cycles.or(all(), all(6))` call.

Fixed the SAH (S-adenosylhomocysteine) test SMILES which had the ribose
connected through an ester oxygen instead of through sulfur. Also relaxed
formal charge matching in reaction-aware MCS since reactions inherently
change charges (e.g. SAM S+ becomes SAH S).

### API naming

Renamed all `Mcs*` types to `MCS*` for consistent casing. Renamed the
`tanimoto` field to `overlapCoefficient` since the metric is actually
Szymkiewicz-Simpson (size/min), not Tanimoto (intersection/union).
Old names still work but are deprecated.

### Performance

Reworked several hot paths in the MCS pipeline:

- Bucketed compatibility graph construction, so we only compare atoms
  that could actually match (same element/aromaticity)
- Flat `int[]` arrays through the pipeline instead of `std::map` /
  `HashMap` at every stage
- Ring-system symmetry detection (union-find + hash signatures) to cut
  backtracking on fused ring systems
- Flat-bucket domain init in VF2++ replacing hash maps with prefix-sum arrays
- Fused `popcountAnd`/`anyBitAnd` bitset operations with ARM NEON paths
  for Apple Silicon
- Stage-aware routing that skips unnecessary pipeline stages based on
  intermediate results
- Added `MCSStageTimers` profiling API for diagnosing pipeline bottlenecks

New Python APIs: `TargetCorpus` for parse-once/query-many batch workflows,
`batch_find_substructure()` returning atom mappings.

### Other

- VF3 switch is now exhaustive
- Pre-allocated domain candidate buffers (capped at 4096)
- `standardiseSafe` logs exceptions instead of swallowing them
- CLI SMILES parser is `ThreadLocal` now
- SDF batch loading capped at 100K molecules


## Compatibility

- Java 25+, C++17, Python 3.9+
- GPU: Metal (Apple Silicon), CUDA (Volta+)


## Copyright

Copyright (c) 2018-2026 BioInception PVT LTD
Algorithm copyright (c) 2009-2026 Syed Asad Rahman
