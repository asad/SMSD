# Changelog

All notable changes to SMSD Pro are documented in this file.

## [6.11.2] - 2026-04-07

### Summary
Correctness, performance, and API cleanup release.

### Included
- Fixed memory leak in SearchEngine cache (`IdentityHashMap` replaced with `WeakHashMap`)
- Fixed BK color-bound overflow when product graph exceeds 64 colors (multi-word array)
- Added missing CIP Rule 3 (Z > E) per IUPAC 2013 in Java and C++
- Thread-safety: `volatile` on lazy-init fields in MolGraph
- Updated tautomer weights: nitroso-oxime 0.95, nitro-aci 0.95, pyridone 0.95
- Added selenium to tautomer compatibility, iodine to reaction-aware scoring
- Corrected SAH test SMILES (thioether, not ester connectivity)
- Relaxed formal charge matching in reaction-aware MCS
- Renamed `Mcs*` types to `MCS*`, `tanimoto` to `overlapCoefficient`
- Bucketed compatibility graph builder for faster MCS construction
- Flat array MCS pipeline replacing `std::map` / `HashMap` at each stage
- Ring-system symmetry breaking via union-find for reduced backtracking
- VF2++ flat-bucket domain init with prefix-sum arrays
- Fused bitset operations with ARM NEON paths for Apple Silicon
- Stage-aware pipeline routing to skip unnecessary MCS stages
- `MCSStageTimers` profiling API for pipeline diagnostics
- `TargetCorpus` and `batch_find_substructure()` Python APIs
- Exhaustive VF3 switch, pre-allocated buffers, `ThreadLocal` SMILES parser
- SDF batch cap at 100K molecules

## [6.11.1] - 2026-04-04

### Bug Fixes
- **ECFP initial invariants**: added missing Rogers & Hahn invariants #3 (bond order
  sum / valence) and #4 (atomic mass number) to both binary and count ECFP variants
  in C++ and Java — fingerprints now encode the full 7-invariant set per the 2010 paper
- **Path fingerprint canonical hash**: replaced double-bit-setting (forward + reverse)
  with `min(fwd, rev)` single canonical hash, correcting bit density inflation
- **FCFP pyrrole-N misclassification**: aromatic nitrogen acceptor classification now
  uses direct hydrogen count instead of bond-sum heuristic, fixing incorrect pyridine-N
  classification as non-acceptor (pyridine N has a free lone pair; pyrrole N does not)
- **Thread safety**: `prewarmGraph()` now calls `ensurePatternFP()` before entering
  OpenMP parallel regions, preventing data races on pattern fingerprint lazy init
- **Dead code removal**: removed unused `seenHashes` set from binary ECFP (count
  variant correctly uses hash-space convergence; binary variant uses bit-space)

## [6.11.0] - 2026-04-04

### Summary
Performance, precision, and depiction release: cache-optimal data structures,
pre-indexed McGregor DFS, publication-quality SVG renderer (ACS 1996 standard),
comprehensive layout engine, 35+ new Python bindings.

### Included
- Core engine: converted hot-path `vector<bool>` to `vector<uint8_t>` for 15-25%
  cache performance improvement across McGregor DFS, BK partition bound, seed-extend
- Pre-indexed candidate sets in McGregor DFS using `compatTargets_[]` — eliminates
  O(n^2) linear scan per frontier atom
- Publication-quality SVG depiction engine (ACS 1996 standard):
  - Jmol/CPK element colors, asymmetric double bonds, wedge/dash stereo bonds
  - Bond-to-label clipping, H-count subscripts, charge superscripts
  - Full customization via DepictOptions (bond_length, colors, fonts, sizes)
  - Side-by-side MCS pair rendering with atom-atom mapping numbers
- 8-phase 2D layout pipeline: template match, ring-first, chain zig-zag, force
  refinement, overlap resolution, crossing reduction, canonical orientation,
  bond-length normalisation
- Distance geometry 3D layout with power iteration eigendecomposition and
  force-field refinement
- 40+ ring scaffold templates (pharmaceutical scaffolds, PAH, spiro, bridged)
- Full 2D/3D coordinate transform suite (translate, rotate, scale, mirror,
  center, align, project, lift)
- 35+ new Python bindings with GIL release for thread safety
- Java: explicit per-atom type matching for robust handling of exotic valence states
- 9 precision chemistry tests (azulene, pyrene, pyridinium, cyclopentadienyl,
  boron, sulfoxide, phosphate, E/Z stereo)
- 27 new layout engine tests (2D/3D generation, transforms, overlaps)
- Comprehensive Python documentation with examples and cautions

## [6.10.2] - 2026-04-03

### Summary
Correctness release: fixed MCS connectivity filter for non-induced mode,
added regression tests for challenging molecule pairs.

### Included
- Corrected connected-component filter to enforce common-bond reachability
  in both query and target molecules (non-induced MCS mode)
- Added GOLDEN_843 regression tests in Python and Java (timeout and size)
- Version bump to 6.10.2

## [6.10.1] - 2026-04-03

### Summary
Stability and correctness release: hardened MCS repair pipeline, deterministic
tests, CI/CD fixes.

### Included
- Rewrote MCS mapping repair to iterative bounded loop — eliminates unbounded
  recursion on large molecules (vancomycin, CoA, paclitaxel)
- Correct duplicate-target handling in mapping repair
- Removed all timing-dependent test assertions — algorithmic correctness is
  now fully deterministic and machine-speed independent
- Removed `forkedProcessTimeoutInSeconds` from Surefire (was killing fork JVM)
- Switched Python publish to manual dispatch only (no auto-publish on release)
- Fixed GitHub Actions artifact version references and Node.js 24 opt-in
- Fixed `atomWeights` array length for benzene queries
- Adjusted MCS thresholds and completeRingsOnly tests for edge cases

## [6.9.0] - 2026-04-02

### Summary
Core chemistry correctness, native I/O hardening, and benchmark alignment release.

### Included
- Direction-stable native/public MCS handling for hard asymmetric pairs
- Symmetric `ringMatchesRingOnly` semantics across C++, Python, and Java
- Mode-matched benchmark leaderboards with explicit `defaults`, `strict`, and `fmcs` comparison modes
- Release-documentation cleanup and benchmark/report alignment
- Native MDL MOL V2000 metadata preservation for molecule name, program line, comment, and SDF properties
- Native MDL MOL V3000 reader and writer for core graph round-trip
- Native patent-style R-group molfile support via `R#` pseudo-atoms and `M  RGP`
- Stronger SMILES and SMARTS attachment-point handling, including labeled placeholders such as `[R1]`
- Additional stereo and CIP round-trip coverage across SMILES and molfile paths
- Python bindings for native mol block read/write APIs
- Java MolGraph metadata parity fields for release alignment

## [6.8.1] - 2026-04-02

### Summary
Release alignment and parity cleanup.

### Included
- Repo-wide version bump
- Python packaging alignment
- Java and C++ release metadata alignment
- Benchmark and documentation refresh

## [6.8.0] - 2026-04-01

### Summary
Open-source release of the SMSD Pro cheminformatics toolkit.

### Included
- Substructure search engine
- Maximum common subgraph (MCS) computation
- Circular fingerprints (ECFP/FCFP) with tautomer awareness
- SMARTS pattern matching
- Molecular similarity and screening
- CIP R/S and E/Z stereodescriptor assignment
- Batch processing with optional GPU acceleration
- Java 21+ / C++17 / Python 3.8+ support
- CLI, SDF batch, and JSON export
