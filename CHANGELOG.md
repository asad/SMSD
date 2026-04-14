# Changelog

All notable changes to SMSD Pro are documented in this file.

## [7.1.0] - 2026-04-14

### Added
- Comprehensive test suite: 597 tests across 9 test files
  - `test_api_coverage.py` (160 tests): MCS utils, batch ops, TargetCorpus,
    file I/O, scaffolds, coordinate transforms, depiction, SMARTS, enums
  - `test_fingerprints.py` (128 tests): all FP types, similarity metrics,
    edge cases (single-atom, disjoint, empty), challenging molecules (taxol,
    C60, cubane, morphine/codeine, enantiomers), mathematical properties
  - Dalke, stress, Ehrlich-Rarey, and Tautobase benchmark tests retained

### Fixed
- `tanimoto_coefficient` / `overlap_coefficient` crashed with sparse count
  input from `circular_fingerprint_counts()` — now auto-detects format
- `counts_to_array()` only accepted `dict` — now accepts `list[tuple]` too
- `fingerprint()` raised ValueError for `kind='ecfp'/'fcfp'/'torsion'`
- Doc crashes: `result.tanimoto` → `result.overlap` (AttributeError),
  `r['rgroups']` → R-group dict keys (KeyError)
- Doc silent bugs: `bond_order_mode=` → `match_bond_order=`,
  removed invalid `solvent=/pH=/chem=` kwargs from `find_mcs` examples
- `decompose_rgroups` → `decompose_r_groups` in PYTHON.md and EXAMPLES.md
- All 196 `smsd.xxx()` references across docs verified against actual exports
- Replaced all stale API names (ecfp_counts, dice_similarity, smarts_search, etc.)

### Removed
- `test_reaction_aware.py` — reaction-aware MCS is not part of the public build

## [7.0.0] - 2026-04-13

### Summary
Major release: unified API, clean break from legacy aliases, full Java parity.

### Breaking Changes
- Removed `smsd.mcs()` — use `smsd.find_mcs()`
- Removed `smsd.substructure_search()` — use `smsd.find_substructure()`
- Removed `smsd.all_mcs()` — use `smsd.find_mcs(mol1, mol2, max_results=N)`
- Removed camelCase aliases: `overlapCoefficient`, `tanimoto`, `count_overlap_coefficient`, `count_tanimoto`

### Added
- Unified Python API: `find_mcs(mol1, mol2, max_results=1)` and `find_substructure(query, target, max_results=1)`
- Java convenience methods: `SearchEngine.findMCS(g1, g2)` and `SearchEngine.findSubstructure(query, target)` with MolGraph and IAtomContainer overloads
- Raw C++ bindings renamed to `_native_*` prefix (clearly internal)

### Changed
- All internal calls updated to unified API names
- `mcs_from_smiles()`, `mcs_rdkit()`, `substructure_rdkit()`, `depict_mcs()`, `depict_substructure()` use new API
- `__all__` cleaned of all deprecated entries
- `overlapCoefficient([], [])` returns 1.0 (trivially identical empty sets)

### Platforms
- macOS (arm64 Apple Silicon, x86_64), Linux (x86_64, aarch64), Windows (AMD64)
- GPU: Metal (Apple Silicon), CUDA (Volta+)
- Java 25+, C++17, Python 3.10-3.13

## [6.12.0] - 2026-04-07

### Summary
Correctness, performance, and API cleanup release.

### Included
- Fixed memory leak in SearchEngine cache
- Fixed overflow in graph-bound computation for large molecular graphs
- Added missing CIP Rule 3 (Z > E) per IUPAC 2013 in Java and C++
- Thread-safety: `volatile` on lazy-init fields in MolGraph
- Updated tautomer weights: nitroso-oxime 0.95, nitro-aci 0.95, pyridone 0.95
- Added selenium to tautomer compatibility, iodine to reaction-aware scoring
- Corrected SAH test SMILES (thioether, not ester connectivity)
- Relaxed formal charge matching in reaction-aware MCS
- Renamed `Mcs*` types to `MCS*`, `tanimoto` to `overlapCoefficient`
- Improved MCS construction throughput via faster compatibility graph traversal
- Reduced allocation pressure throughout the MCS pipeline
- Faster convergence on symmetric ring systems
- Faster substructure search domain initialisation
- Improved throughput on Apple Silicon with native vector operations
- Stage-aware pipeline routing to skip unnecessary MCS stages
- `MCSStageTimers` profiling API for pipeline diagnostics
- `TargetCorpus` and `batch_find_substructure()` Python APIs
- Reduced allocations per query; thread-local SMILES parsing
- SDF batch cap at 100K molecules

## [6.11.1] - 2026-04-04

### Bug Fixes
- **ECFP initial invariants**: corrected circular fingerprint atom invariants to
  include bond-order and mass contributions in both binary and count ECFP variants
  (C++ and Java)
- **Path fingerprint canonical hash**: corrected path fingerprint to use a single
  canonical hash direction, fixing bit density inflation
- **FCFP pyrrole-N misclassification**: aromatic nitrogen acceptor classification
  now uses direct hydrogen count, fixing incorrect non-acceptor assignment for
  pyridine-N (pyridine N has a free lone pair; pyrrole N does not)
- **Thread safety**: `prewarmGraph()` now initialises the pattern fingerprint before
  entering parallel regions, preventing data races on lazy-init fields
- **Dead code removal**: removed unused internal accumulator from binary ECFP path

## [6.11.0] - 2026-04-04

### Summary
Performance, precision, and depiction release: faster MCS engine, publication-quality
SVG renderer (ACS 1996 standard), comprehensive layout engine, 35+ new Python bindings.

### Included
- Core engine: cache-performance improvements on hot MCS computation paths (15-25%)
- Pre-indexed candidate sets in the MCS solver — eliminates repeated linear scans per
  frontier atom
- Publication-quality SVG depiction engine (ACS 1996 standard):
  - Jmol/CPK element colors, asymmetric double bonds, wedge/dash stereo bonds
  - Bond-to-label clipping, H-count subscripts, charge superscripts
  - Full customization via DepictOptions (bond_length, colors, fonts, sizes)
  - Side-by-side MCS pair rendering with atom-atom mapping numbers
- Multi-phase 2D layout pipeline: template match, ring-first, chain zig-zag, force
  refinement, overlap resolution, crossing reduction, canonical orientation,
  bond-length normalisation
- Distance-geometry 3D coordinate generation with iterative coordinate refinement
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
- Mode-matched benchmark leaderboards with explicit `defaults`, `strict`, and `ring-only` comparison modes
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
