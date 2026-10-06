# Changelog

All notable changes to SMSD Pro are documented in this file.

## [7.2.1] - Unreleased

### Fixed
- Native MCS seed searches check the shared deadline before each candidate
  extension. Already-expired seeds and orientation probes return before setup.
  This fixes a Linux regression without relaxing its timeout assertion.

### Repository and release packaging
- Moved Java sources, resources and launchers into `java/src/`, with its Maven
  module at `java/pom.xml` and build artifacts under `java/target/`. The root
  Maven aggregator supports `mvn verify`; publishing uses the Java module.
- Kept C++ under `cpp/` and Python under `python/`, with shared release scripts,
  documentation and licenses at the root. The root `pyproject.toml` remains
  the single manifest for the Python package and C++ extension.
- Restricted the Docker build context to Maven manifests, Java sources and
  legal files, excluding generated API pages, test reports and build artifacts.
- Updated current version examples and release artifacts to 7.2.1. The compact
  release targets portable Java 25 packages and CPython 3.14 CPU/OpenMP wheels
  for Linux x86_64, macOS arm64 and Windows x86_64, plus a source distribution.
- Frozen-source local macOS and emulated Linux wheels each pass all 12 native
  Debug suites and 691 Python tests with 8 optional skips. Native Windows
  execution, full three-wheel collection and publication are pending. Track
  results in `docs/VALIDATION_7.2.1.md`.
- Adjusted Java test wall-clock guards to 35 seconds around a 30-second drug
  pair search and 12 seconds around a 10-second pharmacophore search. Search
  budgets and result assertions are unchanged.
- Carried forward the reviewed 7.2.0 search and chemistry fixes. The benchmark
  report, fingerprints, archive names and measured numbers remain 7.2.0
  evidence. The deadline regression check is separate from those benchmarks;
  no new cross-solver performance ranking is claimed.

## [7.2.0] - 2026-10-05

Published on GitHub. Maven Central and PyPI remain at 7.1.1. The historical
validation and benchmark report retain their 7.2.0 source scope; native Windows
builds were checked separately after publication. Version 7.2.1 is in preparation.

### Fixed
- Tautomer matching preserves element identity and other requested chemistry
  constraints. Relative tetrahedral tags are compared using R/S and mapped
  ligand parity, including SMILES ring-opening neighbor order.
- Signed and bond objectives drive bounds, fragment selection, enumeration and
  constrained target choice. Native non-induced domain refinement preserves
  valid partial mappings and shared search deadlines.
- Native recursive extensions keep undo lists and candidate buffers local to
  each call, preventing inconsistent assignments and invalid bond indices.
- Coverage repair applies the requested bond policy consistently, preserving
  valid aromatic/Kekule matches under flexible aromaticity.
- The optional C++ RDKit adapter now exposes a public header and an installed
  CMake target. Its build enables the implementation and propagates C++20 for
  the current RDKit headers.
- Java library, CLI, source and Javadoc JARs include SMSD's LICENSE and NOTICE
  in project-specific directories, alongside dependency notices.
- Exact symmetry canonicalization uses molecular automorphism generators;
  incomplete generators, orbit caps or expiry now raise explicit errors.
  Enumeration retains raw keys when symmetry proof is unavailable.
- Python RDKit wrappers, multi-result searches and batches preserve original
  atom indices and caller options. Conversion caches detect molecule edits,
  use weak keys and retain metadata for live converted graphs.
- Python batch substructure passes its timeout independently of thread count.
  Invalid weighted batch options raise before OpenMP workers start.
- Python progress reporting performs one native search and forwards the final
  result. It does not currently report intermediate stages.
- The standalone C++ depiction header uses a portable C++17 pi constant rather
  than depending on a platform-specific `M_PI` macro.
- Native MOL/SDF file APIs interpret filenames as UTF-8 on Windows. MSVC
  builds use UTF-8 source and executable character sets.
- Release wheel collection validates the Windows DLL loader added by repair,
  while retaining exact source checks for the application code.
- Windows wheel repair selects compatible Microsoft runtimes explicitly,
  avoiding older DLLs from unrelated applications in the build host's PATH.
  Wheels include the Microsoft runtime license documents.

### Optimised
- Core Python batch and compiled SMARTS multi-target bindings retain graph
  references and their Python owners instead of copying graph caches.
- Batch prewarming includes ring systems before parallel workers start.
- Bounded connected seeds and clique-stage allocation reserve time for native
  recovery, including valid statin and taxane lower-bound fixtures.
- Source build options can override Metal/CUDA auto-detection without duplicate
  CMake arguments. Local release comparisons build CPU-only wheels explicitly.

### Validation and documentation
- Added independent small-graph objective/enumeration oracles, molecular
  symmetry checks, lifetime/index tests and randomized stereo traversals.
- Reworked benchmark protocols around explicit chemistry, objective, budget
  and witness validation. The current comparison uses RDKit 2026.09.1.
- Removed unsupported README speed, quality and dataset provenance claims.
  See `benchmarks/RESULTS_7.2.0.md` for measurements and limits.
- Limited Python wheel builds to CPython 3.14 across Linux x86_64, macOS arm64
  and Windows x86_64, plus a source distribution. Added installed-wheel checks
  and local PyPI and Maven Central publishing commands.
- All 12 native Debug suites and 691 installed-wheel Python tests pass on
  macOS, emulated Linux and native Windows Server 2022, with 8 Python skips
  for optional Open Babel bindings and opt-in external benchmarks.

### Earlier source fixes included in this release
- Standardized `MCS` capitalization in Java/C++ APIs, classes, helper names,
  benchmarks and tests. Java callers should use `setMCSTimeoutMs`,
  `findMCSSmiles`, `findMCSSmarts` and `batchMCSConstrained`; tuning fields use
  `nearMCSDelta` and `nearMCSCandidates`. Python snake_case APIs are unchanged.
- Java benchmark sources use the current timeout field and public graph APIs,
  restoring compilation while retaining their scoring calculations.
- Python bindings use CMake's modern `FindPython` module discovery, removing
  pybind11's CMP0148 deprecation warnings and unnecessary embedding-library
  discovery on Unix.
- MCS fast paths enforce induced, connectivity and fragment constraints;
  weighted searches preserve query atom indices. Full-graph degree sequences
  no longer undercut bounds for partial induced matches.
- Reversed Java and C++ MCS results reapply the original query's complete-ring
  and fragment filters before being compared with the incumbent.
- Java McSplit and product-graph pruning preserve valid partial mappings;
  recursive extensions retain their candidate buffers. MCS enumeration
  applies defaults and validity filters, and similarity bounds honor options.
- Java and C++ connected MCS components respect weighted/bond objectives;
  candidates are postfiltered before promoting incumbents. C++ MCS enumeration
  preserves induced constraints on containment candidates.
- Expired search budgets stay expired during recursive unwinding. Java clique
  search checks deadlines around coloring and pivot work; orientation probes,
  recovery and retries share the caller's total budget.
- Python automatic MCS selection routes unsupported lightweight constraints
  to the native solver; explicit lightweight requests reject unsupported
  options instead of silently ignoring them.
- C++ substructure enumeration retains symmetric self mappings; disconnected
  cycles plus paths cannot enter the linear-path shortcut.
- C++ unrestricted bond order preserves strict aromaticity. CPU/GPU candidate
  domains consistently enforce isotope, chirality, ring and tautomer options.
- Native clique helper seeds and extensions preserve query bonds, select the
  largest connected component and honor result caps and expired deadlines.
- Rectangular assignment rejects ragged and nonfinite inputs instead of
  reading outside rows or failing to terminate.
- C++ weighted MCS rejects nonfinite and out-of-range millipoint scores,
  including partial-score overflow hidden by cancellation.

### Earlier source optimisations included in this release
- C++ matcher setup reuses sorted query neighbors and skips unused target
  canonicalization. Connected-component postprocessing traverses adjacency.
- Rectangular assignment avoids square padding, using
  `O(min(m,n)² max(m,n))` time and `O(m+n)` auxiliary space.
- General matching seeds a valid maximal matching before blossom augmentation,
  preserving maximum cardinality while reducing augmenting-search work.
- Valid connected MCS seeds can extend despite different full neighborhoods;
  Java adds a bounded anchor stage for medium graphs before expensive search.

## [7.1.2] - 2026-10-04

### Updated
- Java dependency: CDK 2.12 to the latest stable CDK 2.13; standardisation
  now uses `Aromaticity.Model.Daylight`.
- Jackson databind 2.20.0 to 2.21.7, covering the patched-version requirements
  of all 11 current repository dependency advisories.
- Version metadata aligned across Java, C++, Python, CLI, and citation.

### Fixed
- Java domain cache isolation across mutable chemical matching options,
  permissive/tautomer pruning, strict aromatic bond matching, and enumeration
  of targets with more than 4,096 candidates.
- Java telemetry defaults, mapping validation bounds, and timeout overflow.
- Java and C++ multi-hop substructure pruning under extra target edges.
- C++ maximum-clique result cap and incumbent tie handling, disconnected-query
  target reuse, policy-incompatible fingerprint pruning, and small-matcher budgets.
- C++ MCS discarding a larger validated directional result.
- C++ non-bipartite aromatic kekulization (including azulene) and undefined
  signed overflow in graph hashing; existing canonical output is preserved.
- Broken Java cage test fixtures and impossible historical C++ MCS thresholds;
  replacement lower bounds have explicit conserved-substructure witnesses.
- CMake package discovery, C++17/OpenMP propagation, and license installation.
- Incomplete nested Python source packaging, PowerShell launcher path, GPU
  build script portability, and skipped CPU batch validation.
- Source launchers and isolated generated distributions; C++ assertions
  remain active in every test build configuration.

### Optimised
- C++ clique pivot scanning removes temporary allocations and adds a safe
  branch bound. Neighborhood construction traverses only the local frontier.

### Release preparation
- Added a local build/test/packaging script with source-distribution builds,
  installed-wheel validation, macOS library repair, and SHA-256 checksums.
- Hosted workflows now require manual dispatch; release artifacts are built
  locally before preparing the GitHub draft.

## [7.1.1] - 2026-04-14

Bug-fix patch on top of v7.1.0.  No new features, no public API breakage.

### Fixed
- Cross-language ECFP / FCFP fingerprint parity between Java, C++, and
  Python.  Two long-standing Java drifts at radius ≥ 1 (signed vs
  unsigned neighbour-hash sort, and pharmacophore implicit-H count
  for pyrrole-type nitrogen) now produce bits byte-identical to the
  C++ / Python reference.
- Java canonical SMILES writer: bond symbol for aromatic-adjacent
  single bonds, and implicit H count inside stereo brackets.
- Python `smsd.canonical_smiles(smi)` / `smsd.to_smiles(smi)` raised
  `TypeError` on string input.  Both now accept `str` or `MolGraph`.
- `MatchResult.overlapCoefficient` returned the wrong similarity metric.
  Both `MatchResult.overlap` and `MatchResult.overlapCoefficient` now
  return Szymkiewicz-Simpson overlap as documented; the new
  `MatchResult.tanimoto` attribute exposes the Jaccard value.
  `__repr__` shows both.
- Canonical SMILES writer now emits `[nH]` for pyrrole-type aromatic
  nitrogen, so SMSD output kekulizes cleanly in downstream readers
  (pyrrole, indole, carbazole, fused benzo-pyrrole).
- FP-level `smsd.overlapCoefficient` and `smsd.count_overlapCoefficient`
  camelCase aliases now return Simpson overlap as documented;
  `smsd.tanimoto_coefficient` / `smsd.count_tanimoto_coefficient`
  expose Jaccard.  Empty-vs-empty convention aligned at 1.0 (trivially
  identical) for all similarity helpers.

### Removed
- Dead C++ fingerprint shim headers under `cpp/include/fp/`.  The
  canonical C++ fingerprint API is `smsd::batch::detail::*` in
  `cpp/include/smsd/batch.hpp`, documented in `docs/CPP.md`.

### Verified
- Python pytest: 603 passed, 6 skipped, 0 failed.
- Java JUnit: 581 tests, 0 failures, 0 errors.
- Canonical SMILES reference set matches character-for-character
  across Java, C++, and Python.

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
- Added selenium to tautomer compatibility, iodine to scoring
- Corrected SAH test SMILES (thioether, not ester connectivity)
- Relaxed formal charge matching in the default MCS profile
- Standardized maximum common substructure types as `MCS*`, and renamed `tanimoto` to `overlapCoefficient`
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
