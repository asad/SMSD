# Changelog

All notable changes to SMSD Pro are documented in this file.

## [6.12.2] - 2026-04-12

### Summary
Performance, quality, and API refinement release.

### Included
- Improved MCS result quality across a broader range of matching configurations
- New named matching profiles: `default`, `strict`, `pharma`, `reaction`
- `ringMatchesRingOnly` and `matchFormalCharge` default values updated for
  wider out-of-the-box applicability
- Fully configurable `mcs()` API — all matching flags available as keyword
  arguments without constructing a ChemOptions object
- `max_stage` parameter for explicit speed/quality trade-off control
- Similarity metrics: `overlap_coefficient` and `tanimoto_coefficient` added
- Java 25 language modernisation: internal codebase quality pass
- Build system fully aligned: Maven, CMake, Python, Docker (JDK 25)
- Security: Docker runtime image now runs as non-root user

## [6.12.0] - 2026-04-07

### Summary
Correctness, performance, and API cleanup release.

### Included
- Corrected memory management in search engine caches
- Fixed integer overflow in graph colouring bound for large molecules
- Added missing CIP Rule 3 (Z > E) per IUPAC 2013 in Java and C++
- Thread-safety improvements on lazily-initialised fields
- Updated tautomer weights for nitroso-oxime, nitro-aci, and pyridone forms
- Extended tautomer compatibility to additional heteroatom classes
- Corrected test molecule SMILES for SAH (thioether connectivity)
- Relaxed formal charge matching in reaction-aware MCS
- API naming alignment: `MCS*` types and `overlapCoefficient` metric
- Significant MCS throughput improvements on large molecules
- Reduced memory allocation in hot paths
- Performance improvements on Apple Silicon
- `TargetCorpus` and `batch_find_substructure()` Python APIs
- SDF batch processing with molecule count limit
- `ThreadLocal` SMILES parser for concurrent use

## [6.11.1] - 2026-04-04

### Bug Fixes
- **Circular fingerprint invariants**: added missing bond-order and mass-number
  invariants to ECFP variants in C++ and Java, aligning with the original
  published specification
- **Path fingerprint hashing**: corrected canonical hash to be
  direction-independent, eliminating spurious bit-density differences between
  forward and reverse path traversals
- **Aromatic nitrogen classification**: fixed hydrogen-count check used to
  distinguish acceptor from non-acceptor aromatic nitrogen in FCFP features
- **Thread safety**: resolved potential data race in fingerprint lazy-init
  cache when called from concurrent contexts

## [6.11.0] - 2026-04-04

### Summary
Performance, precision, and depiction release.

### Included
- Substantial MCS throughput improvements on medium and large molecules
- Publication-quality 2D structure depiction (ACS 1996 style):
  element colours, asymmetric double bonds, wedge/dash stereo bonds,
  H-count subscripts, charge superscripts, full DepictOptions customisation
- Side-by-side MCS pair rendering with atom-mapping annotations
- Full 2D layout engine: ring-first placement, chain extension, overlap
  resolution, canonical orientation, bond-length normalisation
- Distance geometry 3D coordinate generation with force-field refinement
- 40+ ring scaffold templates covering pharmaceutical and PAH cores
- Full 2D/3D coordinate transform suite
- 35+ new Python bindings
- Java: explicit per-atom type matching for exotic valence states
- Expanded chemistry and layout test coverage

## [6.10.2] - 2026-04-03

### Summary
Correctness release.

### Included
- Fixed connected-component filter for non-induced MCS mode
- Regression tests for challenging molecule pairs

## [6.10.1] - 2026-04-03

### Summary
Stability and correctness release.

### Included
- Hardened search for large molecules — iterative bounds replace recursive
  approaches susceptible to stack overflow on complex scaffolds
- Correct handling of duplicate target atoms in result repair
- All test assertions are now fully deterministic and machine-speed independent
- CI/CD stability fixes
- Python publish switched to manual dispatch only

## [6.9.0] - 2026-04-02

### Summary
Core chemistry correctness and I/O hardening release.

### Included
- Deterministic, direction-stable MCS results for asymmetric molecule pairs
- Consistent ring-matching semantics across C++, Python, and Java
- Native MDL MOL V2000 metadata preservation (name, comment, SDF properties)
- Native MDL MOL V3000 reader and writer
- R-group molfile support
- Improved SMILES and SMARTS attachment-point handling
- Stereo and CIP round-trip coverage improvements
- Python bindings for native mol-block read/write

## [6.8.1] - 2026-04-02

### Summary
Release alignment and packaging cleanup.

### Included
- Version alignment across Java, C++, and Python
- Documentation refresh

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
- Java 25 / C++17 / Python 3.10+ support
- CLI, SDF batch, and JSON export
