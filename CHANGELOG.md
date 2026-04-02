# Changelog

All notable changes to SMSD Pro are documented in this file.

## [6.9.1] - 2026-04-02

### Summary
Python binding improvements, stronger aromaticity perception, and release
alignment across code, packaging, and documentation.

### Included
- Improved Python bindings for native graph preparation and aromaticity workflows
- Added explicit aromaticity controls with `AromaticityModel`, `perceive_aromaticity()`, `kekulize()`, and `dearomatize()`
- Improved native aromaticity perception for Kekule inputs, fused aromatic systems, and ring-sensitive matching paths
- Added shared Java/Python aromaticity parity coverage for key aromatic systems and edge cases
- Fixed MCS timeout budgeting and strengthened self-match fast paths
- Updated release-facing metadata, install snippets, and documentation to `6.9.1`

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
