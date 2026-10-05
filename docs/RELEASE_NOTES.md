# SMSD Pro 7.2.0 — Unreleased

This source targets 7.2.0. The reviewed baseline is the 7.1.2 source snapshot
`6807f31`; an earlier 7.1.2 tag or artifact does not contain these changes.
Java requires JDK 25 and uses CDK 2.13. C++ requires C++17.

## Correctness and resource behavior

- Tautomer matching preserves elements when atom-type matching is enabled.
  Tetrahedral matching uses normalized R/S configuration and mapped ligand
  parity, including traversal and CDK ligand-order changes.
- Signed-weight and bond objectives use objective-aware candidate, component,
  fragment and target selection. Java keeps double weight precision; native
  scoring retains its documented integer-millipoint range.
- Small objective searches use admissible score bounds. Bounded large searches
  can return valid incumbents without proving a global optimum.
- Canonical mapping selection computes orbit closure and reports incomplete
  generators or exhausted resources explicitly. Captured permutations preserve
  chemical properties, and skipped/capped generator searches are incomplete.
  Enumeration retains distinct
  raw mappings when symmetry work cannot complete.
- Java molecule-cache values are weak, preventing the cached graph from
  retaining its weak CDK key. Target exclusions are scoped to each call and
  cannot contaminate later cached matching domains.
- Constrained batches retain original target indices and chemistry, support
  programmatic Builder targets, select by the requested objective and respect
  per-pair timeouts without modifying caller options.
- Native recursive extensions retain frame-local assignments and branch
  buffers. Bounded connected seeds reserve search time and retain validated
  statin and taxane lower bounds.
- Coverage validation and recovery share the requested bond policy, including
  aromatic/Kekule matches under flexible aromaticity.

See [the changelog](../CHANGELOG.md) for API behavior changes and
[the algorithm review](ALGORITHM_REVIEW.md) for regression contracts and
reproduction commands.

## Benchmark reporting

The [current benchmark report](../benchmarks/RESULTS_7.2.0.md) records source
versions, dependencies, policies, input hashes, validity, mapping quality and
cancellation observations. Controlled curated comparisons use 10-second
budgets; full-corpus quality runs use common 1-second budgets. Conversion,
startup, SMARTS compilation and tautomer feature measurements are
reported separately. Different chemistry or mapping quality cannot support a
headline speedup.

Derived MoleculeNet pair collections are identified as such rather than being
presented as the original Dalke/Hastings corpus. The checked-in neighbor
collection contains self-source IDs, duplicate pairs and variable similarity;
its measurements remain separate from regenerated input.

## Build and release preparation

Use [the local preparation script](../scripts/prepare-release.sh) to validate
and assemble artifacts for this source. Platform wheels require their own
installed-package checks. These notes do not claim that current artifacts have
been published. The compact set uses Java 25/CDK 2.13 packages shared across
Linux, macOS and Windows, and CPU/OpenMP Python 3.14 wheels for Linux x86_64,
macOS arm64 and Windows x86_64, plus a source distribution. Intel macOS and
Linux arm64 use source builds. See [publishing commands](PUBLISHING.md)
for PyPI, Maven Central and GitHub.
The optional C++ RDKit adapter has a public header and exported CMake target;
the current RDKit headers require C++20. The core remains C++17.
All four Java JARs include SMSD's LICENSE and NOTICE: `META-INF/smsd` for
library, CLI and source JARs, and `doc-files/smsd` for Javadoc.
The platform checks pass all 12 native suites on macOS, emulated Linux and
Windows Server 2022, with 691 installed-wheel Python tests passed and 8 skips
on each. The Windows build uses MSVC 19.44 and CPython 3.14.7; its repaired
wheel bundles the C++ and OpenMP runtime DLLs. These execution checks do not
extend the macOS benchmark comparisons to other platforms.
The first Windows repair selected an older C++ runtime from the runner's
Java installation. Runtime selection is now explicit; the corrected wheel
requires a validation rerun before publication.

Earlier 7.1.2 test counts and primitive measurements are retained as a
[historical validation record](VALIDATION_7.1.2.md), rather than current results.
