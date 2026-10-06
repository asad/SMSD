# SMSD Pro 7.2.1 — Unreleased

Version 7.2.1 separates Java, C++ and Python source modules and updates release
packaging, and fixes native MCS seed deadline checks. Publication is pending.
It carries forward the reviewed 7.2.0
search fixes against the 7.1.2 source snapshot `6807f31`; an earlier 7.1.2 tag
or artifact does not contain those changes. Java requires JDK 25 and uses
CDK 2.13. The native core requires C++17.

## Repository layout

- Java sources, resources and launchers move into `java/src/`; Maven uses
  `java/pom.xml` and writes artifacts to `java/target/`. The root aggregator
  supports `mvn verify`; direct module builds and publication use
  `mvn -f java/pom.xml`.
- C++ remains under `cpp/`, and Python sources/tests remain under `python/`.
  Shared scripts, documentation and licenses stay at the root.
- The root `pyproject.toml` remains the single Python manifest, building the
  native extension from `cpp/` and packaging `python/smsd/`.

Java test wall-clock guards allow 35 seconds for a 30-second drug pair search
and 12 seconds for a 10-second pharmacophore search. Search budgets and result
assertions are unchanged; fresh execution results belong to the release gates
below.

Native seeds check the shared deadline before each candidate extension, and
expired seeds and orientation probes skip setup. The existing Linux deadline
regression retains its 100 ms wall-clock assertion for a 5 ms search budget.

## Search fixes carried forward from 7.2.0

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

The [7.2.0 benchmark report](../benchmarks/RESULTS_7.2.0.md) retains its original
source versions, measured numbers, fingerprints, input hashes and archive
names. The 7.2.1 deadline regression check is recorded separately in the
validation record; no new cross-solver benchmark ranking is claimed.
That report records validity, mapping quality and
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
Java packages include SMSD's LICENSE and NOTICE: `META-INF/smsd` for library,
CLI and source JARs, and `doc-files/smsd` for Javadoc.

The release plan uses local macOS and Linux builds, a native GitHub Windows
build, and verified collection of three wheels from the same source. Each
7.2.1 platform passed all 12 native Debug suites and 691 installed-wheel
Python tests with 8 optional skips. Strict three-wheel collection passed.
Fresh platform checks are tracked in
[7.2.1 validation](VALIDATION_7.2.1.md). The corrected 7.2.0 Windows runtime
build passed, but remains [historical evidence](VALIDATION_7.2.0.md), not a
7.2.1 result. Platform execution checks do not extend macOS benchmark timings
to other operating systems.

Earlier 7.1.2 test counts and primitive measurements are retained as a
[historical validation record](VALIDATION_7.1.2.md), rather than current results.
