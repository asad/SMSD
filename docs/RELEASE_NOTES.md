# SMSD Pro 7.1.2

Java now uses [CDK 2.13](https://github.com/cdk/cdk/releases/tag/cdk-2.13),
the latest stable Chemistry Development Kit release verified on 4 October 2026.
This patch fixes search correctness and release packaging without changing
the public search APIs. Java requires JDK 25; C++ requires C++17.

## Search fixes

- Java cached atom domains now include the selected chemical matching options.
  Changing charge, isotope, aromaticity, or other constraints cannot reuse an
  incompatible cached domain.
- Java permissive atom, aromaticity, and tautomer searches use conservative
  screening. Enumerating targets larger than 4,096 atoms no longer truncates
  candidates; unrestricted bond order still respects strict aromaticity.
- Java and C++ multi-hop pruning counts cumulative neighborhoods so extra
  target edges cannot incorrectly reject a valid substructure.
- C++ result limits cap stored maximum cliques rather than prematurely ending
  the search. Incumbent-size ties remain discoverable; duplicate product edges
  and self-loops are ignored.
- C++ disconnected query components now share one injective mapping instead
  of reusing target atoms. Fingerprint screening only applies features that
  remain valid under the selected matching policy.
- C++ MCS keeps the best validated result from either search direction rather
  than shrinking a larger result to force agreement.
- Small C++ matchers honor time budgets. Java and C++ accept very large
  timeouts without overflow.
- Java telemetry overloads handle default options consistently, and invalid
  mapping indices return validation errors instead of throwing.
- Java standardisation uses CDK's current Daylight aromaticity API.
- C++ kekulization supports non-bipartite aromatic demand graphs such as
  azulene, using general maximum matching only when the faster bipartite path
  does not apply. Graph hashes use defined unsigned arithmetic while retaining
  their previous bit patterns and canonical ordering.

## Optimisation

C++ clique pivot selection scans bitsets directly, removing temporary pivot
arrays. An admissible branch bound skips searches that cannot improve the
result. Cumulative neighborhood construction visits the local frontier instead
of repeatedly scanning the entire graph.

## Packaging

- Java, C++, Python, CLI, and citation versions are aligned at 7.1.2.
- Installed CMake packages propagate C++17 and discover OpenMP when required;
  package configuration remains valid after moving the installation.
- Python builds use the root `pyproject.toml`; the redundant nested metadata
  that produced incomplete source distributions was removed.
- C++ archives and CMake installations include `LICENSE` and `NOTICE`.
- Portable GPU build scripts and the PowerShell repository path are corrected.
- Source launchers locate the current fat JAR, while generated launcher
  distributions are isolated under `target/appassembler/`.
- C++ test assertions remain enabled in Release builds.
- `scripts/prepare-release.sh` tests locally, builds the Python wheel from its
  source distribution, repairs macOS dependencies, and generates checksums.
- GitHub build and release workflows require manual dispatch. Tag pushes
  no longer launch hosted release matrices automatically.

This release is prepared as a GitHub draft. Assets include Java library, CLI,
source and Javadoc JARs, a launcher archive, C++ headers, Python source distribution,
a CPython 3.13 macOS 26+ arm64 wheel with bundled OpenMP and its license,
and SHA-256 checksums. Additional platform wheels and native
installers require separate builds; PyPI and container publishing are separate.

## Local validation

- Java: 1,192 passed, 15 opt-in benchmark cases skipped, no failures/errors.
- C++: all six CPU suites passed with assertions active; the Metal batch suite
  also passed on an Apple M5 GPU.
- Python: 603 passed, 6 optional cases skipped, no failures, testing the installed
  repaired wheel built from its source distribution.
- New general matching passed exhaustive graph oracles, and focused native
  AddressSanitizer/UndefinedBehaviorSanitizer checks passed.
- Relocated CMake consumers passed with OpenMP enabled and disabled.

See `docs/VALIDATION_7.1.2.md` for reproduction commands and measured primitive
optimisations. CUDA, native Windows execution, other Python wheel platforms,
and native installers were not validated in this local release preparation.

Apache 2.0 — see `NOTICE`.
