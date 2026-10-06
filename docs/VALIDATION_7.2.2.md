# SMSD 7.2.2 validation

Release preparation is in progress. Previous release results remain in
[7.2.1 validation](VALIDATION_7.2.1.md); they are not new 7.2.2 results.

Version 7.2.2 now targets Java 8 or later, with Java 25 LTS preferred for
builds and bundled installers. Both Java compatibility suites passed.

| Check | Status |
|---|---|
| Full Java 8 suite | 1,276 passed, 15 opt-in skips, zero failures/errors; 5m9s |
| Full Java 25 suite | 1,276 passed, 15 opt-in skips, zero failures/errors; 5m9s |
| Public value-type and file-format contracts | Passed on both runtimes, including 17 new value-contract cases |
| Four Java JARs | Passed Java 8 bytecode, source and licence checks |
| Frozen source archive | All 224 source files match the clean source commit; one generated metadata file |
| Shaded CLI JAR on Java 8 | Five startup/search cases passed |
| Unix source and generated launchers | Passed on both runtimes; Windows execution pending |
| macOS arm64 DMG | Passed all ten CLI checks, native installation/removal and 58 arm64 binary checks |
| macOS arm64 Python wheel | 12 native Debug suites; 691 Python tests passed, 8 optional skips |
| Linux x86_64 DEB | Passed all ten CLI checks, extraction, installation and removal under QEMU |
| Linux x86_64 Python wheel | Pending frozen-source build |
| Windows x86_64 Python wheel and MSI | Pending native Windows execution |
| Three-wheel and three-installer collection | Pending |
| GitHub publication and download verification | Pending |
| PyPI and Maven Central | Pending |

The installer runtime is Eclipse Temurin 25.0.4.1+1. Local Linux execution uses
QEMU on an arm64 macOS host. macOS and Windows installers are unsigned;
macOS notarisation is unavailable. The macOS application passes ad hoc signature
checks, but Gatekeeper rejects it. Quarantined downloads remain untested.
The macOS Java runtime declares macOS 11
as its minimum version; this does not establish execution on that version.
The macOS Python wheel requires macOS 26 or later. Its fresh installed-wheel
checks use CPython 3.14.8 and RDKit 2026.03.6 on macOS 27.0.1 arm64; native
Debug checks take 99.09 seconds and Python checks take 7.75 seconds.
Metal and CUDA are disabled; bundled OpenMP is required.

The Java tests use the same compiled classes on macOS arm64 with Temurin
25.0.4.1+1 and Corretto 8.504.04.1 (`1.8.0_504-b04`). All 53 production classes
and 167 test classes target class-file version 52. The packaged sources match
the checkout, and all four JARs retain the project's LICENSE and NOTICE.

No new cross-solver benchmark claim is made. The
[7.2.0 benchmark report](../benchmarks/RESULTS_7.2.0.md) retains its measured scope.

The NAD redox fixture validates a 35-atom common-core witness separately from
the five-second search result. This fixes a machine-dependent size expectation;
search algorithms and budgets are unchanged.
