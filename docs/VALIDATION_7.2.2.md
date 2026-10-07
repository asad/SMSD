# SMSD 7.2.2 validation

All three platform builds and strict collection passed. GitHub, PyPI and Maven
Central publication are tracked separately. Previous release results remain in
[7.2.1 validation](VALIDATION_7.2.1.md); they are not new 7.2.2 results.

Version 7.2.2 now targets Java 8 or later, with Java 25 LTS preferred for
builds and bundled installers. Both Java compatibility suites passed.

| Check | Status |
|---|---|
| Full Java 8 suite | 1,278 passed, 15 opt-in skips, zero failures/errors; 5m16s |
| Full Java 25 suite | 1,278 passed, 15 opt-in skips, zero failures/errors; 5m15s |
| Public value-type and file-format contracts | Passed on both runtimes, including 17 new value-contract cases |
| Four Java JARs | Passed Java 8 bytecode, source and licence checks |
| Frozen source archive | All 224 source files match the clean source commit; one generated metadata file |
| Shaded CLI JAR on Java 8 | Seven startup/search cases passed, including compact and pretty Unicode SDF JSON |
| CLI JSON under Windows-style encodings | Ten Java 25 cases passed without a JVM encoding workaround; native MSI cases also passed |
| Unix source and generated launchers | Passed on both runtimes; Windows execution pending |
| macOS arm64 DMG | Passed all ten CLI checks, native installation/removal and 58 arm64 binary checks |
| macOS arm64 Python wheel | 12 native Debug suites on unchanged inputs (reused evidence); 691 fresh Python tests passed, 8 optional skips |
| Linux x86_64 DEB | Passed all ten CLI checks, extraction, installation and removal under QEMU |
| Linux x86_64 Python wheel | 12 native Debug suites on unchanged inputs (reused evidence); 691 fresh Python tests passed, 8 optional skips |
| Docker Linux arm64 and x86_64 images | Ten CLI checks per image, exact release JAR, non-root execution and licence checks passed locally; x86_64 uses QEMU |
| Windows x86_64 MSI | Passed all ten CLI checks, native installation and removal; unsigned |
| Windows x86_64 Python wheel | 12 native Debug suites and 691 installed-package tests passed, 8 optional skips; proof reused for byte-identical payloads after metadata refresh |
| Three-wheel and three-installer collection | Passed for the complete platform set |
| GitHub publication and download verification | Published; all 25 assets downloaded anonymously and matched their checksums |
| PyPI and Maven Central | Pending |

The installer runtime is Eclipse Temurin 25.0.4.1+1. Local Linux execution uses
QEMU on an arm64 macOS host. macOS and Windows installers are unsigned;
macOS notarisation is unavailable. The macOS application passes ad hoc signature
checks, but Gatekeeper rejects it. Quarantined downloads remain untested.
The macOS Java runtime declares macOS 11
as its minimum version; this does not establish execution on that version.
The macOS Python wheel requires macOS 26 or later. Its fresh installed-wheel
checks use CPython 3.14.8 and RDKit 2026.03.6 on macOS 27.0.1 arm64; native
Debug checks take 99.09 seconds; fresh tests of the metadata-refreshed installed
wheel take 12.01 seconds.
Metal and CUDA are disabled; bundled OpenMP is required. Native Debug evidence
is reused after independently verifying identical C++/Python/build inputs across
the packaging and Java updates. Final distributions refresh documentation and
package metadata only; all tested binaries, Python wrappers, headers and licences
remain byte-identical. The final source archive differs only in documentation and generated package
metadata.
Manifests retain the original test identities and make reused evidence explicit.

Docker images provide the Java CLI with Temurin 25.0.4.1+1. Both Linux
architectures run in a local arm64 Colima VM; x86_64 uses QEMU. Each passes
the same ten semantic cases with the exact final shaded JAR. Tests run with
UID/GID 10001, no network and a read-only container filesystem. Image archives
retain project licences, 184 runtime legal files and the vendor notice;
the existing pinned OpenJDK source archive matches both runtime identities.
Container registry publication remains pending.

Linux wheel checks use AlmaLinux 8.10/glibc 2.28, GCC 14.2.1, CPython 3.14.5
and RDKit 2026.03.6 under local QEMU x86_64. The unchanged native Debug suite
took 335.40 seconds; fresh tests of the metadata-refreshed installed wheel
took 34.69 seconds. Release
compilation retains `-O3` and OpenMP. Link-time optimisation is disabled after
GCC crashes under emulation; compile and link commands confirm the setting.

Native Windows checks ran on Windows Server 2022 AMD64 with MSVC
19.44.35229.0, CPython 3.14.7 and RDKit 2026.03.6. All 12 Debug suites passed
in 235.87 seconds; installed-wheel tests passed 691 cases with 8 optional skips
in 8.60 seconds. The final Windows wheel refreshes only its description and
RECORD; native and installed-package proof is reused after verifying all 65 other
members, including their compressed payloads, remain identical. No native Windows
execution of the metadata-refreshed archive is claimed. Repair bundles the checked
Microsoft 14.44.35211.0 release runtimes; all three PE binaries are AMD64. The MSI passed installation, all ten
semantic cases and removal. See [the native run](https://github.com/asad/SMSD/actions/runs/37557118227).

The Java tests use the same compiled classes on macOS arm64 with Temurin
25.0.4.1+1 and Corretto 8.504.04.1 (`1.8.0_504-b04`). All 53 production classes
and 169 test classes target class-file version 52. The packaged sources match
the checkout, and all four JARs retain the project's LICENSE and NOTICE.

Packaged Java classes match the full dual-runtime test build. A separate consumer
declaring only SMSD passes six SMARTS cases on both runtimes and receives CDK
SMARTS transitively. The updated CLI passes compact and pretty JSON regressions
with Unicode file paths and a Windows-1252 output stream; standard output remains
open after writing. The bundled-runtime verifier reads JVM diagnostics as UTF-8.

No new cross-solver benchmark claim is made. The
[7.2.0 benchmark report](../benchmarks/RESULTS_7.2.0.md) retains its measured scope.

The NAD redox fixture validates a 35-atom common-core witness separately from
the five-second search result. This fixes a machine-dependent size expectation;
search algorithms and budgets are unchanged.
