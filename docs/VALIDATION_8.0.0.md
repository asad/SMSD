# SMSD Pro 8.0.0 validation

This report lists the checks required for the SMSD Pro 8.0.0 release packages
and records their results. An entry is completed only with a result tied to the
final source commit and the checksums of the published files. A pending entry
is not a passing result.

## Scope

The report covers the source distribution and its source manifest; the Java
library, shaded CLI, sources and Javadoc JARs; the portable Java CLI archive;
the C++ header archive; the three CPython 3.14 wheels; and the Windows, macOS
and Linux Java CLI installers. Release wheels run batch operations on the CPU
with OpenMP. Metal and CUDA are optional source-build backends and are outside
the scope of this report. Minimum platform versions are listed in the
[installation guide](https://github.com/asad/SMSD/blob/master/docs/HOWTO-INSTALL.md#platform-requirements).

## Status

| Check | Status |
|---|---|
| Frozen source distribution and source manifest | Pending |
| Full Java suite on Java 25 | Pending |
| Java 8 runtime compatibility suite and CLI checks | Pending |
| Java library, shaded CLI, sources and Javadoc JARs | Pending assembly and identity checks |
| Portable Java CLI archive and launchers | Pending |
| C++ headers, native Debug suites and installed CMake consumer | Pending |
| Windows x86_64 (AMD64) CPython 3.14 wheel | Pending native Debug and installed-wheel checks |
| Linux x86_64 (AMD64) CPython 3.14 wheel | Pending native Debug and installed-wheel checks |
| macOS arm64 CPython 3.14 wheel | Pending native Debug and installed-wheel checks |
| Agreement of the Python README, source PKG-INFO and wheel METADATA | Pending |
| Complete wheel collection and build provenance | Pending |
| Windows x86_64 (AMD64) MSI | Pending build, signing status, installation, CLI and removal checks |
| macOS arm64 DMG | Pending build, signing status, installation, CLI and removal checks |
| Linux x86_64 (AMD64) DEB | Pending build, installation, CLI and removal checks |
| Complete installer collection and bundled-runtime licences | Pending |
| Java CLI container image archives, if included | Pending architecture, reload and CLI checks |
| Documented examples against the final distributions | Pending |
| Final asset checksums | Pending |
| GitHub publication and public download verification | Not published |
| PyPI publication and clean installation | Not published |
| Maven Central publication and clean dependency resolution | Not published |

## Java compatibility

The Java module targets class-file version 52 (Java 8). Compilation with
`--release 8` and execution on a Java 8 runtime are separate checks, and both
are recorded. The installers bundle Eclipse Temurin Java 25.0.4.1+1. For each
installer, the report records the runtime and runtime-source archives,
licences, architecture, signing status and execution environment.

## Recording rules

Native and emulated execution are recorded separately, together with the
operating-system versions tested. A minimum platform version that was not
tested is identified as untested. A result on one operating system, Python ABI
or installer does not establish the behaviour of another.
