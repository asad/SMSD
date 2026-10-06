# SMSD 7.2.2 validation

Release preparation is in progress. Previous release results remain in
[7.2.1 validation](VALIDATION_7.2.1.md); they are not new 7.2.2 results.

| Check | Status |
|---|---|
| Java CML/PDB and existing IO regressions | 41 passed, including 16 new file-format cases |
| Full Java 25 suite | 1,258 passed, 15 opt-in skips, zero failures/errors; 5m11s |
| Four Java JARs and launchers | Passed source, class-version, licence and version checks |
| macOS arm64 DMG | Passed native mount/copy/run/remove/detach; ten CLI cases |
| macOS arm64 Python wheel | Pending frozen-source build |
| Linux x86_64 DEB | Passed install/run/remove under local QEMU; ten CLI cases |
| Linux x86_64 Python wheel | Pending frozen-source build |
| Windows x86_64 Python wheel and MSI | Pending native Windows execution |
| Three-wheel and three-installer collection | Pending |
| GitHub publication and download verification | Pending |
| PyPI and Maven Central | Pending |

The installer runtime is Eclipse Temurin 25.0.4.1+1. Local Linux execution uses
QEMU on an arm64 macOS host. macOS and Windows installers are unsigned;
macOS notarisation is unavailable. The macOS Java runtime declares macOS 11
as its minimum version; this does not establish execution on that version.
The macOS Python wheel requires macOS 26 or later.

No new cross-solver benchmark claim is made. The
[7.2.0 benchmark report](../benchmarks/RESULTS_7.2.0.md) retains its measured scope.

The NAD redox fixture validates a 35-atom common-core witness separately from
the five-second search result. This fixes a machine-dependent size expectation;
production search code and budgets are unchanged. The final full suite passed.
