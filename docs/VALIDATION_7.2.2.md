# SMSD 7.2.2 validation

Release preparation is in progress. Previous release results remain in
[7.2.1 validation](VALIDATION_7.2.1.md); they are not new 7.2.2 results.

| Check | Status |
|---|---|
| Java CML/PDB and existing IO regressions | 41 passed, including 16 new file-format cases |
| Java 25 suite | 1,258 distinct passes and 15 skips across full and focused runs; final full run pending |
| Four Java JARs and launchers | Passed source, class-version, licence and version checks |
| macOS arm64 Python wheel and DMG | Pending frozen-source build |
| Linux x86_64 Python wheel and DEB | Pending frozen-source build |
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

The first full run passed 1,257 cases, skipped 15 opt-in cases and failed one
bounded-search size expectation. Released 7.2.1 and current code produced the
same valid result under the same budget. The revised fixture preserves a
validated 35-atom common-core witness and separately checks the timed search
result. All 13 affected cases passed; no production search code or budget changed.
