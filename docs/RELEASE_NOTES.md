# SMSD 7.2.2

- Fix Java CLI CML and PDB input.
- Write CLI JSON output as UTF-8, including Unicode file paths on Windows.
- Include the SMARTS dependency for Maven library users.
- Target Java 8 or later; prefer Java 25 LTS for builds and execution.
- Use Java 8-compatible value classes; see the [changelog](https://github.com/asad/SMSD/blob/master/CHANGELOG.md) for record API changes.
- Reject empty files and files containing multiple molecules/models; use SDF for batch targets.
- Add Windows MSI, macOS DMG and Linux DEB installers with bundled Java 25.
- Refresh Java, C++ and Python examples.

Python wheels target CPython 3.14 on Windows x86_64, Linux x86_64 and macOS arm64.
Java uses CDK 2.13; the C++ core requires C++17.
Java CLI Docker archives target Linux x86_64 and arm64.

Use `SHA256SUMS` to verify downloads. macOS and Windows installers are unsigned.

All three platform checks passed. See the GitHub release and package indexes
for publication status.
