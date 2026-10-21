# SMSD Pro 8.0.0 validation

This record lists the checks run for the SMSD Pro 8.0.0 release, whose
packages are built from the source tag `v8.0.0`.

| Check | Result |
| --- | --- |
| Source package builds and passes the native test suites | Passed |
| C++ native suites, release and sanitiser builds on macOS arm64 | Passed |
| C++ native suites with GCC 14 on Linux | Passed |
| Java library and CLI suites on Java 25 | Passed |
| Java library and CLI suites on a Java 8 runtime | Passed |
| macOS arm64 wheel: native suites from the source package and Python suite on the installed wheel | Passed |
| Linux x86_64 and Windows x86_64 wheels: native suites from the source package and Python suite on the installed wheel | Passed |
| Every wheel built from the published source package `smsd-8.0.0.tar.gz` | Passed |
| Documented examples | Passed |
| Release checksums in `SHA256SUMS` | Published with the release |
