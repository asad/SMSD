# Preparing release assets

Use a dedicated CPython 3.14 environment. Release tooling pins pip 26.2.2
before wheel builds and installed-wheel checks, including CPython 3.10.
These are build environment pins, not requirements for package users.

Freeze a clean commit after all changes and checks are complete. Prepare the
canonical source distribution once, then retain its archive and manifest:

```bash
RELEASE_COMMIT="$(git rev-parse HEAD)"
CANONICAL_SOURCE="/path/to/new-source-directory"
python scripts/build_python_wheels.py prepare \
  --output-dir "$CANONICAL_SOURCE" --release-ref "$RELEASE_COMMIT" --no-isolation
```

The preparation environment needs `pip==26.2.2`, `build==1.6.1` and `twine==7.0.0`; wheel
matrix drivers use `cibuildwheel==4.2.1`. Install dependencies before using
`--no-isolation`. Do not regenerate the archive for individual targets.
Build each declared wheel with `scripts/build_python_wheels.py build`, using
the same `--sdist` and `--source-manifest`. Its declared matrix contains
CPython 3.10–3.14 on Linux x86_64/aarch64, macOS arm64 and Windows x86_64.

For local Java and native checks and final wheel collection, use the same
clean checkout and canonical source:

```bash
SMSD_RELEASE_COMMIT="$RELEASE_COMMIT" \
SMSD_RELEASE_SDIST="$CANONICAL_SOURCE/smsd-VERSION.tar.gz" \
SMSD_RELEASE_SOURCE_MANIFEST="$CANONICAL_SOURCE/SOURCE_MANIFEST.json" \
SMSD_PLATFORM_WHEELS_DIR="/path/to/validated-matrix-wheels" \
SMSD_RELEASE_JAVA8_HOME="/path/to/java8-jdk" \
bash scripts/prepare-release.sh
```

Replace `VERSION` with the frozen version. The local script requires JDK 25
for compilation and a Java 8 JDK for compatibility checks. It checks the
commit, source archive and manifest before building and before collection.
Missing matrix targets fail collection; local results are retained without
being labelled a complete release. Build directories, unrepaired wheels and
failed staging directories are retained and their paths are printed.

Windows MSI, macOS DMG and Linux DEB validation remain separate requirements.
Preparing or collecting assets does not publish, tag or upload them.

The default Java test run includes `AlgorithmTest` and `StressTest`. Tests
tagged `SlowTest` are excluded unless explicitly selected; see the
[Java build instructions](../java/README.md#build-and-test).

## Publishing validated Python distributions

The manual `python-publish.yml` workflow has separate publishing paths.
For `smsd`, select `distribution=smsd`, `platform=all` and `publish=true`.
For the Linux CUDA add-on, select `distribution=smsd-cuda`, `platform=linux`
and `publish=true`. Both use the `pypi` environment and trusted publishing;
each PyPI project needs its own GitHub publisher registration.

The CUDA path uploads only the unchanged published
`smsd_cuda-VERSION-py3-none-manylinux_2_34_x86_64.whl`. Its release must also
contain the canonical source distribution, `SOURCE_MANIFEST.json`, `SHA256SUMS` and
`CUDA_WHEEL_MANIFEST.json`. The latter contains exactly:

- `schema_version`: `1`;
- `source_commit` and `sdist_sha256`: the canonical source identities;
- `smsd_version`: the matching release version;
- `wheel` and `wheel_sha256`: the validated add-on wheel identities;
- `installed_wheel_tests` and `native_gpu_tests`: `"passed"`, taken from
  completed checks of that exact build.

The workflow checks hashes, frozen source, metadata, licences and the build
manifest before staging the wheel, then rechecks the staged bytes immediately
before upload. The core source distribution excludes `packaging/**`; CUDA package metadata,
Python sources and licences are checked separately against Git blobs from the
same frozen commit and the release checkout. The CUDA wheel proof still binds
that commit and the canonical source distribution hash. Missing checks or an
incomplete CPU matrix are not repaired by the CUDA publishing path. PyPI
publisher registration must be complete before requesting a publishing run.
