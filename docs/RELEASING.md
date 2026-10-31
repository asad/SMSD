# Preparing release assets

Before freezing, set the package, native, Java, CUDA and citation versions together:

```bash
python scripts/release_metadata.py set-version VERSION --date YYYY-MM-DD
```

This also sets the citation's release date without changing the cited article date.
Set the same date in the version's `CHANGELOG.md` heading, in place of `Unreleased`.
Release assembly extracts the selected version's text from `docs/RELEASE_NOTES.md`
into the asset `RELEASE_NOTES.md`, suitable for the GitHub release description.
Version 8.1.0 is published on GitHub and PyPI; its Java library and command-line
tool stay at 8.0.3 on Maven Central. PyPI uploads need the full schema 2 wheel
collection described below. Releases in the 8.x line keep the deprecated APIs.

Use a dedicated CPython 3.14 environment. Release tooling pins pip 26.2.2
before wheel builds and installed-wheel checks, including CPython 3.12.
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
CPython 3.12–3.15 on Linux x86_64/aarch64, macOS arm64 and Windows x86_64.

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

macOS wheels default to deployment target 14.0. The compiler, Python runtime
and every bundled library must support that target; changing the wheel tag
alone does not establish compatibility. Delocate checks bundled binaries.

CPU collection writes `WHEEL_MANIFEST.json` from checked build results.
The default schema 2 mode requires full native Debug and installed Python tests
for each wheel; PyPI publication, including 8.1.0, needs this mode. For
GitHub-only collection, the optional lean mode instead shares
the completed frozen-source checks and runs each platform's build, repair, audit
and installed-wheel smoke checks. It does not repeat the full suites per wheel.

After those shared checks pass, prepare a JSON summary containing exactly:

```json
{
  "source_commit": "FULL_FROZEN_COMMIT",
  "sdist_sha256": "CANONICAL_SDIST_SHA256",
  "cpp_release": "passed",
  "cpp_asan": "passed",
  "java8": "passed",
  "java25": "passed",
  "python": "passed",
  "examples": "passed",
  "check_source": "passed"
}
```

Use the exact commit and canonical archive hash, and record only completed checks.
Pass `--lean-validation frozen-validation.json` to each `build` or `verify-local`
command and to `collect-release-wheels.py`. When using `prepare-release.sh`, set
`SMSD_RELEASE_LEAN_VALIDATION` to that same completed summary. The script validates
its frozen-source identity before building and passes it to verification and
collection. Leaving this variable unset preserves the default full mode. Local
verification checks the native host, interpreter, source bytes, installed smoke
and exact installed wheel bytes; its proof records `installed_wheel_identity`
and `installed_wheel_smoke` only.
Retain the local build, repair and audit results separately. The summary is embedded
in schema 3, and every wheel proof binds its digest. Successful cibuildwheel builds
record `platform_recipe_checks` and `installed_wheel_smoke`. Each mode rejects the
other mode's check field and legacy per-wheel full-suite fields. Collection still
requires all 16 targets unless explicitly collecting an incomplete set.
`--check-only` reads either schema without the original summary file. The current
PyPI publication checker accepts schema 2 only; schema 3 is for GitHub collection
and is rejected by that checker.

Final assembly with `prepare-release.sh` requires both `SMSD_RELEASE_CUDA_WHEEL`
and `SMSD_RELEASE_CUDA_PROVENANCE` before building, so a completed collection
contains both manifests. Use the individual collectors for partial preparations.

Collect a CUDA wheel and its completed checks with:

```bash
python scripts/collect_cuda_wheels.py --release-dir RELEASE_DIRECTORY \
  --wheel CUDA_WHEEL --provenance CUDA_BUILD_PROOF
```

The collector checks its source, package contents, licences, recorded checks
and hashes, writes `CUDA_WHEEL_MANIFEST.json`, and updates `SHA256SUMS`.
Set `SMSD_RELEASE_CUDA_WHEEL` and `SMSD_RELEASE_CUDA_PROVENANCE` to include
this step in `prepare-release.sh`.

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

## Publishing to other channels

Publish each channel from the frozen commit and the validated release assets.
Signing keys, tokens and sign-in details stay in the publisher's local
configuration and are never added to the repository or the release. Each
hosted workflow run and each upload needs its own approval.

### Release description

The `RELEASE_NOTES.md` asset links to other files by paths relative to `docs/`.
Before using it as the GitHub release description, replace each relative link
with an absolute link to the same file at the release tag, such as
`https://github.com/asad/SMSD/blob/vVERSION/docs/CPP.md`, so that the links
work on the release page.

### Maven Central

1. Attach the Java assets written by `prepare-release.sh` to the GitHub release:
   `smsd-VERSION.jar`, `smsd-VERSION-sources.jar`, `smsd-VERSION-javadoc.jar`
   and `smsd-cli-VERSION-jar-with-dependencies.jar`.
2. From the same commit and toolchain, deploy the parent, library and
   command-line modules:

   ```bash
   mvn -f java/pom.xml -Prelease deploy
   ```

   The `release` profile attaches sources and Javadoc, signs each file and
   uploads one deployment to the Central Portal, where it is validated and then
   held. Nothing is published automatically.
3. Compare the SHA-256 of each of the four JARs with the file of the same name
   in `java/target/central-publishing/central-bundle.zip`. Publish only if all
   four match.
4. Publish the deployment in the Central Portal. Once it is available, check
   that a new project resolves `com.bioinceptionlabs:smsd:VERSION` and
   `com.bioinceptionlabs:smsd-cli:VERSION`.

### Conda

1. Build and test each package as described in the
   [conda packages guide](../packaging/conda/README.md).
2. Upload the `smsd` and `smsd-cuda` packages to the `test` label of the
   `bioinception` channel on anaconda.org, signed in with an account that may
   upload to that channel:

   ```bash
   anaconda upload --user bioinception --label test \
     /path/to/local-channel/*/smsd-VERSION-* \
     /path/to/local-channel/*/smsd-cuda-VERSION-*
   ```

3. On each platform, install `smsd` from that label in a new environment and
   check the version:

   ```bash
   conda create -n smsd-check -c bioinception/label/test -c conda-forge \
     "python=3.12" "smsd=VERSION"
   conda run -n smsd-check python -c "import smsd; print(smsd.__version__)"
   ```

   Check `smsd-cuda` from the same label on a Linux x86-64 host with an NVIDIA
   driver, as the conda packages guide describes.
4. Only when both packages pass, move them to the `main` label, which
   `conda install -c bioinception -c conda-forge smsd` uses:

   ```bash
   anaconda move --from-label test --to-label main bioinception/smsd/VERSION
   anaconda move --from-label test --to-label main \
     bioinception/smsd-cuda/VERSION
   ```

### Docker image archives

1. Build the CPU images for amd64 and arm64 and the CUDA image for amd64 as
   described in the [Docker section](HOWTO-INSTALL.md#docker) of the
   installation guide.
2. Save each image with `scripts/archive_docker_image.py` as
   `smsd-VERSION-cpu-linux-amd64.tar.xz`, `smsd-VERSION-cpu-linux-arm64.tar.xz`
   and `smsd-VERSION-cuda-linux-amd64.tar.xz`.
3. In the release directory, add the three archives to `SHA256SUMS`, once only:

   ```bash
   shasum -a 256 smsd-VERSION-*-linux-*.tar.xz >> SHA256SUMS
   ```

4. On a host that did not build the images, check the archives against
   `SHA256SUMS`, then load each one and run the commands in the Docker section.
   Run the arm64 image on an arm64 host and the CUDA image on a host with an
   NVIDIA GPU.

   ```bash
   shasum -a 256 --check --ignore-missing SHA256SUMS
   ```

5. Attach the three archives and the updated `SHA256SUMS` to the GitHub release.

### Native installers

`scripts/build_java_installer.py` builds the installer for the host it runs on:
the DMG on macOS arm64 and the DEB on Linux amd64. The manual `installers.yml`
workflow builds the native Windows packages for an existing release.
