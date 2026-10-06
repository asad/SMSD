# Publishing SMSD 7.2.1

Build and validate before publishing. The release uses Java 25 and CDK 2.13,
one CPU/OpenMP wheel per operating system for CPython 3.14, and a Python source
distribution. Keep one Python version across the wheel set:

| Wheel target | Runtime requirement |
|---|---|
| Linux x86_64 | CPython 3.14, glibc 2.28+ |
| macOS arm64 | CPython 3.14, macOS 26+ |
| Windows x86_64 | CPython 3.14, 64-bit Windows |

All three prepared wheels passed target-platform tests and strict collection
against one source archive. Linux execution used local emulation. GitHub 7.2.1
is released; PyPI and Maven publication remain pending. See
[7.2.1 validation](VALIDATION_7.2.1.md) for the complete
build and artifact record.
The [7.2.0 validation](VALIDATION_7.2.0.md) records historical execution results. Other
interpreters and architectures can build from source; they are outside this
wheel set. Each new wheel must pass its own installed-package tests; a
minimum deployment tag does not prove execution on that minimum OS version.

The GitHub asset set also includes the Java library, CLI, sources and Javadoc,
and C++ headers. The portable Java 25 package is shared across the three
operating systems. Java sources and artifacts live under `java/src/` and
`java/target/`; publish Maven artifacts with `mvn -f java/pom.xml deploy`.
The root Maven aggregator supports `mvn verify`. The root `pyproject.toml`
remains the single Python manifest for the C++ extension and Python package.

The retained 7.2.0 benchmarks use Python 3.13.14 with RDKit 2026.09.1. The
7.2.1 deadline regression check is recorded separately; no new cross-solver
benchmark ranking is claimed.
Keep `RESULTS_7.2.0.md` and the original benchmark archive name/hash when
including that historical evidence. RDKit is optional at runtime.

## Prepare locally

Use JDK 25, Maven, CMake, a C++17 compiler and Python 3.14. On macOS, OpenMP
requires an available `libomp`; the wheel repair bundles that library.
The preparation script defaults to a macOS 26 deployment target. If overriding
`MACOSX_DEPLOYMENT_TARGET`, all bundled libraries must support the chosen target.

Run from the SMSD checkout root:

```bash
python3.14 -m venv .venv-release
source .venv-release/bin/activate
python -m pip install --upgrade pip build scikit-build-core pybind11 \
  pytest pytest-timeout twine delocate rdkit==2026.03.6
SMSD_RELEASE_PYTHON="$PWD/.venv-release/bin/python" \
  SMSD_RELEASE_JOBS=2 bash scripts/prepare-release.sh

cd dist/release-7.2.1
shasum -a 256 -c SHA256SUMS
cd ../..
python -m twine check --strict \
  dist/release-7.2.1/smsd-7.2.1.tar.gz \
  dist/release-7.2.1/smsd-7.2.1-*.whl
```

The [benchmark guide](../benchmarks/README.md) gives reproduction commands.
The [7.2.0 report](../benchmarks/RESULTS_7.2.0.md) remains unchanged; do not
relabel its timings or quality cohorts as new 7.2.1 measurements. Record fresh
release checks, source identity and asset hashes in
[7.2.1 validation](VALIDATION_7.2.1.md).
The preparation script tests and assembles artifacts for the local platform;
it does not publish. It can include other tested wheels from a directory
specified by `SMSD_PLATFORM_WHEELS_DIR`. Build the Linux wheel locally in a
manylinux container and the Windows wheel on a native Windows host. Each
platform must pass native Debug tests and the full installed-wheel Python
suite. The local release preflight requires RDKit 2026.03.6 for interoperability
checks and active OpenMP; install both in the selected build/test environment.

The local preparation script builds the macOS wheel from the release source
distribution. Reuse that exact archive for the additional local Linux build;
Docker or Podman provides the manylinux environment:

```bash
python -m pip install cibuildwheel==4.2.1
python scripts/build_python_wheels.py build \
  --sdist dist/release-7.2.1/smsd-7.2.1.tar.gz \
  --platform linux --arch x86_64 \
  --output-dir build/platform-release/linux
```

Wheel output directories must be empty. Reuse that source distribution on a
native Windows machine, from Developer PowerShell:

```powershell
py -3.14 -m venv .venv-release
.\.venv-release\Scripts\python.exe -m pip install cibuildwheel==4.2.1
.\.venv-release\Scripts\python.exe scripts/build_python_wheels.py build --sdist dist/release-7.2.1/smsd-7.2.1.tar.gz --platform windows --arch AMD64 --output-dir build/platform-release/windows
```

The helper runs dependency repair, all native Debug suites and installed-wheel
tests. Keep its source and wheel provenance JSON alongside the build logs.
These native Windows commands are an alternative to the manual hosted workflow.
Windows repair uses `scripts/repair_windows_wheel.py` to select AMD64 Microsoft
release runtimes from Visual Studio redistributables or Windows System32.
It checks the DLL versions against the extension's linker family and verifies
that wheel repair bundled those exact files. Microsoft runtime terms are
included in the wheel's license metadata.

The release plan builds macOS and Linux locally and uses the manual
`python-publish.yml` workflow for native Windows validation. It defaults to
Windows only and does not publish unless explicitly requested. Use the same
exact 7.2.1 release commit for every build and retain the source manifests,
wheel provenance and test logs. The corrected 7.2.0 Windows run passed and is separate
historical evidence; it does not satisfy the 7.2.1 gate.

After downloading the tested Windows artifact, collect all three wheels and
check them against the prepared source distribution:

```bash
python scripts/collect-release-wheels.py \
  --release-dir dist/release-7.2.1 \
  --wheel-dir build/platform-release/linux \
  --wheel-dir build/platform-release/windows
python scripts/collect-release-wheels.py \
  --release-dir dist/release-7.2.1 --check-only
python -m twine check --strict \
  dist/release-7.2.1/smsd-7.2.1.tar.gz \
  dist/release-7.2.1/smsd-7.2.1-*.whl
```

Collection verifies package versions, CPython/ABI/platform tags, native binary
format and architecture, all wheel RECORD hashes, Python wrappers, installed
C++ headers and license copies against the release checkout and source distribution.
For Windows, the supported delvewheel loader is checked separately; the
application code must still match the source exactly. Collection
requires all three target wheels by default and rejects inconsistent builds.
This artifact check complements target-platform tests; it cannot establish
that a Windows wheel executes by inspecting it on macOS. Use
`--allow-incomplete` only during preparation, then run the strict check above
before publishing.

## PyPI

Use an account with permission to publish `smsd`. Create a project-scoped
token under [PyPI account settings](https://pypi.org/manage/account/#api-tokens)
and store it locally. Twine prompts for the token with hidden input; do not
put credentials in a command, repository file or commit.

Upload only the Python source distribution and the three tested wheels:

```bash
source .venv-release/bin/activate
python scripts/collect-release-wheels.py \
  --release-dir dist/release-7.2.1 --check-only
python -m twine upload --username __token__ \
  dist/release-7.2.1/smsd-7.2.1.tar.gz \
  dist/release-7.2.1/smsd-7.2.1-*.whl
```

Confirm [the PyPI release](https://pypi.org/project/smsd/7.2.1/), then install
in a separate environment:

```bash
SMSD_PYPI_CHECK="$(mktemp -d /tmp/smsd-pypi-check.XXXXXX)"
python3.14 -m venv "$SMSD_PYPI_CHECK"
"$SMSD_PYPI_CHECK/bin/python" -m pip install --no-cache-dir smsd==7.2.1
"$SMSD_PYPI_CHECK/bin/python" -c \
  'import smsd; print(smsd.__version__); assert smsd.is_substructure("CC", "CCC")'
```

An existing published filename cannot be replaced. If an upload is interrupted,
check the release files and upload only files that are missing. The
[PyPA publishing guide](https://packaging.python.org/en/latest/tutorials/packaging-projects/)
describes authentication and installation checks.

## Maven Central

The coordinates are `com.bioinceptionlabs:smsd:7.2.1`. The `release` profile
attaches sources and Javadoc, signs the artifacts, uploads through the Central
Publishing plugin, and waits for publication. Its `autoPublish=true` setting
publishes after Central validation succeeds.

Configure the verified `com.bioinceptionlabs` namespace in the
[Central Portal](https://central.sonatype.com/publishing/namespaces).
Store the Portal user-token name and password in the `central` server entry
of `~/.m2/settings.xml`, outside the repository:

```xml
<settings>
  <servers>
    <server>
      <id>central</id>
      <username>PORTAL_USER_TOKEN_NAME</username>
      <password>PORTAL_USER_TOKEN_PASSWORD</password>
    </server>
  </servers>
</settings>
```

Use an existing GPG signing key whose public key is discoverable by Central.
Run signing from your terminal so GPG can prompt for its passphrase:

```bash
chmod 600 ~/.m2/settings.xml
export GPG_TTY="$(tty)"
gpg --list-secret-keys --keyid-format LONG
mvn -f java/pom.xml -Prelease -Dslow.tests.exclude=nothing clean deploy
```

For multiple signing keys, select the intended key with
`-Dgpg.keyname=YOUR_SIGNING_KEY_ID`. Never place a passphrase on the command
line. Follow [Sonatype's Maven publishing instructions](https://central.sonatype.org/publish/publish-portal-maven/)
and [signing requirements](https://central.sonatype.org/publish/requirements/gpg/).

After publication, verify a clean Maven download:

```bash
SMSD_CENTRAL_CHECK="$(mktemp -d /tmp/smsd-central-check.XXXXXX)"
mvn -f java/pom.xml -B -Dmaven.repo.local="$SMSD_CENTRAL_CHECK" dependency:get \
  -Dartifact=com.bioinceptionlabs:smsd:7.2.1 \
  -DremoteRepositories=central::default::https://repo.maven.apache.org/maven2
```

## GitHub

Create the tag from the tested release commit on `master`, then upload the
prepared asset set. This does not require a hosted build:

```bash
git switch master
git pull --ff-only origin master
git tag -a v7.2.1 -m "SMSD 7.2.1"
git push origin v7.2.1
gh release create v7.2.1 --repo asad/SMSD --verify-tag \
  --title "SMSD 7.2.1" --notes-file docs/RELEASE_NOTES.md \
  dist/release-7.2.1/*
```

The GitHub 7.2.1 release is already published; the commands above record the
release process. GitHub publication can precede PyPI and Maven. For future
releases, run them after the fixes are merged and local checks pass. Compare
the tag's production inputs with the source used for validation.
Hosted release workflows remain manual; local publishing does not dispatch
them.
