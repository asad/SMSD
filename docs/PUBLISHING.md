# Publishing SMSD 7.2.0

Build and validate before publishing. The release uses Java 25 and CDK 2.13,
one CPU/OpenMP wheel per operating system for CPython 3.14, and a Python source
distribution. Keep one Python version across the wheel set:

| Wheel target | Runtime requirement |
|---|---|
| Linux x86_64 | CPython 3.14, glibc 2.28+ |
| macOS arm64 | CPython 3.14, macOS 26+ |
| Windows x86_64 | CPython 3.14, 64-bit Windows |

See [validation](VALIDATION_7.2.0.md) for actual execution results. Other
interpreters and architectures can build from source; they are outside this
wheel set. The macOS wheel was executed on macOS 27.0.1, rather than separately
on its minimum deployment target.

The GitHub asset set also includes the Java library, CLI, sources and Javadoc,
and C++ headers. Benchmarks use Python 3.13.14 with RDKit 2026.09.1 to keep the
baseline and candidate on the same interpreter. RDKit is optional at runtime.

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
  pytest pytest-timeout twine delocate
SMSD_RELEASE_PYTHON="$PWD/.venv-release/bin/python" \
  SMSD_RELEASE_JOBS=2 bash scripts/prepare-release.sh

cd dist/release-7.2.0
shasum -a 256 -c SHA256SUMS
cd ../..
python -m twine check --strict \
  dist/release-7.2.0/smsd-7.2.0.tar.gz \
  dist/release-7.2.0/smsd-7.2.0-*.whl
```

Run the complete benchmark protocol in [the benchmark guide](../benchmarks/README.md).
Record the interpreter, dependency versions, input hashes, mapping validity
and budgets in [the results](../benchmarks/RESULTS_7.2.0.md) before release.
The preparation script tests and assembles artifacts for the local platform;
it does not publish. It can include other tested wheels from a directory
specified by `SMSD_PLATFORM_WHEELS_DIR`. Build the Linux wheel locally in a
manylinux container and the Windows wheel on a native Windows host. Each
platform must pass native Debug tests and the full installed-wheel Python
suite. Optional RDKit checks require RDKit in that test environment.

Build the source distribution once, then use it for each wheel. From a clean
release checkout on Linux or macOS with Docker available:

```bash
python -m pip install cibuildwheel==4.2.1
python scripts/build_python_wheels.py prepare \
  --output-dir build/platform-release/source --release-ref "$(git rev-parse HEAD)"
python scripts/build_python_wheels.py build \
  --sdist build/platform-release/source/smsd-7.2.0.tar.gz \
  --source-manifest build/platform-release/source/source-manifest.json \
  --platform linux --arch x86_64 \
  --output-dir build/platform-release/linux
```

Preparation and wheel output directories must be empty. Reuse that source
distribution on a native Windows machine, from Developer PowerShell:

```powershell
py -3.14 -m venv .venv-release
.\.venv-release\Scripts\python.exe -m pip install cibuildwheel==4.2.1
.\.venv-release\Scripts\python.exe scripts/build_python_wheels.py build --sdist build/platform-release/source/smsd-7.2.0.tar.gz --source-manifest build/platform-release/source/source-manifest.json --platform windows --arch AMD64 --output-dir build/platform-release/windows
```

The helper runs dependency repair, all native Debug suites and installed-wheel
tests. Keep its source and wheel provenance JSON alongside the build logs.
These native Windows commands are an alternative to the manual hosted workflow.

The manual `python-publish.yml` workflow defaults to Windows only and does not
publish unless explicitly requested. Its platform selector avoids rebuilding
locally validated Linux and macOS wheels. Use an exact release commit or tag;
keep its build and test logs with the artifacts. Confirm any hosted execution
under the repository's local execution policy before dispatching it.

After downloading the tested Windows artifact, collect all three wheels and
check them against the prepared source distribution:

```bash
python scripts/collect-release-wheels.py \
  --release-dir dist/release-7.2.0 \
  --wheel-dir build/platform-release/linux \
  --wheel-dir build/platform-release/windows
python scripts/collect-release-wheels.py \
  --release-dir dist/release-7.2.0 --check-only
python -m twine check --strict \
  dist/release-7.2.0/smsd-7.2.0.tar.gz \
  dist/release-7.2.0/smsd-7.2.0-*.whl
```

Collection verifies package versions, CPython/ABI/platform tags, native binary
format and architecture, all wheel RECORD hashes, Python wrappers, installed
C++ headers and license copies against the release checkout and source distribution. It
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
  --release-dir dist/release-7.2.0 --check-only
python -m twine upload --username __token__ \
  dist/release-7.2.0/smsd-7.2.0.tar.gz \
  dist/release-7.2.0/smsd-7.2.0-*.whl
```

Confirm [the PyPI release](https://pypi.org/project/smsd/7.2.0/), then install
in a separate environment:

```bash
SMSD_PYPI_CHECK="$(mktemp -d /tmp/smsd-pypi-check.XXXXXX)"
python3.14 -m venv "$SMSD_PYPI_CHECK"
"$SMSD_PYPI_CHECK/bin/python" -m pip install --no-cache-dir smsd==7.2.0
"$SMSD_PYPI_CHECK/bin/python" -c \
  'import smsd; print(smsd.__version__); assert smsd.is_substructure("CC", "CCC")'
```

An existing published filename cannot be replaced. If an upload is interrupted,
check the release files and upload only files that are missing. The
[PyPA publishing guide](https://packaging.python.org/en/latest/tutorials/packaging-projects/)
describes authentication and installation checks.

## Maven Central

The coordinates are `com.bioinceptionlabs:smsd:7.2.0`. The `release` profile
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
mvn -Prelease -Dslow.tests.exclude=nothing clean deploy
```

For multiple signing keys, select the intended key with
`-Dgpg.keyname=YOUR_SIGNING_KEY_ID`. Never place a passphrase on the command
line. Follow [Sonatype's Maven publishing instructions](https://central.sonatype.org/publish/publish-portal-maven/)
and [signing requirements](https://central.sonatype.org/publish/requirements/gpg/).

After publication, verify a clean Maven download:

```bash
SMSD_CENTRAL_CHECK="$(mktemp -d /tmp/smsd-central-check.XXXXXX)"
mvn -B -Dmaven.repo.local="$SMSD_CENTRAL_CHECK" dependency:get \
  -Dartifact=com.bioinceptionlabs:smsd:7.2.0 \
  -DremoteRepositories=central::default::https://repo.maven.apache.org/maven2
```

## GitHub

Create the tag from the tested release commit on `master`, then upload the
prepared asset set. This does not require a hosted build:

```bash
git switch master
git pull --ff-only origin master
git tag -a v7.2.0 -m "SMSD 7.2.0"
git push origin v7.2.0
gh release create v7.2.0 --repo asad/SMSD --verify-tag \
  --title "SMSD 7.2.0" --notes-file docs/RELEASE_NOTES.md \
  dist/release-7.2.0/*
```

Run these commands after the fixes are merged and the package publication
checks succeed. Compare the tag commit with the commit used for validation.
Hosted release workflows remain manual; local publishing does not dispatch
them.
