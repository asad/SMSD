# Publishing SMSD 7.2.2

Build and validate before publishing. The Java JAR targets Java 8 with CDK 2.13;
builds and native installers use Java 25 LTS. Validate on both runtimes. Prepare
one CPU/OpenMP wheel per operating system for CPython 3.14, and a Python source
distribution. Keep one Python version across the wheel set:

| Wheel target | Runtime requirement |
|---|---|
| Linux x86_64 | CPython 3.14, glibc 2.28+ |
| macOS arm64 | CPython 3.14, macOS 26+ |
| Windows x86_64 | CPython 3.14, 64-bit Windows |

Version 7.2.2 is in preparation. Use
[7.2.2 validation](VALIDATION_7.2.2.md) for current checks;
[7.2.1 validation](VALIDATION_7.2.1.md) is historical evidence.

Each release also includes the Java library, CLI, sources and Javadoc,
C++ headers, and one MSI, DMG and DEB with bundled Java. Follow
[installer preparation](INSTALLERS.md) for packaging and installation checks.
Java publication uses `mvn -f java/pom.xml deploy`; the root Maven POM is an
aggregator. The root `pyproject.toml` builds the C++ extension and Python package.

The retained 7.2.0 benchmarks use Python 3.13.14 with RDKit 2026.09.1. The 7.2.1 deadline fix is retained; no new cross-solver benchmark ranking is claimed.
Keep `RESULTS_7.2.0.md` and the original benchmark archive name/hash when
including that historical evidence. RDKit is optional at runtime.

## Prepare locally

Use JDK 25, Maven, CMake, a C++17 compiler and Python 3.14. On macOS, OpenMP
requires an available `libomp`; the wheel repair bundles that library.
Set `SMSD_RELEASE_JAVA8_HOME` to an installed Java 8 JDK for the compatibility
test run. Compilation uses Java 25 with `--release 8`.
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

cd dist/release-7.2.2
shasum -a 256 -c SHA256SUMS
cd ../..
python -m twine check --strict \
  dist/release-7.2.2/smsd-7.2.2.tar.gz \
  dist/release-7.2.2/smsd-7.2.2-*.whl
```

The [benchmark guide](../benchmarks/README.md) gives reproduction commands.
The [7.2.0 report](../benchmarks/RESULTS_7.2.0.md) remains unchanged; do not
relabel its timings or quality cohorts as new 7.2.2 measurements. Record fresh
release checks, source identity and asset hashes in
[7.2.2 validation](VALIDATION_7.2.2.md).
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
  --sdist dist/release-7.2.2/smsd-7.2.2.tar.gz \
  --platform linux --arch x86_64 \
  --output-dir build/platform-release/linux
```

Wheel output directories must be empty. Reuse that source distribution on a
native Windows machine, from Developer PowerShell:

```powershell
py -3.14 -m venv .venv-release
.\.venv-release\Scripts\python.exe -m pip install cibuildwheel==4.2.1
.\.venv-release\Scripts\python.exe scripts/build_python_wheels.py build --sdist dist/release-7.2.2/smsd-7.2.2.tar.gz --platform windows --arch AMD64 --output-dir build/platform-release/windows
```

The helper runs dependency repair, all native Debug suites and installed-wheel
tests. Keep its source and wheel provenance JSON alongside the build logs.
These native Windows commands are an alternative to the manual hosted workflow.
Windows repair uses `scripts/repair_windows_wheel.py` to select AMD64 Microsoft
release runtimes from Visual Studio redistributables or Windows System32.
It checks the DLL versions against the extension's linker family and verifies
that wheel repair bundled those exact files. Microsoft runtime terms are
included in the wheel's license metadata.

Build macOS and Linux locally. The manual `installers.yml` workflow builds
and tests the Windows MSI and Python wheel together from the existing draft
release inputs. The `python-publish.yml` workflow remains available for
Python-only builds and explicit PyPI publication.

After downloading the tested Windows artifact, collect all three wheels and
check them against the prepared source distribution:

```bash
python scripts/collect-release-wheels.py \
  --release-dir dist/release-7.2.2 \
  --wheel-dir build/platform-release/linux \
  --wheel-dir build/platform-release/windows
python scripts/collect-release-wheels.py \
  --release-dir dist/release-7.2.2 --check-only
python -m twine check --strict \
  dist/release-7.2.2/smsd-7.2.2.tar.gz \
  dist/release-7.2.2/smsd-7.2.2-*.whl
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

### Publish through GitHub

Configure the `smsd` project's Trusted Publisher on PyPI:

| Field | Value |
|---|---|
| Owner | `asad` |
| Repository | `SMSD` |
| Workflow filename | `python-publish.yml` |
| Environment | `pypi` |

Open [Build and validate Python wheels](https://github.com/asad/SMSD/actions/workflows/python-publish.yml).
Choose **Run workflow**, branch `master`, `release_tag=v7.2.2`,
`platform=all` and `publish=true`. This rebuilds and tests all three platforms
before publishing. Configure the publisher and finish the release checks first.
See [PyPI's Trusted Publisher setup](https://docs.pypi.org/trusted-publishers/adding-a-publisher/).

Equivalent shell command, once ready:

```bash
gh workflow run python-publish.yml --repo asad/SMSD --ref master \
  -f release_tag=v7.2.2 -f platform=all -f publish=true
gh run list --repo asad/SMSD --workflow python-publish.yml --limit 3
```

### Upload the locally validated files

Use an account with permission to publish `smsd`. Create a project-scoped
token under [PyPI account settings](https://pypi.org/manage/account/#api-tokens)
and store it locally. Twine prompts for the token with hidden input; do not
put credentials in a command, repository file or commit.

Upload only the Python source distribution and the three tested wheels:

```bash
source .venv-release/bin/activate
python scripts/collect-release-wheels.py \
  --release-dir dist/release-7.2.2 --check-only
python -m twine upload --username __token__ \
  dist/release-7.2.2/smsd-7.2.2.tar.gz \
  dist/release-7.2.2/smsd-7.2.2-*.whl
```

Confirm [the PyPI release](https://pypi.org/project/smsd/7.2.2/), then install
in a separate environment:

```bash
SMSD_PYPI_CHECK="$(mktemp -d /tmp/smsd-pypi-check.XXXXXX)"
python3.14 -m venv "$SMSD_PYPI_CHECK"
"$SMSD_PYPI_CHECK/bin/python" -m pip install --no-cache-dir smsd==7.2.2
"$SMSD_PYPI_CHECK/bin/python" -c \
  'import smsd; print(smsd.__version__); assert smsd.is_substructure("CC", "CCC")'
```

An existing published filename cannot be replaced. If an upload is interrupted,
check the release files and upload only files that are missing. The
[PyPA publishing guide](https://packaging.python.org/en/latest/tutorials/packaging-projects/)
describes authentication and installation checks.

## Maven Central

The coordinates are `com.bioinceptionlabs:smsd:7.2.2`. The `release` profile
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
  -Dartifact=com.bioinceptionlabs:smsd:7.2.2 \
  -DremoteRepositories=central::default::https://repo.maven.apache.org/maven2
```

## GitHub

Prepare an authenticated draft release after local checks and source freeze.
Upload the validated JAR, source archive, source manifest and checksum list
so the Windows job uses the same inputs. Keep the release private while native
Windows checks remain pending.

After all three wheels and all three installers pass collection, regenerate
`SHA256SUMS`, upload the complete asset set and publish the draft. Verify the
public downloads against local checksums. Release notes should contain only
changes, downloads and platform requirements.
