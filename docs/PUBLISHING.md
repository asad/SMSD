# Publishing SMSD 7.2.0

Build and validate locally before publishing. The release uses Java 25 and
CDK 2.13, one CPU/OpenMP wheel for CPython 3.14 on macOS arm64, and a Python
source distribution. Other interpreters and platforms can build from source;
this local validation does not certify them.
The macOS wheel targets macOS 26 or later and was executed on macOS 27.0.1.

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
The preparation script tests and assembles artifacts; it does not publish.

## PyPI

Use an account with permission to publish `smsd`. Create a project-scoped
token under [PyPI account settings](https://pypi.org/manage/account/#api-tokens)
and store it locally. Twine prompts for the token with hidden input; do not
put credentials in a command, repository file or commit.

Upload only the Python source distribution and wheel:

```bash
source .venv-release/bin/activate
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
