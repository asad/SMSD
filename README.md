<p align="center">
  <a href="https://github.com/asad/SMSD" aria-label="SMSD Pro">
    <img src="icons/icon.svg" alt="SMSD Pro" width="180"/>
  </a>
</p>

<h1 align="center">SMSD Pro</h1>
<p align="center"><strong>Substructure and MCS Search for Chemical Graphs</strong></p>

<p align="center">
  <a href="https://central.sonatype.com/artifact/com.bioinceptionlabs/smsd"><img src="https://img.shields.io/maven-central/v/com.bioinceptionlabs/smsd" alt="Maven Central"/></a>
  <a href="https://pypi.org/project/smsd/"><img src="https://img.shields.io/pypi/v/smsd" alt="PyPI"/></a>
  <a href="LICENSE"><img src="https://img.shields.io/badge/License-Apache%202.0-blue.svg" alt="Licence"/></a>
  <a href="https://github.com/asad/SMSD/releases"><img src="https://img.shields.io/github/v/release/asad/SMSD" alt="Release"/></a>
</p>

SMSD provides substructure and maximum common substructure (MCS) search for
**Java**, **C++** and **Python**. It also includes fingerprints, molecular
standardisation, stereo handling, depiction and R-group decomposition.
Java integrates with **CDK 2.13**; C++ is a header-only C++17 library; Python
binds the native library and supports optional RDKit integration.

The current [GitHub release is 7.2.2](https://github.com/asad/SMSD/releases/tag/v7.2.2).
It supports Java 8 or later, fixes CML/PDB input and preserves UTF-8 CLI JSON.
Downloads include Java packages, C++ headers, Python wheels and native installers
for Windows, macOS and Linux. [Python 7.2.2 is also on PyPI](https://pypi.org/project/smsd/7.2.2/);
Maven Central publication remains pending.

## Installation

| Package | Requirements |
| --- | --- |
| Java library or portable CLI | Java 8 or later; Java 25 LTS preferred |
| Windows MSI | Windows AMD64; includes Java 25 LTS |
| macOS DMG | macOS arm64; includes Java 25 LTS |
| Linux DEB | Debian-compatible Linux AMD64; includes Java 25 LTS |
| Windows Python wheel | CPython 3.14, AMD64, Windows 10 or later |
| macOS Python wheel | CPython 3.14, arm64, macOS 26 or later |
| Linux Python wheel | CPython 3.14, x86_64, glibc 2.28 or later |
| C++ library | C++17 compiler |

The MSI, DMG and DEB install the **terminal Java CLI**. Python wheels install the
Python package separately. Windows and macOS installers are unsigned; the macOS
application has an ad-hoc signature. See [installer instructions](docs/INSTALLERS.md)
for installation, launcher paths and platform details.

### Java CLI

Download the JAR with its dependencies, then run it:

```bash
curl -LO https://github.com/asad/SMSD/releases/download/v7.2.2/smsd-7.2.2-jar-with-dependencies.jar
java -jar smsd-7.2.2-jar-with-dependencies.jar --version
java -jar smsd-7.2.2-jar-with-dependencies.jar \
  --Q SMI --q "c1ccccc1" --T SMI --t "c1ccc(O)cc1" --json -
```

This searches for benzene in phenol. Add `--mode mcs` for MCS search, or use
`--help` for all options. The CLI also reads MOL, MOL2, CML and PDB files and
supports SDF batch targets. See the [Java guide](docs/JAVA.md).

### Java library

Use these Maven coordinates for 7.2.2:

```xml
<dependency>
  <groupId>com.bioinceptionlabs</groupId>
  <artifactId>smsd</artifactId>
  <version>7.2.2</version>
</dependency>
```

Until 7.2.2 is available on Maven Central, install it from this source checkout
with `mvn -f java/pom.xml install` using Java 25. The downloadable JAR also
provides the complete classpath for the following example.

Save this as `MCSExample.java`:

```java
import com.bioinception.smsd.core.ChemOptions;
import com.bioinception.smsd.core.SMSD;
import java.util.Map;

public class MCSExample {
    public static void main(String[] args) throws Exception {
        SMSD matcher = new SMSD("c1ccccc1", "c1ccc(O)cc1", new ChemOptions());
        Map<Integer, Integer> mapping = matcher.findMCS(false, true, 1000);
        if (!matcher.isSubstructure() || mapping.size() != 6) {
            throw new IllegalStateException("Expected a six-atom benzene match");
        }
        System.out.println("MCS atoms: " + mapping.size());
    }
}
```

Compile with Java 25 and run on Java 8 or later:

```bash
javac --release 8 -cp smsd-7.2.2-jar-with-dependencies.jar MCSExample.java
java -cp "smsd-7.2.2-jar-with-dependencies.jar:." MCSExample
```

When compiling with Java 8, omit `--release 8`. On Windows, use
`"smsd-7.2.2-jar-with-dependencies.jar;."` as the classpath.
See the [Java module](java/README.md) for CDK container examples.

### Python

Install Python 7.2.2 from PyPI:

```bash
python -m pip install smsd==7.2.2
python -c "import smsd; print(smsd.__version__)"
```

You can also download the matching wheel from the
[GitHub release](https://github.com/asad/SMSD/releases/tag/v7.2.2) and install
that file with `python -m pip install`.

```python
import smsd

query = "c1ccccc1"
target = "c1ccc(O)cc1"
assert smsd.is_substructure(query, target)
mapping = smsd.find_mcs(query, target)
assert len(mapping) == 6
print(len(mapping))
```

Release wheels use CPU and OpenMP. Optional CUDA and Apple Metal support
requires a suitable source build. See [Python installation and examples](python/README.md)
and the [Python API guide](docs/PYTHON.md).

### C++

Save this as `example.cpp` in a source checkout:

```cpp
#include <smsd/smsd.hpp>

int main() {
    auto query = smsd::parseSMILES("c1ccccc1");
    auto target = smsd::parseSMILES("c1ccc(O)cc1");
    auto mapping = smsd::findMCS(query, target, smsd::ChemOptions{}, smsd::MCSOptions{});
    return mapping.size() == 6 ? 0 : 1;
}
```

```bash
c++ -std=c++17 -Icpp/include example.cpp -o example
./example
```

See the [C++ module](cpp/README.md) and [C++ API guide](docs/CPP.md) for
chemistry options, file formats and integration.

## Build from source

Use Java 25, a C++17 compiler and a supported Python environment:

```bash
git clone https://github.com/asad/SMSD.git
cd SMSD
mvn -f java/pom.xml clean verify
cmake -S cpp -B build/cpp -DCMAKE_BUILD_TYPE=Release \
  -DSMSD_BUILD_METAL=OFF -DSMSD_BUILD_CUDA=OFF
cmake --build build/cpp --config Release --parallel 4
python -m pip install -e .
```

Each module has its own instructions: [Java](java/README.md),
[C++](cpp/README.md) and [Python](python/README.md).
The [installation guide](docs/HOWTO-INSTALL.md) covers prerequisites.

## Docker CLI

The release includes Linux arm64 and AMD64 Docker image archives. Select the
architecture matching your Docker host:

```bash
SMSD_DOCKER_ARCH=arm64 # use amd64 for an x86_64 Docker host
gh release download v7.2.2 --repo asad/SMSD --dir downloads \
  --pattern "smsd-7.2.2-docker-linux-${SMSD_DOCKER_ARCH}.tar.gz"
docker load --input "downloads/smsd-7.2.2-docker-linux-${SMSD_DOCKER_ARCH}.tar.gz"
docker run --rm "smsd:7.2.2-linux-${SMSD_DOCKER_ARCH}" --version
```

See the release [Docker guide](https://github.com/asad/SMSD/releases/download/v7.2.2/DOCKER.md)
for checksums and search commands.

## Documentation

| Guide | Contents |
| --- | --- |
| [Java API](docs/JAVA.md) | CDK containers, search, fingerprints and CLI |
| [Python API](docs/PYTHON.md) | Native graphs, RDKit integration and bindings |
| [C++ API](docs/CPP.md) | Header-only library and chemistry options |
| [Examples](docs/EXAMPLES.md) | Feature examples and search contracts |
| [Installation](docs/HOWTO-INSTALL.md) | Platform requirements and source builds |
| [Native installers](docs/INSTALLERS.md) | MSI, DMG and DEB installation |
| [Release notes](docs/RELEASE_NOTES.md) | Changes in 7.2.2 |
| [Changelog](CHANGELOG.md) | Versioned API and behaviour changes |
| [Publishing](docs/PUBLISHING.md) | GitHub, PyPI and Maven Central publication |

The release also includes a [complete Java API reference](https://github.com/asad/SMSD/releases/download/v7.2.2/smsd-7.2.2-javadoc.jar).
Search budgets limit work; a returned MCS mapping can be a valid result without
proving global optimality. Atom mappings use zero-based input indices.

## Validation and benchmarks

Java 7.2.2 passed its full suite on Java 8 and 25; Python tests passed on macOS,
Linux and Windows. See [validation results](docs/VALIDATION_7.2.2.md) for exact
platforms, counts and installer checks.

[Benchmark results](benchmarks/RESULTS_7.2.0.md) retain their measured versions,
hardware and matching policies. They were not rerun for 7.2.2.
See [benchmark instructions](benchmarks/README.md) to reproduce them.

## Licence and attribution

SMSD is distributed under the [Apache License 2.0](LICENSE).
Retain the [NOTICE](NOTICE), licence and copyright notices when redistributing.
Bundled dependencies include their own licence notices.

If you use SMSD in research, please cite:

> Rahman SA. *SMSD Pro: Coverage-Driven, Tautomer-Aware Maximum Common Substructure Search.*
> ChemRxiv, 2026.
> DOI: [10.26434/chemrxiv.15001534/v1](https://doi.org/10.26434/chemrxiv.15001534/v1)

> Rahman SA, Bashton M, Holliday GL, Schrader R, Thornton JM.
> *Small Molecule Subgraph Detector (SMSD) toolkit.*
> Journal of Cheminformatics, 1:12, 2009.
> DOI: [10.1186/1758-2946-1-12](https://doi.org/10.1186/1758-2946-1-12)

Citation metadata is available in [CITATION.cff](CITATION.cff).

**Syed Asad Rahman — BioInception PVT LTD**

Copyright (c) 2018-2026 BioInception PVT LTD.
Algorithm Copyright (c) 2009-2026 Syed Asad Rahman.
