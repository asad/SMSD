# SMSD Pro for Java

SMSD Pro for Java searches CDK 2.13 molecules for substructures and maximum
common substructures (MCS). It has two Maven modules, so an application takes
only what it uses:

| Module | Artefact | Contents |
| --- | --- | --- |
| `java/core` | `com.bioinceptionlabs:smsd` | The library: substructure, MCS and SMARTS search, fingerprints, stereo descriptors and standardisation |
| `java/cli` | `com.bioinceptionlabs:smsd-cli` | The command-line interface (CLI), its file readers and JSON output, and the self-contained JAR |

Both run on Java 8 and later; Java 25 LTS is recommended.
The [supported-platforms table](../docs/HOWTO-INSTALL.md#supported-platforms)
also covers the CLI installers.

## Installation

### Library

Maven Central has version 8.0.3, without the 8.1.0 Java changes; to use those,
[build from source](#build-from-source). Maven:

```xml
<dependency>
  <groupId>com.bioinceptionlabs</groupId>
  <artifactId>smsd</artifactId>
  <version>8.0.3</version>
</dependency>
```

Gradle (Kotlin DSL; in the Groovy DSL, write
`implementation 'com.bioinceptionlabs:smsd:8.0.3'`):

```kotlin
repositories {
    mavenCentral()
}

dependencies {
    implementation("com.bioinceptionlabs:smsd:8.0.3")
}
```

The library depends only on `cdk-interfaces`, `cdk-core`, `cdk-silent`,
`cdk-standard`, `cdk-isomorphism`, `cdk-smarts`, `cdk-valencycheck`,
`cdk-atomtype`, `vecmath` and `jackson-annotations`. Add what your application
needs besides:

| To | Add |
| --- | --- |
| Parse SMILES, as the quick start does | `org.openscience.cdk:cdk-smiles:2.13` |
| Read MOL and SDF files | `org.openscience.cdk:cdk-ctab:2.13` |
| Read CML or Mol2 files | `org.openscience.cdk:cdk-io:2.13` |
| Read PDB files | `org.openscience.cdk:cdk-pdb:2.13` |
| Call `Standardiser.PredicateRegistry.loadJsonFile`, or read and write `SearchResult` as JSON | `com.fasterxml.jackson.core:jackson-databind:2.21.7` |

Without `jackson-databind`, `loadJsonFile` throws `IllegalStateException`; the
rest of the library does not need it. On the module path, the library is the
automatic module `com.bioinception.smsd`.

### Command-line tool

The self-contained `smsd-cli-8.0.3-jar-with-dependencies.jar` runs with
`java -jar`, in the same way as the source-built 8.1.0 JAR under
[Command-line interface](#command-line-interface).
The standalone command-line tool stays at 8.0.3 in this release: download it
from the [8.0.3 release](https://github.com/asad/SMSD/releases/tag/v8.0.3), or
copy it from Maven Central:

```bash
mvn dependency:copy \
  -Dartifact=com.bioinceptionlabs:smsd-cli:8.0.3:jar:jar-with-dependencies \
  -DoutputDirectory=.
```

The CPU [Docker images](../docs/HOWTO-INSTALL.md#docker) of this release carry
the 8.1.0 tool, and a [source build](#build-from-source) writes it to
`java/cli/target`.

To call the CLI classes, such as `com.bioinception.smsd.cli.SMSDcli`, from an
application, depend on the CLI module, which brings the library, picocli,
Jackson and the CDK file readers. Maven:

```xml
<dependency>
  <groupId>com.bioinceptionlabs</groupId>
  <artifactId>smsd-cli</artifactId>
  <version>8.0.3</version>
</dependency>
```

Gradle (Kotlin DSL):

```kotlin
repositories {
    mavenCentral()
}

dependencies {
    implementation("com.bioinceptionlabs:smsd-cli:8.0.3")
}
```

On the module path, the CLI JAR is the automatic module
`com.bioinception.smsd.cli`. Packaged installers for Windows, macOS and Linux
are described under
[Java CLI installers](../docs/HOWTO-INSTALL.md#portable-java-cli).

### Build from source

From the repository root, build both modules and install them, with their
parent `com.bioinceptionlabs:smsd-parent` (`java/pom.xml`), into your local
Maven repository:

```bash
mvn -f java/pom.xml install
```

The build writes:

| File | Use |
| --- | --- |
| `java/core/target/smsd-8.1.0.jar` | Library; Maven resolves its dependencies |
| `java/cli/target/smsd-cli-8.1.0.jar` | CLI classes; Maven resolves their dependencies |
| `java/cli/target/smsd-cli-8.1.0-jar-with-dependencies.jar` | Standalone CLI, and a self-contained class path for the examples |
| `java/cli/target/appassembler/` | Portable `bin/smsd` launchers and the JARs they use, in `repo/` |

## SMARTS export

`MCS.matchWithSmarts(target)` returns the ordinary result and a pattern written
under the prepared chemistry. `SearchResult.Match.toSmarts(query, chemistry)`
and `SmartsWriter.toSmarts(graph, chemistry)` also support explicit export.
See [SMARTS export](../docs/SMARTS_EXPORT.md) for an MCS-to-library example and
the unsupported stereo, tautomer and ring constraints. Portable SMARTS does not
encode the induced MCIS objective or establish an optimality proof.

## Quick start

Save this as `CDKExample.java`. It checks that benzene is a substructure of
phenol. It then confirms that the MCS of N-methylacetamide and propanamide,
neither of which contains the other, is a proved four-atom optimum, and prints
its only maximum mapping, query atoms 0 to 3 onto target atoms 1 to 4, in the
original CDK atom indices. A Maven or Gradle project needs `cdk-smiles` beside
`smsd` to run it:

```java
import com.bioinception.smsd.core.MCS;
import com.bioinception.smsd.core.SearchResult;
import com.bioinception.smsd.core.Substructure;
import java.util.Map;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;

public class CDKExample {
    public static void main(String[] args) throws Exception {
        SmilesParser parser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer benzene = parser.parseSmiles("c1ccccc1");
        IAtomContainer phenol = parser.parseSmiles("c1ccc(O)cc1");
        if (!Substructure.find(benzene).match(phenol).found()) {
            throw new IllegalStateException("Expected benzene in phenol");
        }
        IAtomContainer query = parser.parseSmiles("CC(=O)NC");   // N-methylacetamide
        IAtomContainer target = parser.parseSmiles("CCC(=O)N");  // propanamide
        SearchResult result = MCS.find(query).withTimeout(1000).match(target);
        if (!result.optimal() || !result.found()
                || result.matches().get(0).atomCount() != 4
                || MCS.find(query).withTimeout(1000).matchAll(target).matches().size() != 1) {
            throw new IllegalStateException("Expected one proved four-atom MCS");
        }
        SearchResult.Match match = result.matches().get(0);
        for (Map.Entry<Integer, Integer> pair : match.atomMapping().entrySet()) {
            System.out.println(query.getAtom(pair.getKey()).getSymbol() + pair.getKey()
                    + " -> " + target.getAtom(pair.getValue()).getSymbol() + pair.getValue());
        }
        System.out.println("Common fragment: " + match.toSmiles(query));
    }
}
```

After the build above, compile and run it from the repository root against the
self-contained JAR, which includes the library, CDK and `cdk-smiles`:

```bash
javac --release 8 -cp java/cli/target/smsd-cli-8.1.0-jar-with-dependencies.jar CDKExample.java
java -cp "java/cli/target/smsd-cli-8.1.0-jar-with-dependencies.jar:." CDKExample
```

With a Java 8 compiler, omit `--release 8`. On Windows, separate the classpath
entries with a semicolon: `"java/cli/target/smsd-cli-8.1.0-jar-with-dependencies.jar;."`.

By default, `MCS` finds a connected maximum common induced subgraph (MCIS),
maximising atoms and then bonds, and `Substructure` finds a non-induced
embedding of the whole query; both provide `match`, `matchAll` and
`matchMany`. Check `status()`, `optimal()` and `complete()` before relying on
a result: a search stopped by its timeout can return a valid match without
proving it optimal. The [Java guide](../docs/JAVA.md) covers the maximum
common edge subgraph (MCES) objective, disconnected fragments, atom weights,
timeouts and the full result contract.

## Command-line interface

From the repository root:

```bash
java -jar java/cli/target/smsd-cli-8.1.0-jar-with-dependencies.jar --version
java -jar java/cli/target/smsd-cli-8.1.0-jar-with-dependencies.jar \
  --Q SMI --q "c1ccccc1" --T SMI --t "c1ccc(O)cc1" --json -
java -jar java/cli/target/smsd-cli-8.1.0-jar-with-dependencies.jar \
  --Q SMI --q "CC(=O)NC" --T SMI --t "CCC(=O)N" --mode mcs --json -
```

The first search tests whether benzene is a substructure of phenol; the second
finds the four-atom MCS of N-methylacetamide and propanamide
(`mcs_size` is 4). `--help` lists every option, and the
[Java guide](../docs/JAVA.md#command-line-interface) describes input formats,
JSON output and exit status.

## Build and test

From the repository root, with JDK 25:

```bash
mvn -f java/pom.xml clean verify
```

This builds and tests the library, then the CLI; `mvn clean verify` from the
repository root runs the same build. Add `-pl core` for the library alone, or
`-pl cli -am` for the CLI and the library it needs. `AlgorithmTest` and
`StressTest` run by default; only tests tagged `SlowTest` are excluded. Add
`-Dslow.tests.exclude=nothing` to include the tests tagged `SlowTest`; with it,
`-Dbenchmark=true` and `-Dsmsd.benchmark=true` also run the benchmark tests.
Tests named with `-Dtest`, or chosen by tag with `-Dgroups`, run whatever their
tags. To leave out other tags as well, list them in `-Dslow.tests.exclude`, for
example `-Dslow.tests.exclude=SlowTest,MyTag`; Surefire's own `-DexcludedGroups`
is not used.

To run the tests on a Java 8 runtime while Maven runs on Java 25, set
`SMSD_JAVA8_HOME` to a Java 8 JDK and run:

```bash
mvn -f java/pom.xml -Dslow.tests.exclude=nothing \
  "-Djvm=$SMSD_JAVA8_HOME/bin/java" test
```

The library sources and tests are in `java/src/main/java` and
`java/src/test/java`, built by `java/core/pom.xml`; the CLI sources, tests and
launcher scripts are in `java/cli/src/main/java`, `java/cli/src/test/java` and
`java/cli/src/scripts`. The library and CLI JARs carry the repository's
`LICENSE` and `NOTICE` files under `META-INF/smsd/`, and the CDK licence under
`META-INF/smsd/licenses/cdk/`.

## Documentation

- [Java guide](../docs/JAVA.md): API reference, options, results, fingerprints,
  stereo, R-groups and the CLI.
- [Installation guide](../docs/HOWTO-INSTALL.md): platform requirements and
  source builds.
- [Benchmark guide](../benchmarks/README.md) and [changelog](../CHANGELOG.md).

## Licence

Distributed under the [Apache License 2.0](../LICENSE). See [NOTICE](../NOTICE)
for attribution and bundled dependency notices.
