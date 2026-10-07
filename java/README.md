# SMSD Java

The Java library and CLI use CDK 2.13. Version 7.2.2 supports Java 8 or later;
Java 25 LTS is preferred. The same JAR runs on both.

## Installation

Download the [7.2.2 release](https://github.com/asad/SMSD/releases/tag/v7.2.2):

| Package | Use |
| --- | --- |
| `smsd-7.2.2-jar-with-dependencies.jar` | Standalone CLI or complete example classpath |
| `smsd-7.2.2.jar` | Library with dependencies resolved by Maven |
| `smsd-7.2.2-cli.tar.gz` | Portable CLI and shell launchers |
| `smsd-7.2.2-sources.jar` | Java source |
| `smsd-7.2.2-javadoc.jar` | Complete API reference |

Windows AMD64 MSI, macOS arm64 DMG and Linux AMD64 DEB packages install the
terminal CLI with Java 25 LTS included. See [native installers](../docs/INSTALLERS.md)
for installation and platform requirements.

Use these Maven coordinates for 7.2.2:

```xml
<dependency>
  <groupId>com.bioinceptionlabs</groupId>
  <artifactId>smsd</artifactId>
  <version>7.2.2</version>
</dependency>
```

Until 7.2.2 is published to Maven Central, install it from a source checkout
using Java 25:

```bash
mvn -f java/pom.xml install
```

The Maven dependency includes the required CDK modules, including SMARTS support.

## CDK example

Save this as `CDKExample.java`. The mapping refers to atom indices in the
original CDK containers:

```java
import com.bioinception.smsd.core.ChemOptions;
import com.bioinception.smsd.core.MolGraph;
import com.bioinception.smsd.core.SearchEngine;
import java.util.Map;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;

public class CDKExample {
    public static void main(String[] args) throws Exception {
        SmilesParser parser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer query = parser.parseSmiles("c1ccccc1");
        IAtomContainer target = parser.parseSmiles("c1ccc(O)cc1");
        ChemOptions chemistry = new ChemOptions();
        MolGraph queryGraph = new MolGraph(query);
        MolGraph targetGraph = new MolGraph(target);
        SearchEngine.MCSOptions options = new SearchEngine.MCSOptions();
        options.timeoutMs = 1000;
        options.connectedOnly = true;
        Map<Integer, Integer> mapping =
                SearchEngine.findMCS(queryGraph, targetGraph, chemistry, options);
        if (mapping.size() != 6 || !SearchEngine.validateMapping(
                queryGraph, targetGraph, mapping, chemistry).isEmpty()) {
            throw new IllegalStateException("Expected a valid six-atom benzene match");
        }
        for (Map.Entry<Integer, Integer> pair : mapping.entrySet()) {
            System.out.println(query.getAtom(pair.getKey()).getSymbol() + " -> "
                    + target.getAtom(pair.getValue()).getSymbol());
        }
    }
}
```

With the downloaded dependency JAR in the current directory, compile using
Java 25 and run using Java 8 or later:

```bash
javac --release 8 -cp smsd-7.2.2-jar-with-dependencies.jar CDKExample.java
java -cp "smsd-7.2.2-jar-with-dependencies.jar:." CDKExample
```

When compiling with Java 8, omit `--release 8`. On Windows, use
`"smsd-7.2.2-jar-with-dependencies.jar;."` as the classpath.

## CLI

```bash
java -jar smsd-7.2.2-jar-with-dependencies.jar --version
java -jar smsd-7.2.2-jar-with-dependencies.jar --help
java -jar smsd-7.2.2-jar-with-dependencies.jar \
  --Q SMI --q "c1ccccc1" --T SMI --t "c1ccc(O)cc1" --json -
java -jar smsd-7.2.2-jar-with-dependencies.jar \
  --Q SMI --q "c1ccccc1" --T SMI --t "c1ccc(O)cc1" --mode mcs --json -
```

File input supports MOL, MOL2 (`ML2`), CML and PDB. A single-molecule/model
input is required; use SDF for batch targets. JSON output uses UTF-8.

## Build and test

From the repository root, using Java 25:

```bash
mvn -f java/pom.xml clean verify
java -jar java/target/smsd-7.2.2-jar-with-dependencies.jar --version
```

The root `pom.xml` also runs this module with `mvn clean verify`. Include the
algorithm and stress suites with `-Dslow.tests.exclude=nothing`; benchmarks
remain opt-in. To test the same compiled classes on Java 8 while Maven uses
Java 25, set `SMSD_JAVA8_HOME` to an installed Java 8 JDK, then run:

```bash
mvn -f java/pom.xml -Dslow.tests.exclude=nothing \
  "-Djvm=$SMSD_JAVA8_HOME/bin/java" surefire:test
```

Sources and shell launchers are under `java/src/`; build output is under
`java/target/`. Keep the root `LICENSE` and `NOTICE` when building from a checkout.
Use `java/pom.xml` to publish the library; the root POM aggregates local builds.

## Reference

- [Java API guide](../docs/JAVA.md): search options, fingerprints, stereo and R-groups.
- [Examples](../docs/EXAMPLES.md): feature examples and matching contracts.
- [Complete API reference](https://github.com/asad/SMSD/releases/download/v7.2.2/smsd-7.2.2-javadoc.jar).
- [Validation](../docs/VALIDATION_7.2.2.md): Java 8 and 25 test results.
- [Changelog](../CHANGELOG.md): API and behaviour changes.

Distributed under the [Apache License 2.0](../LICENSE); see [NOTICE](../NOTICE)
for attribution and bundled dependency notices.
