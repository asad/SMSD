# Java module

This module contains the Java implementation and CLI, using CDK 2.13.
Java 8 is the minimum target; Java 25 LTS is preferred for builds and execution.
The same JAR runs on both. Release validation covers both runtimes.
Its Maven coordinates remain `com.bioinceptionlabs:smsd`; this source targets
version `7.2.2`. The current published GitHub release is 7.2.1;
7.2.2 preparation and Maven Central publication are pending.

From the repository root:

```bash
mvn -f java/pom.xml clean verify
java -jar java/target/smsd-7.2.2-jar-with-dependencies.jar --version
```

The root `pom.xml` also runs this module with `mvn clean verify`. Include
the normally excluded algorithm and stress suites with
`-Dslow.tests.exclude=nothing`. Benchmark suites remain opt-in.

To run the same compiled tests on a Java 8 runtime while Maven uses Java 25:

```bash
mvn -f java/pom.xml -Dslow.tests.exclude=nothing \
  "-Djvm=$SMSD_JAVA8_HOME/bin/java" test
```

Set `SMSD_JAVA8_HOME` to an installed Java 8 JDK. Published 7.2.1 JARs require
Java 25; the Java 8 target starts with 7.2.2.

Artifacts and generated launchers are under `java/target/`. The source
launchers are under `java/src/scripts/`. The module packages the repository's
root LICENSE and NOTICE; keep those files when building from a checkout.
Use `java/pom.xml` for publishing the Java artifacts; the root POM only
aggregates local builds.

The 7.2.2 code passed 1,276 tests with 15 opt-in skips on both Java 8 and
Java 25, using the same compiled classes. All four JARs target Java 8 and
retain source and licence copies; both Unix launchers ran on both runtimes.
Maven publication remains pending. See [current validation](../docs/VALIDATION_7.2.2.md)
and [historical 7.2.1 results](../docs/VALIDATION_7.2.1.md).

See the [Java guide](../docs/JAVA.md) for CDK examples and the
[installation guide](../docs/HOWTO-INSTALL.md) for platform requirements.

Version 7.2.2 fixes CML/PDB input and rejects empty or multi-molecule/model files.
Use SDF for batch targets. See [current validation](../docs/VALIDATION_7.2.2.md)
and [native installers](../docs/INSTALLERS.md).
