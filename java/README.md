# Java module

This module contains the Java 25 implementation and CLI, using CDK 2.13.
Its Maven coordinates remain `com.bioinceptionlabs:smsd`; this source targets
the unreleased version `7.2.1`.

From the repository root:

```bash
mvn -f java/pom.xml clean verify
java -jar java/target/smsd-7.2.1-jar-with-dependencies.jar --version
```

The root `pom.xml` also runs this module with `mvn clean verify`. Include
the normally excluded algorithm and stress suites with
`-Dslow.tests.exclude=nothing`. Benchmark suites remain opt-in.

Artifacts and generated launchers are under `java/target/`. The source
launchers are under `java/src/scripts/`. The module packages the repository's
root LICENSE and NOTICE; keep those files when building from a checkout.
Use `java/pom.xml` for publishing the Java artifacts; the root POM only
aggregates local builds.

The 7.2.1 module has 1,242 distinct passing cases and 15 opt-in skips across
the full and focused runs. All four JARs and both Unix launchers were checked;
publication remains pending. See [release validation](../docs/VALIDATION_7.2.1.md)
for the combined-run scope and artifact checks.

See the [Java guide](../docs/JAVA.md) for CDK examples and the
[installation guide](../docs/HOWTO-INSTALL.md) for platform requirements.
