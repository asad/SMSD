# Chemistry Development Kit (CDK)

The SMSD Pro™ Java library and command-line tool use the Chemistry Development
Kit (CDK) 2.13, which is licensed under the GNU Lesser General Public License,
version 2.1 or (at your option) any later version. The unmodified licence text
is in [LICENSE](LICENSE). CDK is Copyright (c) The CDK Project.

## Where CDK is bundled

In the file names below, `<version>` is the SMSD Pro version.

- `smsd-cli-<version>-jar-with-dependencies.jar` contains the classes of the
  CDK 2.13 modules that the SMSD Pro library and CLI use, unmodified.
- The Java CLI installers are built from that JAR and carry the same CDK
  classes, as do the CPU Docker images.
- The portable CLI archive, `smsd-<version>-cli.tar.gz`, carries the unmodified
  CDK 2.13 JARs beside `smsd-<version>.jar` and `smsd-cli-<version>.jar`.
- `smsd-<version>.jar` and `smsd-cli-<version>.jar` on Maven Central contain no
  CDK classes; Maven resolves the CDK modules as ordinary dependencies.

## Replacing CDK

SMSD Pro uses CDK only through its public API. To run SMSD Pro with another
CDK build, put `smsd-<version>.jar` (and `smsd-cli-<version>.jar` for the CLI),
that CDK build and the other dependencies listed in `java/core/pom.xml` (and
`java/cli/pom.xml`) on the class path in place of the JAR with dependencies, or
replace the `org/openscience/cdk` classes inside that JAR.

## Source

The CDK source code, including release 2.13, is published at
<https://github.com/cdk/cdk>, and the CDK modules are published on Maven
Central under the group `org.openscience.cdk`.

SMSD Pro's own code remains licensed under the Apache License, Version 2.0.
