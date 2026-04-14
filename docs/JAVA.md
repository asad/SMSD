# SMSD Pro Java Guide

SMSD Pro’s Java API is the full-featured chemistry layer built on CDK, with native
SMSD Pro algorithms for substructure, MCS, fingerprints, stereo/CIP, layout, and
R-group decomposition.

## Install

Maven:

```xml
<dependency>
  <groupId>com.bioinceptionlabs</groupId>
  <artifactId>smsd</artifactId>
  <version>7.1.0</version>
</dependency>
```

CLI:

```bash
java -jar smsd-7.1.0-jar-with-dependencies.jar --Q SMI --q "c1ccccc1" --T SMI --t "c1ccc(O)cc1" --json -
```

## Core API

```java
import com.bioinception.smsd.core.ChemOptions;
import com.bioinception.smsd.core.MolGraph;
import com.bioinception.smsd.core.SearchEngine;

MolGraph g1 = new MolGraph(mol1);
MolGraph g2 = new MolGraph(mol2);

boolean hit = SearchEngine.isSubstructure(g1, g2, new ChemOptions(), 10_000);
var mcs = SearchEngine.findMCS(g1, g2, new ChemOptions(), new SearchEngine.MCSOptions());
```

## Fingerprints

```java
long[] pathFp = SearchEngine.pathFingerprint(mol1, 7, 2048);
long[] mcsFp = SearchEngine.mcsFingerprint(mol1, new ChemOptions(), 7, 2048);
boolean subset = SearchEngine.fingerprintSubset(pathFp, pathFp);
double sim = SearchEngine.mcsFingerprintSimilarity(mcsFp, mcsFp);
```

## Stereo, Layout, and R-groups

```java
var rs = com.bioinception.smsd.core.CIPAssigner.assignRS(g1);
var ez = com.bioinception.smsd.core.CIPAssigner.assignEZ(g1);
var rgroups = SearchEngine.decomposeRGroups(core, molecules, new ChemOptions(), 10_000);
```

Java remains the strongest high-level chemistry integration layer in SMSD Pro for
CDK-based applications, while staying aligned with the native C++ behavior on
core graph-matching semantics.
