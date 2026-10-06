# SMSD Pro Java Guide

The Java API uses CDK 2.13 for molecular input and standardisation, with SMSD
algorithms for substructure, MCS, fingerprints, stereo/CIP, layout and R-group
decomposition. This checkout targets version 7.2.2, currently in preparation.
GitHub's published release is 7.2.1; Maven Central remains at 7.1.1.
Java sources and launchers are in `java/src/`; Maven output is in
`java/target/`. The root aggregator supports `mvn verify`, while direct module
builds use `mvn -f java/pom.xml`. Version 7.2.2 targets Java 8 or later;
Java 25 LTS is preferred. The same JAR is tested on Java 8 and Java 25.
They accept CDK `IAtomContainer` inputs;
the native C++/Python graph layer uses its own molecule representation.

Current measurements and their matching policies are in the
[benchmark report](../benchmarks/RESULTS_7.2.0.md). They do not establish a
universal ranking between Java, native SMSD or CDK.

## Install

The current published Maven Central version is 7.1.1. It does not include the
7.2.2 changes documented here:

```xml
<dependency>
  <groupId>com.bioinceptionlabs</groupId>
  <artifactId>smsd</artifactId>
  <version>7.1.1</version>
</dependency>
```

To use these changes, install this checkout into your local Maven
repository and set the dependency version to `7.2.2`:

```sh
mvn -f java/pom.xml install
```

The 7.2.2 coordinate is available from that local build until it is published
to Maven Central. GitHub's 7.2.1 release is separate from the registry version.

Run the locally built CLI:

```bash
java -jar java/target/smsd-7.2.2-jar-with-dependencies.jar --Q SMI --q "c1ccccc1" --T SMI --t "c1ccc(O)cc1" --json -
```

## Core API

This example follows CDK's parser/container workflow. It returns indices into
the two containers passed to `MolGraph`, so each pair can be used directly with
CDK's `getAtom(index)`:

```java
import java.util.Map;
import com.bioinception.smsd.core.ChemOptions;
import com.bioinception.smsd.core.MolGraph;
import com.bioinception.smsd.core.SearchEngine;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;

public class MCSExample {
  public static void main(String[] args) throws Exception {
    SmilesParser parser = new SmilesParser(SilentChemObjectBuilder.getInstance());
    IAtomContainer query = parser.parseSmiles("NCCO");
    IAtomContainer target = parser.parseSmiles("CC(O)CN");
    MolGraph g1 = new MolGraph(query);
    MolGraph g2 = new MolGraph(target);

    ChemOptions chemistry = new ChemOptions();
    SearchEngine.MCSOptions options = new SearchEngine.MCSOptions();
    options.timeoutMs = 1000;
    options.connectedOnly = true;
    options.induced = false;
    Map<Integer, Integer> mapping = SearchEngine.findMCS(g1, g2, chemistry, options);
    for (Map.Entry<Integer, Integer> pair : mapping.entrySet()) {
      System.out.printf("query[%d] %s -> target[%d] %s%n",
          pair.getKey(), query.getAtom(pair.getKey()).getSymbol(),
          pair.getValue(), target.getAtom(pair.getValue()).getSymbol());
    }
    if (!SearchEngine.validateMapping(g1, g2, mapping, chemistry).isEmpty())
      throw new IllegalStateException("Invalid mapping");
  }
}
```

Save it as `MCSExample.java`, then compile and run against the shaded JAR:

```bash
javac -cp java/target/smsd-7.2.2-jar-with-dependencies.jar MCSExample.java
java -cp java/target/smsd-7.2.2-jar-with-dependencies.jar:. MCSExample
```

On Windows, use a semicolon classpath separator and quote the classpath:

```powershell
java -cp "java/target/smsd-7.2.2-jar-with-dependencies.jar;." MCSExample
```

`MolGraph` retains CDK atom order. Weights and target exclusions use those
zero-based indices, rather than atom-map numbers or canonical-SMILES positions.
If preprocessing changes a container's atoms, retain its correspondence to the
original molecule before searching. Prepare aromaticity, hydrogens and stereo
consistently when importing containers from other CDK readers.

For substructure matching, use the same containers and chemistry options:

```java
Map<Integer, Integer> first = SearchEngine.findSubstructure(query, target, chemistry, 1000);
java.util.List<Map<Integer, Integer>> all = SearchEngine.findAllSubstructures(query, target, chemistry, 100, 1000);
boolean hit = SearchEngine.isSubstructure(query, target, chemistry, 1000);
```

The `SMSD` facade also accepts CDK containers and standardises them by default:

```java
com.bioinception.smsd.core.SMSD matcher = new com.bioinception.smsd.core.SMSD(query, target, chemistry);
Map<Integer, Integer> mapping = matcher.findMCS(false, true, 1000); // non-induced, connected
```

Pass `standardise=false` to the four-argument constructor only for inputs you
have already prepared. Reuse graphs for repeated matching. CDK `DfPattern` and
SMSD have different aromaticity/chemistry policies; equal input containers do
not by themselves establish equivalent matching rules.

## Chemistry and search options

Chemistry belongs to `ChemOptions`; the objective, connectivity and budget
belong to `SearchEngine.MCSOptions`:

```java
ChemOptions chemistry = new ChemOptions();
chemistry.matchFormalCharge = true;
chemistry.matchIsotope = true;
chemistry.useChirality = true;
chemistry.matchBondOrder = ChemOptions.BondOrderMode.STRICT;
chemistry.aromaticityMode = ChemOptions.AromaticityMode.STRICT;

SearchEngine.MCSOptions options = new SearchEngine.MCSOptions();
options.timeoutMs = 1000;           // milliseconds across all search stages
options.connectedOnly = true;
options.maximizeBonds = false;      // atom objective
options.atomWeights = new double[] {10, -30, 1}; // for a three-atom query only
```

Use `ChemOptions.tautomerProfile()` explicitly for relaxed tautomer chemistry.
Specifying a bond objective or weights changes the selected mapping; its atom
count alone cannot be compared with an atom-maximizing result.

## Fingerprints

```java
long[] pathFp = SearchEngine.pathFingerprint(mol1, 7, 2048);
long[] mcsFp = SearchEngine.mcsFingerprint(mol1, new ChemOptions(), 7, 2048);
boolean subset = SearchEngine.fingerprintSubset(pathFp, pathFp);
double sim = SearchEngine.mcsFingerprintSimilarity(mcsFp, mcsFp);
```

## Stereo, Layout, and R-groups

```java
Map<Integer, Character> rs = com.bioinception.smsd.core.CIPAssigner.assignRS(g1);
Map<Long, Character> ez = com.bioinception.smsd.core.CIPAssigner.assignEZ(g1);
java.util.List<Map<String, IAtomContainer>> rgroups = SearchEngine.decomposeRGroups(core, molecules, new ChemOptions(), 10_000);
```

## Matching contracts

Non-induced MCS retains the caller's query direction: mapped query edges must
exist in the target, while extra target edges are allowed. Atom weights apply
to query indices and must be finite. Java compares their double-precision
sums; the native C++ scoring API uses integer millipoints, so sub-millipoint
weights can have different tie behavior between these APIs.

Tautomer-aware matching relaxes applicable bonds while preserving element
identity when atom-type matching is enabled. Chirality compares R/S descriptors
and checks fully mapped tetrahedral ligand permutations. Unspecified stereo
retains its existing wildcard behavior.

`excludedTargetAtoms` excludes original target indices without changing the
target's chemistry. Constrained batches retain those indices, prevent reuse,
rank targets by the requested objective, and use the per-pair timeout argument.
Timeouts cover orientation, recovery and retry stages together; large searches
can still return a valid incumbent without proving global optimality.

`canonicalizeMapping` returns the exact minimum of the explored automorphism
orbit only after closure completes. It throws `IllegalStateException` when
generators are incomplete or the orbit exceeds its time/storage bounds, rather
than presenting a partial result as canonical. Internal enumeration retains
raw mapping keys when symmetry work cannot complete. Captured generators
preserve full chemical properties; graph orbits use only those proven
permutations.
