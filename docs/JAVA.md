# SMSD Java guide

SMSD 7.2.2 integrates with CDK 2.13 and accepts CDK `IAtomContainer` inputs.
It supports substructure and MCS search, fingerprints, stereo handling,
standardisation and R-group decomposition. Java 8 is the minimum runtime;
Java 25 LTS is preferred. Both run the same JAR.

The prepared-query API described here is an unreleased breaking change. Build
from this source checkout; the released 7.2.2 JAR retains the previous API.

## Installation

The [GitHub release](https://github.com/asad/SMSD/releases/tag/v7.2.2) includes
the library, a JAR with dependencies, sources, Javadoc and portable CLI launchers.
Windows AMD64 MSI, macOS arm64 DMG and Linux AMD64 DEB installers include
Java 25 LTS. See [installation](HOWTO-INSTALL.md) and [native installers](INSTALLERS.md).

Download the complete example classpath:

```bash
curl -LO https://github.com/asad/SMSD/releases/download/v7.2.2/smsd-7.2.2-jar-with-dependencies.jar
```

For Maven projects, use:

```xml
<dependency>
  <groupId>com.bioinceptionlabs</groupId>
  <artifactId>smsd</artifactId>
  <version>7.2.2</version>
</dependency>
```

Maven Central publication is separate. Until 7.2.2 is available there, use
Java 25 to install the library from a source checkout:

```bash
mvn -f java/pom.xml install
```

Required CDK dependencies, including SMARTS support, are resolved transitively.
The [Java module](../java/README.md) describes source builds and tests.

## MCS and substructure search

Save this as `SearchExample.java`. It builds CDK containers, searches their
chemical graphs and checks that the returned mapping is valid:

```java
import com.bioinception.smsd.core.MCS;
import com.bioinception.smsd.core.SearchResult;
import com.bioinception.smsd.core.Substructure;
import java.util.Map;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;

public class SearchExample {
    public static void main(String[] args) throws Exception {
        SmilesParser parser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer query = parser.parseSmiles("c1ccccc1");
        IAtomContainer target = parser.parseSmiles("c1ccc(O)cc1");
        SearchResult result = MCS.find(query).withTimeout(1000).match(target);
        if (!result.optimal() || !result.found()
                || result.matches().get(0).atomCount() != 6
                || !Substructure.find(query).match(target).found()) {
            throw new IllegalStateException("Expected a proven six-atom benzene match");
        }
        SearchResult.Match match = result.matches().get(0);
        for (Map.Entry<Integer, Integer> pair : match.atomMapping().entrySet()) {
            System.out.println(query.getAtom(pair.getKey()).getSymbol() + " -> "
                    + target.getAtom(pair.getValue()).getSymbol());
        }
        System.out.println("Common fragment: " + match.toSmiles(query));
    }
}
```

Build the current source API, then compile with Java 25 and run on Java 8 or later:

```bash
mvn -f java/pom.xml package -DskipTests
javac --release 8 -cp java/target/smsd-7.2.2-jar-with-dependencies.jar SearchExample.java
java -cp "java/target/smsd-7.2.2-jar-with-dependencies.jar:." SearchExample
```

When compiling with Java 8, omit `--release 8`. On Windows, use a semicolon
classpath separator:

```powershell
java -cp "java/target/smsd-7.2.2-jar-with-dependencies.jar;." SearchExample
```

`MolGraph` retains the supplied CDK atom order. Mappings, weights and exclusions
use zero-based atom indices, rather than atom-map numbers or canonical-SMILES
positions. Prepare aromaticity, hydrogens and stereo consistently when loading
containers through other CDK readers. If preprocessing changes atom order,
retain the correspondence to the original molecule.

`MCS.find(query)` and `Substructure.find(query)` accept either `MolGraph` or
CDK containers. `Substructure.smarts(expression)` compiles SMARTS, including
recursive and stereo predicates. Prepare chemistry explicitly with
`Standardiser` when the input reader has not supplied it. Prepared queries
own snapshots of their query and options and can be reused concurrently.
Keep target graphs unchanged while a call uses them.

## Chemistry, objectives and completion

`MCS.Options` contains `chemistry`, `objective`, `connectivity`, `atomWeights`,
`minFragmentAtoms`, `maxFragments`, `excludedTargetAtoms`, `timeoutMillis`,
`nodeLimit` and `cancellation`. `Substructure.Options` contains `chemistry`,
`induced` and the same resource limits.

| Choice | Meaning |
| --- | --- |
| `MCS.Objective.MCIS` | Preserve mapped edges and non-edges; maximise atoms, then bonds |
| `MCS.Objective.MCES` | Select compatible common edges; maximise bonds, then incident atoms; exclude isolated padding |
| `MCS.Connectivity.Connected` | Require one common fragment; the default |
| `MCS.Connectivity.Fragments` | Permit disconnected common fragments |
| `atomWeights` | MCIS only; maximise the query-weight sum, then atoms and bonds |
| `minFragmentAtoms`, `maxFragments` | Constrain selected common-edge components during optimisation |
| `excludedTargetAtoms` | Exclude original target indices while retaining the target's chemistry |
| `timeoutMillis` | `-1` chooses an adaptive budget; `0` stops immediately |
| `nodeLimit` | `0` means unlimited nodes |
| `cancellation` | A `BooleanSupplier` checked during searching; interruption also cancels |

Java retains supplied finite double-weight precision. Candidate sums use query
index order and comparisons use `Double.compare`; weights are not quantised.
C++ and Python retain their bounded millipoint weight contract.

`ChemOptions` controls element, charge, isotope, aromaticity, bond order, ring
and stereo comparisons. `ChemOptions.tautomerProfile()` explicitly relaxes
tautomer chemistry. The explicit
`MCS.computeTautomerConfidence(queryGraph, targetGraph, match, chemistry)` helper
reports a heuristic confidence for selected tautomeric chemistry relaxations.
`chemistry.pH` and `chemistry.solvent` apply only to this read-time helper; they do
not change MCIS or MCES optimisation. Identical selected chemistry has confidence
1, and the helper leaves both graphs unchanged.

With `completeRingsOnly`, each retained ring atom or bond
must be covered by a wholly selected relevant ring in both inputs. A shared
fusion atom may survive through one retained complete ring.

Both `match(target)` and `matchAll(target, limit)` return `SearchResult`.
`matchAll` defaults to 1000 distinct atom-and-selected-bond embeddings;
`MCS.canonicalizeMapping` can identify symmetry-equivalent correspondences.
`Match.atomMapping()` and `bondMapping()` use original input atom indices.
`fragments()` follows selected common edges. `toMolecule(query)` and
`toSmiles(query)` cap cut valences; `toGraph(query)` preserves original
chemical properties for further matching.

Inspect `status()` (`Completed`, `TimedOut`, `Cancelled`, `LimitReached`),
`optimal()` and `complete()` separately. A budget-stopped incumbent may be
valid without proving an optimum. A capped MCS enumeration can have
`optimal() == true` and `complete() == false`. Substructure `match` proves
existence when found; `matchAll` reports whether its requested enumeration
finished. Elapsed milliseconds and visited nodes are available on the result.

`MCS.matchMany` and `Substructure.matchMany` search a list of pre-parsed targets
sequentially and return one result per target in input order, including stopped
and negative results. Each target receives the configured per-call budget. `MCS.matchConstrained` searches
successive queries with separate per-pair budgets and excludes claimed atoms
in the original target graph. `MCS.consensus` uses one whole-call budget and
requires `ceil(threshold * moleculeCount)` support. Its progressive result
reports supported mappings and always has `optimal() == false`.

Stereo matching compares retained ligand orientations. Partial common graphs
constrain stereo once selected ligands determine it; unspecified stereo can
act as a wildcard.

## Prepared queries and fragments

Save this as `PreparedSearchExample.java`:

```java
import com.bioinception.smsd.core.MCS;
import com.bioinception.smsd.core.MolGraph;
import com.bioinception.smsd.core.SearchResult;
import com.bioinception.smsd.core.Substructure;
import java.util.Arrays;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;

public class PreparedSearchExample {
    public static void main(String[] args) throws Exception {
        SmilesParser parser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        MolGraph ring = new MolGraph(parser.parseSmiles("C1CC1"));
        MolGraph chain = new MolGraph(parser.parseSmiles("CCC"));
        MCS.Options options = new MCS.Options();
        options.objective = MCS.Objective.MCES;
        SearchResult edges = MCS.find(ring).withOptions(options).match(chain);
        if (!edges.optimal() || edges.matches().get(0).bondCount() != 2) {
            throw new IllegalStateException("Expected a proven two-bond common path");
        }

        options.connectivity = MCS.Connectivity.Fragments;
        options.minFragmentAtoms = 1;
        MolGraph pieces = new MolGraph(parser.parseSmiles("CO.CN"));
        SearchResult fragments = MCS.find(pieces).withOptions(options).match(pieces);
        if (fragments.matches().get(0).fragments().size() != 2) {
            throw new IllegalStateException("Expected two common fragments");
        }

        Substructure query = Substructure.find(new MolGraph(parser.parseSmiles("CO")));
        MolGraph target = new MolGraph(parser.parseSmiles("OCO"));
        if (query.matchAll(target, 10).matches().size() != 2
                || query.matchMany(Arrays.asList(target, chain)).get(1).found()) {
            throw new IllegalStateException("Unexpected enumeration or batch result");
        }
    }
}
```

Compile and run it with the same source-built JAR and classpath as
`SearchExample.java` above.

## Fingerprints, stereo and R-groups

Save this as `FeaturesExample.java`. Fingerprint similarities below compare
a molecule with itself; the R-group example decomposes phenol and toluene
against a benzene core:

```java
import com.bioinception.smsd.core.ChemOptions;
import com.bioinception.smsd.core.CIPAssigner;
import com.bioinception.smsd.core.FingerprintEngine;
import com.bioinception.smsd.core.MolGraph;
import com.bioinception.smsd.core.Substructure;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;

public class FeaturesExample {
    public static void main(String[] args) throws Exception {
        SmilesParser parser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer phenol = parser.parseSmiles("c1ccc(O)cc1");
        long[] path = FingerprintEngine.pathFingerprint(phenol, 7, 2048);
        long[] common = FingerprintEngine.mcsFingerprint(phenol, new ChemOptions(), 7, 2048);
        long[] circular = FingerprintEngine.ecfp(phenol, 2, 2048);
        if (!FingerprintEngine.fingerprintSubset(path, path)
                || FingerprintEngine.tanimoto(common, common) != 1.0
                || FingerprintEngine.tanimoto(circular, circular) != 1.0) {
            throw new IllegalStateException("Expected identical fingerprints to match");
        }

        MolGraph chiral = new MolGraph(parser.parseSmiles("N[C@@H](C)C(=O)O"));
        MolGraph alkene = new MolGraph(parser.parseSmiles("F/C=C/F"));
        Map<Integer, Character> rs = CIPAssigner.assignRS(chiral);
        Map<Long, Character> ez = CIPAssigner.assignEZ(alkene);
        if (rs.isEmpty() || ez.isEmpty()) {
            throw new IllegalStateException("Expected explicit stereo descriptors");
        }
        System.out.println("Atom descriptors: " + rs);
        System.out.println("Bond descriptors: " + ez);

        IAtomContainer core = parser.parseSmiles("c1ccccc1");
        List<IAtomContainer> molecules = Arrays.asList(
                phenol, parser.parseSmiles("Cc1ccccc1"));
        List<Substructure.Decomposition> groups =
                Substructure.decomposeRGroups(core, molecules, new ChemOptions(), 1000);
        if (groups.size() != 2) {
            throw new IllegalStateException("Expected two decomposition rows");
        }
        for (Substructure.Decomposition row : groups) {
            if (!row.search().complete() || row.core() == null) {
                throw new IllegalStateException("Core search did not complete with a match");
            }
            System.out.println("Attached groups: " + row.rgroups().size());
        }
    }
}
```

```bash
javac --release 8 -cp java/target/smsd-7.2.2-jar-with-dependencies.jar FeaturesExample.java
java -cp "java/target/smsd-7.2.2-jar-with-dependencies.jar:." FeaturesExample
```

A fingerprint match is a screening result, not a substitute for graph matching.
Fingerprint similarity and MCS overlap are different measures. CIP assignment
has bounded traversal; inspect the returned descriptors for the molecule being
used rather than assuming every stereo case is assigned.

## CLI

Check the version and available options:

```bash
java -jar smsd-7.2.2-jar-with-dependencies.jar --version
java -jar smsd-7.2.2-jar-with-dependencies.jar --help
```

Search for a benzene substructure, then its MCS with phenol:

```bash
java -jar smsd-7.2.2-jar-with-dependencies.jar \
  --Q SMI --q "c1ccccc1" --T SMI --t "c1ccc(O)cc1" --json -
java -jar smsd-7.2.2-jar-with-dependencies.jar \
  --Q SMI --q "c1ccccc1" --T SMI --t "c1ccc(O)cc1" --mode mcs --json -
```

Use query type `SIG` for a SMARTS pattern requiring a carbonyl group:

```bash
java -jar smsd-7.2.2-jar-with-dependencies.jar \
  --Q SIG --q "C=O" --T SMI --t "CC(=O)O" --json -
java -jar smsd-7.2.2-jar-with-dependencies.jar \
  --Q SIG --q "C=O" --T SMI --t "CCO" --json -
```

The first SMARTS search returns `"exists": true` and exit status 0. The second
returns `"exists": false` and exit status 1; a missing match is a normal search
result. Current source also uses 1 for completed MCS and enumeration requests
without a match, and 2 for invalid input, errors or incomplete searches.
In an SDF batch, 0 means at least one match and 1 means no matches; any failed
or incomplete row makes the exit status 2. Unknown modes are rejected.
File types are `MOL`,
`ML2` (MOL2), `CML` and `PDB`. Each single-file input must contain exactly one
molecule or model. Empty and ambiguous inputs are rejected. Use `SDF` for batch
targets. File paths and JSON output support UTF-8.

## API reference and validation

- [Complete Javadoc](https://github.com/asad/SMSD/releases/download/v7.2.2/smsd-7.2.2-javadoc.jar).
- [Feature examples](EXAMPLES.md).
- [Validation results](VALIDATION_7.2.2.md).
- [Benchmark results](../benchmarks/RESULTS_7.2.0.md), with measured versions and matching policies.
- [Changelog](../CHANGELOG.md).

`canonicalizeMapping` requires a complete automorphism orbit. It throws
`IllegalStateException` when the generators are incomplete or its time/storage
bounds are exceeded. Constrained batch APIs preserve target indices and prevent
reuse; consult the Javadoc for their result types and per-pair budgets.

The library is distributed under the [Apache License 2.0](../LICENSE).
Retain [NOTICE](../NOTICE) and the bundled dependency notices when redistributing.
