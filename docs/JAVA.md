# SMSD Java guide

SMSD 7.2.2 integrates with CDK 2.13 and accepts CDK `IAtomContainer` inputs.
It supports substructure and MCS search, fingerprints, stereo handling,
standardisation and R-group decomposition. Java 8 is the minimum runtime;
Java 25 LTS is preferred. Both run the same JAR.

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
import com.bioinception.smsd.core.ChemOptions;
import com.bioinception.smsd.core.MolGraph;
import com.bioinception.smsd.core.SearchEngine;
import java.util.List;
import java.util.Map;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;

public class SearchExample {
    public static void main(String[] args) throws Exception {
        SmilesParser parser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer query = parser.parseSmiles("c1ccccc1");
        IAtomContainer target = parser.parseSmiles("c1ccc(O)cc1");
        MolGraph queryGraph = new MolGraph(query);
        MolGraph targetGraph = new MolGraph(target);
        ChemOptions chemistry = new ChemOptions();
        SearchEngine.MCSOptions options = new SearchEngine.MCSOptions();
        options.timeoutMs = 1000;
        options.connectedOnly = true;
        options.induced = false;

        Map<Integer, Integer> mapping =
                SearchEngine.findMCS(queryGraph, targetGraph, chemistry, options);
        if (mapping.size() != 6 || !SearchEngine.validateMapping(
                queryGraph, targetGraph, mapping, chemistry).isEmpty()) {
            throw new IllegalStateException("Expected a valid six-atom MCS");
        }
        for (Map.Entry<Integer, Integer> pair : mapping.entrySet()) {
            System.out.printf("query[%d] %s -> target[%d] %s%n",
                    pair.getKey(), query.getAtom(pair.getKey()).getSymbol(),
                    pair.getValue(), target.getAtom(pair.getValue()).getSymbol());
        }

        Map<Integer, Integer> first =
                SearchEngine.findSubstructure(query, target, chemistry, 1000);
        List<Map<Integer, Integer>> all =
                SearchEngine.findAllSubstructures(query, target, chemistry, 100, 1000);
        boolean hit = SearchEngine.isSubstructure(query, target, chemistry, 1000);
        if (!hit || first.size() != 6 || all.isEmpty()) {
            throw new IllegalStateException("Expected a benzene substructure match");
        }
        System.out.println("MCS atoms: " + mapping.size());
    }
}
```

Compile with Java 25, then run on Java 8 or later:

```bash
javac --release 8 -cp smsd-7.2.2-jar-with-dependencies.jar SearchExample.java
java -cp "smsd-7.2.2-jar-with-dependencies.jar:." SearchExample
```

When compiling with Java 8, omit `--release 8`. On Windows, use a semicolon
classpath separator:

```powershell
java -cp "smsd-7.2.2-jar-with-dependencies.jar;." SearchExample
```

`MolGraph` retains the supplied CDK atom order. Mappings, weights and exclusions
use zero-based atom indices, rather than atom-map numbers or canonical-SMILES
positions. Prepare aromaticity, hydrogens and stereo consistently when loading
containers through other CDK readers. If preprocessing changes atom order,
retain the correspondence to the original molecule.

The `SMSD` facade accepts SMILES or CDK containers and standardises inputs by
default. See the complete facade example in the [main README](../README.md).
Its four-argument CDK constructor accepts `standardise=false` for inputs you
have already prepared.

## Chemistry and search options

`ChemOptions` controls molecular matching; `SearchEngine.MCSOptions` controls
the objective, connectivity and search budget.

| Option | Behaviour |
| --- | --- |
| `matchAtomType` | Preserve element identity; enabled by default |
| `matchFormalCharge`, `matchIsotope` | Require matching charges or isotopes when enabled |
| `useChirality`, `useBondStereo` | Compare atom or bond stereo when enabled |
| `matchBondOrder` | `STRICT`, `LOOSE` or `ANY` |
| `aromaticityMode` | `STRICT` or `FLEXIBLE` |
| `ringMatchesRingOnly`, `completeRingsOnly` | Restrict ring matching |
| `timeoutMs` | MCS budget in milliseconds |
| `connectedOnly` | Require a connected common substructure |
| `induced` | Also require mapped target edges to exist in the query |
| `maximizeBonds` | Select the bond objective instead of the atom objective |
| `atomWeights` | Finite weights, one per query atom |
| `excludedTargetAtoms` | Target atom indices unavailable for mapping |

Use `ChemOptions.tautomerProfile()` explicitly for relaxed tautomer chemistry.
Changing the objective or weights changes the selected mapping; an atom count
alone does not describe a weighted or bond-maximising result.

Non-induced MCS preserves query direction: mapped query edges must exist in
the target, while extra target edges are allowed. Target exclusions retain
original indices and do not alter the molecule's chemistry.

Timeouts cover the search stages together. A bounded search may return a valid
mapping without proving global optimality; its deadline does not guarantee a
minimum MCS size. Stereo matching checks R/S descriptors and fully mapped
ligand permutations. Unspecified stereo can act as a wildcard.

## Fingerprints, stereo and R-groups

Save this as `FeaturesExample.java`. Fingerprint similarities below compare
a molecule with itself; the R-group example decomposes phenol and toluene
against a benzene core:

```java
import com.bioinception.smsd.core.ChemOptions;
import com.bioinception.smsd.core.CIPAssigner;
import com.bioinception.smsd.core.FingerprintEngine;
import com.bioinception.smsd.core.MolGraph;
import com.bioinception.smsd.core.SearchEngine;
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
        long[] path = SearchEngine.pathFingerprint(phenol, 7, 2048);
        long[] common = SearchEngine.mcsFingerprint(phenol, new ChemOptions(), 7, 2048);
        long[] circular = FingerprintEngine.ecfp(phenol, 2, 2048);
        if (!SearchEngine.fingerprintSubset(path, path)
                || SearchEngine.mcsFingerprintSimilarity(common, common) != 1.0
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
        List<Map<String, IAtomContainer>> groups =
                SearchEngine.decomposeRGroups(core, molecules, new ChemOptions(), 1000);
        if (groups.size() != 2) {
            throw new IllegalStateException("Expected two decomposition rows");
        }
        System.out.println("R-group rows: " + groups.size());
    }
}
```

```bash
javac --release 8 -cp smsd-7.2.2-jar-with-dependencies.jar FeaturesExample.java
java -cp "smsd-7.2.2-jar-with-dependencies.jar:." FeaturesExample
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
result. File types are `MOL`,
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
