# SMSD Pro Java guide

This guide describes the public Java API of SMSD Pro 8.0.0 in the package
`com.bioinception.smsd.core` (Maven artefact `com.bioinceptionlabs:smsd`), and
the command-line interface (CLI) in `com.bioinception.smsd.cli.SMSDcli`
(`com.bioinceptionlabs:smsd-cli`). The API provides substructure search,
maximum common substructure (MCS) search, SMARTS queries, fingerprints, stereo
descriptors, standardisation and R-group decomposition. It accepts CDK 2.13
`IAtomContainer` molecules or SMSD `MolGraph` graphs. The generated
[Javadoc](#api-reference) lists every public member.

## Requirements and installation

SMSD Pro is compiled with `--release 8` (class-file version 52), so the same
JARs run on Java 8 and later. Java 25 LTS is recommended for building and
running it. Applications that call the API depend on the library alone. Maven:

```xml
<dependency>
  <groupId>com.bioinceptionlabs</groupId>
  <artifactId>smsd</artifactId>
  <version>8.0.0</version>
</dependency>
```

Gradle (Kotlin DSL):

```kotlin
dependencies {
    implementation("com.bioinceptionlabs:smsd:8.0.0")
}
```

The library brings the CDK modules that its search, fingerprint and
standardisation code uses, and no file readers, SMILES parser, JSON library or
CLI. The examples below parse SMILES, so add `org.openscience.cdk:cdk-smiles:2.13`
to run them in your own project; the
[Java module README](../java/README.md#library) lists the other optional
additions, including `jackson-databind` for
`Standardiser.PredicateRegistry.loadJsonFile`. The
[command-line interface](#command-line-interface) is the separate artefact
`com.bioinceptionlabs:smsd-cli`. To build both from source, see
[Build from source](../java/README.md#build-from-source); the
[quick start](../java/README.md#quick-start) shows a complete substructure and
MCS search on CDK molecules.

The examples in this guide use the self-contained
`java/cli/target/smsd-cli-8.0.0-jar-with-dependencies.jar` from that build,
which includes the library, CDK and `cdk-smiles`. Save each example under its
class name, then compile and run it from the repository root:

```bash
javac --release 8 -cp java/cli/target/smsd-cli-8.0.0-jar-with-dependencies.jar PreparedSearchExample.java
java -cp "java/cli/target/smsd-cli-8.0.0-jar-with-dependencies.jar:." PreparedSearchExample
```

With a Java 8 compiler, omit `--release 8`. On Windows, separate the classpath
entries with a semicolon instead of a colon.

## Input graphs and atom indices

`MCS.find`, `Substructure.find`, `match` and `matchAll` accept either a CDK
`IAtomContainer` or a `MolGraph`. `new MolGraph(container)` keeps the CDK atom
order. Mappings, atom weights and exclusions all use these zero-based atom
indices, not atom-map numbers or canonical SMILES positions. If preprocessing
changes the atom order, keep the correspondence to the original molecule.

Searches use the aromaticity, hydrogen counts and stereochemistry supplied
with the input, so prepare them consistently. For example,
`Standardiser.standardise(container, Standardiser.TautomerMode.NONE)` returns a
standardised copy with perceived atom types, implicit hydrogens, rings and
Daylight aromaticity. `MolGraph.Builder` creates a graph without CDK: it
requires `atomCount`, `atomicNumbers` and `neighbors`, and optionally accepts
`formalCharges`, `massNumbers`, `hydrogenCounts`, `aromaticFlags`, `ringFlags`,
`bondOrders`, `bondAromaticFlags`, `bondRingFlags`, `tetrahedralChirality`,
`doubleBondStereo`, `atomIds`, `name`, `programLine`, `comment` and
`properties`; `build()` returns the graph. A graph exposes its per-atom arrays
`atomicNum`, `formalCharge`, `massNumber`, `label`, `tautomerClass` and
`atomId`, the methods `atomCount()`, `hydrogenCount(atom)` and `hasBond(i, j)`,
and the MOL header fields `name`, `programLine`, `comment` and `properties`.

Prepared queries hold their own copies of the query and its options, so they
can be reused, including from several threads at once. Keep a target graph
unchanged while a search is using it.

## Prepared searches

| Call | Returns |
| --- | --- |
| `MCS.find(query)`, `Substructure.find(query)` | A prepared query with default options |
| `Substructure.smarts(expression)` | A prepared SMARTS query, including recursive and stereo predicates |
| `withOptions(options)`, `withTimeout(milliseconds)` | A new prepared query; the options are copied |
| `match(target)` | One `SearchResult` with at most one match: the best MCS mapping found, or one embedding of the whole query |
| `matchAll(target)`, `matchAll(target, limit)` | All tied MCS optima, or all substructure embeddings, up to the limit (default 1,000) |
| `matchMany(targets)`, `matchMany(targets, numThreads)` | One `SearchResult` per `MolGraph` target, in input order |

MCS `matchAll` returns the mappings that tie on the whole objective (for the
default objective: atoms, then bonds). Mappings are distinct by their atom
correspondences and selected bonds, so symmetry-related mappings are reported
separately. An explicit limit must be positive.

`matchMany(targets)` searches the targets in turn. `matchMany(targets,
numThreads)` uses up to `numThreads` threads: 0 uses one thread per available
processor and 1 searches in turn, as in C++ and Python. Each target receives
the full per-call timeout, and stopped and negative results are included.

`Substructure.smarts` raises `IllegalArgumentException` for invalid SMARTS and
for expressions beyond these limits:

| Constant | Limit |
| --- | --- |
| `Substructure.MAX_SMARTS_LENGTH` | 100,000 characters |
| `Substructure.MAX_SMARTS_ATOMS` | 10,000 atoms, including those of recursive SMARTS |
| `Substructure.MAX_SMARTS_EXPRESSION_NODES` | 100,000 atom and bond expression nodes |
| `Substructure.MAX_SMARTS_NODES_PER_EXPRESSION` | 10,000 nodes in one expression |
| `Substructure.MAX_SMARTS_RECURSION_DEPTH` | 20 levels of recursive SMARTS |

SMARTS atom and bond predicates, including any charge, aromaticity and stereo
constraints, come from the expression and follow Daylight SMARTS. Targets are
perceived with the Daylight aromaticity model before SMARTS matching.

`Standardiser.expandNamedPredicates(smarts, registry)` replaces `$name` tokens
with recursive SMARTS from a `Standardiser.PredicateRegistry`, which starts
with ten named predicates such as `isHalogen`, `isHBDonor` and `isHBAcceptor`;
`put(name, smarts)` adds a predicate, `contains(name)` and `get(name)` read
them, and `loadJsonFile(path)` reads a JSON object of names and SMARTS.
`Standardiser.matchAll(smarts, target)` is deprecated in favour of
`Substructure.smarts(smarts).matchAll(target)`.

## Search options

Set fields on `MCS.Options` or `Substructure.Options`, then pass the object to
`withOptions`. Invalid values raise `IllegalArgumentException` when the options
are attached.

| Field | Default | Meaning |
| --- | --- | --- |
| `chemistry` | `new ChemOptions()` | Atom and bond comparison; see [chemistry options](#chemistry-options) |
| `MCS.Options.objective` | `MCS.Objective.MCIS` | `MCIS`: maximum common induced subgraph; mapped atoms are bonded in the query exactly when they are bonded in the target, and the search maximises atoms, then bonds. `MCES`: maximum common edge subgraph; the search maximises common bonds, then their atoms, and never includes isolated atoms |
| `MCS.Options.connectivity` | `MCS.Connectivity.Connected` | `Connected` requires one common fragment; `Fragments` allows several |
| `MCS.Options.minFragmentAtoms`, `maxFragments` | 1, `Integer.MAX_VALUE` | Smallest fragment size and largest number of fragments; both must be at least 1 |
| `MCS.Options.atomWeights` | `null` | MCIS only: one finite weight per query atom; the search maximises the weight sum, then atoms, then bonds. `null` or an empty array selects the unweighted objective |
| `MCS.Options.excludedTargetAtoms` | Empty | Target atom indices that must not be mapped |
| `Substructure.Options.induced` | `false` | Also require unbonded query atom pairs to map to unbonded target atoms |
| `timeoutMillis` | `-1` | Timeout in milliseconds; a negative value selects the adaptive timeout, see below |
| `nodeLimit` | `0` | Maximum number of search nodes; 0 means unlimited |
| `cancellation` | `null` | A `BooleanSupplier` polled during the search; returning `true` cancels it |

With a negative `timeoutMillis`, MCS uses an adaptive timeout of
min(30,000, 500 + 2 × query atoms × target atoms) milliseconds, and substructure
search uses 10,000 milliseconds. A value of 0 stops before searching; a
positive value is the timeout in milliseconds. The timeout and the node limit
together form the search budget of a call. The budget covers the search only:
input parsing, graph construction and query preparation are not counted.
Interrupting the calling thread also cancels a search.

Java uses the supplied `double` weights without rounding. The C++ and Python
APIs score weights in units of 0.001, so weights that differ by less than
0.001 can tie there but not in Java.

## Chemistry options

`ChemOptions` controls how atoms and bonds are compared. Its fields are public;
the setters `withBondOrderMode`, `withAromaticityMode`, `withCompleteRingsOnly`,
`withMatchIsotope`, `withRingFusionMode`, `withTautomerAware`, `withPH` and
`withSolvent` set a field and return the same instance, and
`ChemOptions.copyOf(options)` makes an independent copy.

| Field | Default | Meaning |
| --- | --- | --- |
| `matchAtomType` | `true` | Compare elements |
| `matchFormalCharge` | `false` | Compare formal charges |
| `matchIsotope` | `false` | Compare isotope labels |
| `matchBondOrder` | `BondOrderMode.STRICT` | `STRICT` requires equal bond orders; `LOOSE` relaxes them, subject to aromaticity and ring constraints; `ANY` accepts any order |
| `aromaticityMode` | `AromaticityMode.FLEXIBLE` | `FLEXIBLE` lets aromatic bonds match single or double bonds; `STRICT` matches aromatic bonds only with aromatic bonds |
| `ringMatchesRingOnly` | `false` | Ring atoms and bonds match only ring atoms and bonds |
| `completeRingsOnly` | `false` | MCS keeps a ring atom or bond only when a wholly selected relevant ring covers it in both molecules; a shared fusion atom may be kept through one complete ring |
| `ringFusionMode` | `RingFusionMode.IGNORE` | `IGNORE`, `PERMISSIVE` or `STRICT` comparison of fused ring systems |
| `useChirality` | `false` | Compare tetrahedral stereochemistry |
| `useBondStereo` | `false` | Compare double-bond (E/Z) stereochemistry |
| `tautomerAware` | `false` | Treat tautomeric bonds as equivalent |
| `pH`, `solvent` | `7.4`, `Solvent.AQUEOUS` | Context for `MCS.computeTautomerConfidence` only; `Solvent` is `AQUEOUS`, `DMSO`, `METHANOL`, `CHLOROFORM`, `ACETONITRILE` or `DIETHYL_ETHER` |

The fields `useTwoHopNLF`, `useThreeHopNLF` and `useBitParallelFeasibility`
and their `with` setters are deprecated; they tune the Java substructure
search alone and never change a result.

The named profiles return new instances. `ChemOptions.profile(name)` accepts
case-insensitive names and raises `IllegalArgumentException` for an unknown name.

| Profile | Settings |
| --- | --- |
| `new ChemOptions()`, `ChemOptions.profile("default")` | The defaults above |
| `ChemOptions.profile("strict")` | Strict bond order and strict aromaticity; charge, isotope and stereo matching remain off |
| `ChemOptions.profile("compat-fmcs")`, `profile("compat-substruct")`, `profile("reaction")` | Loose bond order; the deprecated `ChemOptions.fmcsProfile()` returns `profile("compat-fmcs")` |
| `ChemOptions.profile("pharma")` | Formal-charge, ring-only and complete-ring matching |
| `ChemOptions.tautomerProfile()` | Tautomer-aware matching with loose bond order |

### Stereochemistry

With `useChirality`, an atom is a tetrahedral centre when it has four
neighbours and no implicit hydrogen, or three neighbours and either one
implicit hydrogen or, for N, P, S, As, Se, Sb and Te, a lone pair. The C++ core
applies the same rule. The tags passed to `MolGraph.Builder.tetrahedralChirality`
(1 for `@`, 2 for `@@`) wind the neighbours in their listed order, with the
implicit hydrogen or lone pair second when the first neighbour has a lower
index than the centre, and first otherwise.

Stereo matching compares the orientation of mapped ligands. A `Substructure`
search embeds the query, so with `useChirality` a tagged query centre needs a
tagged target centre of the same winding over the mapped ligands, and with
`useBondStereo` a configured query double bond needs a configured target bond
of the same relation, as a SMARTS chirality primitive or pair of directional
bonds does; an untagged or unconfigured query atom or bond matches any target.
In a partial common graph, a stereocentre is constrained once enough of its
ligands are selected to determine it, and an unspecified target atom or bond
acts as a wildcard. The C++ core applies the same rules.

### Tautomers

`ChemOptions.tautomerProfile()` relaxes tautomeric chemistry during the search.
`MCS.computeTautomerConfidence(query, target, match, chemistry)` then reports a
heuristic confidence for the tautomeric relaxations in a selected match, using
`chemistry.pH` (0 to 14) and `chemistry.solvent`. These two fields affect only
this helper, not MCIS or MCES optimisation. The value is 1 when the match needs
no relaxation or `tautomerAware` is off; it is neither a probability nor an MCS
score, and both graphs are left unchanged.

The confidence is the geometric mean, over the relaxed atoms with a tautomer
class in either molecule, of the smaller of the two atom weights. Each atom of
a class takes the weight (`MolGraph.TW_*`) of the first motif in the table
below whose pattern covers the atom, and the ionisable site and aqueous pKa
(`TPKA_*`) that the matches of that motif covering the atom carry. Each motif
takes every atom of all its matches before the next is tried, so the motif
chosen for an atom does not depend on the order of the atoms or of their
neighbours. The atoms of one match share a class, and matches that share an
atom join their classes, so atoms that interconvert through any chain of motifs
form one class, in the manner of the mobile-hydrogen groups of InChI but
without requiring a mobile hydrogen, and the classes hold exactly the atoms
that some motif covers. The partition depends on the graph alone: `ONC(C)=O`
and `CC(=O)NO` both hold one class of all five atoms, succinimide one class of
all seven, and 2-hydroxypyridine, written `Oc1ccccn1` or `n1ccccc1O`, one class
of the oxygen, the ring nitrogen, the carbon between them and its para carbon.
`MolGraph.tautomerClass` numbers the classes from 0 in increasing order of a
key, the sorted colours of their atoms, refined in the manner of Weisfeiler and
Leman from element, charge, mass number, hydrogen count, aromaticity, ring
membership, degree and motif cover; of classes that still share a key, which
the refinement cannot tell apart (in practice those a symmetry of the graph
exchanges), the one with the lowest atom index takes the next number and the
refinement resumes. The numbers are therefore canonical up to Weisfeiler–Leman
refinement, the same for every atom order except among classes that the
refinement cannot tell apart: the two carbonyl classes of `CC(=O)CCC(C)=O` are
numbered 0 and 1, the class of the lower atoms first. An atom outside every
class has class −1, weight 1 and no site. The weight is the relevance of the
dominant neutral tautomer of the motif and does not depend on pH.

| Order | Motif | Weight | Site | pKa |
| --- | --- | --- | --- | --- |
| 1 | Hydroxamic acid | 0.85 | acid | 9 |
| 2 | Imide | 0.92 | acid | 9.5 |
| 3 | Barbituric: a ring N between two ring carbonyls, whose N and carbons the imide motif takes first, so it takes the carbonyl oxygens | 0.80 | — | — |
| 4 | Urea | 0.96 | acid | 14 |
| 5 | 1,3-Diketone | 0.72 | acid | 9.5 |
| 6 | β-Enaminone | 0.85 | — | — |
| 7 | 1,5-Keto-enol | 0.75 | — | — |
| 8 | Pyridone/hydroxypyridine | 0.95 | acid | 11.7, or 8.6 for a pyrimidin-4(3H)-one |
| 9 | Pyridinone | 0.95 | acid | 11.7, or 8.6 for a pyrimidin-4(3H)-one |
| 10 | Hydroxypyrimidine, the pattern of row 8, so it takes no atom | 0.90 | — | — |
| 11 | Amide/imidic acid | 0.97 | acid | 16 |
| 12 | Keto/enol | 0.98 | acid | 20 |
| 13 | Thioamide | 0.94 | acid | 13 |
| 14 | Thione/thiol | 0.94 | acid | 10.5 |
| 15 | Thioamide/iminothiol | 0.95 | acid | 13 |
| 16 | Nitroso/oxime | 0.95 | — | — |
| 17 | Phenol/quinone methide | 0.92 | — | — |
| 18 | Phenol/quinone | 0.88 | acid | 10 |
| 19 | Tetrazole, 1,2,3-triazole, benzimidazole, pyrazole and imidazole N–H shifts | 0.50 | — | — |
| 20 | Guanidine, amidine | 0.50 | — | — |
| 21 | Phosphonate | 0.50 | — | — |
| 22 | Nitro/aci-nitro | 0.95 | — | — |
| 23 | Furanose/pyranose ring–chain | 0.60 | — | — |
| 24 | Lactol ring–chain, the pattern of row 23, so it takes no atom | 0.65 | — | — |
| 25 | Allyl shift | 0.70 | — | — |
| 26 | Selenoketone/selenol | 0.94 | — | — |

`MolGraph.tautomerWeight` holds that weight in every context (C++
`MolGraph::tautomerWeight` holds it in the graph's own solvent and pH),
`MolGraph.tautomerSite` the site (`TAUTOMER_SITE_ACID`, `TAUTOMER_SITE_BASE` or
`TAUTOMER_SITE_NONE`) and `MolGraph.tautomerPKa` the pKa. The helper applies
the solvent shift first; the pH then scales the weight of every atom with a
site by the neutral fraction of its class. With pKa_a the lowest acid pKa and
pKa_b the highest base pKa (that of BH+) among the sites of the class:

```text
f0 = 1 / (1 + 10^(pH - pKa_a) + 10^(pKa_b - pH))
weight = min(0.99, max(0.1, solvent-corrected weight × f0))
```

A class with acid sites only drops the base term, which leaves
1 / (1 + 10^(pH − pKa_a)), the fraction of the neutral acid HA; a class with
base sites only drops the acid term, which leaves 1 / (1 + 10^(pKa_b − pH)),
the fraction of the neutral base B. The library assigns no pKa to the motifs
marked — in the table, and the pH does not change their weights.

A match of a motif with a site carries it only when the atom of the match that
would release the proton bears a hydrogen: the N–O oxygen of a hydroxamic acid,
so an N-alkoxy amide such as a Weinreb amide has none; the nitrogen of an
imide, urea, amide, thioamide or thioamide/iminothiol; the alpha carbon of a
ketone or thione and the bridging carbon of a 1,3-diketone; the ring nitrogen
or the oxygen of a pyridone/hydroxypyridine or pyridinone; and a phenol oxygen
with no second heavy-atom neighbour, so anisoles, aryl esters and lactam or
lactone carbonyls have none. When that atom belongs to an earlier motif, the
match carries the site of that motif instead, and an atom that several matches
of its motif cover takes the lowest pKa among them. The
pyridone/hydroxypyridine and pyridinone motifs take the N–H pKa of 2-pyridone,
11.7, or 8.6 when the other ring neighbour of the nitrogen is a carbon bonded
to a second aromatic ring nitrogen, as in pyrimidin-4(3H)-one. No motif in the
table has a base site. The pKa values are aqueous and are applied in every
solvent, and a class with several sites uses its lowest acid pKa and its
highest base pKa.

For example, against its iminol `ON=C(C)O` the hydroxamic acid `ONC(C)=O`
scores 0.85 at pH 2, 0.829 at pH 7.4, 0.773 at pH 8, 0.425 at pH 9 and the 0.1
floor at pH 10 and 12, and the same with its atoms written `CC(=O)NO`, while
the N-alkoxy amide `CONC(C)=O`, which has no O–H, scores 0.85 at every pH
against `CON=C(C)O`. The 2-pyridone `[nH]1ccccc1=O`, whose class also holds the
carbonyl and para carbons that the phenol motif reaches without an O–H, takes
the N–H pKa of 11.7 and scores 0.914 at pH 2 and 7.4, 0.905 at pH 10 and 0.317
at pH 12 against `n1ccccc1O`; `O=c1cccc[nH]1` scores the same against
`Oc1ccccn1`, whose para carbon takes that pKa through the O–H of the
pyridone/hydroxypyridine motif. 4(3H)-Pyrimidinone `O=c1ccnc[nH]1`, of pKa 8.6,
scores 0.95 at pH 2, 0.894 at pH 7.4, 0.271 at pH 9 and the 0.1 floor at pH 10
and 12 against `Oc1ccncn1`. The amidine `CC(=N)N` against `CC(N)=N` scores 0.5
at every pH. The weights follow the single-precision steps of the native core,
so the confidences of Java, C++ and Python agree within 1e-12.

The three-argument forms `MCS.computeTautomerConfidence(query, target, mapping)`
are deprecated. They select every mapped bond present in both molecules and
return the four-argument value in water at pH 7.4.
`MCS.validateTautomerConsistency(query, target, mapping)` checks that mobile
N–H, O–H and S–H protons are conserved within each mapped tautomeric group.

## Results

`match`, `matchAll` and `matchMany` return `SearchResult` objects.

| Method | Meaning |
| --- | --- |
| `status()` | `Completed`, `TimedOut`, `Cancelled` or `LimitReached` |
| `optimal()`, `complete()` | Whether the result is proved, and whether the requested search finished; see below |
| `found()` | Whether at least one match was retained |
| `matches()` | The retained `SearchResult.Match` objects |
| `score()` | The score of the best match |
| `nodesVisited()`, `elapsedMillis()` | Search effort |

Each `SearchResult.Match` describes one mapping in the original atom indices
of the query and target.

| Method | Meaning |
| --- | --- |
| `atomMapping()` | Query atom index to target atom index |
| `bondMapping()` | Selected query bond to target bond; each `SearchResult.Bond` holds its lower atom index in `first()` and the higher in `second()` |
| `fragments()` | Sets of query atoms joined by selected bonds; MCIS can also keep isolated atoms when the fragment options allow them |
| `atomCount()`, `bondCount()`, `score()` | Size and score of the mapping |
| `toMolecule(query)`, `toSmiles(query)` | The common fragment as a CDK molecule or canonical SMILES, built from the selected bonds with cut valences capped |
| `toGraph(query)` | The common fragment as a `MolGraph` that keeps the original chemical properties, for further matching |

Pass the same query, in the same atom order, that was used to prepare the
search. `score()` is the atom count for unweighted MCIS, the weight sum for
weighted MCIS and the bond count for MCES. Secondary tie-breaks are reflected
in the returned matches, not in the score.

Read `status()`, `optimal()` and `complete()` separately; `found()` reports only
whether a match was retained. A completed MCS search proves its optimum,
including a negative result. A search stopped by its budget can return a valid
match without proving it optimal. `complete()` means that the requested search
finished: a successful `match` does not enumerate all equal matches, and a
capped MCS enumeration can have `optimal() == true` and `complete() == false`.
`LimitReached` means that an enumeration limit or the node limit was reached;
check the flags and the configured limits. A substructure `match` that finds
an embedding proves existence, and a completed search without one proves
absence. A substructure `matchAll` can prove existence without completing its
enumeration.

## Example: prepared queries and fragments

Save this as `PreparedSearchExample.java`. It finds an MCES and disconnected
common fragments, enumerates substructure embeddings with and without a cap,
searches a batch of targets and matches a SMARTS pattern:

```java
import com.bioinception.smsd.core.MCS;
import com.bioinception.smsd.core.MolGraph;
import com.bioinception.smsd.core.SearchResult;
import com.bioinception.smsd.core.Substructure;
import java.util.Arrays;
import java.util.List;
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
        if (!edges.optimal() || !edges.found()
                || edges.matches().get(0).bondCount() != 2) {
            throw new IllegalStateException("Expected a proved two-bond common path");
        }

        options.connectivity = MCS.Connectivity.Fragments;
        options.minFragmentAtoms = 1;
        MolGraph pieces = new MolGraph(parser.parseSmiles("CO.CN"));
        SearchResult fragments = MCS.find(pieces).withOptions(options).match(pieces);
        if (!fragments.optimal() || !fragments.found()
                || fragments.matches().get(0).fragments().size() != 2) {
            throw new IllegalStateException("Expected two common fragments");
        }

        Substructure query = Substructure.find(new MolGraph(parser.parseSmiles("CO")));
        MolGraph target = new MolGraph(parser.parseSmiles("OCO"));
        SearchResult all = query.matchAll(target, 10);
        List<SearchResult> batch = query.matchMany(Arrays.asList(target, chain));
        if (!all.complete() || all.matches().size() != 2
                || !batch.get(0).found() || batch.get(1).found()) {
            throw new IllegalStateException("Unexpected enumeration or batch result");
        }
        SearchResult capped = query.matchAll(target, 1);
        if (capped.status() != SearchResult.Status.LimitReached || capped.complete()) {
            throw new IllegalStateException("Expected an incomplete capped enumeration");
        }
        if (!Substructure.smarts("[OX2H]").match(
                parser.parseSmiles("CCO")).found()) {
            throw new IllegalStateException("Expected an alcohol SMARTS match");
        }
    }
}
```

## Batch, consensus and mapping utilities

| Method | Purpose |
| --- | --- |
| `MCS.matchConstrained(queries, target, options)` | Searches each query in turn, each with its own budget, excluding target atoms claimed by earlier matches, and returns one result per query; once the cancellation token is cancelled, the remaining results report `Cancelled` |
| `MCS.consensus(molecules, options, threshold)` | Finds a common fragment supported by at least `ceil(threshold × number of molecules)` inputs, with one budget for the whole call; `threshold` must lie in (0, 1], and atom weights are rejected |
| `MCS.validateMapping(query, target, match, options)` | Lists the ways in which a match violates the options; an empty list means that it is valid. `validateMapping(query, target, mapping, options)` checks an atom mapping with the bonds `mappingMatch` selects, as in C++ and Python; the `ChemOptions` overload is deprecated |
| `MCS.mappingMatch(query, target, mapping, options)` | The match of an atom mapping: every query bond between mapped atoms that the target bonds compatibly |
| `MCS.isMappingMaximal(query, target, match, options)` | Returns `true` when no valid extension of the match exists and `false` for an invalid match; raises `IllegalStateException` if the budget stops the proof. An overload takes an atom mapping |
| `MCS.fragment(query, atoms, bonds)` | Builds a `Match` from an explicit atom mapping and selected query bonds |
| `MCS.canonicalizeMapping(query, target, mapping)` | Returns a canonical representative of an atom mapping under the symmetries of both graphs |
| `MCS.similarityUpperBound(query, target, chemistry)` | An upper bound on atom-based MCS similarity from the atom labels alone: element, aromaticity under strict aromaticity and ring membership when ring atoms match ring atoms only, as in C++ and Python |
| `MCS.murckoScaffold(molecule)` | The Bemis-Murcko scaffold: the ring systems, the linkers between them and every atom joined to either by a double or triple bond, as a new molecule with each stripped substituent replaced by hydrogen; a stereocentre keeps its configuration when all its neighbours stay, and a double bond when each end keeps a ligand, as in the C++ core |

`MCS.consensus` returns `MCS.ConsensusResult`. `fragment()` holds the fragment
in the atom indices of input `sourceIndex()`, and `graph()` and `toSmiles()`
return it as a compact graph and as SMILES. `mappings()` is keyed by the index
of each supporting input; each match maps `graph()` atom indices to that
input's atom indices. `supportCount()`, `status()`, `complete()`,
`nodesVisited()` and `elapsedMillis()` describe the search. Each input seeds a
candidate in turn, as in the C++ core: the input itself, then its common
subgraph with each other input in input order, found with the search of
`match`; each candidate is embedded in every input to count its support, and
the best supported candidate by the objective is kept. The fragment is reduced
progressively until enough inputs support it, so `optimal()` is always `false`.
When nothing is found, `sourceIndex()` is -1, `fragment()` and `graph()` are
empty and `toSmiles()` is an empty string.

`MCS.canonicalizeMapping` canonicalises the atom correspondence only, not an
MCES bond selection. It raises `IllegalArgumentException` for a mapped atom
outside either graph, and `IllegalStateException` when the automorphism
generators of either graph are incomplete or when its time or state limits are
exceeded, with the messages of the native core.

## Fingerprints, stereo descriptors and R-groups

Save this as `FeaturesExample.java`. It compares fingerprints of a molecule
with themselves, assigns stereo descriptors, and decomposes phenol and toluene
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

`FingerprintEngine` provides these static methods. The `IAtomContainer` forms
standardise a copy of the molecule first. To compute several fingerprints of
one molecule, convert it once with `FingerprintEngine.toMolGraph(molecule)` and
pass the graph to the `MolGraph` overloads.

| Methods | Purpose |
| --- | --- |
| `pathFingerprint(molecule, pathLength, size)` | Path fingerprint |
| `mcsFingerprint(molecule, chemistry, pathLength, size)` | Path fingerprint for MCS screening; honours `tautomerAware` |
| `ecfp`, `fcfp`, `circularFingerprint(..., Mode)` | Circular fingerprints, with `FingerprintEngine.Mode.ECFP` or `FCFP`; radius 2 corresponds to ECFP4 |
| `ecfpCounts`, `fcfpCounts`, `circularFingerprintCounts` | Circular count vectors |
| `topologicalTorsion`, `topologicalTorsionCounts` | Topological torsion fingerprints |
| `tanimoto`, `dice`, `cosine`, `soergel`, `overlapCoefficient` | Bit-vector similarity; `soergel` is the distance 1 − Tanimoto |
| `tanimotoCounts`, `diceCounts`, `cosineCounts` | Count-vector similarity |
| `fingerprintSubset(query, target)` | Whether every query bit is set in the target |
| `toHex`, `fromHex`, `toBitSet`, `fromBitSet`, `toBinaryString`, `countsToArray` | Conversion |
| `analyzeFpQuality`, `suggestFingerprintSize` | Bit density and size diagnostics; `analyzeFingerprintQuality` is deprecated |

Two empty bit fingerprints have a Tanimoto similarity and an overlap
coefficient of 1.0; `dice` and `cosine` return 0.0 for them. A fingerprint
match is a screen, not a substitute for graph matching, and fingerprint
similarity and MCS overlap are different measures.

### Configurable fingerprints

`FingerprintOptions` is a specification shared with C++ and Python: the same
options give the same bits in all three. `FingerprintEngine.fingerprint(graph, options)`
returns a `long[]` bitset and `fingerprintCounts(graph, options)` an `int[]` of
bins; the `IAtomContainer` forms standardise the molecule first. The factories
`FingerprintOptions.ecfp(radius, size)`, `fcfp(radius, size)`,
`path(pathLength, size)` and `torsion(size)`, and their no-argument forms with
radius 2, path length 7 and 2,048 bits, give the classic configuration of each
family, which reproduces the bits of `ecfp`, `fcfp`, `pathFingerprint` and
`topologicalTorsion` exactly; `new FingerprintOptions()` is the classic ECFP4.
Every further switch extends the hashed terms without moving the classic ones.
Save this as `FingerprintOptionsExample.java`:

```java
import com.bioinception.smsd.core.FingerprintEngine;
import com.bioinception.smsd.core.FingerprintOptions;
import com.bioinception.smsd.core.MolGraph;
import java.util.Arrays;
import java.util.BitSet;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;

public class FingerprintOptionsExample {
    public static void main(String[] args) throws Exception {
        SmilesParser parser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        FingerprintOptions options = FingerprintOptions.ecfp(2, 2048);
        options.cipAtoms = true;   // R/S joins the atom environments after the first iteration
        options.ringSize = true;   // the smallest ring through each atom, as a class up to 8
        long[] left = FingerprintEngine.fingerprint(parser.parseSmiles("C[C@H](N)C(=O)O"), options);
        long[] right = FingerprintEngine.fingerprint(parser.parseSmiles("C[C@@H](N)C(=O)O"), options);
        if (Arrays.equals(left, right)) {
            throw new IllegalStateException("Enantiomers must differ with the CIP atom term");
        }

        FingerprintOptions framework = FingerprintOptions.ecfp();
        framework.scaffold = FingerprintOptions.Scaffold.MURCKO;
        long[] toluene = FingerprintEngine.fingerprint(parser.parseSmiles("Cc1ccccc1"), framework);
        long[] benzene = FingerprintEngine.fingerprint(parser.parseSmiles("c1ccccc1"), FingerprintOptions.ecfp());
        if (!Arrays.equals(toluene, benzene)) {
            throw new IllegalStateException("The framework of toluene is benzene");
        }

        MolGraph aspirin = FingerprintEngine.toMolGraph(parser.parseSmiles("CC(=O)Oc1ccccc1C(=O)O"));
        BitSet keys = BitSet.valueOf(FingerprintEngine.fingerprint(aspirin, FingerprintOptions.maccs()));
        System.out.println("MACCS keys set: " + keys.cardinality() + ", six-membered ring: " + keys.get(162));
    }
}
```

| Field | Meaning | ECFP | FCFP | PATH | TORSION |
| --- | --- | --- | --- | --- | --- |
| `family` | `ECFP`, `FCFP`, `PATH`, `TORSION` or `KEYS` | | | | |
| `radius`, `pathLength`, `size` | Radius (−1 for the whole component), path length in bonds, and bits or bins | 2, –, 2048 | 2, –, 2048 | –, 7, 2048 | –, –, 2048 |
| `element` | Atomic number | on | off | on | on |
| `degree` | Number of explicit neighbours | on | on | off | on |
| `valence` | Bond order sum; pi-electron count for `TORSION` | on | off | off | on |
| `formalCharge`, `isotope`, `hydrogenCount`, `tautomerClass`, `aromaticity` | Formal charge, mass number, attached hydrogens, tautomer class when assigned, aromatic atom | on | off | off | off |
| `bondOrder` | Bond order of the neighbour term, or of each path bond | on | on | off | off |
| `aromaticBond`, `ringBond` | Aromatic and ring bond terms | off | off | off | off |
| `ringMembership` | Ring atom | on | on | off | on |
| `ringSize`, `ringSizeCap` | Smallest SSSR ring through the atom, as a class; rings of the cap size or more share one class | off, 8 | off, 8 | off, 8 | off, 8 |
| `cipAtoms`, `cipBonds` | CIP R/S (and r/s) of tetrahedral centres and E/Z of double bonds, from `CIPAssigner` | off | off | off | off |
| `fragments` | `WHOLE`, `LARGEST` (the component with the most heavy atoms, ties resolved without regard to atom order) or `PER_FRAGMENT` (every atom term carries its component's class) | `WHOLE` | | | |
| `scaffold` | `NONE`, or `MURCKO`: the Bemis-Murcko framework, that is ring atoms, the linkers between rings and the atoms double- or triple-bonded to them, capped with hydrogen where side chains were; a molecule with no ring is used whole | `NONE` | | | |
| `keys` | The keys of the `KEYS` family, one bit each: `new FingerprintOptions.Key(smarts, minimum)` is set when the molecule holds at least `minimum` distinct matched atom sets | | | | |

FCFP always hashes the pharmacophore feature class first. The stereo terms are
the CIP descriptors of the whole molecule, so a remote substituent that changes
a centre's priorities changes its label. They join the circular hashes after
the first iteration, so radius 0 does not separate stereoisomers and the bits
of a stereoisomer pair differ only within radius − 1 bonds of the centre or
bond; the path family folds R/S into its atom terms, so single-atom paths
already differ, and E/Z into its bond terms, and the torsion family folds both
into the atom types. Double bonds in rings of fewer than eight atoms carry no
E/Z term. The framework keeps exocyclic carbonyls and imines, as the usual
convention does, whereas `MCS.murckoScaffold` removes every side chain; a
tetrahedral centre or double bond that stays stereogenic after capping keeps
its configuration. The scaffold is taken before the fragment choice. `FingerprintOptions.maccs()` is the key
family of the 166 MACCS structural keys, key k at bit k − 1, where keys 1, 125
and 166 count isotopes, aromatic rings and fragments; `maccsKeys()` lists the
definitions and `withKeys(keys)` builds a custom key family. `copy()` returns an
independent copy. Inconsistent options, and a SMARTS key that does not parse,
throw `IllegalArgumentException`; a key search that does not complete within its
time budget, or that reaches 100,000 embeddings before its count is settled,
throws `IllegalStateException` rather than returning bits that depend on the
machine. The MCS screening fingerprint is a separate method that these options
do not configure.

`CIPAssigner.assignRS(graph)` returns R or S, or r or s for pseudoasymmetric
centres, keyed by atom index. Duplicate atoms follow IUPAC 2013 Rule 1a
literally: at a phosphorus centre the oxo oxygen, which carries a duplicate
phosphorus, ranks above a hydroxy or alkoxy oxygen, so phosphinates and
phosphonates may receive the opposite descriptor from toolkits that rank it
lower. `CIPAssigner.assignEZ(graph)` returns E or Z for
each annotated double bond, keyed by a `long` bond key. The traversal depth is
bounded by `CIPAssigner.MAX_DEPTH`, so check the returned descriptors rather
than assuming that every stereocentre is assigned.

`Substructure.decomposeRGroups(core, molecules, options)`, or
`decomposeRGroups(core, molecules, chemistry, timeoutMillis)`, returns one
`Substructure.Decomposition` per molecule. `search()` holds the core search
result, including no-match and stopped outcomes; inspect it before using
`core()`, the matched atoms with the images of the query bonds only (so a chain
core matched in a ring gives a chain, as in C++ and Python), which is `null`
without a match, or `rgroups()`, the attached groups keyed `R1`, `R2` and so on.

## Graph utilities

| Methods | Purpose |
| --- | --- |
| `toCanonicalSmiles()`, `toCanonicalSmiles(atomMaps)` | Canonical SMILES that retains the supplied tetrahedral and double-bond stereochemistry. The overload adds positive atom-map labels keyed by zero-based atom index; the labels change the output only, not the canonical ranking or the graph |
| `computeRings()`, `numRings()`, `MolGraph.computeSSSR(graph)`, `MolGraph.layoutSSSR(graph)` | Smallest set of smallest rings, and the same rings ordered for layout |
| `computeRelevantCycles()`, `computeURFs()` | Relevant cycles and unique ring families |
| `ensureRingCounts()`, `getRingCounts()` | Per-atom counts of the relevant cycles that the ring search finds, computed on first use |
| `ringSystemOf(atom)`, `ringSystemCount()`, `ringSystemSig(id)` | Ring systems and their signatures |
| `MolGraph.sanitizeSmiles(smiles)` | Balances parentheses and removes unclosed ring openings before CDK parsing |
| `getCanonicalLabeling()`, `getCanonicalHash()`, `getOrbit()` | Canonical labelling, graph hash and symmetry classes, each returned as a copy; `getOrbits()` is deprecated |
| `getAutomorphismGenerators()`, `automorphismGeneratorsTruncated()` | Symmetry generators, and whether they are incomplete |
| `MolGraph.matchTemplate`, `MolGraph.forceDirectedLayout`, `MolGraph.stressMajorisation`, `MolGraph.reduceCrossings` | Two-dimensional layout with `MolGraph.Point2D` coordinates |
| `MolGraph.translateToAtomIds(mapping, query, target)` | Translates a mapping to the external atom IDs set with `Builder.atomIds` |

## Command-line interface

The CLI is the module `com.bioinceptionlabs:smsd-cli`. Its self-contained
`smsd-cli-8.0.0-jar-with-dependencies.jar` runs with `java -jar`; the
[Java module README](../java/README.md#command-line-tool) shows how to obtain it,
and the Maven and Gradle dependencies for calling `SMSDcli` from an
application. From the repository root, after a source build:

```bash
java -jar java/cli/target/smsd-cli-8.0.0-jar-with-dependencies.jar --version
java -jar java/cli/target/smsd-cli-8.0.0-jar-with-dependencies.jar --help
```

Search for benzene as a substructure of phenol, then find their MCS:

```bash
java -jar java/cli/target/smsd-cli-8.0.0-jar-with-dependencies.jar \
  --Q SMI --q "c1ccccc1" --T SMI --t "c1ccc(O)cc1" --json -
java -jar java/cli/target/smsd-cli-8.0.0-jar-with-dependencies.jar \
  --Q SMI --q "c1ccccc1" --T SMI --t "c1ccc(O)cc1" --mode mcs --json -
```

MCS mode defaults to a connected MCIS. Select MCES and allow disconnected
fragments explicitly:

```bash
java -jar java/cli/target/smsd-cli-8.0.0-jar-with-dependencies.jar \
  --Q SMI --q "CO.CN" --T SMI --t "CO.CN" \
  --mode mcs --objective MCES --fragments --timeout 1000 --json -
```

SMARTS queries use the query type `SIG` and require substructure mode. Here
the query requires a carbonyl group:

```bash
java -jar java/cli/target/smsd-cli-8.0.0-jar-with-dependencies.jar \
  --Q SIG --q "C=O" --T SMI --t "CC(=O)O" --json -
java -jar java/cli/target/smsd-cli-8.0.0-jar-with-dependencies.jar \
  --Q SIG --q "C=O" --T SMI --t "CCO" --json -
```

The first search reports `"exists": true` with exit status 0, and the second
reports `"exists": false` with exit status 1.

| Option | Meaning |
| --- | --- |
| `--Q`, `--q` | Query type and value: `SMI` (SMILES) or `SIG` (SMARTS) with the string itself, or `MOL`, `ML2` (MOL2), `CML` or `PDB` with a file path |
| `--T`, `--t` | Target type and value: `SMI`, `MOL`, `ML2`, `CML` or `PDB`, or `SDF` for a batch file |
| `--mode` | `sub` (default) for substructure search or `mcs` for MCS search |
| `--objective` | `MCIS` (default) or `MCES`, in MCS mode |
| `--fragments` | Allow disconnected common fragments, in MCS mode |
| `--induced` | Require an induced substructure embedding |
| `-m`, `--mappings` | Enumerate up to 10,000 substructure embeddings |
| `--timeout` | Timeout per pair in milliseconds; default 10,000 |
| `-a`, `--add-hydrogens`; `-r`, `--remove-hydrogens` | Add or remove explicit hydrogens before matching |
| `--json`, `--json-pretty` | Write JSON to a file, or to standard output with `-`; pretty-print it |
| `--map-out`, `--map-format`, `--map-pretty` | Export mappings; see below |
| `-o`, `--output` | Write the MCS fragment as SMILES to a file, or to standard output with `-`, in MCS mode |
| `--threads` | Worker threads for an SDF batch; 0 (default) uses all available processors |
| `--benchmark` | Time substructure and MCS search on one pair and print the average times |

The CLI standardises each molecule (atom typing, implicit hydrogens, ring
perception and Daylight aromaticity) before searching. CML and PDB files must
contain exactly one molecule or model, empty inputs are rejected, and an SDF
batch may contain up to 100,000 molecules. The CLI timeout defaults to
10,000 milliseconds per pair, whereas the library's MCS default is adaptive.

JSON output is written as UTF-8. It includes `status`, `optimal`, `complete`,
`score`, `nodes_visited` and `elapsed_ms`, and the `bonds` and `fragments` of
the first match when one is found. Substructure mode adds `exists`; MCS mode
adds `mcs_size` and `pairs`, and `mcs_smiles` when a match is found. An SDF
batch writes one row per target under `results`.

`--map-out` writes the mappings of a single-pair search to a file, or to
standard output with `-`: in MCS mode, and in substructure mode with
`--mappings`. `--map-format` selects `json` (default), `smiles`, `smarts` or
`mdl`. SMILES export writes the query and the target with matching atom-map
labels and retains their stereochemistry. In MCS mode, `mdl` writes the
selected common atoms and bonds of the first match as a MOL block.

| Exit status | Single search | SDF batch |
| --- | --- | --- |
| 0 | The search completed with a match | At least one target matched, and every search completed |
| 1 | The search completed without a match | No target matched, and every search completed |
| 2 | Invalid input, an error or an incomplete search | Invalid input, or any target failed or did not complete |

Unknown modes are rejected with exit status 2.

## API reference

Generate the Javadoc for the current source from the repository root:

```bash
mvn -f java/pom.xml -pl core javadoc:javadoc
```

Then open `docs/api/index.html`. The generated files are not committed.

## Further reading

- [Java module README](../java/README.md): installation, quick start, build and tests.
- [Installation guide](HOWTO-INSTALL.md): platform requirements and
  [Java CLI installers](HOWTO-INSTALL.md#java-cli-installers).
- [Python guide](PYTHON.md) and [C++ guide](CPP.md): the same prepared-query API
  in the other languages.
- [Validation record](VALIDATION_8.0.0.md): release checks for 8.0.0.
- [Benchmark guide](../benchmarks/README.md): workloads and mapping validation.
- [Changelog](../CHANGELOG.md): what the release provides.

SMSD Pro is distributed under the [Apache License 2.0](../LICENSE). Retain
[NOTICE](../NOTICE) and the bundled dependency notices when redistributing it.
