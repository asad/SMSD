# SMSD Pro Java guide

Java prepared queries retain their CPU `matchMany` API. The optional native
corpus prefilters exposed in C++ and Python do not add a Java GPU search API.

The package `com.bioinception.smsd.core` (artefact `com.bioinceptionlabs:smsd`)
provides substructure, maximum common substructure (MCS) and SMARTS search,
fingerprints, stereo descriptors, standardisation and R-group decomposition
for CDK 2.13 `IAtomContainer` molecules and `MolGraph` graphs. The
command-line interface (CLI) is `com.bioinception.smsd.cli.SMSDcli`
(`com.bioinceptionlabs:smsd-cli`). The [Javadoc](#api-reference) lists every
public member.

In this release the Java library and command-line tool stay at 8.0.3 on Maven
Central. This guide covers the Java code in the 8.1.0 source, which a
[Java source build](HOWTO-INSTALL.md#java-source-build) turns into the library
and command-line tool; the 8.0.3 artefacts do not have its 8.1.0 changes. The
CPU [Docker images](HOWTO-INSTALL.md#docker) also carry the 8.1.0 command-line
tool.

## Installation

The JARs run on Java 8 and later; Java 25 LTS is recommended. See the
[supported-platforms table](HOWTO-INSTALL.md#supported-platforms) for the CLI
installers. Maven:

```xml
<dependency>
  <groupId>com.bioinceptionlabs</groupId>
  <artifactId>smsd</artifactId>
  <version>8.0.3</version>
</dependency>
```

Gradle (Kotlin DSL):

```kotlin
repositories {
    mavenCentral()
}

dependencies {
    implementation("com.bioinceptionlabs:smsd:8.0.3")
}
```

To parse SMILES, as the examples do, add `org.openscience.cdk:cdk-smiles:2.13`.
The [Java module README](../java/README.md) lists the other
[optional dependencies](../java/README.md#library), such as `jackson-databind`
for `Standardiser.PredicateRegistry.loadJsonFile`, and covers
[building from source](../java/README.md#build-from-source) and a
[quick start](../java/README.md#quick-start). To run an example, save it under
its class name and compile it from the repository root against the
self-contained JAR of the source build:

```bash
javac --release 8 -cp java/cli/target/smsd-cli-8.1.0-jar-with-dependencies.jar PreparedSearchExample.java
java -cp "java/cli/target/smsd-cli-8.1.0-jar-with-dependencies.jar:." PreparedSearchExample
```

With a Java 8 compiler, omit `--release 8`; on Windows, separate the classpath
entries with a semicolon.

## SMARTS export in 8.1.0

`MCS.matchWithSmarts(target)` returns `MCS.SmartsResult`; `result()` retains the
ordinary search status and mappings, and `smarts()` is the best match's portable
pattern under the prepared chemistry. `SearchResult.Match.toSmarts(query,
chemistry)` and `SmartsWriter.toSmarts(graph, chemistry)` provide explicit export.
See [SMARTS export](SMARTS_EXPORT.md) for a library-search example, selected bond
semantics and rejected options. Stereo, tautomer-aware, complete-ring and strict
ring-fusion constraints cannot be silently discarded. Portable SMARTS does not
encode the induced MCIS objective, a resource budget or an optimality proof.

## Input graphs

Searches accept an `IAtomContainer` or a `MolGraph`. `new MolGraph(container)`
keeps the CDK atom order; mappings, atom weights and exclusions use these
zero-based indices. Aromaticity, hydrogen counts and stereochemistry are used
as given; `Standardiser.standardise(container, Standardiser.TautomerMode.NONE)`
returns a copy with atom types, implicit hydrogens, rings and Daylight
aromaticity. The CDK `SmilesParser` gives atoms written in capitals but joined
by `:` ring bonds, as in `C1:C:C:C:C:C:1`, the hydrogens of single bonds;
`MolGraph` and `Standardiser.standardise` read such a ring as aromatic, as the
native reader does, but the CDK molecule does not record which nitrogen a
bracket hydrogen was written on, so where several could carry it, one that
gives a Kekulé structure does, whatever the atom order. The CLI reads such
SMILES as written.

`MolGraph.Builder` builds a graph without CDK from `atomCount`,
`atomicNumbers` and `neighbors`, optionally with `formalCharges`,
`massNumbers`, `hydrogenCounts`, `aromaticFlags`, `ringFlags`, `bondOrders`,
`bondAromaticFlags`, `bondRingFlags`, `tetrahedralChirality`,
`doubleBondStereo`, `atomIds`, `name`, `programLine`, `comment` and
`properties`; `build()` returns the graph. A graph exposes the per-atom arrays
`atomicNum`, `formalCharge`, `massNumber`, `label`, `tautomerClass` and
`atomId`, the methods `atomCount()`, `hydrogenCount(atom)` and `hasBond(i, j)`,
and the MOL header fields `name`, `programLine`, `comment` and `properties`.

Prepared queries copy their query and options, so they can be reused and
shared between threads; do not modify a target during a search.

## Prepared searches

- `MCS.find(query)`, `Substructure.find(query)`: a prepared query with default
  options; `Substructure.smarts(expression)`: a SMARTS query, including
  recursive and stereo predicates.
- `withOptions(options)`, `withTimeout(milliseconds)`: a new prepared query
  with copied options.
- `match(target)`: a `SearchResult` with the best MCS or one embedding.
- `matchAll(target)`, `matchAll(target, limit)`: all tied MCS optima,
  symmetric ones included, or all embeddings, up to a positive limit (default
  1,000).
- `matchMany(targets)`, `matchMany(targets, numThreads)`: one `SearchResult`
  per `MolGraph` target, in order, each with the full timeout; 0 threads uses
  every processor, 1 searches in turn.

`Substructure.smarts` follows Daylight SMARTS and aromaticity and raises
`IllegalArgumentException` for invalid SMARTS or SMARTS beyond
`Substructure.MAX_SMARTS_LENGTH`, `MAX_SMARTS_ATOMS`,
`MAX_SMARTS_EXPRESSION_NODES`, `MAX_SMARTS_NODES_PER_EXPRESSION` or
`MAX_SMARTS_RECURSION_DEPTH`.

`Standardiser.expandNamedPredicates(smarts, registry)` replaces `$name` tokens,
where a name is ASCII letters, digits and underscores, inside recursive SMARTS
as well, with SMARTS from a
`Standardiser.PredicateRegistry`, which offers
`put(name, smarts)`, `contains(name)`, `get(name)` and `loadJsonFile(path)`
(a JSON object of names and SMARTS). A definition is expanded as recursive
SMARTS, so it describes the atom written first. The registry starts with ten
predicates, defined as in the native core:

- `isAmideN`: N with three connections bonded to a carbonyl carbon (amides,
  lactams, imides, ureas and carbamates), aromatic rings such as uracil and
  2-pyridone included.
- `isCarboxylC`: carbonyl carbon bearing OH or O⁻ (carboxylic acids and
  carboxylates, carbamic acids, carbonic acid monoesters and bicarbonate).
- `isSulfonamideN`: N with three connections bonded to a sulfonyl S.
- `isEster`: carbonyl carbon of a carboxylic ester, formates and aromatic
  lactones such as coumarin included; anhydrides are excluded.
- `isKetone`: carbonyl carbon bonded to two carbons, ring carbonyls such as
  chromone and 4-pyridone included.
- `isHalogen`: F, Cl, Br or I.
- `isPositive`, `isNegative`: a charge of +1 to +15, or −1 to −15.
- `isHBDonor`: N, O or S bearing hydrogen, aromatic NH included.
- `isHBAcceptor`: N, O or S of charge 0 to −3 with a free lone pair, aniline,
  nitrile, thioamide, amidine and sulfonamide N included; amide N, aromatic N
  with three connections (pyrrole type), N of valence 5 (nitro and N-oxide
  written without charges) and S of valence 4 or 6 (sulfoxides, sulfones) are
  excluded.

`Standardiser.matchAll(smarts, target)` is deprecated in favour of
`Substructure.smarts(smarts).matchAll(target)`.

## Search options

Set fields on `MCS.Options` or `Substructure.Options` and pass the object to
`withOptions`; invalid values raise `IllegalArgumentException`. Defaults are
in brackets.

- `chemistry` (`new ChemOptions()`): see [chemistry options](#chemistry-options).
- `MCS.Options.objective` (`MCS.Objective.MCIS`): `MCIS`, the maximum common
  induced subgraph, maximises atoms, then bonds; `MCES`, the maximum common
  edge subgraph, maximises bonds, then atoms, with no isolated atoms.
- `MCS.Options.connectivity` (`MCS.Connectivity.Connected`): one common
  fragment, or several with `Fragments`; `minFragmentAtoms` (1) and
  `maxFragments` (`Integer.MAX_VALUE`), both at least 1, bound their size and
  number.
- `MCS.Options.atomWeights` (`null`): MCIS only; one finite weight per query
  atom, maximising the weight sum first. `null` or empty is unweighted.
- `MCS.Options.excludedTargetAtoms` (empty): target atoms never mapped.
- `Substructure.Options.induced` (`false`): unbonded query atom pairs must map
  to unbonded target atoms.
- `timeoutMillis` (`-1`): negative selects the default, which for MCS is
  500 ms plus 2 ms per atom pair (query atoms × target atoms), at most
  30,000 ms: 1,300 ms for two 20-atom molecules, 5,500 ms for two 50-atom
  molecules and 30,000 ms from 122 atoms each. For substructure search it is
  10,000 ms. 0 stops at once; positive is milliseconds.
- `nodeLimit` (`0`, unlimited): the maximum number of search nodes.
- `cancellation` (`null`): a `BooleanSupplier`; returning `true` cancels.

The timeout and node limit cover the search, not parsing or preparation;
interrupting the calling thread also cancels it. Java uses atom weights
unrounded, whereas C++ and Python score them in units of 0.001.

## Chemistry options

`ChemOptions` fields are public. `withBondOrderMode`, `withAromaticityMode`,
`withCompleteRingsOnly`, `withMatchIsotope`, `withRingFusionMode`,
`withTautomerAware`, `withPH` and `withSolvent` set a field and return the same
instance; `ChemOptions.copyOf(options)` returns an independent copy.

- `matchAtomType` (`true`), `matchFormalCharge`, `matchIsotope` (`false`):
  compare elements, formal charges and isotope labels.
- `matchBondOrder` (`BondOrderMode.STRICT`): `STRICT` needs equal orders,
  `LOOSE` relaxes them within the aromaticity and ring settings, `ANY` accepts
  any.
- `aromaticityMode` (`AromaticityMode.FLEXIBLE`): aromatic bonds also match
  single and double bonds; `STRICT` matches them only with aromatic bonds.
- `ringMatchesRingOnly` (`false`): ring atoms and bonds match only ring atoms
  and bonds. `completeRingsOnly` (`false`): MCS keeps ring atoms and bonds only
  within relevant cycles selected whole in both molecules.
- `ringFusionMode` (`RingFusionMode.IGNORE`): `IGNORE`, `PERMISSIVE` or
  `STRICT`; strict matching requires matched ring atoms to have the same
  per-atom ring count.
- `useChirality`, `useBondStereo` (`false`): mapped tetrahedral centres and
  double bonds must agree in configuration, compared by the parity of their
  mapped neighbours, not by CIP labels (R/S, E/Z); see
  [stereochemistry](#stereochemistry).
- `tautomerAware` (`false`): treat tautomeric bonds as equivalent.
- `pH` (`7.4`), `solvent` (`Solvent.AQUEOUS`; also `DMSO`, `METHANOL`,
  `CHLOROFORM`, `ACETONITRILE`, `DIETHYL_ETHER`): used only by
  `MCS.computeTautomerConfidence`.

The deprecated fields `useTwoHopNLF`, `useThreeHopNLF` and
`useBitParallelFeasibility`, and their setters `withTwoHopNLF`,
`withThreeHopNLF` and `withBitParallelFeasibility`, never change a result.

Profiles return new instances; `ChemOptions.profile(name)` ignores case and
raises `IllegalArgumentException` for an unknown name:

- `new ChemOptions()`, `ChemOptions.profile("default")`: the defaults above.
- `profile("strict")`: strict bond order and aromaticity.
- `profile("compat-fmcs")`, `profile("compat-substruct")`, `profile("reaction")`:
  loose bond order. The deprecated `ChemOptions.fmcsProfile()` returns
  `profile("compat-fmcs")`.
- `profile("pharma")`: formal-charge, ring-only and complete-ring matching.
- `ChemOptions.tautomerProfile()`: tautomer-aware matching, loose bond order.

### Stereochemistry

With `useChirality`, a tetrahedral centre has four neighbours and no implicit
hydrogen, or three neighbours and one implicit hydrogen or, for N, P, S, As,
Se, Sb and Te, a lone pair. `MolGraph.Builder.tetrahedralChirality` takes 1 for
`@` and 2 for `@@` over the listed neighbours, with the implicit hydrogen or
lone pair second when the first neighbour's index is below the centre's, and
first otherwise. A configured query centre or double bond needs the same
configuration in the target (in MCS, once enough neighbours are mapped to
define it); unconfigured query atoms and bonds, and in MCS unspecified target
stereochemistry, match any. The test compares the parity of the mapped
neighbours, the `@` or `@@` order around a centre and the cis or trans
placement across a double bond, not CIP labels, whose priorities depend on
atoms outside the mapping.

In MCS, stereochemistry is a hard constraint. An MCIS is the largest common
subgraph whose mapped centres and double bonds agree, so it can be smaller
than the MCIS found without these options; an MCES instead leaves out the
fewest common bonds that resolve each disagreement.

### Tautomers

`MCS.computeTautomerConfidence(query, target, match, chemistry)` scores the
tautomeric relaxations of a match found with `ChemOptions.tautomerProfile()`,
at `chemistry.pH` (0 to 14) in `chemistry.solvent`: 1 when none is needed or
`tautomerAware` is off, lower for less favoured ones. It is not a probability
or an MCS score, changes neither graph and agrees across Java, C++ and Python.
The deprecated three-argument form takes an atom mapping, selects every
mapped bond common to both molecules and assumes water at pH 7.4.
`MCS.validateTautomerConsistency(query, target, mapping)` checks that mobile
N–H, O–H and S–H protons are conserved within each mapped tautomeric group.

Per atom, `MolGraph.tautomerClass` holds the tautomer class (from 0,
independent of atom order except among symmetric classes; −1 outside every
class), `tautomerWeight` its `MolGraph.TW_*` weight (1 outside every class),
`tautomerSite` its ionisable site (`TAUTOMER_SITE_ACID`, `TAUTOMER_SITE_BASE`
or `TAUTOMER_SITE_NONE`) and `tautomerPKa` its aqueous `MolGraph.TPKA_*` pKa.

## Results

A `SearchResult` reports:

- `status()`: `Completed` (the search ran to its end), `TimedOut` (the time
  limit stopped it), `Cancelled` (cancellation or an interrupt stopped it) or
  `LimitReached` (the node limit stopped it, an enumeration reached its match
  limit, or an MCES pair beyond the [size limits](#size-limits) was declined
  before any search, which returns no match and proves nothing).
- `optimal()`: the result is proved. `complete()`: the requested search
  finished. `found()`: a match was retained.
- `matches()`, `score()`: the `SearchResult.Match` objects and the best score.
- `nodesVisited()`, `elapsedMillis()`: search effort.

`completeRingsOnly` requires complete relevant-cycle perception;
`RingFusionMode.STRICT` requires complete perception for its per-atom ring
counts. If that perception is incomplete in either input, the selected mappings
remain, but the result is `LimitReached`, with `optimal()` and `complete()` both
`false`; `TimedOut` and `Cancelled` take precedence. Positive numeric SMARTS
`R<n>` (`n > 0`) uses one deterministic minimum cycle basis, the same as C++
and Python, so a result that 8.0.3 computed from a different minimum basis of
equal size can change. If ring perception cannot complete the basis, a search
that cannot be decided without it returns `LimitReached` without matches,
including for negated, combined and recursive predicates; `found()` is then
false, which does not prove absence. Literal `R` and `R0` remain independent
of that basis. Other functions that need the basis throw
`RingBasisUnavailableException`, an `IllegalStateException`, naming the cause;
`computeSSSR` and `computeRings` may return incomplete rings with truncation
flags. `relevantCyclesTruncated` alone does not make the basis incomplete. Ring
perception works within one fixed, deterministic work limit; see
[Results and budgets](ALGORITHMS.md#results-and-budgets).

A `SearchResult.Match` uses the original atom indices; pass it the query that
prepared the search:

- `atomMapping()`: query atom to target atom.
- `bondMapping()`: query bond to target bond, as `SearchResult.Bond` pairs
  with the lower atom index in `first()` and the higher in `second()`.
- `fragments()`: sets of query atoms joined by selected bonds.
- `atomCount()`, `bondCount()`, `score()`: the score is the atom count for
  MCIS, the weight sum for weighted MCIS and the bond count for MCES.
- `toMolecule(query)`, `toSmiles(query)`, `toGraph(query)`: the common
  fragment as a CDK molecule, as canonical SMILES (cut valences capped), or as
  a `MolGraph` with the original chemical properties.

A completed search proves its result, including absence, and a substructure
match proves existence; a stopped search may return a valid but unproved
match. A capped MCS enumeration can be optimal yet incomplete. The C++ and
Python MCIS search changed in 8.1.0 and Java's did not, so a stopped search can
return a different result in Java; see the
[upgrade notes](RELEASE_NOTES.md#upgrading-from-803). In rare cases Java
reports `Cancelled`, with `optimal()` true and `complete()` false, for a search
whose proof had already completed; the score and the proof are unaffected.

### Size limits

An MCES search declines a pair before searching it, with `LimitReached`, no
match, and `optimal()` and `complete()` both `false`, when:

- query atoms × target atoms exceed 2³⁰ (1,073,741,824), as for two molecules
  of more than 32,768 atoms each;
- the two molecules are too large in bonds, as for two molecules of more than
  16,384 bonds each;
- rarely, a molecule has an atom with 512 or more bonds and the two molecules
  share many kinds of bond.

MCIS searches are never declined; the time and node limits bound them. The
C++ and Python libraries use the same limits.

## Example: prepared queries and fragments

`PreparedSearchExample.java` compares an MCIS with an MCES, finds common
fragments, enumerates embeddings with and without a cap, searches a batch and
matches SMARTS. Cyclopropane and isobutane share two atoms as an MCIS, as any
three ring atoms form a triangle, but a two-bond path as an MCES; diethyl ether
and diethylamine share two separate ethyl groups:

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
        MolGraph ring = new MolGraph(parser.parseSmiles("C1CC1"));       // cyclopropane
        MolGraph branched = new MolGraph(parser.parseSmiles("CC(C)C"));  // isobutane
        SearchResult induced = MCS.find(ring).match(branched);
        MCS.Options options = new MCS.Options();
        options.objective = MCS.Objective.MCES;
        SearchResult edges = MCS.find(ring).withOptions(options).match(branched);
        if (!induced.optimal() || !induced.found()
                || induced.matches().get(0).atomCount() != 2
                || !edges.optimal() || !edges.found()
                || edges.matches().get(0).bondCount() != 2) {
            throw new IllegalStateException("Expected a two-atom MCIS and a two-bond MCES");
        }

        options.connectivity = MCS.Connectivity.Fragments;
        options.minFragmentAtoms = 1;
        MolGraph ether = new MolGraph(parser.parseSmiles("CCOCC"));
        MolGraph amine = new MolGraph(parser.parseSmiles("CCNCC"));
        SearchResult fragments = MCS.find(ether).withOptions(options).match(amine);
        if (!fragments.optimal() || !fragments.found()
                || fragments.matches().get(0).fragments().size() != 2) {
            throw new IllegalStateException("Expected two common ethyl fragments");
        }

        Substructure query = Substructure.find(new MolGraph(parser.parseSmiles("CO")));
        MolGraph target = new MolGraph(parser.parseSmiles("OCO"));
        SearchResult all = query.matchAll(target, 10);
        List<SearchResult> batch = query.matchMany(Arrays.asList(target, branched));
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

## Mapping and consensus utilities

- `MCS.matchConstrained(queries, target, options)`: one result per query,
  searched in turn with its own budget, avoiding target atoms that earlier
  matches claimed; `Cancelled` for the rest after cancellation.
- `MCS.consensus(molecules, options, threshold)`: a fragment shared by at least
  `ceil(threshold × number of molecules)` inputs, under one budget; `threshold`
  lies in (0, 1], and atom weights are rejected. Progressive reduction is an
  approximation; it does not prove a global multi-molecule optimum.
- `MCS.validateMapping(query, target, match, options)`: the ways a match
  violates the options, empty if valid. The atom-mapping overload checks the
  bonds `mappingMatch` selects; the `ChemOptions` overload is deprecated.
- `MCS.mappingMatch(query, target, mapping, options)`: an atom mapping as a
  match, with every compatible query bond between mapped atoms.
- `MCS.isMappingMaximal(query, target, match, options)`, also for an atom
  mapping: `true` if no valid extension exists, `false` for an invalid match,
  `IllegalStateException` if the budget runs out or required ring perception
  was truncated.
- `MCS.fragment(query, atoms, bonds)`: a `Match` from an atom mapping and query
  bonds.
- `MCS.canonicalizeMapping(query, target, mapping)`: the canonical form of an
  atom mapping under the symmetries of both graphs, ignoring bond selections;
  `IllegalArgumentException` for an atom outside either graph,
  `IllegalStateException` if a graph's symmetry cannot be fully determined.
- `MCS.similarityUpperBound(query, target, chemistry)`: an upper bound on
  atom-based MCS similarity from atom labels.
- `MCS.murckoScaffold(molecule)`: the Bemis-Murcko scaffold (ring systems,
  linkers and atoms double- or triple-bonded to them) as a new molecule,
  hydrogen-capped, keeping stereochemistry that stays defined.

`MCS.ConsensusResult` holds `fragment()` in the atom indices of input
`sourceIndex()`, the same fragment as `graph()` and `toSmiles()`, and
`mappings()` from `graph()` atoms to each supporting input; `supportCount()`,
`status()`, `complete()`, `nodesVisited()` and `elapsedMillis()` describe the
search. `optimal()` is always `false`, including empty and single-input calls.
`complete()` reports whether the progressive pass finished. Neither it nor
`Completed` proves a global optimum or absence of a qualifying common graph.
With no result, `sourceIndex()` is -1 and the fragment, graph and SMILES are
empty.
With complete rings or strict ring fusion, truncated ring perception in any
input also makes consensus `LimitReached` and incomplete, unless its status
is already `TimedOut` or `Cancelled`. Its fragment and supporting mappings are
retained.

## Fingerprints, stereo descriptors and R-groups

`FeaturesExample.java` compares fingerprints, assigns stereo descriptors and
decomposes phenol and toluene against benzene:

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

`FingerprintEngine` methods are static; the `IAtomContainer` forms standardise
a copy first, so for several fingerprints convert once with
`FingerprintEngine.toMolGraph(molecule)`:

- `pathFingerprint(molecule, pathLength, size)`; for MCS screening,
  `mcsFingerprint(molecule, chemistry, pathLength, size)`, which honours
  `tautomerAware`.
- `ecfp`, `fcfp`, `circularFingerprint(..., Mode)` (`FingerprintEngine.Mode.ECFP`
  or `FCFP`; radius 2 is ECFP4), `ecfpCounts`, `fcfpCounts`,
  `circularFingerprintCounts`, `topologicalTorsion`, `topologicalTorsionCounts`.
- `tanimoto`, `dice`, `cosine`, `soergel` (1 − Tanimoto), `overlapCoefficient`,
  and for counts `tanimotoCounts`, `diceCounts`, `cosineCounts`. Two empty bit
  fingerprints give 1.0 from `tanimoto` and `overlapCoefficient`, 0.0 from
  `dice` and `cosine`.
- `fingerprintSubset(query, target)`: every query bit is set in the target; a
  screen, not a substitute for graph matching.
- `toHex`, `fromHex`, `toBitSet`, `fromBitSet`, `toBinaryString`,
  `countsToArray`: conversion.
- `analyzeFpQuality(fp, size)`: `set_bits`, `fp_size`, `density`,
  `expected_collision_rate` and `is_saturated` (density over 0.75), with the
  definitions of Python `analyze_fp_quality`: s of N bits set suggest
  n = ln(1 − s/N) / ln(1 − 1/N) uniformly hashed features, of which the share
  1 − s/n fell on a bit already set; path fingerprints of a power-of-two size
  collide more often than uniform hashing predicts.
  `analyzeFingerprintQuality` is deprecated.
- `suggestFingerprintSize(atoms, librarySize, rate)`: the path fingerprint size
  at which a path absent from a molecule finds its bit set with probability of
  about `rate`, allowing seven distinct paths per atom, more than most
  drug-like molecules have, and one bit per path: m = −n / ln(1 − rate),
  rounded up to a multiple of 64; `atoms` must be at least 1, `rate` must lie
  strictly between 0 and 1, and a size beyond the `int` range is rejected.

### Configurable fingerprints

`FingerprintOptions` gives the same bits in Java, C++ and Python.
`FingerprintEngine.fingerprint(graph, options)` returns a `long[]` bitset and
`fingerprintCounts(graph, options)` an `int[]` of bins. The factories
`FingerprintOptions.ecfp(radius, size)`, `fcfp(radius, size)`,
`path(pathLength, size)` and `torsion(size)`, and their no-argument forms,
reproduce `ecfp`, `fcfp`, `pathFingerprint` and `topologicalTorsion` exactly;
`new FingerprintOptions()` is ECFP4. Extra switches add terms without changing
the default ones, as `FingerprintOptionsExample.java` shows:

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
        options.cipAtoms = true;   // add the R/S descriptor of each centre
        options.ringSize = true;   // add the smallest ring size of each atom
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

Fields, with defaults and the families whose terms are on by default:

- `family`: `ECFP`, `FCFP`, `PATH`, `TORSION` or `KEYS`; `radius` (2; −1 for
  the whole component), `pathLength` in bonds (7), `size` in bits or bins
  (2,048).
- Atom terms: `element` (ECFP, PATH, TORSION), `degree` and `ringMembership`
  (ECFP, FCFP, TORSION), `valence`, the bond order sum or for TORSION the pi
  electrons (ECFP, TORSION), and `formalCharge`, `isotope`, `hydrogenCount`,
  `tautomerClass` and `aromaticity` (ECFP).
- Bond terms: `bondOrder` (ECFP, FCFP); `aromaticBond`, `ringBond` (off).
- `ringSize` (off): the smallest ring through each atom; rings of
  `ringSizeCap` (8) atoms or more share one class.
- `cipAtoms`, `cipBonds` (off): R/S (and r/s) of centres and E/Z of double
  bonds, from `CIPAssigner`. Circular fingerprints of radius 0 do not separate
  stereoisomers, and double bonds in rings below eight atoms carry no E/Z term.
- `fragments`: `WHOLE` (default), `LARGEST` (most heavy atoms) or
  `PER_FRAGMENT` (each component tagged). `scaffold`: `NONE` (default) or
  `MURCKO` (the scaffold of `MCS.murckoScaffold`; a molecule without rings is
  used whole).
- `keys`: for `KEYS`, one bit per `new FingerprintOptions.Key(smarts, minimum)`,
  set when at least `minimum` distinct atom sets match.

`FingerprintOptions.maccs()` gives the 166 MACCS keys, key k at bit k − 1;
`maccsKeys()` lists them, `withKeys(keys)` builds a custom key family and
`copy()` copies the options. Inconsistent options and unparsable SMARTS keys
raise `IllegalArgumentException`, and a key count not settled within its
budget raises `IllegalStateException`. `mcsFingerprint` ignores these options.

Ring-size terms, double-bond stereo terms and aromatic-ring keys use the same
minimum cycle basis as C++ and Python, so a Java value that came from a
different minimum basis of equal size can change. They throw
`RingBasisUnavailableException` naming the cause if ring perception cannot
complete that basis, as do SMARTS keys whose search needs it; a separate
relevant-cycle truncation does not prevent their use. The ring refusals replace values that
8.0.3 could compute from incomplete rings. A key whose minimum number of
matches is found before the embedding limit is set; a key search that runs out
of time throws, as in 8.0.3.

`CIPAssigner.assignRS(graph)` returns R or S (r or s if pseudoasymmetric) by
atom index, and `CIPAssigner.assignEZ(graph)` E or Z for each annotated double
bond by a `long` bond key. Rule 1a of IUPAC 2013 is applied literally, so
phosphinate and phosphonate centres may differ from some toolkits. Ranking is
bounded by `CIPAssigner.MAX_DEPTH`, so check which centres were assigned.

`Substructure.decomposeRGroups(core, molecules, options)`, or
`decomposeRGroups(core, molecules, chemistry, timeoutMillis)`, returns one
`Substructure.Decomposition` per molecule. Check `search()` before using
`core()` (the matched atoms and query-bond images; `null` without a match) or
`rgroups()` (attached groups keyed `R1`, `R2` and so on).

## Graph utilities

- `toCanonicalSmiles()`, `toCanonicalSmiles(atomMaps)`: canonical SMILES with
  the supplied stereochemistry; the overload adds positive atom-map labels,
  keyed by atom index, without changing the ranking or the graph. `/` and `\`
  go on single bonds wherever a configured double bond has one at an end; an
  end with only aromatic ring bonds carries them on an aromatic bond, as RDKit
  writes it, which strict OpenSMILES readers take as a single bond; the CDK
  2.13 `SmilesParser` does not keep that configuration.
  For a graph of more than 200 atoms, or when canonical labelling reaches
  its search limit, the output can depend on the input atom order;
  `automorphismGeneratorsTruncated()` is then true.
- `computeRings()`, `numRings()`, `MolGraph.computeSSSR(graph)`,
  `MolGraph.layoutSSSR(graph)`: the smallest set of smallest rings, also
  ordered for layout. `computeRelevantCycles()`, `computeURFs()`: relevant
  cycles and unique ring families. `ringsTruncated()`: whether the ring search
  stopped before finding a complete set of rings. `relevantCyclesTruncated()`:
  whether full relevant-cycle enumeration was incomplete. `ensureRingCounts()`,
  `getRingCounts()`: per-atom ring counts, computed on first use;
  they can depend on the input atom order.
- `computeRelevantCycleFamilies()`: the relevant cycles grouped into families,
  as `RelevantCycleFamilies` with `families()`, `relevantCycleCount()` and
  `truncated()`. Each `RelevantCycleFamily` gives the ring `size()`, the
  `root()` (the highest-numbered atom of every member), `first()`, `second()`
  and `opposite()`, which locate the bond or atom opposite the root, one member
  as `prototype()` and the exact `memberCount()`, so the relevant cycles can be
  counted without listing them; counts are `BigInteger` values.
  `truncated()` reports an incomplete enumeration, whose families are relevant
  but may not be all of them. `relevantCycleFamilyOf(families, ring)` returns
  the index of the family that has the ring as a member, or -1, and
  `relevantCycleFamilyMembers(family, limit)` returns up to `limit` members,
  starting with the prototype. When the enumeration is complete, the relevant
  cycles and their number do not depend on the input atom order; the grouping
  into families can.
- `ringSystemOf(atom)`, `ringSystemCount()`, `ringSystemSig(id)`: ring systems
  and their signatures.
- `MolGraph.sanitizeSmiles(smiles)`: balances parentheses and removes unclosed
  ring openings before CDK parsing.
- `getCanonicalLabeling()`, `getCanonicalHash()`, `getOrbit()`: canonical
  labelling, graph hash and symmetry classes, each a copy (`getOrbits()` is
  deprecated). `getAutomorphismGenerators()`,
  `automorphismGeneratorsTruncated()`: symmetry generators, and whether they
  are incomplete.
- `MolGraph.matchTemplate`, `MolGraph.forceDirectedLayout`,
  `MolGraph.stressMajorisation`, `MolGraph.reduceCrossings`: two-dimensional
  layout with `MolGraph.Point2D` coordinates.
- `MolGraph.translateToAtomIds(mapping, query, target)`: a mapping in the atom
  IDs set with `Builder.atomIds`.

## Command-line interface

A source build of the module `com.bioinceptionlabs:smsd-cli` writes the
self-contained `smsd-cli-8.1.0-jar-with-dependencies.jar`; the
[Java module README](../java/README.md#command-line-tool) shows how to obtain
the 8.0.3 JAR or depend on `SMSDcli`. From the repository root, after a
source build,
print the version and help, test whether benzene is a substructure of phenol,
and find the MCS of N-methylacetamide and propanamide (`mcs_size` is 4):

```bash
java -jar java/cli/target/smsd-cli-8.1.0-jar-with-dependencies.jar --version
java -jar java/cli/target/smsd-cli-8.1.0-jar-with-dependencies.jar --help
java -jar java/cli/target/smsd-cli-8.1.0-jar-with-dependencies.jar \
  --Q SMI --q "c1ccccc1" --T SMI --t "c1ccc(O)cc1" --json -
java -jar java/cli/target/smsd-cli-8.1.0-jar-with-dependencies.jar \
  --Q SMI --q "CC(=O)NC" --T SMI --t "CCC(=O)N" --mode mcs --json -
```

MCS mode defaults to a connected MCIS; request MCES and fragments explicitly,
here for the two ethyl groups of diethyl ether and diethylamine:

```bash
java -jar java/cli/target/smsd-cli-8.1.0-jar-with-dependencies.jar \
  --Q SMI --q "CCOCC" --T SMI --t "CCNCC" \
  --mode mcs --objective MCES --fragments --timeout 1000 --json -
```

SMARTS queries use the type `SIG` in substructure mode; here the first search
reports `"exists": true` (exit status 0) and the second `"exists": false` (1):

```bash
java -jar java/cli/target/smsd-cli-8.1.0-jar-with-dependencies.jar \
  --Q SIG --q "C=O" --T SMI --t "CC(=O)O" --json -
java -jar java/cli/target/smsd-cli-8.1.0-jar-with-dependencies.jar \
  --Q SIG --q "C=O" --T SMI --t "CCO" --json -
```

Options:

- `--Q`, `--q`: query type and value: `SMI` or `SIG` with the string, or
  `MOL`, `ML2` (MOL2), `CML` or `PDB` with a file path.
- `--T`, `--t`: target type and value: `SMI`, `MOL`, `ML2`, `CML`, `PDB`, or
  `SDF` for a batch.
- `--mode`: `sub` (default) or `mcs`, with `--objective` `MCIS` (default) or
  `MCES` and `--fragments` for disconnected fragments.
- `--induced`: require an induced embedding. `-m`, `--mappings`: enumerate up
  to 10,000 embeddings.
- `--timeout`: milliseconds per pair (default 10,000). `--threads`: SDF batch
  threads (default 0, every processor).
- `-a`, `--add-hydrogens`; `-r`, `--remove-hydrogens`: add or remove explicit
  hydrogens first.
- `--json`, `--json-pretty`: (pretty-printed) JSON to a file, or standard
  output with `-`. `-o`, `--output`: the MCS fragment as SMILES, likewise.
- `--map-out`, `--map-format`, `--map-pretty`: mapping export, below.
- `--benchmark`: time substructure and MCS search on one pair and print the
  averages.

Each molecule is standardised (atom typing, implicit hydrogens, rings and
Daylight aromaticity) first. CML and PDB files must hold one molecule or
model, empty inputs are rejected, and an SDF batch may hold up to 100,000
molecules.

JSON output (UTF-8) includes `status`, `optimal`, `complete`, `score`,
`nodes_visited`, `elapsed_ms` and the first match's `bonds` and `fragments`;
substructure mode adds `exists` (with `--mappings`, `mappings` instead), MCS
mode `mcs_size`, `pairs` and, with a match, `mcs_smiles`. An SDF batch has
one row per target under `results`.

`--map-out` writes the mappings of a single-pair search (in substructure mode,
with `--mappings`) to a file or standard output. `--map-format` is `json`
(default), `smiles` (query and target with matching atom maps and
stereochemistry), `smarts` or `mdl` (in MCS mode, the common atoms and bonds of
the first match as a MOL block).

Exit status 0 means a completed search with a match (SDF: at least one match,
every search completed), 1 a completed search without one (SDF: no match,
every search completed), and 2 invalid input, an unknown mode, an error or an
incomplete search (SDF: any target failed or incomplete).

## API reference

Generate the Javadoc from the repository root, then open
`docs/api/index.html`:

```bash
mvn -f java/pom.xml -pl core javadoc:javadoc
```

## Further reading

- [Installation guide](HOWTO-INSTALL.md), including
  [Java CLI installers](HOWTO-INSTALL.md#portable-java-cli).
- [Python guide](PYTHON.md), [C++ guide](CPP.md),
  [benchmark guide](../benchmarks/README.md) and [changelog](../CHANGELOG.md).

SMSD Pro is distributed under the [Apache License 2.0](../LICENSE). Retain
[NOTICE](../NOTICE) and the bundled dependency notices when redistributing it.
