# SMSD Pro C++ guide

A header-only C++17 library in namespace `smsd` for substructure search,
maximum common substructure (MCS) search and related chemistry. Atom indices
are zero-based. Each example is a complete program; build it as in the [C++
quick start](../cpp/README.md#quick-start). See also the [installation
guide](HOWTO-INSTALL.md), the [Java](JAVA.md) and [Python](PYTHON.md) guides
and the [changelog](../CHANGELOG.md).

## Headers

Include `smsd/smsd.hpp`, which includes every header below except
`depictor.hpp`, `bond_energies.hpp`, `color_palette.hpp` and
`rdkit_adapter.hpp`, or only the headers you need: each compiles on its own.
Headers under `smsd/detail/` and namespaces ending in `detail` are
implementation code.

| Area | Headers under `smsd/` |
| --- | --- |
| Search | `mcs.hpp` (`MCS`, mapping utilities); `substructure.hpp` (`Substructure`, `RGroupResult`); `search_types.hpp` (`SearchOptions`, `MCSOptions`, `SubstructureOptions`, `Match`, `SearchResult`, `SearchStatus`, `CancellationToken`); `chem_options.hpp` (`ChemOptions`) |
| Molecules | `mol_graph.hpp` (`MolGraph`, `MolGraph::Builder`); `smiles_parser.hpp` (`parseSMILES`, `toSMILES`, `toSMARTS`); `smarts_parser.hpp` (`parseSMARTS`, `SmartsPredicateRegistry`); `parser_utils.hpp` (input limits); `mol_reader.hpp` (MOL and SDF) |
| Fingerprints and screening | `batch.hpp` (fingerprints, similarity, `TargetCorpus`) and `fingerprint_options.hpp` (`FingerprintOptions`), in `smsd::batch`; `similarity_screen.hpp` (`similarityUpperBound`, `screenTargets`); `gpu.hpp` (Metal or CUDA screening with a CPU fallback, in `smsd::gpu`) |
| Rings, stereo and drawing | `ring_finder.hpp` (`computeSSSR`, `computeRelevantCycles`, `computeURFs`, `computeRings`, and `ringsTruncated`, which says whether the ring search stopped at its work bound); `cip.hpp` (`smsd::cip`); `layout.hpp` (2D and 3D coordinates, `transform`); `depict.hpp` (SVG); `depictor.hpp` (2D coordinates constrained to fixed atoms or a reference) |
| Data and utilities | `hungarian.hpp` (`optimalAssign`); `periodic_table.hpp` (`getElement`, `atomicWeight`, `paulingElectronegativity` and more); `scaffold_library.hpp` (reference ring systems, scaffolds, pharmacophore and functional-group SMARTS, drug-likeness limits); `bond_energies.hpp` (`bondDissociationEnergy`, `mappingBondChangeCost`, in kJ/mol); `color_palette.hpp` (`mapNumberToPaletteColor`, a fixed 20-colour cycle for atom-map numbers); `rdkit_adapter.hpp` ([optional](#optional-native-adapter)) |

## Prepared queries

`MCS::find` and `Substructure::find` take a `MolGraph` or SMILES;
`Substructure::smarts` takes SMARTS text or a `SmartsQuery`. A prepared query
owns copies of its query and options and may be shared between threads.
`withOptions(options)` and `withTimeout(milliseconds)` return new queries;
`MCS::withOptions` throws `std::invalid_argument` for invalid fragment limits
or atom weights. `match` and `matchAll` take a `MolGraph` or SMILES target,
which must not change during a search.

- `match(target)`: one best common subgraph (`MCS`) or the first embedding
  (`Substructure`).
- `matchAll(target, limit = 1000)`: up to `limit` optimal mappings or
  embeddings, symmetric ones included; zero throws `std::invalid_argument`.
- `matchMany(targets, numThreads = 0)`: one `match` result per target, in
  order, each with its own budget; `0` means the hardware concurrency, `1` the
  calling thread, and a negative count throws `std::invalid_argument`. The
  exception of the lowest failing index is rethrown.
- `resultsTimeMillis(target)`: when, in milliseconds from the search start,
  matches found within the timeout are assembled (the effective timeout plus
  one eighth).
- `cancellation()`: the token of the options, or null.

```cpp
#include <cassert>
#include "smsd/smsd.hpp"

int main() {
    const auto query = smsd::parseSMILES("c1ccccc1");
    const auto target = smsd::parseSMILES("c1ccc(O)cc1");
    const auto result = smsd::MCS::find(query).withTimeout(1000).match(target);
    assert(result.status == smsd::SearchStatus::Completed && result.optimal && result.complete);
    assert(result.matches.size() == 1 && result.matches.front().atomMapping.size() == 6);
    assert(smsd::MCS::validateMapping(query, target, result.matches.front()).empty());
    assert(!smsd::Substructure::find(query).match(target).matches.empty());
    assert(!smsd::Substructure::smarts("[OH]").match(target).matches.empty());
    const auto all = smsd::Substructure::find("CO").matchAll("OCO", 10);
    assert(all.complete && all.matches.size() == 2);
}
```

## SMARTS

SMARTS follows Daylight semantics, as in Java: `[H]`, `[2H]` and `[H+]` are
hydrogen atoms, `[CH]` and `[H1]` count hydrogens, `r<n>` is the smallest ring
size, `(C.C)` requires one target component and `(C).(C)` different ones. `/`
and `\` are single, non-aromatic bonds; a double bond between them must have
the written configuration unless one is `/?` or `\?`, and an unpaired one
matches either way.

`parseSMARTS(text, maxRecursionDepth = 20)` throws `std::invalid_argument` for
malformed SMARTS and for input over the `parser_utils.hpp` limits
`MAX_SMARTS_LENGTH`, `MAX_ATOMS` (also applied by `parseSMILES`),
`MAX_SMARTS_EXPRESSION_NODES` and `MAX_SMARTS_NODES_PER_EXPRESSION`.
A `SmartsPredicateRegistry` expands `$name` tokens into recursive SMARTS with
`expandPredicates(text)`; the static `builtinRegistry()` supplies ten, such as
`isHalogen`, `isHBDonor` and `isHBAcceptor`, and
`registerPredicate(name, smarts)` adds more.

## Results

`SearchResult` holds `matches`, `status` (`Completed`, `TimedOut`, `Cancelled`,
or `LimitReached` for a result, node or other limit), `optimal`, `complete`,
`score`, `nodesVisited` and `elapsedMillis`. A stopped search may return an
unproved mapping, and a capped enumeration may be optimal yet incomplete, so
check both flags, and check `matches` before `front()`. `score` counts mapped
atoms for the maximum common induced subgraph (MCIS), selected bonds for the
maximum common edge subgraph (MCES), the summed weight for weighted MCIS, and
query atoms for a `Substructure` match (zero without one).

A `Match` holds `atomMapping` (a `std::map<int, int>` from query to target
atom), `bondMapping` (a `std::array<int, 4>` per selected bond: two query
atoms, then two target atoms; keep it with MCES results, as atom pairs do not
show which bonds were left out) and `fragments` (the components of the selected
bonds, as query atoms). Explicit hydrogen atoms are vertices; bracket hydrogen
counts are not. `Match::toSmiles(query)` writes the selection with valences
capped at cut bonds; `Match::toGraph(query)` returns it, compactly renumbered,
with its original chemical context.

## Objectives and fragments

MCIS, the default, maximises mapped atoms, then common bonds, and preserves
edges and non-edges between mapped atoms. MCES maximises selected common bonds,
then their atoms; isolated atoms do not count. Connectivity and fragment limits
apply to the selected bonds.

`MCS::Options` fields, with defaults:

- `objective` (`Objective::MCIS`): `MCIS` or `MCES`.
- `connectivity` (`Connectivity::Connected`): `Fragments` allows a disconnected
  result.
- `minFragmentAtoms` (`1`) and `maxFragments` (`INT_MAX`): positive limits on
  fragment size and count.
- `excludedTargetAtoms` (empty): target atoms the search may not use.
- `atomWeights` (empty): MCIS only, one finite, signed weight per query atom.

With `atomWeights`, MCIS maximises the summed weight, then atoms, then bonds,
comparing whole thousandths truncated towards zero. Weights with MCES, or whose
positive or negative sum in thousandths exceeds the 32-bit range, throw
`std::invalid_argument`.

```cpp
#include <cassert>
#include "smsd/smsd.hpp"

int main() {
    const auto query = smsd::parseSMILES("C1CC1");
    const auto target = smsd::parseSMILES("CCC");
    const auto induced = smsd::MCS::find(query).match(target);
    assert(induced.matches.front().atomMapping.size() == 2);

    smsd::MCS::Options options;
    options.objective = smsd::MCS::Objective::MCES;
    const auto edges = smsd::MCS::find(query).withOptions(options).match(target);
    const auto& match = edges.matches.front();
    assert(match.atomMapping.size() == 3 && match.bondMapping.size() == 2);
    const auto output = smsd::parseSMILES(match.toSmiles(query));
    std::size_t degreeSum = 0;
    for (const auto& neighbors : output.neighbors) degreeSum += neighbors.size();
    assert(output.n == 3 && degreeSum == 4);
    assert(smsd::MCS::validateMapping(query, target, match, options).empty());

    options.objective = smsd::MCS::Objective::MCIS;
    options.connectivity = smsd::MCS::Connectivity::Fragments;
    const auto fragments = smsd::MCS::find("CO.CN").withOptions(options).match("CO.CN");
    assert(fragments.matches.front().fragments.size() == 2);
}
```

## Search options

`MCS::Options` and `Substructure::Options` share the `SearchOptions` fields:

- `chemistry` (`ChemOptions{}`): atom and bond compatibility rules.
- `timeoutMs` (`-1`): milliseconds; negative selects the default,
  `min(30000, 500 + 2 × query atoms × target atoms)` for `MCS` and 10,000 for
  `Substructure`, and zero stops at once.
- `nodeLimit` (`0`): most search nodes; zero is unlimited.
- `cancellation` (null): a `std::shared_ptr<CancellationToken>`; `cancel()`
  stops the search.
- `induced` (`false`, `Substructure::Options` only): also preserve non-edges.

The budget starts once the query is built and string inputs are parsed, and one
step can overrun it.

```cpp
#include <cassert>
#include <memory>
#include <vector>
#include "smsd/smsd.hpp"

int main() {
    const auto query = smsd::Substructure::smarts("[#6]-[#8]");
    const std::vector<smsd::MolGraph> targets = {
        smsd::parseSMILES("CO"), smsd::parseSMILES("CC")
    };
    const auto results = query.withTimeout(1000).matchMany(targets, 2);
    assert(results.size() == 2);
    assert(results[0].complete && results[0].matches.size() == 1);
    assert(results[1].complete && results[1].matches.empty());

    auto token = std::make_shared<smsd::CancellationToken>();
    smsd::Substructure::Options options;
    options.cancellation = token;
    token->cancel();
    const auto stopped = query.withOptions(options).match(targets.front());
    assert(stopped.status == smsd::SearchStatus::Cancelled);
    assert(stopped.matches.empty() && !stopped.optimal && !stopped.complete);

    const auto immediate = smsd::MCS::find("CC").withTimeout(0).match("CCC");
    assert(immediate.status == smsd::SearchStatus::TimedOut);
    assert(immediate.matches.empty() && !immediate.optimal && !immediate.complete);
}
```

## Chemistry options

`ChemOptions` fields, with defaults:

- `matchAtomType` (`true`), `matchFormalCharge` (`false`) and `matchIsotope`
  (`false`): mapped atoms share their element, their formal charge, or, when
  both specify one, their mass number.
- `matchBondOrder` (`BondOrderMode::STRICT`): `LOOSE` and `ANY` ignore
  bond-order differences.
- `aromaticityMode` (`AromaticityMode::FLEXIBLE`): aromatic bonds may match
  single or double bonds; `STRICT` requires equal aromatic flags.
- `ringMatchesRingOnly` (`false`): ring atoms and bonds map only to ring atoms
  and bonds.
- `completeRingsOnly` (`false`): retained ring atoms and bonds lie in fully
  selected relevant rings of their own input; a fusion atom needs only one.
- `ringFusionMode` (`RingFusionMode::IGNORE`): `STRICT` requires mapped ring
  atoms to lie in equally many relevant cycles.
- `useChirality` and `useBondStereo` (`false`): compare tetrahedral and
  double-bond stereo.
- `tautomerAware` (`false`): relax bond orders where tautomer rules apply.
- `pH` (`7.4`) and `solvent` (`Solvent::AQUEOUS`): the tautomer confidence
  context.

The fluent setters are `withBondOrderMode`, `withAromaticityMode`,
`withCompleteRingsOnly`, `withMatchIsotope`, `withRingFusionMode`,
`withTautomerAware`, `withPH` (clamped to 0–14; throws `std::invalid_argument`
if not finite) and `withSolvent`. `ChemOptions::profile(name)`, ignoring case,
returns `default`, `strict` (strict aromaticity), `pharma` (charge, ring-only
and complete-ring matching), `compat-fmcs`, `compat-substruct` or `reaction`
(loose bond orders), and throws `std::invalid_argument` otherwise;
`ChemOptions::tautomerProfile()` adds tautomer-aware matching with loose bond
orders. Deprecated: `ChemOptions::fmcsProfile()` (use `profile("compat-fmcs")`)
and the unread `aromaticityModel` field (`MolGraph::Builder::build()` and
`perceiveAromaticity()` take the model).

Stereo: tetrahedral centres have four neighbours, or three plus an implicit
hydrogen or, for N, P, S, As, Se, Sb and Te, a lone pair; the SMILES parser
drops other tags, and `MolGraph::Builder` tags on other atoms are compared as
given. A Builder tag (1 for `@`, 2 for `@@`) winds the stored neighbours, the
implicit hydrogen or lone pair coming second if the first neighbour's index is
below the atom's and first otherwise, as in Java. `Substructure` requires
tagged query stereo to recur in the target, and untagged query stereo matches
anything; `MCS` compares the stereo the mapping fixes, treating unspecified
target stereo as a wildcard.

## Mapping utilities

- `MCS::validateMapping(query, target, match, options = {})`: a message per
  rule a `Match` or atom map violates; empty when valid.
- `MCS::mappingMatch(query, target, atoms, options = {})`: a `Match` from an
  atom map, with the compatible bonds between mapped atoms.
- `MCS::isMappingMaximal(query, target, match, options = {})`: whether the
  mapping is valid with no valid strict extension, which is weaker than
  optimal; throws `std::runtime_error` if it cannot decide.
- `MCS::canonicalizeMapping(query, target, mapping)`: one representative under
  the symmetry of both molecules; throws if the symmetry cannot be resolved.
- `smsd::validateTautomerConsistency(query, target, atomMapping)`: whether
  mobile N–H, O–H and S–H protons are conserved in each mapped tautomeric
  group.
- `MCS::computeTautomerConfidence(query, target, match, chemistry)`: a
  heuristic confidence, not a probability, in the selected tautomeric
  relaxations at the solvent and pH of `chemistry`; identical chemistry scores
  1.

Per atom, `MolGraph` exposes `tautomerClass` (from 0, independent of atom order
except among symmetric classes; −1 outside every class), `tautomerWeight` (for
water at pH 7.4 unless `MolGraph::Builder::solvent` and `pH` set another
context), `tautomerSite` (`TAUTOMER_SITE_ACID`, `TAUTOMER_SITE_BASE` or
`TAUTOMER_SITE_NONE`) and `tautomerPKa`; the `MolGraph::TW_*` and `TPKA_*`
constants hold the motif weights and pKa values. Java `MolGraph.tautomerWeight`
is unscaled, so the arrays can differ, but the confidences agree to within
rounding.

## Several molecules

- `MCS::matchConstrained(queries, target, options = {})` runs the queries in
  order against one target, one result each; atoms used by an earlier query are
  unavailable to later ones, and indices refer to the original target.
- `MCS::consensus(molecules, options = {}, threshold = 1.0)` returns a
  `ConsensusResult`: a common graph (`graph`, `fragment`, `sourceIndex`,
  `toSmiles()`) in at least the fraction `threshold` of the inputs, with
  `supportCount` and verified `mappings` keyed by input index; `found` is false
  when none qualifies. `threshold` must lie in (0, 1]; atom weights and
  excluded target atoms are rejected. Consensus is heuristic, shares one budget
  and always reports `optimal == false`.
- `Substructure::decomposeRGroups(core, molecules, options = {})` returns an
  `RGroupResult` per molecule: the matched `core`, substituents in `rgroups`
  (`R1`, `R2`, …) and the `search` result.

## Graphs, SMILES, MOL and SDF

`parseSMILES(text)` throws `std::invalid_argument` for malformed SMILES,
including bracket atoms beyond the OpenSMILES limits;
`parseSMILES(text, options)` with `ParseOptions::lenient = true` recovers where
it can, such as from unclosed rings or unknown elements. `toSMILES(graph)`
writes canonical SMILES, with `/` and `\` on single bonds wherever a configured
double bond has one at an end; an end with only aromatic ring bonds, as at an
exocyclic C=N on an aromatic ring, carries them on an aromatic bond, as RDKit
writes it, which strict OpenSMILES readers take as a single bond.
`toSMARTS(graph, options)` writes canonical SMARTS;
`SmartsWriteOptions` writes aromaticity and charge by default, and ring
membership, hydrogen count and isotope on request.

`MolGraph::Builder` sets atoms (`atomCount`, `atomicNumbers`, `formalCharges`,
`massNumbers`, `hydrogenCounts`, `atomClasses`, `ringFlags`, `aromaticFlags`),
bonds (`setNeighbors`, `setBondOrders`, `bondRingFlags`, `bondAromaticFlags`),
stereo (`tetrahedralChirality`, `doubleBondStereo`, `doubleBondStereoSparse`)
and metadata (`atomIds`, `name`, `programLine`, `comment`, `properties`,
`solvent`, `pH`). `build()` perceives aromaticity, giving supplied aromatic
systems a Kekulé structure and perceiving them again, and `build(false)` keeps
the supplied flags. When perception cannot read every candidate cycle of a
molecule, as with a large bowl-shaped polyarene, the aromatic systems supplied
on its fused ring systems keep their flags.

- `extractSubgraph(graph, atoms)`: the induced subgraph, in the given order.
- `murckoScaffold(graph)`: the Bemis–Murcko scaffold (ring atoms, linkers and
  atoms double- or triple-bonded to them), capped with hydrogen, keeping
  surviving stereo; a molecule without rings is unchanged.
- `translateToAtomIds(mapping, query, target)`: a mapping in external atom IDs.
- `readMolBlock(text)` and `readMolFile(path)` read V2000 or V3000 and throw
  `std::invalid_argument` for a block they cannot read, such as one with a bond
  to a missing atom, a V3000 integer field outside its range, a V3000
  coordinate that is not a decimal number in the range of a double or more than
  `MAX_ATOMS` atoms; V2000 numeric fields are read leniently. Coordinates are
  read by one rule whatever the locale or platform, and one below 1e-307 in
  magnitude reads as zero. `readSDF(path)` reads every record, a
  malformed one giving an empty graph.
- `writeMolBlock(graph)` writes V2000, or V3000 above 999 atoms or bonds;
  `writeMolBlockV3000(graph)` writes V3000; `writeSDFRecord(graph)` writes one
  record ending in `$$$$`, and `writeSDF(graphs, path)` a file.

MOL round trips do not imply support for MDL query features.

```cpp
#include <cassert>
#include <string>
#include "smsd/smsd.hpp"

int main() {
    const auto query = smsd::MolGraph::Builder().atomCount(4)
        .atomicNumbers({7, 6, 6, 8})
        .setNeighbors({{1}, {0, 2}, {1, 3}, {2}})
        .setBondOrders({{1}, {1, 1}, {1, 1}, {1}}).build();
    assert(smsd::MCS::find(query).match("CC(O)CN").matches.front().atomMapping.size() == 4);
    const auto text = smsd::writeMolBlockV3000(query);
    const auto restored = smsd::readMolBlock(text);
    assert(!smsd::Substructure::find(query).match(restored).matches.empty());
    assert(smsd::writeSDFRecord(query).find("$$$$") != std::string::npos);
}
```

## Optional native adapter

Configure with `SMSD_WITH_RDKIT=ON` and include `smsd/rdkit_adapter.hpp`; the
static library is exported as `smsd::smsd_rdkit` and needs an RDKit development
installation and C++20. `smsd::fromRDKit(molecule, removeHs = true)` removes
explicit hydrogens by default, which can change atom indices, and copies
element, charge, isotope, ring, aromaticity and tetrahedral data, but not
hydrogen counts or double-bond stereo. `smsd::fromSmiles(text)` parses with
RDKit and removes hydrogens likewise; `smsd::parseSMILES(text)` is the native
parser.

## Fingerprints, screening and assignment

```cpp
#include <cassert>
#include "smsd/smsd.hpp"

int main() {
    const auto fingerprint = smsd::batch::fingerprintFromSmiles("c1ccc(O)cc1", 2, 2048);
    assert(fingerprint.size() == 32);
    assert(smsd::batch::fingerprintTanimoto(fingerprint, fingerprint) == 1.0);
    const auto assignment = smsd::optimalAssign({{1.0, 3.0}, {4.0, 2.0}});
    assert(assignment.totalCost == 3.0);
}
```

Packed fingerprints are vectors of 64-bit words, count fingerprints vectors of
integer bins; compare only fingerprints of the same family, radius and size.
In `smsd::batch`:

- `fingerprintFromSmiles(smiles, radius = 2, fpSize = 2048)`: ECFP; radius 2 is
  ECFP4.
- `batchFingerprint(molecules, pathLength = 7, fpSize = 1024, numThreads = 0)`:
  path fingerprints over all atoms, explicit hydrogens included, `fpSize`
  rounded up to a multiple of 64; Python `path_fingerprint` uses heavy atoms
  only and differs.
- `fingerprintTanimoto`, `fingerprintOverlap`, `fingerprintDice`,
  `fingerprintCosine` and `fingerprintSoergel` compare packed fingerprints; for
  two empty ones, Tanimoto and overlap give 1 and the others 0, as in Java and
  Python. `countTanimoto`, `countDice` and `countCosine` compare equal-length
  count vectors, giving 0 for empty or unequal ones.
- `fingerprintSubset(query, target)` tests whether every query bit is set in
  the target; `batchFingerprintScreen(query, targets, numThreads = 0)` returns
  the passing target indices.
- `fingerprintQuality(fingerprint, fpSize)`, for packed words or set positions,
  returns a `FingerprintQuality` with `setBits`, `fpSize`, `density`,
  `expectedCollisionRate` and `saturated` (density over
  `FINGERPRINT_SATURATION_DENSITY`, 0.75), as Python `analyze_fp_quality` and
  Java `analyzeFpQuality` report them.
- `toHex`, `fromHex`, `toBinaryString` and `countsToArray` convert formats.
- `TargetCorpus` holds reusable targets: `addTarget`, `addTargets`,
  `prewarm(numThreads = 0)`, `size`, `isPrewarmed`, `targets`,
  `screen(query, threshold, numThreads = 0)`, which returns the indices of the
  targets with ECFP4 Tanimoto similarity of at least `threshold`, and none
  before `prewarm`, and `screenTargets(query, threshold, chemistry,
  numThreads = 0)`, the bound screen below for the held targets, which needs
  no `prewarm`.

`similarityUpperBound(query, target, chemistry)` bounds the atom-count Tanimoto
similarity of the MCS from above;
`screenTargets(query, targets, threshold, chemistry)` returns the indices of
the targets whose bound reaches `threshold`, on the calling thread, or with a
fifth argument on up to `numThreads` threads (0 = the OpenMP thread limit),
giving the same indices. `gpu::batchSimilarityScreenAuto`
computes the bound on CUDA or Metal when built with either and a device is
present, else on the CPU, returning a `ScreenHit{index, score}` per passing
target.
`gpu::isAvailable()` and `gpu::deviceInfo()` report the backend; a non-empty
`SMSD_DISABLE_GPU` other than `0` makes `isAvailable()` false. Batch
fingerprints, screens and `TargetCorpus` run in parallel when compiled with
OpenMP.

`optimalAssign(cost, penalty = 1.0)` solves rectangular minimum-cost assignment
over every row or column of the smaller side, returning the (row, column) pairs
and `totalCost`, which excludes unmatched penalties; ragged or non-finite input
throws `std::invalid_argument`.

### Configurable fingerprints

`FingerprintOptions` gives the same bits in C++, Java and Python for
`computeFingerprint(graph, options)` (packed bits),
`computeFingerprintCounts(graph, options)` (count bins),
`fingerprintFromSmiles` and
`batchFingerprint(molecules, options, numThreads = 0)`.
`FingerprintOptions::ecfp(radius = 2, size = 2048)`,
`fcfp(radius = 2, size = 2048)`, `path(pathLength = 7, size = 2048)` and
`torsion(size = 2048)` reproduce `computeCircularFingerprintECFP`,
`computeCircularFingerprintFCFP`, Python `path_fingerprint` and
`computeTopologicalTorsion`; a default-constructed `FingerprintOptions` is
ECFP4. Turning on another switch adds terms, which can change the classic bits,
so compare only fingerprints built with the same options.

```cpp
#include <cassert>
#include "smsd/smsd.hpp"

int main() {
    using smsd::batch::FingerprintOptions;
    auto options = FingerprintOptions::ecfp(2, 2048);
    options.cipAtoms = true;   // adds CIP R/S terms
    options.ringSize = true;   // adds smallest-ring-size terms
    const auto l = smsd::batch::fingerprintFromSmiles("C[C@H](N)C(=O)O", options);
    const auto d = smsd::batch::fingerprintFromSmiles("C[C@@H](N)C(=O)O", options);
    assert(l != d);
    assert(smsd::batch::fingerprintFromSmiles("C[C@H](N)C(=O)O", FingerprintOptions::ecfp())
           == smsd::batch::fingerprintFromSmiles("C[C@@H](N)C(=O)O", FingerprintOptions::ecfp()));

    auto framework = FingerprintOptions::ecfp();
    framework.scaffold = smsd::batch::ScaffoldMode::MURCKO;
    assert(smsd::batch::fingerprintFromSmiles("Cc1ccccc1", framework)
           == smsd::batch::fingerprintFromSmiles("c1ccccc1", FingerprintOptions::ecfp()));

    const auto keys = smsd::batch::fingerprintFromSmiles("CC(=O)Oc1ccccc1C(=O)O", FingerprintOptions::maccs());
    const auto set = smsd::bitPositions(keys, 166);
    assert(std::find(set.begin(), set.end(), 162) != set.end());   // MACCS key 163, a six-membered ring
}
```

`FingerprintOptions` members:

- `family` (`ECFP`, `FCFP`, `PATH`, `TORSION` or `KEYS`), `radius` (2, or −1
  for the whole component), `pathLength` in bonds (7) and `size` in bits or
  bins (2048).
- Atom terms: `element`, `degree` (explicit neighbours), `valence` (bond order
  sum; pi electrons for `TORSION`), `formalCharge`, `isotope`, `hydrogenCount`,
  `tautomerClass`, `aromaticity`, `ringMembership`, and `ringSize`, the
  smallest SSSR ring, with rings of `ringSizeCap` (8) atoms or more sharing a
  class. Bond terms: `bondOrder`, `aromaticBond` and `ringBond`.
- `cipAtoms` and `cipBonds`: CIP R/S (and r/s) and E/Z from `smsd/cip.hpp`.
  Circular families separate stereoisomers from radius 1; double bonds in rings
  under eight atoms have no E/Z term.
- `fragments`: `WHOLE` (default), `LARGEST` (most heavy atoms) or
  `PER_FRAGMENT` (each atom term carries its component). `scaffold`: `NONE`
  (default) or `MURCKO`, the `murckoScaffold` scaffold, applied first.
- `keys`: the `FingerprintKey(smarts, minimum = 1)` entries of a `KEYS` family;
  a key is set when at least `minimum` distinct atom sets match.
  `FingerprintOptions::maccs()` gives the 166 MACCS keys, key k at bit k − 1,
  listed by `maccsKeys()`; `FingerprintOptions::withKeys(keys)` builds a custom
  family.

Defaults: ECFP sets every atom term except `ringSize`, and `bondOrder`; FCFP
sets `degree`, `ringMembership` and `bondOrder` and always adds the
pharmacophore class; PATH sets `element`; TORSION sets `element`, `degree`,
`valence` and `ringMembership`. Other switches are off. Inconsistent options
and unparsable SMARTS keys throw `std::invalid_argument`; a key count that
cannot be settled within its budget throws `std::runtime_error`. These options
do not configure the MCS screening fingerprint.

## Stereochemistry and depiction

- `cip::assignAll(graph)` returns `CIPDescriptors`, with per-atom R/S labels in
  `rsLabels` and E/Z labels in `ezBonds`; `cip::assignRS` and `cip::assignEZ`
  label one centre or bond, `cip::assignFromSMILES` takes SMILES and
  `cip::to_string` formats a label. IUPAC 2013 Rule 1a applies literally to
  duplicate atoms, so phosphinates and phosphonates may get the opposite
  descriptor from software that ranks the oxo oxygen lower.
- `generateCoords2D(graph, targetBondLength = 1.5)` and
  `generateCoords3D(graph, targetBondLength = 1.5)` generate coordinates;
  `layoutQuality` scores a 2D layout; `transform` translates, rotates, scales,
  mirrors and aligns coordinates.
- `depict(graph, options)` returns SVG;
  `depictWithMapping(graph, atomMapping, options)` highlights mapped atoms;
  `depictPair(query, target, mapping, options)` draws a mapped pair;
  `depictSMILES` and `depictMCS` take SMILES; `DepictOptions` sets the style.
- `depictor.hpp` provides `compute2DCoords(graph, options)`, an overload
  keeping the atoms of a `CoordMap2D` fixed, and
  `generateDepictionMatching2DStructure`, which aligns a molecule to a
  reference by matched atom pairs.

## Deprecated

- `smsd::computeTautomerConfidence(query, target, mapping)`: use
  `MCS::computeTautomerConfidence(query, target, match, chemistry)`; the old
  form selects every mapped bond in both molecules and assumes water at pH 7.4.
- `smsd::implicitHydrogens(z, bondOrderSum, charge)`: use
  `MolGraph::implicitHydrogenCount(atom)`.
- `smsd::ChemOps`: use `MCS::validateMapping`.
- Static `MolGraph::refinePartition`, `comparePerm` and
  `computeCanonicalLabeling`: use `MolGraph::getCanonicalLabeling()`.
- `MolGraph::TW_NITRO_ACI_NITRO`: use `MolGraph::TW_NITRO_ACI`.
- `hungarianSolve` and `HungarianResult`: use `optimalAssign` and
  `AssignmentResult`.
- `runHortonPipeline`: use `computeRings` or `computeRelevantCycles`.

Deprecated with no replacement: `MolGraph::ensurePatternFP()`,
`getNeighborsByDegDesc()`, `getNLF1()`, `getNLF2()`, `getNLF3()`, `nlfCheckOk`,
`TW_SULFOXIDE` and `TW_NITRILE_ISONITRILE`; `popcountAnd`, `popcountN`,
`anyBitAnd` and `popcountAndNot`; and `transform::translate3D`, `rotate3DX`,
`rotate3DY`, `rotate3DZ`, `scale3D`, `center3D` and `rmsd2D`.
