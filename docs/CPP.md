# SMSD Pro C++ guide

The C++ API is a header-only C++17 library in namespace `smsd` for substructure
search, maximum common substructure (MCS) search and related chemistry. Include
`smsd/smsd.hpp` for the main API, or the individual headers listed below. Atom
indices are zero-based throughout. Each example is a complete program; compile
it as shown in the [C++ quick start](../cpp/README.md#quick-start), which also
covers CMake integration and build options.

## Headers

| Header | Provides |
| --- | --- |
| `smsd/smsd.hpp` | Includes every header in this table except `depictor.hpp`, `bond_energies.hpp`, `color_palette.hpp` and `rdkit_adapter.hpp` |
| `smsd/mcs.hpp` | `MCS` prepared queries and mapping utilities |
| `smsd/substructure.hpp` | `Substructure` prepared queries and `RGroupResult` |
| `smsd/search_types.hpp` | `SearchOptions`, `MCSOptions`, `SubstructureOptions`, `Match`, `SearchResult`, `SearchStatus` and `CancellationToken` |
| `smsd/chem_options.hpp` | `ChemOptions` and its named profiles |
| `smsd/mol_graph.hpp` | `MolGraph` and `MolGraph::Builder` |
| `smsd/smiles_parser.hpp` | `parseSMILES`, `toSMILES` and `toSMARTS` |
| `smsd/smarts_parser.hpp` | `parseSMARTS` and `SmartsPredicateRegistry` |
| `smsd/parser_utils.hpp` | Input limits shared by the SMILES and SMARTS parsers |
| `smsd/mol_reader.hpp` | MOL and SDF reading and writing |
| `smsd/batch.hpp` | Fingerprints, similarity metrics and `TargetCorpus`, in `smsd::batch` |
| `smsd/similarity_screen.hpp` | `similarityUpperBound` and `screenTargets` |
| `smsd/gpu.hpp` | Optional Metal or CUDA screening with a CPU fallback, in `smsd::gpu` |
| `smsd/ring_finder.hpp` | `computeSSSR`, `computeRelevantCycles`, `computeURFs` and `computeRings` |
| `smsd/cip.hpp` | CIP descriptors, in `smsd::cip` |
| `smsd/layout.hpp` | 2D and 3D coordinate generation and coordinate transforms |
| `smsd/depict.hpp` | SVG depiction |
| `smsd/depictor.hpp` | 2D coordinates constrained to fixed atoms or to a reference structure |
| `smsd/hungarian.hpp` | `optimalAssign` |
| `smsd/periodic_table.hpp` | Element data such as `getElement`, `atomicWeight` and `paulingElectronegativity` |
| `smsd/scaffold_library.hpp` | Reference tables of ring systems, scaffolds, pharmacophore and functional-group SMARTS, and drug-likeness limits |
| `smsd/bond_energies.hpp` | Bond dissociation energies in kJ/mol, such as `bondDissociationEnergy` and `mappingBondChangeCost` |
| `smsd/color_palette.hpp` | `mapNumberToPaletteColor`, a fixed 20-colour cycle for atom-map numbers |
| `smsd/rdkit_adapter.hpp` | The optional RDKit adapter; see [Optional native adapter](#optional-native-adapter) |

Headers under `smsd/detail/`, and namespaces whose names end in `detail`, hold
implementation code that the headers above include.

Every header directly under `smsd/` compiles on its own, so a program can
include only the headers it uses instead of `smsd/smsd.hpp`. The
`smsd_header_self_containment_tests` CTest entry, built with
`SMSD_BUILD_TESTS=ON`, compiles each one alone as C++17 with `-Wall -Wextra`.

## Prepared queries

`MCS::find` and `Substructure::find` accept a `MolGraph` or a SMILES string;
`Substructure::smarts` accepts SMARTS text or a parsed `SmartsQuery`. A prepared
query owns a copy of its query and options. `withOptions(options)` and
`withTimeout(milliseconds)` return new queries; `MCS::withOptions` validates the
options and throws `std::invalid_argument` for invalid fragment limits or atom
weights.

| Method | Result |
| --- | --- |
| `match(target)` | `MCS`: one best common subgraph. `Substructure`: the first full-query embedding. |
| `matchAll(target, limit = 1000)` | `MCS`: optimal mappings. `Substructure`: full-query embeddings. At most `limit` matches; a limit of zero throws `std::invalid_argument`. |
| `matchMany(targets, numThreads = 0)` | One `match` result per target, in input order |
| `resultsTimeMillis(target)` | The time, in milliseconds from the start of a search of `target`, by which matches found within the timeout are assembled: the effective timeout plus one eighth of it |
| `cancellation()` | The cancellation token of the options, or null |

`match` and `matchAll` accept a `MolGraph` or a SMILES string as the target.

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

Each call keeps its own search state, so one prepared query may be used from
several threads at once. Do not modify a target while it is being searched.
`matchAll` reports symmetry-related embeddings as distinct matches, and a
bounded call can stop before it proves the optimum or completes the
enumeration.

`matchMany` gives each target its own timeout and node limit and runs on
standard C++ threads, not OpenMP. `numThreads = 0` uses the hardware
concurrency, `1` runs on the calling thread and a negative value throws
`std::invalid_argument`. If a target throws, `matchMany` rethrows the exception
of the lowest failing input index.

## SMARTS

SMARTS follows Daylight semantics, as in the Java library. `[H]`, `[2H]` and
`[H+]` are hydrogen atoms, while `[CH]` and `[H1]` count attached hydrogens.
`r<n>` is the size of the smallest ring through the atom. `(C.C)` requires
both atoms to lie in one target component and `(C).(C)` requires them to lie in
different components. A directional bond, `/` or `\`, is a single, non-aromatic
bond; a double bond between directional bonds must have the written
configuration unless one of them is written `/?` or `\?`, and a directional
bond with no configured double bond beside it matches in either direction.

`parseSMARTS(text, maxRecursionDepth = 20)` throws `std::invalid_argument` for
malformed SMARTS and for patterns beyond the limits in `parser_utils.hpp`:
`MAX_SMARTS_LENGTH` (100,000 characters), `MAX_ATOMS` (10,000 atoms, counting
those of recursive SMARTS), `MAX_SMARTS_EXPRESSION_NODES` (100,000) and
`MAX_SMARTS_NODES_PER_EXPRESSION` (10,000). `MAX_ATOMS` also applies to
`parseSMILES`.

`SmartsPredicateRegistry` expands `$name` tokens into recursive SMARTS before
parsing: call `expandPredicates(text)` and pass the result to
`Substructure::smarts`. `SmartsPredicateRegistry::builtinRegistry()` supplies
ten predicates, including `isHalogen`, `isHBDonor` and `isHBAcceptor`, and
`registerPredicate(name, smarts)` adds others.

## Results

`SearchResult` contains `matches`, `status`, `optimal`, `complete`, `score`,
`nodesVisited` and `elapsedMillis`.

| `status` | Meaning |
| --- | --- |
| `Completed` | No limit stopped the search. |
| `TimedOut` | The timeout expired. |
| `Cancelled` | The cancellation token was cancelled. |
| `LimitReached` | The result limit, the node limit or another search limit was reached. |

Because `LimitReached` has several causes, inspect `optimal` and `complete`
separately. A stopped call may return a valid mapping without proving the
optimum. A capped enumeration may prove the objective while leaving the
enumeration incomplete. A completed search can return an empty `matches`
vector, so check it before calling `front()`.

`score` is the primary objective: the atom count for an unweighted maximum
common induced subgraph (MCIS), the selected-bond count for a maximum common
edge subgraph (MCES), or the summed weight for a weighted MCIS.
Substructure results score the query atom count when a match is found and zero
otherwise.

`Match` contains:

- `atomMapping`: a `std::map<int, int>` from query atom to target atom;
- `bondMapping`: one `std::array<int, 4>` per selected bond, holding its two
  query atoms followed by the two matching target atoms;
- `fragments`: the connected components of the selected bonds, as query atom
  indices.

Explicit hydrogen atoms are graph vertices; bracket hydrogen counts do not add
vertices. For an MCES result, keep the bond mapping as well as the atom
mapping, because the atom pairs alone do not identify which input bonds were
left out.

`Match::toSmiles(query)` writes the selected atoms and bonds as a molecule, with
the valences at cut bonds capped. `Match::toGraph(query)` returns the selected
graph with its original chemical context, such as ring perception and hydrogen
counts, for further matching; it is not the same as the capped output molecule.
The extracted graph has compact atom indices, while the mapping and fragment
indices still refer to the original query.

## Objectives and fragments

MCIS is the default objective: it maximises mapped atoms, then common bonds,
and preserves both edges and non-edges between mapped atoms. MCES maximises
selected compatible common bonds, then their incident atoms; isolated
atoms cannot pad the result, and input bonds outside the selection may be left
out. Connectivity and fragment limits apply to the selected common bonds.

| `MCS::Options` field | Default | Effect |
| --- | --- | --- |
| `objective` | `Objective::MCIS` | `MCIS` or `MCES` |
| `connectivity` | `Connectivity::Connected` | `Fragments` allows a disconnected result |
| `minFragmentAtoms` | `1` | Smallest fragment retained; must be positive |
| `maxFragments` | `INT_MAX` | Largest number of fragments; must be positive |
| `excludedTargetAtoms` | empty | Original target atom indices the search may not use |
| `atomWeights` | empty | MCIS only: one finite, signed weight per query atom |

With `atomWeights`, MCIS maximises the summed weight, then atoms, then bonds.
Summed weights are compared in whole thousandths of a unit, truncated towards
zero. If the sum of the positive weights or of the negative weights, counted in
thousandths, falls outside the 32-bit integer range, the options are rejected
with `std::invalid_argument`. MCES rejects atom weights.

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

`MCS::Options` and `Substructure::Options` share the fields of `SearchOptions`.

| Field | Default | Effect |
| --- | --- | --- |
| `chemistry` | `ChemOptions{}` | Atom and bond compatibility rules |
| `timeoutMs` | `-1` | Timeout in milliseconds; a negative value selects the default and zero requests an immediate stop |
| `nodeLimit` | `0` | Largest number of search nodes; zero leaves it unlimited |
| `cancellation` | null | A shared `CancellationToken` |
| `induced` | `false` | `Substructure::Options` only: also require non-edges to be preserved |

`MCS` has no `induced` field; its objective determines how bonds are treated.

The search budget of a call is its timeout together with its node limit. The
default timeout is `min(30000, 500 + 2 × query atoms × target atoms)`
milliseconds for `MCS` and 10,000 milliseconds for `Substructure`. Graph
property preparation, search stages and recursive SMARTS evaluation share the
budget of the call. Constructing a prepared query and parsing string inputs
happen before the budget starts, and an individual step can run past the
deadline. To stop a search cooperatively, set `cancellation` to a
`std::shared_ptr<CancellationToken>` and call `cancel()`; a cancelled token
also stops the assembly of results.

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

| `ChemOptions` field | Default | Effect |
| --- | --- | --- |
| `matchAtomType` | `true` | Mapped atoms must have the same element |
| `matchFormalCharge` | `false` | Mapped atoms must have the same formal charge |
| `matchIsotope` | `false` | Mapped atoms must have the same mass number when both specify one |
| `matchBondOrder` | `BondOrderMode::STRICT` | `LOOSE` and `ANY` ignore bond-order differences once the other bond rules pass |
| `aromaticityMode` | `AromaticityMode::FLEXIBLE` | `FLEXIBLE` lets aromatic bonds match single or double bonds; `STRICT` also requires matching aromatic flags |
| `ringMatchesRingOnly` | `false` | Ring atoms and bonds map only to ring atoms and bonds |
| `completeRingsOnly` | `false` | Retained ring atoms and bonds must belong to complete rings, as described below |
| `ringFusionMode` | `RingFusionMode::IGNORE` | `STRICT` also requires mapped ring atoms to have the same ring count, the number of relevant cycles that the ring search finds through the atom |
| `useChirality` | `false` | Compare mapped tetrahedral stereo |
| `useBondStereo` | `false` | Compare mapped double-bond stereo |
| `tautomerAware` | `false` | Relax bond orders where the tautomer rules apply, without changing `matchBondOrder` |
| `pH` | `7.4` | Context for tautomer confidence; `withPH` clamps it to the range 0 to 14 |
| `solvent` | `Solvent::AQUEOUS` | Context for tautomer confidence |

The fluent setters are `withBondOrderMode`, `withAromaticityMode`,
`withCompleteRingsOnly`, `withMatchIsotope`, `withRingFusionMode`,
`withTautomerAware`, `withPH` and `withSolvent`; `withPH` throws
`std::invalid_argument` for a non-finite value.

`ChemOptions::profile(name)` returns a named profile, ignoring case:
`default`, `strict` (strict aromaticity), `pharma` (charge, ring-only and
complete-ring matching), `compat-fmcs`, `compat-substruct` and `reaction`
(loose bond orders). An unknown name throws `std::invalid_argument`.
`ChemOptions::tautomerProfile()` enables tautomer-aware matching with loose bond
orders. `ChemOptions::fmcsProfile()` is deprecated; it returns
`profile("compat-fmcs")`. The `aromaticityModel` field is deprecated as well:
no search reads it, and `MolGraph::Builder::build()` and `perceiveAromaticity()`
take the model directly.

With `completeRingsOnly`, both inputs are checked: a retained ring atom or bond
must belong to a fully selected relevant ring in its own input, and a fusion
atom may be retained through one complete ring.

With `useChirality`, a tagged atom is a tetrahedral centre when it has four
neighbours and no implicit hydrogen, or three neighbours and either one implicit
hydrogen or, for N, P, S, As, Se, Sb and Te, a lone pair. Tags on other atoms
are compared as they are, and the SMILES parser drops them. A tag (1 for `@`,
2 for `@@`) winds the neighbours in their stored order, with the implicit
hydrogen or lone pair second when the first neighbour has a lower index than the
atom and first otherwise, so `MolGraph::Builder` tags have the same meaning in
C++ and Java.

A `Substructure` search embeds the query, so with `useChirality` a tagged query
centre needs a tagged target centre of the same winding over the mapped
ligands, and with `useBondStereo` a configured query double bond needs a
configured target bond of the same relation, as a SMARTS chirality primitive or
pair of directional bonds does; an untagged or unconfigured query atom or bond
matches any target. An `MCS` search compares the stereo that the selected
mapping determines and keeps an unspecified target atom or bond as a wildcard.

## Mapping utilities

| Function | Purpose |
| --- | --- |
| `MCS::validateMapping(query, target, match, options = {})` | Returns the violated rules of a `Match` or an atom map as messages; an empty vector means the mapping is valid |
| `MCS::mappingMatch(query, target, atoms, options = {})` | Builds a `Match` from an atom map, selecting the compatible bonds between mapped atoms |
| `MCS::isMappingMaximal(query, target, match, options = {})` | Whether the mapping is valid and has no valid strict extension; throws `std::runtime_error` if that search stops before a proof |
| `MCS::canonicalizeMapping(query, target, mapping)` | One representative atom mapping under the symmetry of both molecules; throws if the symmetry calculation cannot complete within its internal limits |
| `MCS::computeTautomerConfidence(query, target, match, chemistry)` | A heuristic confidence for the selected tautomeric relaxations |
| `smsd::validateTautomerConsistency(query, target, atomMapping)` | Whether a mapping conserves the mobile N–H, O–H and S–H protons within each mapped tautomeric group |

Maximality is weaker than a proof of optimality. Tautomer confidence uses the
pH and solvent of `chemistry` and changes neither the molecules nor the search
objective; identical selected chemistry scores 1, and the value is not a
probability.

### Tautomer weights, solvent and pH

The confidence is the geometric mean, over the relaxed atoms with a tautomer
class in either molecule, of the smaller of the two atom weights. Each atom of
a class takes the weight (`MolGraph::TW_*`) of the first motif in the table
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
`MolGraph::tautomerClass` numbers the classes from 0 in increasing order of a
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

The solvent shift is applied first. The pH then scales the weight of every
atom with a site by the neutral fraction of its class. With pKa_a the lowest
acid pKa and pKa_b the highest base pKa (that of BH+) among the sites of the
class:

```text
f0 = 1 / (1 + 10^(pH - pKa_a) + 10^(pKa_b - pH))
weight = min(0.99, max(0.1, solvent-corrected weight × f0))
```

A class with acid sites only drops the base term, which leaves
1 / (1 + 10^(pH − pKa_a)), the fraction of the neutral acid HA; a class with
base sites only drops the acid term, which leaves 1 / (1 + 10^(pKa_b − pH)),
the fraction of the neutral base B. The library assigns no pKa to the motifs
marked — in the table, and the pH does not change their weights.
`MolGraph::tautomerSite` holds each atom's site (`TAUTOMER_SITE_ACID`,
`TAUTOMER_SITE_BASE` or `TAUTOMER_SITE_NONE`) and `tautomerPKa` its pKa.

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
at every pH.

A graph stores `tautomerWeight` in its own context, water at pH 7.4 unless
`MolGraph::Builder::solvent` and `pH` set another;
`MCS::computeTautomerConfidence` recomputes the weights from the baseline in
the context of `chemistry`, so a graph's own context is never applied twice.
Java `MolGraph.tautomerWeight` holds the motif weight in every context, so
the C++ and Java arrays differ wherever a context scales a weight, while the
confidences agree within 1e-12.

## Several molecules

`MCS::matchConstrained(queries, target, options = {})` runs the queries in
order against one target and returns one result per query. Target atoms
matched by an earlier query are unavailable to later ones; indices refer to the
original target, which is not re-perceived.

`MCS::consensus(molecules, options = {}, threshold = 1.0)` returns a
`ConsensusResult`: a common graph (`graph`, `fragment`, `sourceIndex` and
`toSmiles()`) supported by at least the fraction `threshold` of the inputs,
with `supportCount` and verified `mappings` from the compact common graph to
each supporting input, keyed by input index; `found` is false when no common
graph reaches the threshold. The threshold must lie in (0, 1];
atom weights and excluded target atoms are rejected. The whole call shares one
search budget. Consensus is heuristic and always reports `optimal == false`.

`Substructure::decomposeRGroups(core, molecules, options = {})` returns one
`RGroupResult` per molecule with the matched `core`, the substituents in
`rgroups` (keyed `R1`, `R2` and so on) and the underlying `search` result.

## Graphs, SMILES, MOL and SDF

`parseSMILES(text)` throws `std::invalid_argument` for malformed SMILES,
including a bracket atom outside the OpenSMILES limits: a hydrogen count has one
digit, a charge two digits or repeated signs within -15 to +15, and an isotope
three digits (`[CH9]`, `[C+15]` and `[999C]` are accepted; `[CH10]`, `[C+16]`
and `[1000C]` are rejected). `parseSMILES(text, options)` with
`ParseOptions::lenient = true` recovers where it can, for example by ignoring
unclosed rings and skipping unknown elements, but keeps these limits.
`toSMILES(graph)` writes canonical SMILES. `toSMARTS(graph, options)` writes
canonical SMARTS; `SmartsWriteOptions` includes aromaticity and charge by
default, and ring membership, hydrogen count and isotope on request.

`MolGraph::Builder` sets atoms (`atomCount`, `atomicNumbers`, `formalCharges`,
`massNumbers`, `hydrogenCounts`, `atomClasses`, `ringFlags`, `aromaticFlags`),
bonds (`setNeighbors`, `setBondOrders`, `bondRingFlags`, `bondAromaticFlags`),
stereo (`tetrahedralChirality`, `doubleBondStereo`, `doubleBondStereoSparse`)
and metadata (`atomIds`, `name`, `programLine`, `comment`, `properties`,
`solvent`, `pH`). `build()` perceives aromaticity unless called as
`build(false)`. `extractSubgraph(graph, atoms)` returns the induced subgraph in
the given atom order, `murckoScaffold(graph)` keeps ring atoms, the linkers
between them and every atom joined to either by a double or triple bond, so
that exocyclic C=O, C=N and C=S groups stay as in the Bemis and Murcko scaffold
(a molecule without rings is returned unchanged), and
`translateToAtomIds(mapping, query, target)` converts a mapping to external
atom IDs. A scaffold replaces each stripped substituent by hydrogen: an atom
that lost neighbours takes the valence model's hydrogen count for its remaining
bonds, and an aromatic N or P its count plus one per stripped neighbour, as the
Java scaffold and `Match::toSmiles()` do. It keeps the ring and aromaticity
perception of the input, the configuration of each stereocentre whose
neighbours all stay and that of each double bond whose ends each keep a ligand.

| Function | Purpose |
| --- | --- |
| `readMolBlock(text)` | Reads a V2000 or V3000 MOL block |
| `readMolFile(path)` | Reads a MOL file |
| `readSDF(path)` | Reads every record of an SDF file; a malformed record yields an empty graph in its place |
| `writeMolBlock(graph)` | Writes V2000, or V3000 above 999 atoms or bonds |
| `writeMolBlockV3000(graph)` | Writes V3000 |
| `writeSDFRecord(graph)` | Writes one SDF record, ending in `$$$$` |
| `writeSDF(graphs, path)` | Writes an SDF file |

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

An importer that builds graphs with `MolGraph::Builder` should also supply
charge, isotope, hydrogen-count, ring, aromaticity and stereo data. MOL round
trips do not imply support for MDL query features.

## Optional native adapter

Configure with `SMSD_WITH_RDKIT=ON` and include `smsd/rdkit_adapter.hpp`. The
adapter is a compiled static library, `smsd_rdkit`, which the installed CMake
package exports as `smsd::smsd_rdkit`. It requires an RDKit development
installation and C++20, and is separate from the C++17 header-only core.

`smsd::fromRDKit(molecule, removeHs = true)` removes explicit hydrogen atoms by
default, which can change atom indices. It copies element, charge, isotope,
ring, aromaticity and tetrahedral data, but not hydrogen counts or double-bond
stereo. `smsd::fromSmiles(text)` parses with RDKit and removes explicit
hydrogens in the same way; `smsd::parseSMILES(text)` uses the native parser.

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

Packed fingerprints are vectors of 64-bit words; count fingerprints are vectors
of integer bins. In `smsd::batch`:

- `fingerprintFromSmiles(smiles, radius = 2, fpSize = 2048)` returns an ECFP
  fingerprint (radius 2 corresponds to ECFP4).
- `batchFingerprint(molecules, pathLength = 7, fpSize = 1024, numThreads = 0)`
  returns path fingerprints over all atoms, including explicit hydrogen atoms,
  hashing atom labels and bond orders; `fpSize` is rounded up to a multiple of
  64.
- `fingerprintTanimoto`, `fingerprintOverlap`, `fingerprintDice`,
  `fingerprintCosine` and `fingerprintSoergel` compare packed fingerprints, and
  `countTanimoto`, `countDice` and `countCosine` compare count vectors of equal
  length, returning 0 for empty or unequal vectors.
- `fingerprintSubset(query, target)` tests whether every bit of the query is
  set in the target, and `batchFingerprintScreen(query, targets, numThreads = 0)`
  returns the indices of the targets that pass.
- `toHex`, `fromHex`, `toBinaryString` and `countsToArray` convert between
  representations.
- `TargetCorpus` holds reusable targets, with `addTarget`, `addTargets`,
  `prewarm(numThreads = 0)`, `size`, `isPrewarmed`, `targets` and
  `screen(query, threshold, numThreads = 0)`. `screen` returns the indices of
  the targets whose ECFP4 Tanimoto similarity to the query is at least
  `threshold`, and returns no indices until `prewarm` has run.

Two empty fingerprints are identical sets: `fingerprintTanimoto` and
`fingerprintOverlap` return 1 and `fingerprintSoergel` returns 0, as in Java and
Python, while `fingerprintDice` and `fingerprintCosine` return 0. Python
`path_fingerprint` uses a separate element path fingerprint, which hashes
atomic-number paths over heavy atoms only and is not interchangeable with
`batchFingerprint`. Compare fingerprints of the same family, radius and size.
Screening scores and canonical hashes do not prove molecular identity.

### Configurable fingerprints

`smsd/fingerprint_options.hpp` (included by `smsd/smsd.hpp`) adds
`FingerprintOptions`, a specification shared with Java and Python: the same
options give the same bits in all three. `computeFingerprint(graph, options)`
returns packed bits and `computeFingerprintCounts(graph, options)` count bins;
`fingerprintFromSmiles(smiles, options)` and
`batchFingerprint(molecules, options, numThreads = 0)` are overloads of the
fixed functions. The factories `FingerprintOptions::ecfp(radius = 2, size = 2048)`,
`fcfp(radius = 2, size = 2048)`, `path(pathLength = 7, size = 2048)` and
`torsion(size = 2048)` give the classic configuration of each family, which
reproduces the bits of `computeCircularFingerprintECFP`, `computeCircularFingerprintFCFP`,
the Python `path_fingerprint` and `computeTopologicalTorsion` exactly; a
default-constructed `FingerprintOptions` is the classic ECFP4. Every further
switch extends the hashed terms without moving the classic ones.

```cpp
#include <cassert>
#include "smsd/smsd.hpp"

int main() {
    using smsd::batch::FingerprintOptions;
    auto options = FingerprintOptions::ecfp(2, 2048);
    options.cipAtoms = true;   // R/S joins the atom environments after the first iteration
    options.ringSize = true;   // the smallest ring through each atom, as a class up to 8
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

| Member | Meaning | ECFP | FCFP | PATH | TORSION |
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
| `cipAtoms`, `cipBonds` | CIP R/S (and r/s) of tetrahedral centres and E/Z of double bonds, from `smsd/cip.hpp` | off | off | off | off |
| `fragments` | `WHOLE`, `LARGEST` (the component with the most heavy atoms, ties resolved without regard to atom order) or `PER_FRAGMENT` (every atom term carries its component's class) | `WHOLE` | | | |
| `scaffold` | `NONE`, or `MURCKO`: the Bemis-Murcko framework, that is ring atoms, the linkers between rings and the atoms double- or triple-bonded to them, capped with hydrogen where side chains were; a molecule with no ring is used whole | `NONE` | | | |
| `keys` | The keys of the `KEYS` family, one bit each: `FingerprintKey(smarts, minimum = 1)` is set when the molecule holds at least `minimum` distinct matched atom sets | | | | |

FCFP always hashes the pharmacophore feature class first. The stereo terms are
the CIP descriptors of the whole molecule, so a remote substituent that changes
a centre's priorities changes its label. They join the circular hashes after
the first iteration, so radius 0 does not separate stereoisomers and the bits
of a stereoisomer pair differ only within radius − 1 bonds of the centre or
bond; the path family folds R/S into its atom terms, so single-atom paths
already differ, and E/Z into its bond terms, and the torsion family folds both
into the atom types. Double bonds in rings of fewer than eight atoms carry no
E/Z term. The framework keeps exocyclic carbonyls and imines, as the usual
convention does, whereas `murckoScaffold` removes every side chain; a
tetrahedral centre or double bond that stays stereogenic after capping keeps
its configuration. The scaffold is taken before the fragment choice. `FingerprintOptions::maccs()` is the key
family of the 166 MACCS structural keys, key k at bit k − 1, where keys 1, 125
and 166 count isotopes, aromatic rings and fragments; `maccsKeys()` lists the
definitions. `FingerprintOptions::withKeys(keys)` builds a custom key family.
Inconsistent options, and a SMARTS key that does not parse, throw
`std::invalid_argument`; a key search that does not complete within its time
budget, or that reaches 100,000 embeddings before its count is settled, throws
`std::runtime_error` rather than returning bits that depend on the machine. The
MCS screening fingerprint is a separate function that these options do not
configure.

`similarityUpperBound(query, target, chemistry)` returns an upper bound on the
atom-count Tanimoto similarity of the MCS, derived from atom-label frequencies,
and `screenTargets(query, targets, threshold, chemistry)` returns the indices of
targets whose bound is at least `threshold`. `gpu::batchSimilarityScreenAuto`
computes the same bound on CUDA or Metal when the library was built with
either and a device is present, and on the CPU otherwise, returning a
`ScreenHit{index, score}` entry for each target at or above the threshold.
`gpu::isAvailable()` and `gpu::deviceInfo()` report the backend; setting the
environment variable `SMSD_DISABLE_GPU` to a non-empty value other than `0`
makes `isAvailable()` return false. `batchFingerprint`,
`batchFingerprintScreen`, `TargetCorpus` and the CPU path of
`batchSimilarityScreenAuto` run in parallel when the program is compiled with
OpenMP, and sequentially otherwise.

`optimalAssign(cost, penalty = 1.0)` solves the rectangular minimum-cost
assignment problem. It assigns every row or column of the smaller side and
returns the (row, column) pairs with `totalCost`, the sum of their costs,
excluding unmatched penalties. Ragged matrices and non-finite values throw
`std::invalid_argument`.

## Stereochemistry and depiction

- `cip::assignAll(graph)` returns `CIPDescriptors`, with an R/S label per atom
  in `rsLabels` and E/Z labels for stereogenic double bonds in `ezBonds`.
  `cip::assignRS` and `cip::assignEZ` label an individual centre or double
  bond, `cip::assignFromSMILES` accepts SMILES input and `cip::to_string`
  formats a label. Duplicate atoms follow IUPAC 2013 Rule 1a literally: at a
  phosphorus centre the oxo oxygen, which carries a duplicate phosphorus,
  ranks above a hydroxy or alkoxy oxygen, so phosphinates and phosphonates
  may receive the opposite descriptor from toolkits that rank it lower.
- `generateCoords2D(graph, targetBondLength = 1.5)` and
  `generateCoords3D(graph, targetBondLength = 1.5)` generate coordinates;
  `layoutQuality` scores a 2D layout, and the `transform` namespace translates,
  rotates, scales, mirrors and aligns coordinate sets.
- `depict(graph, options)` returns an SVG string.
  `depictWithMapping(graph, atomMapping, options)` highlights the atoms that
  are keys of the mapping, `depictPair(query, target, mapping, options)` draws
  two molecules side by side with their mapping, and `depictSMILES` and
  `depictMCS` accept SMILES input. `DepictOptions` controls the drawing style.
- `depictor.hpp` provides `compute2DCoords(graph, options)`, an overload that
  keeps the atoms of a `CoordMap2D` fixed, and
  `generateDepictionMatching2DStructure`, which aligns a molecule to a
  reference structure through matched atom pairs.

## Deprecated

| Deprecated | Replacement |
| --- | --- |
| `smsd::computeTautomerConfidence(query, target, mapping)` with three arguments | `MCS::computeTautomerConfidence(query, target, match, chemistry)`; the deprecated form selects every mapped bond present in both molecules and assumes water at pH 7.4 |
| `smsd::implicitHydrogens(z, bondOrderSum, charge)` | `MolGraph::implicitHydrogenCount(atom)` gives the count for an atom of a graph |
| `smsd::ChemOps` | `MCS::validateMapping` checks a mapping against the search rules |
| `MolGraph::ensurePatternFP()` and `MolGraph::getNeighborsByDegDesc()` | None; the library does not read these caches |
| `MolGraph::getNLF1()`, `getNLF2()`, `getNLF3()` and `nlfCheckOk` | None; the searches do not read the neighbour-label tables |
| The static `MolGraph::refinePartition`, `comparePerm` and `computeCanonicalLabeling` | `MolGraph::getCanonicalLabeling()` on the graph |
| `MolGraph::TW_NITRO_ACI_NITRO` | `MolGraph::TW_NITRO_ACI` |
| `MolGraph::TW_SULFOXIDE` and `MolGraph::TW_NITRILE_ISONITRILE` | None; no tautomer class uses these weights |
| `hungarianSolve` and `HungarianResult` | `optimalAssign` and `AssignmentResult` |
| `runHortonPipeline` | `computeRings` or `computeRelevantCycles` |
| `popcountAnd`, `popcountN`, `anyBitAnd` and `popcountAndNot` | None; the fused bitset operations are unused |
| `transform::translate3D`, `rotate3DX`, `rotate3DY`, `rotate3DZ`, `scale3D`, `center3D` and `rmsd2D` | None; the library does not use the 3-D transforms |

## Further reading

See the [C++ module README](../cpp/README.md) for building, testing and CMake
integration, the [installation guide](HOWTO-INSTALL.md), the
[Java guide](JAVA.md), the [Python guide](PYTHON.md) and the
[changelog](../CHANGELOG.md).
