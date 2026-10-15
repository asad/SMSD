# SMSD C++ guide

The C++17 header-only source API provides `MCS` and `Substructure` queries.
It replaces the search entry points in release 7.2.2. Include `smsd/smsd.hpp`
for the main API, or include `smsd/mcs.hpp` and `smsd/substructure.hpp` directly.

Save each C++ example as `example.cpp` and compile from the repository root:

```bash
c++ -std=c++17 -O2 -Icpp/include example.cpp -o example
./example
```

An installed CMake package provides the `smsd::smsd` target:

```cmake
find_package(smsd 7.2 CONFIG REQUIRED)
target_link_libraries(my_program PRIVATE smsd::smsd)
```

## Molecule and SMARTS search

```cpp
#include <cassert>
#include "smsd/smsd.hpp"

int main() {
    const auto query = smsd::parseSMILES("c1ccccc1");
    const auto target = smsd::parseSMILES("c1ccc(O)cc1");
    const auto result = smsd::MCS::find(query).withTimeout(1000).match(target);
    assert(result.optimal && result.matches.front().atomMapping.size() == 6);
    assert(smsd::MCS::validateMapping(query, target, result.matches.front()).empty());
    assert(!smsd::Substructure::find(query).match(target).matches.empty());
    assert(!smsd::Substructure::smarts("[OH]").match(target).matches.empty());
    const auto all = smsd::Substructure::find("CO").matchAll("OCO", 10);
    assert(all.complete && all.matches.size() == 2);
}
```

Prepared queries own their query and options. Fluent methods return new
queries. Each call owns its search state, so a prepared query may be reused
concurrently. Do not mutate a target while it is being searched.

`match` requests one best match; `matchAll(target, limit)` enumerates optimal
atom and selected-bond embeddings up to a positive limit. `matchMany(targets,
numThreads)` returns one result per target in input order. Each pair has its
own timeout. Zero workers uses the available hardware concurrency; one worker
runs sequentially. Native matching uses standard C++ threads.

## Results, objectives and fragments

`SearchResult` contains `matches`, `status`, `optimal`, `complete`, `score`,
`nodesVisited` and `elapsedMillis`. Status is `Completed`, `TimedOut`,
`Cancelled` or `LimitReached`. A stopped call may return a valid incumbent
without proving the optimum. A capped enumeration may prove the objective
while leaving enumeration incomplete.

`Match` contains `atomMapping`, `bondMapping` entries
`{queryFirst, querySecond, targetFirst, targetSecond}`, and `fragments`
containing query atom indices. All indices are zero-based. Explicit hydrogen
vertices remain vertices; bracket hydrogen counts do not add vertices.

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

MCIS is the default: maximise atoms, then common bonds, preserving edges and
nonedges between mapped vertices. MCES maximises selected compatible common
edges, then incident atoms; isolated atoms cannot pad the result. Extra input
edges can be omitted. Connectivity and fragment limits apply to selected
common edges. `Match::toSmiles(query)` writes selected atoms and bonds with
cut valences capped. `toGraph(query)` preserves the original chemical context
for further matching and is distinct from a capped output molecule.

## Chemistry and controls

`MCS::computeTautomerConfidence(query, target, match, chemistry)` reports a
heuristic confidence for selected tautomeric chemistry relaxations. The
supplied pH and solvent apply when called, without changing molecules or
the search objective. Identical selected chemistry returns 1; this value
is not a probability.

`MCS::Options` and `Substructure::Options` contain `chemistry`, `timeoutMs`,
`nodeLimit` and `cancellation`. `Substructure::Options::induced` optionally
requires nonedges to be preserved. MCS uses its explicit objective instead.

`ChemOptions` controls element, bond order, aromaticity, ring, isotope, charge,
mapped tetrahedral and double-bond stereo, and tautomer predicates. Defaults
match elements and strict bond order. Charge, isotope and stereo matching are
opt-in. `tautomerAware` does not change `matchBondOrder`.
`completeRingsOnly` checks both inputs: retained ring atoms and bonds must
belong to a fully selected relevant ring in their own input. A fusion atom
may survive through one complete ring.

`timeoutMs=-1` selects an adaptive budget; zero requests an immediate stop.
`nodeLimit=0` leaves the node count unlimited. Set `cancellation` to a shared
`CancellationToken` and call `cancel()` to stop cooperatively. Preparation and
individual steps can extend elapsed time beyond the deadline.

MCIS accepts finite signed query `atomWeights`, with weight before atom and
bond tie-breaks. Native weights use bounded integer millipoints. MCES rejects
atom weights. `minFragmentAtoms`, `maxFragments` and `excludedTargetAtoms`
refer to selected fragments and original target atom indices.

`MCS::matchConstrained(queries, target, options)` returns per-query results
without reusing claimed target atoms or re-perceiving a residual target.
`MCS::consensus(molecules, options, threshold)` returns a progressive common
graph, verified supporting mappings and support count for a fraction of the
inputs. Consensus is heuristic and reports `optimal=false`.

`MCS::canonicalizeMapping` chooses an atom mapping representative under
molecular symmetry and raises if the required symmetry calculation cannot
complete. `MCS::isMappingMaximal` checks for an extension retaining the supplied
correspondence; maximality is weaker than an optimum proof.

## Graphs, MOL and SDF

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

A complete importer must also supply charge, isotope, hydrogen, ring,
aromatic and stereo metadata. MOL round trips do not imply full MDL query
chemistry support. The optional native adapter requires its external SDK;
Python conversion tracks original external atom indices separately.

## Fingerprints and assignment

```cpp
#include <cassert>
#include "smsd/smsd.hpp"

int main() {
    const auto molecule = smsd::parseSMILES("c1ccc(O)cc1");
    const auto fingerprint = smsd::batch::detail::computeCircularFingerprintECFP(molecule, 2, 2048);
    assert(fingerprint.size() == 32);
    assert(smsd::batch::fingerprintTanimoto(fingerprint, fingerprint) == 1.0);
    const auto assignment = smsd::optimalAssign({{1.0, 3.0}, {4.0, 2.0}});
    assert(assignment.totalCost == 3.0);
}
```

Packed fingerprints contain 64-bit words; count fingerprints contain integer
bins. Compare using the same family, radius and size. Screening scores and
canonical hashes do not prove molecular identity. See the
[C++ README](../cpp/README.md) for build configuration and the
[examples](EXAMPLES.md) for other utilities.
