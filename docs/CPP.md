# SMSD C++ guide

SMSD 7.2.2 provides a C++17 header-only core for
MolGraph construction, substructure search, MCS, fingerprints, SMARTS matching,
molfile I/O, stereo/CIP assignment, and layout utilities.
C++ stays under `cpp/`; its CMake package is independent of the Java Maven
module. Python extension builds use this C++ tree through the root
`pyproject.toml`. See [the C++ README](../cpp/README.md) for build, test and
installation commands, and [release validation](VALIDATION_7.2.2.md) for
tested platforms and evidence. The 7.2.0 measurements remain historical.

The examples below are standalone C++17 programs. Save a block as
`example.cpp`, then compile it from the repository root on macOS/Linux:

```bash
c++ -std=c++17 -O2 -I cpp/include example.cpp -o example
./example
```

For Windows or an installed CMake package, use the consumer project in
[the C++ README](../cpp/README.md#build-test-and-install).

## Include

Include `smsd/smsd.hpp` for the main API. SMARTS requires its own
`smsd/smarts_parser.hpp` header.

Installed CMake packages can be consumed without manually setting language or
OpenMP flags:

```cmake
find_package(smsd 7.2 CONFIG REQUIRED)
target_link_libraries(my_program PRIVATE smsd::smsd)
```

## Core Use

```cpp
#include <cassert>
#include "smsd/smsd.hpp"

int main() {
    const auto query = smsd::parseSMILES("c1ccccc1");
    const auto target = smsd::parseSMILES("c1ccc(O)cc1");
    const smsd::ChemOptions chemistry;
    smsd::MCSOptions options;
    options.timeoutMs = 1000;

    assert(smsd::isSubstructure(query, target, chemistry, 1000));
    const auto mapping = smsd::findMCS(query, target, chemistry, options);
    assert(mapping.size() == 6);
    assert(smsd::validateMapping(query, target, mapping, chemistry).empty());
}
```

## Workflow for RDKit users

SMSD uses a query/target workflow, chemistry options and zero-based atom pairs,
but its native types are `smsd::MolGraph` and `smsd::MCSOptions`. Parsing and
search results use the supplied graph's order:

| Operation | Native SMSD API |
|---|---|
| Parse a molecule | `smsd::parseSMILES(smiles)` |
| Test a query against a target | `smsd::isSubstructure(query, target, chemistry, timeoutMs)` |
| Return one substructure embedding | `smsd::findSubstructure(query, target, chemistry, timeoutMs)` |
| Enumerate embeddings | `smsd::findAllSubstructures(query, target, chemistry, timeoutMs)` |
| Compute MCS | `smsd::findMCS(query, target, chemistry, searchOptions)` |

Substructure mappings are vectors of `(query_atom, target_atom)` pairs; MCS
returns `std::map<int, int>`. These are mappings of graph vertices. RDKit FMCS
instead exposes a common query pattern, whose embeddings need a separate
mapping step and may omit edges that SMSD requires between mapped vertices.

Set chemistry separately from the search objective and budget:

```cpp
#include <cassert>
#include "smsd/smsd.hpp"

int main() {
    const auto query = smsd::parseSMILES("c1ccccc1");
    const auto target = smsd::parseSMILES("c1ccc(O)cc1");
    smsd::ChemOptions chemistry;
    chemistry.matchFormalCharge = true;
    chemistry.matchIsotope = true;
    chemistry.useChirality = true;
    chemistry.matchBondOrder = smsd::ChemOptions::BondOrderMode::STRICT;
    chemistry.aromaticityMode = smsd::ChemOptions::AromaticityMode::STRICT;

    smsd::MCSOptions options;
    options.timeoutMs = 1000;
    options.connectedOnly = true;
    options.induced = false;
    const auto mapping = smsd::findMCS(query, target, chemistry, options);
    assert(mapping.size() == 6);
    assert(smsd::validateMapping(query, target, mapping, chemistry).empty());
}
```

For an external molecule representation, `MolGraph::Builder` preserves the
order of the property arrays you supply. A minimal heavy-atom connectivity
example is:

```cpp
#include <cassert>
#include "smsd/smsd.hpp"

int main() {
    const auto query = smsd::MolGraph::Builder()
        .atomCount(4)
        .atomicNumbers({7, 6, 6, 8})
        .setNeighbors({{1}, {0, 2}, {1, 3}, {2}})
        .setBondOrders({{1}, {1, 1}, {1, 1}, {1}})
        .build();
    const auto target = smsd::parseSMILES("CC(O)CN");
    const auto mapping = smsd::findMCS(
        query, target, smsd::ChemOptions{}, smsd::MCSOptions{});
    assert(mapping.size() == 4);
    assert(smsd::validateMapping(query, target, mapping, smsd::ChemOptions{}).empty());
}
```

A complete importer must also preserve charge, isotope, hydrogen, ring,
aromatic and stereo metadata. Connectivity alone does not establish equivalent
chemistry.

## RDKit integration

The optional `smsd/rdkit_adapter.hpp` API and
`smsd::smsd_rdkit` CMake target require a compatible RDKit development
installation and C++20. Enable them with `SMSD_WITH_RDKIT=ON`. `fromRDKit` imports atom
numbers, charge, isotope and ring/aromatic flags, but does not import hydrogen
counts or double-bond stereo. Explicit hydrogens are removed by default,
changing indices; unsupported bond types become single bonds. Aromaticity is
re-perceived, and tetrahedral tags lack full ligand-order normalisation. Python's
`smsd.from_rdkit` is the maintained RDKit conversion entry point and tracks
original RDKit atom indices; see [the Python guide](PYTHON.md).

## Batch matching and build configuration

The batch namespace retains target order and applies the search budget to
each query/target pair. Parse once when performing repeated queries:

```cpp
#include <cassert>
#include <vector>
#include "smsd/smsd.hpp"

int main() {
    const auto query = smsd::parseSMILES("CCO");
    const std::vector<smsd::MolGraph> targets{
        smsd::parseSMILES("CCO"), smsd::parseSMILES("CCCC")};
    smsd::MCSOptions options;
    options.timeoutMs = 1000;
    const auto mappings = smsd::batch::batchMCS(
        query, targets, smsd::ChemOptions{}, options, 1);
    const auto counts = smsd::batch::batchMCSSize(
        query, targets, smsd::ChemOptions{}, options, 1);
    assert((counts == std::vector<int>{3, 2}));
    assert(mappings.size() == targets.size());
    for (std::size_t i = 0; i < targets.size(); ++i) {
        assert(mappings[i].size() == static_cast<std::size_t>(counts[i]));
        assert(smsd::validateMapping(
            query, targets[i], mappings[i], smsd::ChemOptions{}).empty());
    }
}
```

The final argument is the worker count: `1` uses one worker and `0` uses OpenMP
defaults. Without OpenMP, operations run sequentially. Core matching uses
CPU/OpenMP; GPU domain kernels require a source build and their own hardware
validation. Published Python wheels use CPU/OpenMP.
For a reproducible CPU build:

```bash
cmake -S cpp -B build/native -DCMAKE_BUILD_TYPE=Release \
  -DSMSD_BUILD_TESTS=OFF -DSMSD_BUILD_PYTHON=OFF \
  -DSMSD_BUILD_METAL=OFF -DSMSD_BUILD_CUDA=OFF -DSMSD_WITH_RDKIT=OFF
cmake --build build/native --parallel 4
cmake --install build/native --prefix "$PWD/build/native-install"
```

Mappings are oriented from query to target and validate every mapped query
bond. In non-induced mode, extra target bonds are allowed, so the best valid
mapping can have a different size when the arguments are reversed. Non-induced
searches retain the caller's query orientation. Atom weights refer to query indices.
Weights must be finite and cover every query atom. The native scoring API uses
integer millipoints: positive and negative weight sums must each fit that score
range after multiplication by 1,000 and truncation. Unsupported inputs raise
`std::invalid_argument` before search, including on identity fast paths.

## Fingerprints

The shared fingerprint functions used by the Python binding live in
`smsd::batch::detail`. Regression fixtures compare their outputs with Java;
this does not establish identical results for every molecule representation,
option or radius.

```cpp
#include <cassert>
#include <cmath>
#include <cstdint>
#include <vector>
#include "smsd/smsd.hpp"

int main() {
    const auto molecule = smsd::parseSMILES("c1ccc(O)cc1");
    const auto ecfp = smsd::batch::detail::computeCircularFingerprintECFP(molecule, 2, 2048);
    const auto fcfp = smsd::batch::detail::computeCircularFingerprintFCFP(molecule, 2, 2048);
    const auto ecfpCounts = smsd::batch::detail::computeCircularFingerprintECFPCounts(molecule, 2, 2048);
    const auto fcfpCounts = smsd::batch::detail::computeCircularFingerprintFCFPCounts(molecule, 2, 2048);
    const auto path = smsd::batch::detail::computePathFingerprint(molecule, 7, 2048);
    const auto torsion = smsd::batch::detail::computeTopologicalTorsion(molecule, 2048);
    assert(ecfp.size() == 32 && fcfp.size() == 32);
    assert(path.size() == 32 && torsion.size() == 32);
    assert(ecfpCounts.size() == 2048 && fcfpCounts.size() == 2048);
    assert(smsd::batch::fingerprintTanimoto(ecfp, ecfp) == 1.0);
    assert(smsd::batch::countTanimoto(ecfpCounts, ecfpCounts) == 1.0);

    const std::vector<std::uint64_t> a{0b11}, b{0b10};
    assert(smsd::batch::fingerprintTanimoto(a, b) == 0.5);
    assert(std::abs(smsd::batch::fingerprintDice(a, b) - 2.0 / 3.0) < 1e-12);
}
```

Binary fingerprints contain packed 64-bit words; count fingerprints contain
integer bins. Compare molecules using the same fingerprint family, radius
and size. `fingerprintTanimoto` computes intersection divided by union;
`fingerprintDice` computes twice the intersection divided by the sum of set
bits. Count Tanimoto uses bin-wise minima and maxima. These measures differ
from the overlap coefficient.

The removed `fp/mol/*.hpp` shims are not part of the current API. There is
no native `computeMACCSKeys` entry point in 7.2.2.

## Public MCS / Substructure Entry Points

The high-level entry points are `smsd::findMCS()`, `smsd::findSubstructure()`
and `smsd::isSubstructure()` declared in `smsd/mcs.hpp` and `smsd/vf2pp.hpp`.
The internal solver headers (`smsd/clique_solver.hpp`, the partition-refinement
backtracker, edge-growth refinement) are private implementation details whose
signatures may change between minor releases — do not depend on them in
out-of-tree code.

```cpp
#include <cassert>
#include "smsd/smiles_parser.hpp"
#include "smsd/mcs.hpp"
#include "smsd/vf2pp.hpp"

int main() {
    const auto query = smsd::parseSMILES("c1ccccc1");
    const auto target = smsd::parseSMILES("c1ccc(O)cc1");
    const auto mapping = smsd::findMCS(query, target, smsd::ChemOptions{}, smsd::MCSOptions{});
    const auto embedding = smsd::findSubstructure(query, target, smsd::ChemOptions{});
    const auto allEmbeddings = smsd::findAllSubstructures(query, target, smsd::ChemOptions{});
    assert(smsd::isSubstructure(query, target, smsd::ChemOptions{}));
    assert(mapping.size() == 6 && embedding.size() == 6);
    assert(!allEmbeddings.empty());
}
```

`findAllSubstructures` includes distinct atom mappings related by molecular
symmetry, including self matches. It returns up to 10,000 mappings within the
requested time budget; a timeout can return a partial enumeration.

See [the algorithm review](ALGORITHM_REVIEW.md) for regression oracles and local
validation commands. Historical timing, quality and cancellation observations are
in the [benchmark report](../benchmarks/RESULTS_7.2.0.md); timings with different
matching policies are not pooled into a headline speedup.

Weighted objectives and bond maximisation use objective-aware component and
fragment selection. Signed weights can favour a smaller subgraph than an
identity mapping. Small graphs use bounded exact exploration; reaching a time
or node limit does not prove an optimum.

`canonicalizeMapping` computes automorphism-orbit closure and reports incomplete
generators or exceeded storage bounds with `std::length_error`, and deadline
expiry with `std::runtime_error`. Internal
mapping deduplication retains raw keys when canonicalisation cannot complete.
Weighted mapping deduplication keeps query weights attached to their vertices.

## Scaffold extraction

`smsd::murckoScaffold` retains ring systems and connecting paths and removes
side chains. It returns an acyclic input unchanged; callers expecting an
empty scaffold for acyclic molecules must handle that case separately.

```cpp
#include <cassert>
#include "smsd/smsd.hpp"

int main() {
    const auto molecule = smsd::parseSMILES("Cc1ccccc1");
    const auto scaffold = smsd::murckoScaffold(molecule);
    assert(scaffold.n == 6);
    assert(smsd::isSubstructure(scaffold, molecule, smsd::ChemOptions{}));
}
```

`smsd/scaffold_library.hpp` separately supplies reference scaffold data;
it does not contain the extraction function.

## Optimal assignment

Optimal assignment solver for atom matching cost matrices.

```cpp
#include <cassert>
#include <vector>
#include "smsd/hungarian.hpp"

int main() {
    const std::vector<std::vector<double>> costs{{1.0, 3.0}, {4.0, 2.0}};
    const auto result = smsd::optimalAssign(costs);
    assert(result.assignment.size() == 2);
    assert(result.totalCost == 3.0);
}
```

For an `m × n` matrix, the solver assigns `min(m,n)` pairs using
`O(min(m,n)² max(m,n))` time and `O(m+n)` auxiliary space. A uniform unmatched
penalty does not change the selected pairs. `totalCost` sums real assigned costs
and excludes unmatched penalties. Ragged matrices and nonfinite inputs raise
`std::invalid_argument`; reduced-cost arithmetic overflow raises `std::overflow_error`.

## SMARTS and CIP

```cpp
#include <cassert>
#include <tuple>
#include "smsd/smiles_parser.hpp"
#include "smsd/smarts_parser.hpp"
#include "smsd/cip.hpp"

int main() {
    const auto query = smsd::parseSMARTS("[#6]~[#7]");
    const auto target = smsd::parseSMILES("CCN");
    const auto matches = query.findAll(target, 100);
    assert(!matches.empty() && matches.front().size() == 2);

    const auto alanine = smsd::parseSMILES("N[C@@H](C)C(=O)O");
    const auto alanineDescriptors = smsd::cip::assignAll(alanine);
    assert(alanineDescriptors.rsLabels[1] == smsd::cip::RSLabel::S);
    const auto butene = smsd::parseSMILES("C/C=C/C");
    const auto buteneDescriptors = smsd::cip::assignAll(butene);
    assert(buteneDescriptors.ezBonds.size() == 1);
    assert(std::get<2>(buteneDescriptors.ezBonds.front()) == smsd::cip::EZLabel::E);
}
```

`smsd::cip::assignAll` returns per-atom R/S labels in `rsLabels` and
stereogenic double-bond atom pairs with E/Z labels in `ezBonds`.

## Native MOL/SDF I/O

```cpp
#include <cassert>
#include <string>
#include "smsd/smsd.hpp"

int main() {
    const auto molecule = smsd::parseSMILES("CCO");
    const std::string v2000 = smsd::writeMolBlock(molecule);
    const std::string v3000 = smsd::writeMolBlockV3000(molecule);
    const std::string sdf = smsd::writeSDFRecord(molecule);
    const auto fromV2000 = smsd::readMolBlock(v2000);
    const auto fromV3000 = smsd::readMolBlock(v3000);
    assert(fromV2000.n == 3 && fromV3000.n == 3);
    assert(smsd::isSubstructure(molecule, fromV2000, smsd::ChemOptions{}));
    assert(smsd::isSubstructure(molecule, fromV3000, smsd::ChemOptions{}));
    assert(sdf.find("$$$$") != std::string::npos);
}
```

The native I/O path covers practical V2000/V3000 graph round-trip,
metadata, SDF properties, atom maps/classes, and patent-style `R#` handling.
Graph round-trip support does not imply complete MDL query chemistry support.
