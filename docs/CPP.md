# SMSD Pro C++ Guide

SMSD Pro’s C++ layer is header-only and provides the native implementations for
MolGraph construction, substructure search, MCS, fingerprints, SMARTS matching,
molfile I/O, stereo/CIP assignment, and layout utilities. The checkout targets
version 7.2.1, with current changes under Unreleased. Build this source for those
changes; a version label alone does not identify the reviewed source snapshot.
C++ stays under `cpp/`; its CMake package is independent of the Java Maven
module. Python extension builds use this C++ tree through the root
`pyproject.toml`. The 7.2.0 measurements remain historical; new release checks
are tracked in [7.2.1 validation](VALIDATION_7.2.1.md).

## Include

```cpp
#include "smsd/smsd.hpp"
```

Installed CMake packages can be consumed without manually setting language or
OpenMP flags:

```cmake
find_package(smsd 7.2 CONFIG REQUIRED)
target_link_libraries(my_program PRIVATE smsd::smsd)
```

## Core Use

```cpp
#include "smsd/smsd.hpp"

auto q = smsd::parseSMILES("c1ccccc1");
auto t = smsd::parseSMILES("c1ccc(O)cc1");

smsd::ChemOptions chem;
smsd::MCSOptions mcsOpts;
mcsOpts.timeoutMs = 10000;

bool hit = smsd::isSubstructure(q, t, chem, 10000);
auto mcs = smsd::findMCS(q, t, chem, mcsOpts);
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
smsd::ChemOptions chemistry;
chemistry.matchFormalCharge = true;
chemistry.matchIsotope = true;
chemistry.useChirality = true;
chemistry.matchBondOrder = smsd::ChemOptions::BondOrderMode::STRICT;
chemistry.aromaticityMode = smsd::ChemOptions::AromaticityMode::STRICT;

smsd::MCSOptions options;
options.timeoutMs = 1000;       // milliseconds, not seconds
options.connectedOnly = true;
options.induced = false;
auto mapping = smsd::findMCS(q, t, chemistry, options);
auto errors = smsd::validateMapping(q, t, mapping, chemistry);
```

For an external molecule representation, `MolGraph::Builder` preserves the
order of the property arrays you supply. A minimal heavy-atom connectivity
example is:

```cpp
auto query = smsd::MolGraph::Builder()
    .atomCount(4)
    .atomicNumbers({7, 6, 6, 8})
    .setNeighbors({{1}, {0, 2}, {1, 3}, {2}})
    .setBondOrders({{1}, {1, 1}, {1, 1}, {1}})
    .build();
auto target = smsd::parseSMILES("CC(O)CN");
auto mapping = smsd::findMCS(query, target, smsd::ChemOptions{}, smsd::MCSOptions{});
```

A complete importer must also preserve charge, isotope, hydrogen, ring,
aromatic and stereo metadata. The optional `smsd/rdkit_adapter.hpp` API and
`smsd::smsd_rdkit` CMake target require a compatible RDKit development
installation and C++20. Enable them with `SMSD_WITH_RDKIT=ON`. `fromRDKit` imports atom
numbers, charge, isotope and ring/aromatic flags, but does not import hydrogen
counts or double-bond stereo. Explicit hydrogens are removed by default,
changing indices; unsupported bond types become single bonds. Aromaticity is
re-perceived, and tetrahedral tags lack full ligand-order normalization. Python's
`smsd.from_rdkit` is the maintained RDKit conversion entry point and tracks
original RDKit atom indices; see [the Python guide](PYTHON.md).

## Batch matching and build configuration

The batch namespace retains target order and applies the search budget to
each query/target pair. Parse once when performing repeated queries:

```cpp
std::vector<smsd::MolGraph> targets{
    smsd::parseSMILES("CC(O)CN"), smsd::parseSMILES("CCCC")};
smsd::MCSOptions options;
options.timeoutMs = 1000;
auto mappings = smsd::batch::batchMCS(q, targets, smsd::ChemOptions{}, options, 1);
auto counts = smsd::batch::batchMCSSize(q, targets, smsd::ChemOptions{}, options, 1);
```

The final argument is the worker count: `1` is sequential and `0` uses OpenMP
defaults. Core matching uses CPU/OpenMP; GPU domain kernels are optional.
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
#include "smsd/batch.hpp"

auto q = smsd::parseSMILES("c1ccc(O)cc1");

// Circular (ECFP / FCFP) — binary and count-based
auto ecfp  = smsd::batch::detail::computeCircularFingerprintECFP(q, 2, 2048);
auto fcfp  = smsd::batch::detail::computeCircularFingerprintFCFP(q, 2, 2048);
auto ecfpc = smsd::batch::detail::computeCircularFingerprintECFPCount(q, 2, 2048);
auto fcfpc = smsd::batch::detail::computeCircularFingerprintFCFPCount(q, 2, 2048);

// Path / topological torsion / MACCS
auto pathFP  = smsd::batch::detail::computePathFingerprint(q, 7, 2048);
auto torsion = smsd::batch::detail::computeTopologicalTorsion(q, 2048);
auto maccs   = smsd::batch::detail::computeMACCSKeys(q);

// Similarity
double tani = smsd::batch::detail::tanimoto(ecfp, fcfp);
double dice = smsd::batch::detail::dice(ecfp, fcfp);
```

> **Note.** The pre-7.1.1 `fp/mol/circular.hpp`, `fp/mol/path.hpp`,
> `fp/mol/pharmacophore.hpp`, and `fp/mol/torsion.hpp` headers were
> unmaintained shims that drifted from the real Python/Java bit pattern.
> They have been removed. Use the `smsd::batch::detail::*` entry points
> documented above — these are the exact functions the Python binding and
> the Java `FingerprintEngine` are tested against.

## Public MCS / Substructure Entry Points

The high-level entry points are `smsd::findMCS()`, `smsd::findSubstructure()`
and `smsd::isSubstructure()` declared in `smsd/mcs.hpp` and `smsd/vf2pp.hpp`.
The internal solver headers (`smsd/clique_solver.hpp`, the partition-refinement
backtracker, edge-growth refinement) are private implementation details whose
signatures may change between minor releases — do not depend on them in
out-of-tree code.

```cpp
#include "smsd/mcs.hpp"
#include "smsd/vf2pp.hpp"

auto mapping     = smsd::findMCS(g1, g2, smsd::ChemOptions{}, smsd::MCSOptions{});
auto sub_mapping = smsd::findSubstructure(query, target, smsd::ChemOptions{});
bool contained   = smsd::isSubstructure(query, target, smsd::ChemOptions{});
auto all_maps    = smsd::findAllSubstructures(query, target, smsd::ChemOptions{});
```

`findAllSubstructures` includes distinct atom mappings related by molecular
symmetry, including self matches. It returns up to 10,000 mappings within the
requested time budget; a timeout can return a partial enumeration.

See [the algorithm review](ALGORITHM_REVIEW.md) for regression oracles and local
validation commands. Current timing, quality and cancellation observations are
in the [benchmark report](../benchmarks/RESULTS_7.2.0.md); timings with different
matching policies are not pooled into a headline speedup.

Weighted objectives and bond maximization use objective-aware component and
fragment selection. Signed weights can favor a smaller subgraph than an
identity mapping. Small graphs use bounded exact exploration; reaching a time
or node limit does not prove an optimum.

`canonicalizeMapping` computes automorphism-orbit closure and reports incomplete
generators or exceeded storage bounds with `std::length_error`, and deadline
expiry with `std::runtime_error`. Internal
mapping deduplication retains raw keys when canonicalization cannot complete.
Weighted mapping deduplication keeps query weights attached to their vertices.

## Scaffold Library (7.1.0)

```cpp
#include "smsd/scaffold_library.hpp"
auto scaffold = smsd::scaffold::murckoScaffold(mol);
```

## Hungarian Algorithm (7.1.0)

Optimal assignment solver for atom matching cost matrices.

```cpp
#include "smsd/hungarian.hpp"
auto assignment = smsd::optimalAssign(costMatrix);
```

For an `m × n` matrix, the solver assigns `min(m,n)` pairs using
`O(min(m,n)² max(m,n))` time and `O(m+n)` auxiliary space. A uniform unmatched
penalty does not change the selected pairs. `totalCost` sums real assigned costs
and excludes unmatched penalties. Ragged matrices and nonfinite inputs raise
`std::invalid_argument`; reduced-cost arithmetic overflow raises `std::overflow_error`.

## SMARTS and CIP

```cpp
auto query = smsd::parseSMARTS("[#6]~[#7]");
auto rs = smsd::cip::assignRSAll(q);
auto ez = smsd::cip::assignEZAll(q);
```

## Native MOL/SDF I/O

```cpp
auto mol = smsd::readMolBlock(molBlockText);
std::string v2000 = smsd::writeMolBlock(mol);
std::string v3000 = smsd::writeMolBlockV3000(mol);
std::string sdf = smsd::writeSDFRecord(mol);
```

The native I/O path in `7.1.0` covers practical V2000/V3000 graph round-trip,
metadata, SDF properties, atom maps/classes, and patent-style `R#` handling.
The most exotic MDL query chemistry features are still intentionally documented
as out of scope until they are implemented natively.
