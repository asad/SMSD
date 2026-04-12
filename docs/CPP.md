# SMSD Pro C++ Guide

SMSD Pro’s C++ layer is header-only and provides the native implementations for
MolGraph construction, substructure search, MCS, fingerprints, SMARTS matching,
molfile I/O, stereo/CIP assignment, and layout utilities.

## Include

```cpp
#include "smsd/smsd.hpp"
```

## Core Use

```cpp
#include "smsd/smsd.hpp"

auto q = smsd::parseSMILES("c1ccccc1");
auto t = smsd::parseSMILES("c1ccc(O)cc1");

smsd::ChemOptions chem;
smsd::MCSOptions mcsOpts;
mcsOpts.timeout_ms = 10000;

bool hit = smsd::isSubstructure(q, t, chem, 10000);
auto mcs = smsd::findMCS(q, t, chem, mcsOpts);
```

## Fingerprints

```cpp
// Standalone fingerprint modules (6.12.1)
#include "fp/mol/circular.hpp"
#include "fp/mol/path.hpp"
#include "fp/mol/pharmacophore.hpp"
#include "fp/mol/torsion.hpp"
#include "fp/similarity.hpp"

auto ecfp = fp::mol::computeCircularFingerprintECFP(q, 2, 2048);
auto fcfp = fp::mol::computeCircularFingerprintFCFP(q, 2, 2048);
auto pathFp = fp::mol::computePathFingerprint(q, 7, 1024);
auto torsion = fp::mol::computeTopologicalTorsion(q, 2048);
double tani = fp::fingerprintTanimoto(ecfp, pathFp);
double dice = fp::fingerprintDice(ecfp, pathFp);
bool subset = fp::fingerprintSubset(ecfp, pathFp);

// Batch API (via batch.hpp, delegates to fp/)
auto pathFp = smsd::batch::detail::computePathFingerprint(q, 7, 2048);
auto ecfp = smsd::batch::detail::computeCircularFingerprintECFP(q, 2, 2048);
```

## Clique Solver (6.12.1)

Lightweight maximum clique finder for Python/RDKit integration where
chemistry stays in the caller and only the combinatorial search runs
in C++.

```cpp
#include "smsd/clique_solver.hpp"

smsd::clique::ProductGraph pg;
pg.build(vertices, edges);

auto result = smsd::clique::findMaxCliques(pg, 8, 1000);
// result.cliques, result.max_size, result.timed_out

auto mcs = smsd::clique::findMCSPipeline(
    compat, bonds_a, bonds_b, n_a, n_b);
// mcs.candidates, mcs.best_size, mcs.lfub

auto sub = smsd::clique::substructureMatch(
    n_query, n_target, compat, bonds_query, bonds_target);
```

## Scaffold Library (6.12.1)

```cpp
#include "smsd/scaffold_library.hpp"
auto scaffold = smsd::scaffold::murckoScaffold(mol);
```

## Hungarian Algorithm (6.12.1)

Optimal O(n^3) assignment for atom matching cost matrices.

```cpp
#include "smsd/hungarian.hpp"
auto assignment = smsd::hungarian::solve(costMatrix);
```

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

The native I/O path in `6.12.1` covers practical V2000/V3000 graph round-trip,
metadata, SDF properties, atom maps/classes, and patent-style `R#` handling.
The most exotic MDL query chemistry features are still intentionally documented
as out of scope until they are implemented natively.
