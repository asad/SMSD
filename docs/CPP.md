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

The canonical C++ fingerprint API lives in `smsd::batch::detail` and is
verified byte-identical with the Java and Python engines in 7.1.1 (same
bits for the same molecules at every radius).

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
and `smsd::isSubstructure()` declared in `smsd/mcs.hpp` and `smsd/substructure.hpp`.
The internal solver headers (`smsd/clique_solver.hpp`, the partition-refinement
backtracker, edge-growth refinement) are private implementation details whose
signatures may change between minor releases — do not depend on them in
out-of-tree code.

```cpp
#include "smsd/mcs.hpp"
#include "smsd/substructure.hpp"

auto mapping     = smsd::findMCS(g1, g2, smsd::ChemOptions{}, smsd::MCSOptions{});
auto sub_mapping = smsd::findSubstructure(query, target, smsd::ChemOptions{});
bool contained   = smsd::isSubstructure(query, target, smsd::ChemOptions{});
```

## Scaffold Library (7.1.0)

```cpp
#include "smsd/scaffold_library.hpp"
auto scaffold = smsd::scaffold::murckoScaffold(mol);
```

## Hungarian Algorithm (7.1.0)

Optimal assignment solver for atom matching cost matrices.

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

The native I/O path in `7.1.0` covers practical V2000/V3000 graph round-trip,
metadata, SDF properties, atom maps/classes, and patent-style `R#` handling.
The most exotic MDL query chemistry features are still intentionally documented
as out of scope until they are implemented natively.
