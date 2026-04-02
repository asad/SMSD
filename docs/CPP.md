# SMSD C++ Guide

SMSD’s C++ layer is header-only and provides the native implementations for
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
smsd::McsOptions mcsOpts;
mcsOpts.timeout_ms = 10000;

bool hit = smsd::isSubstructure(q, t, chem, 10000);
auto mcs = smsd::findMCS(q, t, chem, mcsOpts);
```

## Fingerprints

```cpp
auto pathFp = smsd::batch::detail::computePathFingerprint(q, 7, 2048);
auto ecfp = smsd::batch::detail::computeCircularFingerprintECFP(q, 2, 2048);
bool subset = smsd::batch::fingerprintSubset(pathFp, pathFp);
double tanimoto = smsd::batch::fingerprintTanimoto(pathFp, pathFp);
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

The native I/O path in `6.9.1` covers practical V2000/V3000 graph round-trip,
metadata, SDF properties, atom maps/classes, and patent-style `R#` handling.
The most exotic MDL query chemistry features are still intentionally documented
as out of scope until they are implemented natively.
