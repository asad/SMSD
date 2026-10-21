# SMSD Pro 8.0.0

SMSD Pro 8.0.0 provides substructure and maximum common substructure (MCS)
search for Java, C++ and Python through prepared `MCS` and `Substructure`
queries. A query is prepared once and searches one or more targets; every
result carries its matches with their atom and bond mappings, the search
status and whether the optimum was proved. The
[changelog](https://github.com/asad/SMSD/blob/master/CHANGELOG.md) describes
the release in full.

## Highlights

- One search model in the three languages: connected and disconnected
  maximum common induced subgraphs (MCIS), maximum common edge subgraphs
  (MCES), weighted MCIS and substructure search, with multi-target,
  constrained, consensus and R-group searches in the same query APIs.
- Results that report completion, stopping status and optimality separately
  from the matches, under timeouts, node limits and cancellation tokens.
- Daylight SMARTS semantics in C++ and Java, published SMARTS size limits in
  every language, and the OpenSMILES limits for SMILES bracket atoms in C++
  and Python.
- One valence table, one tetrahedral-centre rule and CIP descriptors that
  follow the IUPAC 2013 Sequence Rules 1 to 5, with Rule 6 as proposed by
  Hanson et al. (2018), for R/S, shared by C++, Java and Python.
- Tautomer classes that do not depend on the atom order, with pH-dependent
  weights from the Henderson–Hasselbalch relation and a tautomer confidence
  for the selected bonds of a match.
- `FingerprintOptions`, one configurable fingerprint specification for
  circular, path, torsion and substructure-key fingerprints, including the 166
  MACCS keys, with identical bits in C++, Java and Python.
- The Java library `com.bioinceptionlabs:smsd`, which depends only on the CDK
  modules it uses, and the separate command-line tool
  `com.bioinceptionlabs:smsd-cli`.

## Downloads and platform requirements

- **Python:** CPython 3.14 CPU/OpenMP wheels for Windows 10 or later on x86_64
  (AMD64), Linux x86_64 (AMD64) with glibc 2.28 or later, and macOS 26 or later
  on arm64. The source distribution `smsd-8.0.0.tar.gz` builds with Python 3.9
  or later, a C++17 compiler and CMake 3.18 or later.
- **Java:** the library `smsd-8.0.0.jar` with its sources and Javadoc JARs,
  the self-contained CLI `smsd-cli-8.0.0-jar-with-dependencies.jar`, and the
  portable CLI `smsd-8.0.0-cli.tar.gz`. They run on Java 8 or later (Java 25
  LTS preferred) with CDK 2.13. On Maven Central, the library is
  `com.bioinceptionlabs:smsd` and the CLI `com.bioinceptionlabs:smsd-cli`.
- **Java CLI installers:** Windows MSI on x86_64 (AMD64), macOS DMG for macOS
  11.0 or later on arm64, and Linux DEB on x86_64 (AMD64), each with a bundled
  Eclipse Temurin Java 25.0.4.1+1 runtime.
- **C++:** headers in `smsd-cpp-8.0.0-headers.tar.gz`, for a C++17 compiler.

Each release file is listed with its SHA-256 checksum in `SHA256SUMS`, and the
release's validation record is published with it as `VALIDATION.md`; release
checks are recorded there as they complete. The
[installation guide](https://github.com/asad/SMSD/blob/master/docs/HOWTO-INSTALL.md)
explains how to verify downloads, install each package and build from source.

## Upgrading from 7.x

Code written for 7.x must be updated: the 7.2.2 search entry points below are
removed. Prepare a query once, then search one or more targets with it.
`withOptions(options)` (Python `with_options`) takes the search options, whose
`chemistry` field holds the chemistry settings that the 7.2.2 calls took as
arguments. Java `matchMany` and `MCS.matchConstrained` take `MolGraph` lists;
`new MolGraph(container)` converts a CDK molecule.

| 7.2.2 | 8.0.0 |
| --- | --- |
| Python `smsd.find_mcs(a, b)`, `smsd.mcs_result(a, b)`, `smsd.mcs_from_smiles(a, b)`, `smsd.mcs_rdkit(a, b)`, `smsd.mcs_rdkit_native(a, b)` | `smsd.MCS.find(a).match(b)`, which takes SMILES strings, `MolGraph` objects and RDKit molecules |
| Python `smsd.find_substructure(q, t)`, `smsd.is_substructure(q, t)`, `smsd.substructure_rdkit(q, t)` | `smsd.Substructure.find(q).match(t)`, then `.found` or `.best` |
| Python `smsd.batch_mcs(q, targets)`, `smsd.batch_mcs_size(q, targets)`, `smsd.batch_mcs_rdkit(q, targets)`, `smsd.batch_substructure(q, targets)`, `smsd.batch_find_substructure(q, targets)` | `smsd.MCS.find(q).match_many(targets)`, `smsd.Substructure.find(q).match_many(targets)` |
| Python `smsd.screen_and_match(q, targets, threshold)`, `smsd.screen_and_mcs_size(q, targets, threshold)` | `smsd.screen_targets(q, targets, threshold)`, then `smsd.MCS.find(q).match_many` on the targets it selects |
| Python `TargetCorpus.substructure(q)`, `find_substructure(q)` and `mcs_size(q)` | `smsd.Substructure.find(q).match_many(corpus)`, `smsd.MCS.find(q).match_many(corpus)` |
| Python `smsd.mcs_smiles(a, b)` | `smsd.MCS.find(a).match(b)`, then `.best.to_smiles(a)` |
| Python `smsd.find_mcs_smarts(smarts, target)` | `smsd.Substructure.smarts(smarts).match(target)`, then `.best.atom_mapping` |
| Python `smsd.smarts_match(smarts, t)`, `smsd.smarts_find_all(smarts, t)`, `smsd.compile_smarts(smarts)` | `smsd.Substructure.smarts(smarts)`, then `.match(t).found` or `.match_all(t, limit)` |
| Python `smsd.batch_mcs_constrained(queries, targets)` | `smsd.MCS.match_constrained(queries, target)`, for one target |
| Python `smsd.decompose_r_groups(core, molecules)` | `smsd.Substructure.decompose_r_groups(core, molecules)` |
| Python `timeout_ms=` keyword arguments | `.with_timeout(ms)` on the prepared query |
| Java `SearchEngine.findMCS(a, b)`, `new SMSD(q, t, chem).findMCS()` | `MCS.find(a).match(b)` |
| Java `SearchEngine.findMCSSmiles(a, b, chem, options)`, `new SMSD(q, t, chem).findMCSSmiles()` | `MCS.find(a).match(b)`, then `toSmiles(a)` on the first match |
| Java `SearchEngine.findSubstructure(q, t)`, `SearchEngine.isSubstructure(q, t, chem, ms)`, `new SMSD(q, t, chem).isSubstructure()` | `Substructure.find(q).match(t)`, then `found()`, or `atomMapping()` of the first match |
| Java `SearchEngine.findMCSSmarts(smarts, target, chem, ms)` | `Substructure.smarts(pattern).withTimeout(ms).match(target)` |
| Java `SearchEngine.findAllMCS`, `SearchEngine.findAllSubstructures` | `matchAll(target, limit)` on the prepared query |
| Java `SearchEngine.batchSubstructure(q, targets, chem, ms)` | `Substructure.find(q).matchMany(targets, numThreads)`, then `found()` on each result |
| Java `SearchEngine.batchMCS(molecules, chem, options, threshold)`, every pair of a list | For each molecule, `MCS.find(molecule).matchMany(laterMolecules, numThreads)`, after `MCS.similarityUpperBound(a, b, chem)` has removed the pairs below the threshold |
| Java `SearchEngine.batchMCSConstrained(queries, targets, chem, options, ms)` | `MCS.matchConstrained(queries, target, options)`, for one target |
| Java `new SMSD(smarts, target, chem)` | `Substructure.smarts(pattern).match(target)` |
| Java `SMSD.setMCSTimeoutMs(ms)`, `setSubstructureTimeoutMs(ms)` | `withTimeout(ms)` on the prepared query |
| Java `SearchEngine.findDisconnectedMCS(a, b, chem, options)` | `MCS.find(a).withOptions(options).match(b)`, with `options.connectivity = MCS.Connectivity.Fragments` |
| Java `SearchEngine.findNMCS` and `findNMCSMolecule`, and their `SMSD` forms | `MCS.consensus(graphs, options, threshold)` on `MolGraph` copies of the molecules |
| Java `SearchEngine.decomposeRGroups(core, molecules, chem, ms)`, `SMSD.decomposeRGroups` | `Substructure.decomposeRGroups(core, molecules, chem, ms)`, whose results hold `core()`, `rgroups()` and `search()` |
| Java `SearchEngine` and `SMSD` helpers `similarityUpperBound`, `validateMapping`, `isMappingMaximal`, `canonicalizeMapping`, `computeTautomerConfidence`, `validateTautomerConsistency` and `murckoScaffold` | The `MCS` methods of the same names |
| Java `SearchEngine` and `SMSD` fingerprint methods | `FingerprintEngine`, for example `ecfp`, `fcfp`, `pathFingerprint`, `mcsFingerprint`, `topologicalTorsion` and `tanimoto` |
| C++ `smsd::findMCS(a, b, ...)` | `smsd::MCS::find(a).match(b)` |
| C++ `smsd::isSubstructure(q, t, ...)`, `smsd::findSubstructure(q, t, ...)` | `smsd::Substructure::find(q).match(t)` |
| C++ `smsd::findAllMCS`, `smsd::findAllSubstructures` | `matchAll(target, limit)` on the prepared query |
| C++ `smsd::batch::batchMCS`, `smsd::batch::batchSubstructure` | `matchMany(targets, numThreads)` on the prepared query |
| C++ `smsd::findMCSSmiles(a, b, ...)`, `smsd::mcsToSmiles(g, mapping)` | `smsd::MCS::find(a).match(b)`, then `toSmiles(a)` on the first match |
| C++ `smsd::findDisconnectedMCS(a, b, ...)` | `smsd::MCS::find(a).withOptions(options).match(b)`, with `options.connectivity = smsd::MCS::Connectivity::Fragments` |
| C++ `smsd::findNMCS(molecules, ...)` | `smsd::MCS::consensus(molecules, options, threshold)` |
| C++ `smsd::batchMCSConstrained(queries, targets, ...)` | `smsd::MCS::matchConstrained(queries, target, options)`, for one target |
| C++ `smsd::validateMapping`, `smsd::isMappingMaximal`, `smsd::canonicalizeMapping` and `smsd::decomposeRGroups` | The static members of the same names of `smsd::MCS`, and `smsd::Substructure::decomposeRGroups` |
| C++ `smsd::batch::batchFindSubstructure`, `batchMCSSize`, `screenAndMatch` and `screenAndMCSSize`, and the `TargetCorpus` members `substructure`, `findSubstructure` and `mcsSize` | `matchMany(targets, numThreads)` on the prepared query, after `smsd::screenTargets` or `TargetCorpus::screen` for the screening forms |

### Breaking changes

- Prepared `MCS` and `Substructure` queries replace the search facades and
  algorithm selectors in Java, C++ and Python, and results report completion,
  stopping status and optimality separately from the matches.
- The C++ similarity-bound screening header `smsd/rascal.hpp` is
  `smsd/similarity_screen.hpp`, and the GPU entry points `batchRascalScreen`
  and `gpu::batchRascalScreenAuto` of `smsd/gpu.hpp` are
  `batchSimilarityScreen` and `gpu::batchSimilarityScreenAuto`.
- The Java library and command-line tool are separate Maven artefacts.
  `com.bioinceptionlabs:smsd` depends only on the CDK modules that the library
  uses, 10 artefacts at run time instead of the 21 of 7.2.2, and
  `com.bioinceptionlabs:smsd-cli` holds `com.bioinception.smsd.cli` with
  picocli, Jackson and the CDK file readers. Projects that used picocli,
  `cdk-data`, `cdk-io`, `cdk-ioformats`, `cdk-ctab`, `cdk-pdb`, `cdk-smiles` or
  `jackson-databind` through `smsd` declare them. The self-contained CLI JAR is
  `smsd-cli-8.0.0-jar-with-dependencies.jar`.
- SMARTS: C++ bracket atoms with an uppercase element symbol, such as `[C]` or
  `[NX3]`, match aliphatic atoms only, and lowercase symbols aromatic atoms
  only, as in Java; use `[#6]` or `[C,c]` to match both. C++ reads `[H]`,
  `[2H]`, `[H+]` and `[H:1]` as hydrogen atoms rather than hydrogen counts.
  Java directional bonds no longer match the single bonds of an aromatic ring,
  and Java `r<n>` is the smallest ring through the atom. SMARTS beyond the
  published size limits are rejected in every language, and C++ and Python
  reject SMILES bracket atoms outside the OpenSMILES limits, such as `[CH10]`.
- MCES: Java and C++ run every MCES search through one bounded search. Proved
  optima are unchanged, but a search may return a different, equally optimal
  mapping, and Java `matchAll` results are sorted as in C++. A node limit counts
  mapped atom pairs; a search that reaches it returns status `LimitReached` and
  is never reported optimal, and a C++ `matchAll` that stops early is not
  reported optimal.
- Java disconnected MCIS searches may return a different one of several tied
  optimal mappings, and the default Java `Substructure` budget is that of C++
  and Python.
- `Substructure` searches with `useChirality` or `useBondStereo` require a
  specified target centre or bond for a specified query one; an unspecified
  target centre or bond matched it in 7.2.2.
- Building results counts towards the timeout, and Java molecule substructure
  pairs that element and degree counts rule out report `Completed` with no
  nodes visited.
- Stereochemistry and CIP: C++ and Java share one tetrahedral-centre rule, and
  Java `MolGraph.Builder` tags place the implicit hydrogen or lone pair in the
  order C++ uses rather than always first. R/S labels apply Sequence Rules 4 to
  6 and rank mancude ring duplicates by mean atomic number, so some
  pseudoasymmetric, ring and isotope-labelled centres take different labels;
  a direct C++ `assignRS` call ignores a supplied descriptor map.
- Valence and rings: C++ and Java share one implicit-hydrogen valence table,
  which changes SMARTS `H`, `X` and `v` on graphs built without hydrogen counts
  and Java canonical SMILES of thiophene-type sulfur. C++ counts each atom's
  rings over its relevant cycles and keeps the input order of equal-size rings,
  which can change ring order, smallest-ring choices, `R<n>` counts and 2D
  coordinates where rings tie.
- Tautomers: classes no longer depend on the atom order and can hold more
  atoms, and weights follow the Henderson–Hasselbalch relation instead of a
  linear damping by pH. Tautomer confidences, tautomer-aware search results and
  fingerprints that hash the tautomer class change for molecules whose classes
  change. C++ and Python tautomer consistency counts stored hydrogens.
- Murcko scaffolds keep atoms double or triple bonded to the scaffold, replace
  stripped substituents with hydrogen in C++ and Python, and keep
  stereochemistry where its neighbours stay.
- Fingerprints: two empty fingerprints score as identical sets (Tanimoto and
  overlap 1.0, Soergel distance 0.0); Python `circular_fingerprint` rejects
  modes other than `ecfp` and `fcfp`, and `from_hex` rejects non-hexadecimal
  characters.
- Java result details follow the native core: `MCS.consensus` candidates and
  exclusions, `decomposeRGroups` cores, `similarityUpperBound`,
  `matchConstrained` after a cancellation, `isMappingMaximal` for invalid
  input, `validateMapping` messages and order, `canonicalizeMapping`
  exceptions, and any negative `timeoutMillis` as the adaptive timeout. A
  `MolGraph` reads CDK implicit hydrogen counts when it is built. A call that
  passes a bare `null` as the options of `MCS.validateMapping`, `mappingMatch`
  or `isMappingMaximal` needs a cast to `MCS.Options` or `ChemOptions`.
- Python: `MolGraph.canonical_label` and `orbit` compute on first access,
  `TargetCorpus.add_target` and `add_targets` leave the invariants to
  `prewarm()` or the first search, `translate_mapping` raises `IndexError` for
  out-of-range indices, and every `ChemOptions` setter, including
  `with_solvent` and `with_ph`, returns the options it is called on rather than
  a copy.
- The Java CLI exits with status 0 for a match, 1 for no match, and 2 for
  invalid input, an error or an incomplete search, and rejects unknown modes.

### Deprecated names

Each deprecated name keeps working; the Python names warn with
`DeprecationWarning`.

| Deprecated | Replacement |
| --- | --- |
| C++ and Java `ChemOptions.fmcsProfile()` | `ChemOptions.profile("compat-fmcs")`, which returns the same options |
| C++ and Java three-argument `computeTautomerConfidence(query, target, mapping)` | `MCS.computeTautomerConfidence(query, target, match, chemistry)` |
| C++ `smsd::ChemOps` | `MCS::validateMapping` |
| C++ `smsd::implicitHydrogens(z, bondOrderSum, charge)` | `MolGraph::implicitHydrogenCount(atom)`; `allowedValences(z)` gives the per-element valence data |
| C++ `ChemOptions::aromaticityModel`, Python `ChemOptions.aromaticity_model` | The model argument of `MolGraph::Builder::build()` or `perceiveAromaticity()` |
| C++ `MolGraph::ensurePatternFP()`, `getNeighborsByDegDesc()`, `getNLF1()`, `getNLF2()`, `getNLF3()` and `nlfCheckOk` | None; no search reads these caches |
| C++ static `MolGraph::refinePartition`, `comparePerm` and `computeCanonicalLabeling` | `MolGraph::getCanonicalLabeling()` on the graph |
| C++ and Java `MolGraph.TW_NITRO_ACI_NITRO` | `MolGraph.TW_NITRO_ACI` |
| C++ `MolGraph::TW_SULFOXIDE` and `TW_NITRILE_ISONITRILE` | None; no tautomer class uses these weights |
| C++ `hungarianSolve` and `HungarianResult` | `optimalAssign` and `AssignmentResult` |
| C++ `runHortonPipeline` | `computeRings` or `computeRelevantCycles` |
| C++ `popcountAnd`, `popcountN`, `anyBitAnd` and `popcountAndNot` | None |
| C++ `transform::translate3D`, `rotate3DX`, `rotate3DY`, `rotate3DZ`, `scale3D`, `center3D` and `rmsd2D` | None |
| Java `ChemOptions.useTwoHopNLF`, `useThreeHopNLF` and `useBitParallelFeasibility`, with their `with` setters | None; they never change a result |
| Java `FingerprintEngine.clearFPCache()` | None; each `MolGraph` keeps its own default fingerprints |
| Java `FingerprintEngine.circularFingerprint(molecule, radius, size)` | `ecfp` or the `Mode` overload |
| Java `FingerprintEngine.analyzeFingerprintQuality` | `analyzeFpQuality` |
| Java `MolGraph.getOrbits()` | `getOrbit()`, which returns a copy |
| Java `Standardiser.matchAll(smarts, target)` | `Substructure.smarts(smarts).matchAll(target)` |
| Java `ChemOptions` overloads of `MCS.validateMapping` and `MCS.isMappingMaximal` | The `Match` and `MCS.Options` forms |
| Python `smsd.dearomatize()` | `smsd.kekulize()` |
| Python `canonical_smiles`, `overlapCoefficient`, `count_overlapCoefficient`, `prewarm_graph`, `read_molfile`, `export_sdf` and `to_svg` | `to_smiles`, `overlap_coefficient`, `count_overlap_coefficient`, `prewarm`, `read_mol_file`, `write_sdf` and `depict_svg` |
| Python `chem` keyword of `similarity_upper_bound` and `screen_targets` | `chemistry` |
| Python `highlight_color` argument of `depict_mcs` | None; it has no effect |

### Removed entry points

- The 7.2.2 search entry points in the table above, including the Java
  `SearchEngine` and `SMSD` classes.
- The near-MCS, post-filter and progress-callback options, including Python
  `find_mcs_progressive`; SMARTS queries use `Substructure`.
- The C++ headers `smsd/vf2pp.hpp` and `smsd/clique_solver.hpp`;
  `smsd/substructure.hpp` and `smsd/mcs.hpp` declare the searches, and
  `smsd/smsd.hpp` includes both.
- The undocumented in-class C++ SMILES writer members of `MolGraph`; use
  `smsd::toSMILES`.
- The private native entry points `depict_smiles`, `depict_mcs_svg`,
  `batch_fingerprint_screen`, `overlapCoefficient`, `count_overlapCoefficient`
  and `prewarm_graph` of `smsd._smsd`; the package-level functions remain.
- The superseded benchmark drivers, the `diverse_molecules.txt` pool and the
  SIDER collection `benchmarks/data/sider.csv`, which is distributed under a
  non-commercial Creative Commons licence. `chembl_mcs_benchmark.smi` and the
  Dalke-style pair files are rebuilt without SIDER structures.
