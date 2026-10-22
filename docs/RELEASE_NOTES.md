# SMSD Pro 8.0.1

SMSD Pro 8.0.1 corrects parsing, aromaticity, CIP descriptors and
stereochemistry, gives fingerprint quality one definition in C++, Java and
Python, and screens targets on every core in Python. The prepared-query API is
unchanged. The [changelog](https://github.com/asad/SMSD/blob/master/CHANGELOG.md)
lists the changes in full.

## Highlights

- SMILES, SMARTS and V3000 MOL input outside the supported characters or
  ranges is rejected with the same error on every platform, and SMILES atoms
  joined by `:` ring bonds are aromatic, as RDKit reads them.
- CIP R/S and E/Z labels that only the later sequence rules decide, and C++
  SMARTS chirality, ring closures and MCES stereochemistry that match Java.
- Python `screen_targets` takes `num_threads` and a prepared `TargetCorpus`,
  and `rings_truncated` reports when a bounded ring search is cut short.

## Downloads

| File | Requirements |
| --- | --- |
| C++ headers: `smsd-cpp-8.0.1-headers.tar.gz` | A C++17 compiler |
| Java library: `smsd-8.0.1.jar`, with `smsd-8.0.1-sources.jar` and `smsd-8.0.1-javadoc.jar` | Java 8 or later and CDK 2.13 |
| Command-line tool: `smsd-cli-8.0.1-jar-with-dependencies.jar` and `smsd-8.0.1-cli.tar.gz` | Java 8 or later |
| Python source distribution: `smsd-8.0.1.tar.gz` | Python 3.9 or later, a C++17 compiler and CMake 3.18 or later |
| Python wheels for macOS arm64, Linux x86_64 and Windows x86_64 | CPython 3.14, as for 8.0.0 |
| `SHA256SUMS` | Lists the SHA-256 checksum of each file |

## Upgrading from 8.0.0

- `expected_collision_rate` (Python `analyze_fp_quality`, Java
  `analyzeFpQuality`) is now the estimated share of hashed features that fell
  on a bit already set, and Java `suggestFingerprintSize` returns larger
  sizes; review thresholds based on either.
- In C++ and Python, rings written aromatic that are not, such as
  cyclobutadiene, are no longer aromatic.
- A long fragment MCES search for one best mapping can return a different
  mapping of the same size.

# SMSD Pro 8.0.0

SMSD Pro 8.0.0 provides substructure and maximum common substructure (MCS)
search for C++, Java and Python through prepared `MCS` and `Substructure`
queries. A query is prepared once and searches one or more targets; each
result holds its matches with their atom and bond mappings, the search status
and whether the optimum was proved. The
[changelog](https://github.com/asad/SMSD/blob/master/CHANGELOG.md) lists the
changes in full.

## Highlights

- One search API in three languages: connected and disconnected maximum
  common induced subgraphs (MCIS), maximum common edge subgraphs (MCES),
  atom-weighted MCIS and substructure search, through `match`, `matchAll` and
  `matchMany` (Python `match`, `match_all` and `match_many`).
- Constrained, consensus and R-group searches: `MCS.matchConstrained`,
  `MCS.consensus` and `Substructure.decomposeRGroups` (Python
  `match_constrained`, `consensus` and `decompose_r_groups`).
- Results that report status, proved optimality and completeness separately
  from the matches, under timeouts, node limits and cancellation.
- Daylight SMARTS semantics in C++ and Java through `Substructure.smarts`,
  documented SMARTS size limits in every language, and the OpenSMILES limits
  for SMILES bracket atoms in C++ and Python.
- The same valence and tetrahedral-centre rules in C++, Java and Python, and
  CIP descriptors for R/S that follow the IUPAC 2013 Sequence Rules 1 to 5,
  with Rule 6 as proposed by Hanson et al. (2018).
- Tautomer classes that do not depend on the atom order, pH-dependent tautomer
  weights and a tautomer confidence for the selected bonds of a match.
- `FingerprintOptions`, one fingerprint specification for circular, path,
  torsion and structural-key fingerprints, including the 166 MACCS keys, with
  identical bits in C++, Java and Python.
- The Java library `com.bioinceptionlabs:smsd`, which depends only on the CDK
  modules it uses, and the separate command-line tool
  `com.bioinceptionlabs:smsd-cli`.

## Downloads

| File | Requirements |
| --- | --- |
| C++ headers: `smsd-cpp-8.0.0-headers.tar.gz` | A C++17 compiler |
| Java library: `smsd-8.0.0.jar`, with `smsd-8.0.0-sources.jar` and `smsd-8.0.0-javadoc.jar` | Java 8 or later (Java 25 LTS preferred) and CDK 2.13 |
| Self-contained CLI: `smsd-cli-8.0.0-jar-with-dependencies.jar` | Java 8 or later |
| Portable CLI with Unix and Windows launchers: `smsd-8.0.0-cli.tar.gz` | Java 8 or later |
| Python source distribution: `smsd-8.0.0.tar.gz` | Python 3.9 or later, a C++17 compiler and CMake 3.18 or later |
| Python wheel for macOS arm64 | macOS 26 or later and CPython 3.14 |
| Python wheels for Linux x86_64 and Windows x86_64 | glibc 2.28 or later, or Windows 10 or later, and CPython 3.14 |
| `SHA256SUMS` | Lists the SHA-256 checksum of each file |

The wheels run batch operations on the CPU with OpenMP; the optional Metal and
CUDA screening backends need a source build. The
[installation guide](https://github.com/asad/SMSD/blob/master/docs/HOWTO-INSTALL.md)
explains how to verify the downloads and install each package, and the
[validation record](https://github.com/asad/SMSD/blob/master/docs/VALIDATION_8.0.0.md)
lists the release checks.

## Upgrading from 7.x

The 7.x search functions are removed. Prepare a query with `find` (or
`Substructure.smarts` for SMARTS), pass search options with `withOptions`
(Python `with_options`), whose `chemistry` field holds the chemistry
settings, then search with `match`, `matchAll` or `matchMany`. Java
`matchMany` takes a `List<MolGraph>`; `new MolGraph(container)` converts a CDK
molecule.

| 7.x | 8.0.0 |
| --- | --- |
| Python `smsd.find_mcs(a, b)`, `smsd.mcs_result(a, b)`, `smsd.mcs_rdkit(a, b)` | `smsd.MCS.find(a).match(b)`, which accepts SMILES, `MolGraph` objects and RDKit molecules |
| Python `smsd.find_substructure(q, t)`, `smsd.is_substructure(q, t)` | `smsd.Substructure.find(q).match(t)`, then `.found` or `.best` |
| Python `smsd.batch_mcs(q, targets)`, `smsd.batch_substructure(q, targets)` | `smsd.MCS.find(q).match_many(targets)`, `smsd.Substructure.find(q).match_many(targets)` |
| Python `smsd.smarts_match(smarts, t)`, `smsd.smarts_find_all(smarts, t)`, `smsd.compile_smarts(smarts)` | `smsd.Substructure.smarts(smarts)`, then `.match(t).found` or `.match_all(t, limit)` |
| Java `SearchEngine.findMCS(a, b)`, `new SMSD(q, t, chem).findMCS()` | `MCS.find(a).match(b)` |
| Java `SearchEngine.findSubstructure(q, t)`, `SearchEngine.isSubstructure(q, t, chem, ms)`, `new SMSD(q, t, chem).isSubstructure()` | `Substructure.find(q).match(t)`, then `found()` |
| Java `SearchEngine.batchSubstructure(q, targets, chem, ms)`, `SearchEngine.findAllMCS`, `SearchEngine.findAllSubstructures` | `matchMany(targets, numThreads)` or `matchAll(target, limit)` on the prepared query |
| Java `SearchEngine.findMCSSmarts(smarts, target, chem, ms)`, `new SMSD(smarts, target, chem)` | `Substructure.smarts(pattern).match(target)` |
| C++ `smsd::findMCS(a, b, ...)`, `smsd::isSubstructure(q, t, ...)`, `smsd::findSubstructure(q, t, ...)` | `smsd::MCS::find(a).match(b)`, `smsd::Substructure::find(q).match(t)` |
| C++ `smsd::batch::batchMCS`, `smsd::batch::batchSubstructure`, `smsd::findAllMCS`, `smsd::findAllSubstructures` | `matchMany(targets, numThreads)` or `matchAll(target, limit)` on the prepared query |
| Python `timeout_ms=` keyword arguments, Java `SMSD.setMCSTimeoutMs(ms)` and `setSubstructureTimeoutMs(ms)` | `withTimeout(ms)` (Python `with_timeout(ms)`) on the prepared query |
| Java artefact `com.bioinceptionlabs:smsd` with the CLI | `com.bioinceptionlabs:smsd` for the library and `com.bioinceptionlabs:smsd-cli` for the command-line tool |

- **Status and optimality:** each result reports `status` (`Completed`,
  `TimedOut`, `Cancelled` or `LimitReached`; in Python `SearchStatus.COMPLETED`,
  `TIMED_OUT`, `CANCELLED` or `LIMIT_REACHED`), `optimal` and `complete`. A
  stopped search may return a valid mapping without proving the optimum, and a
  capped `matchAll` (Python `match_all`) may prove it without listing every
  match, so check `optimal` and `complete` separately.
- **Java library and CLI:** the library no longer brings picocli,
  `jackson-databind` or CDK modules such as `cdk-smiles`, `cdk-data` and the
  CDK file readers; projects that used them through `smsd` declare them. The
  self-contained CLI is `smsd-cli-8.0.0-jar-with-dependencies.jar`.
- **SMARTS:** in C++, as in Java, uppercase bracket atoms such as `[C]` match
  aliphatic atoms only and lowercase ones aromatic atoms only; use `[#6]` or
  `[C,c]` for both. C++ reads `[H]`, `[2H]` and `[H+]` as hydrogen atoms, and
  SMARTS beyond the size limits is rejected in every language.
- **Deprecated names** keep working in 8.0.0, and Python names warn with
  `DeprecationWarning`. The warning, or the Java and C++ API documentation of
  the name, gives its replacement where one exists.
