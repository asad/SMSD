# Changelog

All notable changes to SMSD Pro are documented in this file. The format is
based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

## [8.0.0] - 2026-10-21

Substructure and maximum common substructure (MCS) search for C++, Java and
Python through prepared `MCS` and `Substructure` queries, with one search
model, one set of chemistry rules and the same fingerprints in the three
languages. The [release notes](docs/RELEASE_NOTES.md) list the downloads and
the upgrade notes.

### Added

- Prepared `MCS` and `Substructure` queries with `match`, `matchAll` and `matchMany` (Python `match`, `match_all` and `match_many`), and `Substructure.smarts` for SMARTS queries.
- MCS objectives: connected and disconnected maximum common induced subgraphs (MCIS), maximum common edge subgraphs (MCES) and atom-weighted MCIS, with bond mappings, fragment results and fragment limits.
- Constrained, consensus and R-group searches: `MCS.matchConstrained`, `MCS.consensus` and `Substructure.decomposeRGroups`.
- Results that report status (`Completed`, `TimedOut`, `Cancelled` or `LimitReached`), proved optimality and completeness separately from the matches.
- Search budgets: timeouts, node limits and cancellation.
- Mapping tools: `MCS.validateMapping`, `MCS.isMappingMaximal`, `MCS.mappingMatch` and `MCS.canonicalizeMapping`.
- Tautomer confidence for the selected bonds of a match: `MCS.computeTautomerConfidence` (Python `MCS.compute_tautomer_confidence`).
- `FingerprintOptions`, a configurable fingerprint specification shared by C++, Java and Python.
- Similarity bounds for screening: `similarityUpperBound` and `screenTargets` (Python `similarity_upper_bound` and `screen_targets`).
- Python: `parse_smiles(smiles, lenient=True)`, `MCS.validate_tautomer_consistency`, and the full set of `ChemOptions` setters, with copy support.
- Java: `MolGraph.toCanonicalSmiles(Map<Integer, Integer>)` and `MolGraph.Builder.hydrogenCounts(int[])`.
- Shared benchmark workloads for C++, Java and Python, with independent checks of every returned mapping.

### Search behaviour

- A search is reported optimal only when it runs to completion; a search stopped by a limit returns the best valid mapping it found.
- Shared golden fixtures pin the same search and chemistry results in C++, Java and Python, and the C++ and Java MCES searches visit the same nodes in the same order.
- Stereo-aware substructure searches require a specified, matching target centre or double bond for each specified query centre or double bond.
- Batch searches keep input order, and their results do not depend on the number of threads unless a search reaches its time limit.
- The Java command-line tool exits with status 0 for a match, 1 for no match and 2 for an error or an incomplete search.

### Chemistry

- SMARTS follows Daylight semantics in C++ and Java, including directional bonds, and SMARTS input has documented size limits.
- SMILES bracket atoms are checked against the OpenSMILES limits.
- One implicit-hydrogen valence table and one tetrahedral-centre rule are shared by C++ and Java.
- CIP descriptors follow the IUPAC 2013 Recommendations, including pseudoasymmetric centres (r and s), identically in C++ and Java.
- Ring perception has bounded work and stack use on very large ring systems, and C++ and Java count rings by the same rule.
- Tautomer classes do not depend on the atom order, and tautomer weights depend on pH through the Henderson–Hasselbalch relation.
- Murcko scaffolds keep exocyclic multiply bonded atoms and stereochemistry.
- MOL and SDF input and output keep tetrahedral and double-bond stereochemistry.

### Fingerprints

- `FingerprintOptions` selects circular (ECFP or FCFP), path, topological-torsion or structural-key fingerprints, including the 166 MACCS keys, and switches atom, bond, ring, ring-size, stereo and fragment terms, with an optional Murcko framework.
- The default configuration of each family reproduces the bits of the corresponding fingerprint function.
- Two empty fingerprints are identical sets: Tanimoto and overlap are 1.0 in every language.

### Packaging

- Java: the library `com.bioinceptionlabs:smsd`, with a lean set of CDK dependencies, and the command-line tool `com.bioinceptionlabs:smsd-cli`, built from one Maven project.
- C++: a header-only C++17 library whose public headers each compile on their own.
- Python: the `smsd` package on the native library, with optional RDKit, NumPy and Open Babel support loaded only when used.
- The CDK licence ships with every Java artefact that bundles CDK.
- The guides document the prepared-query API, and `scripts/check_examples.py` runs every example in them.
