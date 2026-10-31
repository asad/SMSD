# Changelog

All notable changes to SMSD Pro are documented in this file. The format is
based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

## [Unreleased]

### Added

- A task-based [user guide](docs/guide/README.md) in `docs/guide` and five
  topic notebooks in `python/examples`.

### Changed

- Clearer CUDA, conda, Python and SMARTS export guides: the 256 MiB limit
  applies to each retained device workspace, CUDA Graph replay needs the
  call's buffers to fit a free retained workspace, and identical fingerprint
  bits hold for the same SMSD version.

## [8.1.0] - 2026-10-31

Python data and machine-learning adapters, local MCP tools, packed fingerprint
corpora, corpus threshold searches and SMARTS export of MCS results, and a
correction to the C++ RDKit adapter. Deprecated 8.x APIs remain available and
supported throughout the 8.x line. The Java library and command-line tool stay
at 8.0.3 on Maven Central. The CPU Docker images carry the 8.1.0 command-line
tool; otherwise the Java changes below are available only from a source build.
The [release notes](docs/RELEASE_NOTES.md) list the downloads and upgrade notes.

### Added

- `relevantCyclesTruncated` (C++ and Java) and Python
  `relevant_cycles_truncated` report incomplete full relevant-cycle enumeration
  separately from the existing `ringsTruncated`/`rings_truncated` basis flag.
- `RingBasisUnavailable` identifies a result refused because ring perception
  reached its work limit before the minimum cycle basis was complete: C++
  `smsd::RingBasisUnavailable` (a `std::runtime_error`), Python
  `smsd.RingBasisUnavailable` (a `RuntimeError`) and Java
  `RingBasisUnavailableException` (an `IllegalStateException`), so existing
  handlers still catch it. C++ `MolGraph::ringCountsTruncated()` reports when
  the per-atom ring counts may be incomplete for that reason, and
  `MolGraph::computeRelevantCyclesFromBasis` finds relevant cycles from a given
  basis.
- Relevant-cycle families: C++ `computeRelevantCycleFamilies`
  (`relevant_cycle_families.hpp`), Python `compute_relevant_cycle_families` and
  Java `MolGraph.computeRelevantCycleFamilies()` group the relevant cycles into
  families, each with one member and its exact number of members, so the
  relevant cycles can be counted without listing them; `truncated` reports an
  incomplete enumeration. `relevantCycleFamilyOf` (Python
  `relevant_cycle_family_of`) finds the family of a ring, and
  `relevantCycleFamilyMembers` (Python `relevant_cycle_family_members`) lists
  members up to a given limit. When the enumeration is complete, the relevant
  cycles and their number do not depend on the input atom order; the grouping
  into families can.
- `Substructure.match_corpus` and `MCS.match_corpus` (C++ `matchCorpus`)
  decide, for each target in a corpus, whether the query is contained or shares
  at least a given number of atoms, with an optional GPU pre-screen
  (`gpu_prefilter`). Matching runs on the CPU; ordinary batch methods are
  unchanged.
- Optional Python feature adapters for NumPy, SciPy CSR, pandas, Polars,
  Arrow, PyTorch and PyTorch Geometric, plus scikit-learn fingerprint and
  similarity transformers. Arrow records buffer ownership and required copies.
- Local molecular tools for MCS, substructure, SMARTS, fingerprints, similarity,
  screening, standardisation and depiction; one registry for MCP, function
  definitions and optional LangChain/LlamaIndex adapters. Results retain limits
  and native proof/status where applicable.
- MCP container and Hugging Face demonstration source recipes using supplied
  wheels and offline installation.
- A consolidated API reference, an `llms.txt` documentation index and a guide
  to connecting AI assistants to the local tools.
- Python examples for NVIDIA screening checked against the CPU, an offline
  chain of the molecular tools and an MCP client.
- A C++ example, `cpp/examples/metal_fingerprint_screen.cpp`, screens packed
  fingerprints on a Metal GPU through `metal::fingerprintAccelerator()`, falls
  back to the CPU screen when the provider declines or fails, and reports which
  path ran; build it with `-DSMSD_BUILD_EXAMPLES=ON`.
- Native SVG highlights of existing MCS results.
- Memory-mapped packed binary fingerprints with streamed construction, CPU/CUDA
  threshold and exact top-k scans, Butina and MaxMin. Explicit backend selection
  and `cuda_calls` expose completed packed device work; clustering bookkeeping
  remains deterministic CPU work.
- CUDA kernel preload, retained workspace, opt-in CUDA Graph replay and on-device
  top-k selection.
- `SMSD_CUDA_MIN_TARGETS` lets automatic ECFP4 screening of a prepared
  `TargetCorpus` use any NVIDIA GPU with compute capability 7.5 or later from a
  chosen corpus size; on a Tesla T4 it replaces the 100,000-target default.
- Streamed native corpus screening and `warm_start` for an existing target corpus.
- `smsd-cuda` includes native code for GPUs of compute capability 10.0 and
  12.0, such as data-centre Blackwell GPUs, the GeForce RTX 50 series and
  RTX PRO Blackwell cards.
- `MCS::matchWithSmarts`, Java `MCS.matchWithSmarts` and Python
  `MCS.match_with_smarts`, returning the ordinary result beside a SMARTS pattern
  under the prepared chemistry, in Java as an `MCS.SmartsResult`. Java
  `SearchResult.Match.toSmarts(query, chemistry)`, with the query as a
  `MolGraph` or `IAtomContainer`, exports one match, and `SmartsWriter` exports
  molecular graphs.
- Conda recipes for `smsd` (CPU) and `smsd-cuda` (Linux CUDA), with conda,
  mamba and Pixi install routes.
- Release tooling updates citation metadata, extracts release notes and
  collects wheel manifests.

### Changed

- `MCS.consensus` remains an approximation: `optimal` is always false,
  including for empty and single-input calls, and `complete` means the
  computation finished, not that no larger qualifying common graph exists.
- Ring perception works within one fixed, deterministic work limit, the same
  in C++, Java and Python and independent of time, thread count and hardware.
  Very large fused-ring systems, such as synthetic cages, can reach it. Some
  ring systems that 8.0.3 left incomplete now complete, so their ring-dependent
  results can change. Where rings stay incomplete, searches that need them
  report `LimitReached` and other calls raise an error naming the cause, with
  the zero-based input index in C++ and Python batch calls, where 8.0.3
  computed such results from incomplete rings.
- Complete-ring matching uses all relevant cycles, so a relevant cycle outside
  the smallest set of smallest rings now counts as a complete ring; 8.0.3 could
  reject it and report a smaller result, even an empty one, as proved.
  `relevantCyclesTruncated` reports when the enumeration is incomplete. With
  `completeRingsOnly`, Java substructure matching now rejects embeddings that
  cover only part of a ring, as C++ and Python do, and reports `LimitReached`
  when the ring perception it needs is incomplete.
- Per-atom ring counts, and with them strict ring-fusion results, can differ
  from 8.0.3 and can still depend on the input atom order.
- Positive numeric SMARTS `R<n>` (`n > 0`) uses one deterministic minimum cycle
  basis, the same in C++ and Java, so a Java result that came from a different
  minimum basis of equal size can change. If ring perception cannot complete
  the basis, a search that cannot be decided without it returns `LimitReached`
  without matches, other functions that need the basis raise an error, and
  `computeSSSR` and `computeRings` may return incomplete rings with truncation
  flags. Literal `R` and `R0` do not depend on the basis, and
  `relevantCyclesTruncated` alone does not make it incomplete.
- Ring-size and double-bond stereo fingerprint terms, aromatic-ring keys (as in
  MACCS) and SMARTS keys whose search needs the minimum cycle basis raise an
  error naming the cause when ring perception cannot complete that basis, where
  8.0.3 could return values computed from incomplete rings. Java now uses the
  same basis for these terms as C++ and Python, so a Java value that came from
  a different minimum basis of equal size can change. Time and embedding limits
  on key searches behave as in 8.0.3, but in C++ and Python batch calls their
  errors now also give the zero-based input index, as ring refusals do.
- Compiled SMARTS evaluate distinct recursive predicates separately.
  Common-graph search rejects inconsistent bond properties.
- On Linux and macOS x86-64, fingerprint similarity and screening use the
  processor's population-count instruction when it is present, detected at run
  time; Windows builds already used it. With NumPy 2.0 or later, packed corpora
  count bits with NumPy's `bitwise_count`. Results are unchanged.
- C++ and Python MCIS search (the default MCS objective) has changed: when
  8.0.3 or 8.1.0 stops a search at its node or time limit, including the
  default time limit, the two results can differ, and a result of the same size
  can come with a different atom mapping. Java does not have this change, so a
  stopped search can return a different result in Java than in C++ and Python.
- macOS arm64 wheels now target macOS 14 or later; the 8.0.3 wheels targeted
  macOS 26.
- The Python wheels and conda recipes support CPython 3.12–3.15. CPython 3.10
  and 3.11 are no longer supported; on those versions, pip continues to install
  8.0.3.
- Python source builds, including installation from the source distribution,
  leave Metal and CUDA out by default; set `SMSD_BUILD_METAL` or
  `SMSD_BUILD_CUDA` to `ON` or `AUTO` with `-Ccmake.define` to include them.
  This changes only what `gpu_is_available()` and `gpu_device_info()` report;
  the published wheels are unchanged.

### Fixed

- The RDKit adapter (`fromRDKit` and `fromSmiles`, C++) now keeps each atom's
  hydrogen count and the E/Z configuration of double bonds where it is specified,
  bringing its graphs into line with SMSD's own SMILES parser. Previously,
  aromatic N–H atoms, such as those in pyrrole, indole and imidazole, were imported
  without their hydrogen.
- NOTICE and the packaged licence files now credit the RDKit for the MACCS key
  definitions and pybind11 for the Python extension module.
- The `smsd-cuda` package now ships the CUB and CCCL licences, and its NOTICE
  points to them.
- SMILES writers reuse closure labels and traverse iteratively, and the C++
  SMARTS writer, including `toSMARTS` (Python `to_smarts`), now reuses closure
  labels as Java `SmartsWriter` does. They reject graphs needing more than 99
  simultaneously open closures.
- SMARTS parsing in C++ and Python now rejects a pattern that ends with a bond,
  such as `CC-` or `CC=`, as Java already did; 8.0.3 C++ and Python accepted
  it and dropped the bond.
- Native corpus prewarming, native corpus screening on the CPU, the
  path-fingerprint form of C++ `batchFingerprint` and `batchFingerprintScreen`
  pass worker exceptions to the caller. On several threads the exception of the
  first input that fails is raised, as on one thread, instead of ending the
  process.
- A CUDA device below compute capability 7.5 is reported as unavailable, with
  the reason, rather than as available with failing device calls. Automatic
  screening then uses the CPU.
- Bond edits refresh cached topology and ring perception for later reads and
  searches.
- The MOL and SDF readers no longer record E/Z on double bonds that cannot carry
  it, such as those in rings of fewer than eight atoms; such files now give the
  same canonical SMILES as the equivalent SMILES and can be written as SMILES.
- `extractSubgraph` (Python `extract_subgraph`), and with it R-group
  decomposition, now keeps tetrahedral and E/Z stereo where the kept atoms
  determine it, as Java already did.

### Limitations

- Packed all-against-all scans take quadratic time, and keeping every result
  can need quadratic memory.
- Explicit packed CUDA failures propagate; native corpus screening retains its
  separate dispatch/fallback contract. Device availability alone does not prove
  that a native screen used the GPU.
- Portable SMARTS rejects unsupported stereo, tautomer, complete-ring and strict
  ring-fusion rules; it does not encode the induced MCIS objective or a proof.
- The PyTorch adapters reject detected bundled OpenMP conflicts; on macOS, an
  OpenMP-free SMSD source build may be needed.
- Canonical SMILES and canonical hashes of a molecule with more than 5,000
  atoms (200 in Java), or whose canonical labelling reaches its search limit,
  can depend on the input atom order; `automorphismGeneratorsTruncated()`
  (Python `automorphism_generators_truncated()`) is then true.
- In C++ and Python, canonical SMILES of molecules with very large rings, of
  several hundred ring atoms, can take minutes to compute.
- With complete-ring constraints, an MCIS search on a large polycyclic aromatic
  system or a cyclodextrin can, for some atom orders, reach its time limit
  without a match; its status reports the stop.
- In rare cases Java reports `Cancelled`, with `optimal()` true and `complete()`
  false, for a search whose proof had already completed; the score and the
  proof are unaffected.
- E/Z marks on a double bond with two identical substituents at one end are
  kept as written, so such a molecule written with and without the marks has
  two canonical SMILES.

## [8.0.3] - 2026-10-25

Corrections to complete-ring MCS matching, to the status of MCS results that
depend on truncated ring perception and to the built-in named SMARTS
predicates; opt-in Metal fingerprint screening for C++ source builds on macOS;
and Python `fingerprints` for many molecules in one call. The
[release notes](docs/RELEASE_NOTES.md) list the downloads and upgrade notes.

### Added

- C++: `smsd::metal::fingerprintAccelerator()` provides Metal fingerprint screening in C++ source builds on macOS. It keeps target fingerprints on the GPU between screens, returns the same hits as the CPU screen and runs only when passed explicitly; it is never selected automatically.
- Python: `fingerprints` computes the fingerprints of many molecules or SMILES strings in one call, on `num_threads` threads (0 uses the OpenMP default), each identical to `fingerprint` and in input order.

### Changed

- C++ and Java MCIS searches with stereochemistry find the largest mapping that keeps the configuration of mapped stereocentres and double bonds.
- C++ and Python: a single-result connected MCIS search without stereochemistry or ring options can return a different mapping of the same size when every atom of the target is bonded to every other.
- Java MCIS searches can use less memory on large molecules; results are unchanged.
- The guides describe each MCS status, the size limits and the default time budgets. They state that an MCES search beyond the size limits reports `LimitReached` without searching, and that stereochemistry is compared by the parity of mapped neighbours, not by CIP labels.
- MCS examples use N-methylacetamide against propanamide, where neither molecule contains the other.

### Fixed

- With complete rings or strict ring fusion, an MCS search that depends on truncated ring perception keeps its mappings but reports `LimitReached`, with `optimal` and `complete` false. `isMappingMaximal` (Python `is_mapping_maximal`) then raises its proof-unavailable error rather than returning a value.
- Complete-ring MCS matching recognises rings in graphs built without ring flags.
- C++ MCIS searches return the maximum mapping whatever the atom and bond order of the target.
- GPU similarity-bound screening returns the same hits as the CPU screen, including scores at the threshold and molecules with very large atom counts; the float-valued Metal similarity API is unchanged.
- C++: `batchFingerprint` accepts any sequence of graphs and, on several threads, throws the exception of the first molecule that fails, as on one thread.
- Built-in named SMARTS predicates, in C++ and Java:
  - `isHBAcceptor` matches N, O and S with a charge of 0 to −3. It no longer requires a charge, and it excludes positive atoms, amide N, pyrrole-type aromatic N, nitro and N-oxide N, and sulfoxide and sulfone S.
  - `isHBDonor` includes aromatic NH.
  - `isAmideN` no longer matches amines; it matches amide, imide and lactam N, aromatic ones included.
  - `isKetone` and `isEster` mark the carbonyl carbon, not the carbon next to it, and now include carbonyls in aromatic rings such as chromone and coumarin. `isEster` includes formates and excludes anhydrides.
  - `isPositive` and `isNegative` match charges of +1 to +15 and −1 to −15, not only +1 and −1.
- Java expands named predicates inside recursive SMARTS, as C++ does. As in C++, a name consists of ASCII letters, digits and underscores, so a registered name with other characters is no longer expanded.

## [8.0.2] - 2026-10-23

Optional GPU screening through the separate `smsd-cuda` package, ordered
batches of corpus screening queries, and CPU wheels for CPython 3.10 to 3.14 on
four platforms. The [release notes](docs/RELEASE_NOTES.md) list the downloads
and the upgrade notes.

### Added

- Optional GPU screening with the separate `smsd-cuda` package (Linux x86_64 with an NVIDIA GPU): corpus fingerprints stay on the GPU between screens, scores use the CPU expressions and hits keep target order. MCS and substructure matching stay on the CPU.
- `TargetCorpus.screen_many` and `TargetCorpus.screen_targets_many` screen ordered batches of queries (C++ `screenMany` and `screenTargetsMany`).
- Substructure searches over a `TargetCorpus` first skip targets whose element paths rule out a match, then match the rest on the CPU.
- CPU wheels for CPython 3.10–3.14 on macOS arm64, Linux x86_64 and aarch64, and Windows x86_64; source builds need Python 3.10 or later.
- CPU Docker images for Linux amd64 and arm64, and a GPU image for Linux amd64.

## [8.0.1] - 2026-10-22

Corrections to SMILES, SMARTS and MOL parsing, aromaticity, CIP descriptors
and stereochemistry, one definition of fingerprint quality in C++, Java and
Python, and target screening on every core in Python. The
[release notes](docs/RELEASE_NOTES.md) list the downloads and the upgrade
notes.

### Added

- C++: `smsd::batch::fingerprintQuality`, the fingerprint density and collision estimates that Python `analyze_fp_quality` reports.
- Python: `screen_targets` takes `num_threads` and, as the other batch calls do, runs on every core unless `OMP_NUM_THREADS` sets fewer; it also accepts a prepared `TargetCorpus` for repeated screens. The indices are unchanged. C++ `screenTargets` takes a thread count, and `TargetCorpus::screenTargets` screens a corpus.
- `ringsTruncated` (C++ `smsd::ringsTruncated(graph)`, Java `MolGraph.ringsTruncated()`, Python `rings_truncated(mol)`) and the C++ `MolGraph::computeRelevantCyclesReportingTruncation(truncated)` report when a bounded ring search is cut short; the rings are unchanged.

### Changed

- A long fragment MCES search for one best mapping, without stereochemistry, complete rings, fragment-size limits, a required mapping or excluded atoms, also searches from the other molecule, so that more pairs of molecules are proved within a node limit. Both directions count towards the node limit and the nodes visited, so a long search can return a different mapping of the same size, and one stopped at the limit a different best mapping; shorter searches return the same mappings as before.
- Fingerprint quality has one definition in C++, Java and Python: `expected_collision_rate` (Python `analyze_fp_quality`, Java `analyzeFpQuality`) is now the estimated share of hashed features that fell on a bit already set, assuming uniform hashing, so its values change in both languages, and Java reports saturation above a density of 0.75, as Python does.
- Java `suggestFingerprintSize` and the deprecated `analyzeFingerprintQuality` size a fingerprint for one hashed bit per feature, which gives larger sizes (about tenfold at a 1% rate), and `suggestFingerprintSize` rejects an atom count below 1, a target rate that is not strictly between 0 and 1 and a size that needs more bits than an `int` holds.
- Python: importing `smsd` raises `ImportError` when the native extension, which now reports its `__version__`, was built for another version.
- Python: `fingerprint` and `fingerprint_from_smiles` with `options`, `MCS.validate_mapping`, `MCS.is_mapping_maximal`, `read_mol_file`, `read_sdf` and `write_sdf` release the global interpreter lock while they work, as searches do.
- The CLI's `--version` reports the version recorded in its JAR manifest.
- Python source builds need scikit-build-core below 2 and pybind11 below 4, and the `dev` and `rdkit` extras have upper bounds.
- The Docker image is built from base images pinned by digest and carries the standard OCI labels.

### Fixed

- V3000 MOL blocks with an integer field that is not an integer in its range (a charge from −99 to 99, for example), a coordinate that is not a decimal number in the range of a double or more than 10,000 atoms raise `std::invalid_argument` (Python `ValueError`). MOL coordinates are read by one rule whatever the locale or platform, and one below 1e-307 in magnitude reads as zero.
- C++: SMILES output and path fingerprints take any formal charge without overflow, and MOL output rejects a hydrogen count that the valence field cannot hold.
- SMILES and SMARTS parsing rejects characters outside ASCII safely and with the same message on every platform; in Python the error is a `ValueError` that shows the byte as an escape.
- With MSVC, the C++ bit operations use the x64 intrinsics only on x64, the ARM64 intrinsics on ARM64 and portable code on other targets; defining `SMSD_PORTABLE_BITOPS` selects the portable code on any compiler.
- C++ MCES searches with bond stereochemistry treat a double bond in a three-membered ring as Java does, so they find the same mappings as Java and those mappings pass `validateMapping`.
- With complete rings, C++ treats an atom or bond on a relevant ring as a ring member even when its ring flag is clear, as Java does.
- CIP: ligands that atomic number alone leaves tied are ranked by Rules 1b, 2 and 3 read over the whole hierarchical digraph, so a ring closure nearer the stereocentre outranks an isotope or a double-bond configuration further out (R/S and E/Z).
- CIP: a double bond whose ends are told apart only by another double bond, as in 1,3-dialkylidenecyclobutanes, or only by like and unlike descriptor pairs (Rule 4), now gets its E/Z descriptor.
- CIP: in a mancude ring the duplicate atom of a ring double bond takes the mean atomic number of its ring neighbours whenever that differs from its partner's, so a carbon between boron and nitrogen no longer makes labels depend on the Kekulé structure.
- C++: `cip::assignRS` applies every sequence rule, as `cip::assignAll` does, so it also returns r and s.
- SMARTS: `@?` and `@@?` match an atom with that configuration or with none, as Daylight defines them; C++ no longer reads `@?` as either configuration or rejects `@@?`.
- SMARTS: in C++, a chiral query atom whose pattern neighbours do not define a configuration, as in `[C@]` or `F[C@H]Cl`, matches either configuration of a configured atom, as in Java and RDKit, rather than comparing the stored tags.
- SMARTS: C++ rejects a ring closure written with different bond expressions at its two ends, as Java does, and both reject one written with a bond at both ends and a direction (`/` or `\`) at either; Java, as C++, reads a direction written at the closing digit alone from the closing atom, so `F/C=C1.Cl/1` matches only the cis isomer.
- SMILES: an atom written in capitals but joined to its ring by `:` bonds is aromatic and takes the hydrogens of an aromatic atom, as RDKit reads it, so `C1:C:C:C:C:C:1` is benzene rather than six CH2 groups; an atom in brackets keeps the hydrogens written for it. In Java, `MolGraph` and `Standardiser.standardise` give the atoms CDK reads from such SMILES the same hydrogens, placing a hydrogen that the CDK molecule leaves open on a nitrogen that gives a Kekulé structure, whatever the atom order, and the CLI reads the hydrogens as written.
- SMILES output writes `/` and `\` on single bonds wherever a configured double bond has one at an end, so the configuration of a double bond in an aromatic ring is written on its substituent bonds rather than its ring bonds.
- Aromaticity (C++, Python): atoms written aromatic and joined by aromatic bonds are perceived again rather than kept, as RDKit and CDK read them, so a ring written aromatic that is not, such as cyclobutadiene or p-benzoquinone, is no longer aromatic. A ring carbanion gives two electrons only with three σ bonds, whether written with aromatic or Kekulé bonds. When perception cannot read every candidate cycle of a molecule, as with a large bowl-shaped polyarene, the aromatic systems supplied on its fused ring systems keep their flags, so that their aromaticity does not depend on the atom order.
- The CLI JAR with dependencies, the installers, the portable CLI archive and the Docker image carry the licence texts of the bundled Beam and vecmath libraries, and the image's licence label names every bundled licence.

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
