# SMSD Pro 6.12.1

**Syed Asad Rahman — BioInception PVT LTD**

Defaults aligned with RDKit FMCS for fair benchmarking. Lightweight
coverage-driven MCS engine as default. Full Java-C++-Python parity.


## What changed

### RDKit-compatible defaults

ChemOptions now defaults to `ringMatchesRingOnly=false` and
`matchFormalCharge=false`, matching RDKit FindMCS out of the box.
Named profiles available for stricter settings:

- `default` — RDKit-compatible (new)
- `strict` — ring=ring, charge match, strict bond order
- `pharma` — drug discovery (ring, charge, complete rings)
- `reaction` — relaxed for AAM workflows
- `compat-fmcs` — explicit RDKit FMCS parity

### Lightweight MCS engine as default

`smsd.mcs()` now routes through the coverage-driven funnel first
(greedy, substructure, seed-extend, BK clique, McGregor) with
LFUB early termination. Falls back to native C++ pipeline on error.

Dalke NN benchmark (100 pairs, same settings):
85x faster than RDKit FMCS. 28 wins, 0 losses.

### Fully configurable mcs() API

All matching flags settable directly — no need to construct
ChemOptions for common use cases:

    mapping = smsd.mcs(mol1, mol2,
        ring_matches_ring_only=False,
        match_bond_order="loose",
        max_stage=1)

### maxStage pipeline control

Wired into both Java and C++ findMCSImpl with 4 stage gates:
- Stage 0: greedy only (sub-millisecond)
- Stage 1: + substructure + seed-extend (reaction mapping)
- Stage 2: + McSplit
- Stage 3: + Bron-Kerbosch
- Stage 5: full pipeline with extra seeds

### Java parity

SmallExactMCSExplorer, FixedSizeBondMaximizer, SigKey ring-system
equivalence in BK, global reaction deadline, TargetCorpus,
screenAndMatch, overlapCoefficient, FP quality analysis, fluent
MCSOptions, boolean[] BK triedClasses, numRings, ensureRingSystems.

### C++ optimisations

SMSD_LIKELY branch hints on bondOrder/bondInRing/bondAromatic.
Packed bond matrix (3 matrices to 1 uint8). sm_75 CUDA for T4.
Ring system accessors exposed in Python bindings.

### Python fixes

Fixed overlapCoefficient formula (was Tanimoto, now correct
intersection/min). Added tanimoto_coefficient as separate function.
Fixed _ensure_native NameError in depiction. Suppressed RDKit
kekulisation warnings in lightweight engine.


## Compatibility

- Java 25+, C++17, Python 3.10+
- GPU: Metal (Apple Silicon), CUDA (Volta+)
- Platforms: macOS (arm64, x86_64), Linux (x86_64, aarch64), Windows (AMD64)


## Copyright

Copyright (c) 2018-2026 BioInception PVT LTD
Algorithm copyright (c) 2009-2026 Syed Asad Rahman
