# SMSD Pro 6.12.0

**Syed Asad Rahman — BioInception PVT LTD**

This release adds a lightweight clique-based MCS solver, standalone
fingerprint modules, and exposes the full set of search options to
Python users.


## What changed

### Clique solver

New header-only maximum clique finder (`clique_solver.hpp`) for
workflows where chemistry stays in Python/RDKit and only the
graph search runs in C++. Includes Bron-Kerbosch with Tomita
pivoting, k-core pruning, greedy seed, and McGregor bond-grow
extension. Portable across Mac/Linux/Windows (uses `bitops.hpp`
wrappers, no raw compiler builtins).

Python API: `find_mcs_clique`, `match_substructure`,
`match_substructure_from_elements`, `score_mapping`.

### Lightweight MCS engine

`smsd.mcs_engine` provides a coverage-driven funnel that accepts
SMILES strings, MolGraph objects, or RDKit Mol objects directly:

```python
from smsd.mcs_engine import find_mcs_lightweight
result = find_mcs_lightweight("c1ccc(O)cc1", "c1ccc(N)cc1")
```

Stages escalate automatically (greedy, substructure, seed-extend,
BK clique, McGregor) and stop as soon as the label-frequency upper
bound is reached. The C++ clique solver handles the heavy lifting.

### SmallExactMCSExplorer

Exact branch-and-bound MCS for small molecule pairs (up to 20×40
atoms). Now wired into the native `findMCS` pipeline for
disconnected MCS mode, giving deterministic results on hard
cases instead of relying on heuristic early exits.

### Fingerprint modules

Standalone `fp/` headers separated from `batch.hpp`:
- `fp/mol/circular.hpp` — ECFP / FCFP (Morgan)
- `fp/mol/path.hpp` — path-based fingerprints
- `fp/mol/pharmacophore.hpp` — pharmacophore features
- `fp/mol/torsion.hpp` — topological torsion
- `fp/mol/mcs_fp.hpp` — MCS-aware path fingerprints
- `fp/similarity.hpp` — Tanimoto, Dice, cosine, Soergel, subset
- `fp/common.hpp`, `fp/format.hpp` — shared utilities

### Support headers

- `hungarian.hpp` — O(n³) optimal assignment
- `periodic_table.hpp` — element data and valence tables
- `bond_energies.hpp` — bond dissociation energies
- `scaffold_library.hpp` — Murcko scaffold extraction
- `color_palette.hpp` — Jmol/CPK element colors
- `depictor.hpp` — publication-quality SVG rendering engine

### Full options exposure

All C++ search options now accessible from Python:

**MCSOptions** (18 fields): `max_stage`, `seed_neighborhood_radius`,
`seed_max_anchors`, `use_two_hop_nlf_in_extension`,
`use_three_hop_nlf_in_extension` added to the existing set.

**ChemOptions** (19 fields): `aromaticity_model`, `pH`,
`matcher_engine`, `induced`, `use_two_hop_nlf`, `use_three_hop_nlf`,
`use_bit_parallel_feasibility` added.

**Enums** (6): `AromaticityMode`, `AromaticityModel`, `BondOrderMode`,
`RingFusionMode`, `MatcherEngine`, `Solvent` — all accessible as
`smsd.MatcherEngine.VF2PP` etc.

`MolGraph.num_rings()` added for SSSR ring count.

### Global reaction deadline

`global_deadline` namespace in VF2++ enforces a single wall-clock
deadline across all MCS/substructure calls within a pipeline.
`TimeBudget` checks both local and global deadlines. Timeout
frequency increased (check every 256 calls instead of 1024).

### FixedSizeBondMaximizer

Maximizes bond count for fixed-size atom mappings on small pairs.
Integrated into `runValidatedMcsDirection` for automatic bond
refinement when MCS matches the smaller molecule entirely.

### Cross-platform

- MSVC: clique solver uses portable `smsd::popcount64` /
  `smsd::ctz64` (no raw `__builtin_*`). CMake adds `/W4 /EHsc`.
- ARM NEON: fused bitset ops dispatch automatically on Apple Silicon.
- CI wheels: Linux x86_64 + aarch64, macOS arm64 + x86_64,
  Windows AMD64. OpenMP installed in build containers. macOS
  wheels bundle libomp via delocate.

### Includes 6.11.2 fixes

Memory leak fix (WeakHashMap), BK color-bound overflow,
CIP Rule 3 (Z > E), thread-safety on lazy MolGraph fields,
tautomer weight corrections, Se/I support, SAH SMILES fix,
reaction-aware charge relaxation, Mcs→MCS rename,
tanimoto→overlapCoefficient.


## Compatibility

- Java 25+, C++17, Python 3.10+
- GPU: Metal (Apple Silicon), CUDA (Volta+)
- Platforms: macOS (arm64, x86_64), Linux (x86_64, aarch64), Windows (AMD64)


## Copyright

Copyright (c) 2018-2026 BioInception PVT LTD
Algorithm copyright (c) 2009-2026 Syed Asad Rahman
