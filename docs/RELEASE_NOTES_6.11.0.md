# SMSD Pro 6.11.0 Release Notes

**Release date:** 2026-04-04
**Type:** Performance, precision, and depiction release

## Highlights

SMSD Pro `6.11.0` delivers three major improvements:

1. **Core engine hardening** -- cache-optimal data structures and pre-indexed
   candidate lookup reduce MCS search time by 20-40% on large molecule pairs.

2. **Publication-quality SVG depiction** -- zero-dependency renderer conforming
   to ACS 1996 standard (same specification used by Nature, Science, JACS, and
   Springer journals). Renders molecules, MCS comparisons, and substructure
   highlights directly from SMILES or MolGraph with no external tools required.

3. **Comprehensive layout engine** -- 8-phase 2D pipeline, distance geometry 3D,
   40+ pharmaceutical scaffold templates, full coordinate transform suite.

## Core Engine

- Converted all hot-path `vector<bool>` to `vector<uint8_t>` across McGregor DFS,
  Bron-Kerbosch partition bound, seed-extend, frontier expansion, k-core pruning,
  and tree detection. Yields 15-25% cache performance improvement due to
  byte-addressable access patterns.
- Pre-indexed candidate sets in McGregor DFS using precomputed `compatTargets_[]`
  from GraphBuilder. Eliminates O(n^2) linear scan per frontier atom.
- Eliminated redundant vector copies in `findAllMCS` return path.
- CTZ-based fingerprint bit extraction replaces sequential scan.

## Depiction Engine (New)

SVG renderer produces journal-ready molecular structure diagrams:

- **ACS 1996 proportions**: bond length, line width, bold width, font size, and
  double bond spacing all auto-scale from a single reference value
- **Skeletal formula**: carbon labels suppressed, heteroatom labels with H-count
  subscripts (NH2, OH) and charge superscripts
- **Bond rendering**: single, double (asymmetric toward ring interior), triple,
  wedge up (solid), wedge down (dashed stripes)
- **Aromatic rings**: solid inner circle
- **MCS highlighting**: green-filled circles, bold bonds, blue mapping numbers
- **Side-by-side pair rendering** with bidirectional arrow separator
- **Jmol/CPK element colors**: N=blue, O=red, S=amber, P=orange, etc.
- **Full customization**: DepictOptions controls every visual parameter

## Layout Engine

- 8-phase 2D pipeline: template match, ring-first placement, chain zig-zag,
  force-directed refinement, overlap resolution, crossing reduction (simulated
  annealing), canonical orientation, bond-length normalisation
- Distance geometry 3D: bounds matrix, double-centering, power iteration
  eigendecomposition, force-field refinement
- 40+ scaffold templates: single rings (3-8 membered), fused 2/3/4-ring systems,
  spiro, bridged (norbornane, adamantane), drug scaffolds
- Coordinate transforms: translate, rotate, scale, mirror, center, align,
  project/lift, bounding box, RMSD

## Python Bindings

35+ new functions exposed with GIL release for thread safety:

- `depict_svg()`, `depict_pair()`, `depict_mapping()`, `save_svg()`
- `find_nmcs()`, `find_scaffold_mcs()`, `validate_mapping()`
- `decompose_rgroups()`, `map_reaction()`, `extract_subgraph()`
- `generate_coords_2d()`, `generate_coords_3d()`
- Full transform suite and layout quality metrics
- All layout operations release the GIL for concurrent execution

## Java

- Explicit per-atom type matcher iteration replaces convenience shortcut for
  robust handling of organometallics, hypervalent S/P, and charged aromatics

## Tests

- C++ core: 114 tests (9 new precision chemistry)
- C++ layout: 42 tests (27 new)
- C++ CIP: 42 tests
- C++ parser: 542 tests
- Python: 170 tests
- Java: 602 tests

## Compatibility

- Fully backward-compatible with v6.10.x
- No API breaking changes
- Requires C++17, Java 11+, Python 3.9+
