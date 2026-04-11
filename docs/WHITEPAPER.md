# SMSD Pro 6.12.0 — Algorithms & Design

**Author**: Syed Asad Rahman, BioInception PVT LTD
**Version**: 6.12.0 | **Toolkit**: CDK 2.12 | **Language**: Java 25+ / C++17 / Python 3.9+ (recommended 3.12)

## Abstract

SMSD Pro provides exact substructure and maximum common substructure (MCS) search
for chemical graphs using:

- **VF2++ backtracking** with NLF pruning and bit-parallel domain filtering for substructure search
- **Modular product + maximum clique** (Bron-Kerbosch with RRSplit pivoting) for MCS
- **RASCAL screening** for fast similarity upper-bound filtering
- **FASTiso ordering** for accelerated candidate selection during backtracking
- **Greedy probing** to seed and bound clique search
- **McGregor extension** to grow clique seeds into larger common subgraphs
- **MolGraph precomputation** eliminating all CDK hot-path overhead
- **Canonical labeling** via Bliss-style individualization-refinement
- **Canonical SMILES** generation (deterministic, CDK-free)
- **Orbit-based automorphism** for atom equivalence classes
- Recursive SMARTS predicates with named registry
- Stereo-aware matching (E/Z double bonds, tetrahedral chirality)

## 1. Background

The original SMSD toolkit (Rahman et al., *J. Cheminformatics*, 2009)
established a multi-strategy approach to MCS by integrating VF-family graph
matching with CDK routines. Version 3.0 is a ground-up rewrite that retains
the algorithmic foundations while modernizing the implementation for
performance and correctness on CDK 2.11.

Key algorithmic references:

- **VF2/VF2++**: Cordella et al. (2004); Juttner & Madarasi (2018)
- **Maximum clique**: Bron & Kerbosch (1973); Tomita et al. (2006)
- **RRSplit pivoting**: Sun & Luo (2025)
- **RASCAL screening**: Raymond & Willett (2002)
- **FASTiso ordering**: Vu & Ahn (2025)
- **McGregor extension**: McGregor (1982)
- **CDK**: Willighagen et al., *J. Cheminformatics* (2017)

## 2. MolGraph Precomputation

All CDK object access is performed once at construction time. The `MolGraph`
class pre-builds:

| Data | Type | Purpose |
|------|------|---------|
| `atomicNum[]` | `int[]` | Atomic number per atom |
| `formalCharge[]` | `int[]` | Formal charge per atom |
| `label[]` | `int[]` | Encoded label: `(Z<<2)\|(aromatic<<1)\|ring` |
| `ring[]`, `aromatic[]` | `boolean[]` | Ring/aromaticity flags |
| `degree[]` | `int[]` | Vertex degree |
| `neighbors[][]` | `int[][]` | Adjacency lists (primitive arrays) |
| `bonds[][]` | `IBond[][]` | O(1) bond lookup matrix |
| `adj[]` | `BitSet[]` | Adjacency bitsets for pruning |

All subsequent operations (NLF building, feasibility checks, compatibility
tests, MCS seed construction) operate on these primitive arrays with zero
CDK method calls in inner loops.

## 3. Chemistry Constraints

`ChemOptions` governs matching behaviour:

- **Atom matching**: element type, formal charge, aromaticity mode
  (strict/flexible), ring-only constraint
- **Bond matching**: strict/loose/any bond order, ring-size guard,
  geometric (E/Z) stereo via `DoubleBondStereochemistry`
- **Chirality**: tetrahedral stereo via `ITetrahedralChirality`
- **Pruning**: 1/2/3-hop NLF (neighbor label frequency), bit-parallel
  adjacency feasibility

Constraints are checked incrementally during VF2++ backtracking.

## 4. Substructure Search (VF2++)

The `VF2PPMatcher` implements VF2++ with enhancements:

1. **Atom ordering**: BFS from highest-scoring atom (aromatic > ring > degree)
   with degree-sorted neighbor expansion
2. **Frontier pruning**: target candidates restricted to frontier set of
   already-matched neighbors
3. **Multi-level NLF**: 1-hop (direct neighbors), 2-hop, and 3-hop label
   frequency comparison eliminates incompatible candidates early
4. **Bit-parallel feasibility**: BitSet adjacency subset test
5. **O(1) bond lookup**: pre-built `IBond[][]` matrix
6. **Reusable scratch**: BitSets cleared and reused instead of allocated

The `VF2Matcher` (baseline) uses simpler neighbor-based candidate selection
without frontier tracking.

## 5. MCS via Modular Product + Maximum Clique

### Modular Product Construction

Each node represents a compatible `(query_atom, target_atom)` pair.
Edges connect nodes whose pairwise bond relationships are consistent:

- **MCCS** (connected): edges only when both pairs have compatible bonds
- **MCIS** (induced): edges also when neither pair has a bond

### Maximum Clique (Bron-Kerbosch with RRSplit)

The clique finder uses Bron-Kerbosch with RRSplit pivoting (Sun & Luo, 2025),
which applies a reduced-recursion split strategy to minimize redundant
branching. This is combined with Tomita-style pivoting and degree-ordered
start vertices. Multiple seed strategies run in parallel under time budget:

1. **Clique seed**: modular product maximum clique
2. **Ring anchor seed**: match ring atoms by aromatic/degree score
3. **Label-degree anchor seed**: match by rarest label, highest degree
4. **VF2++ ring skeleton seed**: substructure match on ring skeletons
5. **VF2++ core seed**: substructure match on heavy-atom core

### McGregor Extension

Each seed is extended by bounded BFS augmentation with NLF filtering.
Best result across all seeds is retained. Post-processing applies:

- Induced subgraph pruning (if requested)
- Largest connected component extraction
- Ring-anchor guard (prevents off-ring inflation)

## 6. Optimizations (v5.2 onwards)

### 6.1 Bit-Parallel Domain Filtering

Candidate domains for each query atom are represented as `BitSet` vectors.
Adjacency constraints are enforced via bitwise AND operations, pruning
impossible assignments before backtracking begins. This eliminates large
portions of the search space in O(n/w) time per constraint check, where
w is the machine word size.

### 6.2 Bron-Kerbosch with RRSplit Pivoting

The RRSplit strategy (Sun & Luo, SIGMOD 2025) partitions the candidate set
at each recursion level to reduce redundant subtree exploration. By splitting
candidates into reduced and recursive subsets, this approach achieves tighter
pruning than standard Tomita pivoting, particularly on dense modular product
graphs common in MCS computation.

### 6.3 RASCAL Screening

Before committing to full MCS computation, SMSD applies RASCAL-based
upper-bound estimation (Raymond & Willett, 2002). This computes a fast
similarity ceiling from the modular product graph structure. Pairs that
cannot exceed the current best similarity are skipped entirely, providing
significant speedups on batch MCS workflows. Exposed via the
`SearchEngine.similarityUpperBound()` API.

### 6.4 FASTiso Ordering

The FASTiso heuristic (Vu & Ahn, 2025) provides an improved candidate
ordering for subgraph isomorphism. It assigns priorities based on
neighborhood label diversity and connectivity patterns, producing orderings
that lead to earlier pruning and fewer backtracks compared to simple
degree-based orderings.

### 6.5 Greedy Probing

Before launching the full Bron-Kerbosch search, a greedy probing phase
constructs an initial clique by iteratively selecting the highest-degree
compatible vertex. This provides a non-trivial lower bound that prunes
unproductive branches early in the exact search, and serves as a fallback
result if the time budget expires.

### 6.6 Canonical Labeling (Bliss-Style Individualization-Refinement)

SMSD computes a canonical labeling for each `MolGraph` using a Bliss-style
individualization-refinement algorithm (McKay & Piperno, 2014). The procedure:

1. **Color refinement**: iteratively partition atoms by their neighborhood
   label multisets until the partition stabilises (equitable colouring).
2. **Individualization**: when the partition is not discrete, select the
   first non-singleton cell, individualize each atom in it, and recurse.
3. **Canonical leaf selection**: among all discrete partitions reached,
   select the lexicographically smallest as the canonical form.
4. **Automorphism pruning**: detect automorphisms during the search tree
   traversal and prune equivalent branches.

The result is a canonical relabeling `canonicalLabel[i]` giving the canonical
index for each atom `i`, and a `canonicalHash` (64-bit) derived from the
canonical adjacency encoding. Two molecules are isomorphic if and only if
their canonical hashes match, enabling **O(1) graph isomorphism** checks.

Exposed via `MolGraph.getCanonicalLabeling()` and `MolGraph.getCanonicalHash()`.

### 6.7 Canonical SMILES Generation

`MolGraph.toCanonicalSmiles()` produces a deterministic canonical SMILES
string directly from the MolGraph representation, with no CDK dependency.
The algorithm:

1. Applies the canonical labeling to determine atom visit order.
2. Performs a DFS traversal in canonical order, emitting atoms, branches
   (parentheses), and ring closures.
3. Handles disconnected components via dot-separated notation.
4. Encodes aromaticity (lowercase letters), charges, and hydrogen counts.

Because the traversal follows the canonical labeling, the output is identical
for all isomorphic input representations.

### 6.8 Orbit Computation and Automorphism Pruning

During canonical labeling, the algorithm accumulates automorphism generators.
Each generator permutes atom indices while preserving adjacency and labels.
The **orbit** of an atom is the set of all atoms reachable from it under the
automorphism group.

`MolGraph.getOrbits()` returns an array where `orbit[i]` is the smallest
atom index in the equivalence class of atom `i`. Atoms with the same orbit
value are symmetrically equivalent. This information is useful for:

- **Symmetry-aware enumeration**: avoid redundant atom mappings
- **Reaction site prediction**: identify equivalent functional groups
- **Canonicalization verification**: atoms in the same orbit are interchangeable

## 7. Public API

SMSD provides three API tiers, each with different trade-offs:

### Tier 1: SMSD facade (convenience, auto-standardises)
- `new SMSD(query, target, opts)` — standardises molecules automatically
- `new SMSD(smiles1, smiles2, opts)` — parses SMILES + standardises
- `smsd.isSubstructure()`, `smsd.findMCS()`, `smsd.similarityUpperBound()`

### Tier 2: SMSD facade (pre-standardised, skip preprocessing)
- `new SMSD(query, target, opts, false)` — molecules used as-is
- Same API as Tier 1, but 2-5x faster construction

### Tier 3: SearchEngine (zero overhead, full control)
- `SearchEngine.isSubstructure(q, t, opts, timeoutMs)` — substructure
- `SearchEngine.findMCS(q, t, opts, mcsOpts)` — MCS with full options
- `SearchEngine.similarityUpperBound(q, t, opts)` — RASCAL screening
- `SearchEngine.batchMCS(molecules, threshold, mcsOpts)` — batch MCS
- `SearchEngine.findDisconnectedMCS(q, t, opts, mcsOpts)` — disconnected MCS
- No standardisation, no object creation overhead

Standardisation (`Standardiser.standardise()`) is a CDK preprocessing step
independent of the matching algorithms. Users with their own preprocessing
pipelines should use Tier 2 or 3.

### Tier 4: MolGraph.Builder (CDK-free, toolkit-agnostic)
- `MolGraph.Builder` constructs molecules from raw primitive arrays
- No CDK dependency required at runtime
- Enables adapters for OpenBabel, RDKit, Schrödinger, or any toolkit
- `SearchEngine.isSubstructure(MolGraph, MolGraph, opts, timeout)`
- `SearchEngine.findMCS(MolGraph, MolGraph, opts, mcsOpts)`

### MolGraph Canonical Operations
- `MolGraph.getCanonicalLabeling()` — canonical atom permutation (Bliss-style)
- `MolGraph.getCanonicalHash()` — 64-bit canonical hash for O(1) isomorphism
- `MolGraph.toCanonicalSmiles()` — deterministic canonical SMILES (CDK-free)
- `MolGraph.getOrbits()` — automorphism orbit equivalence classes

Architecture: CDK is one adapter via `new MolGraph(IAtomContainer)`.
Users can write their own adapter using `MolGraph.Builder` for any toolkit.

## 8. SMARTS Support

- `Standardiser.matchAll(smarts, target)` matches SMARTS via CDK's
  SmartsPattern (direct API, no reflection)
- `SmartPredicateProcessor` expands named predicates (`$isKetone`,
  `$isAmideN`, etc.) from `PredicateRegistry` before matching
- Recursive `$(SMARTS)` blocks are anchored to candidate atoms

## 9. Performance

Benchmarks (best-of-500 substructure / best-of-100 MCS, JIT-warmed):

| Pair | Substructure | MCS |
|------|:-----------:|:---:|
| Benzene / Naphthalene | 97 us | 471 us |
| Aspirin / Acetaminophen | 90 us | 608 us |
| Morphine / Codeine | 153 us | 1,745 us |
| C15 / C20 chain | 104 us | 1,286 us |
| Triphenylene / Chrysene | 191 us | 2,496 us |

Both substructure and MCS are NP-hard in the worst case. Practical
performance is achieved through strong NLF pruning, bit-parallel
feasibility, RASCAL screening, FASTiso ordering, greedy probing,
RRSplit pivoting, frontier-based candidate selection, and time budgeting.

## 10. Test Suite

The project includes 420 JUnit 5 tests covering heterocycles, reactions,
drug pairs, stereo configurations, SMARTS predicates, batch workflows,
adversarial edge cases, large molecules (100-1000+ atoms), nucleotide MCS pairs,
and Kekule/aromatic equivalence across all supported molecule formats.

## 11. References

1. Rahman SA, Bashton M, Holliday GL, Schrader R, Thornton JM.
   *Small Molecule Subgraph Detector (SMSD) toolkit.*
   J. Cheminformatics, 1:12, 2009.

2. Willighagen EL, Mayfield JW, Alvarsson J, et al.
   *The Chemistry Development Kit (CDK) v2.0.*
   J. Cheminformatics, 9:33, 2017.

3. Cordella LP, Foggia P, Sansone C, Vento M.
   *A (sub)graph isomorphism algorithm for matching large graphs.*
   IEEE TPAMI, 26(10):1367-1372, 2004.

4. Tomita E, Tanaka A, Takahashi H.
   *The worst-case time complexity for generating all maximal cliques.*
   Theoretical Computer Science, 363(1):28-42, 2006.

5. McGregor JJ.
   *Backtrack search algorithms and the maximal common subgraph problem.*
   Software Practice and Experience, 12(1):23-34, 1982.

6. Raymond JW, Willett P.
   *Maximum common subgraph isomorphism algorithms for the matching
   of chemical structures.*
   J. Computer-Aided Molecular Design, 16(7):521-533, 2002.

7. Sun S, Luo S.
   *RRSplit: Efficient Maximum Clique Finding via Reduced-Recursion Splitting.*
   Proc. ACM SIGMOD, 2025.

8. Vu XT, Ahn HK.
   *FASTiso: A Faster Algorithm for Subgraph Isomorphism.*
   Proc. VLDB Endowment, 2025.

9. Sun S, Luo S.
   *SIGMo: Subgraph Isomorphism on GPU via Maximum Clique.*
   Proc. SC (Supercomputing), 2025.

10. McKay BD, Piperno A.
    *Practical Graph Isomorphism, II.*
    J. Symbolic Computation, 60:94-112, 2014.
    (Nauty/Bliss canonical labeling foundations)

11. Thakkar A, et al.
    *SynKit: A Comprehensive Toolkit for Synthesis Planning.*
    J. Chemical Information and Modeling (JCIM), 2025.

12. Thakkar A, et al.
    *Rxn-INSIGHT: Fast Chemical Reaction Analysis Using Bond-Electron Matrices.*
    J. Cheminformatics, 16:37, 2024.

13. Zhong Z, et al.
    *FlowER: Flow-based Enzymatic Retrosynthesis.*
    arXiv preprint, 2025.

14. Zhong Z, et al.
    *SynRXN: Synthesizability-Aware Retrosynthesis with Reaction Networks.*
    arXiv preprint, 2026.

15. McKay BD, Piperno A.
    *Nauty and Traces, version 2.9.3.*
    2026. https://pallini.di.uniroma1.it/
