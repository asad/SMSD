# SMSD — Algorithms & Design Whitepaper

**Version**: 3.0.0  	**Toolkit**: Chemistry Development Kit (CDK)  	**Language**: Pure Java

## Abstract
SMSD (CDK‑compat) provides exact substructure and maximum common substructure (MCS) search using a pragmatic combination of:
- VF2‑style backtracking with terminality/degree look‑ahead for substructure;
- Modular‑product construction and maximum clique (Bron–Kerbosch/BBMC lineage) for MCS;
- Optional McGregor extension to grow clique seeds into chemically larger common subgraphs;
- Advanced recursive SMARTS predicates ($name and $(SMARTS)) anchored to candidate atoms;
- Chemically meaningful invariants (ring size/system, aromaticity, charge, valence, donor/acceptor, stereo).

The implementation is CDK‑only. Prior work on SMSD and related algorithms is discussed in the references.

## 1. Background & Prior Work
- The original **SMSD toolkit** (J. Cheminf., 2009) established a multi‑strategy approach to MCS, integrating VF family methods with CDK routines and heuristics. citeturn0search0turn0search5turn0search4
- **CDK** provides the core chemical model, ring analysis, SMILES/SMARTS parsing and substructure checks used by our code. citeturn0search1turn0search19
- **Maximum clique** algorithms of the **Bron–Kerbosch** family and **Tomita** improvements underpin many exact MCS approaches. citeturn0search2turn0search20turn0search16
- **McGregor (1982)** proposed an extension scheme that incrementally augments a seed mapping; we employ a similar bounded augmentation after clique finding. citeturn0search3turn0search10
- Reviews of **MCS algorithms** in cheminformatics provide useful context and variants. citeturn0search21

## 2. Standardisation & Pre‑processing
- Inputs: SMILES, SMARTS or `IAtomContainer`.
- Standardisation stages (configurable): largest fragment, normalisation/re‑ionisation, uncharge, optional kekulisation; **explicit hydrogens are removed** for matching. Ring info and ring system IDs are cached once.
- Rationale: canonicalise inputs and expose ring/system features to early filters; improve determinism.

## 3. Chemistry Constraints
A single options object (`ChemOptions`) governs:
- **Atoms**: element (wildcards allowed), isotope, formal charge, valence windows, aromaticity mode (strict/flexible), ring‑only matching, ring‑system and bridgehead enforcement, donor/acceptor compatibility, optional tetrahedral chirality.
- **Bonds**: strict vs loose order (aromatic compatibility), ring‑size guard (exact/subset with tolerance), double‑bond stereo (off/defined/exact).
- Constraints are checked **incrementally** during search.

## 4. Recursive SMARTS
`SmartPredicateProcessor` evaluates both `$name` and `$(SMARTS)`:
- Named predicates `$name` are resolved via a registry (either mapped to SMARTS or Java functions).
- `$(SMARTS)` blocks are **anchored**: the first atom in the SMARTS must map to the candidate atom—implemented via a CDK substructure check on the target with the anchor index pre‑bound.
- This enables domain‑specific chemistry predicates without leaving CDK.

## 5. Substructure (VF2‑style) Search
State keeps partial maps `q→t` and `t→q`, terminal sets, degrees, and a pre‑filtered compatibility matrix.
Expansion:
1. Pick the next unmapped query atom with **smallest candidate set**, favouring high degree and ring membership for tie‑breaks.
2. For each candidate target atom: enforce atom compatibility, neighbour bond compatibility to already mapped neighbours, terminality/degree look‑ahead, and optional **induced** constraint (non‑bond in `q` must map to non‑bond in `t`).
Success: on full depth, a mapping is reported; enumeration is supported with de‑duplication by mapping or by image. Time limits stop cleanly.

## 6. MCS via Modular Product + Maximum Clique
**Modular product**: each node is a compatible `(qi, tj)` atom pair. Two nodes connect if their **pairwise relation is consistent**:
- **MCCS (default)**: edges only when both `q(i,k)` and `t(j,l)` are bonds and **bond pair** is compatible.
- **MCIS (induced)**: edges present when both pairs are bonds and compatible **or both non‑bonds**.
**Maximum clique**: we use a Bron–Kerbosch/BBMC‑style backtracker with Tomita pivoting and a greedy colouring upper bound for pruning. citeturn0search2
**McGregor extension**: the resulting clique is treated as a seed mapping and **augmented** by a bounded backtracker to add further atom pairs while preserving feasibility. citeturn0search3
Connectivity: for MCCS we retain the largest connected component in the query image.

## 7. Multi‑way (N‑MCS)
For multiple targets, we employ a **sequential reduction** strategy (and optional best‑first ordering by size). This is simple, deterministic and scales to moderate N; it matches common practice for multi‑set MCS. citeturn0search7

## 8. Stereochemistry & Rings
- **Tetrahedral** chirality, when enabled, must agree (like‑for‑like tags).
- **Cis/trans** double bonds: “defined” requires both sides defined; “exact” requires equality.
- **Rings**: ring‑only matching, ring **size** guards on bonds, optional **complete‑ring** validation for final mappings; ring‑system and bridgehead heuristics available.

## 9. Output & CLI
- **JSON** payloads include: query/target atom indices, target bond indices, `(query,target)` pairs, optional subgraph SMILES, and Tanimoto/Jaccard on atoms and bonds.
- CLI supports jar‑first usage and generated launchers; `--induced` toggles MCIS, `-m` enumerates mappings, `--json` emits machine‑readable output.

## 10. Complexity & Performance
Both substructure and MCS are exponential in the worst case; practical performance is achieved through strong domain filters, early bond feasibility, tight next‑vertex selection, and bounding in the clique stage. Tomita‑style bounds are particularly effective. citeturn0search2

## 11. References
- Rahman et al., **Small Molecule Subgraph Detector (SMSD) toolkit**, J. Cheminf., 2009. citeturn0search0turn0search5
- Willighagen et al., **The Chemistry Development Kit (CDK) v2.0**, J. Cheminf., 2017. citeturn0search1
- Tomita et al., **Worst‑case complexity for generating all maximal cliques**, 2006. citeturn0search2
- McGregor, **Backtrack search algorithms and the maximal common subgraph problem**, 1982. citeturn0search3
- Raymond & Willett, **MCS algorithms for cheminformatics**, 2002. citeturn0search21
- CDK site and documentation. citeturn0search19
