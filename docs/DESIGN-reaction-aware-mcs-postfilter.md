# Design Document: Reaction-Aware MCS Post-Filter

**Author:** Syed Asad Rahman
**Date:** 2026-03-30
**Status:** IMPLEMENTED (v6.5.0+) -- archived design document
**Implemented in:** v6.5.0 (ReactionAwareScorer), v6.6.0 (BondChangeScorer)

---

## 1. Problem Statement

SMSD's MCS algorithm is mathematically correct: it returns the largest common
subgraph.  For general cheminformatics (similarity, scaffold extraction) this is
the right answer.  For **reaction atom-atom mapping**, however, the chemically
most relevant mapping may be slightly smaller than the mathematical maximum.

**Motivating example -- SAM methyltransferase:**

| Candidate | Size | Heteroatom types covered |
|-----------|------|--------------------------|
| A         | 7    | {C, N, O}                |
| B         | 6    | {C, N, O, **S**}         |

Candidate A is mathematically larger, but Candidate B captures the sulfur atom
that participates in the methyl-transfer reaction center.  A reaction-mapping
tool should prefer B.

### Design Constraint

The core MCS algorithm (`SearchEngine.findMCS`) is **immutable**.  It stays
pure and mathematically correct.  The reaction-aware logic is a separate
**post-filter** that re-ranks candidates produced by the existing enumeration
pipeline.

---

## 2. Architecture Overview

```
                 findAllMCS(g1, g2, ..., maxResults)
                          |
          +---------------+---------------+
          |  existing pipeline (unchanged) |
          |  returns List<Map<Int,Int>>    |
          |  all of size K (maximum)       |
          +-------------------------------+
                          |
                          v
           +-----------------------------+
           | NEW: findNearMCS(g1, g2, K) |
           | generates candidates of     |
           | sizes K, K-1, K-2           |
           +-----------------------------+
                          |
                          v
          +-------------------------------+
          | NEW: ReactionAwareScorer      |
          | scores each candidate         |
          | returns ranked list           |
          +-------------------------------+
                          |
                          v
               best-scoring candidate
```

**Key insight:** The current `findAllMCS` only returns size-K candidates (all
maximum-size).  For reaction relevance we also need near-miss candidates of size
K-1 and K-2 that may cover more heteroatom types.  This requires a new
`findNearMCS` helper (Section 5).

---

## 3. Scoring Function

### 3.1 Composite Score

Each candidate mapping `m` (a `Map<Integer, Integer>` from g1 indices to g2
indices) receives a composite score:

```
score(m) = w_size * S_size(m)
         + w_hetero * S_hetero(m)
         + w_rare * S_rare(m)
         + w_conn * S_conn(m)
```

All sub-scores are normalised to [0, 1].  Default weights:

| Weight     | Default | Rationale |
|------------|---------|-----------|
| `w_size`   | 0.40    | Size still matters -- avoid collapsing to trivial mappings |
| `w_hetero` | 0.35    | Heteroatom coverage is the primary reaction signal |
| `w_rare`   | 0.15    | Rare elements (S, P, Se, B) are disproportionately important |
| `w_conn`   | 0.10    | Prefer connected mappings when scores are close |

### 3.2 Sub-Score Definitions

**S_size (size preservation):**

```
S_size(m) = |m| / K
```

where K is the maximum MCS size.  A candidate of size K scores 1.0; size K-1
scores (K-1)/K, etc.  This penalises smaller candidates proportionally but
allows them to win if heteroatom coverage is sufficiently better.

**S_hetero (heteroatom type coverage):**

Let `H_query` = set of non-carbon element symbols present in g1.
Let `H_target` = set of non-carbon element symbols present in g2.
Let `H_universe` = `H_query` intersection `H_target` (heteroatoms present in both
molecules -- only these are mappable).
Let `H_mapped(m)` = set of non-carbon element symbols that appear among the
mapped atoms in candidate m.

```
S_hetero(m) = |H_mapped(m)| / |H_universe|      (if |H_universe| > 0)
            = 1.0                                 (if |H_universe| == 0)
```

When both molecules are pure hydrocarbons, all candidates score equally on
heteroatoms, and size dominates.

**S_rare (rare-element bonus):**

Some elements are disproportionately important as reaction centers.  The scoring
assigns rarity tiers:

| Tier | Elements         | Tier weight |
|------|------------------|-------------|
| 3    | S, P, Se, B, Si  | 3.0         |
| 2    | N, O, F, Cl, Br  | 1.5         |
| 1    | all others (non-C)| 1.0         |

```
R_universe = sum of tier_weight(e) for e in H_universe
R_mapped(m) = sum of tier_weight(e) for e in H_mapped(m)
S_rare(m) = R_mapped(m) / R_universe       (if R_universe > 0)
          = 1.0                              (if R_universe == 0)
```

This makes capturing a sulfur atom worth twice as much as capturing an extra
nitrogen, reflecting sulfur's rarity and its prevalence as a reaction center
(thiol-ester, methyltransferase, disulfide).

**S_conn (connectivity):**

```
S_conn(m) = 1.0   if mapped atoms in g1 form a single connected component
          = 1 / (number of connected components)   otherwise
```

Reactions operate on contiguous molecular fragments; disconnected mappings
usually indicate artifacts.

### 3.3 Tie-Breaking

When two candidates have identical composite scores (within epsilon = 1e-9):

1. Prefer the candidate with more mapped bonds (not just atoms).
2. Prefer the candidate whose mapped subgraph has lower graph edit distance
   to the original g1 substructure (preserves topology).
3. If still tied, prefer the candidate discovered first (stable sort).

### 3.4 SAM Methyltransferase Walkthrough

Given: g1 = SAM cofactor, g2 = SAH product.  `H_universe` = {N, O, S}.

| Candidate | |m| | H_mapped  | S_size | S_hetero | S_rare     | S_conn | **Total** |
|-----------|----|-----------|--------|----------|------------|--------|-----------|
| A (7-no-S)| 7  | {N, O}    | 1.00   | 0.667    | 0.462*     | 1.0    | **0.757** |
| B (6+S)   | 6  | {N, O, S} | 0.857  | 1.000    | 1.000      | 1.0    | **0.943** |

*S_rare(A): (1.5+1.5)/(1.5+1.5+3.0) = 3.0/6.0 = 0.500

Candidate B wins decisively despite being one atom smaller.

---

## 4. API Design

### 4.1 Approach: Callback on McsOptions (Option C from requirements)

After analysing the three options, **option (c) -- a callback** is the most
extensible.  However, a bare callback imposes too much burden on users.  The
design combines (b) and (c):

- A boolean flag `McsOptions.reactionAware` for the common case.
- A pluggable `McsPostFilter` interface for advanced users.
- The flag simply installs the built-in `ReactionAwareScorer`.

### 4.2 New Types

```java
// In SearchEngine.java (or a new file ReactionAwareScorer.java)

/**
 * Post-filter that re-ranks MCS candidates for reaction relevance.
 * Implementations receive the full candidate list and both molecule
 * graphs, and must return the candidates in preferred order.
 */
@FunctionalInterface
public interface McsPostFilter {
    /**
     * @param candidates  unmodifiable list of MCS candidates (size K down to K-delta)
     * @param g1          first molecule graph
     * @param g2          second molecule graph
     * @return            candidates sorted best-first; may be a subset
     */
    List<Map<Integer, Integer>> rank(
        List<Map<Integer, Integer>> candidates, MolGraph g1, MolGraph g2);
}
```

### 4.3 McsOptions Additions

```java
public static final class McsOptions {
    // ... existing fields unchanged ...

    /**
     * Enable built-in reaction-aware post-filtering.  Off by default.
     * When true, findMCS and mapReaction will enumerate near-MCS
     * candidates (sizes K, K-1, K-2) and re-rank by heteroatom
     * coverage and reaction-center proximity.
     */
    public boolean reactionAware = false;

    /**
     * Maximum size deficit from the mathematical maximum K to consider.
     * Only used when reactionAware=true or postFilter!=null.
     * Default: 2 (consider candidates of size K, K-1, K-2).
     */
    public int nearMcsDelta = 2;

    /**
     * Maximum number of near-MCS candidates to generate before scoring.
     * Default: 20.  Higher values improve result quality at the cost
     * of enumeration time.
     */
    public int nearMcsCandidates = 20;

    /**
     * Custom post-filter.  If non-null, overrides the built-in
     * reaction-aware scorer even when reactionAware=true.
     */
    public McsPostFilter postFilter = null;
}
```

### 4.4 SMSD.java Convenience Method

```java
/**
 * Reaction atom-atom mapping with reaction-aware MCS post-filtering.
 * Enumerates near-MCS candidates and re-ranks by heteroatom coverage,
 * rare-element importance, and connectivity.
 */
public static Map<Integer, Integer> mapReactionAware(
    IAtomContainer reactants, IAtomContainer products,
    ChemOptions chem, long timeoutMs) { ... }
```

### 4.5 Python API Surface

```python
# New function
smsd.map_reaction_aware(reactants_smi, products_smi, timeout_ms=10000)

# Or via options
smsd.find_mcs(smi1, smi2, reaction_aware=True)

# Custom scorer (advanced)
smsd.find_mcs(smi1, smi2, post_filter=my_scorer_fn)
```

### 4.6 ChemOptions: No Changes

The `ChemOptions` class controls atom/bond matching constraints.  The
reaction-aware post-filter is orthogonal -- it re-ranks results after matching,
not during.  Therefore `ChemOptions` stays unchanged.

---

## 5. Near-MCS Candidate Generation

### 5.1 The Gap in findAllMCS

The current `findAllMCS` returns only candidates of the maximum size K.  For the
post-filter we need candidates of size K-1 and K-2 as well.  There are two
strategies:

**Strategy A: Deletion from K-sized mappings (preferred)**

For each K-sized mapping, systematically remove one mapped pair and collect the
resulting (K-1)-sized mappings.  Then for each (K-1)-sized mapping, check if a
different atom can be added that covers additional heteroatom types.  This is
efficient because it reuses existing K-sized solutions.

Steps:
1. Run `findAllMCS(g1, g2, C, M, N)` to get up to N size-K candidates.
2. For each size-K candidate, generate K single-deletion variants (size K-1).
3. For each deletion variant, attempt greedy re-extension with the constraint
   that the re-added atom must be a different element type than the deleted one.
4. Deduplicate by canonical SMILES of the mapped subgraph.
5. Optionally repeat step 2-4 on size-(K-1) candidates to get size-(K-2).

**Strategy B: Iterative deepening with reduced K**

Temporarily set a ceiling `maxMcsSize = K-1` and re-run the MCS pipeline.  This
is clean but expensive (runs the full MCS algorithm again).

**Decision: Strategy A** -- deletion and re-extension is O(K * N) additional
work on top of the already-computed candidates, vs. a full re-run.  It also
naturally generates candidates that differ from the maximum by exactly the atoms
that matter.

### 5.2 Candidate Budget

| Molecule size (atoms) | Recommended `nearMcsCandidates` | Rationale |
|-----------------------|---------------------------------|-----------|
| < 15                  | 10                              | Small molecules have few alternatives |
| 15 - 40               | 20 (default)                    | Typical drug-like molecules |
| > 40                  | 30                              | Larger molecules, more room for alternatives |

The budget is a soft cap.  If deletion generates fewer candidates than the
budget, no further work is done.  The time cost is bounded by the existing
`McsOptions.timeoutMillis`.

### 5.3 New Internal Method

```java
// Package-private, inside SearchEngine
static List<Map<Integer, Integer>> findNearMCS(
    MolGraph g1, MolGraph g2, ChemOptions C, McsOptions M,
    List<Map<Integer, Integer>> exactMCS, int delta, int maxCandidates)
```

This method is NOT part of the public API.  It is called internally when
`reactionAware=true` or `postFilter != null`.

---

## 6. Spectator Heteroatoms

### 6.1 Problem

Not every heteroatom is reaction-relevant.  A peripheral fluorine on a large
hydrocarbon scaffold, or a distant ether oxygen, may be a spectator.

### 6.2 Mitigation: Rarity Weighting Handles Most Cases

The S_rare scoring already addresses this partially.  A peripheral F has tier
weight 1.5, while S has weight 3.0.  So a mapping that captures F at the cost of
one carbon atom is less rewarded than one that captures S.

### 6.3 Optional Distance Heuristic (Future Enhancement)

For molecules where multiple heteroatoms compete, a proximity bonus could be
added:

```
S_proximity(m) = average(1 / (1 + dist(a, nearest_changed_bond)))
                 for each mapped heteroatom a
```

where `dist` is the graph-shortest-path to the nearest bond that differs between
g1 and g2.  This requires knowing which bonds change, which is only available
after an initial mapping.  This creates a chicken-and-egg problem.

**Deferred to v6.5.0.**  The current scoring (size + heteroatom coverage +
rarity) handles the SAM case and similar scenarios without needing proximity
information.

### 6.4 User Override via Custom Post-Filter

Users who need domain-specific logic (e.g., enzyme active site knowledge) can
supply a custom `McsPostFilter` that encodes arbitrary scoring.

---

## 7. Integration with mapReaction

The existing `SearchEngine.mapReaction` is a thin wrapper:

```java
public static Map<Integer, Integer> mapReaction(...) {
    McsOptions opts = new McsOptions();
    opts.disconnectedMCS = true;
    opts.connectedOnly = false;
    opts.timeoutMillis = timeoutMs;
    return findMCS(reactants, products, chem, opts);
}
```

The new `mapReactionAware` differs only in:

```java
opts.reactionAware = true;
opts.nearMcsDelta = 2;
opts.nearMcsCandidates = 20;
```

The existing `mapReaction` is unchanged.  Users opt in explicitly.

---

## 8. Performance Impact

### 8.1 When reactionAware = false (default)

**Zero impact.**  No new code paths are executed.  The post-filter is gated
behind `if (M.reactionAware || M.postFilter != null)`.

### 8.2 When reactionAware = true

Additional cost beyond standard findMCS:

| Phase | Cost |
|-------|------|
| findAllMCS (N=5 candidates) | ~2-3x single findMCS |
| Deletion variants (K * N) | O(K * N) subgraph extractions |
| Scoring (composite, all candidates) | O(candidates * K) -- negligible |
| **Total overhead** | **~3-5x single findMCS** |

For typical drug-sized molecules (20-50 atoms) with a 10-second timeout, the
additional time is 2-4 seconds.  Acceptable for reaction-mapping workflows which
are not latency-sensitive.

### 8.3 Memory

Each candidate is a `Map<Integer, Integer>` of at most K entries.  With 20
candidates and K=30, that is 600 map entries -- negligible.

---

## 9. Edge Cases

### 9.1 Pure Hydrocarbons

When `H_universe` is empty, `S_hetero` = 1.0 for all candidates.  The scorer
degrades gracefully to pure size ranking (same as standard MCS).

### 9.2 Single Candidate

If findAllMCS returns only one size-K candidate and no near-MCS candidates are
generated (all deletions produce duplicates), the single candidate is returned.
No performance penalty beyond the attempted enumeration.

### 9.3 All Candidates Score Identically

Falls through to tie-breaking: bond count, then topology, then insertion order.

### 9.4 K <= 2

For tiny MCS (1-2 atoms), the near-MCS delta is clamped:
`effectiveDelta = min(delta, K - 1)`.  This avoids generating empty mappings.

### 9.5 Symmetric Molecules

Symmetric molecules may produce many structurally identical mappings that
differ only in atom indices.  The existing SMILES-based deduplication in
`findAllMCS` handles this.  The near-MCS pipeline reuses the same deduplication.

### 9.6 Timeout Exhaustion

If the time budget is exhausted during near-MCS generation, the scorer works
with whatever candidates have been collected so far.  At minimum, the standard
size-K findMCS result is always available as a fallback.

---

## 10. Testing Strategy

### 10.1 Unit Tests

1. **SAM methyltransferase:** Verify 6-atom-with-S is ranked above 7-atom-no-S.
2. **Pure hydrocarbon pair:** Verify scorer returns size-K candidate (no regression).
3. **Single heteroatom (peripheral F):** Verify that a mapping capturing F but
   losing 2 carbons does NOT win over the size-K all-carbon mapping (S_size
   penalty outweighs S_hetero gain for a single low-tier heteroatom).
4. **Phosphoryl transfer (kinase):** Verify that P is captured.
5. **Disulfide bond formation:** Verify both S atoms are captured.
6. **Custom post-filter callback:** Verify the interface works with a lambda.
7. **Timeout edge case:** Verify graceful fallback to size-K under tight timeout.

### 10.2 Regression Tests

All existing MCS tests must pass unchanged when `reactionAware = false`
(the default).

### 10.3 Performance Benchmarks

Measure wall-clock time for `mapReaction` vs `mapReactionAware` on 10 known
reaction pairs of varying sizes (10-atom to 80-atom molecules).  The overhead
must stay under 5x.

---

## 11. File Layout (Proposed)

```
src/main/java/com/bioinception/smsd/core/
  McsPostFilter.java          # @FunctionalInterface (Section 4.2)
  ReactionAwareScorer.java    # Built-in implementation (Section 3)
  SearchEngine.java           # Modified: findNearMCS (pkg-private), gating logic
  SMSD.java                   # Modified: mapReactionAware convenience method
  ChemOptions.java            # UNCHANGED
  MolGraph.java               # UNCHANGED

src/test/java/com/bioinception/smsd/core/
  ReactionAwareScorerTest.java
```

---

## 12. Open Questions for Review

1. **Delta ceiling:** Is K-2 sufficient, or should reactions with large
   scaffold rearrangements allow K-3?  The default of 2 can be overridden via
   `nearMcsDelta`, so this is a question of defaults only.

2. **Weight tuning:** The default weights (0.40 / 0.35 / 0.15 / 0.10) are
   initial estimates.  Should we run a grid search over the SAM + kinase +
   disulfide test set to optimise them before release?

3. **Rarity tiers:** The current tiers {S,P,Se,B,Si}=3.0, {N,O,F,Cl,Br}=1.5
   are biochemistry-oriented.  Should we offer an alternative tier table for
   organic synthesis (where halogens are more important)?

4. **mapReaction default:** Should the existing `mapReaction` method turn on
   `reactionAware` by default in a future major version?  Or should it stay
   purely mathematical forever, with `mapReactionAware` as the opt-in?

---

## 13. Summary

| Decision               | Choice                                       |
|------------------------|----------------------------------------------|
| Core MCS algorithm     | UNCHANGED                                    |
| Post-filter trigger    | `McsOptions.reactionAware = true`            |
| Extensibility          | `McsPostFilter` callback interface           |
| Candidate generation   | Deletion from K-sized + greedy re-extension  |
| Scoring formula        | Weighted: size + heteroatom + rarity + conn  |
| Default candidate cap  | 20 candidates, delta = 2                     |
| Performance impact     | Zero when off; 3-5x findMCS when on          |
| New public API surface | 1 interface, 1 class, 1 convenience method   |
