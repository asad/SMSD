/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms. */
package com.bioinception.smsd.core;

import java.util.*;

/**
 * Built-in reaction-aware MCS post-filter that re-ranks candidates by
 * a composite score of size preservation, heteroatom coverage, rare-element
 * importance, and connectivity.
 *
 * <p>Scoring formula:
 * <pre>
 * score(m) = 0.25 * S_size + 0.40 * S_hetero + 0.25 * S_rare + 0.10 * S_conn
 * </pre>
 *
 * @author Syed Asad Rahman
 * @since 6.4.0
 * @see McsPostFilter
 */
public final class ReactionAwareScorer implements McsPostFilter {

  // Scoring weights tuned for reaction mapping: heteroatom coverage is critical.
  // A 6-atom mapping with S should beat a 7-atom mapping without S.
  private static final double W_SIZE   = 0.25;
  private static final double W_HETERO = 0.40;
  private static final double W_RARE   = 0.25;
  private static final double W_CONN   = 0.10;

  private static final double EPSILON = 1e-9;

  // Rarity tiers by atomic number
  // Tier 3 (weight 3.0): S(16), P(15), Se(34), B(5), Si(14)
  // Tier 2 (weight 1.5): N(7), O(8), F(9), Cl(17), Br(35)
  // Tier 1 (weight 1.0): all other non-carbon
  private static double rarityWeight(int atomicNum) {
    switch (atomicNum) {
      case 16: case 15: case 34: case 5: case 14: return 3.0;
      case 7: case 8: case 9: case 17: case 35:   return 1.5;
      default: return 1.0;
    }
  }

  @Override
  public List<Map<Integer, Integer>> rank(
      List<Map<Integer, Integer>> candidates, MolGraph g1, MolGraph g2) {
    if (candidates == null || candidates.isEmpty()) return candidates;

    // Determine K (max candidate size)
    int K = 0;
    for (Map<Integer, Integer> m : candidates) {
      if (m.size() > K) K = m.size();
    }
    if (K == 0) return candidates;

    // Build H_universe: heteroatom types present in both molecules
    Set<Integer> heteroG1 = heteroAtomTypes(g1);
    Set<Integer> heteroG2 = heteroAtomTypes(g2);
    Set<Integer> hUniverse = new HashSet<>(heteroG1);
    hUniverse.retainAll(heteroG2);

    // Compute R_universe (sum of rarity weights for the universe)
    double rUniverse = 0.0;
    for (int z : hUniverse) rUniverse += rarityWeight(z);

    // Score each candidate
    final int maxK = K;
    final double rUniv = rUniverse;
    List<ScoredCandidate> scored = new ArrayList<>(candidates.size());
    for (int idx = 0; idx < candidates.size(); idx++) {
      Map<Integer, Integer> m = candidates.get(idx);
      double sSize  = (double) m.size() / maxK;
      double sHetero = scoreHetero(m, g1, hUniverse);
      double sRare   = scoreRare(m, g1, hUniverse, rUniv);
      double sConn   = scoreConnectivity(m, g1);

      double total = W_SIZE * sSize + W_HETERO * sHetero + W_RARE * sRare + W_CONN * sConn;

      int bondCount = countMappedBonds(m, g1);
      scored.add(new ScoredCandidate(m, total, bondCount, idx));
    }

    // Sort best-first with tie-breaking
    scored.sort((a, b) -> {
      double diff = b.score - a.score;
      if (Math.abs(diff) > EPSILON) return diff > 0 ? 1 : -1;
      // Tie-break 1: more mapped bonds
      if (a.bondCount != b.bondCount) return b.bondCount - a.bondCount;
      // Tie-break 2: insertion order (stable sort)
      return a.insertionOrder - b.insertionOrder;
    });

    List<Map<Integer, Integer>> result = new ArrayList<>(scored.size());
    for (ScoredCandidate sc : scored) result.add(sc.mapping);
    return result;
  }

  // ---- Sub-score helpers ----

  private static Set<Integer> heteroAtomTypes(MolGraph g) {
    Set<Integer> types = new HashSet<>();
    for (int i = 0; i < g.n; i++) {
      int z = g.atomicNum[i];
      if (z != 6 && z != 1) types.add(z); // non-C, non-H
    }
    return types;
  }

  private static double scoreHetero(Map<Integer, Integer> m, MolGraph g1,
      Set<Integer> hUniverse) {
    if (hUniverse.isEmpty()) return 1.0;
    Set<Integer> hMapped = new HashSet<>();
    for (int qi : m.keySet()) {
      int z = g1.atomicNum[qi];
      if (z != 6 && z != 1 && hUniverse.contains(z)) hMapped.add(z);
    }
    return (double) hMapped.size() / hUniverse.size();
  }

  private static double scoreRare(Map<Integer, Integer> m, MolGraph g1,
      Set<Integer> hUniverse, double rUniverse) {
    if (rUniverse <= 0.0) return 1.0;
    Set<Integer> hMapped = new HashSet<>();
    for (int qi : m.keySet()) {
      int z = g1.atomicNum[qi];
      if (z != 6 && z != 1 && hUniverse.contains(z)) hMapped.add(z);
    }
    double rMapped = 0.0;
    for (int z : hMapped) rMapped += rarityWeight(z);
    return rMapped / rUniverse;
  }

  private static double scoreConnectivity(Map<Integer, Integer> m, MolGraph g1) {
    if (m.size() <= 1) return 1.0;
    // BFS to count connected components among mapped atoms in g1
    Set<Integer> mapped = m.keySet();
    Set<Integer> visited = new HashSet<>();
    int components = 0;
    for (int start : mapped) {
      if (visited.contains(start)) continue;
      components++;
      Deque<Integer> queue = new ArrayDeque<>();
      queue.add(start);
      visited.add(start);
      while (!queue.isEmpty()) {
        int u = queue.poll();
        for (int v : g1.neighbors[u]) {
          if (mapped.contains(v) && !visited.contains(v)) {
            visited.add(v);
            queue.add(v);
          }
        }
      }
    }
    return 1.0 / components;
  }

  private static int countMappedBonds(Map<Integer, Integer> m, MolGraph g1) {
    Set<Integer> mapped = m.keySet();
    int bonds = 0;
    for (int u : mapped) {
      for (int v : g1.neighbors[u]) {
        if (v > u && mapped.contains(v)) bonds++;
      }
    }
    return bonds;
  }

  private static final class ScoredCandidate {
    final Map<Integer, Integer> mapping;
    final double score;
    final int bondCount;
    final int insertionOrder;

    ScoredCandidate(Map<Integer, Integer> mapping, double score,
        int bondCount, int insertionOrder) {
      this.mapping = mapping;
      this.score = score;
      this.bondCount = bondCount;
      this.insertionOrder = insertionOrder;
    }
  }
}
