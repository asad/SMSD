/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms. */
package com.bioinception.smsd.core;

import java.util.*;

/**
 * Bond-change-aware MCS post-filter that re-ranks candidates by the
 * plausibility of the implied bond transformations.
 *
 * <p>For reaction mapping, an MCS that implies many C-C bond breakings is
 * chemically less plausible than one that preserves most C-C bonds and
 * breaks only at expected reaction centres (heteroatom bonds).
 *
 * <p>Scoring: lower bond-change penalty = better mapping.
 * <pre>
 * penalty(m) = sum_{mapped_bonds} bondChangeCost(queryBondOrder, targetBondOrder)
 * </pre>
 *
 * Bond-change costs (biochemistry-informed):
 * <ul>
 *   <li>C-C bond broken/formed: 3.0 (rare in most biochemistry)</li>
 *   <li>C-N, C-O bond change: 1.5 (amide hydrolysis, esterification)</li>
 *   <li>C-S, C-P, heteroatom bonds: 0.5 (common reaction centres)</li>
 *   <li>Bond order change (single↔double): 1.0</li>
 *   <li>No change: 0.0</li>
 * </ul>
 *
 * @author Syed Asad Rahman
 * @since 6.5.3
 * @see MCSPostFilter
 * @see ReactionAwareScorer
 */
public final class BondChangeScorer implements MCSPostFilter {

  // Bond-change cost weights
  private static final double COST_CC_BREAK = 3.0;
  private static final double COST_CN_CO_CHANGE = 1.5;
  private static final double COST_HETERO_CHANGE = 0.5;
  private static final double COST_ORDER_CHANGE = 1.0;

  @Override
  public List<Map<Integer, Integer>> rank(
      List<Map<Integer, Integer>> candidates, MolGraph g1, MolGraph g2) {
    if (candidates == null || candidates.isEmpty()) return candidates;

    // Find the maximum candidate size K for normalization
    int K = 0;
    for (Map<Integer, Integer> m : candidates) {
      if (m.size() > K) K = m.size();
    }
    if (K == 0) return candidates;

    List<ScoredCandidate> scored = new ArrayList<>(candidates.size());
    for (int idx = 0; idx < candidates.size(); idx++) {
      Map<Integer, Integer> m = candidates.get(idx);
      double penalty = computePenalty(m, g1, g2);
      int size = m.size();

      // Combined score: reward size, penalise bond changes.
      // A K-1 mapping with 0 bond changes scores higher than a K mapping
      // with 3+ bond changes.  This allows near-MCS candidates to win
      // when they have significantly better bond alignment.
      //
      // score = size/K - penalty_weight * penalty/max_possible_penalty
      // penalty_weight = 0.15 (a K-1 with 0 penalty beats K with penalty > K*0.15/3)
      double sizeScore = (double) size / K;
      double maxPenalty = Math.max(1.0, K * COST_CC_BREAK); // theoretical max using global K
      double penaltyScore = penalty / maxPenalty;
      double combinedScore = sizeScore - 0.15 * penaltyScore;

      scored.add(new ScoredCandidate(m, combinedScore, size, idx));
    }

    // Sort: highest combined score first. Tie-break: larger MCS, then insertion order.
    scored.sort((a, b) -> {
      if (Math.abs(a.combinedScore() - b.combinedScore()) > 1e-9)
        return Double.compare(b.combinedScore(), a.combinedScore()); // higher combined score better
      if (a.size() != b.size()) return b.size() - a.size(); // larger MCS better
      return a.insertionOrder() - b.insertionOrder();
    });

    List<Map<Integer, Integer>> result = new ArrayList<>(scored.size());
    for (ScoredCandidate sc : scored) result.add(sc.mapping());
    return result;
  }

  /**
   * Compute bond-change penalty for a mapping.
   *
   * For each bond in g1 between two mapped atoms, check if the corresponding
   * bond in g2 exists and has the same order. Penalise changes based on the
   * element types involved.
   */
  static double computePenalty(Map<Integer, Integer> mapping,
                               MolGraph g1, MolGraph g2) {
    double penalty = 0.0;
    Set<Integer> mapped = mapping.keySet();

    for (int qi : mapped) {
      int ti = mapping.get(qi);
      for (int qj : g1.neighbors[qi]) {
        if (qj <= qi) continue; // count each bond once
        if (!mapping.containsKey(qj)) continue; // unmapped neighbor
        int tj = mapping.get(qj);

        int boQuery = g1.bondOrder(qi, qj);
        int boTarget = g2.bondOrder(ti, tj);

        if (boQuery == boTarget && boTarget > 0) continue; // no change

        // Determine element types at the bond
        int z1 = g1.atomicNum[qi], z2 = g1.atomicNum[qj];

        if (boTarget == 0) {
          // Bond in query, no bond in target — bond broken
          penalty += bondBreakCost(z1, z2);
        } else if (boQuery != boTarget) {
          // Bond order changed
          penalty += COST_ORDER_CHANGE;
        }
      }
    }

    // Also check for bonds in target between mapped atoms that
    // don't exist in query (bond formation)
    for (int qi : mapped) {
      int ti = mapping.get(qi);
      for (int tj : g2.neighbors[ti]) {
        // Find the query atom mapped to tj
        int qj = -1;
        for (Map.Entry<Integer, Integer> e : mapping.entrySet()) {
          if (e.getValue() == tj && e.getKey() > qi) { qj = e.getKey(); break; }
        }
        if (qj < 0) continue;
        int boQuery = g1.bondOrder(qi, qj);
        if (boQuery == 0) {
          // Bond in target but not in query — bond formed
          penalty += bondBreakCost(g1.atomicNum[qi], g1.atomicNum[qj]);
        }
      }
    }

    return penalty;
  }

  /** Cost of breaking/forming a bond between atoms of given elements. */
  private static double bondBreakCost(int z1, int z2) {
    // C-C: very rare breakage in biochemistry
    if (z1 == 6 && z2 == 6) return COST_CC_BREAK;
    // C-N or C-O: moderate (amide, ester)
    if ((z1 == 6 && (z2 == 7 || z2 == 8)) || (z2 == 6 && (z1 == 7 || z1 == 8)))
      return COST_CN_CO_CHANGE;
    // Heteroatom bonds (S, P, Se, etc.)
    return COST_HETERO_CHANGE;
  }

  private record ScoredCandidate(Map<Integer, Integer> mapping, double combinedScore,
      int size, int insertionOrder) {}
}
