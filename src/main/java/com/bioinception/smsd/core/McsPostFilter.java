/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms. */
package com.bioinception.smsd.core;

import java.util.List;
import java.util.Map;

/**
 * Post-filter that re-ranks MCS candidates for reaction relevance.
 * Implementations receive the full candidate list and both molecule
 * graphs, and must return the candidates in preferred order.
 *
 * @author Syed Asad Rahman
 * @since 6.4.0
 */
@FunctionalInterface
public interface McsPostFilter {
  /**
   * Rank MCS candidates by domain-specific relevance.
   *
   * @param candidates unmodifiable list of MCS candidates (size K down to K-delta)
   * @param g1         first molecule graph
   * @param g2         second molecule graph
   * @return candidates sorted best-first; may be a subset
   */
  List<Map<Integer, Integer>> rank(
      List<Map<Integer, Integer>> candidates, MolGraph g1, MolGraph g2);
}
