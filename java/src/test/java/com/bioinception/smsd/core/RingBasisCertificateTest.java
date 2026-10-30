/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
package com.bioinception.smsd.core;

import static com.bioinception.smsd.TestSupport.rawGraph;
import static org.junit.jupiter.api.Assertions.*;

import org.junit.jupiter.api.Test;

/** The minimum-basis certificate and the ring-perception work limit, as in the native core. */
class RingBasisCertificateTest {
  private static final long AMPLE = 1L << 34;

  private static String cyclodextrin(int units) {
    StringBuilder smiles = new StringBuilder("C%101C(O)C(O)C(OC1CO)O");
    for (int unit = 1; unit < units; unit++) smiles.append("C1C(O)C(O)C(OC1CO)O");
    return smiles.append("%10").toString();
  }

  private static int weight(int[][] cycles) {
    int total = 0;
    for (int[] cycle : cycles) total += cycle.length;
    return total;
  }

  @Test void documentedLimitAndMessage() {
    assertEquals(1L << 29, RingPerception.RING_SEARCH_WORK_LIMIT);
    assertEquals("Exact ring basis unavailable: ring perception reached its work limit before the minimum cycle basis was complete",
                 RingBasisUnavailableException.MESSAGE);
    assertTrue(new RingBasisUnavailableException() instanceof IllegalStateException);
  }

  @Test void fallbackBasisIsCertifiedExactlyWhenMinimum() throws Exception {
    for (String smiles : new String[] {"c1ccc2ccccc2c1", "C1C2CC3CC1CC(C2)C3", "C12C3C4C1C5C2C3C45",
                                       "C12C3C4C5C6C1C7C2C3C4C5C67", cyclodextrin(8)}) {
      MolGraph graph = rawGraph(smiles);
      RingPerception.Rings exact = RingPerception.perceive(graph.n, graph.neighbors, false, AMPLE);
      int rank = exact.smallest.length;
      assertFalse(exact.basisTruncated, smiles);
      assertTrue(RingPerception.certifiedMinimumBasis(graph.n, graph.neighbors, exact.smallest, rank, AMPLE), smiles);
      RingPerception.Rings bounded = RingPerception.perceive(graph.n, graph.neighbors, false, 0);
      if (bounded.truncated) {
        assertTrue(bounded.basisTruncated, smiles);
        boolean minimum = bounded.smallest.length == rank && weight(bounded.smallest) == weight(exact.smallest);
        assertEquals(minimum, RingPerception.certifiedMinimumBasis(graph.n, graph.neighbors, bounded.smallest, rank, AMPLE), smiles);
      }
      if (rank > 1) {
        int[][] dependent = exact.smallest.clone();
        dependent[rank - 1] = dependent[0];
        assertFalse(RingPerception.certifiedMinimumBasis(graph.n, graph.neighbors, dependent, rank, AMPLE), smiles);
      }
    }
  }

  @Test void largeCyclodextrinFallbackIsExact() throws Exception {
    MolGraph graph = rawGraph(cyclodextrin(34));
    RingPerception.Rings rings = RingPerception.perceive(graph.n, graph.neighbors, false, 1L << 20);
    assertTrue(rings.truncated);
    assertFalse(rings.basisTruncated);
    assertEquals(35, rings.smallest.length);
    assertEquals(34 * 11, weight(rings.smallest));
    assertFalse(graph.ringsTruncated());
  }
}
