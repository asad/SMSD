/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
package com.bioinception.smsd.core;

import static com.bioinception.smsd.TestSupport.rawGraph;
import static org.junit.jupiter.api.Assertions.*;

import java.util.Arrays;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.CsvSource;

class RelevantCycleCompletenessTest {
  @ParameterizedTest
  @CsvSource({"C1CCCCC1,1", "C12C3C4C1C5C4C3C25,6", "C1C2CC3CC1CC(C2)C3,4"})
  void ordinaryRelevantCyclesRemainComplete(String text, int count) throws Exception {
    MolGraph graph = rawGraph(text);
    assertEquals(count, graph.computeRelevantCycles().length);
    assertFalse(graph.ringsTruncated());
  }

  @Test
  void missingRelevantAlternativesReportIncompletePerception() {
    MolGraph exact = cyclobutaneMacrocycle(5);
    assertFalse(exact.relevantCyclesTruncated());
    assertEquals(5 + (1 << 5), exact.computeRelevantCycles().length);
    MolGraph graph = cyclobutaneMacrocycle(13);
    int count = graph.computeRelevantCycles().length;
    assertFalse(graph.ringsTruncated());
    assertTrue(graph.relevantCyclesTruncated());
    assertTrue(count <= 13 + (1 << 13), "the reported large relevant set is bounded");
  }

  private static MolGraph cyclobutaneMacrocycle(int units) {
    int atoms = 5 * units + 1;
    boolean[][] bonds = new boolean[atoms][atoms];
    for (int unit = 0; unit < units; unit++) {
      int corner = 3 * unit, opposite = 3 * (unit + 1);
      bond(bonds, corner, corner + 1);
      bond(bonds, corner + 1, opposite);
      bond(bonds, corner, corner + 2);
      bond(bonds, corner + 2, opposite);
    }
    int previous = 3 * units;
    for (int atom = 3 * units + 1; atom < atoms; atom++) {
      bond(bonds, previous, atom);
      previous = atom;
    }
    bond(bonds, previous, 0);
    int[] elements = new int[atoms];
    Arrays.fill(elements, 6);
    int[][] neighbours = new int[atoms][], orders = new int[atoms][];
    for (int atom = 0; atom < atoms; atom++) {
      int degree = 0;
      for (boolean bonded : bonds[atom]) if (bonded) degree++;
      neighbours[atom] = new int[degree];
      orders[atom] = new int[degree];
      Arrays.fill(orders[atom], 1);
      int index = 0;
      for (int other = 0; other < atoms; other++)
        if (bonds[atom][other]) neighbours[atom][index++] = other;
    }
    return new MolGraph.Builder().atomCount(atoms).atomicNumbers(elements)
        .neighbors(neighbours).bondOrders(orders).build();
  }

  private static void bond(boolean[][] bonds, int first, int second) {
    bonds[first][second] = bonds[second][first] = true;
  }
}
