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
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.CsvSource;
class RingPerceptionCompatibilityTest {
  private static String hoop(int units) {
    StringBuilder text = new StringBuilder("c%101ccc(cc1)");
    for (int unit = 1; unit + 1 < units; unit++) text.append("c1ccc(cc1)");
    return text.append("c1ccc%10cc1").toString();
  }
  @ParameterizedTest
  @CsvSource({"C1CCCC1,C12C3C4C5C1C6C2C3C4C56,5", "C1CCCCC1,C12C3C4C5C6C1C7C2C3C4C5C67,6", "C1CCCC1,C12C3C4C1C5C3C4C25,5"})
  void completeRingPrismsAndCage(String query, String target, int expected) throws Exception {
    for (MCS.Objective objective : MCS.Objective.values()) {
      MCS.Options options = new MCS.Options();
      options.timeoutMillis = 10000;
      options.objective = objective;
      options.chemistry.completeRingsOnly = true;
      SearchResult result = MCS.find(rawGraph(query)).withOptions(options).match(rawGraph(target));
      assertEquals(SearchResult.Status.Completed, result.status());
      assertTrue(result.optimal()); assertTrue(result.complete()); assertEquals(expected, result.score());
    }
  }
  @Test void basisAndRelevantFlagsAreIndependent() throws Exception {
    MolGraph graph = rawGraph(hoop(10));
    int[][] basis = graph.computeRings();
    assertFalse(graph.ringsTruncated());
    graph.computeRelevantCycles();
    assertSame(basis, graph.computeRings());
    assertFalse(graph.ringsTruncated());
    assertTrue(graph.relevantCyclesTruncated());
    assertFalse(graph.ringCountsIncomplete());
    SearchResult found = Substructure.smarts("[R2]").withTimeout(10000).match(graph);
    assertEquals(SearchResult.Status.Completed, found.status());
    assertFalse(found.matches().isEmpty());
  }
  @Test void numericPredicatesRefuseBeforeNegationOrRecursion() throws Exception {
    MolGraph graph = rawGraph("C12C3C4C5C6C1C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C1C6C5C4C3C21");
    assertTrue(graph.ringsTruncated());
    for (String text : new String[] {"[R1]", "[!R1]", "[R1,#6]", "[$([R1])]", "[!$([R1])]"}) {
      SearchResult result = Substructure.smarts(text).withTimeout(10000).match(graph);
      assertEquals(SearchResult.Status.LimitReached, result.status(), text);
      assertFalse(result.optimal()); assertFalse(result.complete()); assertTrue(result.matches().isEmpty());
      result = Substructure.smarts(text).withTimeout(0).match(graph);
      assertEquals(SearchResult.Status.TimedOut, result.status(), text);
    }
    for (String text : new String[] {"[R]", "[R0]", "[!R0]", "[$([R0])]", "[#6]"}) {
      SearchResult result = Substructure.smarts(text).withTimeout(10000).match(graph);
      assertEquals(SearchResult.Status.Completed, result.status(), text);
      assertTrue(result.optimal()); assertTrue(result.complete());
    }
  }
}
