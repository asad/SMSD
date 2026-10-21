/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
package com.bioinception.smsd;

import static com.bioinception.smsd.TestSupport.list;
import static com.bioinception.smsd.TestSupport.mapping;
import static org.junit.jupiter.api.Assertions.*;

import com.bioinception.smsd.core.*;
import java.util.List;
import java.util.Map;
import org.junit.jupiter.api.Test;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;

class SearchApiRegressionTest {
  private static MolGraph ethane() {
    return new MolGraph.Builder()
        .atomCount(2)
        .atomicNumbers(new int[] {6, 6})
        .neighbors(new int[][] {{1}, {0}})
        .build();
  }

  @Test
  void graphResultsUseTheSameDefaultsForExistenceAndEnumeration() {
    Substructure query = Substructure.find(ethane());
    SearchResult first = query.match(ethane()), all = query.matchAll(ethane());
    assertEquals(first.found(), all.found());
    assertTrue(first.found());
    assertEquals(2, all.matches().size());
    assertTrue(all.complete());
  }
  @Test
  void cdkResultsAcceptDefaultChemicalOptions() throws Exception {
    SmilesParser parser = new SmilesParser(SilentChemObjectBuilder.getInstance());
    IAtomContainer query = parser.parseSmiles("CC"), target = parser.parseSmiles("CCC");
    assertTrue(Substructure.find(query).match(target).found());
    assertEquals(4, Substructure.find(query).matchAll(target).matches().size());
  }
  @Test
  void missingGraphsAreRejectedConsistently() {
    assertThrows(NullPointerException.class, () -> Substructure.find((MolGraph) null));
    assertThrows(
        NullPointerException.class, () -> Substructure.find(ethane()).match((MolGraph) null));
    assertThrows(NullPointerException.class, () -> MCS.find((MolGraph) null));
    assertThrows(NullPointerException.class, () -> MCS.find(ethane()).match((MolGraph) null));
  }
  @Test
  void invalidMappingIndicesAreReportedWithoutDereferencingThem() {
    MolGraph query = ethane(), target = ethane();
    for (Map<Integer, Integer> mapping :
        list(mapping(-1, 0), mapping(2, 0), mapping(0, -1), mapping(0, 2), mapping(0, 0, 1, 2))) {
      List<String> errors =
          assertDoesNotThrow(() -> MCS.validateMapping(query, target, mapping, new ChemOptions()));
      assertTrue(errors.stream().anyMatch(error -> error.contains("out of range")));
    }
    assertTrue(MCS.validateMapping(query, target, mapping(0, 0, 1, 1), (ChemOptions) null).isEmpty());
  }
  @Test
  void anyNegativeTimeoutSelectsTheAdaptiveBudget() {
    MolGraph query = ethane(), target = ethane();
    for (long timeout : new long[] {-1L, -5L, Long.MIN_VALUE}) {
      SearchResult mcs = MCS.find(query).withTimeout(timeout).match(target);
      assertTrue(mcs.found() && mcs.optimal(), "timeout " + timeout);
      SearchResult embedding = Substructure.find(query).withTimeout(timeout).match(target);
      assertEquals(SearchResult.Status.Completed, embedding.status(), "timeout " + timeout);
      MCS.Options options = new MCS.Options();
      options.timeoutMillis = timeout;
      assertTrue(MCS.find(query).withOptions(options).match(target).found(), "timeout " + timeout);
    }
  }

  @Test
  void hugeTimeoutsDoNotOverflowIntoExpiredBudgets() {
    SearchResult result = Substructure.find(ethane()).withTimeout(Long.MAX_VALUE).match(ethane());
    assertTrue(result.found());
    assertEquals(SearchResult.Status.Completed, result.status());
    assertTrue(MCS.find(ethane()).withTimeout(Long.MAX_VALUE).match(ethane()).optimal());
  }
  @Test
  void immediateTimeoutDoesNotAffectConcurrentPreparedQueryReuse() {
    MCS query = MCS.find(ethane());
    assertEquals(SearchResult.Status.TimedOut, query.withTimeout(0).match(ethane()).status());
    assertTrue(query.match(ethane()).optimal());
  }
}
