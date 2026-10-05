/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
package com.bioinception.smsd;

import com.bioinception.smsd.core.ChemOptions;
import com.bioinception.smsd.core.MolGraph;
import com.bioinception.smsd.core.SearchEngine;
import java.util.List;
import java.util.Map;
import org.junit.jupiter.api.Test;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;
import static org.junit.jupiter.api.Assertions.*;

class SearchApiRegressionTest {
  private static MolGraph ethane() {
    return new MolGraph.Builder().atomCount(2).atomicNumbers(new int[] {6, 6})
        .neighbors(new int[][] {{1}, {0}}).build();
  }

  @Test
  void graphTelemetryUsesTheSameDefaultsAsPlainSearch() {
    MolGraph query = ethane(), target = ethane();
    assertEquals(SearchEngine.isSubstructure(query, target, null, 0),
        SearchEngine.isSubstructureWithStats(query, target, null, 0).exists());
    assertEquals(SearchEngine.findAllSubstructures(query, target, null, 0, 0),
        SearchEngine.findAllSubstructuresWithStats(query, target, null, 0, 0).mappings());
  }

  @Test
  void cdkTelemetryAcceptsDefaultChemicalOptions() throws Exception {
    SmilesParser parser = new SmilesParser(SilentChemObjectBuilder.getInstance());
    var query = parser.parseSmiles("CC");
    var target = parser.parseSmiles("CCC");
    assertTrue(SearchEngine.isSubstructureWithStats(query, target, null, 1_000).exists());
    assertEquals(SearchEngine.findAllSubstructures(query, target, null, 0, 1_000),
        SearchEngine.findAllSubstructuresWithStats(query, target, null, 0, 1_000).mappings());
  }

  @Test
  void graphTelemetryHandlesMissingGraphsConsistently() {
    MolGraph target = ethane();
    assertFalse(SearchEngine.isSubstructureWithStats((MolGraph) null, target, null, 1_000).exists());
    assertTrue(SearchEngine.findAllSubstructuresWithStats(target, (MolGraph) null, null, 10, 1_000)
        .mappings().isEmpty());
  }

  @Test
  void invalidMappingIndicesAreReportedWithoutDereferencingThem() {
    MolGraph query = ethane(), target = ethane();
    for (Map<Integer, Integer> mapping : List.of(
        Map.of(-1, 0), Map.of(2, 0), Map.of(0, -1), Map.of(0, 2), Map.of(0, 0, 1, 2))) {
      List<String> errors = assertDoesNotThrow(
          () -> SearchEngine.validateMapping(query, target, mapping, new ChemOptions()));
      assertTrue(errors.stream().anyMatch(error -> error.contains("out of range")));
    }
    assertTrue(SearchEngine.validateMapping(query, target, Map.of(0, 0, 1, 1), null).isEmpty());
  }

  @Test
  void hugeTimeoutsDoNotOverflowIntoExpiredBudgets() {
    SearchEngine.TimeBudget budget = new SearchEngine.TimeBudget(Long.MAX_VALUE);
    assertFalse(budget.expiredNow());
    for (int i = 0; i < 2_048; i++) assertFalse(budget.expired());
  }
}
