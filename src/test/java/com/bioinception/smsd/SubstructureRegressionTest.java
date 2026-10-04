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
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.EnumSource;
import static org.junit.jupiter.api.Assertions.*;

/** Regression coverage for pruning, cache isolation, and complete enumeration. */
class SubstructureRegressionTest {
  private static final long TIMEOUT_MS = 5_000;

  private static ChemOptions options(ChemOptions.MatcherEngine engine) {
    ChemOptions options = new ChemOptions();
    options.matcherEngine = engine;
    return options;
  }

  private static MolGraph graph(int[] elements, int[][] neighbors) {
    return new MolGraph.Builder().atomCount(elements.length)
        .atomicNumbers(elements).neighbors(neighbors).build();
  }

  private static MolGraph isolatedCarbons(int count) {
    int[] elements = new int[count];
    Arrays.fill(elements, 6);
    return graph(elements, new int[count][]);
  }

  @ParameterizedTest
  @EnumSource(ChemOptions.MatcherEngine.class)
  void mutableOptionsDoNotReuseIncompatibleDomains(ChemOptions.MatcherEngine engine) {
    MolGraph query = new MolGraph.Builder().atomCount(1).atomicNumbers(new int[] {6})
        .neighbors(new int[][] {{}}).formalCharges(new int[] {1})
        .massNumbers(new int[] {13}).tetrahedralChirality(new int[] {1}).build();
    MolGraph target = new MolGraph.Builder().atomCount(1).atomicNumbers(new int[] {6})
        .neighbors(new int[][] {{}}).formalCharges(new int[] {0})
        .massNumbers(new int[] {14}).tetrahedralChirality(new int[] {2})
        .ringFlags(new boolean[] {true}).aromaticFlags(new boolean[] {true}).build();
    ChemOptions options = options(engine);

    assertTrue(SearchEngine.isSubstructure(query, target, options, TIMEOUT_MS));
    options.matchFormalCharge = true;
    assertFalse(SearchEngine.isSubstructure(query, target, options, TIMEOUT_MS));
    options.matchFormalCharge = false;
    options.matchIsotope = true;
    assertFalse(SearchEngine.isSubstructure(query, target, options, TIMEOUT_MS));
    options.matchIsotope = false;
    options.useChirality = true;
    assertFalse(SearchEngine.isSubstructure(query, target, options, TIMEOUT_MS));
    options.useChirality = false;
    options.ringMatchesRingOnly = true;
    assertFalse(SearchEngine.isSubstructure(query, target, options, TIMEOUT_MS));
    options.ringMatchesRingOnly = false;
    options.aromaticityMode = ChemOptions.AromaticityMode.STRICT;
    assertFalse(SearchEngine.isSubstructure(query, target, options, TIMEOUT_MS));
    options.aromaticityMode = ChemOptions.AromaticityMode.FLEXIBLE;
    assertTrue(SearchEngine.isSubstructure(query, target, options, TIMEOUT_MS));
  }

  @ParameterizedTest
  @EnumSource(ChemOptions.MatcherEngine.class)
  void relaxedAtomTypesAlsoRelaxNeighborLabels(ChemOptions.MatcherEngine engine) {
    MolGraph query = graph(new int[] {6, 8}, new int[][] {{1}, {0}});
    MolGraph target = graph(new int[] {6, 7}, new int[][] {{1}, {0}});
    ChemOptions options = options(engine);
    assertFalse(SearchEngine.isSubstructure(query, target, options, TIMEOUT_MS));
    options.matchAtomType = false;
    assertTrue(SearchEngine.isSubstructure(query, target, options, TIMEOUT_MS));
    assertEquals(2, SearchEngine.findAllSubstructures(query, target, options, 10, TIMEOUT_MS).size());
    options.matchAtomType = true;
    assertFalse(SearchEngine.isSubstructure(query, target, options, TIMEOUT_MS));
  }

  @ParameterizedTest
  @EnumSource(ChemOptions.MatcherEngine.class)
  void flexibleAromaticityAlsoRelaxesNeighborLabels(ChemOptions.MatcherEngine engine) {
    int[][] neighbors = {{1}, {0}};
    MolGraph query = graph(new int[] {6, 6}, neighbors);
    MolGraph target = new MolGraph.Builder().atomCount(2).atomicNumbers(new int[] {6, 6})
        .neighbors(neighbors).aromaticFlags(new boolean[] {true, true}).build();
    ChemOptions options = options(engine);
    assertTrue(SearchEngine.isSubstructure(query, target, options, TIMEOUT_MS));
    options.aromaticityMode = ChemOptions.AromaticityMode.STRICT;
    assertFalse(SearchEngine.isSubstructure(query, target, options, TIMEOUT_MS));
  }

  @ParameterizedTest
  @EnumSource(ChemOptions.MatcherEngine.class)
  void extraTargetEdgesCannotInvalidateMultiHopPruning(ChemOptions.MatcherEngine engine) {
    // In a clique all mapped path neighbors are distance 1. Exact distance
    // shells at hops 2 and 3 therefore cannot be used as rejection tests.
    int count = 22;
    int[] elements = new int[count];
    Arrays.fill(elements, 6);
    int[][] path = new int[count][], clique = new int[count][];
    for (int i = 0; i < count; i++) {
      path[i] = i == 0 ? new int[] {1} : i == count - 1 ? new int[] {i - 1}
          : new int[] {i - 1, i + 1};
      clique[i] = new int[count - 1];
      for (int j = 0, k = 0; j < count; j++) if (j != i) clique[i][k++] = j;
    }
    MolGraph query = graph(elements, path), target = graph(elements, clique);
    ChemOptions options = options(engine);
    options.useThreeHopNLF = true;
    assertTrue(SearchEngine.isSubstructure(query, target, options, TIMEOUT_MS));
    List<Map<Integer, Integer>> maps = SearchEngine.findAllSubstructures(query, target, options, 1, TIMEOUT_MS);
    assertEquals(1, maps.size());
    assertEquals(count, maps.get(0).size());
    assertTrue(SearchEngine.validateMapping(query, target, maps.get(0), options).isEmpty());
    options.useTwoHopNLF = false;
    assertTrue(SearchEngine.isSubstructure(query, target, options, TIMEOUT_MS));
  }

  @ParameterizedTest
  @EnumSource(ChemOptions.MatcherEngine.class)
  void enumerationIncludesTargetsPast4096(ChemOptions.MatcherEngine engine) {
    List<Map<Integer, Integer>> maps = SearchEngine.findAllSubstructures(
        isolatedCarbons(1), isolatedCarbons(4_200), options(engine), 5_000, TIMEOUT_MS);
    assertEquals(4_200, maps.size());
    assertEquals(4_199, maps.getLast().get(0));
  }

  @ParameterizedTest
  @EnumSource(ChemOptions.MatcherEngine.class)
  void anyBondOrderStillChecksStrictAromaticity(ChemOptions.MatcherEngine engine) {
    MolGraph query = new MolGraph.Builder().atomCount(2).atomicNumbers(new int[] {6, 6})
        .neighbors(new int[][] {{1}, {0}}).aromaticFlags(new boolean[] {true, true})
        .bondAromaticFlags(new boolean[][] {{true}, {true}}).build();
    MolGraph target = new MolGraph.Builder().atomCount(2).atomicNumbers(new int[] {6, 6})
        .neighbors(new int[][] {{1}, {0}}).aromaticFlags(new boolean[] {true, true})
        .bondAromaticFlags(new boolean[][] {{false}, {false}}).build();
    ChemOptions options = options(engine);
    options.matchBondOrder = ChemOptions.BondOrderMode.ANY;
    options.aromaticityMode = ChemOptions.AromaticityMode.STRICT;
    assertFalse(SearchEngine.isSubstructure(query, target, options, TIMEOUT_MS));
    assertTrue(SearchEngine.findAllSubstructures(query, target, options, 10, TIMEOUT_MS).isEmpty());
  }
}
