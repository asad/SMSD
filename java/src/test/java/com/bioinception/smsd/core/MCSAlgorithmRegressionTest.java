/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
package com.bioinception.smsd.core;

import com.bioinception.smsd.core.MCS;

import com.bioinception.smsd.TestSupport;

import static com.bioinception.smsd.TestSupport.mapping;

import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.List;
import java.util.Random;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.CsvSource;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;
import static org.junit.jupiter.api.Assertions.*;

class MCSAlgorithmRegressionTest {
  @ParameterizedTest
  @CsvSource({
      "CCC(C)CC, CCCCCC, 5", "CCCCCC, CCC(C)CC, 5",
      "C#CC, CC(=O)O, 2", "CC(=O)O, C#CC, 2",
      "CCCCO, OC=O, 2", "OC=O, CCCCO, 2",
      "C#CC, NCC(=O)O, 2", "NCC(=O)O, C#CC, 2"
  })
  void connectedMCSSearchFindsTheOptimumInEitherArgumentOrder(
      String querySmiles, String targetSmiles, int expectedSize) throws Exception {
    SmilesParser parser = new SmilesParser(SilentChemObjectBuilder.getInstance());
    MolGraph query = new MolGraph(parser.parseSmiles(querySmiles));
    MolGraph target = new MolGraph(parser.parseSmiles(targetSmiles));
    ChemOptions chemistry = new ChemOptions();
    MCS.Options options = new MCS.Options();
    options.timeoutMillis = 1_000;
    Map<Integer, Integer> mapping = TestSupport.mcsMapping(query, target, chemistry, options);
    assertEquals(expectedSize, mapping.size());
    assertTrue(MCS.validateMapping(query, target, mapping, chemistry).isEmpty());
    assertEquals(mapping, TestSupport.largestConnected(query, mapping));
  }

  @Test
  void connectedMCSRetainsTheSharedAminoAlcoholAndRingAttachment() throws Exception {
    SmilesParser parser = new SmilesParser(SilentChemObjectBuilder.getInstance());
    MolGraph query = new MolGraph(Standardiser.standardise(parser.parseSmiles(
        "CC(C)NCC(O)COc1cccc2ccccc12"), Standardiser.TautomerMode.NONE));
    MolGraph target = new MolGraph(Standardiser.standardise(parser.parseSmiles(
        "CC(C)NCC(O)c1ccc(O)c(O)c1"), Standardiser.TautomerMode.NONE));
    ChemOptions chemistry = new ChemOptions();
    MCS.Options options = new MCS.Options();
    options.timeoutMillis = 2_000;
    Map<Integer, Integer> mapping = TestSupport.mcsMapping(query, target, chemistry, options);
    assertTrue(mapping.size() >= 9);
    assertTrue(MCS.validateMapping(query, target, mapping, chemistry).isEmpty());
    assertEquals(mapping, TestSupport.largestConnected(query, mapping));
  }

  private static MolGraph graph(int[] elements, int[][] neighbors) {
    return new MolGraph.Builder().atomCount(elements.length).atomicNumbers(elements).neighbors(neighbors).build();
  }

  private static MCS.Options disconnected(boolean induced) {
    MCS.Options options = new MCS.Options();
    options.connectivity = (false) ? MCS.Connectivity.Connected : MCS.Connectivity.Fragments;
    options.timeoutMillis = 2_000;
    return options;
  }

  @Test
  void inducedMCSCannotUseFullGraphDegreesAsAnUpperBound() {
    MolGraph query = graph(new int[] {6, 8, 8}, new int[][] {{1}, {0}, {}});
    MolGraph target = graph(new int[] {6, 7, 8}, new int[][] {{}, {2}, {1}});
    Map<Integer, Integer> mapping = TestSupport.mcsMapping(query, target, new ChemOptions(), disconnected(true));
    assertEquals(2, mapping.size(), "The isolated C/O pair is a valid induced common subgraph");
    assertTrue(TestSupport.atomBound(query, target, new ChemOptions()) >= 2);
  }

  @Test
  void upperBoundIsAdmissibleForEveryPairOfFourVertexGraphs() {
    MolGraph[] graphs = new MolGraph[64];
    for (int mask = 0; mask < graphs.length; mask++) graphs[mask] = fourVertexGraph(mask);
    for (int q = 0; q < graphs.length; q++) {
      for (int t = 0; t < graphs.length; t++) {
        int exact = exhaustiveInducedSize(graphs[q], graphs[t]);
        int bound = TestSupport.atomBound(graphs[q], graphs[t], new ChemOptions());
        assertTrue(bound >= exact, "Inadmissible bound for graph masks " + q + ", " + t);
      }
    }
  }

  @Test
  void mcsEnumerationAppliesInducedConstraintsToContainmentCandidates() {
    MolGraph query = graph(new int[] {6, 6}, new int[][] {{}, {}});
    MolGraph target = graph(new int[] {6, 6, 6}, new int[][] {{1}, {0, 2}, {1}});
    List<Map<Integer, Integer>> mappings = TestSupport.mcsMappings(query, target, new ChemOptions(), disconnected(true), 10);
    assertFalse(mappings.isEmpty());
    for (Map<Integer, Integer> mapping : mappings) {
      assertEquals(2, mapping.size());
      assertFalse(target.hasBond(mapping.get(0), mapping.get(1)),
          "An induced mapping of two isolated atoms cannot map to an adjacent pair");
    }
  }

  @Test
  void mcsEnumerationUsesDefaultOptionsConsistently() {
    MolGraph query = graph(new int[] {6, 6}, new int[][] {{1}, {0}});
    MolGraph target = graph(new int[] {6, 6, 6}, new int[][] {{1}, {0, 2}, {1}});
    List<Map<Integer, Integer>> expected = TestSupport.mcsMappings(query, target, new ChemOptions(), new MCS.Options(), 10);
    assertEquals(expected, TestSupport.mcsMappings(query, target, null, null, 10));
    assertEquals(expected, TestSupport.mcsMappings(query, target, null, new MCS.Options(), 10));
  }

  @Test
  void queryAtomWeightsNeverMoveToTheTargetDuringOrientation() {
    MolGraph query = graph(new int[] {6, 8, 6}, new int[][] {{2}, {2}, {0, 1}});
    MolGraph target = graph(new int[] {6, 8, 7, 7}, new int[][] {{2}, {}, {0}, {}});
    MCS.Options options = disconnected(false);
    options.atomWeights = new double[] {1, 10, 1};
    Map<Integer, Integer> mapping = assertDoesNotThrow(
        () -> TestSupport.mcsMapping(query, target, new ChemOptions(), options));
    assertTrue(mapping.containsKey(1), "The query oxygen has the highest weight");
    assertTrue(MCS.validateMapping(query, target, mapping, new ChemOptions()).isEmpty());
  }

  @Test
  void greedyExtensionPreservesEveryMappedQueryBond() {
    MolGraph query=graph(new int[]{6,6,8},new int[][]{{1},{0,2},{1}});
    MolGraph target=graph(new int[]{6,6,8},new int[][]{{},{},{}});
    Map<Integer,Integer> result=TestSupport.mcsMapping(query,target,new ChemOptions(),disconnected(true));
    assertEquals(2,result.size());
    assertTrue(MCS.validateMapping(query,target,result,new ChemOptions()).isEmpty());
  }

  @Test
  void greedyInducedExtensionChecksMappedNonNeighbors() {
    MolGraph query=graph(new int[]{6,6,8},new int[][]{{1},{0,2},{1}});
    MolGraph target=graph(new int[]{6,6,8},new int[][]{{1,2},{0,2},{0,1}});
    SearchResult result=MCS.find(query).withOptions(disconnected(true)).match(target);
    assertEquals(2,result.matches().get(0).atomCount());
    assertTrue(result.optimal());
  }

  @Test
  void similarityScreenUsesTheRequestedChemicalConstraints() {
    MolGraph query = graph(new int[] {6, 6}, new int[][] {{1}, {0}});
    MolGraph target = new MolGraph.Builder().atomCount(3).atomicNumbers(new int[] {6, 6, 6})
        .neighbors(new int[][] {{1, 2}, {0, 2}, {0, 1}}).ringFlags(new boolean[] {true, true, true}).build();
    assertTrue(MCS.similarityUpperBound(query, target, new ChemOptions()) >= 2.0 / 3.0,
        "Default matching permits ring atoms to match chain atoms");
    MolGraph oxygen = graph(new int[] {8}, new int[][] {{}});
    MolGraph nitrogen = graph(new int[] {7}, new int[][] {{}});
    ChemOptions untyped = new ChemOptions();
    untyped.matchAtomType = false;
    assertEquals(1.0, MCS.similarityUpperBound(oxygen, nitrogen, untyped));
    ChemOptions fused = new ChemOptions();
    fused.ringFusionMode = ChemOptions.RingFusionMode.STRICT;
    assertDoesNotThrow(() -> MCS.similarityUpperBound(query, target, fused));
  }

  @Test
  void connectedFilteringPreservesFirstLargestComponent() {
    MolGraph graph = graph(new int[] {6, 6, 6, 6, 6, 6},
        new int[][] {{1}, {0, 2}, {1}, {4}, {3, 5}, {4}});
    Map<Integer, Integer> mapping = new LinkedHashMap<>();
    for (int atom : new int[] {4, 0, 3, 1, 5, 2}) mapping.put(atom, atom + 10);
    assertEquals(mapping(3, 13, 4, 14, 5, 15), TestSupport.largestConnected(graph, mapping));
  }

  @Test
  void nonInducedMCSKeepsCallerQueryOrientation() {
    MolGraph query = graph(new int[] {6, 6, 6}, new int[][] {{}, {}, {}});
    MolGraph target = graph(new int[] {6, 6, 6}, new int[][] {{1}, {0}, {}});
    Map<Integer, Integer> mapping = TestSupport.mcsMapping(query, target, new ChemOptions(), disconnected(false));
    assertEquals(2, mapping.size(), "MCIS preserves mapped non-edges in both inputs");
  }

  @Test
  void productGraphIncludesCompatibleNonEdgesAndPartialVertices() {
    MolGraph query=graph(new int[]{6,6,6},new int[][]{{},{},{}});
    MolGraph target=graph(new int[]{6,6,6},new int[][]{{1},{0},{}});
    SearchResult forward=MCS.find(query).withOptions(disconnected(true)).match(target);
    SearchResult reverse=MCS.find(target).withOptions(disconnected(true)).match(query);
    assertEquals(2,forward.matches().get(0).atomCount());
    assertEquals(2,reverse.matches().get(0).atomCount());
    assertTrue(forward.optimal());assertTrue(reverse.optimal());
  }

  @Test
  void invalidExtensionCannotDisplaceALargerValidSeed() {
    MolGraph query = graph(new int[] {8, 7, 7, 6, 8}, new int[][] {{3}, {2, 3}, {1, 3}, {0, 1, 2}, {}});
    MolGraph target = graph(new int[] {6, 8, 7, 8, 7, 7}, new int[][] {{}, {}, {4}, {5}, {2}, {3}});
    Map<Integer, Integer> mapping = TestSupport.mcsMapping(query, target, new ChemOptions(), disconnected(false));
    assertEquals(4, mapping.size());
    assertTrue(MCS.validateMapping(query, target, mapping, new ChemOptions()).isEmpty());
  }

  @Test
  void emptyFallbackCannotEraseAValidChainSeed() {
    MolGraph query = graph(new int[] {6, 7, 8}, new int[][] {{2}, {2}, {0, 1}});
    MolGraph target = graph(new int[] {6, 8, 6}, new int[][] {{}, {}, {}});
    assertEquals(1, TestSupport.mcsMapping(query, target, new ChemOptions(), disconnected(false)).size());
  }

  @Test
  void smallGraphSearchMatchesAnIndependentExhaustiveOracle() {
    Random random = new Random(194);
    for (int sample = 0; sample < 200; sample++) {
      MolGraph query = randomGraph(random, 3 + random.nextInt(4));
      MolGraph target = randomGraph(random, 3 + random.nextInt(4));
      for (boolean induced : new boolean[] {true}) {
        int expected = exhaustiveSize(query, target, induced);
        Map<Integer, Integer> mapping = TestSupport.mcsMapping(query, target, new ChemOptions(), disconnected(induced));
        String context = "sample=" + sample + ", induced=" + induced;
        assertEquals(expected, mapping.size(), context);
        assertTrue(MCS.validateMapping(query, target, mapping, new ChemOptions()).isEmpty(), context);
        if (induced) {
          for (int qi : mapping.keySet())
            for (int qk : mapping.keySet())
              assertEquals(query.hasBond(qi, qk), target.hasBond(mapping.get(qi), mapping.get(qk)), context);
        }
      }
    }
  }

  @Test
  void observedDeadlineExpirationRemainsExpired() throws InterruptedException {
    SearchResult.Control budget = new SearchResult.Control(1,null);
    Thread.sleep(5);
    assertTrue(budget.stopped());
    for (int i = 0; i < 32; i++) assertTrue(budget.stopped(), "Expiration must remain true between clock checks");
    SearchResult.Control unchecked = new SearchResult.Control(1,null);
    Thread.sleep(5);
    assertTrue(unchecked.stopped(), "The first call must read the clock");
  }

  @Test
  void connectedFilteringUsesQueryWeights() {
    MolGraph graph = graph(new int[] {6, 6, 6, 6, 6, 6},
        new int[][] {{1}, {0, 2}, {1, 3}, {2}, {5}, {4}});
    MCS.Options options = new MCS.Options();
    options.atomWeights = new double[] {1, 1, 1, 1, 5, 5};
    Map<Integer,Integer> result=TestSupport.mcsMapping(graph,graph,new ChemOptions(),options);
    assertEquals(com.bioinception.smsd.TestSupport.set(4,5),result.keySet());
    assertTrue(MCS.validateMapping(graph,graph,result,new ChemOptions()).isEmpty());
  }

  @Test
  void connectedFilteringUsesMappedBondCount() {
    MolGraph graph = graph(new int[] {6, 6, 6, 6, 6, 6, 6, 6, 6},
        new int[][] {{1}, {0, 2}, {1, 3}, {2, 4}, {3}, {6, 7, 8}, {5, 7, 8}, {5, 6, 8}, {5, 6, 7}});
    MCS.Options options = new MCS.Options();
    options.objective = (true) ? MCS.Objective.MCES : MCS.Objective.MCIS;
    Map<Integer, Integer> mapping = TestSupport.mcsMapping(graph, graph, new ChemOptions(), options);
    assertEquals(mapping(5, 5, 6, 6, 7, 7, 8, 8), mapping);
    assertEquals(6, MCS.find(graph).withOptions(options).match(graph).score());
  }

  @Test
  void augmentationCannotEraseAnInducedIncumbent() {
    MolGraph query = graph(new int[] {8, 7, 6, 7}, new int[][] {{}, {3}, {3}, {1, 2}});
    MolGraph target = graph(new int[] {8, 7, 7, 8, 8, 7, 6},
        new int[][] {{3, 5}, {5, 6}, {5}, {0, 6}, {5}, {0, 1, 2, 4, 6}, {1, 3, 5}});
    MCS.Options options = disconnected(true);
    Map<Integer, Integer> incumbent = new LinkedHashMap<>();
    incumbent.put(3, 1);
    incumbent.put(2, 6);
    incumbent.put(0, 0);
    Map<Integer, Integer> augmented = new LinkedHashMap<>(incumbent);
    augmented.put(1, 5);
    assertTrue(MCS.validateMapping(query,target,incumbent,new ChemOptions()).isEmpty());
    assertFalse(query.hasBond(1,2));
    assertTrue(target.hasBond(augmented.get(1),augmented.get(2)),"The raw augmentation violates induced non-edge matching");
    assertTrue(TestSupport.mcsMapping(query, target, new ChemOptions(), options).size() >= incumbent.size());
  }

  @Test
  void qualityRetriesShareOneCallerTimeBudget() throws Exception {
    SmilesParser parser = new SmilesParser(SilentChemObjectBuilder.getInstance());
    MolGraph query = new MolGraph(Standardiser.standardise(parser.parseSmiles(
        "c1cc(c(c(c1)Cl)N2c3cc(cc(c3CNC2=O)c4ccc(cc4F)F)N5CCNCC5)Cl"), Standardiser.TautomerMode.NONE));
    MolGraph target = new MolGraph(Standardiser.standardise(parser.parseSmiles(
        "CCNc1cc(c2c(c1)N(C(=O)NC2)c3ccc(cc3)n4ccc-5ncnc5c4)c6ccnnc6"), Standardiser.TautomerMode.NONE));
    query.ensureCanonical();
    target.ensureCanonical();
    MCS.Options options = new MCS.Options();
    options.timeoutMillis = 300;
    long started = System.nanoTime();
    Map<Integer, Integer> mapping = TestSupport.mcsMapping(query, target, new ChemOptions(), options);
    long elapsedMillis = (System.nanoTime() - started) / 1_000_000;
    assertTrue(elapsedMillis < 500, "Quality retries must share the 300ms budget; elapsed=" + elapsedMillis);
    assertTrue(MCS.validateMapping(query, target, mapping, new ChemOptions()).isEmpty());
    assertEquals(300, options.timeoutMillis, "The caller's options must remain unchanged");
    assertFalse(new SearchResult.Control(1_000,null).stopped(), "The MCS budget scope must be cleared afterward");
  }

  @Test
  void rejectedMCSClearsItsBudgetScope() throws InterruptedException {
    MolGraph graph = graph(new int[] {6, 6}, new int[][] {{1}, {0}});
    MCS.Options invalid = new MCS.Options();
    invalid.timeoutMillis = 1;
    invalid.atomWeights = new double[] {1};
    assertThrows(IllegalArgumentException.class, () -> TestSupport.mcsMapping(graph, graph, new ChemOptions(), invalid));
    Thread.sleep(5);
    assertFalse(new SearchResult.Control(1_000,null).stopped(), "A rejected call must restore the budget scope");
  }

  @Test
  void partialConnectedSeedsExtendDespiteFullGraphNeighborhoodDifferences() throws Exception {
    SmilesParser parser = new SmilesParser(SilentChemObjectBuilder.getInstance());
    MolGraph query = new MolGraph(Standardiser.standardise(parser.parseSmiles("c1ccc(-c2ccccc2)cc1"),
        Standardiser.TautomerMode.NONE));
    MolGraph target = new MolGraph(Standardiser.standardise(parser.parseSmiles("c1ccc(Cc2ccccc2)cc1"),
        Standardiser.TautomerMode.NONE));
    ChemOptions chemistry = new ChemOptions();
    // Complete first ring plus its attachment. Default FLEXIBLE aromaticity
    // permits the query's aromatic junction to match the target methylene.
    Map<Integer, Integer> witness = mapping(0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 10, 11, 11, 12);
    assertTrue(MCS.validateMapping(query, target, witness, chemistry).isEmpty());
    MCS.Options options = new MCS.Options();
    options.timeoutMillis = 2_000;
    Map<Integer, Integer> mapping = TestSupport.mcsMapping(query, target, chemistry, options);
    assertTrue(mapping.size() >= witness.size(), "The valid seven-atom witness must not be missed");
    assertTrue(MCS.validateMapping(query, target, mapping, chemistry).isEmpty());
  }

  @Test
  void mediumMacrolideAnchorsRetainALargeConnectedMatch() throws Exception {
    SmilesParser parser = new SmilesParser(SilentChemObjectBuilder.getInstance());
    MolGraph query = new MolGraph(Standardiser.standardise(parser.parseSmiles(
        "CCC1OC(=O)C(C)C(OC2CC(C)(OC)C(O)C(C)O2)C(C)C(OC2OC(C)CC(C2O)N(C)C)C(C)(O)CC(C)C(=O)C(C)C(O)C1(C)O"), Standardiser.TautomerMode.NONE));
    MolGraph target = new MolGraph(Standardiser.standardise(parser.parseSmiles(
        "CCC1C(C(C(N(CC(CC(C(C(C(C(C(=O)O1)C)OC2CC(C(C(O2)C)O)(C)OC)C)OC3C(C(CC(O3)C)N(C)C)O)(C)O)C)C)C)O)(C)O"), Standardiser.TautomerMode.NONE));
    ChemOptions chemistry = new ChemOptions();
    // This connected 25-atom subset comes from an independently computed
    // 49-atom literal-graph witness, rather than the search being tested.
    Map<Integer, Integer> witness = mapping(
        0, 0,
        1, 1,
        2, 2,
        3, 17,
        4, 15,
        5, 16,
        6, 14,
        7, 18,
        8, 13,
        9, 19,
        10, 20,
        11, 21,
        12, 22,
        13, 28,
        14, 29,
        15, 30,
        16, 23,
        17, 27,
        18, 24,
        19, 26,
        20, 25,
        21, 12,
        22, 31,
        23, 11,
        24, 32);
    assertTrue(MCS.validateMapping(query, target, witness, chemistry).isEmpty());
    assertEquals(witness, TestSupport.largestConnected(query, witness));
    MCS.Options options = new MCS.Options();
    options.timeoutMillis = 1_000;
    Map<Integer, Integer> mapping = TestSupport.mcsMapping(query, target, chemistry, options);
    assertTrue(mapping.size() >= witness.size(), "The valid 25-atom common subgraph must not be missed");
    assertTrue(MCS.validateMapping(query, target, mapping, chemistry).isEmpty());
    assertEquals(mapping, TestSupport.largestConnected(query, mapping));
  }

  @Test
  void inducedOrientationPreservesCallerQueryRingCompleteness() {
    MolGraph query = graph(new int[] {6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7},
        new int[][] {{1, 9, 10}, {0, 2}, {1, 3}, {2, 8, 4}, {3, 5}, {4, 6},
            {5, 7}, {6, 8}, {3, 9, 7}, {8, 0}, {0, 11}, {10}});
    MolGraph target = graph(new int[] {6, 6, 6, 6, 6, 6, 7, 7},
        new int[][] {{1, 5, 6}, {0, 2}, {1, 3}, {2, 4}, {3, 5}, {4, 0}, {0, 7}, {6}});
    ChemOptions chemistry = new ChemOptions();
    chemistry.completeRingsOnly = true;
    MCS.Options options = new MCS.Options();
    options.timeoutMillis = 1_000;
    // A complete mapped fused query system would require all ten carbons,
    // while the target has only six. The pendant N-N bond is a valid witness.
    Map<Integer, Integer> witness = mapping(10, 6, 11, 7);
    assertTrue(MCS.validateMapping(query, target, witness, chemistry).isEmpty());
    assertEquals(2,witness.size());
    Map<Integer, Integer> mapping = TestSupport.mcsMapping(query, target, chemistry, options);
    assertEquals(8,mapping.size(),"A fusion atom may be retained through one complete relevant ring");
    assertTrue(MCS.validateMapping(query,target,mapping,chemistry).isEmpty());
  }

  private static MolGraph randomGraph(Random random, int size) {
    int[] elements = new int[size];
    for (int i = 0; i < size; i++) elements[i] = 6 + random.nextInt(3);
    boolean[][] edges = new boolean[size][size];
    for (int i = 0; i < size; i++)
      for (int j = i + 1; j < size; j++) edges[i][j] = edges[j][i] = random.nextInt(3) == 0;
    int[][] neighbors = new int[size][];
    for (int i = 0; i < size; i++) {
      int[] row = new int[size - 1];
      int count = 0;
      for (int j = 0; j < size; j++) if (edges[i][j]) row[count++] = j;
      neighbors[i] = Arrays.copyOf(row, count);
    }
    return graph(elements, neighbors);
  }

  private static int exhaustiveSize(MolGraph query, MolGraph target, boolean induced) {
    int[] mapping = new int[query.n];
    Arrays.fill(mapping, -1);
    return exhaustiveSize(query, target, induced, mapping, new boolean[target.n], 0);
  }

  private static int exhaustiveSize(MolGraph query, MolGraph target, boolean induced,
      int[] mapping, boolean[] used, int qi) {
    if (qi == query.n) return 0;
    int best = exhaustiveSize(query, target, induced, mapping, used, qi + 1);
    for (int tj = 0; tj < target.n; tj++) {
      if (used[tj] || query.atomicNum[qi] != target.atomicNum[tj]) continue;
      boolean valid = true;
      for (int qk = 0; qk < qi; qk++) {
        if (mapping[qk] < 0) continue;
        boolean queryBond = query.hasBond(qi, qk), targetBond = target.hasBond(tj, mapping[qk]);
        if ((queryBond && !targetBond) || (induced && queryBond != targetBond)) { valid = false; break; }
      }
      if (!valid) continue;
      mapping[qi] = tj;
      used[tj] = true;
      best = Math.max(best, 1 + exhaustiveSize(query, target, induced, mapping, used, qi + 1));
      used[tj] = false;
      mapping[qi] = -1;
    }
    return best;
  }

  private static MolGraph fourVertexGraph(int mask) {
    int[][] neighbors = new int[4][];
    boolean[][] edges = new boolean[4][4];
    for (int i = 0, bit = 0; i < 4; i++)
      for (int j = i + 1; j < 4; j++, bit++) edges[i][j] = edges[j][i] = (mask & (1 << bit)) != 0;
    for (int i = 0; i < 4; i++) {
      int[] row = new int[3];
      int count = 0;
      for (int j = 0; j < 4; j++) if (edges[i][j]) row[count++] = j;
      neighbors[i] = Arrays.copyOf(row, count);
    }
    return graph(new int[] {6, 6, 6, 6}, neighbors);
  }

  private static int exhaustiveInducedSize(MolGraph query, MolGraph target) {
    int[] mapping = new int[query.n];
    Arrays.fill(mapping, -1);
    return exhaustiveInducedSize(query, target, mapping, new boolean[target.n], 0);
  }

  private static int exhaustiveInducedSize(MolGraph query, MolGraph target, int[] mapping, boolean[] used, int qi) {
    if (qi == query.n) return 0;
    int best = exhaustiveInducedSize(query, target, mapping, used, qi + 1);
    for (int tj = 0; tj < target.n; tj++) {
      if (used[tj]) continue;
      boolean valid = true;
      for (int qk = 0; qk < qi; qk++) {
        if (mapping[qk] >= 0 && query.hasBond(qi, qk) != target.hasBond(tj, mapping[qk])) { valid = false; break; }
      }
      if (!valid) continue;
      mapping[qi] = tj;
      used[tj] = true;
      best = Math.max(best, 1 + exhaustiveInducedSize(query, target, mapping, used, qi + 1));
      used[tj] = false;
      mapping[qi] = -1;
    }
    return best;
  }
}
