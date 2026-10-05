/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
package com.bioinception.smsd.core;

import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.Random;
import org.junit.jupiter.api.Test;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;
import static org.junit.jupiter.api.Assertions.*;

class MCSAlgorithmRegressionTest {
  private static MolGraph graph(int[] elements, int[][] neighbors) {
    return new MolGraph.Builder().atomCount(elements.length).atomicNumbers(elements).neighbors(neighbors).build();
  }

  private static SearchEngine.MCSOptions disconnected(boolean induced) {
    SearchEngine.MCSOptions options = new SearchEngine.MCSOptions();
    options.connectedOnly = false;
    options.induced = induced;
    options.timeoutMs = 2_000;
    return options;
  }

  @Test
  void inducedMCSCannotUseFullGraphDegreesAsAnUpperBound() {
    MolGraph query = graph(new int[] {6, 8, 8}, new int[][] {{1}, {0}, {}});
    MolGraph target = graph(new int[] {6, 7, 8}, new int[][] {{}, {2}, {1}});
    Map<Integer, Integer> mapping = SearchEngine.findMCS(query, target, new ChemOptions(), disconnected(true));
    assertEquals(2, mapping.size(), "The isolated C/O pair is a valid induced common subgraph");
    assertTrue(SearchEngine.degreeSequenceUpperBound(query, target, new ChemOptions()) >= 2);
  }

  @Test
  void upperBoundIsAdmissibleForEveryPairOfFourVertexGraphs() {
    MolGraph[] graphs = new MolGraph[64];
    for (int mask = 0; mask < graphs.length; mask++) graphs[mask] = fourVertexGraph(mask);
    for (int q = 0; q < graphs.length; q++) {
      for (int t = 0; t < graphs.length; t++) {
        int exact = exhaustiveInducedSize(graphs[q], graphs[t]);
        int bound = SearchEngine.degreeSequenceUpperBound(graphs[q], graphs[t], new ChemOptions());
        assertTrue(bound >= exact, "Inadmissible bound for graph masks " + q + ", " + t);
      }
    }
  }

  @Test
  void mcsEnumerationAppliesInducedConstraintsToContainmentCandidates() {
    MolGraph query = graph(new int[] {6, 6}, new int[][] {{}, {}});
    MolGraph target = graph(new int[] {6, 6, 6}, new int[][] {{1}, {0, 2}, {1}});
    var mappings = SearchEngine.findAllMCS(query, target, new ChemOptions(), disconnected(true), 10);
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
    var expected = SearchEngine.findAllMCS(query, target, new ChemOptions(), new SearchEngine.MCSOptions(), 10);
    assertEquals(expected, SearchEngine.findAllMCS(query, target, null, null, 10));
    assertEquals(expected, SearchEngine.findAllMCS(query, target, null, new SearchEngine.MCSOptions(), 10));
  }

  @Test
  void queryAtomWeightsNeverMoveToTheTargetDuringOrientation() {
    MolGraph query = graph(new int[] {6, 8, 6}, new int[][] {{2}, {2}, {0, 1}});
    MolGraph target = graph(new int[] {6, 8, 7, 7}, new int[][] {{2}, {}, {0}, {}});
    SearchEngine.MCSOptions options = disconnected(false);
    options.atomWeights = new double[] {1, 10, 1};
    Map<Integer, Integer> mapping = assertDoesNotThrow(
        () -> SearchEngine.findMCS(query, target, new ChemOptions(), options));
    assertTrue(mapping.containsKey(1), "The query oxygen has the highest weight");
    assertTrue(SearchEngine.validateMapping(query, target, mapping, new ChemOptions()).isEmpty());
  }

  @Test
  void greedyExtensionPreservesEveryMappedQueryBond() {
    MolGraph query = graph(new int[] {6, 6, 8}, new int[][] {{1}, {0, 2}, {1}});
    MolGraph target = graph(new int[] {6, 6, 8}, new int[][] {{}, {}, {}});
    Map<Integer, Integer> extended = SearchEngine.greedyAtomExtend(
        query, target, Map.of(0, 0), new ChemOptions(), disconnected(false));
    assertEquals(Map.of(0, 0), extended);
    assertTrue(SearchEngine.validateMapping(query, target, extended, new ChemOptions()).isEmpty());
  }

  @Test
  void greedyInducedExtensionChecksMappedNonNeighbors() {
    MolGraph query = graph(new int[] {6, 6, 8}, new int[][] {{1}, {0, 2}, {1}});
    MolGraph target = graph(new int[] {6, 6, 8}, new int[][] {{1, 2}, {0, 2}, {0, 1}});
    Map<Integer, Integer> seed = Map.of(0, 0, 1, 1);
    assertEquals(seed, SearchEngine.greedyAtomExtend(query, target, seed, new ChemOptions(), disconnected(true)));
  }

  @Test
  void similarityScreenUsesTheRequestedChemicalConstraints() {
    MolGraph query = graph(new int[] {6, 6}, new int[][] {{1}, {0}});
    MolGraph target = new MolGraph.Builder().atomCount(3).atomicNumbers(new int[] {6, 6, 6})
        .neighbors(new int[][] {{1, 2}, {0, 2}, {0, 1}}).ringFlags(new boolean[] {true, true, true}).build();
    assertTrue(SearchEngine.similarityUpperBound(query, target, new ChemOptions()) >= 2.0 / 3.0,
        "Default matching permits ring atoms to match chain atoms");
    MolGraph oxygen = graph(new int[] {8}, new int[][] {{}});
    MolGraph nitrogen = graph(new int[] {7}, new int[][] {{}});
    ChemOptions untyped = new ChemOptions();
    untyped.matchAtomType = false;
    assertEquals(1.0, SearchEngine.similarityUpperBound(oxygen, nitrogen, untyped));
    ChemOptions fused = new ChemOptions();
    fused.ringFusionMode = ChemOptions.RingFusionMode.STRICT;
    assertDoesNotThrow(() -> SearchEngine.similarityUpperBound(query, target, fused));
  }

  @Test
  void connectedFilteringPreservesFirstLargestComponent() {
    MolGraph graph = graph(new int[] {6, 6, 6, 6, 6, 6},
        new int[][] {{1}, {0, 2}, {1}, {4}, {3, 5}, {4}});
    Map<Integer, Integer> mapping = new LinkedHashMap<>();
    for (int atom : new int[] {4, 0, 3, 1, 5, 2}) mapping.put(atom, atom + 10);
    assertEquals(Map.of(3, 13, 4, 14, 5, 15), SearchEngine.largestConnected(graph, mapping));
  }

  @Test
  void nonInducedMCSKeepsCallerQueryOrientation() {
    MolGraph query = graph(new int[] {6, 6, 6}, new int[][] {{}, {}, {}});
    MolGraph target = graph(new int[] {6, 6, 6}, new int[][] {{1}, {0}, {}});
    Map<Integer, Integer> mapping = SearchEngine.findMCS(query, target, new ChemOptions(), disconnected(false));
    assertEquals(3, mapping.size(), "Target edges are permitted when the query has no edges");
  }

  @Test
  void productGraphIncludesCompatibleNonEdgesAndPartialVertices() {
    MolGraph query = graph(new int[] {6, 6, 6}, new int[][] {{}, {}, {}});
    MolGraph target = graph(new int[] {6, 6, 6}, new int[][] {{1}, {0}, {}});
    query.ensureCanonical();
    target.ensureCanonical();
    SearchEngine.GraphBuilder product = new SearchEngine.GraphBuilder(query, target, new ChemOptions(), false);
    assertEquals(3, product.maximumCliqueSeed(new SearchEngine.TimeBudget(2_000)).size());
    assertEquals(3, product.mcSplitSeed(new SearchEngine.TimeBudget(2_000), new long[1]).size());
  }

  @Test
  void invalidExtensionCannotDisplaceALargerValidSeed() {
    MolGraph query = graph(new int[] {8, 7, 7, 6, 8}, new int[][] {{3}, {2, 3}, {1, 3}, {0, 1, 2}, {}});
    MolGraph target = graph(new int[] {6, 8, 7, 8, 7, 7}, new int[][] {{}, {}, {4}, {5}, {2}, {3}});
    Map<Integer, Integer> mapping = SearchEngine.findMCS(query, target, new ChemOptions(), disconnected(false));
    assertEquals(4, mapping.size());
    assertTrue(SearchEngine.validateMapping(query, target, mapping, new ChemOptions()).isEmpty());
  }

  @Test
  void emptyFallbackCannotEraseAValidChainSeed() {
    MolGraph query = graph(new int[] {6, 7, 8}, new int[][] {{2}, {2}, {0, 1}});
    MolGraph target = graph(new int[] {6, 8, 6}, new int[][] {{}, {}, {}});
    assertEquals(1, SearchEngine.findMCS(query, target, new ChemOptions(), disconnected(false)).size());
  }

  @Test
  void smallGraphSearchMatchesAnIndependentExhaustiveOracle() {
    Random random = new Random(194);
    for (int sample = 0; sample < 200; sample++) {
      MolGraph query = randomGraph(random, 3 + random.nextInt(4));
      MolGraph target = randomGraph(random, 3 + random.nextInt(4));
      for (boolean induced : new boolean[] {true, false}) {
        int expected = exhaustiveSize(query, target, induced);
        Map<Integer, Integer> mapping = SearchEngine.findMCS(query, target, new ChemOptions(), disconnected(induced));
        String context = "sample=" + sample + ", induced=" + induced;
        assertEquals(expected, mapping.size(), context);
        assertTrue(SearchEngine.validateMapping(query, target, mapping, new ChemOptions()).isEmpty(), context);
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
    SearchEngine.TimeBudget budget = new SearchEngine.TimeBudget(1);
    Thread.sleep(5);
    assertTrue(budget.expiredNow());
    for (int i = 0; i < 32; i++) assertTrue(budget.expired(), "Expiration must remain true between clock checks");
    assertEquals(0, budget.remainingMillis());
    SearchEngine.TimeBudget unchecked = new SearchEngine.TimeBudget(1);
    Thread.sleep(5);
    assertTrue(unchecked.expired(), "The first call must read the clock");
  }

  @Test
  void connectedFilteringUsesQueryWeights() {
    MolGraph graph = graph(new int[] {6, 6, 6, 6, 6, 6},
        new int[][] {{1}, {0, 2}, {1, 3}, {2}, {5}, {4}});
    SearchEngine.MCSOptions options = new SearchEngine.MCSOptions();
    options.atomWeights = new double[] {1, 1, 1, 1, 5, 5};
    assertEquals(Map.of(4, 4, 5, 5), SearchEngine.findMCS(graph, graph, new ChemOptions(), options));
  }

  @Test
  void connectedFilteringUsesMappedBondCount() {
    MolGraph graph = graph(new int[] {6, 6, 6, 6, 6, 6, 6, 6, 6},
        new int[][] {{1}, {0, 2}, {1, 3}, {2, 4}, {3}, {6, 7, 8}, {5, 7, 8}, {5, 6, 8}, {5, 6, 7}});
    SearchEngine.MCSOptions options = new SearchEngine.MCSOptions();
    options.maximizeBonds = true;
    Map<Integer, Integer> mapping = SearchEngine.findMCS(graph, graph, new ChemOptions(), options);
    assertEquals(Map.of(5, 5, 6, 6, 7, 7, 8, 8), mapping);
    assertEquals(6, SearchEngine.mcsScore(graph, mapping, options));
  }

  @Test
  void augmentationCannotEraseAnInducedIncumbent() {
    MolGraph query = graph(new int[] {8, 7, 6, 7}, new int[][] {{}, {3}, {3}, {1, 2}});
    MolGraph target = graph(new int[] {8, 7, 7, 8, 8, 7, 6},
        new int[][] {{3, 5}, {5, 6}, {5}, {0, 6}, {5}, {0, 1, 2, 4, 6}, {1, 3, 5}});
    SearchEngine.MCSOptions options = disconnected(true);
    Map<Integer, Integer> incumbent = new LinkedHashMap<>();
    incumbent.put(3, 1);
    incumbent.put(2, 6);
    incumbent.put(0, 0);
    Map<Integer, Integer> augmented = new LinkedHashMap<>(incumbent);
    augmented.put(1, 5);
    assertEquals(3, SearchEngine.ppx(query, target, incumbent, new ChemOptions(), options).size());
    assertEquals(2, SearchEngine.ppx(query, target, augmented, new ChemOptions(), options).size(),
        "A larger raw augmentation can shrink after induced filtering");
    assertTrue(SearchEngine.findMCS(query, target, new ChemOptions(), options).size() >= incumbent.size());
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
    SearchEngine.MCSOptions options = new SearchEngine.MCSOptions();
    options.timeoutMs = 300;
    long started = System.nanoTime();
    Map<Integer, Integer> mapping = SearchEngine.findMCS(query, target, new ChemOptions(), options);
    long elapsedMillis = (System.nanoTime() - started) / 1_000_000;
    assertTrue(elapsedMillis < 500, "Quality retries must share the 300ms budget; elapsed=" + elapsedMillis);
    assertTrue(SearchEngine.validateMapping(query, target, mapping, new ChemOptions()).isEmpty());
    assertEquals(300, options.timeoutMs, "The caller's options must remain unchanged");
    assertFalse(new SearchEngine.TimeBudget(1_000).expiredNow(), "The MCS budget scope must be cleared afterward");
  }

  @Test
  void rejectedMCSClearsItsBudgetScope() throws InterruptedException {
    MolGraph graph = graph(new int[] {6, 6}, new int[][] {{1}, {0}});
    SearchEngine.MCSOptions invalid = new SearchEngine.MCSOptions();
    invalid.timeoutMs = 1;
    invalid.atomWeights = new double[] {1};
    assertThrows(IllegalArgumentException.class, () -> SearchEngine.findMCS(graph, graph, new ChemOptions(), invalid));
    Thread.sleep(5);
    assertFalse(new SearchEngine.TimeBudget(1_000).expiredNow(), "A rejected call must restore the budget scope");
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
    Map<Integer, Integer> witness = Map.of(0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 10, 11, 11, 12);
    assertTrue(SearchEngine.validateMapping(query, target, witness, chemistry).isEmpty());
    SearchEngine.MCSOptions options = new SearchEngine.MCSOptions();
    options.timeoutMs = 2_000;
    Map<Integer, Integer> mapping = SearchEngine.findMCS(query, target, chemistry, options);
    assertTrue(mapping.size() >= witness.size(), "The valid seven-atom witness must not be missed");
    assertTrue(SearchEngine.validateMapping(query, target, mapping, chemistry).isEmpty());
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
    Map<Integer, Integer> witness = Map.ofEntries(
        Map.entry(0, 0),
        Map.entry(1, 1),
        Map.entry(2, 2),
        Map.entry(3, 17),
        Map.entry(4, 15),
        Map.entry(5, 16),
        Map.entry(6, 14),
        Map.entry(7, 18),
        Map.entry(8, 13),
        Map.entry(9, 19),
        Map.entry(10, 20),
        Map.entry(11, 21),
        Map.entry(12, 22),
        Map.entry(13, 28),
        Map.entry(14, 29),
        Map.entry(15, 30),
        Map.entry(16, 23),
        Map.entry(17, 27),
        Map.entry(18, 24),
        Map.entry(19, 26),
        Map.entry(20, 25),
        Map.entry(21, 12),
        Map.entry(22, 31),
        Map.entry(23, 11),
        Map.entry(24, 32));
    assertTrue(SearchEngine.validateMapping(query, target, witness, chemistry).isEmpty());
    assertEquals(witness, SearchEngine.largestConnected(query, witness));
    SearchEngine.MCSOptions options = new SearchEngine.MCSOptions();
    options.timeoutMs = 1_000;
    Map<Integer, Integer> mapping = SearchEngine.findMCS(query, target, chemistry, options);
    assertTrue(mapping.size() >= witness.size(), "The valid 25-atom common subgraph must not be missed");
    assertTrue(SearchEngine.validateMapping(query, target, mapping, chemistry).isEmpty());
    assertEquals(mapping, SearchEngine.largestConnected(query, mapping));
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
    SearchEngine.MCSOptions options = new SearchEngine.MCSOptions();
    options.induced = true;
    options.timeoutMs = 1_000;
    // A complete mapped fused query system would require all ten carbons,
    // while the target has only six. The pendant N-N bond is a valid witness.
    Map<Integer, Integer> witness = Map.of(10, 6, 11, 7);
    assertTrue(SearchEngine.validateMapping(query, target, witness, chemistry).isEmpty());
    assertEquals(witness, SearchEngine.ppx(query, target, witness, chemistry, options));
    Map<Integer, Integer> mapping = SearchEngine.findMCS(query, target, chemistry, options);
    assertEquals(2, mapping.size(), "Reversing an induced search must not admit partial caller-query rings");
    assertEquals(witness.keySet(), mapping.keySet());
    assertEquals(mapping, SearchEngine.ppx(query, target, mapping, chemistry, options));
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
