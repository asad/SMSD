/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
package com.bioinception.smsd.core;

import static com.bioinception.smsd.TestSupport.rawGraph;
import static org.junit.jupiter.api.Assertions.*;

import java.util.*;
import java.util.concurrent.*;
import java.util.concurrent.atomic.AtomicBoolean;
import org.junit.jupiter.api.Test;

class SubstructureRingCacheTest {
  private static SearchResult.Control control() {
    return new SearchResult.Control(60_000, null);
  }

  private static boolean original(MolGraph query, int[][] neighbors, MolGraph target, int[] mapping) {
    Set<Integer> queryAtoms = new HashSet<>(), targetAtoms = new HashSet<>();
    Set<SearchResult.Bond> queryBonds = new HashSet<>(), targetBonds = new HashSet<>();
    for (int atom = 0; atom < mapping.length; atom++) {
      queryAtoms.add(atom);
      targetAtoms.add(mapping[atom]);
      for (int neighbor : neighbors[atom]) if (neighbor > atom) {
        queryBonds.add(new SearchResult.Bond(atom, neighbor));
        targetBonds.add(new SearchResult.Bond(mapping[atom], mapping[neighbor]));
      }
    }
    return (query == null || VertexSearch.completeRings(query, queryAtoms, queryBonds))
        && VertexSearch.completeRings(target, targetAtoms, targetBonds);
  }

  private static MolGraph graph(int atoms, int mask, boolean atomFlag, boolean bondFlag) {
    int[] elements = new int[atoms];
    Arrays.fill(elements, 6);
    boolean[] flags = new boolean[atoms];
    if (atomFlag && atoms > 0) flags[0] = true;
    int[][] neighbors = new int[atoms][], orders = new int[atoms][atoms];
    boolean[][] bondFlags = new boolean[atoms][atoms];
    int bit = 0;
    for (int first = 0; first < atoms; first++)
      for (int second = first + 1; second < atoms; second++, bit++)
        if ((mask & (1 << bit)) != 0) {
          orders[first][second] = orders[second][first] = 1;
          bondFlags[first][second] = bondFlags[second][first] = bondFlag;
        }
    for (int atom = 0; atom < atoms; atom++) {
      int count = 0;
      for (int other = 0; other < atoms; other++) if (orders[atom][other] != 0) count++;
      neighbors[atom] = new int[count];
      count = 0;
      for (int other = 0; other < atoms; other++)
        if (orders[atom][other] != 0) neighbors[atom][count++] = other;
    }
    return new MolGraph.Builder().atomCount(atoms).atomicNumbers(elements).ringFlags(flags)
        .neighbors(neighbors).bondOrders(orders).bondRingFlags(bondFlags).build();
  }

  @Test
  void matchesOriginalPredicateForEveryThreeAtomGraphAndFourAtomImage() {
    long comparisons = 0;
    for (int targetMask = 0; targetMask < 64; targetMask++) {
      MolGraph target = graph(4, targetMask, false, false);
      for (int queryMask = 0; queryMask < 8; queryMask++) {
        MolGraph query = graph(3, queryMask, false, false);
        Substructure.CompleteRingCheck molecular = new Substructure.CompleteRingCheck(
            query, query.neighbors, target, control());
        Substructure.CompleteRingCheck smarts = new Substructure.CompleteRingCheck(
            null, query.neighbors, target, control());
        for (int first = 0; first < 4; first++)
          for (int second = 0; second < 4; second++) if (second != first)
            for (int third = 0; third < 4; third++) if (third != first && third != second) {
              int[] mapping = {first, second, third};
              assertEquals(original(query, query.neighbors, target, mapping), molecular.holds(mapping));
              assertEquals(original(null, query.neighbors, target, mapping), smarts.holds(mapping));
              comparisons += 2;
            }
      }
    }
    assertEquals(24_576, comparisons);
  }

  @Test
  void preservesFlaggedAtomsAndBondsWithoutRelevantCycles() {
    int[] mapping = {0, 1};
    for (int flags = 0; flags < 4; flags++) {
      MolGraph marked = graph(2, 1, (flags & 1) != 0, (flags & 2) != 0);
      MolGraph plain = graph(2, 1, false, false);
      for (MolGraph query : new MolGraph[] {marked, plain})
        for (MolGraph target : new MolGraph[] {marked, plain}) {
          Substructure.CompleteRingCheck check = new Substructure.CompleteRingCheck(
              query, query.neighbors, target, control());
          assertEquals(original(query, query.neighbors, target, mapping), check.holds(mapping));
          assertEquals(original(query, query.neighbors, target, mapping), check.holds(mapping));
        }
    }
  }

  @Test
  void targetExtraEdgesDoNotCompleteAQueryChain() throws Exception {
    MolGraph query = rawGraph("CCC"), target = rawGraph("C1CC1");
    assertFalse(new Substructure.CompleteRingCheck(query, query.neighbors, target, control())
        .holds(new int[] {0, 1, 2}));
    assertFalse(new Substructure.CompleteRingCheck(null, query.neighbors, target, control())
        .holds(new int[] {0, 1, 2}));
  }

  @Test
  void retainsFusedBridgedSpiroAndMacrocycleCoverageAcrossEmbeddings() throws Exception {
    for (String smiles : new String[] {"c1ccc2ccccc2c1", "C1CC2CCC1C2", "C1CCC2(CC1)CCCC2",
        "C12C3C4C1C5C2C3C45", "C1CCCCCCCCCCC1", "CC1CC1.CCC"}) {
      MolGraph query = rawGraph(smiles), target = rawGraph(smiles);
      Substructure.CompleteRingCheck check = new Substructure.CompleteRingCheck(
          query, query.neighbors, target, control());
      int[] mapping = new int[query.n];
      for (int shift = 0; shift < query.n; shift++) {
        for (int atom = 0; atom < query.n; atom++) mapping[atom] = (atom + shift) % query.n;
        assertEquals(original(query, query.neighbors, target, mapping), check.holds(mapping), smiles);
      }
    }
  }

  @Test
  void cachedCoverageDoesNotOverrideCancellation() throws Exception {
    MolGraph graph = rawGraph("C1CC1");
    AtomicBoolean cancelled = new AtomicBoolean();
    SearchResult.Control control = new SearchResult.Control(60_000, cancelled::get);
    Substructure.CompleteRingCheck check = new Substructure.CompleteRingCheck(
        graph, graph.neighbors, graph, control);
    assertTrue(check.holds(new int[] {0, 1, 2}));
    cancelled.set(true);
    assertFalse(check.holds(new int[] {0, 1, 2}));
    assertEquals(SearchResult.Status.Cancelled, control.status());
  }

  @Test
  void preparedQueriesKeepCoverageLocalToConcurrentSearches() throws Exception {
    Substructure.Options options = new Substructure.Options();
    options.chemistry.completeRingsOnly = true;
    options.timeoutMillis = 60_000;
    MolGraph queryGraph = rawGraph("C1CC1");
    ExecutorService executor = Executors.newFixedThreadPool(2);
    try {
      for (Substructure query : new Substructure[] {Substructure.find(queryGraph).withOptions(options),
          Substructure.smarts("C1CC1").withOptions(options)}) {
        MolGraph target = rawGraph("C1CC1.C1CC1");
        SearchResult expected = query.matchAll(target, 100);
        assertTrue(expected.found());
        List<Future<SearchResult>> results = new ArrayList<>();
        for (int iteration = 0; iteration < 6; iteration++)
          results.add(executor.submit(() -> query.matchAll(target, 100)));
        for (Future<SearchResult> task : results) {
          SearchResult result = task.get();
          assertEquals(expected.matches().size(), result.matches().size());
          assertEquals(SearchResult.Status.Completed, result.status());
          assertTrue(result.complete());
          assertTrue(result.optimal());
        }
      }
    } finally {
      executor.shutdownNow();
    }
  }
}
