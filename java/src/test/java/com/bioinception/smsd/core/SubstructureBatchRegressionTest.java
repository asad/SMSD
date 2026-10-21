/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
package com.bioinception.smsd.core;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.IdentityHashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;
import org.junit.jupiter.api.Test;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smarts.SmartsPattern;
import org.openscience.cdk.smiles.SmilesParser;
import static com.bioinception.smsd.TestSupport.rawGraph;
import static org.junit.jupiter.api.Assertions.*;

class SubstructureBatchRegressionTest {
  /** Targets with fused, bridged and spiro rings, stereo, charges and several components. */
  private static final String[] TARGETS = {
      "CC(=O)Nc1ccc(O)cc1", "OC(=O)c1ccccc1OC(C)=O", "Cn1cnc2c1c(=O)n(C)c(=O)n2C",
      "C[NH3+].[O-]C(=O)CC(F)(F)F", "c1ccc2[nH]ccc2c1", "NS(=O)(=O)c1ccc(Cl)cc1", "C=CC#N",
      "O=C1NC(=O)C(N1)=O", "CC(C)(C)OC(=O)N[C@@H](Cc1ccccc1)C(=O)O", "N[C@@H](C)C(=O)O",
      "N[C@H](C)C(=O)O", "F/C=C/Cl", "F/C=C\\Cl", "C1CCC2CCCCC2C1", "C1C2CC3CC1CC(C2)C3",
      "C1CC2CCC1C2", "C1CCC2(CC1)CCCC2", "c1ccc2c(c1)ccc1ccccc12", "C1CC2CCC3CCCC4CCC(C1)C2C34",
      "OC1CCC2C(C1)CCC1C2CCC2CCCC12", "CCCCCC", "C.O", "CO", "c1ccccc1CC",
      "O=C(O)CCc1ccccc1.Cl", "CC1=CC(=O)C=CC1=O"};

  /** SMARTS with recursion, ring counts and sizes, a negated stereo test, stereo and components. */
  private static final String[] PATTERNS = {
      "[CX3](=O)[NX3]", "[$(C(=O)[O,N])]", "[R2]", "[r5]", "[x3;R]", "[R]@[R]@[R]",
      "[!R]-[R]", "[r6;!r5]~[R2]", "[C;R2;!@]", "[$([R2]~[R2])]", "c1ccccc1",
      "N[C@@H](C)C(=O)O", "F/C=C/Cl", "(C).(O)", "[#6].[#8]", "*~*~*"};

  /** Every deterministic part of two results: matches, status, completeness, score and nodes. */
  private static void same(SearchResult expected, SearchResult actual, String label) {
    assertEquals(expected.matches(), actual.matches(), label);
    assertEquals(expected.status(), actual.status(), label);
    assertEquals(expected.complete(), actual.complete(), label);
    assertEquals(expected.optimal(), actual.optimal(), label);
    assertEquals(expected.score(), actual.score(), label);
    assertEquals(expected.nodesVisited(), actual.nodesVisited(), label);
  }

  @Test void moleculeBatchRetainsInputOrderAndCanBeReused() throws Exception {
    Substructure query = Substructure.find(rawGraph("CC"));
    List<MolGraph> targets = Arrays.asList(rawGraph("CCC"), rawGraph("N"), rawGraph("CC"));
    for (int repeat = 0; repeat < 2; repeat++) {
      List<SearchResult> results = query.matchMany(targets);
      assertEquals(3, results.size());
      assertTrue(results.get(0).found());
      assertFalse(results.get(1).found());
      assertTrue(results.get(2).found());
      for (SearchResult result : results) {
        assertEquals(SearchResult.Status.Completed, result.status());
        assertTrue(result.complete());
        if (result.found()) assertEquals(2, result.matches().get(0).atomCount());
      }
      assertThrows(UnsupportedOperationException.class, () -> results.clear());
    }
    assertTrue(query.matchMany(Collections.emptyList()).isEmpty());
  }

  @Test void smartsBatchUsesItsCompiledQueryForEveryTarget() throws Exception {
    List<SearchResult> results = Substructure.smarts("[CX3](=O)[NX3]")
        .matchMany(Arrays.asList(rawGraph("CCO"), rawGraph("CC(=O)N"), rawGraph("NC(=O)N")));
    assertFalse(results.get(0).found());
    assertTrue(results.get(1).found());
    assertTrue(results.get(2).found());
    assertEquals(3, results.get(1).matches().get(0).atomCount());
  }

  @Test void stoppedOutcomesRemainAssociatedWithEveryTarget() throws Exception {
    List<MolGraph> targets = Arrays.asList(rawGraph("CCC"), rawGraph("N"));
    for (SearchResult result : Substructure.find(rawGraph("CC")).withTimeout(0).matchMany(targets)) {
      assertEquals(SearchResult.Status.TimedOut, result.status());
      assertFalse(result.complete());
      assertFalse(result.found());
    }
    Substructure.Options options = new Substructure.Options();
    options.cancellation = () -> true;
    for (SearchResult result : Substructure.find(rawGraph("CC")).withOptions(options).matchMany(targets)) {
      assertEquals(SearchResult.Status.Cancelled, result.status());
      assertFalse(result.complete());
    }
  }

  /** Batches on threads return, in input order, what searches of fresh graphs in turn return. */
  @Test void smartsBatchesMatchSearchesInTurnOnEveryThreadCount() throws Exception {
    List<MolGraph> targets = new ArrayList<>();
    // Each graph appears three times, so that threads also read one target at the same time.
    for (int copy = 0; copy < 3; copy++)
      for (String smiles : TARGETS) targets.add(rawGraph(smiles));
    for (int copy = 1; copy < 3; copy++)
      for (int index = 0; index < TARGETS.length; index++)
        targets.set(copy * TARGETS.length + index, targets.get(index));
    for (String pattern : PATTERNS) {
      Substructure query = Substructure.smarts(pattern);
      List<SearchResult> inTurn = new ArrayList<>();
      for (int index = 0; index < targets.size(); index++)
        inTurn.add(Substructure.smarts(pattern).match(rawGraph(TARGETS[index % TARGETS.length])));
      for (int threads : new int[] {1, 2, 4, 0}) {
        List<SearchResult> batch = query.matchMany(targets, threads);
        assertEquals(targets.size(), batch.size());
        for (int index = 0; index < targets.size(); index++)
          same(inTurn.get(index), batch.get(index), pattern + " on " + threads + " threads, target " + index);
      }
    }
  }

  /** One prepared SMARTS query shared by four threads over many targets finds what it finds in turn. */
  @Test void oneSmartsQuerySharedByFourThreadsMatchesAsInTurn() throws Exception {
    String[] patterns = {"[$([R2]~[R2]),$(C(=O)[O,N]),r6]~*", "[R2;!@]~[!R]", "[x3;R]@[R]"};
    List<MolGraph> targets = new ArrayList<>();
    for (int copy = 0; copy < 4; copy++)
      for (String smiles : TARGETS) targets.add(rawGraph(smiles));
    for (String pattern : patterns) {
      Substructure query = Substructure.smarts(pattern);
      List<SearchResult> expected = new ArrayList<>();
      for (String smiles : TARGETS) expected.add(Substructure.smarts(pattern).matchAll(rawGraph(smiles)));
      ExecutorService threads = Executors.newFixedThreadPool(4);
      try {
        CountDownLatch start = new CountDownLatch(1);
        List<Future<?>> running = new ArrayList<>();
        for (int thread = 0; thread < 4; thread++) {
          int offset = thread * 7;
          running.add(threads.submit(() -> {
            start.await();
            for (int round = 0; round < 3; round++)
              for (int step = 0; step < targets.size(); step++) {
                int index = (offset + step) % targets.size();
                same(expected.get(index % TARGETS.length), query.matchAll(targets.get(index)),
                    pattern + " on target " + index);
              }
            return null;
          }));
        }
        start.countDown();
        for (Future<?> thread : running) thread.get(120, TimeUnit.SECONDS);
      } finally {
        threads.shutdownNow();
      }
    }
  }

  /** Later SMARTS searches of a graph reuse what earlier ones prepared and find what a fresh graph gives. */
  @Test void repeatedSmartsSearchesOfOneGraphMatchAFreshGraph() throws Exception {
    for (String smiles : TARGETS) {
      MolGraph reused = rawGraph(smiles);
      for (int round = 0; round < 3; round++)
        for (String pattern : PATTERNS) {
          Substructure query = Substructure.smarts(pattern);
          same(query.matchAll(rawGraph(smiles)), query.matchAll(reused), pattern + " in " + smiles + ", round " + round);
          same(query.match(rawGraph(smiles)), query.match(reused), pattern + " in " + smiles + ", round " + round);
        }
    }
  }

  /**
   * Searches stopped by the node limit, or ended by a failing cancellation check, leave the
   * prepared query and the graph to later searches as a fresh query and graph would be.
   */
  @Test void stoppedSmartsSearchesLeaveLaterSearchesUnchanged() throws Exception {
    String[] patterns = {"[$([R2]~[R2]),$(C(=O)[O,N]),r6]~*", "[R2;!@]~[!R]", "[x3;R]@[R]", "*~*~*"};
    for (String pattern : patterns) {
      Substructure query = Substructure.smarts(pattern);
      for (String smiles : TARGETS) {
        MolGraph reused = rawGraph(smiles);
        SearchResult expected = Substructure.smarts(pattern).matchAll(rawGraph(smiles));
        for (int limit = 1; limit <= 40; limit++) {
          String label = pattern + " in " + smiles + ", limit " + limit;
          Substructure.Options limited = new Substructure.Options();
          limited.nodeLimit = limit;
          same(Substructure.smarts(pattern).withOptions(limited).matchAll(rawGraph(smiles)),
              query.withOptions(limited).matchAll(reused), label);
          int[] polls = {0};
          int failAt = limit;
          Substructure.Options failing = new Substructure.Options();
          failing.cancellation = () -> {
            if (++polls[0] >= failAt) throw new IllegalStateException("cancellation check failed");
            return false;
          };
          try {
            query.withOptions(failing).matchAll(reused);
          } catch (IllegalStateException failure) {
            assertEquals("cancellation check failed", failure.getMessage(), label);
          }
          same(expected, query.matchAll(reused), label);
        }
      }
    }
  }

  /**
   * Ring counts, alone and beside a negated stereo test, find what SmartsPattern finds in fused,
   * bridged and spiro rings.
   */
  @Test void ringCountTestsFindWhatTheReferenceMatcherFinds() throws Exception {
    String[] targets = {"C1CCC2CCCCC2C1", "C1CC2CCC1C2", "C1CCC2(CC1)CCCC2", "C1CCC2CC3CCCCC3CC2C1",
        "c1ccc2c(c1)ccc1ccccc12", "CCCCC"};
    String[] patterns = {"[R1]", "[R2]", "[C;R2;!@]", "[R2]~[R2]", "[!R2;R]", "[$(*~[R2])]"};
    SmilesParser parser = new SmilesParser(SilentChemObjectBuilder.getInstance());
    for (String pattern : patterns)
      for (String smiles : targets) {
        IAtomContainer target = parser.parseSmiles(smiles);
        int expected = SmartsPattern.create(pattern).matchAll(target.clone()).count();
        SearchResult result = Substructure.smarts(pattern).matchAll(new MolGraph(target));
        assertEquals(SearchResult.Status.Completed, result.status(), pattern + " in " + smiles);
        assertEquals(expected, result.matches().size(), pattern + " in " + smiles);
      }
  }

  /** A batch fails with the failure of its lowest failing target, whatever the number of threads. */
  @Test void batchReportsTheLowestFailingTarget() throws Exception {
    List<MolGraph> targets = new ArrayList<>();
    for (int index = 0; index < 40; index++) targets.add(rawGraph(index % 2 == 0 ? "CCO" : "c1ccccc1"));
    Map<MolGraph, Integer> indices = new IdentityHashMap<>();
    for (int index = 0; index < targets.size(); index++) indices.put(targets.get(index), index);
    Substructure query = Substructure.find(rawGraph("CO"));
    for (int threads : new int[] {1, 2, 4, 8}) {
      IllegalStateException failure = assertThrows(IllegalStateException.class,
          () -> SearchResult.searchEach(targets, threads, target -> {
            int index = indices.get(target);
            // Slow enough that other threads take part after the head start.
            for (long until = System.nanoTime() + 300_000L; System.nanoTime() < until;) continue;
            if (index == 31 || index == 7 || index == 23) throw new IllegalStateException("target " + index);
            if (index == 11) throw new AssertionError("target 11");
            return query.match(target);
          }));
      assertEquals("target 7", failure.getMessage(), threads + " threads");
    }
    for (int threads : new int[] {1, 2, 4, 8}) {
      List<SearchResult> results = SearchResult.searchEach(targets, threads, query::match);
      for (int index = 0; index < targets.size(); index++)
        assertEquals(index % 2 == 0, results.get(index).found(), threads + " threads, target " + index);
    }
  }
}
