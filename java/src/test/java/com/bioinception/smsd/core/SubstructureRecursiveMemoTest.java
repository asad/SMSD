/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
package com.bioinception.smsd.core;

import static com.bioinception.smsd.TestSupport.rawGraph;
import static com.bioinception.smsd.TestSupport.rawMolecule;
import static org.junit.jupiter.api.Assertions.*;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.CsvSource;
import org.openscience.cdk.smarts.SmartsPattern;

class SubstructureRecursiveMemoTest {
  @ParameterizedTest
  @CsvSource(value = {
      "[$(C(=O)[O,N])]~*|CC(=O)Nc1ccc(O)cc1",
      "[C;!$(C=O)]~*|CC(=O)NCC(=O)O",
      "[$([C;$(C(=O)N)])]~*|CC(=O)NCC(=O)N",
      "[$([R2]~[R2]),$(C(=O)N)]~*|CC(=O)NC1CCC2CCCCC2C1",
      "[$(N[C@@H](C)C(=O)O)]|N[C@@H](C)C(=O)O",
      "[$(N[C@@H](C)C(=O)O)]|N[C@H](C)C(=O)O",
      "[$(F/C=C/Cl)]|F/C=C/Cl",
      "[$(F/C=C/Cl)]|F/C=C\\Cl"
  }, delimiter = '|')
  void keepsRecursiveBooleanRingAndStereoPredicates(String pattern, String smiles) throws Exception {
    int expected = SmartsPattern.create(pattern).matchAll(rawMolecule(smiles)).count();
    SearchResult actual = Substructure.smarts(pattern).matchAll(rawGraph(smiles));
    assertEquals(SearchResult.Status.Completed, actual.status());
    assertTrue(actual.complete());
    assertEquals(expected, actual.matches().size());
  }

  @Test void scopesAnswersToOneTargetAndOneCall() throws Exception {
    Substructure query = Substructure.smarts("[$(C=O)]~*");
    for (String smiles : Arrays.asList("CC=O", "CCC", "CC=O", "COC", "CC=O")) {
      SearchResult expected = Substructure.smarts("[$(C=O)]~*").matchAll(rawGraph(smiles));
      SearchResult actual = query.matchAll(rawGraph(smiles));
      assertEquals(expected.matches(), actual.matches());
      assertEquals(expected.nodesVisited(), actual.nodesVisited());
      assertTrue(actual.complete());
    }
  }

  @Test void aStoppedCallLeavesThePreparedQueryReusable() throws Exception {
    Substructure query = Substructure.smarts("[$(C(=O)N)]~*~*");
    MolGraph target = rawGraph("CC(=O)NCC(=O)NCC(=O)N");
    SearchResult expected = query.matchAll(target);
    Substructure.Options limited = new Substructure.Options();
    limited.nodeLimit = 3;
    SearchResult stopped = query.withOptions(limited).matchAll(target);
    assertEquals(SearchResult.Status.LimitReached, stopped.status());
    assertFalse(stopped.complete());
    assertEquals(expected.matches(), query.matchAll(target).matches());
    Substructure.Options cancelled = new Substructure.Options();
    cancelled.cancellation = () -> true;
    assertEquals(SearchResult.Status.Cancelled, query.withOptions(cancelled).match(target).status());
    assertEquals(expected.matches(), query.matchAll(target).matches());
  }

  @Test void concurrentCallsKeepTheirOwnRecursiveAnswers() throws Exception {
    Substructure query = Substructure.smarts("[$(C(=O)[N,O])]~*");
    List<MolGraph> targets = Arrays.asList(rawGraph("CC(=O)N"), rawGraph("CCCC"), rawGraph("CC(=O)O"));
    List<SearchResult> expected = new ArrayList<>();
    for (MolGraph target : targets) expected.add(query.matchAll(target));
    ExecutorService threads = Executors.newFixedThreadPool(2);
    try {
      List<Future<SearchResult>> pending = new ArrayList<>();
      for (int call = 0; call < 12; call++) {
        MolGraph target = targets.get(call % targets.size());
        pending.add(threads.submit(() -> query.matchAll(target)));
      }
      for (int call = 0; call < pending.size(); call++) {
        SearchResult actual = pending.get(call).get(10, TimeUnit.SECONDS);
        assertEquals(expected.get(call % targets.size()).matches(), actual.matches());
        assertTrue(actual.complete());
      }
    } finally {
      threads.shutdownNow();
    }
  }
}
