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
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.ValueSource;

class SubstructureEnumerationTest {
  private static Substructure disconnectedQuery(boolean smarts) throws Exception {
    return smarts ? Substructure.smarts("[C]-[C].[O]") : Substructure.find(rawGraph("CC.O"));
  }

  private static Set<List<Integer>> embeddings(SearchResult result) {
    Set<List<Integer>> embeddings = new HashSet<>();
    for (SearchResult.Match match : result.matches())
      embeddings.add(new ArrayList<>(match.atomMapping().values()));
    return embeddings;
  }

  private static Set<List<Integer>> expectedEmbeddings() {
    Set<List<Integer>> expected = new HashSet<>();
    for (int[] bond : new int[][] {{0, 1}, {1, 0}, {3, 4}, {4, 3}})
      for (int oxygen : new int[] {2, 5}) expected.add(Arrays.asList(bond[0], bond[1], oxygen));
    return expected;
  }

  @ParameterizedTest
  @ValueSource(booleans = {false, true})
  void retainsEveryEmbeddingAndItsSelectedGraph(boolean smarts) throws Exception {
    Substructure query = disconnectedQuery(smarts);
    MolGraph target = rawGraph("CC.O.CC.O");
    SearchResult all = query.matchAll(target);
    assertEquals(SearchResult.Status.Completed, all.status());
    assertTrue(all.complete());
    assertEquals(expectedEmbeddings(), embeddings(all));
    assertEquals(8, all.matches().size());
    Set<Set<Integer>> fragments = new HashSet<>(Arrays.asList(
        new TreeSet<>(Arrays.asList(0, 1)), Collections.singleton(2)));
    for (SearchResult.Match match : all.matches()) {
      assertEquals(1, match.bondCount());
      SearchResult.Bond mapped = new SearchResult.Bond(match.atomMapping().get(0), match.atomMapping().get(1));
      assertEquals(Collections.singletonMap(new SearchResult.Bond(0, 1), mapped), match.bondMapping());
      assertEquals(fragments, new HashSet<>(match.fragments()));
      assertThrows(UnsupportedOperationException.class, () -> match.fragments().get(0).add(9));
    }
    SearchResult limited = query.matchAll(target, 3);
    assertEquals(SearchResult.Status.LimitReached, limited.status());
    assertFalse(limited.complete());
    assertEquals(3, limited.matches().size());
    assertTrue(expectedEmbeddings().containsAll(embeddings(limited)));
    assertEquals(expectedEmbeddings(), embeddings(query.matchAll(target)));
  }

  @ParameterizedTest
  @ValueSource(booleans = {false, true})
  void enumerationHasIndependentConcurrentResults(boolean smarts) throws Exception {
    Substructure query = disconnectedQuery(smarts);
    MolGraph target = rawGraph("CC.O.CC.O");
    ExecutorService workers = Executors.newFixedThreadPool(3);
    try {
      List<Future<SearchResult>> calls = new ArrayList<>();
      for (int call = 0; call < 6; call++) calls.add(workers.submit(() -> query.matchAll(target)));
      for (Future<SearchResult> call : calls) {
        SearchResult result = call.get(2, TimeUnit.SECONDS);
        assertEquals(SearchResult.Status.Completed, result.status());
        assertEquals(expectedEmbeddings(), embeddings(result));
        assertEquals(2, result.matches().get(0).fragments().size());
      }
    } finally { workers.shutdownNow(); }
  }
}
