/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
package com.bioinception.smsd.core;

import static com.bioinception.smsd.TestSupport.rawGraph;
import static org.junit.jupiter.api.Assertions.*;
import com.fasterxml.jackson.databind.ObjectMapper;
import java.util.*;
import java.util.concurrent.atomic.AtomicInteger;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.CsvSource;

/** Search order, node counts, match construction and per-search work arrays of molecule queries. */
class SubstructureSearchInternalsTest {
  private static final ObjectMapper JSON = new ObjectMapper();

  /** Node counts and match counts of first and all-match searches, fixed by the search order. */
  @ParameterizedTest
  @CsvSource({
      "CC.O, OCC.CC, 3, 10, 4",
      "C.N.O, NCO.CN.OC, 3, 19, 12",
      "C.CC, CCC.CC, 3, 30, 18",
      "c1ccccc1, c1ccc2ccccc2c1, 13, 129, 24",
      "CC(C)(C)C, CC(C)(C)CC(C)(C)C, 5, 131, 48"
  })
  void searchOrderKeepsItsNodeCounts(String query, String target, long firstNodes, long allNodes, int matches)
      throws Exception {
    Substructure prepared = Substructure.find(rawGraph(query));
    SearchResult first = prepared.match(rawGraph(target)), all = prepared.matchAll(rawGraph(target));
    assertEquals(firstNodes, first.nodesVisited());
    assertEquals(allNodes, all.nodesVisited());
    assertEquals(matches, all.matches().size());
    assertEquals(SearchResult.Status.Completed, all.status());
  }

  @Test
  void firstMatchesKeepTheirMapping() throws Exception {
    SearchResult result = Substructure.find(rawGraph("c1ccccc1")).match(rawGraph("c1ccc2ccccc2c1"));
    SearchResult.Match match = result.matches().get(0);
    assertEquals("{0=0, 1=1, 2=2, 3=3, 4=8, 5=9}", match.atomMapping().toString());
    assertEquals("{0-1=0-1, 0-5=0-9, 1-2=1-2, 2-3=2-3, 3-4=3-8, 4-5=8-9}", match.bondMapping().toString());
    assertEquals("[[0, 1, 2, 3, 4, 5]]", match.fragments().toString());
  }

  /** Lines of the atom mapping and fragments of every match, in result order. */
  private static String describe(SearchResult result) {
    StringBuilder text = new StringBuilder();
    for (SearchResult.Match match : result.matches())
      text.append(match.atomMapping()).append(" ").append(match.fragments()).append("\n");
    return text.toString();
  }

  @Test
  void fragmentsOfEachEmbeddingFollowTheOrderTheSearchReachedThem() throws Exception {
    // A first match lists fragments in query order; each enumerated match, in the order the
    // search mapped them, here the oxygen first since it has the fewest targets.
    Substructure query = Substructure.find(rawGraph("CC.O"));
    MolGraph target = rawGraph("OCC.CC");
    assertEquals("{0=2, 1=1, 2=0} [[0, 1], [2]]\n", describe(query.match(target)));
    assertEquals("{0=2, 1=1, 2=0} [[2], [0, 1]]\n"
        + "{0=3, 1=4, 2=0} [[2], [0, 1]]\n"
        + "{0=4, 1=3, 2=0} [[2], [0, 1]]\n"
        + "{0=1, 1=2, 2=0} [[2], [0, 1]]\n", describe(query.matchAll(target)));
    SearchResult three = Substructure.find(rawGraph("C.N.O")).matchAll(rawGraph("NCO.CN.OC"));
    assertEquals(12, three.matches().size());
    for (SearchResult.Match match : three.matches()) assertEquals("[[1], [2], [0]]", match.fragments().toString());
    String expected = "";
    for (String mapping : Arrays.asList("2, 0, 1", "3, 0, 1", "4, 0, 1", "0, 2, 1", "3, 2, 1", "4, 2, 1",
        "0, 3, 4", "2, 3, 4", "1, 3, 4", "0, 4, 3", "2, 4, 3", "1, 4, 3",
        "2, 1, 0", "3, 1, 0", "4, 1, 0", "0, 1, 2", "3, 1, 2", "4, 1, 2")) {
      String[] targets = mapping.split(", ");
      expected += "{0=" + targets[0] + ", 1=" + targets[1] + ", 2=" + targets[2] + "} [[1, 2], [0]]\n";
    }
    assertEquals(expected, describe(Substructure.find(rawGraph("C.CC")).matchAll(rawGraph("CCC.CC"))));
  }

  @Test
  void ownedMatchesBehaveAsCopiedMatches() throws Exception {
    SearchResult.Bond first = new SearchResult.Bond(0, 1), second = new SearchResult.Bond(1, 2),
        third = new SearchResult.Bond(0, 2);
    SearchResult.QueryBonds bonds = new SearchResult.QueryBonds(Arrays.asList(first, second, third));
    List<Set<Integer>> fragments = Collections.unmodifiableList(Arrays.asList(
        Collections.unmodifiableSet(new TreeSet<>(Arrays.asList(0, 1, 2)))));
    SearchResult.Match owned = SearchResult.Match.of(new int[] {5, 3, 9}, bonds, fragments, 3);
    Map<Integer, Integer> atoms = new LinkedHashMap<>();
    atoms.put(2, 9);
    atoms.put(0, 5);
    atoms.put(1, 3);
    Map<SearchResult.Bond, SearchResult.Bond> images = new LinkedHashMap<>();
    images.put(first, new SearchResult.Bond(5, 3));
    images.put(second, new SearchResult.Bond(3, 9));
    images.put(third, new SearchResult.Bond(5, 9));
    SearchResult.Match copied = new SearchResult.Match(atoms, images, fragments, 3);

    assertEquals(copied, owned);
    assertEquals(owned, copied);
    assertEquals(copied.hashCode(), owned.hashCode());
    assertEquals(copied.atomMapping().hashCode(), owned.atomMapping().hashCode());
    assertEquals(copied.bondMapping().hashCode(), owned.bondMapping().hashCode());
    assertEquals(copied.atomMapping().toString(), owned.atomMapping().toString());
    assertEquals(copied.bondMapping().toString(), owned.bondMapping().toString());
    assertEquals(new ArrayList<>(copied.atomMapping().entrySet()), new ArrayList<>(owned.atomMapping().entrySet()));
    assertEquals(new ArrayList<>(copied.bondMapping().keySet()), new ArrayList<>(owned.bondMapping().keySet()));
    assertEquals(new ArrayList<>(copied.bondMapping().values()), new ArrayList<>(owned.bondMapping().values()));
    assertEquals(JSON.writeValueAsString(copied), JSON.writeValueAsString(owned));
    assertEquals(owned, JSON.readValue(JSON.writeValueAsString(owned), SearchResult.Match.class));
    assertEquals(3, owned.atomCount());
    assertEquals(3, owned.bondCount());

    assertEquals(Integer.valueOf(3), owned.atomMapping().get(1));
    assertNull(owned.atomMapping().get(3));
    assertNull(owned.atomMapping().get(-1));
    assertNull(owned.atomMapping().get("1"));
    assertTrue(owned.atomMapping().containsKey(2));
    assertFalse(owned.atomMapping().containsKey(3));
    assertTrue(owned.atomMapping().containsValue(9));
    assertEquals(new SearchResult.Bond(3, 9), owned.bondMapping().get(new SearchResult.Bond(2, 1)));
    assertNull(owned.bondMapping().get(new SearchResult.Bond(1, 3)));
    assertNull(owned.bondMapping().get("0-1"));
    assertTrue(owned.bondMapping().containsKey(third));
    assertThrows(UnsupportedOperationException.class, () -> owned.atomMapping().put(0, 0));
    assertThrows(UnsupportedOperationException.class, () -> owned.bondMapping().remove(first));
    assertThrows(UnsupportedOperationException.class,
        () -> owned.atomMapping().entrySet().iterator().next().setValue(1));
    Iterator<Map.Entry<Integer, Integer>> entries = owned.atomMapping().entrySet().iterator();
    for (int entry = 0; entry < 3; entry++) entries.next();
    assertThrows(NoSuchElementException.class, entries::next);

    SearchResult.Match empty = SearchResult.Match.of(new int[0],
        new SearchResult.QueryBonds(Collections.emptyList()), Collections.emptyList(), 0);
    assertEquals(new SearchResult.Match(Collections.emptyMap(), Collections.emptyMap(), Collections.emptyList(), 0), empty);
    assertEquals("{}", empty.atomMapping().toString());
  }

  @Test
  void searchResultsSurviveTheirJsonForm() throws Exception {
    SearchResult result = Substructure.find(rawGraph("CC.O")).matchAll(rawGraph("OCC.CC"));
    SearchResult restored = JSON.readValue(JSON.writeValueAsString(result), SearchResult.class);
    assertEquals(result.matches(), restored.matches());
    assertEquals(JSON.writeValueAsString(result), JSON.writeValueAsString(restored));
  }

  @Test
  void neighbourhoodLabelsFromReusedWorkArraysEqualFreshOnes() throws Exception {
    MolGraph query = rawGraph("CCO");
    for (String smiles : Arrays.asList("c1ccc2ccccc2c1", "CC(C)(C)CC(=O)NC1CCCC1", "OCC.CC", "C")) {
      MolGraph target = rawGraph(smiles);
      Substructure.SubgraphSearch search = new Substructure.SubgraphSearch(query,
          new Substructure.QueryTables(query, new ChemOptions()), target, new ChemOptions(),
          new SearchResult.Control(60_000, null));
      for (int pass = 0; pass < 2; pass++)
        for (int atom = target.atomCount() - 1; atom >= 0; atom--) {
          assertArrayEquals(MolGraph.buildNLF3(target, atom), search.neighbourhoodLabels(atom, 3), smiles);
          assertArrayEquals(MolGraph.buildNLF2(target, atom), search.neighbourhoodLabels(atom, 2), smiles);
        }
    }
    // Stamps are compared for equality only, so they may run through negative values.
    MolGraph graph = rawGraph("CC(C)(C)CC(=O)NC1CCCC1");
    int[] seen = new int[graph.atomCount()], reached = new int[graph.atomCount()];
    int stamp = Integer.MAX_VALUE - 3;
    for (int atom = 0; atom < graph.atomCount(); atom++, stamp++)
      assertArrayEquals(MolGraph.buildNLF2(graph, atom), MolGraph.buildNeighborhoodNLF(graph, atom, 2, seen, stamp, reached));
  }

  @Test
  void sampledStopsReadCancellationOnEveryCallAndTheClockOnEverySixteenth() {
    AtomicInteger checks = new AtomicInteger();
    // A zero time limit has passed at once, but only the sixteenth sampled call reads the clock.
    SearchResult.Control expired = new SearchResult.Control(0, () -> checks.incrementAndGet() < 0);
    for (int call = 1; call < 16; call++) assertFalse(expired.stoppedSampled());
    assertTrue(expired.stoppedSampled());
    assertEquals(16, checks.get());
    assertEquals(SearchResult.Status.TimedOut, expired.status());
    SearchResult.Control cancelled = new SearchResult.Control(60_000, () -> true);
    assertTrue(cancelled.stoppedSampled());
    assertEquals(SearchResult.Status.Cancelled, cancelled.status());
    SearchResult.Control unsampled = new SearchResult.Control(0, null);
    assertTrue(unsampled.stopped());
    assertEquals(SearchResult.Status.TimedOut, unsampled.status());
  }
}
