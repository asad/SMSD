/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
package com.bioinception.smsd.core;

import static com.bioinception.smsd.core.RingSearchFixtures.*;
import static org.junit.jupiter.api.Assertions.*;

import com.bioinception.smsd.core.RingCandidateSearch.PairRecord;
import com.bioinception.smsd.core.RingCandidateSearch.RingSearchNotes;
import com.bioinception.smsd.core.RingCandidateSearch.RingSearchOptions;
import com.fasterxml.jackson.databind.JsonNode;
import com.fasterxml.jackson.databind.ObjectMapper;
import com.fasterxml.jackson.databind.node.ArrayNode;
import com.fasterxml.jackson.databind.node.ObjectNode;
import java.io.BufferedReader;
import java.io.InputStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Deque;
import java.util.List;
import java.util.Map;
import java.util.Random;
import org.junit.jupiter.api.Tag;
import org.junit.jupiter.api.Test;

/** The ring candidates read from one tree per atom, against the pairwise search they stand for. */
class RingCandidateSearchTest {
  private static final long LIMIT = RingPerception.RING_SEARCH_WORK_LIMIT;
  private static final RingSearchOptions DEFAULTS = RingSearchOptions.DEFAULTS;
  private static final RingSearchOptions LOOSE = new RingSearchOptions(false, false, true, false);
  private static final RingSearchOptions REPLAY_ALL = new RingSearchOptions(false, false, false, true);
  private static final RingSearchOptions LOOSE_REPLAY = new RingSearchOptions(false, false, true, true);
  private static final RingSearchOptions LIFT = new RingSearchOptions(true, false, false, false);
  private static final RingSearchOptions PROVE = new RingSearchOptions(false, true, false, false);

  private static RingSearchNotes traced() {
    RingSearchNotes notes = new RingSearchNotes();
    notes.trace = new ArrayList<>();
    return notes;
  }

  /** Small graphs for every run: the small fixtures, bespoke graphs and seeded random graphs. */
  private static List<Object[]> corpus() throws Exception {
    List<Object[]> graphs = new ArrayList<>();
    for (Map.Entry<String, String> fixture : smiles().entrySet())
      if (small(fixture.getKey())) graphs.add(new Object[] {fixture.getKey(), neighborsOf(fixture.getValue())});
    graphs.add(new Object[] {"extended class", extendedClassGraph()});
    graphs.add(new Object[] {"eight ring", eightRingGraph()});
    for (int layers = 2; layers <= 8; layers++) graphs.add(new Object[] {"stack " + layers, prismStack(layers)});
    graphs.add(new Object[] {"cyclobutane 5", cyclobutaneChain(5)});
    graphs.add(new Object[] {"double hub 12", doubleHub(12)});
    Random random = new Random(20261030L);
    for (int i = 0; i < 120; i++) graphs.add(new Object[] {"random " + i, randomRingGraph(random)});
    return graphs;
  }

  /** The number of a bond in edge order: lower atom first, then neighbour order. */
  private static int bondIndex(int[][] neighbors, int a, int b) {
    int index = 0;
    for (int i = 0; i < neighbors.length; i++)
      for (int j : neighbors[i])
        if (j > i) {
          if (i == Math.min(a, b) && j == Math.max(a, b)) return index;
          index++;
        }
    return -1;
  }

  @Test void everyPairReadsThePathsAndTheChargeOfThePairwiseSearch() throws Exception {
    long[] seen = new long[6];
    for (Object[] item : corpus()) {
      int[][] graph = (int[][]) item[1];
      for (int[][] neighbors : new int[][][] {graph, sortedNeighbours(graph)}) {
        int n = neighbors.length;
        RingPerceptionReference.Replay replay = RingPerceptionReference.replay(n, neighbors);
        RingSearchNotes notes = traced();
        RingPerception.perceive(n, neighbors, false, Long.MAX_VALUE, DEFAULTS, notes);
        assertEquals(replay.pairs.size(), notes.trace.size(), (String) item[0]);
        for (int k = 0; k < notes.trace.size(); k++) {
          PairRecord pair = notes.trace.get(k);
          String label = item[0] + " pair " + k;
          assertEquals(replay.pairs.get(k)[0], pair.root, label);
          assertEquals(replay.pairs.get(k)[1], pair.bond, label);
          assertEquals(replay.valid.get(k), pair.valid, label);
          long charge = replay.searchCharge.get(k);
          if (pair.bfsLow == pair.bfsHigh) {
            assertEquals(charge, pair.bfsLow, label);
          } else {
            assertEquals(RingCandidateSearch.PAIR_B2, pair.kind, label);
            assertTrue(pair.bfsLow <= charge && charge <= pair.bfsHigh, label);
          }
          if (pair.kind != RingCandidateSearch.PAIR_SIMPLE_RING) {
            assertArrayEquals(replay.firstPaths.get(k), pair.firstPath, label);
            assertArrayEquals(replay.secondPaths.get(k), pair.secondPath, label);
          }
          assertEquals(replay.firstPaths.get(k).length + replay.secondPaths.get(k).length, pair.pathLength, label);
          seen[pair.kind]++;
        }
      }
    }
    for (int kind = 0; kind < seen.length; kind++) assertTrue(seen[kind] > 0, "pair class " + kind);
  }

  @Test void aBypassFromOneLevelDeeperIsChargedExactly() {
    int[][] neighbors = extendedClassGraph();
    int n = neighbors.length, bond = bondIndex(neighbors, 0, 2);
    assertArrayEquals(new int[] {4, 2}, RingCandidateSearch.coversForTesting(n, neighbors, 0)[2]);
    RingSearchNotes notes = traced();
    RingPerception.Rings rings = RingPerception.perceive(n, neighbors, false, LIMIT, DEFAULTS, notes);
    PairRecord pair = null;
    for (PairRecord record : notes.trace) if (record.root == 0 && record.bond == bond) pair = record;
    assertNotNull(pair);
    assertEquals(RingCandidateSearch.PAIR_B1_PLUS, pair.kind);
    assertArrayEquals(new int[] {0}, pair.firstPath);
    assertArrayEquals(new int[] {0, 1, 4, 2}, pair.secondPath);
    assertEquals(15, pair.bfsLow);
    assertEquals(15, pair.bfsHigh);
    RingPerceptionReference.Replay replay = RingPerceptionReference.replay(n, neighbors);
    for (int k = 0; k < replay.pairs.size(); k++)
      if (replay.pairs.get(k)[0] == 0 && replay.pairs.get(k)[1] == bond) assertEquals(15L, replay.searchCharge.get(k));
    boolean found = false;
    for (int[] ring : rings.smallest) found |= Arrays.equals(new int[] {0, 2, 4, 1}, ring);
    assertTrue(found, Arrays.deepToString(rings.smallest));
  }

  @Test void everyAtomTakesTheLeastHalfBondIntoItsSubtree() throws Exception {
    long seed = 11;
    for (Object[] item : corpus())
      for (int[][] neighbors : variants((int[][]) item[1], seed++)) {
        int n = neighbors.length;
        RingPerception.Systems systems = RingPerception.systems(n, neighbors);
        if (systems == null) continue;
        for (int root = 0; root < n; root++) {
          int s = systems.ofAtom[root];
          if (s < 0 || systems.simpleRing[s]) continue;
          assertArrayEquals(leastHalfBonds(neighbors, systems.ofAtom, root),
                            RingCandidateSearch.coversForTesting(n, neighbors, root), item[0] + " root " + root);
        }
      }
  }

  /** The covers by brute force: for each atom, the least closing half-bond into its subtree. */
  private static int[][] leastHalfBonds(int[][] neighbors, int[] ofAtom, int root) {
    int n = neighbors.length, s = ofAtom[root];
    int[] par = new int[n], dep = new int[n];
    Arrays.fill(par, -2);
    List<Integer> order = new ArrayList<>();
    order.add(root);
    par[root] = -1;
    for (int head = 0; head < order.size(); head++) {
      int u = order.get(head);
      for (int w : neighbors[u])
        if (ofAtom[w] == s && par[w] == -2) { par[w] = u; dep[w] = dep[u] + 1; order.add(w); }
    }
    List<int[]> halves = new ArrayList<>();
    Deque<int[]> stack = new ArrayDeque<>();
    stack.push(new int[] {root, 0});
    while (!stack.isEmpty()) {
      int[] top = stack.peek();
      int x = top[0];
      if (top[1] == neighbors[x].length) { stack.pop(); continue; }
      int w = neighbors[x][top[1]++];
      if (ofAtom[w] != s || w == par[x]) continue;
      if (par[w] == x) stack.push(new int[] {w, 0});
      else halves.add(new int[] {x, w});
    }
    int[][] covers = new int[n][];
    for (int c : order) {
      if (c == root) continue;
      int best = Integer.MAX_VALUE;
      for (int[] half : halves) {
        if (!below(par, half[1], c) || below(par, half[0], c)) continue;
        int length = dep[half[0]] + dep[half[1]] + 1;
        if (length < best) { best = length; covers[c] = half; }
      }
    }
    return covers;
  }

  private static boolean below(int[] par, int u, int c) {
    for (int x = u; x >= 0; x = par[x]) if (x == c) return true;
    return false;
  }

  @Test void resultsEqualThePairwiseSearchForEveryOrderAndLimit() throws Exception {
    long seed = 7;
    for (Object[] item : corpus()) {
      List<int[][]> variants = variants((int[][]) item[1], seed++);
      for (int v = 0; v < variants.size(); v++) compareWithReference(variants.get(v), item[0] + " variant " + v);
    }
  }

  private static void compareWithReference(int[][] neighbors, String label) {
    int n = neighbors.length;
    RingPerceptionReference.Replay replay = RingPerceptionReference.replay(n, neighbors);
    long before = replay.beforeLast;
    for (long limit : new long[] {LIMIT, before, before - 1})
      for (boolean enumerate : new boolean[] {false, true}) {
        RingSearchNotes notes = new RingSearchNotes();
        RingPerception.Rings rings = RingPerception.perceive(n, neighbors, enumerate, limit, DEFAULTS, notes);
        String where = label + " limit " + limit + " enumerate " + enumerate;
        assertSameRings(RingPerceptionReference.perceive(n, neighbors, enumerate, limit), rings, where);
        if (replay.pairs.isEmpty()) continue;
        assertEquals(replay.fallsBack(limit), notes.fallback, where);
        assertEquals(replay.pairs.size(), notes.pairs, where);
        if (notes.stoppedAfterPair < 0) assertEquals(replay.pairs.size() - 1, notes.pairsSummed, where);
      }
  }

  @Test void keysReadAlongThePathsGiveTheSameResults() throws Exception {
    for (Object[] item : corpus()) {
      int[][] neighbors = (int[][]) item[1];
      int n = neighbors.length;
      for (boolean enumerate : new boolean[] {false, true}) {
        RingSearchNotes kept = new RingSearchNotes(), walked = new RingSearchNotes();
        RingPerception.Rings expected = RingPerception.perceive(n, neighbors, enumerate, LIMIT, DEFAULTS, kept);
        RingPerception.Rings rings = RingCandidateSearch.perceiveWalkingKeysForTesting(n, neighbors, enumerate, LIMIT, walked);
        assertSameRings(expected, rings, (String) item[0]);
        assertEquals(kept.candidates, walked.candidates, (String) item[0]);
        assertEquals(kept.fallback, walked.fallback, (String) item[0]);
      }
    }
  }

  /** A seeded graph whose last pair, never charged, is a bracketed one. */
  private static int[][] lastPairBracketedGraph() {
    Random random = new Random(42);
    for (int attempt = 0; attempt < 20000; attempt++) {
      int[][] neighbors = randomRingGraph(random);
      RingSearchNotes notes = traced();
      RingPerception.perceive(neighbors.length, neighbors, false, Long.MAX_VALUE, DEFAULTS, notes);
      if (notes.trace.isEmpty()) continue;
      PairRecord last = notes.trace.get(notes.trace.size() - 1);
      if (last.kind == RingCandidateSearch.PAIR_B2 && last.bfsHigh - last.bfsLow >= 4) return neighbors;
    }
    fail("no seeded graph has a bracketed last pair");
    return null;
  }

  @Test void theFallbackFollowsTheWorkBeforeTheLastPair() throws Exception {
    List<Object[]> graphs = new ArrayList<>();
    graphs.add(new Object[] {"CD8", neighborsOf(cyclodextrin(8))});
    graphs.add(new Object[] {"cubane", neighborsOf("C12C3C4C1C5C2C3C45")});
    graphs.add(new Object[] {"adamantane", neighborsOf("C1C2CC3CC1CC(C2)C3")});
    for (int layers = 3; layers <= 8; layers++) graphs.add(new Object[] {"stack " + layers, prismStack(layers)});
    graphs.add(new Object[] {"cyclobutane 5", cyclobutaneChain(5)});
    graphs.add(new Object[] {"bracketed last pair", lastPairBracketedGraph()});
    for (Object[] item : graphs) {
      int[][] neighbors = (int[][]) item[1];
      int n = neighbors.length;
      RingPerceptionReference.Replay replay = RingPerceptionReference.replay(n, neighbors);
      long before = replay.beforeLast;
      for (RingSearchOptions options : new RingSearchOptions[] {DEFAULTS, LOOSE, REPLAY_ALL, LOOSE_REPLAY})
        for (long limit : new long[] {before, before - 1}) {
          String label = item[0] + " limit " + limit + " loose " + options.loosenBracketsForTesting
              + " replay " + options.replayAllForTesting;
          RingSearchNotes notes = new RingSearchNotes();
          RingPerception.Rings rings = RingPerception.perceive(n, neighbors, false, limit, options, notes);
          assertSameRings(RingPerceptionReference.perceive(n, neighbors, false, limit), rings, label);
          assertEquals(limit < before, notes.fallback, label);
          assertTrue(notes.pairsSummed <= replay.pairs.size() - 1, label);
          if (notes.stoppedAfterPair < 0) assertEquals(replay.pairs.size() - 1, notes.pairsSummed, label);
          if (options.replayAllForTesting || (options.loosenBracketsForTesting && !notes.fallback)) {
            assertTrue(notes.baseWorkExact, label);
            assertEquals(before, notes.baseWorkLow, label);
            assertEquals(before, notes.baseWorkHigh, label);
          }
          if (notes.baseWorkExact) assertEquals(before, notes.baseWorkLow, label);
          if (!notes.baseWorkHighUnbounded) assertTrue(notes.baseWorkLow <= before && before <= notes.baseWorkHigh, label);
        }
    }
  }

  @Test void certificateHookValuesKeepTheirResults() throws Exception {
    for (String text : new String[] {"c1ccc2ccccc2c1", "C1C2CC3CC1CC(C2)C3", "C12C3C4C1C5C2C3C45",
                                     "C12C3C4C5C6C1C7C2C3C4C5C67", cyclodextrin(8)}) {
      int[][] neighbors = neighborsOf(text);
      for (long limit : new long[] {0, 1L << 34})
        for (boolean enumerate : new boolean[] {false, true})
          assertSameRings(RingPerceptionReference.perceive(neighbors.length, neighbors, enumerate, limit),
                          RingPerception.perceive(neighbors.length, neighbors, enumerate, limit), text + " " + limit);
    }
    int[][] cd34 = neighborsOf(cyclodextrin(34));
    RingPerception.Rings rings = RingPerception.perceive(cd34.length, cd34, false, 1L << 20);
    assertTrue(rings.truncated);
    assertFalse(rings.basisTruncated);
    assertEquals(35, rings.smallest.length);
    assertEquals(34 * 11, weight(rings.smallest));
    assertSameRings(RingPerceptionReference.perceive(cd34.length, cd34, false, 1L << 20), rings, "CD34");
  }

  @Test void aGraphListingANeighbourTwiceKeepsThePairwiseSearch() {
    int[][] neighbors = {{1, 1, 2}, {0, 0, 2, 3}, {0, 1, 4}, {1, 4, 5}, {2, 3, 5}, {3, 4}};
    assertNull(RingPerception.systems(neighbors.length, neighbors));
    for (long limit : new long[] {LIMIT, 0})
      for (boolean enumerate : new boolean[] {false, true})
        assertSameRings(RingPerceptionReference.perceive(neighbors.length, neighbors, enumerate, limit),
                        RingPerception.perceive(neighbors.length, neighbors, enumerate, limit), "repeated neighbour");
  }

  @Test void theSwitchesAreOffForEveryGraph() {
    assertSame(RingSearchOptions.DEFAULTS, MolGraph.RING_SEARCH_OPTIONS);
    assertFalse(DEFAULTS.liftRefusedBasis);
    assertFalse(DEFAULTS.proveRelevantByFamilies);
    assertFalse(DEFAULTS.loosenBracketsForTesting);
    assertFalse(DEFAULTS.replayAllForTesting);
  }

  @Test void aRefusedBasisIsCompletedOnlyWhenAskedAndProved() {
    int[][] neighbors = prismStack(80);
    int n = neighbors.length, rank = RingPerception.cycleRank(n, neighbors);
    RingSearchNotes plain = new RingSearchNotes();
    RingPerception.Rings refused = RingPerception.perceive(n, neighbors, false, LIMIT, DEFAULTS, plain);
    assertTrue(plain.fallback);
    assertTrue(plain.baseBasisTruncated);
    assertTrue(refused.basisTruncated);
    assertFalse(refused.basisLifted);
    RingSearchNotes notes = new RingSearchNotes();
    RingPerception.Rings lifted = RingPerception.perceive(n, neighbors, false, LIMIT, LIFT, notes);
    assertTrue(notes.baseBasisTruncated);
    assertTrue(notes.basisLifted);
    assertTrue(lifted.basisLifted);
    assertFalse(lifted.basisTruncated);
    assertTrue(notes.searchWork <= LIMIT);
    assertEquals(rank, lifted.smallest.length);
    assertTrue(RingPerception.certifiedMinimumBasis(n, neighbors, lifted.smallest, rank, LIMIT));
    assertArrayEquals(RingPerception.perceive(n, neighbors, false, Long.MAX_VALUE).smallest, lifted.smallest);
    assertArrayEquals(refused.smallest, lifted.baseSmallest);
    assertArrayEquals(refused.relevant, lifted.relevant);
    assertArrayEquals(refused.legacyRelevant, lifted.legacyRelevant);
    assertEquals(refused.relevantTruncated, lifted.relevantTruncated);
    assertEquals(refused.truncated, lifted.truncated);
    RingPerception.Rings relevantRefused = RingPerception.completeRelevant(n, neighbors, refused);
    RingPerception.Rings relevantLifted = RingPerception.completeRelevant(n, neighbors, lifted);
    assertArrayEquals(relevantRefused.relevant, relevantLifted.relevant);
    assertTrue(relevantLifted.relevantTruncated);
    // Without the room to search every pair the refusal stays.
    RingPerception.Rings small = RingPerception.perceive(n, neighbors, false, 1L << 20, LIFT, null);
    assertFalse(small.basisLifted);
    assertSameRings(RingPerception.perceive(n, neighbors, false, 1L << 20), small, "small budget");
  }

  @Test void familiesProveATruncatedRelevantListComplete() {
    int[][] hub = doubleHub(51);
    int n = hub.length;
    RingPerception.Rings basis = RingPerception.perceive(n, hub, false);
    RingPerception.Rings plain = RingPerception.completeRelevant(n, hub, basis);
    assertTrue(plain.relevantTruncated);
    assertEquals(51 * 50 / 2, plain.relevant.length);
    RingPerception.Rings proved = RingPerception.completeRelevant(n, hub, basis, PROVE);
    assertFalse(proved.relevantTruncated);
    assertArrayEquals(plain.relevant, proved.relevant);
    RingPerception.Rings enumerated = RingPerception.perceive(n, hub, true, LIMIT, DEFAULTS, null);
    RingPerception.Rings enumeratedProved = RingPerception.perceive(n, hub, true, LIMIT, PROVE, null);
    assertTrue(enumerated.relevantTruncated);
    assertFalse(enumeratedProved.relevantTruncated);
    assertArrayEquals(enumerated.relevant, enumeratedProved.relevant);
    // A list shorter than the members of all families stays incomplete.
    int[][] chain = cyclobutaneChain(13);
    RingPerception.Rings chainBasis = RingPerception.perceive(chain.length, chain, false);
    RingPerception.Rings chainRelevant = RingPerception.completeRelevant(chain.length, chain, chainBasis, PROVE);
    assertTrue(chainRelevant.relevantTruncated);
    assertTrue(RelevantCycleFamilySearch.compute(chain.length, chain).relevantCycleCount()
                   .compareTo(java.math.BigInteger.valueOf(chainRelevant.relevant.length)) > 0);
  }

  @Test void theKeyHashIsFixed() {
    assertEquals(0L, RingCandidateSearch.hash(new long[0], 0, 0));
    long[] words = {1L, 0x0123456789ABCDEFL, -1L};
    long h = 0;
    for (long word : words) h = (Long.rotateLeft(h, 5) ^ word) * 0x9E3779B97F4A7C15L;
    h ^= h >>> 33;
    h *= 0xFF51AFD7ED558CCDL;
    h ^= h >>> 33;
    assertEquals(h, RingCandidateSearch.hash(words, 0, 3));
    assertEquals(RingCandidateSearch.hash(new long[] {0x0123456789ABCDEFL, -1L}, 0, 2),
                 RingCandidateSearch.hash(words, 1, 2));
  }

  // ---- The golden values shared with the native core ----

  private static final String GOLDEN = "com/bioinception/smsd/ring_search_golden.json";

  /** The graphs of the golden file: the small fixtures in both neighbour orders and bespoke graphs. */
  private static List<Object[]> goldenGraphs() throws Exception {
    List<Object[]> graphs = new ArrayList<>();
    for (Map.Entry<String, String> fixture : smiles().entrySet()) {
      if (!small(fixture.getKey())) continue;
      int[][] neighbors = neighborsOf(fixture.getValue());
      graphs.add(new Object[] {fixture.getKey(), fixture.getValue(), "given", neighbors});
      graphs.add(new Object[] {fixture.getKey(), fixture.getValue(), "sorted", sortedNeighbours(neighbors)});
    }
    graphs.add(new Object[] {"extended class", "bespoke", "given", extendedClassGraph()});
    graphs.add(new Object[] {"eight ring", "bespoke", "given", eightRingGraph()});
    graphs.add(new Object[] {"stack 3", "bespoke", "given", prismStack(3)});
    graphs.add(new Object[] {"stack 8", "bespoke", "given", prismStack(8)});
    graphs.add(new Object[] {"cyclobutane 5", "bespoke", "given", cyclobutaneChain(5)});
    graphs.add(new Object[] {"double hub 12", "bespoke", "given", doubleHub(12)});
    graphs.add(new Object[] {"bracketed last pair", "bespoke", "given", lastPairBracketedGraph()});
    graphs.add(new Object[] {"six ring", "bespoke", "given", simpleRing(6)});
    graphs.add(new Object[] {"five ring", "bespoke", "given", simpleRing(5)});
    graphs.add(new Object[] {"ring beside a fused pair", "bespoke", "given", ringBesideFusedPair()});
    return graphs;
  }

  /** One simple ring of the given size. */
  private static int[][] simpleRing(int size) {
    int[][] bonds = new int[size][];
    for (int i = 0; i < size; i++) bonds[i] = new int[] {i, (i + 1) % size};
    return graph(size, bonds);
  }

  /** A five-membered ring joined by one bond to two fused six-membered rings. */
  private static int[][] ringBesideFusedPair() {
    return graph(15, new int[][] {{0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 0}, {4, 5}, {5, 6}, {6, 7}, {7, 8}, {8, 9},
                                  {9, 10}, {10, 5}, {9, 11}, {11, 12}, {12, 13}, {13, 14}, {14, 10}});
  }

  private static ArrayNode ints(ObjectMapper mapper, int[] values) {
    ArrayNode array = mapper.createArrayNode();
    for (int value : values) array.add(value);
    return array;
  }

  private static ArrayNode lists(ObjectMapper mapper, int[][] values) {
    ArrayNode array = mapper.createArrayNode();
    for (int[] value : values) array.add(ints(mapper, value));
    return array;
  }

  /** Every golden value of one graph, as this implementation finds it. */
  private static ObjectNode goldenEntry(ObjectMapper mapper, String name, String source, String order, int[][] neighbors) {
    int n = neighbors.length;
    ObjectNode entry = mapper.createObjectNode();
    entry.put("name", name);
    entry.put("source", source);
    entry.put("order", order);
    entry.set("neighbors", lists(mapper, neighbors));
    RingPerceptionReference.Replay replay = RingPerceptionReference.replay(n, neighbors);
    entry.put("pairs", replay.pairs.size());
    entry.put("workBeforeLastPair", replay.beforeLast);
    entry.put("lastPairCharge", replay.last);
    ArrayNode decisions = mapper.createArrayNode();
    // The native core has no limit below zero, so a graph without work before its last pair has
    // no limit below that work.
    long[] limits = replay.beforeLast > 0 ? new long[] {LIMIT, replay.beforeLast, replay.beforeLast - 1, 1L << 20, 0}
                                          : new long[] {LIMIT, replay.beforeLast, 1L << 20, 0};
    for (long limit : limits) {
      RingSearchNotes notes = new RingSearchNotes();
      RingPerception.perceive(n, neighbors, false, limit, DEFAULTS, notes);
      ArrayNode decision = mapper.createArrayNode();
      decision.add(limit);
      decision.add(notes.fallback);
      decisions.add(decision);
    }
    entry.set("fallback", decisions);
    RingSearchNotes notes = new RingSearchNotes();
    RingPerception.Rings rings = RingPerception.perceive(n, neighbors, false, LIMIT, DEFAULTS, notes);
    ArrayNode classes = mapper.createArrayNode();
    for (long count : new long[] {notes.nonTreePairs, notes.classAPairs, notes.classB1Pairs, notes.classB1PlusPairs,
                                  notes.classB2Pairs, notes.simpleRingPairs}) classes.add(count);
    entry.set("classes", classes);
    entry.put("baseWorkLow", notes.baseWorkLow);
    entry.put("baseWorkHigh", notes.baseWorkHigh);
    entry.put("replayedPairs", notes.replayedPairs);
    entry.put("searchWork", notes.searchWork);
    entry.put("candidates", notes.candidates);
    entry.set("basis", lists(mapper, rings.smallest));
    entry.put("legacyRelevant", rings.legacyRelevant.length);
    entry.put("truncated", rings.truncated);
    entry.put("basisTruncated", rings.basisTruncated);
    entry.put("relevantTruncated", rings.relevantTruncated);
    // The internal switches, off by default, counted alike: brackets loosened at the default limit,
    // so that every bracketed pair is replayed, and a refused basis completed at half the work.
    RingSearchNotes loosened = new RingSearchNotes();
    RingPerception.perceive(n, neighbors, false, LIMIT, LOOSE, loosened);
    entry.put("loosenedReplayedPairs", loosened.replayedPairs);
    entry.put("loosenedSearchWork", loosened.searchWork);
    RingSearchNotes lift = new RingSearchNotes();
    RingPerception.Rings lifted = RingPerception.perceive(n, neighbors, false, replay.beforeLast / 2, LIFT, lift);
    entry.put("liftedAtHalf", lift.basisLifted);
    entry.put("liftSearchWork", lift.searchWork);
    entry.set("liftBasis", lists(mapper, lifted.smallest));
    long[] work = new long[1];
    RelevantCycleFamilies families = RelevantCycleFamilySearch.compute(n, neighbors, RelevantCycleFamilySearch.FAMILY_WORK_LIMIT, work);
    ArrayNode records = mapper.createArrayNode();
    for (RelevantCycleFamily family : families.families()) {
      ObjectNode record = mapper.createObjectNode();
      record.put("size", family.size());
      record.put("root", family.root());
      record.put("first", family.first());
      record.put("second", family.second());
      record.put("opposite", family.opposite());
      record.set("prototype", ints(mapper, family.prototype()));
      record.put("memberCount", family.memberCount().toString());
      records.add(record);
    }
    entry.set("families", records);
    entry.put("relevantCycleCount", families.relevantCycleCount().toString());
    entry.put("familiesTruncated", families.truncated());
    entry.put("familyWork", work[0]);
    return entry;
  }

  private static ObjectNode goldenHashes(ObjectMapper mapper) {
    long[][] cases = {{}, {1L}, {0x0123456789ABCDEFL, -1L}, {-1L, 0L, 1L << 63, 42L}};
    ArrayNode list = mapper.createArrayNode();
    for (long[] words : cases) {
      ObjectNode item = mapper.createObjectNode();
      ArrayNode values = mapper.createArrayNode();
      for (long word : words) values.add(Long.toUnsignedString(word));
      item.set("words", values);
      item.put("hash", Long.toUnsignedString(RingCandidateSearch.hash(words, 0, words.length)));
      list.add(item);
    }
    ObjectNode node = mapper.createObjectNode();
    node.set("cases", list);
    return node;
  }

  /**
   * The golden values: the decision, charges, classes, counted work, basis and families of every
   * golden graph. With -Dsmsd.ringSearchGolden.write=path the values are written there instead.
   */
  @Test void goldenValuesAreReproduced() throws Exception {
    ObjectMapper mapper = new ObjectMapper();
    String target = System.getProperty("smsd.ringSearchGolden.write");
    List<Object[]> graphs = goldenGraphs();
    if (target != null) {
      StringBuilder text = new StringBuilder("{\n\"hashes\": ");
      text.append(mapper.writeValueAsString(goldenHashes(mapper))).append(",\n\"graphs\": [\n");
      for (int i = 0; i < graphs.size(); i++) {
        Object[] graph = graphs.get(i);
        text.append(mapper.writeValueAsString(goldenEntry(mapper, (String) graph[0], (String) graph[1], (String) graph[2],
                                                          (int[][]) graph[3])));
        text.append(i + 1 < graphs.size() ? ",\n" : "\n");
      }
      text.append("]\n}\n");
      Files.write(Paths.get(target), text.toString().getBytes(StandardCharsets.UTF_8));
      return;
    }
    JsonNode golden;
    try (InputStream in = getClass().getClassLoader().getResourceAsStream(GOLDEN)) {
      assertNotNull(in, GOLDEN);
      golden = mapper.readTree(in);
    }
    assertEquals(mapper.readTree(mapper.writeValueAsString(goldenHashes(mapper))), golden.get("hashes"));
    JsonNode entries = golden.get("graphs");
    assertEquals(graphs.size(), entries.size());
    for (int i = 0; i < entries.size(); i++) {
      JsonNode expected = entries.get(i);
      int[][] neighbors = new int[expected.get("neighbors").size()][];
      for (int atom = 0; atom < neighbors.length; atom++) {
        JsonNode list = expected.get("neighbors").get(atom);
        neighbors[atom] = new int[list.size()];
        for (int k = 0; k < list.size(); k++) neighbors[atom][k] = list.get(k).asInt();
      }
      assertArrayEquals((int[][]) graphs.get(i)[3], neighbors, expected.get("name").asText());
      JsonNode actual = mapper.readTree(mapper.writeValueAsString(
          goldenEntry(mapper, expected.get("name").asText(), expected.get("source").asText(),
                      expected.get("order").asText(), neighbors)));
      java.util.Iterator<String> fields = expected.fieldNames();
      while (fields.hasNext()) {
        String field = fields.next();
        assertEquals(expected.get(field), actual.get(field), expected.get("name").asText() + " "
                     + expected.get("order").asText() + ": " + field);
      }
    }
  }

  @Tag("SlowTest")
  @Test void largeRingSystemsKeepTheirResultsAndDecisions() throws Exception {
    for (Map.Entry<String, String> fixture : smiles().entrySet()) {
      if (small(fixture.getKey()) || fixture.getKey().contains("gamma")) continue;
      int[][] neighbors = neighborsOf(fixture.getValue());
      for (int[][] graph : new int[][][] {neighbors, sortedNeighbours(neighbors)})
        compareAtTheLimit(graph, fixture.getKey());
    }
    for (int units : new int[] {48, 60, 97}) compareAtTheLimit(neighborsOf(cyclodextrin(units)), "CD" + units);
    compareAtTheLimit(prismStack(80), "stack 80");
  }

  private static void compareAtTheLimit(int[][] neighbors, String label) {
    int n = neighbors.length;
    RingSearchNotes notes = new RingSearchNotes();
    RingPerception.Rings rings = RingPerception.perceive(n, neighbors, false, LIMIT, DEFAULTS, notes);
    assertSameRings(RingPerceptionReference.perceive(n, neighbors, false, LIMIT), rings, label);
    if (!notes.fallback) assertTrue(notes.baseWorkHigh <= LIMIT, label);
    else assertTrue(notes.baseWorkLow > LIMIT, label);
  }

  @Tag("SlowTest")
  @Test void benchmarkMoleculesKeepTheirRings() throws Exception {
    Path file = Paths.get("benchmarks/data/chembl_mcs_benchmark.smi");
    assertTrue(Files.exists(file), file.toAbsolutePath().toString());
    long graphs = 0, pairs = 0;
    long[] classes = new long[6];
    long seed = 5;
    try (BufferedReader reader = Files.newBufferedReader(file, StandardCharsets.UTF_8)) {
      for (String line; (line = reader.readLine()) != null; ) {
        String[] fields = line.trim().split("\\s+");
        if (fields.length == 0 || fields[0].isEmpty()) continue;
        int[][] neighbors;
        try {
          neighbors = neighborsOf(fields[0]);
        } catch (Exception unreadable) {
          continue;
        }
        if (RingPerception.cycleRank(neighbors.length, neighbors) <= 0) continue;
        Random random = new Random(seed++);
        for (int[][] graph : new int[][][] {neighbors, sortedNeighbours(neighbors),
                                            shuffledNeighbours(neighbors, random), relabelled(neighbors, random)}) {
          int n = graph.length;
          for (boolean enumerate : new boolean[] {false, true}) {
            RingSearchNotes notes = new RingSearchNotes();
            RingPerception.Rings rings = RingPerception.perceive(n, graph, enumerate, LIMIT, DEFAULTS, notes);
            assertSameRings(RingPerceptionReference.perceive(n, graph, enumerate, LIMIT), rings, fields.length > 1 ? fields[1] : fields[0]);
            if (enumerate) continue;
            graphs++;
            pairs += notes.pairs;
            classes[0] += notes.nonTreePairs;
            classes[1] += notes.classAPairs;
            classes[2] += notes.classB1Pairs;
            classes[3] += notes.classB1PlusPairs;
            classes[4] += notes.classB2Pairs;
            classes[5] += notes.simpleRingPairs;
          }
        }
      }
    }
    System.out.println("RING_SEARCH_BENCHMARK graphs=" + graphs + " pairs=" + pairs + " classes=" + Arrays.toString(classes));
    assertTrue(graphs > 4000);
  }
}
