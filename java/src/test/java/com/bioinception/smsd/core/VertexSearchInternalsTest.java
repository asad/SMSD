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
import org.junit.jupiter.api.Assumptions;
import org.junit.jupiter.api.Test;

/** The complete-ring test and the work arrays of the general MCIS search. */
class VertexSearchInternalsTest {
  /** Fused, bridged, spiro, cage and macrocyclic ring systems, with side chains and a chain. */
  private static final String[] RING_SYSTEMS = {
      "c1ccc2ccccc2c1", "C1CCC2CCCCC2C1", "c1ccc2c(c1)ccc1ccccc12", "C1CCC2C(C1)CCC1C2CCC2CCCC21",
      "C1CC2CCC3CCC4CCC1C1C2C3C41", "C1CC2CC3CCC2CC13", "C1CC2CCC1C2", "C1C2CC3CC1CC(C2)C3",
      "C1CC2(C1)CC2", "C1CCC2(CC1)CCCC2", "C1CCC2(CC1)CCC1(CC2)CCCC1", "C12C3C4C1C5C2C3C45",
      "C12C3C4C1C5C3C4C25", "C1CCCCCCCCCCC1", "C1COCCOCCOCCOCCOCCO1",
      "CC(C)CC1NC(=O)C(CC(=O)O)NC(=O)C(CC(C)C)NC(=O)C(Cc2ccccc2)NC(=O)C(CCCCN)NC(=O)C(CO)NC1=O",
      "OC(=O)CCc1ccc(cc1)C1CC1", "CCCCCCCC"};

  private static SearchResult.Match match(Map<Integer, Integer> mapping) {
    return new SearchResult.Match(
        mapping, Collections.emptyMap(), Collections.emptyList(), mapping.size());
  }

  /** The complete-ring rule as {@link MCS#validateMapping} applies it to both graphs. */
  private static boolean validRings(MolGraph query, MolGraph target, Map<Integer, Integer> mapping) {
    Set<SearchResult.Bond> queryBonds = new HashSet<>(), targetBonds = new HashSet<>();
    for (int atom : mapping.keySet())
      for (int neighbor : query.neighbors[atom])
        if (atom < neighbor && mapping.containsKey(neighbor)) {
          queryBonds.add(new SearchResult.Bond(atom, neighbor));
          targetBonds.add(new SearchResult.Bond(mapping.get(atom), mapping.get(neighbor)));
        }
    return VertexSearch.completeRings(query, mapping.keySet(), queryBonds)
        && VertexSearch.completeRings(target, new HashSet<>(mapping.values()), targetBonds);
  }

  /** A copy of {@code graph} in which each atom is relabelled {@code order[atom]}. */
  private static MolGraph relabelled(MolGraph graph, int[] order) {
    int[] elements = new int[graph.n];
    boolean[] ring = new boolean[graph.n];
    int[][] neighbors = new int[graph.n][], bondOrders = new int[graph.n][];
    boolean[][] bondRings = new boolean[graph.n][];
    for (int atom = 0; atom < graph.n; atom++) {
      int label = order[atom], degree = graph.neighbors[atom].length;
      elements[label] = graph.atomicNum[atom];
      ring[label] = graph.ring[atom];
      neighbors[label] = new int[degree];
      bondOrders[label] = new int[degree];
      bondRings[label] = new boolean[degree];
      for (int slot = 0; slot < degree; slot++) {
        int neighbor = graph.neighbors[atom][slot];
        neighbors[label][slot] = order[neighbor];
        bondOrders[label][slot] = graph.bondOrder(atom, neighbor);
        bondRings[label][slot] = graph.bondInRing(atom, neighbor);
      }
    }
    return new MolGraph.Builder().atomCount(graph.n).atomicNumbers(elements).ringFlags(ring)
        .neighbors(neighbors).bondOrders(bondOrders).bondRingFlags(bondRings).build();
  }

  /**
   * Query atoms chosen in one of three ways: whole relevant cycles, with most atoms off cycles
   * and a few on them; the same less one atom; or each atom at random.
   */
  private static boolean[] selection(MolGraph query, Random random, int kind) {
    boolean[] selected = new boolean[query.n];
    if (kind == 2) {
      int percent = 25 + random.nextInt(70);
      for (int atom = 0; atom < query.n; atom++)
        selected[atom] = random.nextInt(100) < percent;
      return selected;
    }
    boolean[] onCycle = query.ring.clone();
    for (int[] cycle : query.computeRelevantCycles()) {
      boolean whole = random.nextBoolean();
      for (int atom : cycle) {
        onCycle[atom] = true;
        selected[atom] |= whole;
      }
    }
    List<Integer> chosen = new ArrayList<>();
    for (int atom = 0; atom < query.n; atom++) {
      if (!selected[atom] && random.nextInt(onCycle[atom] ? 16 : 2) == 0)
        selected[atom] = true;
      if (selected[atom])
        chosen.add(atom);
    }
    if (kind == 1 && !chosen.isEmpty())
      selected[chosen.get(random.nextInt(chosen.size()))] = false;
    return selected;
  }

  /**
   * Images of the selected atoms: {@code order[atom]} when an order is given; otherwise atoms
   * taken in breadth-first order, each mapped to an unused neighbour of the image of a mapped
   * neighbour where there is one and {@code alongBonds}, and to any unused target atom otherwise.
   */
  private static Map<Integer, Integer> mapping(MolGraph query, MolGraph target, boolean[] selected,
      Random random, int[] order, boolean alongBonds) {
    Map<Integer, Integer> mapping = new TreeMap<>();
    if (order != null) {
      for (int atom = 0; atom < query.n; atom++)
        if (selected[atom])
          mapping.put(atom, order[atom]);
      return mapping;
    }
    boolean[] used = new boolean[target.n], queued = new boolean[query.n];
    int[] queue = new int[query.n];
    int offset = random.nextInt(query.n);
    for (int index = 0; index < query.n; index++) {
      int root = (index + offset) % query.n;
      if (!selected[root] || queued[root])
        continue;
      int head = 0, tail = 0;
      queue[tail++] = root;
      queued[root] = true;
      while (head < tail) {
        int atom = queue[head++];
        List<Integer> candidates = new ArrayList<>();
        if (alongBonds)
          for (int neighbor : query.neighbors[atom]) {
            Integer image = mapping.get(neighbor);
            if (image != null)
              for (int adjacent : target.neighbors[image])
                if (!used[adjacent])
                  candidates.add(adjacent);
          }
        if (candidates.isEmpty())
          for (int other = 0; other < target.n; other++)
            if (!used[other])
              candidates.add(other);
        if (!candidates.isEmpty()) {
          int image = candidates.get(random.nextInt(candidates.size()));
          used[image] = true;
          mapping.put(atom, image);
        }
        for (int neighbor : query.neighbors[atom])
          if (selected[neighbor] && !queued[neighbor]) {
            queued[neighbor] = true;
            queue[tail++] = neighbor;
          }
      }
    }
    return mapping;
  }

  /**
   * Checks 36 mappings on one search against the rule, three in four of them by {@code order}
   * where it is given. The stamps start near their largest value, so the search also passes
   * their reset. Counts the accepted mappings that hold a ring atom, and the rejected ones.
   */
  private static void check(MolGraph query, MolGraph target, int[] order, Random random, int[] counts) {
    MCS.Options options = new MCS.Options();
    options.chemistry.completeRingsOnly = true;
    VertexSearch search = new VertexSearch(query, target, options, 1, false);
    search.ringStamp = Integer.MAX_VALUE - 16;
    for (int trial = 0; trial < 36; trial++) {
      boolean[] selected = selection(query, random, trial % 3);
      Map<Integer, Integer> mapping = mapping(query, target, selected, random,
          trial % 4 != 3 ? order : null, trial % 5 != 4);
      search.seed(match(mapping));
      boolean expected = validRings(query, target, mapping);
      assertEquals(expected, search.ringsComplete(), query.n + " to " + target.n + " " + mapping);
      boolean onRing = false;
      for (Map.Entry<Integer, Integer> pair : mapping.entrySet())
        onRing |= query.ring[pair.getKey()] || target.ring[pair.getValue()];
      if (!expected)
        counts[1]++;
      else if (onRing)
        counts[0]++;
    }
    assertEquals(20, search.ringStamp);
  }

  @Test
  void ringTestAgreesWithTheValidationRule() throws Exception {
    // Random mappings between every pair of ring systems, and from each system to a copy with its
    // atoms relabelled, against the rule that validateMapping applies.
    List<MolGraph> graphs = new ArrayList<>();
    for (String smiles : RING_SYSTEMS)
      graphs.add(rawGraph(smiles));
    Random random = new Random(803);
    int[] counts = new int[2];
    for (MolGraph query : graphs) {
      int[] identity = new int[query.n];
      List<Integer> labels = new ArrayList<>();
      for (int atom = 0; atom < query.n; atom++) {
        identity[atom] = atom;
        labels.add(atom);
      }
      for (MolGraph target : graphs)
        check(query, target, query == target ? identity : null, random, counts);
      Collections.shuffle(labels, random);
      int[] order = new int[query.n];
      for (int atom = 0; atom < query.n; atom++)
        order[atom] = labels.get(atom);
      check(query, relabelled(query, order), order, random, counts);
    }
    assertTrue(counts[0] >= 400 && counts[1] >= 5000,
        counts[0] + " accepted, " + counts[1] + " rejected");
  }

  @Test
  void ringStampsClearEveryMarkWhenTheyWrap() throws Exception {
    // The first test marks both rings of naphthalene with stamp 1, and the last test before the
    // wrap marks only the second ring. After the wrap the stamp is 1 again, so a mark left on the
    // first ring would accept it once it is no longer whole.
    MolGraph naphthalene = rawGraph("c1ccc2ccccc2c1");
    MCS.Options options = new MCS.Options();
    options.chemistry.completeRingsOnly = true;
    VertexSearch search = new VertexSearch(naphthalene, naphthalene, options, 1, false);
    Map<Integer, Integer> whole = new TreeMap<>(), second = new TreeMap<>();
    for (int atom = 0; atom < naphthalene.n; atom++) {
      whole.put(atom, atom);
      if (atom >= 3 && atom <= 8)
        second.put(atom, atom);
    }
    search.seed(match(whole));
    assertTrue(search.ringsComplete());
    assertEquals(1, search.ringStamp);
    search.ringStamp = Integer.MAX_VALUE - 1;
    search.seed(match(second));
    assertTrue(search.ringsComplete());
    assertEquals(Integer.MAX_VALUE, search.ringStamp);
    Map<Integer, Integer> broken = new TreeMap<>(whole);
    broken.remove(0);
    search.seed(match(broken));
    assertFalse(search.ringsComplete());
    assertEquals(1, search.ringStamp);
  }

  /** The allocation counter of each thread; the test is skipped where the virtual machine has none. */
  private static com.sun.management.ThreadMXBean allocations() {
    java.lang.management.ThreadMXBean threads = java.lang.management.ManagementFactory.getThreadMXBean();
    Assumptions.assumeTrue(threads instanceof com.sun.management.ThreadMXBean);
    com.sun.management.ThreadMXBean allocations = (com.sun.management.ThreadMXBean) threads;
    Assumptions.assumeTrue(allocations.isThreadAllocatedMemorySupported()
        && allocations.isThreadAllocatedMemoryEnabled());
    return allocations;
  }

  /** Options whose cancellation stops the search after {@code polls} polls, with no clock. */
  private static MCS.Options polled(MCS.Connectivity connectivity, long polls) {
    MCS.Options options = new MCS.Options();
    options.connectivity = connectivity;
    options.timeoutMillis = 3_600_000L;
    long[] count = {0};
    options.cancellation = () -> ++count[0] > polls;
    return options;
  }

  @Test
  void fragmentSearchAllocatesWorkOnlyForTheDepthsItReaches() throws Exception {
    // Two chains of 4000 carbons in the native fragment order, stopped soon after the domains are
    // built. Each domain takes 500 bytes of bits; an int for every pair of query and target
    // atoms, as domain scratch or as candidates for every depth, would take 64 MB.
    com.sun.management.ThreadMXBean allocations = allocations();
    StringBuilder smiles = new StringBuilder();
    for (int atom = 0; atom < 4000; atom++)
      smiles.append('C');
    MolGraph chain = rawGraph(smiles.toString());
    MCS.Options options = polled(MCS.Connectivity.Fragments, 4040);
    long thread = Thread.currentThread().getId();
    long start = allocations.getThreadAllocatedBytes(thread);
    VertexSearch search = new VertexSearch(chain, chain, options, 1, false);
    long built = allocations.getThreadAllocatedBytes(thread);
    SearchResult result = search.search();
    long searched = allocations.getThreadAllocatedBytes(thread);
    assertEquals(SearchResult.Status.Cancelled, result.status());
    assertEquals(38L, result.nodesVisited());
    assertEquals(37, result.matches().get(0).atomCount());
    assertTrue(built - start < 16L << 20, "construction allocated " + ((built - start) >> 20) + " MiB");
    assertTrue(searched - built < 16L << 20, "search allocated " + ((searched - built) >> 20) + " MiB");
  }

  @Test
  void completeRingTestAllocatesNothingAtEachNode() throws Exception {
    // The connected strychnine/quinine search runs under a poll budget. Domain pruning can
    // finish it before cancellation, so the stamped ring helper also runs 10,000 times on
    // its own; fewer search nodes cannot hide an allocation at every ring test.
    com.sun.management.ThreadMXBean allocations = allocations();
    MolGraph query = rawGraph("C1CN2CC3=CCOC4CC(=O)N5C6C4C3CC2C61C7=CC=CC=C75");
    MolGraph target = rawGraph("COC1=CC2=C(C=CN=C2C=C1)C(C3CC4CCN3CC4C=C)O");
    // Ring perception is cached on the graphs and is not part of the search.
    query.computeRelevantCycles();
    target.computeRelevantCycles();
    MCS.Options options = polled(MCS.Connectivity.Connected, 10_000);
    options.chemistry.completeRingsOnly = true;
    VertexSearch search = new VertexSearch(query, target, options, 1, false);
    long thread = Thread.currentThread().getId();
    long start = allocations.getThreadAllocatedBytes(thread);
    SearchResult result = search.search();
    long allocated = allocations.getThreadAllocatedBytes(thread) - start;
    boolean completed = result.status() == SearchResult.Status.Completed;
    if (!completed)
      assertEquals(SearchResult.Status.Cancelled, result.status());
    assertEquals(completed, result.optimal());
    assertEquals(completed, result.complete());
    assertEquals(6, result.matches().get(0).atomCount());
    assertEquals(Collections.emptyList(),
        MCS.validateMapping(query, target, result.matches().get(0), options));
    assertTrue(allocated < 4L << 20, "search allocated " + (allocated >> 20) + " MiB");

    search.seed(result.matches().get(0));
    assertTrue(search.ringsComplete());
    long ringStart = allocations.getThreadAllocatedBytes(thread);
    int complete = 0;
    for (int check = 0; check < 10_000; check++)
      if (search.ringsComplete())
        complete++;
    long ringAllocated = allocations.getThreadAllocatedBytes(thread) - ringStart;
    assertEquals(10_000, complete);
    assertTrue(ringAllocated < 4L << 20,
        "ring tests allocated " + (ringAllocated >> 20) + " MiB");
  }
}
