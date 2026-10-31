/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
package com.bioinception.smsd.core;

import static com.bioinception.smsd.TestSupport.rawGraph;
import static com.bioinception.smsd.core.RingSearchFixtures.*;
import static org.junit.jupiter.api.Assertions.*;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;
import org.junit.jupiter.api.Test;

/** Relevant cycles grouped into families, against the relevant cycles enumerated one by one. */
class RelevantCycleFamiliesTest {
  private static MolGraph graphOf(int[][] neighbors) {
    int n = neighbors.length;
    int[] elements = new int[n];
    Arrays.fill(elements, 6);
    int[][] orders = new int[n][];
    for (int atom = 0; atom < n; atom++) {
      orders[atom] = new int[neighbors[atom].length];
      Arrays.fill(orders[atom], 1);
    }
    return new MolGraph.Builder().atomCount(n).atomicNumbers(elements).neighbors(copy(neighbors))
        .bondOrders(orders).build();
  }

  private static List<Long> key(int[] cycle) {
    List<Long> bonds = new ArrayList<>();
    for (int i = 0; i < cycle.length; i++) bonds.add(MolGraph.bondKey(cycle[i], cycle[(i + 1) % cycle.length]));
    bonds.sort(null);
    return bonds;
  }

  /** Every member of every family, checked to be a cycle of the graph that names its own family. */
  private static Set<List<Long>> members(MolGraph graph, RelevantCycleFamilies families) {
    Set<List<Long>> found = new HashSet<>();
    List<RelevantCycleFamily> list = families.families();
    for (int index = 0; index < list.size(); index++) {
      RelevantCycleFamily family = list.get(index);
      int[][] members = graph.relevantCycleFamilyMembers(family, Integer.MAX_VALUE);
      assertEquals(family.memberCount(), BigInteger.valueOf(members.length));
      assertArrayEquals(family.prototype(), members[0]);
      for (int[] member : members) {
        assertEquals(family.size(), member.length);
        assertEquals(family.root(), member[0]);
        assertEquals(index, graph.relevantCycleFamilyOf(families, member), Arrays.toString(member));
        assertTrue(found.add(key(member)), "members are distinct");
      }
    }
    return found;
  }

  /** Whether the relevant cycles were complete, so that the families could be compared with them. */
  private static boolean assertFamiliesAreTheRelevantCycles(MolGraph graph, String label) {
    int[][] relevant = graph.computeRelevantCycles();
    if (graph.relevantCyclesTruncated()) return false;
    RelevantCycleFamilies families = graph.computeRelevantCycleFamilies();
    assertFalse(families.truncated(), label);
    assertEquals(BigInteger.valueOf(relevant.length), families.relevantCycleCount(), label);
    Set<List<Long>> expected = new HashSet<>();
    for (int[] cycle : relevant) expected.add(key(cycle));
    assertEquals(expected, members(graph, families), label);
    for (int i = 1; i < families.families().size(); i++) {
      RelevantCycleFamily before = families.families().get(i - 1), after = families.families().get(i);
      assertTrue(before.size() < after.size() || (before.size() == after.size() && before.root() <= after.root()), label);
    }
    return true;
  }

  @Test void theFamiliesHoldExactlyTheRelevantCycles() throws Exception {
    int fixtures = 0, compared = 0;
    for (Map.Entry<String, String> fixture : smiles().entrySet())
      if (small(fixture.getKey()) && assertFamiliesAreTheRelevantCycles(graphOf(neighborsOf(fixture.getValue())), fixture.getKey()))
        fixtures++;
    assertTrue(fixtures >= 20, "fixtures compared: " + fixtures);
    for (int units = 6; units <= 12; units++)
      assertFamiliesAreTheRelevantCycles(rawGraph(cyclodextrin(units)), "CD" + units);
    for (String text : new String[] {"C12C3C4C1C5C2C3C45", "C1C2CC3CC1CC(C2)C3", "C1CC2CCC1C2", "c1ccc2ccccc2c1",
                                     "C1CC1", "C1CCCCC1", "C1CCC2(CC1)CCCC2"})
      assertFamiliesAreTheRelevantCycles(rawGraph(text), text);
    for (int layers = 2; layers <= 8; layers++) assertFamiliesAreTheRelevantCycles(graphOf(prismStack(layers)), "stack " + layers);
    assertFamiliesAreTheRelevantCycles(graphOf(cyclobutaneChain(5)), "cyclobutane 5");
    assertFamiliesAreTheRelevantCycles(graphOf(doubleHub(12)), "double hub 12");
    Random random = new Random(20261030L);
    for (int i = 0; i < 150; i++) {
      int[][] neighbors = randomRingGraph(random);
      for (int[][] graph : variants(neighbors, i))
        if (assertFamiliesAreTheRelevantCycles(graphOf(graph), "random " + i)) compared++;
    }
    assertTrue(compared >= 1000, "random graphs compared: " + compared);
  }

  @Test void cyclodextrinTotalsAreExact() throws Exception {
    assertEquals("17179869218", rawGraph(cyclodextrin(34)).computeRelevantCycleFamilies().relevantCycleCount().toString());
    assertEquals("1099511627816", rawGraph(cyclodextrin(40)).computeRelevantCycleFamilies().relevantCycleCount().toString());
    RelevantCycleFamilies cd97 = rawGraph(cyclodextrin(97)).computeRelevantCycleFamilies();
    assertFalse(cd97.truncated());
    assertEquals("158456325028528675187087900769", cd97.relevantCycleCount().toString());
    assertEquals(BigInteger.ONE.shiftLeft(97).add(BigInteger.valueOf(97)), cd97.relevantCycleCount());
    // The 97 rings of the units, each a family of one; the macrocycles share the rest.
    BigInteger macrocycles = BigInteger.ZERO;
    int units = 0;
    for (RelevantCycleFamily family : cd97.families()) {
      if (family.size() == 6) { units++; assertEquals(BigInteger.ONE, family.memberCount()); }
      else macrocycles = macrocycles.add(family.memberCount());
    }
    assertEquals(97, units);
    assertEquals(BigInteger.ONE.shiftLeft(97), macrocycles);
    assertEquals("442", graphOf(prismStack(64)).computeRelevantCycleFamilies().relevantCycleCount().toString());
    assertEquals("554", graphOf(prismStack(80)).computeRelevantCycleFamilies().relevantCycleCount().toString());
  }

  @Test void cubaneHasSixFamiliesOfOneMember() throws Exception {
    MolGraph cubane = rawGraph("C12C3C4C1C5C2C3C45");
    RelevantCycleFamilies families = cubane.computeRelevantCycleFamilies();
    assertEquals(6, families.families().size());
    for (RelevantCycleFamily family : families.families()) {
      assertEquals(4, family.size());
      assertEquals(BigInteger.ONE, family.memberCount());
    }
    assertEquals(BigInteger.valueOf(6), families.relevantCycleCount());
  }

  @Test void aSimpleRingNamesItsHighestAtomAndTheAtomsOppositeIt() {
    MolGraph six = graphOf(graph(6, new int[][] {{0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 5}, {5, 0}}));
    RelevantCycleFamily even = six.computeRelevantCycleFamilies().families().get(0);
    assertEquals(6, even.size());
    assertEquals(5, even.root());
    assertEquals(2, even.opposite());
    assertEquals(1, even.first());
    assertEquals(3, even.second());
    assertArrayEquals(new int[] {5, 0, 1, 2, 3, 4}, even.prototype());
    assertEquals(BigInteger.ONE, even.memberCount());
    MolGraph five = graphOf(graph(5, new int[][] {{4, 3}, {3, 2}, {2, 1}, {1, 0}, {0, 4}}));
    RelevantCycleFamily odd = five.computeRelevantCycleFamilies().families().get(0);
    assertEquals(5, odd.size());
    assertEquals(4, odd.root());
    assertEquals(-1, odd.opposite());
    assertEquals(1, odd.first());
    assertEquals(2, odd.second());
    assertArrayEquals(new int[] {4, 0, 1, 2, 3}, odd.prototype());
    assertEquals(0, five.relevantCycleFamilyOf(five.computeRelevantCycleFamilies(), new int[] {2, 3, 4, 0, 1}));
    assertEquals(0, five.relevantCycleFamilyOf(five.computeRelevantCycleFamilies(), new int[] {1, 0, 4, 3, 2}));
  }

  @Test void membersComeInOneFixedOrderStartingWithThePrototype() throws Exception {
    MolGraph cd8 = rawGraph(cyclodextrin(8));
    RelevantCycleFamilies families = cd8.computeRelevantCycleFamilies();
    RelevantCycleFamily macrocycle = families.families().get(families.families().size() - 1);
    assertEquals(BigInteger.valueOf(256), macrocycle.memberCount());
    int[][] all = cd8.relevantCycleFamilyMembers(macrocycle, 1000);
    assertEquals(256, all.length);
    assertArrayEquals(macrocycle.prototype(), all[0]);
    assertArrayEquals(all, cd8.relevantCycleFamilyMembers(macrocycle, 1000));
    assertArrayEquals(Arrays.copyOf(all, 5), cd8.relevantCycleFamilyMembers(macrocycle, 5));
    assertEquals(0, cd8.relevantCycleFamilyMembers(macrocycle, 0).length);
    Set<List<Long>> distinct = new HashSet<>();
    for (int[] member : all) assertTrue(distinct.add(key(member)));
  }

  @Test void aLowerLimitMarksTheFamiliesTruncated() throws Exception {
    MolGraph cubane = rawGraph("C12C3C4C1C5C2C3C45");
    long[] work = new long[1];
    RelevantCycleFamilies complete = RelevantCycleFamilySearch.compute(cubane.n, cubane.neighbors,
                                                                       RelevantCycleFamilySearch.FAMILY_WORK_LIMIT, work);
    assertFalse(complete.truncated());
    RelevantCycleFamilies cut = RelevantCycleFamilySearch.compute(cubane.n, cubane.neighbors, work[0] - 1, null);
    assertTrue(cut.truncated());
    assertTrue(cut.families().isEmpty());
    assertFalse(RelevantCycleFamilySearch.compute(cubane.n, cubane.neighbors, work[0], null).truncated());
    // A ring system that runs out leaves the systems before it.
    int[][] two = graph(9, new int[][] {{0, 1}, {1, 2}, {2, 0}, {2, 3}, {3, 4}, {4, 5}, {5, 6}, {6, 7}, {7, 8}, {8, 3}, {4, 7}});
    long[] twoWork = new long[1];
    RelevantCycleFamilies twoComplete = RelevantCycleFamilySearch.compute(9, two, RelevantCycleFamilySearch.FAMILY_WORK_LIMIT, twoWork);
    assertEquals(3, twoComplete.families().size());
    RelevantCycleFamilies twoCut = RelevantCycleFamilySearch.compute(9, two, twoWork[0] - 1, null);
    assertTrue(twoCut.truncated());
    assertEquals(1, twoCut.families().size());
    assertEquals(3, twoCut.families().get(0).size());
    assertEquals(BigInteger.ONE, twoCut.relevantCycleCount());
  }

  @Test void halvesThatAreNotShortestPathsNameNoFamily() {
    MolGraph graph = graphOf(eightRingGraph());
    RelevantCycleFamilies families = graph.computeRelevantCycleFamilies();
    assertFalse(families.truncated());
    assertEquals(-1, graph.relevantCycleFamilyOf(families, eightRing()));
    assertEquals(-1, graph.relevantCycleFamilyOf(families, new int[] {7, 0, 3}));
    assertEquals(-1, graph.relevantCycleFamilyOf(families, new int[] {7, 0}));
    assertEquals(-1, graph.relevantCycleFamilyOf(families, null));
    assertEquals(-1, graph.relevantCycleFamilyOf(families, new int[] {7, 0, 3, 1, 3, 0}));
    for (int[] cycle : graph.computeRelevantCycles())
      assertTrue(graph.relevantCycleFamilyOf(families, cycle) >= 0, Arrays.toString(cycle));
  }

  @Test void theFamiliesAreKeptAndCannotBeChanged() throws Exception {
    MolGraph graph = rawGraph("C1C2CC3CC1CC(C2)C3");
    RelevantCycleFamilies families = graph.computeRelevantCycleFamilies();
    assertSame(families, graph.computeRelevantCycleFamilies());
    assertThrows(UnsupportedOperationException.class, () -> families.families().clear());
    RelevantCycleFamily family = families.families().get(0);
    int[] prototype = family.prototype();
    prototype[0] = -5;
    assertNotEquals(-5, family.prototype()[0]);
    assertEquals(0, new MolGraph.Builder().atomCount(0).atomicNumbers(new int[0]).neighbors(new int[0][])
                        .bondOrders(new int[0][]).build().computeRelevantCycleFamilies().families().size());
  }
}
