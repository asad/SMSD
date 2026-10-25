/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
package com.bioinception.smsd.core;

import static com.bioinception.smsd.TestSupport.rawGraph;
import static org.junit.jupiter.api.Assertions.*;

import java.util.Arrays;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.CsvSource;

class MCSCompleteRingDomainTest {
  @ParameterizedTest
  @CsvSource({"MCES,Connected,false", "MCES,Connected,true", "MCES,Fragments,false",
      "MCES,Fragments,true", "MCIS,Connected,false", "MCIS,Connected,true",
      "MCIS,Fragments,false", "MCIS,Fragments,true"})
  void completeCycleMatchesWithoutRingFlags(
      MCS.Objective objective, MCS.Connectivity connectivity, boolean reverse) throws Exception {
    MolGraph built = cycle(6).graph(0);
    MolGraph parsed = rawGraph("C1CCCCC1");
    for (int atom = 0; atom < built.n; atom++) {
      assertFalse(built.ring[atom]);
      assertTrue(parsed.ring[atom]);
      for (int other : built.neighbors[atom]) assertFalse(built.bondInRing(atom, other));
    }
    assertFalse(built.ringsTruncated());
    assertFalse(parsed.ringsTruncated());
    MolGraph query = reverse ? parsed : built, target = reverse ? built : parsed;
    MCS.Options options = options(objective, connectivity);
    SearchResult result = MCS.find(query).withOptions(options).match(target);
    assertCompleted(result);
    assertEquals(6, result.score());
    assertTrue(result.found());
    SearchResult.Match match = result.matches().get(0);
    assertEquals(6, match.atomCount());
    assertEquals(6, match.bondCount());
    assertTrue(MCS.validateMapping(query, target, match, options).isEmpty());
  }

  @ParameterizedTest
  @CsvSource({"MCIS,Connected,false,false", "MCIS,Connected,false,true",
      "MCIS,Connected,true,false", "MCIS,Connected,true,true",
      "MCIS,Fragments,false,false", "MCIS,Fragments,false,true",
      "MCIS,Fragments,true,false", "MCIS,Fragments,true,true",
      "MCES,Connected,false,false", "MCES,Connected,false,true",
      "MCES,Connected,true,false", "MCES,Connected,true,true",
      "MCES,Fragments,false,false", "MCES,Fragments,false,true",
      "MCES,Fragments,true,false", "MCES,Fragments,true,true"})
  void explicitRingOnlyKeepsItsRawFlagContract(MCS.Objective objective,
      MCS.Connectivity connectivity, boolean reverse, boolean complete) throws Exception {
    MolGraph built = cycle(6).graph(0), parsed = rawGraph("C1CCCCC1");
    MolGraph query = reverse ? parsed : built, target = reverse ? built : parsed;
    MCS.Options options = options(objective, connectivity);
    options.chemistry.completeRingsOnly = complete;
    options.chemistry.ringMatchesRingOnly = true;
    SearchResult result = MCS.find(query).withOptions(options).match(target);
    assertCompleted(result);
    assertFalse(result.found());
    assertEquals(0, result.score());
  }

  @ParameterizedTest
  @CsvSource({"fused,Connected,false", "fused,Connected,true", "fused,Fragments,false",
      "fused,Fragments,true", "bridged,Connected,false", "bridged,Connected,true",
      "bridged,Fragments,false", "bridged,Fragments,true", "spiro,Connected,false",
      "spiro,Connected,true", "spiro,Fragments,false", "spiro,Fragments,true",
      "mixed,Connected,false", "mixed,Connected,true", "mixed,Fragments,false",
      "mixed,Fragments,true", "incomplete,Connected,false", "incomplete,Connected,true",
      "incomplete,Fragments,false", "incomplete,Fragments,true"})
  void everyMCISMappingAgreesWithSmallGraphOracle(
      String family, MCS.Connectivity connectivity, boolean reverse) {
    Fixture[] pair = pair(family);
    Fixture queryShape = pair[reverse ? 1 : 0], targetShape = pair[reverse ? 0 : 1];
    Oracle expected = oracle(queryShape, targetShape, connectivity);
    MCS.Options options = options(MCS.Objective.MCIS, connectivity);
    for (int flags = 0; flags < 3; flags++) {
      MolGraph query = queryShape.graph(flags == 0 ? 1 : 0);
      MolGraph target = targetShape.graph(flags == 0 ? 1 : flags == 1 ? 2 : 0);
      assertFalse(query.ringsTruncated());
      assertFalse(target.ringsTruncated());
      MCS prepared = MCS.find(query).withOptions(options);
      SearchResult all = prepared.matchAll(target, 1000);
      assertCompleted(all);
      assertEquals(expected.atoms, all.score());
      Set<String> actual = new TreeSet<>();
      for (SearchResult.Match match : all.matches()) {
        assertEquals(expected.atoms, match.atomCount());
        assertEquals(expected.bonds, match.bondCount());
        assertTrue(MCS.validateMapping(query, target, match, options).isEmpty());
        assertTrue(actual.add(key(match.atomMapping())));
      }
      assertEquals(expected.keys, actual, family + " flags=" + flags);
      SearchResult one = prepared.match(target);
      assertCompleted(one);
      assertEquals(expected.atoms, one.score());
      assertEquals(!expected.keys.isEmpty(), one.found());
      if (one.found()) {
        assertEquals(expected.bonds, one.matches().get(0).bondCount());
        assertTrue(expected.keys.contains(key(one.matches().get(0).atomMapping())));
        assertTrue(MCS.validateMapping(query, target, one.matches().get(0), options).isEmpty());
      }
    }
  }

  private static MCS.Options options(MCS.Objective objective, MCS.Connectivity connectivity) {
    MCS.Options options = new MCS.Options();
    options.objective = objective;
    options.connectivity = connectivity;
    options.chemistry.completeRingsOnly = true;
    options.timeoutMillis = 10_000;
    return options;
  }

  private static void assertCompleted(SearchResult result) {
    assertEquals(SearchResult.Status.Completed, result.status());
    assertTrue(result.optimal());
    assertTrue(result.complete());
  }

  private static Fixture cycle(int atoms) {
    int[][] edges = new int[atoms][2];
    int[] cycle = new int[atoms];
    for (int atom = 0; atom < atoms; atom++) {
      edges[atom] = new int[] {atom, (atom + 1) % atoms};
      cycle[atom] = atom;
    }
    return new Fixture(atoms, edges, new int[][] {cycle});
  }

  private static Fixture[] pair(String family) {
    switch (family) {
      case "fused":
        return new Fixture[] {cycle(3), new Fixture(4,
            new int[][] {{0, 1}, {1, 2}, {0, 2}, {1, 3}, {0, 3}},
            new int[][] {{0, 1, 2}, {0, 1, 3}})};
      case "bridged":
        return new Fixture[] {cycle(4), new Fixture(5,
            new int[][] {{0, 2}, {1, 2}, {0, 3}, {1, 3}, {0, 4}, {1, 4}},
            new int[][] {{0, 2, 1, 3}, {0, 2, 1, 4}, {0, 3, 1, 4}})};
      case "spiro":
        return new Fixture[] {cycle(3), new Fixture(5,
            new int[][] {{0, 1}, {1, 2}, {0, 2}, {0, 3}, {3, 4}, {0, 4}},
            new int[][] {{0, 1, 2}, {0, 3, 4}})};
      case "mixed":
        return new Fixture[] {new Fixture(4,
            new int[][] {{0, 1}, {1, 2}, {0, 2}, {2, 3}}, new int[][] {{0, 1, 2}}),
            new Fixture(5, new int[][] {{0, 1}, {1, 2}, {0, 2}, {2, 3}, {3, 4}},
                new int[][] {{0, 1, 2}})};
      case "incomplete":
        return new Fixture[] {cycle(3), new Fixture(3,
            new int[][] {{0, 1}, {1, 2}}, new int[0][])};
      default: throw new IllegalArgumentException(family);
    }
  }

  /** Small single-bond carbon graphs with their explicitly listed relevant cycles. */
  private static final class Fixture {
    final int atoms;
    final int[][] cycles;
    final boolean[][] edges, ringEdges;
    final boolean[] ringAtoms;

    Fixture(int atoms, int[][] bonds, int[][] cycles) {
      this.atoms = atoms;
      this.cycles = cycles;
      edges = new boolean[atoms][atoms];
      ringEdges = new boolean[atoms][atoms];
      ringAtoms = new boolean[atoms];
      for (int[] bond : bonds) edges[bond[0]][bond[1]] = edges[bond[1]][bond[0]] = true;
      for (int[] cycle : cycles)
        for (int index = 0; index < cycle.length; index++) {
          int first = cycle[index], second = cycle[(index + 1) % cycle.length];
          ringAtoms[first] = true;
          ringEdges[first][second] = ringEdges[second][first] = true;
        }
    }

    /** Flags: absent, all topological ring flags, or alternating topological ring flags. */
    MolGraph graph(int flags) {
      int[] elements = new int[atoms];
      Arrays.fill(elements, 6);
      int[][] neighbours = new int[atoms][], orders = new int[atoms][];
      boolean[] atomFlags = new boolean[atoms];
      boolean[][] bondFlags = new boolean[atoms][atoms];
      for (int atom = 0; atom < atoms; atom++) {
        int degree = 0;
        for (boolean bonded : edges[atom]) if (bonded) degree++;
        neighbours[atom] = new int[degree];
        orders[atom] = new int[degree];
        Arrays.fill(orders[atom], 1);
        int offset = 0;
        for (int other = 0; other < atoms; other++) if (edges[atom][other]) {
          neighbours[atom][offset++] = other;
          bondFlags[atom][other] = flags != 0 && ringEdges[atom][other]
              && (flags == 1 || ((atom + other) & 1) == 0);
        }
        atomFlags[atom] = flags != 0 && ringAtoms[atom] && (flags == 1 || (atom & 1) == 0);
      }
      return new MolGraph.Builder().atomCount(atoms).atomicNumbers(elements)
          .neighbors(neighbours).bondOrders(orders).ringFlags(atomFlags)
          .bondRingFlags(bondFlags).build();
    }

    boolean complete(boolean[] selected) {
      boolean[] coveredAtoms = new boolean[atoms];
      boolean[][] coveredEdges = new boolean[atoms][atoms];
      for (int[] cycle : cycles) {
        boolean whole = true;
        for (int atom : cycle) whole &= selected[atom];
        if (!whole) continue;
        for (int index = 0; index < cycle.length; index++) {
          int first = cycle[index], second = cycle[(index + 1) % cycle.length];
          coveredAtoms[first] = true;
          coveredEdges[first][second] = coveredEdges[second][first] = true;
        }
      }
      for (int first = 0; first < atoms; first++) if (selected[first]) {
        if (ringAtoms[first] && !coveredAtoms[first]) return false;
        for (int second = first + 1; second < atoms; second++)
          if (selected[second] && ringEdges[first][second] && !coveredEdges[first][second])
            return false;
      }
      return true;
    }
  }

  private static final class Oracle {
    int atoms, bonds;
    final Set<String> keys = new TreeSet<>();
  }

  private static Oracle oracle(Fixture query, Fixture target, MCS.Connectivity connectivity) {
    Oracle result = new Oracle();
    int[] mapping = new int[query.atoms];
    Arrays.fill(mapping, -1);
    enumerate(query, target, connectivity, 0, mapping, new boolean[target.atoms], result);
    return result;
  }

  private static void enumerate(Fixture query, Fixture target, MCS.Connectivity connectivity,
      int queryAtom, int[] mapping, boolean[] used, Oracle result) {
    if (queryAtom < query.atoms) {
      enumerate(query, target, connectivity, queryAtom + 1, mapping, used, result);
      for (int targetAtom = 0; targetAtom < target.atoms; targetAtom++) {
        if (used[targetAtom]) continue;
        boolean induced = true;
        for (int previous = 0; previous < queryAtom; previous++)
          if (mapping[previous] >= 0
              && query.edges[queryAtom][previous] != target.edges[targetAtom][mapping[previous]])
            induced = false;
        if (!induced) continue;
        mapping[queryAtom] = targetAtom;
        used[targetAtom] = true;
        enumerate(query, target, connectivity, queryAtom + 1, mapping, used, result);
        used[targetAtom] = false;
        mapping[queryAtom] = -1;
      }
      return;
    }
    boolean[] selected = new boolean[query.atoms];
    int atoms = 0, bonds = 0, firstSelected = -1;
    for (int atom = 0; atom < query.atoms; atom++) if (mapping[atom] >= 0) {
      selected[atom] = true;
      firstSelected = atom;
      atoms++;
      for (int other = 0; other < atom; other++)
        if (selected[other] && query.edges[atom][other]) bonds++;
    }
    if (atoms == 0 || !query.complete(selected) || !target.complete(used)) return;
    if (connectivity == MCS.Connectivity.Connected) {
      boolean[] reached = new boolean[query.atoms];
      reached[firstSelected] = true;
      for (int pass = 0; pass < query.atoms; pass++)
        for (int first = 0; first < query.atoms; first++) if (reached[first])
          for (int second = 0; second < query.atoms; second++)
            if (selected[second] && query.edges[first][second]) reached[second] = true;
      for (int atom = 0; atom < query.atoms; atom++) if (selected[atom] && !reached[atom]) return;
    }
    if (atoms < result.atoms || atoms == result.atoms && bonds < result.bonds) return;
    if (atoms > result.atoms || bonds > result.bonds) {
      result.atoms = atoms;
      result.bonds = bonds;
      result.keys.clear();
    }
    StringBuilder key = new StringBuilder();
    for (int atom = 0; atom < mapping.length; atom++)
      if (mapping[atom] >= 0) key.append(atom).append(':').append(mapping[atom]).append(',');
    result.keys.add(key.toString());
  }

  private static String key(Map<Integer, Integer> mapping) {
    StringBuilder key = new StringBuilder();
    for (int atom : new TreeSet<>(mapping.keySet()))
      key.append(atom).append(':').append(mapping.get(atom)).append(',');
    return key.toString();
  }
}
