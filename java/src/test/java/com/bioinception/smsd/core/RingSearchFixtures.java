/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
package com.bioinception.smsd.core;

import static com.bioinception.smsd.TestSupport.rawGraph;
import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;

import java.io.BufferedReader;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;

/** Graphs and comparisons shared by the ring search tests. */
final class RingSearchFixtures {
  private RingSearchFixtures() {}

  /** The ring fixtures by name, as listed in ring_search_fixtures.smi. */
  static Map<String, String> smiles() throws Exception {
    Map<String, String> fixtures = new LinkedHashMap<>();
    try (InputStream in = RingSearchFixtures.class.getClassLoader()
             .getResourceAsStream("com/bioinception/smsd/ring_search_fixtures.smi")) {
      assertNotNull(in, "ring_search_fixtures.smi");
      BufferedReader reader = new BufferedReader(new InputStreamReader(in, StandardCharsets.UTF_8));
      for (String line; (line = reader.readLine()) != null; ) {
        String[] fields = line.trim().split("\\s+");
        if (fields.length >= 2) fixtures.put(fields[1], fields[0]);
      }
    }
    return fixtures;
  }

  /** Fixtures whose pairwise search stays small enough for every test run. */
  static boolean small(String name) {
    return !name.matches("(?i).*(cd3[4-9]|cd[4-9][0-9]|stack_(32|64|80|170)|c294|gamma).*");
  }

  /**
   * The neighbour lists of a graph read from SMILES: those of its MolGraph, or, when CDK cannot
   * assign a Kekulé structure, those of the parsed atoms and bonds in container order.
   */
  static int[][] neighborsOf(String smiles) throws Exception {
    try {
      return copy(rawGraph(smiles).neighbors);
    } catch (InvalidSmilesException unassigned) {
      SmilesParser parser = new SmilesParser(SilentChemObjectBuilder.getInstance());
      parser.kekulise(false);
      IAtomContainer molecule = parser.parseSmiles(smiles);
      int[][] neighbors = new int[molecule.getAtomCount()][];
      for (int atom = 0; atom < neighbors.length; atom++) {
        List<IBond> bonds = molecule.getConnectedBondsList(molecule.getAtom(atom));
        neighbors[atom] = new int[bonds.size()];
        for (int k = 0; k < bonds.size(); k++)
          neighbors[atom][k] = molecule.indexOf(bonds.get(k).getOther(molecule.getAtom(atom)));
      }
      return neighbors;
    }
  }

  static int[][] copy(int[][] neighbors) {
    int[][] copy = new int[neighbors.length][];
    for (int i = 0; i < neighbors.length; i++) copy[i] = neighbors[i].clone();
    return copy;
  }

  /** A graph from its bonds, each atom listing its neighbours in bond order. */
  static int[][] graph(int n, int[][] bonds) {
    List<List<Integer>> lists = new ArrayList<>();
    for (int i = 0; i < n; i++) lists.add(new ArrayList<>());
    for (int[] bond : bonds) {
      lists.get(bond[0]).add(bond[1]);
      lists.get(bond[1]).add(bond[0]);
    }
    return arrays(lists);
  }

  static int[][] arrays(List<List<Integer>> lists) {
    int[][] neighbors = new int[lists.size()][];
    for (int i = 0; i < neighbors.length; i++) {
      neighbors[i] = new int[lists.get(i).size()];
      for (int k = 0; k < neighbors[i].length; k++) neighbors[i][k] = lists.get(i).get(k);
    }
    return neighbors;
  }

  /** A stack of hexagonal prisms: layers of six-membered rings joined atom to atom. */
  static int[][] prismStack(int layers) {
    List<int[]> bonds = new ArrayList<>();
    for (int layer = 0; layer < layers; layer++)
      for (int j = 0; j < 6; j++) {
        int a = 6 * layer + j;
        bonds.add(new int[] {a, 6 * layer + (j + 1) % 6});
        if (layer + 1 < layers) bonds.add(new int[] {a, a + 6});
      }
    return graph(6 * layers, bonds.toArray(new int[0][]));
  }

  /** The cyclobutane macrocycle of RelevantCycleCompletenessTest, neighbours in index order. */
  static int[][] cyclobutaneChain(int units) {
    int atoms = 5 * units + 1;
    boolean[][] bonded = new boolean[atoms][atoms];
    List<int[]> bonds = new ArrayList<>();
    for (int unit = 0; unit < units; unit++) {
      int corner = 3 * unit, opposite = 3 * (unit + 1);
      bonds.add(new int[] {corner, corner + 1});
      bonds.add(new int[] {corner + 1, opposite});
      bonds.add(new int[] {corner, corner + 2});
      bonds.add(new int[] {corner + 2, opposite});
    }
    int previous = 3 * units;
    for (int atom = 3 * units + 1; atom < atoms; atom++) {
      bonds.add(new int[] {previous, atom});
      previous = atom;
    }
    bonds.add(new int[] {previous, 0});
    for (int[] bond : bonds) bonded[bond[0]][bond[1]] = bonded[bond[1]][bond[0]] = true;
    int[][] neighbors = new int[atoms][];
    for (int atom = 0; atom < atoms; atom++) {
      int degree = 0;
      for (boolean b : bonded[atom]) if (b) degree++;
      neighbors[atom] = new int[degree];
      for (int other = 0, k = 0; other < atoms; other++) if (bonded[atom][other]) neighbors[atom][k++] = other;
    }
    return neighbors;
  }

  static String cyclodextrin(int units) {
    StringBuilder smiles = new StringBuilder("C%101C(O)C(O)C(OC1CO)O");
    for (int unit = 1; unit < units; unit++) smiles.append("C1C(O)C(O)C(OC1CO)O");
    return smiles.append("%10").toString();
  }

  /** Six atoms whose root-0 tree bypasses bond (0, 2) through a half-bond from one level deeper. */
  static int[][] extendedClassGraph() {
    return new int[][] {{1, 2, 3}, {0, 3, 4}, {0, 4, 5}, {0, 1}, {1, 2, 5}, {2, 4}};
  }

  /**
   * The eight-ring r-a1-b1-a2-b2-y-b3-a3 whose two halves from r are not shortest paths, with
   * a1, a2, a3, b1, b2, b3, y, r numbered 0 to 7.
   */
  static int[][] eightRingGraph() {
    int a1 = 0, a2 = 1, a3 = 2, b1 = 3, b2 = 4, b3 = 5, y = 6, r = 7;
    return graph(8, new int[][] {{r, a1}, {r, a2}, {r, a3}, {a1, b1}, {a2, b1}, {a2, b2}, {a3, b3},
                                 {b1, y}, {b2, y}, {b3, y}});
  }

  static int[] eightRing() {
    return new int[] {7, 0, 3, 1, 4, 6, 5, 2};
  }

  /** Two hubs joined through count middle atoms, every pair of middle atoms closing a four-ring. */
  static int[][] doubleHub(int count) {
    List<int[]> bonds = new ArrayList<>();
    for (int middle = 2; middle < count + 2; middle++) {
      bonds.add(new int[] {0, middle});
      bonds.add(new int[] {1, middle});
    }
    return graph(count + 2, bonds.toArray(new int[0][]));
  }

  /**
   * A seeded graph of one to three ring systems, each a cycle with chords, joined by short chains
   * and carrying pendant atoms, with degrees of at most four and neighbours in random order.
   */
  static int[][] randomRingGraph(Random random) {
    List<List<Integer>> lists = new ArrayList<>();
    int systems = 1 + random.nextInt(3), anchor = -1;
    for (int s = 0; s < systems; s++) {
      int size = 3 + random.nextInt(18), base = lists.size();
      for (int i = 0; i < size; i++) lists.add(new ArrayList<>());
      for (int i = 0; i < size; i++) bond(lists, base + i, base + (i + 1) % size);
      int chords = random.nextInt(size);
      for (int c = 0; c < chords; c++) bond(lists, base + random.nextInt(size), base + random.nextInt(size));
      if (anchor >= 0) {
        int from = anchor;
        for (int i = random.nextInt(3); i > 0; i--) {
          int atom = lists.size();
          lists.add(new ArrayList<>());
          bond(lists, from, atom);
          from = atom;
        }
        bond(lists, from, base + random.nextInt(size));
      }
      anchor = base + random.nextInt(size);
    }
    for (int p = random.nextInt(4); p > 0; p--) {
      int atom = lists.size();
      lists.add(new ArrayList<>());
      bond(lists, random.nextInt(atom), atom);
    }
    for (List<Integer> list : lists) Collections.shuffle(list, random);
    return arrays(lists);
  }

  private static void bond(List<List<Integer>> lists, int a, int b) {
    if (a == b || lists.get(a).contains(b) || lists.get(a).size() >= 4 || lists.get(b).size() >= 4) return;
    lists.get(a).add(b);
    lists.get(b).add(a);
  }

  static int[][] sortedNeighbours(int[][] neighbors) {
    int[][] sorted = copy(neighbors);
    for (int[] list : sorted) Arrays.sort(list);
    return sorted;
  }

  static int[][] shuffledNeighbours(int[][] neighbors, Random random) {
    int[][] shuffled = copy(neighbors);
    for (int[] list : shuffled)
      for (int i = list.length - 1; i > 0; i--) {
        int j = random.nextInt(i + 1), swap = list[i];
        list[i] = list[j];
        list[j] = swap;
      }
    return shuffled;
  }

  /** The graph with its atoms renumbered at random, each list keeping its order. */
  static int[][] relabelled(int[][] neighbors, Random random) {
    int n = neighbors.length;
    List<Integer> order = new ArrayList<>();
    for (int i = 0; i < n; i++) order.add(i);
    Collections.shuffle(order, random);
    int[] label = new int[n];
    for (int i = 0; i < n; i++) label[i] = order.get(i);
    int[][] relabelled = new int[n][];
    for (int atom = 0; atom < n; atom++) {
      int[] list = new int[neighbors[atom].length];
      for (int k = 0; k < list.length; k++) list[k] = label[neighbors[atom][k]];
      relabelled[label[atom]] = list;
    }
    return relabelled;
  }

  /** The variants of a graph the differential runs: as given, sorted, three shuffles and three relabellings. */
  static List<int[][]> variants(int[][] neighbors, long seed) {
    Random random = new Random(seed);
    List<int[][]> variants = new ArrayList<>();
    variants.add(neighbors);
    variants.add(sortedNeighbours(neighbors));
    for (int i = 0; i < 3; i++) variants.add(shuffledNeighbours(neighbors, random));
    for (int i = 0; i < 3; i++) variants.add(relabelled(neighbors, random));
    return variants;
  }

  /** Field-by-field equality, with the order of every ring list and of the atoms of every ring. */
  static void assertSameRings(RingPerception.Rings expected, RingPerception.Rings actual, String label) {
    assertArrayEquals(expected.smallest, actual.smallest, label + ": basis");
    assertArrayEquals(expected.relevant, actual.relevant, label + ": relevant");
    assertArrayEquals(expected.legacyRelevant, actual.legacyRelevant, label + ": legacy relevant");
    assertArrayEquals(expected.numericBasis, actual.numericBasis, label + ": numeric basis");
    assertEquals(expected.basisTruncated, actual.basisTruncated, label + ": basis truncated");
    assertEquals(expected.relevantTruncated, actual.relevantTruncated, label + ": relevant truncated");
    assertEquals(expected.truncated, actual.truncated, label + ": truncated");
  }

  static int weight(int[][] cycles) {
    int total = 0;
    for (int[] cycle : cycles) total += cycle.length;
    return total;
  }
}
