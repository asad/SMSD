/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
package com.bioinception.smsd.core;

import static org.junit.jupiter.api.Assertions.*;

import java.lang.reflect.Field;
import java.util.Arrays;
import java.util.Set;
import java.util.TreeSet;
import org.junit.jupiter.api.Test;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;

/** The writer refusals and the SMARTS inputs that the native tests check, with the same outcomes. */
class WriterAndSmartsParityTest {
  private static final String[] BONDS = {"-", "=", "#", ":", "~", "@", "/", "\\", "/?", "\\?"};
  /** The message these inputs have always had in Java, kept word for word. */
  private static final String UNFINISHED = "Invalid SMARTS: Unclosed ring, component group, or branch";

  @Test
  void ninetyNineSimultaneouslyOpenSmartsClosuresRoundTrip() {
    MolGraph graph = parallelClosures(99);
    String text = SmartsWriter.toSmarts(graph);
    assertTrue(text.contains("%99"), text);
    // The chain is written in atom order, so the parsed atoms keep the source numbering.
    IAtomContainer query = SmartsInput.parse(text);
    assertEquals(edges(graph), edges(query));
  }

  @Test
  void excessiveSimultaneouslyOpenSmartsClosuresAreRejected() {
    MolGraph graph = parallelClosures(100);
    IllegalArgumentException failure = assertThrows(IllegalArgumentException.class,
        () -> SmartsWriter.toSmarts(graph));
    assertEquals("SMARTS output needs more than 99 simultaneously open ring closures", failure.getMessage());
  }

  @Test
  void excessiveSimultaneouslyOpenSmilesClosuresAreRejected() throws Exception {
    MolGraph graph = parallelClosures(100);
    useIdentityWriterOrder(graph);
    IllegalArgumentException failure = assertThrows(IllegalArgumentException.class,
        graph::toCanonicalSmiles);
    assertEquals("SMILES output needs more than 99 simultaneously open ring closures", failure.getMessage());
  }

  /** More than 99 closures in all, but never more than one open, are written with reused numbers. */
  @Test
  void manySmartsClosuresReuseNumbers() {
    MolGraph graph = pendantRings(101);
    String text = SmartsWriter.toSmarts(graph);
    assertFalse(text.contains("%"), text);
    assertSameShape(graph, text);
    assertTrue(Substructure.smarts(text).match(graph).found());
  }

  @Test
  void trailingBondsAreRejected() {
    for (String bond : BONDS) {
      String dangling = "C" + bond;
      String recursive = "[$(" + dangling + ")]";
      assertRejected(dangling, UNFINISHED);
      assertRejected(dangling + " \t\r\n", UNFINISHED);
      assertRejected("C" + dangling, UNFINISHED);
      assertRejected(recursive, null);
      assertRejected("[$(" + recursive + ")]", null);

      String complete = bond.equals(":") ? "c:c" : "C" + bond + "C";
      String completeRecursive = "[$(" + complete + ")]";
      assertShape(complete, 2, 1);
      assertShape(complete + " \t\r\n", 2, 1);
      assertShape(completeRecursive, 1, 0);
      assertShape("[$(" + completeRecursive + ")]", 1, 0);
    }
    for (String unfinished : new String[] {"CC-,=", "CC!-", "CC-;@"}) assertRejected(unfinished, UNFINISHED);
  }

  @Test
  void otherUnfinishedPatternsKeepTheirMessages() {
    assertRejected("C(", UNFINISHED);
    assertRejected("C1CC", UNFINISHED);
    assertRejected("[C", "Invalid SMARTS: Invalid atom expression");
    assertRejected("C)C", "Invalid SMARTS: Closing unopened component grouping");
  }

  @Test
  void bondsAtRingClosuresAreNotTrailing() {
    for (String closure : new String[] {"C-1CC1", "C1CC-1", "C=1CC1", "C1CC=1"}) assertShape(closure, 3, 3);
    assertShape("C", 1, 0);
    assertShape("", 0, 0);
  }

  private static void assertRejected(String smarts, String message) {
    IllegalArgumentException failure = assertThrows(IllegalArgumentException.class,
        () -> SmartsInput.parse(smarts), smarts);
    if (message != null) assertEquals(message, failure.getMessage(), smarts);
  }

  private static void assertShape(String smarts, int atoms, int bonds) {
    IAtomContainer query = SmartsInput.parse(smarts);
    assertEquals(atoms, query.getAtomCount(), smarts);
    assertEquals(bonds, query.getBondCount(), smarts);
  }

  private static Set<Long> edges(MolGraph graph) {
    Set<Long> edges = new TreeSet<>();
    for (int atom = 0; atom < graph.atomCount(); atom++)
      for (int neighbor : graph.neighbors[atom]) edges.add(MolGraph.bondKey(atom, neighbor));
    return edges;
  }

  private static Set<Long> edges(IAtomContainer query) {
    Set<Long> edges = new TreeSet<>();
    for (IBond bond : query.bonds())
      edges.add(MolGraph.bondKey(query.indexOf(bond.getBegin()), query.indexOf(bond.getEnd())));
    return edges;
  }

  private static void assertSameShape(MolGraph graph, String smarts) {
    int edges = 0;
    for (int atom = 0; atom < graph.atomCount(); atom++) edges += graph.neighbors[atom].length;
    assertShape(smarts, graph.atomCount(), edges / 2);
  }

  /** Fix the writer order without invoking the separately tested canonical labeller. */
  private static void useIdentityWriterOrder(MolGraph graph) throws Exception {
    graph.canonicalLabel = new int[graph.atomCount()];
    for (int atom = 0; atom < graph.atomCount(); atom++) graph.canonicalLabel[atom] = atom;
    Field computed = MolGraph.class.getDeclaredField("canonicalComputed");
    computed.setAccessible(true);
    computed.setBoolean(graph, true);
  }

  /** A chain of 2k atoms with atom i also bonded to atom i + k, so k closures are open at once. */
  private static MolGraph parallelClosures(int ringCount) {
    int[][] neighbors = new int[2 * ringCount][0];
    for (int atom = 1; atom < 2 * ringCount; atom++) addBond(neighbors, atom - 1, atom);
    for (int ring = 0; ring < ringCount; ring++) addBond(neighbors, ring, ring + ringCount);
    return carbonGraph(neighbors);
  }

  /** A chain of k atoms, each carrying its own three-membered ring. */
  private static MolGraph pendantRings(int count) {
    int[][] neighbors = new int[4 * count][0];
    for (int atom = 1; atom < count; atom++) addBond(neighbors, atom - 1, atom);
    for (int ring = 0; ring < count; ring++) {
      int first = count + 3 * ring;
      addBond(neighbors, ring, first);
      addBond(neighbors, first, first + 1);
      addBond(neighbors, first + 1, first + 2);
      addBond(neighbors, first + 2, first);
    }
    return carbonGraph(neighbors);
  }

  private static void addBond(int[][] neighbors, int first, int second) {
    neighbors[first] = Arrays.copyOf(neighbors[first], neighbors[first].length + 1);
    neighbors[first][neighbors[first].length - 1] = second;
    neighbors[second] = Arrays.copyOf(neighbors[second], neighbors[second].length + 1);
    neighbors[second][neighbors[second].length - 1] = first;
  }

  private static MolGraph carbonGraph(int[][] neighbors) {
    int[] atomicNumbers = new int[neighbors.length];
    Arrays.fill(atomicNumbers, 6);
    return new MolGraph.Builder().atomCount(neighbors.length).atomicNumbers(atomicNumbers)
        .neighbors(neighbors).build();
  }
}
