/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
package com.bioinception.smsd.core;

import static org.junit.jupiter.api.Assertions.*;

import java.lang.reflect.Field;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.atomic.AtomicReference;
import org.junit.jupiter.api.Test;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.isomorphism.Pattern;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;

class SmilesWriterSafetyRegressionTest {
  @Test
  void moreThanNinetyNineRingsRoundTrip() throws Exception {
    StringBuilder input = new StringBuilder();
    for (int ring = 0; ring < 101; ring++) input.append("C(C1CC1)");
    IAtomContainer original = molecule(input.toString());
    MolGraph graph = new MolGraph(original);
    assertSameTopology(original, molecule(graph.toCanonicalSmiles()));

    Map<Integer, Integer> mapping = new LinkedHashMap<>();
    for (int atom = 0; atom < graph.atomCount(); atom++) mapping.put(atom, atom);
    Map<SearchResult.Bond, SearchResult.Bond> bonds = MatchGraph.allBonds(graph, mapping);
    SearchResult.Match fragment = new SearchResult.Match(mapping, bonds,
        MatchGraph.fragments(mapping, bonds.keySet()), mapping.size());
    assertSameTopology(original, molecule(fragment.toSmiles(graph)));
  }

  @Test
  void ringLabelsAreReleasedBetweenComponents() throws Exception {
    assertEquals("C1CC1.C2CC2", new MolGraph(molecule("C1CC1.C1CC1")).toCanonicalSmiles());
  }

  @Test
  void ninetyNineSimultaneousRingClosuresRoundTrip() throws Exception {
    MolGraph graph = parallelClosures(99);
    useIdentityWriterOrder(graph);
    String output = graph.toCanonicalSmiles();
    assertTrue(output.contains("%99"));
    assertSameTopology(molecule(parallelClosureSmiles(99)), molecule(output));
  }

  @Test
  void excessiveSimultaneousRingClosuresAreRejected() throws Exception {
    MolGraph graph = parallelClosures(100);
    useIdentityWriterOrder(graph);
    IllegalArgumentException failure = assertThrows(IllegalArgumentException.class,
        graph::toCanonicalSmiles);
    assertEquals("SMILES output needs more than 99 simultaneously open ring closures", failure.getMessage());
  }

  @Test
  void ordinaryChainRoundTripsThroughCanonicalPreparation() throws Exception {
    char[] atoms = new char[128];
    Arrays.fill(atoms, 'C');
    IAtomContainer original = molecule(new String(atoms));
    assertSameTopology(original, molecule(new MolGraph(original).toCanonicalSmiles()));
  }

  @Test
  void writerTraversesLongChainWithSmallThreadStack() throws Exception {
    MolGraph graph = chain(9_999);
    useIdentityWriterOrder(graph);
    AtomicReference<Throwable> failure = new AtomicReference<>();
    AtomicReference<String> result = new AtomicReference<>();
    Thread writer = new Thread(null, () -> {
      try {
        result.set(graph.toCanonicalSmiles());
      } catch (Throwable caught) {
        failure.set(caught);
      }
    }, "smiles-writer", 512L * 1024);
    writer.start();
    writer.join();
    assertNull(failure.get(), () -> String.valueOf(failure.get()));
    char[] expected = new char[graph.atomCount()];
    Arrays.fill(expected, 'C');
    assertEquals(new String(expected), result.get());
  }

  /** Fix the writer order without invoking the separately tested canonical labeller. */
  private static void useIdentityWriterOrder(MolGraph graph) throws Exception {
    graph.canonicalLabel = new int[graph.atomCount()];
    for (int atom = 0; atom < graph.atomCount(); atom++) graph.canonicalLabel[atom] = atom;
    Field computed = MolGraph.class.getDeclaredField("canonicalComputed");
    computed.setAccessible(true);
    computed.setBoolean(graph, true);
  }

  private static MolGraph chain(int atomCount) {
    int[][] neighbors = new int[atomCount][];
    for (int atom = 0; atom < atomCount; atom++) {
      if (atom == 0) neighbors[atom] = new int[] {1};
      else if (atom == atomCount - 1) neighbors[atom] = new int[] {atom - 1};
      else neighbors[atom] = new int[] {atom - 1, atom + 1};
    }
    return carbonGraph(neighbors);
  }

  private static MolGraph parallelClosures(int ringCount) {
    int atomCount = 2 * ringCount;
    @SuppressWarnings("unchecked")
    List<Integer>[] neighbors = new List[atomCount];
    for (int atom = 0; atom < atomCount; atom++) neighbors[atom] = new ArrayList<>();
    for (int atom = 1; atom < atomCount; atom++) addBond(neighbors, atom - 1, atom);
    for (int ring = 0; ring < ringCount; ring++) addBond(neighbors, ring, ring + ringCount);
    int[][] rows = new int[atomCount][];
    for (int atom = 0; atom < atomCount; atom++) {
      rows[atom] = new int[neighbors[atom].size()];
      for (int neighbor = 0; neighbor < rows[atom].length; neighbor++)
        rows[atom][neighbor] = neighbors[atom].get(neighbor);
    }
    return carbonGraph(rows);
  }

  private static void addBond(List<Integer>[] neighbors, int firstAtom, int secondAtom) {
    neighbors[firstAtom].add(secondAtom);
    neighbors[secondAtom].add(firstAtom);
  }

  private static MolGraph carbonGraph(int[][] neighbors) {
    int[] atomicNumbers = new int[neighbors.length];
    Arrays.fill(atomicNumbers, 6);
    return new MolGraph.Builder().atomCount(neighbors.length).atomicNumbers(atomicNumbers)
        .neighbors(neighbors).build();
  }

  private static String parallelClosureSmiles(int ringCount) {
    StringBuilder text = new StringBuilder();
    for (int atom = 0; atom < 2 * ringCount; atom++) {
      text.append('C');
      int label = atom % ringCount + 1;
      if (label >= 10) text.append('%');
      text.append(label);
    }
    return text.toString();
  }

  private static IAtomContainer molecule(String text) throws Exception {
    return new SmilesParser(SilentChemObjectBuilder.getInstance()).parseSmiles(text);
  }

  private static void assertSameTopology(IAtomContainer expected, IAtomContainer actual) {
    assertEquals(expected.getAtomCount(), actual.getAtomCount());
    assertEquals(expected.getBondCount(), actual.getBondCount());
    assertTrue(Pattern.findIdentical(expected).matches(actual));
  }
}
