/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
package com.bioinception.smsd.core;

import static com.bioinception.smsd.TestSupport.rawMolecule;
import static org.junit.jupiter.api.Assertions.*;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.concurrent.atomic.AtomicReference;
import org.junit.jupiter.api.Test;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.silent.SilentChemObjectBuilder;

/** SMARTS input limits shared with the native core, and ring perception of very large rings. */
class InputLimitsRegressionTest {
  private static final long SMALL_STACK = 512 * 1024;

  private static String repeat(String text, int count) {
    return String.join("", Collections.nCopies(count, text));
  }

  /** Runs work on a thread with a 512 KB stack and rethrows what it threw. */
  private static void onSmallStack(ThrowingRunnable work) throws Throwable {
    AtomicReference<Throwable> failure = new AtomicReference<>();
    Thread thread = new Thread(null, () -> {
      try {
        work.run();
      } catch (Throwable thrown) {
        failure.set(thrown);
      }
    }, "small-stack", SMALL_STACK);
    thread.start();
    thread.join();
    if (failure.get() != null) throw failure.get();
  }

  private interface ThrowingRunnable {
    void run() throws Exception;
  }

  /** A chain of n carbons, built without parsing SMILES. */
  private static MolGraph chain(int n) {
    List<int[]> edges = new ArrayList<>();
    for (int i = 1; i < n; i++) edges.add(new int[] {i - 1, i});
    return graph(n, edges);
  }

  /** A ring of n carbons, optionally fused to a benzene ring, built without parsing SMILES. */
  private static MolGraph ring(int n, boolean fused) {
    List<int[]> edges = new ArrayList<>();
    int ringAtoms = fused ? n - 4 : n;
    for (int i = 0; i < ringAtoms; i++) edges.add(new int[] {i, (i + 1) % ringAtoms});
    if (fused) {
      edges.add(new int[] {0, n - 4});
      edges.add(new int[] {n - 4, n - 3});
      edges.add(new int[] {n - 3, n - 2});
      edges.add(new int[] {n - 2, n - 1});
      edges.add(new int[] {n - 1, 1});
    }
    return graph(n, edges);
  }

  private static MolGraph graph(int n, List<int[]> edges) {
    List<List<Integer>> rows = new ArrayList<>();
    for (int i = 0; i < n; i++) rows.add(new ArrayList<>());
    for (int[] edge : edges) {
      rows.get(edge[0]).add(edge[1]);
      rows.get(edge[1]).add(edge[0]);
    }
    int[][] neighbors = new int[n][], orders = new int[n][];
    for (int i = 0; i < n; i++) {
      neighbors[i] = rows.get(i).stream().mapToInt(Integer::intValue).toArray();
      orders[i] = new int[neighbors[i].length];
      Arrays.fill(orders[i], 1);
    }
    int[] elements = new int[n];
    Arrays.fill(elements, 6);
    return new MolGraph.Builder().atomCount(n).atomicNumbers(elements).neighbors(neighbors).bondOrders(orders).build();
  }

  @Test
  void smartsLimitsAreInclusiveAndConsistent() throws Exception {
    assertEquals(100_000, Substructure.MAX_SMARTS_LENGTH);
    assertEquals(10_000, Substructure.MAX_SMARTS_ATOMS);
    assertEquals(100_000, Substructure.MAX_SMARTS_EXPRESSION_NODES);
    assertEquals(10_000, Substructure.MAX_SMARTS_NODES_PER_EXPRESSION);
    assertEquals(20, Substructure.MAX_SMARTS_RECURSION_DEPTH);
    assertTrue(Substructure.smarts(repeat("[$(", 20) + "C" + repeat(")]", 20)).match(rawMolecule("C")).found());
    assertThrows(IllegalArgumentException.class, () -> Substructure.smarts(repeat("[$(", 21) + "C" + repeat(")]", 21)));
    assertTrue(Substructure.smarts("[" + repeat("!", 99_997) + "C]").match(rawMolecule("N")).found());
    IllegalArgumentException length = assertThrows(IllegalArgumentException.class,
        () -> Substructure.smarts("[" + repeat("!", 99_998) + "C]"));
    assertTrue(length.getMessage().contains("maximum length"), length.getMessage());
    Substructure largest = Substructure.smarts(repeat("C", 10_000));
    assertFalse(largest.match(chain(20)).found());
    assertThrows(IllegalArgumentException.class, () -> Substructure.smarts(repeat("C", 10_001)));
    assertThrows(IllegalArgumentException.class, () -> Substructure.smarts(repeat("[Cl]", 10_001)));
    assertDoesNotThrow(() -> Substructure.smarts(repeat("Cl", 5_000) + repeat("Br", 5_000)));
    assertThrows(IllegalArgumentException.class, () -> Substructure.smarts("[$(" + repeat("C", 10_000) + ")]"));
    // [C R R ...] has one primitive per symbol and one implicit AND between neighbours, so
    // 4999 Rs make 9999 nodes, counted as the native parser counts them.
    String ringAtom = "[C" + repeat("R", 4_999) + "]";
    assertTrue(Substructure.smarts(ringAtom).match(rawMolecule("C1CC1")).found());
    IllegalArgumentException expression = assertThrows(IllegalArgumentException.class,
        () -> Substructure.smarts("[C" + repeat("R", 5_000) + "]"));
    assertTrue(expression.getMessage().contains("atom or bond expression"), expression.getMessage());
    assertDoesNotThrow(() -> Substructure.smarts(repeat(ringAtom, 10)));
    IllegalArgumentException nodes = assertThrows(IllegalArgumentException.class,
        () -> Substructure.smarts(repeat(ringAtom, 10) + ".[CR]"));
    assertTrue(nodes.getMessage().contains("expression nodes"), nodes.getMessage());
    StringBuilder nested = new StringBuilder("[");
    for (int part = 0; part < 11; part++) nested.append(part == 0 ? "$(" : ";$(").append(ringAtom).append(')');
    assertThrows(IllegalArgumentException.class, () -> Substructure.smarts(nested + "]"));
  }

  /** Long operator chains within the limits parse and match without deep recursion. */
  @Test
  void longExpressionsMatchOnASmallStack() throws Throwable {
    IAtomContainer chain = rawMolecule("CCCC");
    onSmallStack(() -> {
      // 5,000 operands joined by 4,999 operators make 9,999 nodes.
      for (String operator : new String[] {"&", ";", ",", ""}) {
        String atoms = "[" + String.join(operator, Collections.nCopies(5_000, "C")) + "]";
        assertTrue(Substructure.smarts(atoms).match(chain).found(), operator);
      }
      String bond = "C" + String.join("&", Collections.nCopies(5_000, "-")) + "C";
      assertTrue(Substructure.smarts(bond).match(chain).found());
      assertTrue(Substructure.smarts("[" + repeat("!", 60_001) + "N]").match(chain).found());
      assertFalse(Substructure.smarts("[" + String.join("&", Collections.nCopies(4_999, "C")) + "&N]").match(chain).found());
    });
  }

  /** Ring perception of very large rings needs no deep call stack and bounded work. */
  @Test
  void veryLargeRingsArePerceivedOnASmallStack() throws Throwable {
    onSmallStack(() -> {
      long started = System.nanoTime();
      MolGraph macrocycle = ring(9_999, false);
      assertEquals(1, macrocycle.computeRings().length);
      assertEquals(1, macrocycle.computeRelevantCycles().length);
      MolGraph fused = ring(3_000, true);
      assertEquals(2, fused.computeRings().length);
      assertEquals(2, fused.computeRelevantCycles().length);
      // The work bound keeps rings of thousands of atoms from taking minutes.
      assertTrue(System.nanoTime() - started < 20_000_000_000L, "large-ring perception took too long");
    });
  }

  /** Ring flags of a long chain and a macrocycle are marked without deep recursion, as CDK marks them. */
  @Test
  void ringFlagsOfLargeMoleculesNeedNoDeepRecursion() throws Throwable {
    IAtomContainer chain = SilentChemObjectBuilder.getInstance().newAtomContainer();
    for (int atom = 0; atom < 9_999; atom++) {
      chain.addAtom(SilentChemObjectBuilder.getInstance().newInstance(IAtom.class, "C"));
      chain.getAtom(atom).setImplicitHydrogenCount(atom == 0 || atom == 9_998 ? 3 : 2);
      if (atom > 0) chain.addBond(atom - 1, atom, IBond.Order.SINGLE);
    }
    IAtomContainer macrocycle = chain.clone();
    macrocycle.addBond(0, 9_998, IBond.Order.SINGLE);
    macrocycle.getAtom(0).setImplicitHydrogenCount(2);
    macrocycle.getAtom(9_998).setImplicitHydrogenCount(2);
    onSmallStack(() -> {
      MolGraph open = new MolGraph(chain);
      for (IAtom atom : chain.atoms()) assertFalse(atom.isInRing());
      MolGraph closed = new MolGraph(macrocycle);
      for (IAtom atom : macrocycle.atoms()) assertTrue(atom.isInRing());
      for (IBond bond : macrocycle.bonds()) assertTrue(bond.isInRing());
      Substructure.Options options = new Substructure.Options();
      options.timeoutMillis = 10_000;
      assertTrue(Substructure.smarts("[R]").withOptions(options).match(closed).found());
      assertFalse(Substructure.smarts("[R]").withOptions(options).match(open).found());
    });
  }
}
