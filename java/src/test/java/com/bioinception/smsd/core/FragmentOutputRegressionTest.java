/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
package com.bioinception.smsd.core;

import static org.junit.jupiter.api.Assertions.*;

import java.util.LinkedHashMap;
import java.util.Map;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.ValueSource;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.isomorphism.Pattern;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;

class FragmentOutputRegressionTest {
  @ParameterizedTest
  @ValueSource(strings = {"[nH]1cccc1", "c1ncc[nH]1", "c1ccc2[nH]ccc2c1", "[NH4+]",
      "[NH3+]CC(=O)[O-]", "[nH+]1ccccc1", "[SiH4]", "C[SiH2]C", "[2H]", "[CH3]"})
  void canonicalSmilesPreservesHydrogens(String text) throws Exception {
    SmilesParser parser = new SmilesParser(SilentChemObjectBuilder.getInstance());
    IAtomContainer original = parser.parseSmiles(text);
    String output = new MolGraph(original).toCanonicalSmiles();
    IAtomContainer restored = parser.parseSmiles(output);
    assertEquals(hydrogenCount(original), hydrogenCount(restored), output);
    assertTrue(Pattern.findIdentical(original).matches(restored), output);
  }

  @ParameterizedTest
  @ValueSource(strings = {"c1ccccc1", "c1ncc[nH]1", "c1ccc2[nH]ccc2c1",
      "O=c1cc[nH]c(=O)[nH]1", "Cn1cccc1", "[NH4+]", "[CH3]", "C/C=C/C", "F[C@H](Cl)CCO"})
  void completeMCSKeepsChemistryAndStereo(String text) throws Exception {
    SmilesParser parser = new SmilesParser(SilentChemObjectBuilder.getInstance());
    IAtomContainer original = parser.parseSmiles(text);
    Map<Integer, Integer> mapping = new LinkedHashMap<>();
    for (int atom = 0; atom < original.getAtomCount(); atom++) mapping.put(atom, atom);
    String output = fragmentSmiles(new MolGraph(original), mapping);
    IAtomContainer restored = parser.parseSmiles(output);
    assertEquals(hydrogenCount(original), hydrogenCount(restored), output);
    assertEquals(aromaticCount(original), aromaticCount(restored), output);
    assertTrue(Pattern.findIdentical(original).matches(restored), output);
  }

  @Test
  void retainedPhenylRingIsBenzene() throws Exception {
    SmilesParser parser = new SmilesParser(SilentChemObjectBuilder.getInstance());
    IAtomContainer original = parser.parseSmiles("Cc1ccccc1");
    Map<Integer, Integer> mapping = new LinkedHashMap<>();
    for (int atom = 1; atom <= 6; atom++) mapping.put(atom, atom - 1);
    String output = fragmentSmiles(new MolGraph(original), mapping);
    IAtomContainer restored = parser.parseSmiles(output);
    assertTrue(Pattern.findIdentical(parser.parseSmiles("c1ccccc1")).matches(restored), output);
    assertEquals(6, hydrogenCount(restored), output);
  }

  @Test
  void openedAromaticRingHasCappedValences() throws Exception {
    SmilesParser parser = new SmilesParser(SilentChemObjectBuilder.getInstance());
    Map<Integer, Integer> mapping = new LinkedHashMap<>();
    for (int atom = 0; atom < 5; atom++) mapping.put(atom, atom);
    String output = fragmentSmiles(new MolGraph(parser.parseSmiles("c1ccccc1")), mapping);
    IAtomContainer restored = parser.parseSmiles(output);
    assertEquals(0, aromaticCount(restored), output);
    for (IAtom atom : restored.atoms())
      assertEquals(4, atom.getImplicitHydrogenCount() + (int) restored.getBondOrderSum(atom), output);
  }

  @Test
  void cutAlkylBondIsCappedWithHydrogen() throws Exception {
    SmilesParser parser = new SmilesParser(SilentChemObjectBuilder.getInstance());
    Map<Integer, Integer> mapping = new LinkedHashMap<>();
    for (int atom = 1; atom <= 3; atom++) mapping.put(atom, atom - 1);
    String output = fragmentSmiles(new MolGraph(parser.parseSmiles("CCCO")), mapping);
    IAtomContainer restored = parser.parseSmiles(output);
    assertTrue(Pattern.findIdentical(parser.parseSmiles("CCO")).matches(restored), output);
    assertEquals(6, hydrogenCount(restored), output);
  }

  @Test
  void retainedFusedRingIsBenzene() throws Exception {
    SmilesParser parser = new SmilesParser(SilentChemObjectBuilder.getInstance());
    MolGraph original = new MolGraph(parser.parseSmiles("c1ccc2ccccc2c1"));
    Map<Integer, Integer> mapping = new LinkedHashMap<>();
    for (int atom : original.computeRings()[0]) mapping.put(atom, atom);
    String output = fragmentSmiles(original, mapping);
    IAtomContainer restored = parser.parseSmiles(output);
    assertEquals(6, aromaticCount(restored), output);
    assertEquals(6, hydrogenCount(restored), output);
  }

  @Test
  void builderRetainsSuppliedHydrogensAndDefaultValence() throws Exception {
    MolGraph inferred = new MolGraph.Builder().atomCount(1).atomicNumbers(new int[] {6})
        .neighbors(new int[][] {new int[0]}).build();
    assertEquals("C", inferred.toCanonicalSmiles());
    MolGraph explicitZero = new MolGraph.Builder().atomCount(1).atomicNumbers(new int[] {6})
        .hydrogenCounts(new int[] {0}).neighbors(new int[][] {new int[0]}).build();
    assertEquals("[C]", explicitZero.toCanonicalSmiles());
    MolGraph ammonium = new MolGraph.Builder().atomCount(1).atomicNumbers(new int[] {7})
        .formalCharges(new int[] {1}).neighbors(new int[][] {new int[0]}).build();
    assertEquals("[NH4+]", ammonium.toCanonicalSmiles());
    MolGraph pyrrole = new MolGraph.Builder().atomCount(5).atomicNumbers(new int[] {7, 6, 6, 6, 6})
        .hydrogenCounts(new int[] {1, 1, 1, 1, 1}).aromaticFlags(new boolean[] {true, true, true, true, true})
        .neighbors(new int[][] {{1, 4}, {0, 2}, {1, 3}, {2, 4}, {3, 0}})
        .bondOrders(new int[][] {{4, 4}, {4, 4}, {4, 4}, {4, 4}, {4, 4}})
        .bondAromaticFlags(new boolean[][] {{true, true}, {true, true}, {true, true}, {true, true}, {true, true}})
        .build();
    SmilesParser parser = new SmilesParser(SilentChemObjectBuilder.getInstance());
    assertEquals(5, hydrogenCount(parser.parseSmiles(pyrrole.toCanonicalSmiles())));
    Map<Integer, Integer> mapping = new LinkedHashMap<>();
    for (int atom = 0; atom < 5; atom++) mapping.put(atom, atom);
    assertEquals(5, hydrogenCount(parser.parseSmiles(fragmentSmiles(pyrrole, mapping))));
  }

  @Test
  void builderValidatesSuppliedHydrogenCounts() {
    assertThrows(IllegalArgumentException.class, () -> new MolGraph.Builder().atomCount(1)
        .atomicNumbers(new int[] {6}).hydrogenCounts(new int[] {0, 0}).neighbors(new int[][] {new int[0]}).build());
    assertThrows(IllegalArgumentException.class, () -> new MolGraph.Builder().atomCount(1)
        .atomicNumbers(new int[] {6}).hydrogenCounts(new int[] {-1}).neighbors(new int[][] {new int[0]}).build());
  }

  private static int hydrogenCount(IAtomContainer molecule) {
    int count = 0;
    for (IAtom atom : molecule.atoms()) count += atom.getImplicitHydrogenCount();
    return count;
  }

  private static int aromaticCount(IAtomContainer molecule) {
    int count = 0;
    for (IAtom atom : molecule.atoms()) if (atom.isAromatic()) count++;
    return count;
  }
  private static String fragmentSmiles(MolGraph graph, Map<Integer, Integer> mapping) {
    Map<SearchResult.Bond, SearchResult.Bond> bonds = MCS.allBonds(graph, mapping);
    return new SearchResult.Match(mapping, bonds, MCS.fragments(mapping, bonds.keySet()), mapping.size()).toSmiles(graph);
  }
}
