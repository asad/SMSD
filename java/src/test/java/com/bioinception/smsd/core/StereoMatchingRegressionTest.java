/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
package com.bioinception.smsd.core;

import static org.junit.jupiter.api.Assertions.*;

import java.util.Arrays;
import java.util.Map;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.CsvSource;
import org.junit.jupiter.params.provider.ValueSource;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.isomorphism.Pattern;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;

class StereoMatchingRegressionTest {
  private static IAtomContainer molecule(String smiles) throws Exception {
    return new SmilesParser(SilentChemObjectBuilder.getInstance()).parseSmiles(smiles);
  }

  @ParameterizedTest
  @CsvSource({
      "C[C@H](O)CC, ClC[C@H](O)CC, true",
      "C[C@@H](O)CC, ClC[C@H](O)CC, false",
      "N[C@@H](C)C(=O)O, N[C@@H](CS)C(=O)O, true",
      "N[C@@H](C)C(=O)O, N[C@@](F)(C)C(=O)O, true",
      "N[C@@H](C)C(=O)O, N[C@](F)(C)C(=O)O, false"
  })
  void tetrahedralMatchingUsesMappedLigands(String queryText, String targetText, boolean expected) throws Exception {
    ChemOptions chemistry = new ChemOptions();
    chemistry.useChirality = true;
    assertEntryPoints(queryText, targetText, chemistry, expected);
  }

  @ParameterizedTest
  @CsvSource({
      "C/C=C/C, Cl/C(C)=C/C, false",
      "C/C=C\\C, Cl/C(C)=C/C, true",
      "C/C(Cl)=C/C, Cl/C(C)=C\\C, true",
      "Cl/C(C)=C\\C, C/C(Cl)=C/C, true",
      "F/C=C/Cl, F/C=C\\Cl, false"
  })
  void doubleBondMatchingUsesMappedReferences(String queryText, String targetText, boolean expected) throws Exception {
    ChemOptions chemistry = new ChemOptions();
    chemistry.useBondStereo = true;
    assertEntryPoints(queryText, targetText, chemistry, expected);
  }

  private static void assertEntryPoints(String queryText, String targetText, ChemOptions chemistry,
                                       boolean expected) throws Exception {
    IAtomContainer queryMolecule = molecule(queryText), targetMolecule = molecule(targetText);
    assertEquals(expected, Pattern.findSubstructure(queryMolecule).matches(targetMolecule));
    MolGraph query = new MolGraph(queryMolecule), target = new MolGraph(targetMolecule);
    for (boolean bitParallel : new boolean[] {false, true}) {
      chemistry.useBitParallelFeasibility = bitParallel;
      Substructure matcher = Substructure.find(query).withOptions(substructureOptions(chemistry)).withTimeout(2_000);
      SearchResult first = matcher.match(target);
      assertEquals(expected, first.found());
      assertEquals(expected, !first.matches().isEmpty());
      assertEquals(expected, matcher.matchAll(target, 20).found());
    }
  }

  @ParameterizedTest
  @CsvSource({"F/C=C(Cl)/C, Z", "F/C=C(C)/Cl, E", "OC(=O)/C=C(C)/C(=O)O, E"})
  void ezAssignmentUsesTheStoredReferenceConvention(String text, char expected) throws Exception {
    IAtomContainer original = molecule(text);
    for (int pass = 0; pass < 2; pass++) {
      Map<Long, Character> labels = CIPAssigner.assignEZ(new MolGraph(original));
      assertEquals(1, labels.size());
      assertEquals(expected, labels.values().iterator().next());
      IBond[] reversed = new IBond[original.getBondCount()];
      for (int bond = 0; bond < reversed.length; bond++) reversed[bond] = original.getBond(reversed.length - 1 - bond);
      original.setBonds(reversed);
    }
  }

  @ParameterizedTest
  @ValueSource(ints = {195, 196, 197, 200})
  void largeGraphsRetainDoubleBondStereo(int extraAtoms) throws Exception {
    StringBuilder targetText = new StringBuilder("C/C=C/C");
    for (int atom = 0; atom < extraAtoms; atom++) targetText.append('C');
    MolGraph target = new MolGraph(molecule(targetText.toString()));
    ChemOptions chemistry = new ChemOptions();
    chemistry.useBondStereo = true;
    assertTrue(Substructure.find(new MolGraph(molecule("C/C=C/C"))).withOptions(substructureOptions(chemistry)).withTimeout(2_000).match(target).found());
    assertFalse(Substructure.find(new MolGraph(molecule("C/C=C\\C"))).withOptions(substructureOptions(chemistry)).withTimeout(2_000).match(target).found());
    assertEquals('E', CIPAssigner.assignEZ(target).get(MolGraph.bondKey(1, 2)));
  }

  @Test
  void builderStereoSurvivesSparseStorageAndWriting() throws Exception {
    int atoms = 204;
    int[] elements = new int[atoms];
    Arrays.fill(elements, 6);
    int[][] neighbors = new int[atoms][], bondOrders = new int[atoms][];
    for (int atom = 0; atom < atoms; atom++) {
      neighbors[atom] = atom == 0 ? new int[] {1} : atom == atoms - 1 ? new int[] {atom - 1} : new int[] {atom - 1, atom + 1};
      bondOrders[atom] = new int[neighbors[atom].length];
      for (int position = 0; position < neighbors[atom].length; position++)
        bondOrders[atom][position] = MolGraph.bondKey(atom, neighbors[atom][position]) == MolGraph.bondKey(1, 2) ? 2 : 1;
    }
    int[][] stereo = new int[atoms][atoms];
    stereo[1][2] = stereo[2][1] = 2;
    MolGraph target = new MolGraph.Builder().atomCount(atoms).atomicNumbers(elements)
        .neighbors(neighbors).bondOrders(bondOrders).doubleBondStereo(stereo).build();
    ChemOptions chemistry = new ChemOptions();
    chemistry.useBondStereo = true;
    assertFalse(Substructure.find(new MolGraph(molecule("C/C=C\\C"))).withOptions(substructureOptions(chemistry)).withTimeout(2_000).match(target).found());
    assertEquals('E', CIPAssigner.assignEZ(target).get(MolGraph.bondKey(1, 2)));
    assertEquals('E', CIPAssigner.assignEZ(new MolGraph(molecule(target.toCanonicalSmiles()))).values().iterator().next());
  }

  @Test
  void mcsValidationUsesMappedStereo() throws Exception {
    ChemOptions chemistry = new ChemOptions();
    chemistry.useChirality = true;
    MolGraph query = new MolGraph(molecule("C[C@H](O)CC"));
    MolGraph target = new MolGraph(molecule("ClC[C@H](O)CC"));
    assertEquals(5, firstMapping(MCS.find(query).withOptions(mcsOptions(chemistry)).match(target)).size());
    chemistry.useChirality = false;
    chemistry.useBondStereo = true;
    query = new MolGraph(molecule("C/C=C/C"));
    target = new MolGraph(molecule("Cl/C(C)=C/C"));
    Map<Integer, Integer> mapping = firstMapping(MCS.find(query).withOptions(mcsOptions(chemistry)).match(target));
    assertTrue(mapping.size() < 4);
    int[] indexed = new int[query.n];
    Arrays.fill(indexed, -1);
    for (Map.Entry<Integer, Integer> atom : mapping.entrySet()) indexed[atom.getKey()] = atom.getValue();
    assertEquals(-1, Substructure.stereoMismatchAtom(query, target, indexed, chemistry));
    query = new MolGraph(molecule("C/C=C\\C"));
    assertEquals(4, firstMapping(MCS.find(query).withOptions(mcsOptions(chemistry)).match(target)).size());
  }
  private static Substructure.Options substructureOptions(ChemOptions chemistry) {
    Substructure.Options options = new Substructure.Options();
    options.chemistry = chemistry;
    return options;
  }

  private static MCS.Options mcsOptions(ChemOptions chemistry) {
    MCS.Options options = new MCS.Options();
    options.chemistry = chemistry;
    return options;
  }

  private static Map<Integer, Integer> firstMapping(SearchResult result) {
    return result.found() ? result.matches().get(0).atomMapping() : java.util.Collections.emptyMap();
  }
}
