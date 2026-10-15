/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
package com.bioinception.smsd.core;

import com.bioinception.smsd.core.MCS;

import com.bioinception.smsd.TestSupport;

import static org.junit.jupiter.api.Assertions.*;

import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.CsvSource;
import org.junit.jupiter.params.provider.ValueSource;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;

class GraphIdentityRegressionTest {
  @ParameterizedTest
  @ValueSource(strings = {"C=CCC", "OCC=O", "CC[CH2-]", "[13CH3]CCC", "CC[CH]C", "[nH]1ccnc1",
      "CC(=O)O", "C=CC=C", "O=C1CCC(=O)C1", "[NH3+]CC(=O)[O-]", "F[C@H](Cl)CCO",
      "C[C@@H](N)C(=O)O", "N[C@@H]1CCC[C@H]1O", "C/C=C(/F)C(Cl)=C/C", "C/C=C\\C=C/C",
      "C[S@](=O)c1ccccc1", "CCO.CC=O", "CCO.CCO", "c1ccc2[nH]ccc2c1"})
  void atomAndBondOrderDoNotChangeIdentity(String text) throws Exception {
    SmilesParser parser = new SmilesParser(SilentChemObjectBuilder.getInstance());
    IAtomContainer original = parser.parseSmiles(text);
    MolGraph expected = new MolGraph(original);
    long hash = expected.getCanonicalHash();
    String smiles = expected.toCanonicalSmiles();
    for (int order = 0; order < 3; order++) {
      IAtomContainer reordered = original.clone();
      int atomCount = reordered.getAtomCount();
      IAtom[] atoms = new IAtom[atomCount];
      for (int atom = 0; atom < atomCount; atom++)
        atoms[atom] = reordered.getAtom(order == 0 ? atomCount - atom - 1 : (atom + order) % atomCount);
      reordered.setAtoms(atoms);
      IBond[] bonds = new IBond[reordered.getBondCount()];
      for (int bond = 0; bond < bonds.length; bond++) bonds[bond] = reordered.getBond(bonds.length - bond - 1);
      reordered.setBonds(bonds);
      MolGraph actual = new MolGraph(reordered);
      assertEquals(hash, actual.getCanonicalHash(), text);
      assertEquals(smiles, actual.toCanonicalSmiles(), text);
    }
  }

  @ParameterizedTest
  @CsvSource({"C,[CH3]", "C,[13CH4]", "[NH4+],N", "C=CCC,CC=CC", "C/C=C/C,C/C=C\\C",
      "F[C@H](Cl)Br,F[C@@H](Cl)Br"})
  void distinctChemistryHasDistinctIdentity(String first, String second) throws Exception {
    SmilesParser parser = new SmilesParser(SilentChemObjectBuilder.getInstance());
    MolGraph firstGraph = new MolGraph(parser.parseSmiles(first));
    MolGraph secondGraph = new MolGraph(parser.parseSmiles(second));
    assertNotEquals(firstGraph.getCanonicalHash(), secondGraph.getCanonicalHash());
  }

  @ParameterizedTest
  @ValueSource(strings = {"Cc1ccc(C(=O)O)cc1", "NC(=O)c1ccccc1"})
  void carbonylSubstituentsDoNotPreventSymmetryDeduplication(String text) throws Exception {
    SmilesParser parser = new SmilesParser(SilentChemObjectBuilder.getInstance());
    MolGraph query = new MolGraph(parser.parseSmiles("OC(=O)c1ccccc1"));
    MolGraph target = new MolGraph(parser.parseSmiles(text));
    MCS.Options options = new MCS.Options();
    options.timeoutMillis = 1000;
    java.util.List<java.util.Map<Integer,Integer>> embeddings=TestSupport.mcsMappings(query,target,new ChemOptions(),options,10);
    assertEquals(2,embeddings.size());
    java.util.Set<java.util.Map<Integer,Integer>> representatives=new java.util.HashSet<>();
    for(java.util.Map<Integer,Integer> embedding:embeddings)representatives.add(MCS.canonicalizeMapping(query,target,embedding));
    assertEquals(1,representatives.size());
  }
}
