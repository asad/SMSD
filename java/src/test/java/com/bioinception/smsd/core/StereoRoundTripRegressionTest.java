/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
package com.bioinception.smsd.core;

import static org.junit.jupiter.api.Assertions.*;

import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.ValueSource;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IStereoElement;
import org.openscience.cdk.isomorphism.Pattern;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;

class StereoRoundTripRegressionTest {
  @ParameterizedTest
  @ValueSource(strings = {
      "C[C@H](N)C(=O)O", "C[C@@H](N)C(=O)O", "CC[C@@H](O)C",
      "[C@H](F)(Cl)Br", "F[C@](Cl)(Br)I", "N[C@@H]1CCC[C@H]1O",
      "C/C=C/C", "C/C=C\\C", "C(=C/Cl)\\Cl", "Cl/C(C)=C/C",
      "C/C=C/C=C/C", "C/C=C\\C=C/C", "C/C=C/C=C\\C", "C/C=C\\C=C\\C",
      "C/C=C(/F)\\C(Cl)=C/C", "C/C=C(\\F)/C(Cl)=C/C"
  })
  void canonicalSmilesPreservesStereo(String smiles) throws Exception {
    SmilesParser parser = new SmilesParser(SilentChemObjectBuilder.getInstance());
    IAtomContainer original = parser.parseSmiles(smiles);
    assertRoundTrip(parser, original);

    IAtomContainer reordered = original.clone();
    IBond[] bonds = new IBond[reordered.getBondCount()];
    for (int bondIndex = 0; bondIndex < bonds.length; bondIndex++)
      bonds[bondIndex] = reordered.getBond(bonds.length - bondIndex - 1);
    reordered.setBonds(bonds);
    assertRoundTrip(parser, reordered);
  }

  private static void assertRoundTrip(SmilesParser parser, IAtomContainer original) throws Exception {
    String output = new MolGraph(original).toCanonicalSmiles();
    IAtomContainer reparsed = parser.parseSmiles(output);
    assertEquals(stereoCount(original), stereoCount(reparsed), output);
    assertTrue(Pattern.findIdentical(original).matches(reparsed), output);
  }

  private static int stereoCount(IAtomContainer molecule) {
    int count = 0;
    for (IStereoElement<?, ?> stereo : molecule.stereoElements()) count++;
    return count;
  }
}
