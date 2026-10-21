/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
package com.bioinception.smsd.core;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNull;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.CsvSource;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.silent.SilentChemObjectBuilder;

/** Valence and radical information supplies missing hydrogen counts. */
class MolHydrogenRegressionTest {
  private static IAtomContainer radicalCarbon(Integer hydrogens, Integer valence) {
    IChemObjectBuilder builder = SilentChemObjectBuilder.getInstance();
    IAtomContainer molecule = builder.newInstance(IAtomContainer.class);
    IAtom carbon = builder.newInstance(IAtom.class, "C");
    carbon.setImplicitHydrogenCount(hydrogens);
    carbon.setValency(valence);
    molecule.addAtom(carbon);
    molecule.addAtom(builder.newInstance(IAtom.class, "C"));
    molecule.addBond(0, 1, IBond.Order.SINGLE);
    molecule.addSingleElectron(0);
    molecule.addSingleElectron(0);
    return molecule;
  }

  @ParameterizedTest
  @CsvSource({"0,0", "1,0", "2,1"})
  void missingHydrogensRespectExplicitValence(int valence, int expected) {
    IAtomContainer molecule = radicalCarbon(null, valence);
    assertEquals(expected, new MolGraph(molecule).hydrogenCount(0));
    assertNull(molecule.getAtom(0).getImplicitHydrogenCount());
  }

  @Test
  void missingHydrogensUseRadicalElectronDemand() {
    assertEquals(1, new MolGraph(radicalCarbon(null, null)).hydrogenCount(0));
  }

  @ParameterizedTest
  @CsvSource({"0", "1", "2"})
  void suppliedHydrogensRemainAuthoritative(int hydrogens) {
    IAtomContainer molecule = radicalCarbon(hydrogens, 2);
    assertEquals(hydrogens, new MolGraph(molecule).hydrogenCount(0));
    assertEquals(hydrogens, molecule.getAtom(0).getImplicitHydrogenCount());
    assertEquals(hydrogens, new MolGraph.Builder().atomCount(1).atomicNumbers(new int[] {6})
        .neighbors(new int[][] {new int[0]}).hydrogenCounts(new int[] {hydrogens}).build().hydrogenCount(0));
  }
}
