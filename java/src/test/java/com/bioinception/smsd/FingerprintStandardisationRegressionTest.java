/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
package com.bioinception.smsd;

import com.bioinception.smsd.TestSupport;

import static org.junit.jupiter.api.Assertions.*;

import com.bioinception.smsd.core.ChemOptions;
import com.bioinception.smsd.core.FingerprintEngine;
import com.bioinception.smsd.core.MolGraph;
import com.bioinception.smsd.TestSupport.MoleculePair;
import com.bioinception.smsd.core.MCS;
import java.io.StringReader;
import java.io.StringWriter;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.EnumSource;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.io.MDLV2000Reader;
import org.openscience.cdk.io.MDLV2000Writer;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;

class FingerprintStandardisationRegressionTest {
  private static IAtomContainer phenol(boolean aromatic) throws Exception {
    return new SmilesParser(SilentChemObjectBuilder.getInstance())
        .parseSmiles(aromatic ? "Oc1ccccc1" : "OC1=CC=CC=C1");
  }

  @ParameterizedTest
  @EnumSource(FingerprintEngine.Mode.class)
  void circularBinaryDispatchIgnoresKekuleRepresentation(FingerprintEngine.Mode mode) throws Exception {
    assertArrayEquals(FingerprintEngine.circularFingerprint(phenol(true), 2, 2048, mode),
        FingerprintEngine.circularFingerprint(phenol(false), 2, 2048, mode));
  }

  @ParameterizedTest
  @EnumSource(FingerprintEngine.Mode.class)
  void circularCountDispatchIgnoresKekuleRepresentation(FingerprintEngine.Mode mode) throws Exception {
    assertArrayEquals(FingerprintEngine.circularFingerprintCounts(phenol(true), 2, 2048, mode),
        FingerprintEngine.circularFingerprintCounts(phenol(false), 2, 2048, mode));
  }

  @ParameterizedTest
  @EnumSource(FingerprintEngine.Mode.class)
  void circularBinaryMethodsIgnoreKekuleRepresentation(FingerprintEngine.Mode mode) throws Exception {
    IAtomContainer aromatic = phenol(true), kekule = phenol(false);
    assertArrayEquals(mode == FingerprintEngine.Mode.ECFP ? FingerprintEngine.ecfp(aromatic, 2, 2048)
            : FingerprintEngine.fcfp(aromatic, 2, 2048),
        mode == FingerprintEngine.Mode.ECFP ? FingerprintEngine.ecfp(kekule, 2, 2048)
            : FingerprintEngine.fcfp(kekule, 2, 2048));
  }

  @ParameterizedTest
  @EnumSource(FingerprintEngine.Mode.class)
  void circularCountMethodsIgnoreKekuleRepresentation(FingerprintEngine.Mode mode) throws Exception {
    IAtomContainer aromatic = phenol(true), kekule = phenol(false);
    assertArrayEquals(mode == FingerprintEngine.Mode.ECFP ? FingerprintEngine.ecfpCounts(aromatic, 2, 2048)
            : FingerprintEngine.fcfpCounts(aromatic, 2, 2048),
        mode == FingerprintEngine.Mode.ECFP ? FingerprintEngine.ecfpCounts(kekule, 2, 2048)
            : FingerprintEngine.fcfpCounts(kekule, 2, 2048));
  }

  @ParameterizedTest
  @EnumSource(FingerprintEngine.Mode.class)
  void batchConversionIgnoresKekuleRepresentation(FingerprintEngine.Mode mode) throws Exception {
    MolGraph aromatic = FingerprintEngine.toMolGraph(phenol(true));
    MolGraph kekule = FingerprintEngine.toMolGraph(phenol(false));
    assertArrayEquals(FingerprintEngine.circularFingerprint(aromatic, 2, 2048, mode),
        FingerprintEngine.circularFingerprint(kekule, 2, 2048, mode));
    assertArrayEquals(FingerprintEngine.circularFingerprintCounts(aromatic, 2, 2048, mode),
        FingerprintEngine.circularFingerprintCounts(kekule, 2, 2048, mode));
  }

  @Test
  void topologicalTorsionIgnoresKekuleRepresentation() throws Exception {
    assertArrayEquals(FingerprintEngine.topologicalTorsion(phenol(true), 2048),
        FingerprintEngine.topologicalTorsion(phenol(false), 2048));
  }

  @Test
  void pathFingerprintForMCSIgnoresKekuleRepresentation() throws Exception {
    assertArrayEquals(FingerprintEngine.mcsFingerprint(phenol(true), new ChemOptions(), 7, 2048),
        FingerprintEngine.mcsFingerprint(phenol(false), new ChemOptions(), 7, 2048));
    assertArrayEquals(FingerprintEngine.mcsFingerprint(phenol(true), new ChemOptions(), 7, 2048),
        FingerprintEngine.mcsFingerprint(phenol(false), new ChemOptions(), 7, 2048));
  }

  @Test
  void smilesConvenienceIgnoresKekuleRepresentation() throws Exception {
    assertArrayEquals(TestSupport.fingerprintFromSmiles("Oc1ccccc1", 2, 2048),
        TestSupport.fingerprintFromSmiles("OC1=CC=CC=C1", 2, 2048));
  }

  @Test
  void ecfpSurvivesAMolfileRoundTrip() throws Exception {
    IAtomContainer molecule = phenol(true);
    StringWriter text = new StringWriter();
    try (MDLV2000Writer writer = new MDLV2000Writer(text)) {
      writer.write(molecule);
    }
    IAtomContainer restored;
    try (MDLV2000Reader reader = new MDLV2000Reader(new StringReader(text.toString()))) {
      restored = reader.read(SilentChemObjectBuilder.getInstance().newAtomContainer());
    }
    assertArrayEquals(FingerprintEngine.ecfp(molecule, 2, 2048),
        FingerprintEngine.ecfp(restored, 2, 2048));
  }

  @Test
  void fingerprintingPreservesTheInputAromaticFlags() throws Exception {
    IAtomContainer molecule = phenol(false);
    FingerprintEngine.ecfp(molecule, 2, 2048);
    for (int atom = 0; atom < molecule.getAtomCount(); atom++) assertFalse(molecule.getAtom(atom).isAromatic());
    for (int bond = 0; bond < molecule.getBondCount(); bond++) assertFalse(molecule.getBond(bond).isAromatic());
  }
}
