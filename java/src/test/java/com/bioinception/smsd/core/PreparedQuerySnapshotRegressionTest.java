/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
package com.bioinception.smsd.core;

import java.util.Arrays;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.ValueSource;
import org.openscience.cdk.interfaces.IAtomContainer;
import static com.bioinception.smsd.TestSupport.rawMolecule;
import static org.junit.jupiter.api.Assertions.*;

class PreparedQuerySnapshotRegressionTest {
  private static MolGraph carbon() {
    return new MolGraph.Builder().atomCount(1).atomicNumbers(new int[]{6}).neighbors(new int[][]{{}}).build();
  }

  @Test void preparedGraphQueriesRetainTheOriginalElement() {
    MolGraph query = carbon(), target = carbon();
    MCS maximum = MCS.find(query);
    Substructure substructure = Substructure.find(query);
    query.atomicNum[0] = 7;
    assertTrue(maximum.match(target).found());
    assertTrue(substructure.match(target).found());
    assertFalse(MCS.find(query).match(target).found());
    assertFalse(Substructure.find(query).match(target).found());
  }

  @Test void preparedGraphQueriesRetainChargeAndIsotope() {
    MolGraph query = new MolGraph.Builder().atomCount(1).atomicNumbers(new int[]{6})
        .formalCharges(new int[]{1}).massNumbers(new int[]{13}).neighbors(new int[][]{{}}).build();
    MolGraph target = new MolGraph.Builder().atomCount(1).atomicNumbers(new int[]{6})
        .formalCharges(new int[]{1}).massNumbers(new int[]{13}).neighbors(new int[][]{{}}).build();
    MCS.Options maximumOptions = new MCS.Options();
    maximumOptions.chemistry.matchFormalCharge = true;
    maximumOptions.chemistry.matchIsotope = true;
    Substructure.Options substructureOptions = new Substructure.Options();
    substructureOptions.chemistry = maximumOptions.chemistry;
    MCS maximum = MCS.find(query).withOptions(maximumOptions);
    Substructure substructure = Substructure.find(query).withOptions(substructureOptions);
    query.formalCharge[0] = 0;
    query.massNumber[0] = 12;
    assertTrue(maximum.match(target).found());
    assertTrue(substructure.match(target).found());
    assertFalse(MCS.find(query).withOptions(maximumOptions).match(target).found());
    assertFalse(Substructure.find(query).withOptions(substructureOptions).match(target).found());
  }

  @ParameterizedTest @ValueSource(strings = {"[C@H](F)(Cl)Br", "Cl/C=C(/F)Br"})
  void preparedQueriesRetainStereoAndReferenceLigands(String smiles) throws Exception {
    MolGraph query = new MolGraph(rawMolecule(smiles)), target = new MolGraph(rawMolecule(smiles));
    MCS.Options maximumOptions = new MCS.Options();
    maximumOptions.chemistry.useChirality = true;
    maximumOptions.chemistry.useBondStereo = true;
    Substructure.Options substructureOptions = new Substructure.Options();
    substructureOptions.chemistry = maximumOptions.chemistry;
    MCS maximum = MCS.find(query).withOptions(maximumOptions);
    Substructure substructure = Substructure.find(query).withOptions(substructureOptions);
    if (smiles.startsWith("[")) {
      query.tetraChirality[0] = 3 - query.tetraChirality[0];
    } else {
      for (int atom = 0; atom < query.n; atom++)
        for (int neighbour : query.neighbors[atom])
          if (atom < neighbour && query.dbStereo(atom, neighbour) != 0) {
            int configuration = 3 - query.dbStereo(atom, neighbour);
            query.dbStereoConf[atom][neighbour] = query.dbStereoConf[neighbour][atom] = configuration;
          }
    }
    assertEquals(query.n, maximum.match(target).matches().get(0).atomCount());
    assertTrue(substructure.match(target).found());
    assertFalse(Substructure.find(query).withOptions(substructureOptions).match(target).found());
  }

  @Test void sparseBondStereoIsPartOfTheSnapshot() throws Exception {
    char[] prefix = new char[201];
    Arrays.fill(prefix, 'C');
    String smiles = new String(prefix) + "/C=C/Cl";
    MolGraph query = new MolGraph(rawMolecule(smiles)), target = new MolGraph(rawMolecule(smiles));
    assertTrue(query.n > 200);
    MCS.Options maximumOptions = new MCS.Options();
    maximumOptions.chemistry.useBondStereo = true;
    Substructure.Options substructureOptions = new Substructure.Options();
    substructureOptions.chemistry = maximumOptions.chemistry;
    MCS maximum = MCS.find(query).withOptions(maximumOptions);
    Substructure substructure = Substructure.find(query).withOptions(substructureOptions);
    for (int[] row : query.neighbourDbStereo)
      for (int k = 0; k < row.length; k++)
        if (row[k] != 0) row[k] = 3 - row[k];
    assertEquals(query.n, maximum.match(target).matches().get(0).atomCount());
    assertTrue(substructure.match(target).found());
    assertFalse(Substructure.find(query).withOptions(substructureOptions).match(target).found());
  }

  @Test void containerInputsAlreadyHaveIndependentPreparedState() throws Exception {
    IAtomContainer query = rawMolecule("C");
    MolGraph target = carbon();
    MCS maximum = MCS.find(query);
    Substructure substructure = Substructure.find(query);
    query.getAtom(0).setAtomicNumber(7);
    assertTrue(maximum.match(target).found());
    assertTrue(substructure.match(target).found());
  }
}
