/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
package com.bioinception.smsd.core;

import static org.junit.jupiter.api.Assertions.*;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.CsvSource;
import org.junit.jupiter.api.Test;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IStereoElement;
import org.openscience.cdk.interfaces.ITetrahedralChirality;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.stereo.TetrahedralChirality;

class CIPPriorityRegressionTest {
  @ParameterizedTest
  @CsvSource({
      "F[C@H](C(OC)C)C(O)CO, 1, S",
      "[2H]C([2H])([C@H](O)CCF)C, 3, S",
      "O[C@H](C[13CH2]C)CCCF, 1, R",
      "F[C@H](C(CO)C(CBr)(CBr)CBr)C(C(O)C)C(C)C, 1, R",
      "O[C@@H](c1ccccc1)C(C)(C)C, 1, R",
      "O[C@@H](C1=CC=CC=C1)C(C)(C)C, 1, R",
      "O=C1Nc2ccc(Cl)cc2[C@@](C#CC2CC2)(C(F)(F)F)O1, 10, S",
      "C[S@](=O)c1ccccc1, 1, S",
      "COc1ccc2[nH]c([S@@](=O)Cc3ncc(C)c(OC)c3C)nc2c1, 8, S",
      "C/C=C/[C@@H](F)/C=C\\C, 3, R",
      "C/C([C@H](F)/C(C)=C/C)=C/C, 2, R",
      "C/C=C(C)/[C@@H](F)/C(C)=C\\C, 4, R",
      "C/C(Cl)=C(C)/[C@@H](F)/C(C)=C(Cl)\\C, 5, S",
      "N[C@@H](C)C(=O)O, 1, S",
      "N[C@H](C)C(=O)O, 1, R"
  })
  void assignsPrioritiesByHierarchicalSequenceRules(String smiles, int centre, char expected) throws Exception {
    MolGraph molecule = new MolGraph(new SmilesParser(SilentChemObjectBuilder.getInstance()).parseSmiles(smiles));
    assertEquals(expected, CIPAssigner.assignRS(molecule).get(centre));
  }

  @Test
  void lonePairDescriptorUsesTheSpecifiedLigandOrder() throws Exception {
    IAtomContainer molecule = new SmilesParser(SilentChemObjectBuilder.getInstance())
        .parseSmiles("C[S@](=O)c1ccccc1");
    ITetrahedralChirality original = (ITetrahedralChirality) molecule.stereoElements().iterator().next();
    IAtom[] ligands = original.getLigands().clone();
    IAtom first = ligands[0];
    ligands[0] = ligands[1];
    ligands[1] = first;
    ITetrahedralChirality reordered = new TetrahedralChirality(original.getChiralAtom(), ligands,
        original.getStereo().invert());
    molecule.setStereoElements(java.util.Collections.<IStereoElement>singletonList(reordered));
    assertEquals('S', CIPAssigner.assignRS(new MolGraph(molecule)).get(1));
  }
}
