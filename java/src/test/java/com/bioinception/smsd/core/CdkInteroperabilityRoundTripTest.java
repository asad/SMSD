/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
package com.bioinception.smsd.core;

import static com.bioinception.smsd.TestSupport.rawMolecule;
import static org.junit.jupiter.api.Assertions.*;

import java.util.Map;
import java.util.Objects;
import java.util.TreeMap;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.CsvSource;
import org.junit.jupiter.params.provider.ValueSource;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IStereoElement;
import org.openscience.cdk.isomorphism.Pattern;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;

class CdkInteroperabilityRoundTripTest {
  @ParameterizedTest
  @ValueSource(strings = {"CC(=O)N", "[H]OC([H])([H])[H]", "c1cc[nH]c1",
      "[13CH3][NH3+]", "N[C@@H](C)C(=O)O", "F/C=C/Cl", "[Na+].CC(=O)[O-]"})
  void cdkSubstructureInputsReturnEquivalentCdkMolecules(String smiles) throws Exception {
    IAtomContainer query = rawMolecule(smiles), target = rawMolecule(smiles);
    IAtomContainer queryBefore = query.clone(), targetBefore = target.clone();
    Substructure.Options options = new Substructure.Options();
    options.chemistry = chemistry();
    SearchResult result = Substructure.find(query).withOptions(options).match(target);

    SearchResult.Match match = completeMatch(result);
    assertEquals(query.getAtomCount(), match.atomCount());
    assertEquals(query.getBondCount(), match.bondCount());
    assertExportsEquivalent(query, match);
    assertChemistryUnchanged(queryBefore, query);
    assertChemistryUnchanged(targetBefore, target);
  }

  @ParameterizedTest
  @CsvSource(value = {
      "MCIS|Connected|[13CH3][C@H](F)[NH3+]",
      "MCES|Connected|[13CH3][C@H](F)[NH3+]",
      "MCIS|Fragments|CCO.CN",
      "MCES|Fragments|CCO.CN"
  }, delimiter = '|')
  void cdkMCSExportsRetainWholeMatchedComponents(
      MCS.Objective objective, MCS.Connectivity connectivity, String smiles) throws Exception {
    IAtomContainer query = rawMolecule(smiles), target = rawMolecule(smiles);
    MCS.Options options = new MCS.Options();
    options.objective = objective;
    options.connectivity = connectivity;
    options.chemistry = chemistry();
    SearchResult.Match match = completeMatch(MCS.find(query).withOptions(options).match(target));

    assertEquals(query.getAtomCount(), match.atomCount());
    assertEquals(query.getBondCount(), match.bondCount());
    assertEquals(connectivity == MCS.Connectivity.Fragments ? 2 : 1, match.fragments().size());
    assertExportsEquivalent(query, match);
  }

  @Test
  void mappingUsesOriginalCdkIndicesAndPreparedInputsAndExportsAreIndependent() throws Exception {
    IAtomContainer query = rawMolecule("NCCO"), target = rawMolecule("CC(O)CN");
    IAtom[] reversed = {query.getAtom(3), query.getAtom(2), query.getAtom(1), query.getAtom(0)};
    query.setAtoms(reversed);
    IAtomContainer queryBefore = query.clone(), targetBefore = target.clone();
    IAtom[] targetAtoms = new IAtom[target.getAtomCount()];
    for (int atom = 0; atom < targetAtoms.length; atom++) targetAtoms[atom] = target.getAtom(atom);
    Substructure prepared = Substructure.find(query);
    SearchResult result = prepared.matchAll(target);
    SearchResult.Match match = completeMatch(result);
    assertEquals(1, result.matches().size());
    Map<Integer, Integer> expected = new TreeMap<>();
    expected.put(0, 2);
    expected.put(1, 1);
    expected.put(2, 3);
    expected.put(3, 4);
    assertEquals(expected, match.atomMapping());
    for (int atom = 0; atom < reversed.length; atom++) assertSame(reversed[atom], query.getAtom(atom));
    for (int atom = 0; atom < targetAtoms.length; atom++) assertSame(targetAtoms[atom], target.getAtom(atom));

    IAtomContainer exported = match.toMolecule(query);
    assertEquivalent(query, exported);
    for (int atom = 0; atom < query.getAtomCount(); atom++)
      assertNotSame(query.getAtom(atom), exported.getAtom(atom));
    exported.getAtom(0).setFormalCharge(-1);
    assertChemistryUnchanged(queryBefore, query);
    assertChemistryUnchanged(targetBefore, target);
    query.getAtom(3).setAtomicNumber(8);
    assertEquals(expected, completeMatch(prepared.match(target)).atomMapping());
  }

  @Test
  void cdkFragmentExportCapsCutValences() throws Exception {
    IAtomContainer query = rawMolecule("CCO"), target = rawMolecule("CCN");
    IAtomContainer queryBefore = query.clone();
    SearchResult.Match match = completeMatch(MCS.find(query).match(target));
    assertEquals(2, match.atomCount());
    assertEquals(1, match.bondCount());
    IAtomContainer ethane = rawMolecule("CC");
    assertEquivalent(ethane, match.toMolecule(query));
    assertEquivalent(ethane, rawMolecule(match.toSmiles(query)));
    assertChemistryUnchanged(queryBefore, query);
  }

  @ParameterizedTest
  @CsvSource(value = {"N[C@@H](C)C(=O)O|N[C@H](C)C(=O)O", "F/C=C/Cl|F/C=C\\Cl",
      "[13CH3]O|CO", "C[NH3+]|CN", "c1ccccc1|C1CCCCC1", "CCO|CCN"}, delimiter = '|')
  void cdkIdentityCheckRejectsChangedChemistry(String expected, String changed) throws Exception {
    IAtomContainer original = rawMolecule(expected), different = rawMolecule(changed);
    assertThrows(AssertionError.class, () -> assertEquivalent(original, different));
  }

  private static ChemOptions chemistry() {
    ChemOptions chemistry = new ChemOptions();
    chemistry.matchFormalCharge = true;
    chemistry.matchIsotope = true;
    chemistry.useChirality = true;
    chemistry.useBondStereo = true;
    chemistry.aromaticityMode = ChemOptions.AromaticityMode.STRICT;
    return chemistry;
  }

  private static SearchResult.Match completeMatch(SearchResult result) {
    assertEquals(SearchResult.Status.Completed, result.status());
    assertTrue(result.optimal());
    assertTrue(result.complete());
    assertTrue(result.found());
    return result.matches().get(0);
  }

  private static void assertExportsEquivalent(IAtomContainer query, SearchResult.Match match)
      throws Exception {
    IAtomContainer exported = match.toMolecule(query);
    assertEquivalent(query, exported);
    assertEquivalent(query, rawMolecule(
        new SmilesGenerator(SmiFlavor.Isomeric | SmiFlavor.UseAromaticSymbols).create(exported)));
    assertEquivalent(query, rawMolecule(match.toSmiles(query)));
  }

  private static void assertEquivalent(IAtomContainer expected, IAtomContainer actual) {
    assertEquals(expected.getAtomCount(), actual.getAtomCount());
    assertEquals(expected.getBondCount(), actual.getBondCount());
    assertEquals(stereoCount(expected), stereoCount(actual));
    int[] mapping = Pattern.findIdentical(expected).matchAll(actual).filter(indices -> {
      for (int atom = 0; atom < indices.length; atom++)
        if (!sameAtomChemistry(expected.getAtom(atom), actual.getAtom(indices[atom]))) return false;
      for (IBond bond : expected.bonds()) {
        IBond mapped = actual.getBond(actual.getAtom(indices[expected.indexOf(bond.getBegin())]),
            actual.getAtom(indices[expected.indexOf(bond.getEnd())]));
        if (mapped == null || bond.isAromatic() != mapped.isAromatic()) return false;
      }
      return true;
    }).first();
    assertEquals(expected.getAtomCount(), mapping.length, "CDK chemical identity including stereo");
  }

  private static boolean sameAtomChemistry(IAtom first, IAtom second) {
    return Objects.equals(first.getAtomicNumber(), second.getAtomicNumber())
        && Objects.equals(first.getFormalCharge(), second.getFormalCharge())
        && Objects.equals(first.getMassNumber(), second.getMassNumber())
        && Objects.equals(first.getImplicitHydrogenCount(), second.getImplicitHydrogenCount())
        && first.isAromatic() == second.isAromatic();
  }

  private static void assertChemistryUnchanged(IAtomContainer expected, IAtomContainer actual)
      throws Exception {
    assertEquals(expected.getAtomCount(), actual.getAtomCount());
    assertEquals(expected.getBondCount(), actual.getBondCount());
    for (int atom = 0; atom < expected.getAtomCount(); atom++)
      assertTrue(sameAtomChemistry(expected.getAtom(atom), actual.getAtom(atom)), "atom " + atom);
    for (int bond = 0; bond < expected.getBondCount(); bond++) {
      IBond before = expected.getBond(bond), after = actual.getBond(bond);
      assertEquals(expected.indexOf(before.getBegin()), actual.indexOf(after.getBegin()));
      assertEquals(expected.indexOf(before.getEnd()), actual.indexOf(after.getEnd()));
      assertEquals(before.getOrder(), after.getOrder());
      assertEquals(before.isAromatic(), after.isAromatic());
    }
    assertEquals(SmilesGenerator.isomeric().create(expected), SmilesGenerator.isomeric().create(actual));
  }

  private static int stereoCount(IAtomContainer molecule) {
    int count = 0;
    for (IStereoElement<?, ?> stereo : molecule.stereoElements()) count++;
    return count;
  }
}
