/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
package com.bioinception.smsd.core;

import static com.bioinception.smsd.TestSupport.rawGraph;
import static org.junit.jupiter.api.Assertions.*;

import java.util.Arrays;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.TreeSet;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.CsvSource;
import org.junit.jupiter.params.provider.ValueSource;

class SmartsWriterTest {
  @ParameterizedTest
  @CsvSource(value = {
      "CC|CCC|STRICT|FLEXIBLE|true|false|false|false",
      "CC|C=C|STRICT|FLEXIBLE|true|false|false|false",
      "CC|c1ccccc1|STRICT|FLEXIBLE|true|false|false|false",
      "C=C|c1ccccc1|STRICT|FLEXIBLE|true|false|false|false",
      "C#C|c1ccccc1|STRICT|FLEXIBLE|true|false|false|false",
      "c1ccccc1|C1CCCCC1|STRICT|FLEXIBLE|true|false|false|false",
      "CC|C=C|LOOSE|STRICT|true|false|false|false",
      "CC|c1ccccc1|ANY|STRICT|true|false|false|false",
      "CC|NN|STRICT|FLEXIBLE|false|false|false|false",
      "N|C[NH3+]|STRICT|FLEXIBLE|true|true|false|false",
      "N|CN|STRICT|FLEXIBLE|true|true|false|false",
      "[NH3+]|CN|STRICT|FLEXIBLE|true|true|false|false",
      "[NH3+]|C[NH3+]|STRICT|FLEXIBLE|true|true|false|false",
      "[13CH4]|C|STRICT|FLEXIBLE|true|false|true|false",
      "[13CH4]|[13CH4]|STRICT|FLEXIBLE|true|false|true|false",
      "[13CH4]|[14CH4]|STRICT|FLEXIBLE|true|false|true|false",
      "C|[13CH4]|STRICT|FLEXIBLE|true|false|true|false",
      "[13CH4]|[13NH3]|STRICT|FLEXIBLE|false|false|true|false",
      "C1CC1|C1CC1|STRICT|FLEXIBLE|true|false|false|true",
      "CC|C1CC1|STRICT|FLEXIBLE|true|false|false|true",
      "C1CC1|CCC|ANY|FLEXIBLE|true|false|false|true",
      "C1CC2CCC1C2|C1CC2CCC1C2|STRICT|STRICT|true|false|false|true"
  }, delimiter = '|')
  void portableQueriesPreserveSupportedMoleculePredicates(String querySmiles, String targetSmiles,
      ChemOptions.BondOrderMode order, ChemOptions.AromaticityMode aromaticity,
      boolean elements, boolean charges, boolean isotopes, boolean rings) throws Exception {
    ChemOptions chemistry = new ChemOptions();
    chemistry.matchBondOrder = order;
    chemistry.aromaticityMode = aromaticity;
    chemistry.matchAtomType = elements;
    chemistry.matchFormalCharge = charges;
    chemistry.matchIsotope = isotopes;
    chemistry.ringMatchesRingOnly = rings;
    MolGraph query = rawGraph(querySmiles), target = rawGraph(targetSmiles);
    Substructure.Options options = new Substructure.Options();
    options.chemistry = chemistry;
    SearchResult expected = Substructure.find(query).withOptions(options).matchAll(target);
    String pattern = SmartsWriter.toSmarts(query, chemistry);
    SearchResult actual = Substructure.smarts(pattern).matchAll(target);
    assertEquals(SearchResult.Status.Completed, expected.status());
    assertEquals(SearchResult.Status.Completed, actual.status(), pattern);
    assertEquals(expected.matches().size(), actual.matches().size(), pattern);
  }

  @Test void anMCSCanSearchAnotherLibrary() throws Exception {
    MolGraph query = rawGraph("CC(=O)Nc1ccccc1"), target = rawGraph("CC(=O)Oc1ccccc1");
    SearchResult common = MCS.find(query).match(target);
    assertTrue(common.found());
    String pattern = common.matches().get(0).toSmarts(query, new ChemOptions());
    Substructure libraryQuery = Substructure.smarts(pattern);
    assertTrue(libraryQuery.match(query).found());
    assertTrue(libraryQuery.match(target).found());
    assertFalse(libraryQuery.match(rawGraph("CCO")).found());
  }

  @Test void exportsOnlyTheSelectedMCESEdges() throws Exception {
    MolGraph query = rawGraph("C1CC1"), target = rawGraph("CCC");
    MCS.Options options = new MCS.Options();
    options.objective = MCS.Objective.MCES;
    SearchResult common = MCS.find(query).withOptions(options).match(target);
    SearchResult.Match match = common.matches().get(0);
    assertEquals(3, match.atomCount());
    assertEquals(2, match.bondCount());
    String pattern = match.toSmarts(query, options.chemistry);
    SearchResult roundTrip = Substructure.smarts(pattern).match(target);
    assertTrue(roundTrip.found(), pattern);
    assertEquals(2, roundTrip.matches().get(0).bondCount());
  }

  @Test void keepsDisconnectedSelectedAtomsAndOriginalRingContext() throws Exception {
    MolGraph source = rawGraph("C1CCCCC1.O");
    Map<Integer, Integer> atoms = new LinkedHashMap<>();
    atoms.put(0, 0);
    atoms.put(1, 1);
    atoms.put(6, 6);
    SearchResult.Bond bond = new SearchResult.Bond(0, 1);
    SearchResult.Match match = new SearchResult.Match(atoms, Collections.singletonMap(bond, bond),
        Arrays.asList(new TreeSet<>(Arrays.asList(0, 1)), Collections.singleton(6)), 3.0);
    ChemOptions chemistry = new ChemOptions();
    chemistry.ringMatchesRingOnly = true;
    String pattern = match.toSmarts(source, chemistry);
    assertTrue(pattern.contains("."));
    assertTrue(Substructure.smarts(pattern).match(source).found(), pattern);
    assertFalse(Substructure.smarts(pattern).match(rawGraph("CC.O")).found(), pattern);
  }

  @ParameterizedTest
  @ValueSource(strings = {"chirality", "bondStereo", "tautomers", "completeRings", "ringFusion"})
  void rejectsPredicatesThatPortableSmartsCannotExpress(String option) throws Exception {
    ChemOptions chemistry = new ChemOptions();
    switch (option) {
      case "chirality": chemistry.useChirality = true; break;
      case "bondStereo": chemistry.useBondStereo = true; break;
      case "tautomers": chemistry.tautomerAware = true; break;
      case "completeRings": chemistry.completeRingsOnly = true; break;
      case "ringFusion": chemistry.ringFusionMode = ChemOptions.RingFusionMode.STRICT; break;
      default: fail(option);
    }
    MolGraph source = rawGraph("CCO");
    assertThrows(IllegalArgumentException.class, () -> SmartsWriter.toSmarts(source, chemistry));
  }

  @Test void rejectsTheConditionalAromaticEndpointResonanceRule() throws Exception {
    MolGraph source = rawGraph("c1ccccc1-c1ccccc1");
    assertThrows(IllegalArgumentException.class, () -> SmartsWriter.toSmarts(source));
    ChemOptions strict = ChemOptions.profile("strict");
    assertTrue(Substructure.smarts(SmartsWriter.toSmarts(source, strict)).match(source).found());
  }

  @Test void rejectsStrictKekuleFlagsChangedBySmartsTargetPreparation() throws Exception {
    MolGraph source = rawGraph("C1=CC=CC=C1");
    ChemOptions strict = ChemOptions.profile("strict");
    Substructure.Options options = new Substructure.Options();
    options.chemistry = strict;
    assertTrue(Substructure.find(source).withOptions(options).match(source).found());
    assertFalse(Substructure.smarts("[#6;A]").match(source).found());
    IllegalArgumentException exception = assertThrows(IllegalArgumentException.class,
        () -> SmartsWriter.toSmarts(source, strict));
    assertTrue(exception.getMessage().contains("SMARTS target preparation"));
    MolGraph aromatic = rawGraph("c1ccccc1");
    assertTrue(Substructure.smarts(SmartsWriter.toSmarts(aromatic, strict)).match(aromatic).found());
  }

  @Test void writesAnEmptyGraphAsEmptyText() throws Exception {
    assertEquals("", SmartsWriter.toSmarts(rawGraph("")));
  }
}
