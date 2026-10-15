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

import java.util.Arrays;
import java.util.Map;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.CsvSource;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;

class SearchChemistryRegressionTest {
  private static IAtomContainer molecule(String smiles) throws Exception {
    return new SmilesParser(SilentChemObjectBuilder.getInstance()).parseSmiles(smiles);
  }

  private static ChemOptions strictFusion() {
    ChemOptions chemistry = new ChemOptions();
    chemistry.ringFusionMode = ChemOptions.RingFusionMode.STRICT;
    return chemistry;
  }

  @Test
  void strictFusionBoundAllowsRingToChainMatches() throws Exception {
    MolGraph query = new MolGraph(molecule("C1CCCCC1"));
    MolGraph target = new MolGraph(molecule("CCCCCC"));
    ChemOptions chemistry = strictFusion();
    assertTrue(MCS.similarityUpperBound(query, target, chemistry) >= 5.0 / 7.0);
    chemistry.ringMatchesRingOnly = true;
    assertEquals(0.0, MCS.similarityUpperBound(query, target, chemistry));
  }

  @Test
  void strictFusionScreenRetainsAValidRingToChainPair() throws Exception {
    IAtomContainer query = molecule("C1CCCCC1");
    IAtomContainer target = molecule("CCCCCC");
    Map<Integer, Map<Integer, Integer>> mappings = TestSupport.pairMappings(
        Arrays.asList(query, target), strictFusion(), new MCS.Options(), 0.5, 1);
    assertEquals(1, mappings.size());
    assertEquals(5, mappings.get(1).size());
  }

  @ParameterizedTest
  @CsvSource({
      "OC1CCCCC1O, NCCCCCCN, 5",
      "C1CCC(CC1)CCC1CCCCC1, CCCCCCCCCCCC, 12"
  })
  void strictFusionSearchKeepsValidRingToChainMappings(
      String querySmiles, String targetSmiles, int expectedSize) throws Exception {
    MolGraph query = new MolGraph(molecule(querySmiles));
    MolGraph target = new MolGraph(molecule(targetSmiles));
    ChemOptions chemistry = strictFusion();
    MCS.Options options = new MCS.Options();
    options.timeoutMillis = 2_000;
    Map<Integer, Integer> mapping = TestSupport.mcsMapping(query, target, chemistry, options);
    assertEquals(expectedSize, mapping.size());
    assertTrue(MCS.validateMapping(query, target, mapping, chemistry).isEmpty());
  }

  @ParameterizedTest
  @CsvSource({
      "C1CCCCC1, C1CCC2CCCCC2C1, false, 6",
      "CCC, C1CC1, true, 2",
      "CCC, C1CC1, false, 3"
  })
  void bondGrowthPreservesMappedBondsAndTheRequestedInducedPolicy(
      String querySmiles, String targetSmiles, boolean induced, int expectedSize) throws Exception {
    MolGraph query = new MolGraph(molecule(querySmiles));
    MolGraph target = new MolGraph(molecule(targetSmiles));
    query.ensureCanonical();
    target.ensureCanonical();
    ChemOptions chemistry = new ChemOptions();
    MCS.Options options=new MCS.Options();options.chemistry=chemistry;
    options.objective=induced?MCS.Objective.MCIS:MCS.Objective.MCES;
    SearchResult result=MCS.find(query).withOptions(options).match(target);
    Map<Integer,Integer> mapping=result.matches().get(0).atomMapping();
    assertTrue(result.optimal());
    assertEquals(expectedSize, mapping.size());
    if (induced) {
      for (int queryAtom : mapping.keySet())
        for (int otherQueryAtom : mapping.keySet())
          assertEquals(query.hasBond(queryAtom, otherQueryAtom),
              target.hasBond(mapping.get(queryAtom), mapping.get(otherQueryAtom)));
    }
  }
}
