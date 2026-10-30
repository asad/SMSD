/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
package com.bioinception.smsd.core;

import static org.junit.jupiter.api.Assertions.*;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.CsvSource;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;

class MCSRingViabilityTest {
  private static MolGraph molecule(String smiles) throws Exception {
    return new MolGraph(new SmilesParser(SilentChemObjectBuilder.getInstance()).parseSmiles(smiles));
  }

  private static MCS.Options ringOptions() {
    MCS.Options options = new MCS.Options();
    options.chemistry.completeRingsOnly = true;
    options.timeoutMillis = 1000;
    options.nodeLimit = 100000;
    return options;
  }

  @ParameterizedTest
  @CsvSource({
      "CC1CCCCC1, CC1CCCC1, 1, 0",
      "CCC1CCCC1, CCC1CCCCC1, 2, 1",
      "CCOC1CCCCC1, CCNC1CCCCC1, 6, 6",
      "CC1CCC2CCCCC2C1, CC1CCC2CCCC2C1, 7, 7",
      "CC1CCC2(CC1)CCCC2, CC1CCC2(CC1)CCCCC2, 7, 7"
  })
  void wholeRingsAndTheirBranchesKeepTheMaximum(String first, String second, int atoms, int bonds) throws Exception {
    for (boolean reverse : new boolean[] {false, true}) {
      MolGraph query = molecule(reverse ? second : first), target = molecule(reverse ? first : second);
      MCS.Options options = ringOptions();
      SearchResult result = MCS.find(query).withOptions(options).match(target);
      assertEquals(SearchResult.Status.Completed, result.status());
      assertTrue(result.optimal() && result.complete());
      assertTrue(result.found());
      assertEquals(atoms, result.matches().get(0).atomCount());
      assertEquals(bonds, result.matches().get(0).bondCount());
      assertTrue(MCS.validateMapping(query, target, result.matches().get(0), options).isEmpty());
    }
  }

  @Test
  void ringWeightsAndEnumerationKeepTheirContracts() throws Exception {
    MolGraph query = molecule("CCC1CCCC1"), target = molecule("CCC1CCCCC1");
    MCS.Options options = ringOptions();
    options.atomWeights = new double[] {0.25, 0.5, 100, 100, 100, 100, 100};
    SearchResult weighted = MCS.find(query).withOptions(options).match(target);
    assertTrue(weighted.optimal() && weighted.complete());
    assertEquals(0.75, weighted.score());
    assertEquals(2, weighted.matches().get(0).atomCount());
    MolGraph molecule = molecule("CC1CCCC1");
    options = ringOptions();
    SearchResult all = MCS.find(molecule).withOptions(options).matchAll(molecule, 100);
    assertTrue(all.optimal() && all.complete());
    assertEquals(2, all.matches().size());
    for (SearchResult.Match match : all.matches())
      assertTrue(MCS.validateMapping(molecule, molecule, match, options).isEmpty());
  }

  @Test
  void ringBudgetsKeepTheirStatus() throws Exception {
    MolGraph query = molecule("CC1CCCCC1"), target = molecule("CC1CCCC1");
    MCS.Options options = ringOptions();
    options.nodeLimit = 1;
    SearchResult capped = MCS.find(query).withOptions(options).match(target);
    assertEquals(SearchResult.Status.LimitReached, capped.status());
    assertFalse(capped.optimal() || capped.complete());
    assertTrue(capped.nodesVisited() <= 1);
    options.timeoutMillis = 0;
    SearchResult timed = MCS.find(query).withOptions(options).match(target);
    assertEquals(SearchResult.Status.TimedOut, timed.status());
    assertFalse(timed.optimal() || timed.complete());
    options.timeoutMillis = 1000;
    options.cancellation = () -> true;
    SearchResult cancelled = MCS.find(query).withOptions(options).match(target);
    assertEquals(SearchResult.Status.Cancelled, cancelled.status());
    assertFalse(cancelled.optimal() || cancelled.complete());
  }
}
