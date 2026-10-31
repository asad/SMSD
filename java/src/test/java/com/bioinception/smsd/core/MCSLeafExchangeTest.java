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

class MCSLeafExchangeTest {
  private static MolGraph molecule(String smiles) throws Exception {
    SmilesParser parser = new SmilesParser(SilentChemObjectBuilder.getInstance());
    return new MolGraph(parser.parseSmiles(smiles));
  }

  private static MCS.Options strictOptions() {
    MCS.Options options = new MCS.Options();
    options.chemistry = ChemOptions.profile("strict");
    options.timeoutMillis = 1000;
    options.nodeLimit = 100000;
    return options;
  }

  @ParameterizedTest
  @CsvSource({
      "CCC(C)CC, CCCCCC, 5, 4",
      "CC(C)(C)OCC, CC(C)(C)NCC, 4, 3",
      "CC(C)NCC(C)O, CC(C)NCC(C)N, 7, 6",
      "CC(C)(C)CC(C)(C)O, CC(C)(C)CCC(C)(C)N, 7, 6",
      "OC(C)(C)CC(C)(C)C, NC(C)(C)CCC(C)(C)C, 7, 6"
  })
  void pendantBranchesKeepTheMaximumInEitherArgumentOrder(
      String first, String second, int atoms, int bonds) throws Exception {
    for (boolean reverse : new boolean[] {false, true}) {
      MolGraph query = molecule(reverse ? second : first);
      MolGraph target = molecule(reverse ? first : second);
      MCS.Options options = strictOptions();
      SearchResult result = MCS.find(query).withOptions(options).match(target);
      assertEquals(SearchResult.Status.Completed, result.status());
      assertTrue(result.optimal() && result.complete());
      assertEquals(atoms, result.matches().get(0).atomCount());
      assertEquals(bonds, result.matches().get(0).bondCount());
      assertTrue(MCS.validateMapping(query, target, result.matches().get(0), options).isEmpty());
    }
  }

  @Test
  void isotopeWildcardsKeepAlternativePairs() throws Exception {
    MCS.Options options = strictOptions();
    options.chemistry.matchIsotope = true;
    MolGraph query = molecule("CN[13CH2]O"), target = molecule("CN[12CH2]S");
    SearchResult result = MCS.find(query).withOptions(options).match(target);
    assertTrue(result.optimal() && result.complete());
    assertEquals(3, result.matches().get(0).atomCount());
    assertEquals(2, result.matches().get(0).bondCount());
    assertTrue(MCS.validateMapping(query, target, result.matches().get(0), options).isEmpty());
  }

  @Test
  void weightsAndEnumerationKeepTheirContracts() throws Exception {
    MolGraph graph = molecule("CC(C)C");
    MCS.Options options = strictOptions();
    SearchResult all = MCS.find(graph).withOptions(options).matchAll(graph, 100);
    assertTrue(all.optimal() && all.complete());
    assertEquals(6, all.matches().size());
    options.atomWeights = new double[] {-5, 1, 1, 1};
    SearchResult weighted = MCS.find(graph).withOptions(options).match(graph);
    assertTrue(weighted.optimal() && weighted.complete());
    assertEquals(3, weighted.score());
    assertEquals(3, weighted.matches().get(0).atomCount());
    assertFalse(weighted.matches().get(0).atomMapping().containsKey(0));
  }

  @Test
  void stoppedSearchesKeepTheirStatus() throws Exception {
    MolGraph query = molecule("CC(C)(C)CC(C)(C)O");
    MolGraph target = molecule("CC(C)(C)CCC(C)(C)N");
    MCS.Options options = strictOptions();
    options.nodeLimit = 1;
    SearchResult capped = MCS.find(query).withOptions(options).match(target);
    assertEquals(SearchResult.Status.LimitReached, capped.status());
    assertFalse(capped.optimal() || capped.complete());
    assertTrue(capped.nodesVisited() <= options.nodeLimit);
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
