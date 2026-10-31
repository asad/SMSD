/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
package com.bioinception.smsd.core;

import static org.junit.jupiter.api.Assertions.*;

import java.util.Collections;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.CsvSource;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;

class MCSTargetTwinTest {
  private static MolGraph molecule(String smiles) throws Exception {
    return new MolGraph(new SmilesParser(SilentChemObjectBuilder.getInstance()).parseSmiles(smiles));
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
      "CCCC, CC1CC1, 3, 2, false",
      "CCCCC, CCC1CC1, 4, 3, false",
      "CCCC, C1CCC1, 3, 2, false",
      "CCOC, C1COC1, 3, 2, false",
      "CCCO, CCOC, 3, 2, false",
      "CC=CC, C1=CCC1, 3, 2, false",
      "CCCO, C1CCN1, 3, 2, false",
      "C[CH-]C, [CH-]1CCC1, 3, 2, true"
  })
  void singleMatchesBelongToTheCompleteOptimumSet(
      String first, String second, int atoms, int bonds, boolean charge) throws Exception {
    for (boolean reverse : new boolean[] {false, true}) {
      MolGraph query = molecule(reverse ? second : first), target = molecule(reverse ? first : second);
      MCS.Options options = strictOptions();
      options.chemistry.matchFormalCharge = charge;
      MCS prepared = MCS.find(query).withOptions(options);
      SearchResult all = prepared.matchAll(target, 10000), single = prepared.match(target);
      assertEquals(SearchResult.Status.Completed, all.status());
      assertTrue(all.optimal() && all.complete());
      assertEquals(SearchResult.Status.Completed, single.status());
      assertTrue(single.optimal() && single.complete());
      assertEquals(1, single.matches().size());
      assertEquals(all.score(), single.score());
      assertEquals(atoms, single.score());
      assertEquals(atoms, single.matches().get(0).atomCount());
      assertEquals(bonds, single.matches().get(0).bondCount());
      Set<Map<Integer, Integer>> mappings = new HashSet<>();
      for (SearchResult.Match match : all.matches()) {
        assertTrue(MCS.validateMapping(query, target, match, options).isEmpty());
        assertTrue(mappings.add(match.atomMapping()));
      }
      assertTrue(mappings.contains(single.matches().get(0).atomMapping()));
      assertTrue(MCS.validateMapping(query, target, single.matches().get(0), options).isEmpty());
    }
  }

  @Test
  void weightsEnumerationAndPinnedAtomsKeepTheirContracts() throws Exception {
    MolGraph graph = molecule("C1CCC1");
    MCS.Options options = strictOptions();
    SearchResult all = MCS.find(graph).withOptions(options).matchAll(graph, 100);
    assertTrue(all.optimal() && all.complete());
    assertEquals(8, all.matches().size());
    options.atomWeights = new double[] {-5, 1, 1, 1};
    SearchResult weighted = MCS.find(graph).withOptions(options).match(graph);
    assertTrue(weighted.optimal() && weighted.complete());
    assertEquals(3, weighted.score());
    assertEquals(3, weighted.matches().get(0).atomCount());
    assertFalse(weighted.matches().get(0).atomMapping().containsKey(0));
    options.atomWeights = null;
    assertFalse(MCS.isMappingMaximal(graph, graph, Collections.singletonMap(0, 2), options));
    options.excludedTargetAtoms.add(0);
    SearchResult excluded = MCS.find(graph).withOptions(options).match(graph);
    assertTrue(excluded.optimal() && excluded.complete());
    assertEquals(3, excluded.score());
    assertFalse(excluded.matches().get(0).atomMapping().containsValue(0));
    assertTrue(MCS.validateMapping(graph, graph, excluded.matches().get(0), options).isEmpty());
  }

  @Test
  void stoppedSearchesKeepTheirStatus() throws Exception {
    MolGraph query = molecule("CCCCC"), target = molecule("C1CCC1");
    MCS.Options options = strictOptions();
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
