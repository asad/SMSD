/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
package com.bioinception.smsd.core;

import static com.bioinception.smsd.TestSupport.rawMolecule;
import static org.junit.jupiter.api.Assertions.*;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.CsvSource;

class SubstructureCompleteRingsTest {
  @ParameterizedTest
  @CsvSource(value = {
      "false|C|C1CCCCC1|false", "true|C|C1CCCCC1|false",
      "false|C1CCCCC1|C1CCCCC1|true", "true|C1CCCCC1|C1CCCCC1|true",
      "false|C|CCC|true", "true|C|CCC|true",
      "false|CC.CC|C1CCC1|false", "true|CC.CC|C1CCC1|false"
  }, delimiter = '|')
  void acceptsOnlyCompleteSelectedRings(boolean smarts, String queryText, String targetText,
                                       boolean expected) throws Exception {
    Substructure.Options options = new Substructure.Options();
    options.chemistry.completeRingsOnly = true;
    options.timeoutMillis = 1000;
    Substructure query = (smarts ? Substructure.smarts(queryText)
        : Substructure.find(new MolGraph(rawMolecule(queryText)))).withOptions(options);
    MolGraph target = new MolGraph(rawMolecule(targetText));
    for (SearchResult result : new SearchResult[] {query.match(target), query.matchAll(target, 100)}) {
      assertEquals(expected, result.found());
      assertEquals(SearchResult.Status.Completed, result.status());
      assertTrue(result.optimal());
      assertTrue(result.complete());
    }
  }

  @Test
  void partialRingEmbeddingsRemainAvailableByDefault() throws Exception {
    MolGraph target = new MolGraph(rawMolecule("C1CCCCC1"));
    assertTrue(Substructure.find(new MolGraph(rawMolecule("C"))).match(target).found());
    assertTrue(Substructure.smarts("C").match(target).found());
  }
}
