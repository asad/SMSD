/* SPDX-License-Identifier: Apache-2.0 */
package com.bioinception.smsd.core;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import org.junit.jupiter.api.Test;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;
import static org.junit.jupiter.api.Assertions.*;

class SubstructureBatchRegressionTest {
  private static MolGraph graph(String smiles) throws Exception {
    return new MolGraph(new SmilesParser(SilentChemObjectBuilder.getInstance()).parseSmiles(smiles));
  }

  @Test void moleculeBatchRetainsInputOrderAndCanBeReused() throws Exception {
    Substructure query = Substructure.find(graph("CC"));
    List<MolGraph> targets = Arrays.asList(graph("CCC"), graph("N"), graph("CC"));
    for (int repeat = 0; repeat < 2; repeat++) {
      List<SearchResult> results = query.matchMany(targets);
      assertEquals(3, results.size());
      assertTrue(results.get(0).found());
      assertFalse(results.get(1).found());
      assertTrue(results.get(2).found());
      for (SearchResult result : results) {
        assertEquals(SearchResult.Status.Completed, result.status());
        assertTrue(result.complete());
        if (result.found()) assertEquals(2, result.matches().get(0).atomCount());
      }
      assertThrows(UnsupportedOperationException.class, () -> results.clear());
    }
    assertTrue(query.matchMany(Collections.emptyList()).isEmpty());
  }

  @Test void smartsBatchUsesItsCompiledQueryForEveryTarget() throws Exception {
    List<SearchResult> results = Substructure.smarts("[CX3](=O)[NX3]")
        .matchMany(Arrays.asList(graph("CCO"), graph("CC(=O)N"), graph("NC(=O)N")));
    assertFalse(results.get(0).found());
    assertTrue(results.get(1).found());
    assertTrue(results.get(2).found());
    assertEquals(3, results.get(1).matches().get(0).atomCount());
  }

  @Test void stoppedOutcomesRemainAssociatedWithEveryTarget() throws Exception {
    List<MolGraph> targets = Arrays.asList(graph("CCC"), graph("N"));
    for (SearchResult result : Substructure.find(graph("CC")).withTimeout(0).matchMany(targets)) {
      assertEquals(SearchResult.Status.TimedOut, result.status());
      assertFalse(result.complete());
      assertFalse(result.found());
    }
    Substructure.Options options = new Substructure.Options();
    options.cancellation = () -> true;
    for (SearchResult result : Substructure.find(graph("CC")).withOptions(options).matchMany(targets)) {
      assertEquals(SearchResult.Status.Cancelled, result.status());
      assertFalse(result.complete());
    }
  }
}
