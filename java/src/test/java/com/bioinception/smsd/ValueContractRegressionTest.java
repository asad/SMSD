/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * See the NOTICE file for attribution.
 */
package com.bioinception.smsd;

import static org.junit.jupiter.api.Assertions.*;

import com.bioinception.smsd.cli.SMSDcli.MolIO.Query;
import com.bioinception.smsd.core.MCS;
import com.bioinception.smsd.core.MolGraph;
import com.bioinception.smsd.core.SearchResult;
import com.fasterxml.jackson.databind.ObjectMapper;
import org.junit.jupiter.api.Test;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.silent.SilentChemObjectBuilder;

/** Shared immutable mapping, completion and JSON contracts on Java 8. */
class ValueContractRegressionTest {
  private static final ObjectMapper JSON = new ObjectMapper();
  private static MolGraph ethane() {
    return new MolGraph.Builder()
        .atomCount(2)
        .atomicNumbers(new int[] {6, 6})
        .neighbors(new int[][] {{1}, {0}})
        .build();
  }

  @Test
  void resultsAndSelectedBondKeysRoundTripThroughJson() throws Exception {
    SearchResult result = MCS.find(ethane()).match(ethane());
    SearchResult restored = JSON.readValue(JSON.writeValueAsString(result), SearchResult.class);
    assertEquals(result, restored);
    assertEquals(restored, result);
    assertEquals(result.hashCode(), restored.hashCode());
    assertEquals(result.toString(), restored.toString());
    assertNotEquals(result, null);
    assertNotEquals(result, result.toString());
    assertEquals(new SearchResult.Bond(0, 1), new SearchResult.Bond(1, 0));
    assertEquals(new SearchResult.Bond(0, 1), SearchResult.Bond.valueOf("0-1"));
  }

  @Test
  void returnedCollectionsAreImmutable() {
    SearchResult result = MCS.find(ethane()).match(ethane());
    SearchResult.Match match = result.matches().get(0);
    assertThrows(UnsupportedOperationException.class, () -> result.matches().clear());
    assertThrows(UnsupportedOperationException.class, () -> match.atomMapping().clear());
    assertThrows(UnsupportedOperationException.class, () -> match.bondMapping().clear());
    assertThrows(UnsupportedOperationException.class, () -> match.fragments().clear());
    assertThrows(UnsupportedOperationException.class, () -> match.fragments().get(0).clear());
    assertEquals(2, match.atomCount());
    assertEquals(1, match.bondCount());
    assertEquals(2.0, match.score());
    assertEquals(SearchResult.Status.Completed, result.status());
    assertTrue(result.optimal());
    assertTrue(result.complete());
    assertTrue(result.nodesVisited() > 0);
    assertTrue(result.elapsedMillis() >= 0);
  }

  @Test
  void completionComponentsParticipateInEquality() throws Exception {
    SearchResult result = MCS.find(ethane()).match(ethane());
    com.fasterxml.jackson.databind.node.ObjectNode node =
        (com.fasterxml.jackson.databind.node.ObjectNode) JSON.valueToTree(result);
    for (String field :
        new String[] {"optimal", "complete", "score", "nodesVisited", "elapsedMillis"}) {
      com.fasterxml.jackson.databind.node.ObjectNode changed = node.deepCopy();
      if (node.get(field).isBoolean())
        changed.put(field, !node.get(field).asBoolean());
      else
        changed.put(field, node.get(field).asDouble() + 1);
      assertNotEquals(result, JSON.treeToValue(changed, SearchResult.class), field);
    }
  }

  @Test
  void queryValueAndJsonContractIsPreserved() throws Exception {
    Query query = new Query(true, "C*N", null);
    assertEquals(JSON.readTree("{\"isSmarts\":true,\"text\":\"C*N\",\"container\":null}"),
        JSON.valueToTree(query));
    assertEquals(query, JSON.readValue(JSON.writeValueAsString(query), Query.class));
    assertEquals(3221768, query.hashCode());
    assertEquals("Query[isSmarts=true, text=C*N, container=null]", query.toString());
    assertTrue(query.isSmarts());
    assertEquals("C*N", query.text());
    assertNull(query.container());
    assertNotEquals(query, new Query(false, "C*N", null));
    assertNotEquals(query, new Query(true, "C*O", null));
    IAtomContainer container = SilentChemObjectBuilder.getInstance().newAtomContainer();
    Query molecule = new Query(true, "C*N", container);
    assertSame(container, molecule.container());
    assertEquals(molecule, new Query(true, "C*N", container));
    assertNotEquals(query, molecule);
    Query empty = new Query(false, null, null);
    assertEquals(1188757, empty.hashCode());
    assertEquals("Query[isSmarts=false, text=null, container=null]", empty.toString());
  }
}
