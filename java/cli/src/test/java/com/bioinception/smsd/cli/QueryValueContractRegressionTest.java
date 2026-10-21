/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
package com.bioinception.smsd.cli;

import static org.junit.jupiter.api.Assertions.*;

import com.bioinception.smsd.cli.SMSDcli.MolIO.Query;
import com.fasterxml.jackson.databind.ObjectMapper;
import org.junit.jupiter.api.Test;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.silent.SilentChemObjectBuilder;

/** Value and JSON contract of the CLI query on Java 8. */
class QueryValueContractRegressionTest {
  private static final ObjectMapper JSON = new ObjectMapper();

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
