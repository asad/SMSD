/* SPDX-License-Identifier: Apache-2.0 */
package com.bioinception.smsd.core;

import static org.junit.jupiter.api.Assertions.*;
import static com.bioinception.smsd.TestSupport.rawGraph;
import org.junit.jupiter.api.Test;

class MCSQuerySmartsResultTest {
  @Test void searchResultCarriesPortableLibraryQuery() throws Exception {
    MolGraph query = rawGraph("CCO");
    MCS.SmartsResult common = MCS.find(query).matchWithSmarts(rawGraph("CCN"));
    assertTrue(common.result().found());
    assertTrue(common.result().optimal());
    assertTrue(Substructure.smarts(common.smarts()).match(rawGraph("CCCl")).found());
  }

  @Test void stoppedSearchHasNoPattern() throws Exception {
    MCS.SmartsResult common = MCS.find(rawGraph("CCO")).withTimeout(0)
        .matchWithSmarts(rawGraph("CCN"));
    assertEquals("", common.smarts());
    assertFalse(common.result().optimal());
  }
}
