/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
package com.bioinception.smsd.core;

import static com.bioinception.smsd.TestSupport.rawGraph;
import static org.junit.jupiter.api.Assertions.*;

import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.CsvSource;

class SubstructureRingProofStatusTest {
  private static final String STACK_SMILES =
      "C12C3C4C5C6C1C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C" +
      "86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17" +
      "C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C4" +
      "5C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C" +
      "12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54" +
      "C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C7" +
      "1C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C" +
      "56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34" +
      "C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C4" +
      "3C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C" +
      "86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17" +
      "C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C4" +
      "5C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C" +
      "12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54" +
      "C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C7" +
      "1C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C" +
      "56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34" +
      "C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C1C6C5C4C3C21";

  private static Substructure query(boolean smarts, String text) throws Exception {
    return smarts ? Substructure.smarts(text) : Substructure.find(rawGraph(text));
  }

  private static Substructure.Options options(boolean strict) {
    Substructure.Options options = new Substructure.Options();
    options.timeoutMillis = 60_000;
    options.chemistry.completeRingsOnly = !strict;
    if (strict) options.chemistry.ringFusionMode = ChemOptions.RingFusionMode.STRICT;
    return options;
  }

  @ParameterizedTest
  @CsvSource({"false,false", "false,true", "true,false", "true,true"})
  void truncatedRingSearchRetainsUnprovedStatus(boolean smarts, boolean strict) throws Exception {
    MolGraph target = rawGraph(STACK_SMILES);
    assertTrue(target.ringsTruncated());
    SearchResult result = query(smarts, "C").withOptions(options(strict)).match(target);
    assertEquals(SearchResult.Status.LimitReached, result.status());
    assertFalse(result.optimal());
    assertFalse(result.complete());
  }

  @ParameterizedTest
  @CsvSource({"false,false", "false,true", "true,false", "true,true"})
  void stoppedRingSearchKeepsItsReason(boolean smarts, boolean strict) throws Exception {
    MolGraph target = rawGraph(STACK_SMILES);
    Substructure.Options options = options(strict);
    options.timeoutMillis = 0;
    SearchResult timedOut = query(smarts, "C").withOptions(options).match(target);
    options.cancellation = () -> true;
    SearchResult cancelled = query(smarts, "C").withOptions(options).match(target);
    assertEquals(SearchResult.Status.TimedOut, timedOut.status());
    assertEquals(SearchResult.Status.Cancelled, cancelled.status());
    assertFalse(timedOut.optimal() || timedOut.complete());
    assertFalse(cancelled.optimal() || cancelled.complete());
  }

  @ParameterizedTest
  @CsvSource({"false,false", "false,true", "true,false", "true,true"})
  void independentElementNegativeRetainsProof(boolean smarts, boolean strict) throws Exception {
    MolGraph target = rawGraph(STACK_SMILES);
    assertTrue(target.ringsTruncated());
    SearchResult result = query(smarts, "N").withOptions(options(strict)).match(target);
    assertFalse(result.found());
    assertEquals(SearchResult.Status.Completed, result.status());
    assertTrue(result.optimal() && result.complete());
  }

  @ParameterizedTest
  @CsvSource({"false,false", "false,true", "true,false", "true,true"})
  void ordinaryWholeRingRetainsProof(boolean smarts, boolean strict) throws Exception {
    MolGraph target = rawGraph("C1CCC1");
    assertFalse(target.ringsTruncated());
    SearchResult result = query(smarts, "C1CCC1").withOptions(options(strict)).match(target);
    assertTrue(result.found());
    assertEquals(SearchResult.Status.Completed, result.status());
    assertTrue(result.optimal() && result.complete());
  }
}
