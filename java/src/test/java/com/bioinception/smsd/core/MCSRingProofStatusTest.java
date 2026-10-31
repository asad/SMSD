/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
package com.bioinception.smsd.core;

import static com.bioinception.smsd.TestSupport.rawGraph;
import static org.junit.jupiter.api.Assertions.*;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.CsvSource;

class MCSRingProofStatusTest {
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
  private static MolGraph stack, cyclobutane;

  @BeforeAll
  static void prepareGraphs() throws Exception {
    stack = rawGraph(STACK_SMILES);
    cyclobutane = rawGraph("C1CCC1");
    assertTrue(stack.ringsTruncated());
    assertFalse(cyclobutane.ringsTruncated());
  }

  @ParameterizedTest
  @CsvSource({"MCIS,Connected,false", "MCIS,Connected,true", "MCIS,Fragments,false",
      "MCIS,Fragments,true", "MCES,Connected,false", "MCES,Connected,true",
      "MCES,Fragments,false", "MCES,Fragments,true"})
  void completeRingSearchDoesNotClaimProof(
      MCS.Objective objective, MCS.Connectivity connectivity, boolean reverse) {
    MCS.Options options = options(objective, connectivity);
    options.chemistry.completeRingsOnly = true;
    MolGraph query = reverse ? cyclobutane : stack;
    MolGraph target = reverse ? stack : cyclobutane;
    SearchResult result = MCS.find(query).withOptions(options).match(target);
    assertUnproved(query, target, options, result);
  }

  @ParameterizedTest
  @CsvSource({"MCIS,Connected,false", "MCIS,Connected,true", "MCIS,Fragments,false",
      "MCIS,Fragments,true", "MCES,Connected,false", "MCES,Connected,true",
      "MCES,Fragments,false", "MCES,Fragments,true"})
  void unconstrainedSearchRemainsComplete(
      MCS.Objective objective, MCS.Connectivity connectivity, boolean reverse) {
    MCS.Options options = options(objective, connectivity);
    MolGraph query = reverse ? cyclobutane : stack;
    MolGraph target = reverse ? stack : cyclobutane;
    SearchResult result = MCS.find(query).withOptions(options).match(target);
    assertCompleted(query, target, options, result);
  }

  @Test
  void strictRingCountsAlsoNeedUntruncatedPerception() {
    MCS.Options options = options(MCS.Objective.MCIS, MCS.Connectivity.Connected);
    options.chemistry.ringFusionMode = ChemOptions.RingFusionMode.STRICT;
    assertFalse(options.chemistry.completeRingsOnly);
    SearchResult result = MCS.find(stack).withOptions(options).match(cyclobutane);
    assertUnproved(stack, cyclobutane, options, result);
  }

  @Test
  void permissiveFusionWithoutCompleteRingsRemainsComplete() {
    MCS.Options options = options(MCS.Objective.MCIS, MCS.Connectivity.Connected);
    options.chemistry.ringFusionMode = ChemOptions.RingFusionMode.PERMISSIVE;
    assertCompleted(stack, cyclobutane, options,
        MCS.find(stack).withOptions(options).match(cyclobutane));
  }

  @Test
  void ordinaryCompleteRingSearchRemainsComplete() {
    MCS.Options options = options(MCS.Objective.MCIS, MCS.Connectivity.Connected);
    options.chemistry.completeRingsOnly = true;
    assertCompleted(cyclobutane, cyclobutane, options,
        MCS.find(cyclobutane).withOptions(options).match(cyclobutane));
  }

  @Test
  void timeoutAndCancellationKeepTheirStatus() {
    MCS.Options options = options(MCS.Objective.MCIS, MCS.Connectivity.Connected);
    options.chemistry.completeRingsOnly = true;
    options.timeoutMillis = 0;
    SearchResult timedOut = MCS.find(stack).withOptions(options).match(cyclobutane);
    options.cancellation = () -> true;
    SearchResult cancelled = MCS.find(stack).withOptions(options).match(cyclobutane);
    assertAll(() -> assertEquals(SearchResult.Status.TimedOut, timedOut.status()),
        () -> assertFalse(timedOut.optimal()), () -> assertFalse(timedOut.complete()),
        () -> assertEquals(SearchResult.Status.Cancelled, cancelled.status()),
        () -> assertFalse(cancelled.optimal()), () -> assertFalse(cancelled.complete()));
  }

  @Test
  void threadedBatchRetainsEveryTarget() {
    MCS.Options options = options(MCS.Objective.MCIS, MCS.Connectivity.Connected);
    options.chemistry.completeRingsOnly = true;
    List<SearchResult> results = MCS.find(cyclobutane).withOptions(options)
        .matchMany(Arrays.asList(stack, cyclobutane), 2);
    assertEquals(2, results.size());
    assertAll(() -> assertUnproved(cyclobutane, stack, options, results.get(0)),
        () -> assertCompleted(cyclobutane, cyclobutane, options, results.get(1)));
  }

  @ParameterizedTest
  @CsvSource({"false", "true"})
  void consensusDoesNotClaimCompleteRingProof(boolean strict) {
    MCS.Options options = options(MCS.Objective.MCIS, MCS.Connectivity.Connected);
    options.chemistry.completeRingsOnly = !strict;
    if (strict) options.chemistry.ringFusionMode = ChemOptions.RingFusionMode.STRICT;
    MCS.ConsensusResult result =
        MCS.consensus(Arrays.asList(cyclobutane, stack), options, 1.0);
    assertAll(() -> assertEquals(SearchResult.Status.LimitReached, result.status()),
        () -> assertFalse(result.optimal()), () -> assertFalse(result.complete()));
  }

  private static MCS.Options options(MCS.Objective objective, MCS.Connectivity connectivity) {
    MCS.Options options = new MCS.Options();
    options.objective = objective;
    options.connectivity = connectivity;
    options.timeoutMillis = 60_000;
    return options;
  }

  private static void assertUnproved(
      MolGraph query, MolGraph target, MCS.Options options, SearchResult result) {
    SearchResult.Match match = result.matches().isEmpty()
        ? new SearchResult.Match(Collections.emptyMap(), Collections.emptyMap(),
            Collections.emptyList(), 0)
        : result.matches().get(0);
    assertAll(() -> assertEquals(SearchResult.Status.LimitReached, result.status()),
        () -> assertFalse(result.optimal()), () -> assertFalse(result.complete()),
        () -> assertTrue(MCS.validateMapping(query, target, match, options).isEmpty()),
        () -> assertThrows(IllegalStateException.class,
            () -> MCS.isMappingMaximal(query, target, match, options)),
        () -> assertThrows(IllegalStateException.class,
            () -> MCS.isMappingMaximal(query, target, match.atomMapping(), options)));
  }

  private static void assertCompleted(
      MolGraph query, MolGraph target, MCS.Options options, SearchResult result) {
    assertAll(() -> assertEquals(SearchResult.Status.Completed, result.status()),
        () -> assertTrue(result.optimal()), () -> assertTrue(result.complete()),
        () -> {
          for (SearchResult.Match match : result.matches())
            assertTrue(MCS.validateMapping(query, target, match, options).isEmpty());
        });
  }
}
