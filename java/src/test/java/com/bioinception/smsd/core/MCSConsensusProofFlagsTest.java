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
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.CsvSource;

class MCSConsensusProofFlagsTest {
  @ParameterizedTest
  @CsvSource({"MCIS,empty,0,0", "MCES,empty,0,0", "MCIS,single,3,1", "MCES,single,3,1",
      "MCIS,supported,2,2", "MCES,supported,2,2", "MCIS,absent,0,0", "MCES,absent,0,0"})
  void completedReductionRetainsItsCompletionFlag(
      MCS.Objective objective, String scenario, int atoms, int support) throws Exception {
    MCS.Options options = new MCS.Options();
    options.objective = objective;
    options.timeoutMillis = 10_000;
    List<MolGraph> molecules;
    double threshold = 1.0;
    if (scenario.equals("empty")) molecules = Collections.emptyList();
    else if (scenario.equals("single")) molecules = Collections.singletonList(rawGraph("CCO"));
    else if (scenario.equals("supported")) {
      molecules = Arrays.asList(rawGraph("N"), rawGraph("CCO"), rawGraph("CCN"));
      threshold = 2.0 / 3.0;
    } else molecules = Arrays.asList(rawGraph("N"), rawGraph("O"));
    MCS.ConsensusResult result = MCS.consensus(molecules, options, threshold);
    assertAll(() -> assertEquals(SearchResult.Status.Completed, result.status()),
        () -> assertFalse(result.optimal()), () -> assertTrue(result.complete()),
        () -> assertEquals(atoms > 0, result.found()),
        () -> assertEquals(atoms, result.graph().n),
        () -> assertEquals(support, result.supportCount()));
    if (scenario.equals("empty")) assertEquals(-1, result.sourceIndex());
    if (scenario.equals("single")) assertEquals(0, result.sourceIndex());
  }

  @ParameterizedTest
  @CsvSource({"MCIS,false", "MCES,false", "MCIS,true", "MCES,true"})
  void stoppedReductionRetainsItsStatus(MCS.Objective objective, boolean cancelled) throws Exception {
    MCS.Options options = new MCS.Options();
    options.objective = objective;
    options.timeoutMillis = cancelled ? 10_000 : 0;
    if (cancelled) options.cancellation = () -> true;
    MCS.ConsensusResult result =
        MCS.consensus(Collections.singletonList(rawGraph("CCO")), options, 1.0);
    assertAll(() -> assertEquals(cancelled ? SearchResult.Status.Cancelled : SearchResult.Status.TimedOut,
        result.status()), () -> assertFalse(result.optimal()), () -> assertFalse(result.complete()));
  }
}
