/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
package com.bioinception.smsd;

import static org.junit.jupiter.api.Assertions.*;
import org.junit.jupiter.api.Test;

class BenchmarkHarnessRegressionTest {
  @Test
  void interruptedBenchmarkDoesNotScheduleAnotherTrial() {
    int completed = 0;
    try {
      BenchmarkRunSettings.checkInterrupted();
      completed++;
      Thread.currentThread().interrupt();
      assertThrows(AssertionError.class, BenchmarkRunSettings::checkInterrupted);
      assertTrue(Thread.currentThread().isInterrupted(), "Stop check must retain interruption");
      assertEquals(1, completed);
    } finally {
      Thread.interrupted();
    }
  }
}
