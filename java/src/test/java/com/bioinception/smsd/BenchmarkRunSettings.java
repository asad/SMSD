/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
package com.bioinception.smsd;

import com.bioinception.smsd.core.ChemOptions;
import com.bioinception.smsd.core.MolGraph;
import com.bioinception.smsd.core.SearchEngine;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardOpenOption;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import org.openscience.cdk.interfaces.IAtomContainer;

/** Shared bounded protocol and flushed checkpoints for opt-in benchmark tests. */
final class BenchmarkRunSettings {
  static final long TIMEOUT_MS = positiveLong("smsd.benchmark.timeoutMs", 1000);
  static final int WARMUP = count("smsd.benchmark.warmup", 0, true);
  static final int ROUNDS = count("smsd.benchmark.rounds", 1, false);
  private static final Set<Path> STARTED = new HashSet<>();

  private BenchmarkRunSettings() {}

  private static long positiveLong(String property, long fallback) {
    long value = Long.parseLong(System.getProperty(property, Long.toString(fallback)));
    if (value <= 0) throw new IllegalArgumentException(property + " must be positive");
    return value;
  }

  private static int count(String property, int fallback, boolean allowZero) {
    int value = Integer.parseInt(System.getProperty(property, Integer.toString(fallback)));
    if (value < (allowZero ? 0 : 1))
      throw new IllegalArgumentException(property + " must be " + (allowZero ? "nonnegative" : "positive"));
    return value;
  }

  static void checkInterrupted() {
    if (Thread.currentThread().isInterrupted())
      throw new AssertionError("Benchmark interrupted; no further pairs or trials will run");
  }

  static void validate(IAtomContainer query, IAtomContainer target,
                       Map<Integer, Integer> mapping, ChemOptions chemistry, String name) {
    validate(new MolGraph(query), new MolGraph(target), mapping, chemistry, name);
  }

  static void validate(MolGraph query, MolGraph target,
                       Map<Integer, Integer> mapping, ChemOptions chemistry, String name) {
    var errors = SearchEngine.validateMapping(query, target, mapping, chemistry);
    if (!errors.isEmpty()) throw new AssertionError(name + " invalid mapping: " + errors);
  }

  static synchronized void checkpoint(String suite, String row, String status, int atoms,
                                      long elapsedNs, String detail) {
    Path output = Path.of(System.getProperty("smsd.benchmark.outputDir", "build/local-benchmarks/java"),
                          suite + ".tsv");
    try {
      Files.createDirectories(output.getParent());
      if (STARTED.add(output)) {
        Files.writeString(output,
            "# timeout_ms=" + TIMEOUT_MS + " warmup=" + WARMUP + " rounds=" + ROUNDS
                + " java=" + System.getProperty("java.version") + "\n"
                + "row\tstatus\tatoms\telapsed_ns\telapsed_budget_crossing\tdetail\n",
            StandardCharsets.UTF_8);
      }
      Files.writeString(output, clean(row) + "\t" + status + "\t" + atoms + "\t" + elapsedNs
          + "\t" + (elapsedNs / 1_000_000.0 > TIMEOUT_MS) + "\t" + clean(detail) + "\n",
          StandardCharsets.UTF_8, StandardOpenOption.APPEND);
    } catch (IOException exception) {
      throw new AssertionError("Cannot write benchmark checkpoint " + output, exception);
    }
    checkInterrupted();
  }

  private static String clean(String text) {
    return text == null ? "" : text.replace('\t', ' ').replace('\n', ' ').replace('\r', ' ');
  }
}
