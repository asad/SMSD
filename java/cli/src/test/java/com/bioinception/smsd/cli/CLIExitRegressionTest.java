/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
package com.bioinception.smsd.cli;

import static org.junit.jupiter.api.Assertions.*;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.ValueSource;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;
import org.openscience.cdk.io.SDFWriter;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;
import picocli.CommandLine;

class CLIExitRegressionTest {
  @TempDir Path directory;
  private int execute(String... arguments) {
    return new CommandLine(new SMSDcli()).setErr(new PrintWriter(new StringWriter()))
        .execute(arguments);
  }

  @ParameterizedTest
  @ValueSource(strings = {"single", "mappings", "mcs"})
  void completedNoMatchHasOneExitCode(String operation) {
    List<String> arguments = new ArrayList<>(Arrays.asList(
        "--Q", "SMI", "--q", "N", "--T", "SMI", "--t", "O"));
    if (operation.equals("mcs")) arguments.addAll(Arrays.asList("--mode", "mcs"));
    if (operation.equals("mappings")) arguments.add("--mappings");
    assertEquals(1, execute(arguments.toArray(new String[0])));
  }

  @ParameterizedTest
  @ValueSource(strings = {"C1CC", "not-a-smiles"})
  void malformedInputIsAnErrorRatherThanANonMatch(String query) {
    assertEquals(2, execute("--Q", "SMI", "--q", query, "--T", "SMI", "--t", "CC"));
  }

  @ParameterizedTest
  @ValueSource(strings = {"substructure", "unknown"})
  void unknownModeIsRejectedBeforeSearching(String mode) {
    assertEquals(2, execute("--Q", "SMI", "--q", "C", "--T", "SMI", "--t", "C",
        "--mode", mode));
  }

  @Test
  void missingFileAndUnsupportedQueryModeAreErrors() {
    assertEquals(2, execute("--Q", "MOL", "--q", directory.resolve("missing.mol").toString(),
        "--T", "SMI", "--t", "CC"));
    assertEquals(2, execute("--Q", "SIG", "--q", "C", "--T", "SMI", "--t", "CC",
        "--mode", "mcs"));
    assertEquals(2, execute("--Q", "SIG", "--q", "C", "--T", "SMI", "--t", "CC",
        "--mode", "mcs", "--benchmark"));
  }

  @ParameterizedTest
  @ValueSource(strings = {"nomatch", "timeout"})
  void batchExitDistinguishesNoMatchFromIncompleteSearch(String outcome) throws Exception {
    StringWriter contents = new StringWriter();
    try (SDFWriter writer = new SDFWriter(contents)) {
      writer.write(new SmilesParser(SilentChemObjectBuilder.getInstance()).parseSmiles("O"));
    }
    Path input = directory.resolve("targets.sdf");
    Files.write(input, contents.toString().getBytes(StandardCharsets.UTF_8));
    assertEquals(outcome.equals("timeout") ? 2 : 1, execute(
        "--Q", "SMI", "--q", "N", "--T", "SDF", "--t", input.toString(),
        "--timeout", outcome.equals("timeout") ? "0" : "1000", "--threads", "1"));
  }

  @ParameterizedTest
  @ValueSource(strings = {"match", "nomatch", "timeout"})
  void benchmarkExitUsesItsTimedSearchOutcomes(String outcome) {
    assertEquals(outcome.equals("timeout") ? 2 : outcome.equals("match") ? 0 : 1, execute(
        "--Q", "SMI", "--q", "N", "--T", "SMI", "--t", outcome.equals("match") ? "N" : "O",
        "--benchmark", "--timeout", outcome.equals("timeout") ? "0" : "1000"));
  }
}
