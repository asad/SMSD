/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * See the NOTICE file for attribution.
 */
package com.bioinception.smsd.cli;

import static org.junit.jupiter.api.Assertions.*;
import com.fasterxml.jackson.databind.JsonNode;
import com.fasterxml.jackson.databind.ObjectMapper;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.*;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;
import org.openscience.cdk.io.SDFWriter;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;
import picocli.CommandLine;

class SearchCLIRegressionTest {
  @TempDir Path directory;
  private final ObjectMapper mapper = new ObjectMapper();

  private JsonNode run(int expectedExit, String... arguments) throws Exception {
    Path json = directory.resolve(UUID.randomUUID() + ".json");
    List<String> args = new ArrayList<>(Arrays.asList(arguments));
    args.addAll(Arrays.asList("--json", json.toString()));
    StringWriter errors = new StringWriter();
    assertEquals(expectedExit, new CommandLine(new SMSDcli()).setErr(new PrintWriter(errors, true))
        .execute(args.toArray(new String[0])), errors.toString());
    return mapper.readTree(json.toFile());
  }

  @Test
  void matchesTheDocumentedSmartsAndEnumeratesMappings() throws Exception {
    JsonNode result = run(0, "--Q", "SIG", "--q", "[CX3](=O)[NX3]", "--T", "SMI",
        "--t", "CC(=O)Nc1ccc(O)cc1", "--mappings");
    assertEquals("Completed", result.path("status").asText());
    assertTrue(result.path("complete").asBoolean());
    assertEquals(1, result.path("mappings").size());
    assertEquals(3, result.path("mappings").get(0).path("pairs").size());
    assertEquals(2, result.path("bonds").size());
  }

  @Test
  void exposesTimeoutInsteadOfReportingACompleteNoMatch() throws Exception {
    JsonNode result = run(2, "--Q", "SIG", "--q", "C", "--T", "SMI", "--t", "C", "--timeout", "0");
    assertEquals("TimedOut", result.path("status").asText());
    assertFalse(result.path("complete").asBoolean());
    assertFalse(result.path("exists").asBoolean());
  }

  @Test
  void distinguishesInducedAndEdgeObjectives() throws Exception {
    JsonNode induced = run(0, "--Q", "SMI", "--q", "CCC", "--T", "SMI", "--t", "C1CC1", "--mode", "mcs");
    JsonNode edges = run(0, "--Q", "SMI", "--q", "CCC", "--T", "SMI", "--t", "C1CC1", "--mode", "mcs",
        "--objective", "MCES");
    assertEquals(2, induced.path("mcs_size").asInt());
    assertEquals(1, induced.path("bonds").size());
    assertEquals(3, edges.path("mcs_size").asInt());
    assertEquals(2, edges.path("bonds").size());
    assertEquals(2, new SmilesParser(SilentChemObjectBuilder.getInstance())
        .parseSmiles(edges.path("mcs_smiles").asText()).getBondCount());
  }

  @Test
  void disconnectedFragmentsAreAnExplicitChoice() throws Exception {
    JsonNode connected = run(0, "--Q", "SMI", "--q", "CC.O", "--T", "SMI", "--t", "CC.N.O", "--mode", "mcs");
    JsonNode fragments = run(0, "--Q", "SMI", "--q", "CC.O", "--T", "SMI", "--t", "CC.N.O", "--mode", "mcs", "--fragments");
    assertEquals(2, connected.path("mcs_size").asInt());
    assertEquals(1, connected.path("fragments").size());
    assertEquals(3, fragments.path("mcs_size").asInt());
    assertEquals(2, fragments.path("fragments").size());
  }

  @Test
  void batchSmartsUsesTheSameCompiledQueryAndOutcome() throws Exception {
    Path input = directory.resolve("targets.sdf");
    StringWriter sdf = new StringWriter();
    SmilesParser parser = new SmilesParser(SilentChemObjectBuilder.getInstance());
    try (SDFWriter writer = new SDFWriter(sdf)) {
      writer.write(parser.parseSmiles("CC(=O)N"));
      writer.write(parser.parseSmiles("CCO"));
    }
    Files.write(input, sdf.toString().getBytes(java.nio.charset.StandardCharsets.UTF_8));
    JsonNode result = run(0, "--Q", "SIG", "--q", "[CX3](=O)[NX3]", "--T", "SDF", "--t", input.toString(), "--threads", "2");
    assertTrue(result.path("results").get(0).path("exists").asBoolean());
    assertFalse(result.path("results").get(1).path("exists").asBoolean());
    for (JsonNode row : result.path("results")) assertEquals("Completed", row.path("status").asText());
  }
}
