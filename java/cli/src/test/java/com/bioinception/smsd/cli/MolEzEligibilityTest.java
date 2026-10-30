/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
package com.bioinception.smsd.cli;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;

import com.bioinception.smsd.cli.SMSDcli.MolIO;
import com.bioinception.smsd.core.MolGraph;
import com.fasterxml.jackson.databind.JsonNode;
import com.fasterxml.jackson.databind.ObjectMapper;
import java.io.InputStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IDoubleBondStereochemistry;
import org.openscience.cdk.interfaces.IStereoElement;

/**
 * MOL and SDF blocks give E/Z only to double bonds that can carry it, so a drawn molecule has
 * the canonical SMILES of the same molecule read as SMILES. The cases are shared with the C++
 * and Python tests.
 */
class MolEzEligibilityTest {
  @TempDir Path temporary;

  private static int configuredBonds(IAtomContainer molecule) {
    int count = 0;
    for (IStereoElement<?, ?> element : molecule.stereoElements())
      if (element instanceof IDoubleBondStereochemistry) count++;
    return count;
  }

  @Test
  void sharedCasesGiveTheSmilesConfiguration() throws Exception {
    JsonNode fixture;
    try (InputStream input = MolEzEligibilityTest.class.getResourceAsStream("mol_ez_eligibility.json")) {
      assertNotNull(input, "mol_ez_eligibility.json must be on the test class path");
      fixture = new ObjectMapper().readTree(input);
    }
    for (JsonNode item : fixture.get("cases")) {
      String name = item.get("name").asText();
      StringBuilder block = new StringBuilder();
      for (JsonNode line : item.get("molblock")) block.append(line.asText()).append('\n');
      Path mol = Files.write(temporary.resolve(name + ".mol"),
          block.toString().getBytes(StandardCharsets.UTF_8));
      Path sdf = Files.write(temporary.resolve(name + ".sdf"),
          (block + "$$$$\n").getBytes(StandardCharsets.UTF_8));
      String expected = new MolGraph(MolIO.loadTarget("SMI", item.get("smiles").asText())).toCanonicalSmiles();
      for (IAtomContainer drawn : new IAtomContainer[] {
          MolIO.loadTarget("MOL", mol.toString()), MolIO.loadTargetsSDF(sdf.toString()).get(0)}) {
        assertEquals(item.get("ez_bonds").asInt(), configuredBonds(drawn), name + ": configured double bonds");
        assertEquals(expected, new MolGraph(drawn).toCanonicalSmiles(), name + ": canonical SMILES");
      }
    }
  }
}
