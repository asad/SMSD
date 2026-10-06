/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * See the NOTICE file for attribution.
 */
package com.bioinception.smsd.cli;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import com.bioinception.smsd.cli.SMSDcli.MolIO;
import com.fasterxml.jackson.databind.JsonNode;
import com.fasterxml.jackson.databind.ObjectMapper;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.ValueSource;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import picocli.CommandLine;

/** Checks file-reader chemistry and CLI results independently of SMILES input. */
class CLIFileFormatRegressionTest {
  @TempDir Path temporary;

  private Path input(String type, String contents) throws IOException {
    Path directory = Files.createDirectories(temporary.resolve("Molécules with spaces"));
    return Files.write(directory.resolve("éthane." + type.toLowerCase()),
        contents.getBytes(StandardCharsets.UTF_8));
  }

  private static String cmlMolecule(String id) {
    return "<molecule id=\"" + id + "\"><atomArray>"
        + "<atom id=\"a1\" elementType=\"C\" hydrogenCount=\"3\"/>"
        + "<atom id=\"a2\" elementType=\"C\" hydrogenCount=\"3\"/>"
        + "</atomArray><bondArray><bond atomRefs2=\"a1 a2\" order=\"1\"/>"
        + "</bondArray></molecule>";
  }

  private static String cml(String molecules) {
    return "<?xml version=\"1.0\" encoding=\"UTF-8\"?>"
        + "<cml xmlns=\"http://www.xml-cml.org/schema\">" + molecules + "</cml>";
  }

  private static String pdbAtoms() {
    return "HETATM    1  C1  ETH A   1       0.000   0.000   0.000  1.00  0.00           C  \n"
        + "HETATM    2  C2  ETH A   1       1.540   0.000   0.000  1.00  0.00           C  \n"
        + "CONECT    1    2\nCONECT    2    1\n";
  }

  private static String ethane(String type) {
    switch (type) {
      case "CML": return cml(cmlMolecule("ethane"));
      case "PDB": return pdbAtoms() + "END\n";
      case "MOL": return "Ethane\n  SMSD\n\n  2  1  0  0  0  0            999 V2000\n"
          + "    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
          + "    1.5400    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
          + "  1  2  1  0  0  0  0\nM  END\n";
      case "ML2": return "@<TRIPOS>MOLECULE\nEthane\n2 1 0 0 0\nSMALL\nNO_CHARGES\n\n"
          + "@<TRIPOS>ATOM\n1 C1 0.0 0.0 0.0 C.3 1 ETH 0.0\n"
          + "2 C2 1.54 0.0 0.0 C.3 1 ETH 0.0\n@<TRIPOS>BOND\n1 1 2 1\n";
      default: throw new IllegalArgumentException("Unknown fixture format: " + type);
    }
  }

  private CommandLine cli(StringWriter errors) {
    return new CommandLine(new SMSDcli()).setErr(new PrintWriter(errors, true))
        .setOut(new PrintWriter(new StringWriter()));
  }

  @ParameterizedTest
  @ValueSource(strings = {"CML", "PDB", "MOL", "ML2"})
  void readsQueryAndTargetWithElementsAndBond(String type) throws Exception {
    Path file = input(type, ethane(type));
    IAtomContainer target = MolIO.loadTarget(type, file.toString());
    IAtomContainer query = MolIO.loadQuery(type, file.toString()).container();
    for (IAtomContainer molecule : new IAtomContainer[] {target, query}) {
      assertEquals(2, molecule.getAtomCount());
      assertEquals(1, molecule.getBondCount());
      assertEquals(6, molecule.getAtom(0).getAtomicNumber());
      assertEquals(6, molecule.getAtom(1).getAtomicNumber());
      assertEquals(IBond.Order.SINGLE, molecule.getBond(0).getOrder());
      assertTrue(molecule.getBond(0).contains(molecule.getAtom(0)));
      assertTrue(molecule.getBond(0).contains(molecule.getAtom(1)));
    }
  }

  @ParameterizedTest
  @ValueSource(strings = {"CML", "PDB", "MOL", "ML2"})
  void cliMatchesFileTargetAndReturnsNegativeExitCode(String type) throws Exception {
    Path file = input(type, ethane(type));
    Path output = temporary.resolve("résultat.json");
    StringWriter errors = new StringWriter();
    assertEquals(0, cli(errors).execute("--Q", "SMI", "--q", "CC", "--T", type,
        "--t", file.toString(), "--json", output.toString()), errors.toString());
    ObjectMapper mapper = new ObjectMapper();
    assertTrue(mapper.readTree(output.toFile()).path("exists").asBoolean());
    assertEquals(1, cli(errors).execute("--Q", "SMI", "--q", "N", "--T", type,
        "--t", file.toString(), "--json", output.toString()), errors.toString());
    assertFalse(mapper.readTree(output.toFile()).path("exists").asBoolean());
  }

  @ParameterizedTest
  @ValueSource(strings = {"CML", "PDB"})
  void cliUsesFileQueryForMCS(String type) throws Exception {
    Path file = input(type, ethane(type));
    Path output = temporary.resolve("MCS résultat.json");
    StringWriter errors = new StringWriter();
    assertEquals(0, cli(errors).execute("--Q", type, "--q", file.toString(), "--T", "SMI",
        "--t", "CCC", "--mode", "mcs", "--json", output.toString()), errors.toString());
    JsonNode result = new ObjectMapper().readTree(output.toFile());
    assertEquals(2, result.path("mcs_size").asInt());
    assertEquals(2, result.path("pairs").size());
    assertTrue(result.path("mcs_smiles").isTextual());
  }

  @ParameterizedTest
  @ValueSource(strings = {"CML", "PDB"})
  void rejectsEmptyFilesWithoutWritingAResult(String type) throws Exception {
    Path file = input(type, "CML".equals(type) ? cml("") : "END\n");
    IOException error = assertThrows(IOException.class, () -> MolIO.loadTarget(type, file.toString()));
    assertTrue(error.getMessage().contains("molecule"));
    Path output = temporary.resolve("empty-result.json");
    StringWriter errors = new StringWriter();
    assertEquals(1, cli(errors).execute("--Q", "SMI", "--q", "C", "--T", type,
        "--t", file.toString(), "--json", output.toString()));
    assertFalse(Files.exists(output));
    assertTrue(errors.toString().contains("molecule"));
  }

  @ParameterizedTest
  @ValueSource(strings = {"CML", "PDB"})
  void rejectsMultipleMoleculesOrModels(String type) throws Exception {
    String contents = "CML".equals(type) ? cml(cmlMolecule("first") + cmlMolecule("second"))
        : "MODEL        1\n" + pdbAtoms() + "ENDMDL\nMODEL        2\n" + pdbAtoms() + "ENDMDL\nEND\n";
    Path file = input(type, contents);
    IOException error = assertThrows(IOException.class, () -> MolIO.loadTarget(type, file.toString()));
    assertTrue(error.getMessage().contains("exactly one"));
    assertTrue(error.getMessage().contains("found 2"));
    assertThrows(IOException.class, () -> MolIO.loadQuery(type, file.toString()));
    Path output = temporary.resolve("multiple-result.json");
    StringWriter errors = new StringWriter();
    assertEquals(1, cli(errors).execute("--Q", "SMI", "--q", "C", "--T", type,
        "--t", file.toString(), "--json", output.toString()));
    assertFalse(Files.exists(output));
    assertTrue(errors.toString().contains("exactly one"));
  }

  @Test
  void rejectsNonEmptyCMLDocumentWithAnEmptyMolecule() throws Exception {
    Path file = input("CML", cml("<molecule id=\"empty\"/>"));
    IOException error = assertThrows(IOException.class, () -> MolIO.loadTarget("CML", file.toString()));
    assertTrue(error.getMessage().contains("non-empty"));
  }

  @Test
  void rejectsMalformedCML() throws Exception {
    Path file = input("CML", "<cml><molecule>");
    assertThrows(CDKException.class, () -> MolIO.loadTarget("CML", file.toString()));
  }
}
