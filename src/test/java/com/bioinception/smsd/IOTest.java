/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
package com.bioinception.smsd;

import static org.junit.jupiter.api.Assertions.*;

import com.bioinception.smsd.cli.SMSDcli.MolIO;
import com.bioinception.smsd.cli.SMSDcli.OutputUtil;
import com.bioinception.smsd.core.ChemOptions;
import com.bioinception.smsd.core.SMSD;
import com.bioinception.smsd.core.Standardiser;
import java.io.ByteArrayOutputStream;
import java.util.*;
import org.junit.jupiter.api.*;
import org.openscience.cdk.interfaces.IAtomContainer;

/**
 * Consolidated I/O tests covering OutputUtil and MolIO.
 *
 * @author Syed Asad Rahman
 */
@DisplayName("I/O Tests")
class IOTest extends TestBase {

  private static IAtomContainer parse(String smi) throws Exception {
    return mol(smi);
  }

  // ======================================================================
  // OutputUtil
  // ======================================================================

  @Nested
  @DisplayName("OutputUtil: JSON Output")
  class JsonOutput {

    @Test
    @DisplayName("JSON contains expected keys")
    void jsonKeys() throws Exception {
      IAtomContainer q = parse("c1ccccc1");
      IAtomContainer t = parse("c1ccc(C)cc1");
      SMSD smsd = new SMSD(q, t, new ChemOptions());
      smsd.isSubstructure();
      List<Map<Integer, Integer>> mappings = smsd.findAllSubstructures(10, 5000);

      ByteArrayOutputStream baos = new ByteArrayOutputStream();
      OutputUtil.writeMappings(q, t, mappings, OutputUtil.OutType.JSON, baos, false);
      String json = baos.toString("UTF-8");

      assertTrue(json.contains("query_size"), "JSON should contain query_size");
      assertTrue(json.contains("target_size"), "JSON should contain target_size");
      assertTrue(json.contains("mappings"), "JSON should contain mappings");
      assertTrue(json.contains("pairs"), "JSON should contain pairs");
      assertTrue(json.contains("query_atoms"), "JSON should contain query_atoms");
      assertTrue(json.contains("target_atoms"), "JSON should contain target_atoms");
    }

    @Test
    @DisplayName("Pretty JSON contains newlines")
    void prettyJson() throws Exception {
      IAtomContainer q = parse("c1ccccc1");
      IAtomContainer t = parse("c1ccc(C)cc1");
      SMSD smsd = new SMSD(q, t, new ChemOptions());
      smsd.isSubstructure();
      List<Map<Integer, Integer>> mappings = smsd.findAllSubstructures(10, 5000);

      ByteArrayOutputStream baos = new ByteArrayOutputStream();
      OutputUtil.writeMappings(q, t, mappings, OutputUtil.OutType.JSON, baos, true);
      String json = baos.toString("UTF-8");

      assertTrue(json.contains("\n"), "Pretty JSON should contain newlines");
      assertTrue(json.contains("  "), "Pretty JSON should contain indentation");
    }

    @Test
    @DisplayName("JSON with similarity upper bound")
    void jsonWithSimilarity() throws Exception {
      IAtomContainer q = parse("c1ccccc1");
      IAtomContainer t = parse("c1ccc(C)cc1");
      List<Map<Integer, Integer>> mappings = new ArrayList<>();
      Map<Integer, Integer> m = new LinkedHashMap<>();
      m.put(0, 0);
      m.put(1, 1);
      m.put(2, 2);
      mappings.add(m);

      ByteArrayOutputStream baos = new ByteArrayOutputStream();
      OutputUtil.writeMappings(q, t, mappings, OutputUtil.OutType.JSON, baos, false, 0.85);
      String json = baos.toString("UTF-8");

      assertTrue(
          json.contains("similarity_upper_bound"), "JSON should contain similarity_upper_bound");
      assertTrue(json.contains("0.85"), "JSON should contain the bound value");
    }

    @Test
    @DisplayName("JSON without similarity (NaN)")
    void jsonWithoutSimilarity() throws Exception {
      IAtomContainer q = parse("c1ccccc1");
      IAtomContainer t = parse("c1ccc(C)cc1");
      List<Map<Integer, Integer>> mappings = Collections.singletonList(Collections.emptyMap());

      ByteArrayOutputStream baos = new ByteArrayOutputStream();
      OutputUtil.writeMappings(q, t, mappings, OutputUtil.OutType.JSON, baos, false, Double.NaN);
      String json = baos.toString("UTF-8");

      assertFalse(
          json.contains("similarity_upper_bound"),
          "JSON should NOT contain similarity_upper_bound when NaN");
    }
  }

  @Nested
  @DisplayName("OutputUtil: SMILES Output")
  class SmilesOutput {

    @Test
    @DisplayName("SMI output produces mapped SMILES")
    void smiOutput() throws Exception {
      IAtomContainer q = parse("c1ccccc1");
      IAtomContainer t = parse("c1ccc(C)cc1");
      SMSD smsd = new SMSD(q, t, new ChemOptions());
      smsd.isSubstructure();
      List<Map<Integer, Integer>> mappings = smsd.findAllSubstructures(1, 5000);

      ByteArrayOutputStream baos = new ByteArrayOutputStream();
      OutputUtil.writeMappings(q, t, mappings, OutputUtil.OutType.SMI, baos, false);
      String output = baos.toString("UTF-8").trim();

      assertFalse(output.isEmpty(), "SMILES output should not be empty");
    }
  }

  @Nested
  @DisplayName("OutputUtil: Edge Cases")
  class OutputEdgeCases {

    @Test
    @DisplayName("Empty mappings -> valid JSON with empty array")
    void emptyMappings() throws Exception {
      IAtomContainer t = parse("c1ccccc1");
      List<Map<Integer, Integer>> mappings = Collections.emptyList();

      ByteArrayOutputStream baos = new ByteArrayOutputStream();
      OutputUtil.writeMappings(null, t, mappings, OutputUtil.OutType.JSON, baos, false);
      String json = baos.toString("UTF-8");

      assertTrue(
          json.contains("\"mappings\":[]") || json.contains("\"mappings\" : [ ]"),
          "Empty mappings should produce empty array");
    }

    @Test
    @DisplayName("Null target -> NullPointerException")
    void nullTarget() {
      assertThrows(
          NullPointerException.class,
          () ->
              OutputUtil.writeMappings(
                  null,
                  null,
                  Collections.emptyList(),
                  OutputUtil.OutType.JSON,
                  new ByteArrayOutputStream(),
                  false));
    }

    @Test
    @DisplayName("Null mappings -> NullPointerException")
    void nullMappings() throws Exception {
      IAtomContainer t = parse("C");
      assertThrows(
          NullPointerException.class,
          () ->
              OutputUtil.writeMappings(
                  null, t, null, OutputUtil.OutType.JSON, new ByteArrayOutputStream(), false));
    }

    @Test
    @DisplayName("writeMolecule SMI format")
    void writeMoleculeSmi() throws Exception {
      IAtomContainer mol = parse("c1ccccc1");
      ByteArrayOutputStream baos = new ByteArrayOutputStream();
      OutputUtil.writeMolecule(mol, OutputUtil.OutType.SMI, baos);
      String output = baos.toString("UTF-8").trim();
      assertFalse(output.isEmpty(), "SMILES output should not be empty");
    }
  }

  // ======================================================================
  // MolIO
  // ======================================================================

  @Nested
  @DisplayName("MolIO: SMILES Loading")
  class MolIOSmilesLoading {

    @Test
    @DisplayName("Load benzene SMILES -> 6 atoms, 6 bonds")
    void loadBenzene() throws Exception {
      MolIO.Query q = MolIO.loadQuery("SMI", "c1ccccc1");
      assertFalse(q.isSmarts());
      assertNotNull(q.container());
      assertEquals(6, q.container().getAtomCount());
      assertEquals(6, q.container().getBondCount());
    }

    @Test
    @DisplayName("Load ethanol SMILES -> 3 heavy atoms")
    void loadEthanol() throws Exception {
      MolIO.Query q = MolIO.loadQuery("SMI", "CCO");
      assertEquals(3, q.container().getAtomCount());
    }

    @Test
    @DisplayName("Load target SMILES -> valid molecule")
    void loadTarget() throws Exception {
      IAtomContainer t = MolIO.loadTarget("SMI", "c1ccccc1");
      assertNotNull(t);
      assertEquals(6, t.getAtomCount());
    }
  }

  @Nested
  @DisplayName("MolIO: SMARTS Loading")
  class MolIOSmartsLoading {

    @Test
    @DisplayName("Load SMARTS -> isSmarts=true, container=null")
    void loadSmarts() throws Exception {
      MolIO.Query q = MolIO.loadQuery("SIG", "[#6]");
      assertTrue(q.isSmarts());
      assertNull(q.container());
      assertEquals("[#6]", q.text());
    }

    @Test
    @DisplayName("Load complex SMARTS")
    void loadComplexSmarts() throws Exception {
      MolIO.Query q = MolIO.loadQuery("SIG", "[#6][#7]");
      assertTrue(q.isSmarts());
      assertEquals("[#6][#7]", q.text());
    }
  }

  @Nested
  @DisplayName("MolIO: Error Handling")
  class MolIOErrorHandling {

    @Test
    @DisplayName("Null query type -> NullPointerException")
    void nullType() {
      assertThrows(NullPointerException.class, () -> MolIO.loadQuery(null, "C"));
    }

    @Test
    @DisplayName("Null query value -> NullPointerException")
    void nullValue() {
      assertThrows(NullPointerException.class, () -> MolIO.loadQuery("SMI", null));
    }

    @Test
    @DisplayName("Null target type -> NullPointerException")
    void nullTargetType() {
      assertThrows(NullPointerException.class, () -> MolIO.loadTarget(null, "C"));
    }

    @Test
    @DisplayName("Null target value -> NullPointerException")
    void nullTargetValue() {
      assertThrows(NullPointerException.class, () -> MolIO.loadTarget("SMI", null));
    }

    @Test
    @DisplayName("Invalid SMILES -> exception")
    void invalidSmiles() {
      assertThrows(Exception.class, () -> MolIO.loadQuery("SMI", "XYZ123!!!"));
    }

    @Test
    @DisplayName("SDF via loadTarget -> IllegalArgumentException")
    void sdfViaLoadTarget() {
      assertThrows(IllegalArgumentException.class, () -> MolIO.loadTarget("SDF", "file.sdf"));
    }

    @Test
    @DisplayName("Unsupported type -> IllegalArgumentException")
    void unsupportedType() {
      assertThrows(Exception.class, () -> MolIO.loadTarget("XYZ", "fake.xyz"));
    }

    @Test
    @DisplayName("Null SDF path -> NullPointerException")
    void nullSdfPath() {
      assertThrows(NullPointerException.class, () -> MolIO.loadTargetsSDF(null));
    }

    @Test
    @DisplayName("Non-existent file -> FileNotFoundException")
    void fileNotFound() {
      assertThrows(
          java.io.FileNotFoundException.class,
          () -> MolIO.loadTarget("MOL", "/nonexistent/path/fake.mol"));
    }
  }

  @Nested
  @DisplayName("MolIO: Roundtrip Tests")
  class MolIORoundtrip {

    @Test
    @DisplayName("Load two SMILES -> substructure check works")
    void roundtripSubstructure() throws Exception {
      MolIO.Query q = MolIO.loadQuery("SMI", "c1ccccc1");
      IAtomContainer target = MolIO.loadTarget("SMI", "c1ccc(C)cc1");
      IAtomContainer qStd = Standardiser.standardise(q.container(), Standardiser.TautomerMode.NONE);
      IAtomContainer tStd = Standardiser.standardise(target, Standardiser.TautomerMode.NONE);
      SMSD smsd = new SMSD(qStd, tStd, new ChemOptions());
      assertTrue(smsd.isSubstructure(), "Benzene should be substructure of toluene");
    }

    @Test
    @DisplayName("SMARTS query roundtrip")
    void smartsRoundtrip() throws Exception {
      MolIO.Query q = MolIO.loadQuery("SIG", "[#6]~[#6]");
      IAtomContainer target = MolIO.loadTarget("SMI", "CC");
      IAtomContainer tStd = Standardiser.standardise(target, Standardiser.TautomerMode.NONE);
      SMSD smsd = new SMSD(q.text(), tStd, new ChemOptions());
      assertTrue(smsd.isSubstructure(), "C~C SMARTS should match ethane");
    }
  }
}
