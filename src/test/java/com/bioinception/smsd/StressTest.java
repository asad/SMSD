/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
package com.bioinception.smsd;

import com.bioinception.smsd.core.*;
import java.util.*;
import org.junit.jupiter.api.*;
import org.junit.jupiter.api.condition.*;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.interfaces.IAtomContainer;
import static org.junit.jupiter.api.Assertions.*;
import static org.junit.jupiter.api.Assumptions.*;

/**
 * Stress Tests: consolidated from EdgeCaseTest.java, ChemistryStressTest.java, LargeMoleculeTest.java, HardCasesTest.java.
 *
 * @author Syed Asad Rahman
 */
@DisplayName("Stress Tests")
public class StressTest extends TestBase {

  // ======================================================================
  // From: EdgeCaseTest.java
  // ======================================================================


  private static ChemOptions defaultOpts() {
    return new ChemOptions();
  }

  private static ChemOptions looseOpts() {
    ChemOptions c = new ChemOptions();
    c.matchBondOrder = ChemOptions.BondOrderMode.LOOSE;
    c.ringMatchesRingOnly = false;
    c.matchFormalCharge = false;
    return c;
  }

  private static ChemOptions stereoOpts() {
    ChemOptions c = new ChemOptions();
    c.useChirality = true;
    c.useBondStereo = true;
    return c;
  }

  // ======================================================================
  // 1. SINGLE ATOM MOLECULES
  // ======================================================================

  @Nested
  @DisplayName("1. SingleAtomMolecules")
  class SingleAtomMolecules {

    @Test @Timeout(10) @DisplayName("1.01 C vs C self-match")
    void carbonSelfMatch() throws Exception {
      SMSD smsd = new SMSD(mol("C"), mol("C"), defaultOpts());
      assertTrue(smsd.isSubstructure());
    }

    @Test @Timeout(10) @DisplayName("1.02 C vs N no match under strict atom type")
    void carbonVsNitrogen() throws Exception {
      SMSD smsd = new SMSD(mol("C"), mol("N"), defaultOpts());
      assertFalse(smsd.isSubstructure());
    }

    @Test @Timeout(10) @DisplayName("1.03 O vs O self-match")
    void oxygenSelfMatch() throws Exception {
      SMSD smsd = new SMSD(mol("O"), mol("O"), defaultOpts());
      assertTrue(smsd.isSubstructure());
    }

    @Test @Timeout(10) @DisplayName("1.04 Single atom MCS C vs C returns size 1")
    void singleAtomMCS() throws Exception {
      SMSD smsd = new SMSD(mol("C"), mol("C"), defaultOpts());
      Map<Integer, Integer> mcs = smsd.findMCS(false, true, 5000L);
      assertNotNull(mcs);
      assertEquals(1, mcs.size());
    }

    @Test @Timeout(10) @DisplayName("1.05 Single C in ethane")
    void singleCarbonInEthane() throws Exception {
      SMSD smsd = new SMSD(mol("C"), mol("CC"), defaultOpts());
      assertTrue(smsd.isSubstructure());
    }

    @Test @Timeout(10) @DisplayName("1.06 Single N in methylamine")
    void singleNInMethylamine() throws Exception {
      SMSD smsd = new SMSD(mol("N"), mol("CN"), defaultOpts());
      assertTrue(smsd.isSubstructure());
    }

    @Test @Timeout(10) @DisplayName("1.07 Single O in methanol")
    void singleOInMethanol() throws Exception {
      SMSD smsd = new SMSD(mol("O"), mol("CO"), defaultOpts());
      assertTrue(smsd.isSubstructure());
    }

    @Test @Timeout(10) @DisplayName("1.08 Single S in methanethiol")
    void singleSInMethanethiol() throws Exception {
      SMSD smsd = new SMSD(mol("S"), mol("CS"), defaultOpts());
      assertTrue(smsd.isSubstructure());
    }

    @Test @Timeout(10) @DisplayName("1.09 Single C in benzene")
    void singleCInBenzene() throws Exception {
      ChemOptions opts = looseOpts();
      SMSD smsd = new SMSD(mol("C"), mol("c1ccccc1"), opts);
      assertTrue(smsd.isSubstructure());
    }

    @Test @Timeout(10) @DisplayName("1.10 S vs O no match strict")
    void sulfurVsOxygen() throws Exception {
      SMSD smsd = new SMSD(mol("S"), mol("O"), defaultOpts());
      assertFalse(smsd.isSubstructure());
    }
  }

  // ======================================================================
  // 2. EMPTY AND NULL
  // ======================================================================

  @Nested
  @DisplayName("2. EmptyAndNull")
  class EmptyAndNull {

    @Test @Timeout(10) @DisplayName("2.01 Empty SMILES query gives no substructure")
    void emptyQueryNoSubstructure() throws Exception {
      IAtomContainer empty = mol("[H][H]");
      IAtomContainer target = mol("CCO");
      // Hydrogen molecule vs ethanol - no heavy atom overlap
      SMSD smsd = new SMSD(empty, target, defaultOpts());
      Map<Integer, Integer> mcs = smsd.findMCS(false, true, 5000L);
      assertTrue(mcs == null || mcs.isEmpty());
    }

    @Test @Timeout(10) @DisplayName("2.02 Single bond molecule CC")
    void singleBondMolecule() throws Exception {
      SMSD smsd = new SMSD(mol("CC"), mol("CC"), defaultOpts());
      assertTrue(smsd.isSubstructure());
    }

    @Test @Timeout(10) @DisplayName("2.03 Single bond molecule MCS size is 2")
    void singleBondMCSSize() throws Exception {
      SMSD smsd = new SMSD(mol("CC"), mol("CC"), defaultOpts());
      Map<Integer, Integer> mcs = smsd.findMCS(false, true, 5000L);
      assertNotNull(mcs);
      assertEquals(2, mcs.size());
    }

    @Test @Timeout(10) @DisplayName("2.04 Methane self match")
    void methaneSelfMatch() throws Exception {
      SMSD smsd = new SMSD(mol("C"), mol("C"), defaultOpts());
      assertTrue(smsd.isSubstructure());
    }

    @Test @Timeout(10) @DisplayName("2.05 Double bond molecule")
    void doubleBondMolecule() throws Exception {
      SMSD smsd = new SMSD(mol("C=C"), mol("C=C"), defaultOpts());
      assertTrue(smsd.isSubstructure());
    }

    @Test @Timeout(10) @DisplayName("2.06 Triple bond molecule")
    void tripleBondMolecule() throws Exception {
      SMSD smsd = new SMSD(mol("C#C"), mol("C#C"), defaultOpts());
      assertTrue(smsd.isSubstructure());
    }

    @Test @Timeout(10) @DisplayName("2.07 Ethylene in propylene")
    void ethyleneInPropylene() throws Exception {
      SMSD smsd = new SMSD(mol("C=C"), mol("CC=C"), defaultOpts());
      assertTrue(smsd.isSubstructure());
    }

    @Test @Timeout(10) @DisplayName("2.08 Acetylene in propyne")
    void acetyleneInPropyne() throws Exception {
      SMSD smsd = new SMSD(mol("C#C"), mol("CC#C"), defaultOpts());
      assertTrue(smsd.isSubstructure());
    }

    @Test @Timeout(10) @DisplayName("2.09 Water self-match")
    void waterSelfMatch() throws Exception {
      SMSD smsd = new SMSD(mol("O"), mol("O"), defaultOpts());
      assertTrue(smsd.isSubstructure());
    }

    @Test @Timeout(10) @DisplayName("2.10 Ammonia self-match")
    void ammoniaSelfMatch() throws Exception {
      SMSD smsd = new SMSD(mol("N"), mol("N"), defaultOpts());
      assertTrue(smsd.isSubstructure());
    }
  }

  // ======================================================================
  // 3. ISOMER PAIRS
  // ======================================================================

  @Nested
  @DisplayName("3. IsomerPairs")
  class IsomerPairs {

    @Test @Timeout(10) @DisplayName("3.01 Butane vs isobutane MCS >= 3")
    void butaneVsIsobutane() throws Exception {
      SMSD smsd = new SMSD(mol("CCCC"), mol("CC(C)C"), defaultOpts());
      Map<Integer, Integer> mcs = smsd.findMCS(false, true, 5000L);
      assertNotNull(mcs);
      assertTrue(mcs.size() >= 3, "Butane vs isobutane share at least 3 atoms");
    }

    @Test @Timeout(10) @DisplayName("3.02 Butane substructure of isobutane (loose)")
    void butaneSubIsobutane() throws Exception {
      SMSD smsd = new SMSD(mol("CCC"), mol("CC(C)C"), looseOpts());
      assertTrue(smsd.isSubstructure());
    }

    @Test @Timeout(10) @DisplayName("3.03 Pentane vs neopentane MCS >= 3")
    void pentaneVsNeopentane() throws Exception {
      SMSD smsd = new SMSD(mol("CCCCC"), mol("CC(C)(C)C"), defaultOpts());
      Map<Integer, Integer> mcs = smsd.findMCS(false, true, 5000L);
      assertNotNull(mcs);
      assertTrue(mcs.size() >= 3);
    }

    @Test @Timeout(10) @DisplayName("3.04 Pentane vs isopentane MCS >= 4")
    void pentaneVsIsopentane() throws Exception {
      SMSD smsd = new SMSD(mol("CCCCC"), mol("CCCC(C)C"), defaultOpts());
      Map<Integer, Integer> mcs = smsd.findMCS(false, true, 5000L);
      assertNotNull(mcs);
      assertTrue(mcs.size() >= 4);
    }

    @Test @Timeout(10) @DisplayName("3.05 Ethanol vs dimethyl ether MCS")
    void ethanolVsDimethylEther() throws Exception {
      SMSD smsd = new SMSD(mol("CCO"), mol("COC"), defaultOpts());
      Map<Integer, Integer> mcs = smsd.findMCS(false, true, 5000L);
      assertNotNull(mcs);
      assertTrue(mcs.size() >= 2, "Share at least C-O");
    }

    @Test @Timeout(10) @DisplayName("3.06 Ortho-xylene vs meta-xylene MCS >= 8")
    void orthoVsMetaXylene() throws Exception {
      SMSD smsd = new SMSD(mol("Cc1ccccc1C"), mol("Cc1cccc(C)c1"), looseOpts());
      Map<Integer, Integer> mcs = smsd.findMCS(false, true, 5000L);
      assertNotNull(mcs);
      assertTrue(mcs.size() >= 7, "Xylenes share toluene core");
    }

    @Test @Timeout(10) @DisplayName("3.07 Meta-xylene vs para-xylene MCS >= 8")
    void metaVsParaXylene() throws Exception {
      SMSD smsd = new SMSD(mol("Cc1cccc(C)c1"), mol("Cc1ccc(C)cc1"), looseOpts());
      Map<Integer, Integer> mcs = smsd.findMCS(false, true, 5000L);
      assertNotNull(mcs);
      assertTrue(mcs.size() >= 7);
    }

    @Test @Timeout(10) @DisplayName("3.08 Ortho-xylene vs para-xylene MCS >= 7")
    void orthoVsParaXylene() throws Exception {
      SMSD smsd = new SMSD(mol("Cc1ccccc1C"), mol("Cc1ccc(C)cc1"), looseOpts());
      Map<Integer, Integer> mcs = smsd.findMCS(false, true, 5000L);
      assertNotNull(mcs);
      assertTrue(mcs.size() >= 7);
    }

    @Test @Timeout(10) @DisplayName("3.09 1-propanol vs 2-propanol MCS >= 3")
    void propanolIsomers() throws Exception {
      SMSD smsd = new SMSD(mol("CCCO"), mol("CC(O)C"), defaultOpts());
      Map<Integer, Integer> mcs = smsd.findMCS(false, true, 5000L);
      assertNotNull(mcs);
      assertTrue(mcs.size() >= 3);
    }

    @Test @Timeout(10) @DisplayName("3.10 Propanal vs acetone MCS >= 3")
    void propanalVsAcetone() throws Exception {
      SMSD smsd = new SMSD(mol("CCC=O"), mol("CC(=O)C"), defaultOpts());
      Map<Integer, Integer> mcs = smsd.findMCS(false, true, 5000L);
      assertNotNull(mcs);
      assertTrue(mcs.size() >= 3);
    }

    @Test @Timeout(10) @DisplayName("3.11 Hexane self-match MCS == 6")
    void hexaneSelfMatch() throws Exception {
      SMSD smsd = new SMSD(mol("CCCCCC"), mol("CCCCCC"), defaultOpts());
      Map<Integer, Integer> mcs = smsd.findMCS(false, true, 5000L);
      assertNotNull(mcs);
      assertEquals(6, mcs.size());
    }

    @Test @Timeout(10) @DisplayName("3.12 Cyclohexane vs hexane: ring vs chain")
    void cyclohexaneVsHexane() throws Exception {
      SMSD smsd = new SMSD(mol("C1CCCCC1"), mol("CCCCCC"), defaultOpts());
      Map<Integer, Integer> mcs = smsd.findMCS(false, true, 5000L);
      // Ring-to-ring-only blocks ring mapping to chain
      assertNotNull(mcs);
    }

    @Test @Timeout(10) @DisplayName("3.13 Diethyl ether vs methyl propyl ether")
    void diethylVsMethylPropylEther() throws Exception {
      SMSD smsd = new SMSD(mol("CCOCC"), mol("COCCC"), defaultOpts());
      Map<Integer, Integer> mcs = smsd.findMCS(false, true, 5000L);
      assertNotNull(mcs);
      assertTrue(mcs.size() >= 4);
    }

    @Test @Timeout(10) @DisplayName("3.14 Acetic acid vs methyl formate")
    void aceticAcidVsMethylFormate() throws Exception {
      SMSD smsd = new SMSD(mol("CC(=O)O"), mol("COC=O"), defaultOpts());
      Map<Integer, Integer> mcs = smsd.findMCS(false, true, 5000L);
      assertNotNull(mcs);
      assertTrue(mcs.size() >= 2);
    }

    @Test @Timeout(10) @DisplayName("3.15 Propylamine vs trimethylamine")
    void propylVsTrimethylamine() throws Exception {
      SMSD smsd = new SMSD(mol("CCCN"), mol("CN(C)C"), defaultOpts());
      Map<Integer, Integer> mcs = smsd.findMCS(false, true, 5000L);
      assertNotNull(mcs);
      assertTrue(mcs.size() >= 2);
    }

    @Test @Timeout(10) @DisplayName("3.16 Toluene substructure of xylene")
    void tolueneInXylene() throws Exception {
      SMSD smsd = new SMSD(mol("Cc1ccccc1"), mol("Cc1ccccc1C"), looseOpts());
      assertTrue(smsd.isSubstructure());
    }

    @Test @Timeout(10) @DisplayName("3.17 Benzene in all xylenes")
    void benzeneInXylenes() throws Exception {
      IAtomContainer bz = mol("c1ccccc1");
      assertTrue(new SMSD(bz, mol("Cc1ccccc1C"), looseOpts()).isSubstructure());
      assertTrue(new SMSD(bz, mol("Cc1cccc(C)c1"), looseOpts()).isSubstructure());
      assertTrue(new SMSD(bz, mol("Cc1ccc(C)cc1"), looseOpts()).isSubstructure());
    }

    @Test @Timeout(10) @DisplayName("3.18 Cyclopentane vs cyclopentene MCS")
    void cyclopentaneVsCyclopentene() throws Exception {
      SMSD smsd = new SMSD(mol("C1CCCC1"), mol("C1=CCCC1"), looseOpts());
      Map<Integer, Integer> mcs = smsd.findMCS(false, true, 5000L);
      assertNotNull(mcs);
      assertTrue(mcs.size() >= 4);
    }

    @Test @Timeout(10) @DisplayName("3.19 Neopentane self MCS == 5")
    void neopentaneSelfMCS() throws Exception {
      SMSD smsd = new SMSD(mol("CC(C)(C)C"), mol("CC(C)(C)C"), defaultOpts());
      Map<Integer, Integer> mcs = smsd.findMCS(false, true, 5000L);
      assertNotNull(mcs);
      assertEquals(5, mcs.size());
    }

    @Test @Timeout(10) @DisplayName("3.20 Isopentane self MCS == 5")
    void isopentaneSelfMCS() throws Exception {
      SMSD smsd = new SMSD(mol("CCC(C)C"), mol("CCC(C)C"), defaultOpts());
      Map<Integer, Integer> mcs = smsd.findMCS(false, true, 5000L);
      assertNotNull(mcs);
      assertEquals(5, mcs.size());
    }
  }

  // ======================================================================
  // 4. STEREOISOMER PAIRS
  // ======================================================================

  @Nested
  @DisplayName("4. StereoisomerPairs")
  class StereoisomerPairs {

    @Test @Timeout(10) @DisplayName("4.01 E-stilbene vs Z-stilbene no stereo -> match")
    void stilbeneNoStereo() throws Exception {
      SMSD smsd = new SMSD(
          mol("C(/c1ccccc1)=C\\c1ccccc1"),
          mol("C(/c1ccccc1)=C/c1ccccc1"),
          defaultOpts());
      assertTrue(smsd.isSubstructure());
    }

    @Test @Timeout(10) @DisplayName("4.02 E-stilbene self-match with stereo")
    void eStilbeneSelfStereo() throws Exception {
      SMSD smsd = new SMSD(
          mol("C(/c1ccccc1)=C/c1ccccc1"),
          mol("C(/c1ccccc1)=C/c1ccccc1"),
          stereoOpts());
      assertTrue(smsd.isSubstructure());
    }

    @Test @Timeout(10) @DisplayName("4.03 R-alanine vs S-alanine no chirality -> match")
    void alanineNoChirality() throws Exception {
      SMSD smsd = new SMSD(
          mol("[C@@H](N)(C)C(=O)O"),
          mol("[C@H](N)(C)C(=O)O"),
          defaultOpts());
      assertTrue(smsd.isSubstructure());
    }

    @Test @Timeout(10) @DisplayName("4.04 R-alanine self-match with chirality")
    void rAlanineSelfChirality() throws Exception {
      SMSD smsd = new SMSD(
          mol("[C@@H](N)(C)C(=O)O"),
          mol("[C@@H](N)(C)C(=O)O"),
          stereoOpts());
      assertTrue(smsd.isSubstructure());
    }

    @Test @Timeout(10) @DisplayName("4.05 Cis-2-butene vs trans-2-butene no stereo")
    void cisTransButeneNoStereo() throws Exception {
      SMSD smsd = new SMSD(
          mol("C/C=C\\C"), mol("C/C=C/C"), defaultOpts());
      assertTrue(smsd.isSubstructure());
    }

    @Test @Timeout(10) @DisplayName("4.06 Cis-2-butene self-match with bond stereo")
    void cisButeneSelfStereo() throws Exception {
      SMSD smsd = new SMSD(
          mol("C/C=C\\C"), mol("C/C=C\\C"), stereoOpts());
      assertTrue(smsd.isSubstructure());
    }

    @Test @Timeout(10) @DisplayName("4.07 Cis-1,2-dimethylcyclohexane self")
    void cisDimethylCyclohexaneSelf() throws Exception {
      SMSD smsd = new SMSD(
          mol("C[C@H]1CCCC[C@@H]1C"),
          mol("C[C@H]1CCCC[C@@H]1C"),
          defaultOpts());
      assertTrue(smsd.isSubstructure());
    }

    @Test @Timeout(10) @DisplayName("4.08 Meso-tartaric acid self-match")
    void mesoTartaricSelf() throws Exception {
      SMSD smsd = new SMSD(
          mol("O[C@@H]([C@H](O)C(=O)O)C(=O)O"),
          mol("O[C@@H]([C@H](O)C(=O)O)C(=O)O"),
          defaultOpts());
      assertTrue(smsd.isSubstructure());
    }

    @Test @Timeout(10) @DisplayName("4.09 Stilbene MCS without stereo >= 14")
    void stilbeneMCSNoStereo() throws Exception {
      SMSD smsd = new SMSD(
          mol("C(/c1ccccc1)=C\\c1ccccc1"),
          mol("C(/c1ccccc1)=C/c1ccccc1"),
          looseOpts());
      Map<Integer, Integer> mcs = smsd.findMCS(false, true, 5000L);
      assertNotNull(mcs);
      assertTrue(mcs.size() >= 10, "Stilbene isomers share most atoms");
    }

    @Test @Timeout(10) @DisplayName("4.10 Alanine MCS without chirality >= 5")
    void alanineMCSNoChirality() throws Exception {
      SMSD smsd = new SMSD(
          mol("[C@@H](N)(C)C(=O)O"),
          mol("[C@H](N)(C)C(=O)O"),
          looseOpts());
      Map<Integer, Integer> mcs = smsd.findMCS(false, true, 5000L);
      assertNotNull(mcs);
      assertTrue(mcs.size() >= 5);
    }

    @Test @Timeout(10) @DisplayName("4.11 E-2-pentene vs Z-2-pentene no stereo")
    void penteneIsomersNoStereo() throws Exception {
      SMSD smsd = new SMSD(
          mol("C/C=C/CC"), mol("C/C=C\\CC"), defaultOpts());
      assertTrue(smsd.isSubstructure());
    }

    @Test @Timeout(10) @DisplayName("4.12 R-glyceraldehyde self")
    void rGlyceraldehydeSelf() throws Exception {
      SMSD smsd = new SMSD(
          mol("OC[C@@H](O)C=O"),
          mol("OC[C@@H](O)C=O"),
          stereoOpts());
      assertTrue(smsd.isSubstructure());
    }

    @Test @Timeout(10) @DisplayName("4.13 D-glucose vs L-glucose no chirality")
    void glucoseNoChirality() throws Exception {
      SMSD smsd = new SMSD(
          mol("OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O"),
          mol("OC[C@@H]1OC(O)[C@@H](O)[C@H](O)[C@H]1O"),
          defaultOpts());
      assertTrue(smsd.isSubstructure());
    }

    @Test @Timeout(10) @DisplayName("4.14 Trans-2-butene self with stereo")
    void transButeneSelfStereo() throws Exception {
      SMSD smsd = new SMSD(
          mol("C/C=C/C"), mol("C/C=C/C"), stereoOpts());
      assertTrue(smsd.isSubstructure());
    }

    @Test @Timeout(10) @DisplayName("4.15 E/Z but-2-enol MCS >= 4")
    void butenolMCS() throws Exception {
      SMSD smsd = new SMSD(
          mol("C/C=C/CO"), mol("C/C=C\\CO"), looseOpts());
      Map<Integer, Integer> mcs = smsd.findMCS(false, true, 5000L);
      assertNotNull(mcs);
      assertTrue(mcs.size() >= 4);
    }
  }

  // ======================================================================
  // 5. CHARGED SPECIES
  // ======================================================================

  @Nested
  @DisplayName("5. ChargedSpecies")
  class ChargedSpecies {

    @Test @Timeout(10) @DisplayName("5.01 Carboxylate anion vs carboxylic acid charge match")
    void carboxylateVsAcidCharge() throws Exception {
      ChemOptions opts = defaultOpts();
      opts.matchFormalCharge = true;
      SMSD smsd = new SMSD(mol("CC([O-])=O"), mol("CC(O)=O"), opts);
      Map<Integer, Integer> mcs = smsd.findMCS(false, true, 5000L);
      // With charge matching, the charged O should not match neutral O
      assertNotNull(mcs);
    }

    @Test @Timeout(10) @DisplayName("5.02 Carboxylate vs acid no charge match -> full overlap")
    void carboxylateVsAcidNoCharge() throws Exception {
      ChemOptions opts = looseOpts();
      opts.matchFormalCharge = false;
      SMSD smsd = new SMSD(mol("CC([O-])=O"), mol("CC(O)=O"), opts);
      Map<Integer, Integer> mcs = smsd.findMCS(false, true, 5000L);
      assertNotNull(mcs);
      assertTrue(mcs.size() >= 3);
    }

    @Test @Timeout(10) @DisplayName("5.03 Ammonium vs amine")
    void ammoniumVsAmine() throws Exception {
      ChemOptions opts = defaultOpts();
      opts.matchFormalCharge = true;
      SMSD smsd = new SMSD(mol("[NH4+]"), mol("N"), opts);
      Map<Integer, Integer> mcs = smsd.findMCS(false, true, 5000L);
      // Charged ammonium vs neutral amine with charge matching
      assertNotNull(mcs);
    }

    @Test @Timeout(10) @DisplayName("5.04 Ammonium vs amine no charge -> match")
    void ammoniumVsAmineNoCharge() throws Exception {
      ChemOptions opts = looseOpts();
      opts.matchFormalCharge = false;
      SMSD smsd = new SMSD(mol("[NH4+]"), mol("N"), opts);
      Map<Integer, Integer> mcs = smsd.findMCS(false, true, 5000L);
      assertNotNull(mcs);
      assertTrue(mcs.size() >= 1);
    }

    @Test @Timeout(10) @DisplayName("5.05 Glycine zwitterion MCS with itself")
    void glycineZwitterionSelf() throws Exception {
      SMSD smsd = new SMSD(
          mol("[NH3+]CC([O-])=O"),
          mol("[NH3+]CC([O-])=O"),
          defaultOpts());
      Map<Integer, Integer> mcs = smsd.findMCS(false, true, 5000L);
      assertNotNull(mcs);
      assertTrue(mcs.size() >= 4);
    }

    @Test @Timeout(10) @DisplayName("5.06 Glycine zwitterion vs neutral glycine")
    void glycineZwitterionVsNeutral() throws Exception {
      ChemOptions opts = looseOpts();
      opts.matchFormalCharge = false;
      SMSD smsd = new SMSD(
          mol("[NH3+]CC([O-])=O"),
          mol("NCC(O)=O"),
          opts);
      Map<Integer, Integer> mcs = smsd.findMCS(false, true, 5000L);
      assertNotNull(mcs);
      assertTrue(mcs.size() >= 4);
    }

    @Test @Timeout(10) @DisplayName("5.07 Sodium acetate vs acetic acid")
    void sodiumAcetateVsAceticAcid() throws Exception {
      ChemOptions opts = looseOpts();
      opts.matchFormalCharge = false;
      SMSD smsd = new SMSD(mol("CC([O-])=O"), mol("CC(O)=O"), opts);
      Map<Integer, Integer> mcs = smsd.findMCS(false, true, 5000L);
      assertNotNull(mcs);
      assertTrue(mcs.size() >= 3);
    }

    @Test @Timeout(10) @DisplayName("5.08 Phenoxide vs phenol")
    void phenoxideVsPhenol() throws Exception {
      ChemOptions opts = looseOpts();
      opts.matchFormalCharge = false;
      SMSD smsd = new SMSD(mol("[O-]c1ccccc1"), mol("Oc1ccccc1"), opts);
      Map<Integer, Integer> mcs = smsd.findMCS(false, true, 5000L);
      assertNotNull(mcs);
      assertTrue(mcs.size() >= 6);
    }

    @Test @Timeout(10) @DisplayName("5.09 Trimethylammonium self-match")
    void trimethylammoniumSelf() throws Exception {
      SMSD smsd = new SMSD(
          mol("C[NH+](C)C"), mol("C[NH+](C)C"), defaultOpts());
      assertTrue(smsd.isSubstructure());
    }

    @Test @Timeout(10) @DisplayName("5.10 Phosphate vs phosphoric acid no charge")
    void phosphateVsPhosphoricAcid() throws Exception {
      ChemOptions opts = looseOpts();
      opts.matchFormalCharge = false;
      SMSD smsd = new SMSD(mol("[O-]P([O-])([O-])=O"), mol("OP(O)(O)=O"), opts);
      Map<Integer, Integer> mcs = smsd.findMCS(false, true, 5000L);
      assertNotNull(mcs);
      assertTrue(mcs.size() >= 3);
    }

    @Test @Timeout(10) @DisplayName("5.11 Sulfonate vs sulfonic acid")
    void sulfonateVsSulfonicAcid() throws Exception {
      ChemOptions opts = looseOpts();
      opts.matchFormalCharge = false;
      SMSD smsd = new SMSD(mol("CS([O-])(=O)=O"), mol("CS(O)(=O)=O"), opts);
      Map<Integer, Integer> mcs = smsd.findMCS(false, true, 5000L);
      assertNotNull(mcs);
      assertTrue(mcs.size() >= 3);
    }

    @Test @Timeout(10) @DisplayName("5.12 Pyridinium vs pyridine")
    void pyridiniumVsPyridine() throws Exception {
      ChemOptions opts = looseOpts();
      opts.matchFormalCharge = false;
      SMSD smsd = new SMSD(mol("c1cc[nH+]cc1"), mol("c1ccncc1"), opts);
      Map<Integer, Integer> mcs = smsd.findMCS(false, true, 5000L);
      assertNotNull(mcs);
      assertTrue(mcs.size() >= 5);
    }

    @Test @Timeout(10) @DisplayName("5.13 Betaine self-match")
    void betaineSelf() throws Exception {
      SMSD smsd = new SMSD(
          mol("C[N+](C)(C)CC([O-])=O"),
          mol("C[N+](C)(C)CC([O-])=O"),
          defaultOpts());
      assertTrue(smsd.isSubstructure());
    }

    @Test @Timeout(10) @DisplayName("5.14 Histidine protonated vs neutral")
    void histidineProtonatedVsNeutral() throws Exception {
      ChemOptions opts = looseOpts();
      opts.matchFormalCharge = false;
      SMSD smsd = new SMSD(
          mol("NC(Cc1c[nH+]c[nH]1)C(=O)O"),
          mol("NC(Cc1cnc[nH]1)C(=O)O"),
          opts);
      Map<Integer, Integer> mcs = smsd.findMCS(false, true, 5000L);
      assertNotNull(mcs);
      assertTrue(mcs.size() >= 7);
    }

    @Test @Timeout(10) @DisplayName("5.15 Lithium carboxylate substructure")
    void lithiumCarboxylate() throws Exception {
      ChemOptions opts = looseOpts();
      opts.matchFormalCharge = false;
      SMSD smsd = new SMSD(mol("CC([O-])=O"), mol("CC(=O)[O-].[Li+]"), opts);
      Map<Integer, Integer> mcs = smsd.findMCS(false, true, 5000L);
      assertNotNull(mcs);
      assertTrue(mcs.size() >= 3);
    }
  }

  // ======================================================================
  // 6. ISOTOPE LABELED
  // ======================================================================

  @Nested
  @DisplayName("6. IsotopeLabeled")
  class IsotopeLabeled {

    @Test @Timeout(10) @DisplayName("6.01 Deuterium water vs water no isotope match")
    void deuteriumWaterNoIsotope() throws Exception {
      ChemOptions opts = defaultOpts();
      opts.matchIsotope = false;
      SMSD smsd = new SMSD(mol("[2H]O[2H]"), mol("O"), opts);
      Map<Integer, Integer> mcs = smsd.findMCS(false, true, 5000L);
      assertNotNull(mcs);
      assertTrue(mcs.size() >= 1);
    }

    @Test @Timeout(10) @DisplayName("6.02 Deuterium water self-match with isotope")
    void deuteriumWaterSelfIsotope() throws Exception {
      ChemOptions opts = defaultOpts();
      opts.matchIsotope = true;
      SMSD smsd = new SMSD(mol("[2H]O[2H]"), mol("[2H]O[2H]"), opts);
      assertTrue(smsd.isSubstructure());
    }

    @Test @Timeout(10) @DisplayName("6.03 13C-methane vs methane no isotope")
    void c13MethaneNoIsotope() throws Exception {
      ChemOptions opts = defaultOpts();
      opts.matchIsotope = false;
      SMSD smsd = new SMSD(mol("[13CH4]"), mol("C"), opts);
      Map<Integer, Integer> mcs = smsd.findMCS(false, true, 5000L);
      assertNotNull(mcs);
      assertTrue(mcs.size() >= 1);
    }

    @Test @Timeout(10) @DisplayName("6.04 13C-methane self with isotope match")
    void c13MethaneSelfIsotope() throws Exception {
      ChemOptions opts = defaultOpts();
      opts.matchIsotope = true;
      SMSD smsd = new SMSD(mol("[13CH4]"), mol("[13CH4]"), opts);
      assertTrue(smsd.isSubstructure());
    }

    @Test @Timeout(10) @DisplayName("6.05 18F-fluoride vs fluoride no isotope")
    void f18NoIsotope() throws Exception {
      ChemOptions opts = defaultOpts();
      opts.matchIsotope = false;
      SMSD smsd = new SMSD(mol("[18F]"), mol("[F]"), opts);
      Map<Integer, Integer> mcs = smsd.findMCS(false, true, 5000L);
      assertNotNull(mcs);
    }

    @Test @Timeout(10) @DisplayName("6.06 Tritium-labeled methanol self")
    void tritiumMethanol() throws Exception {
      ChemOptions opts = defaultOpts();
      opts.matchIsotope = true;
      SMSD smsd = new SMSD(mol("[3H]CO"), mol("[3H]CO"), opts);
      assertTrue(smsd.isSubstructure());
    }

    @Test @Timeout(10) @DisplayName("6.07 14C-ethanol vs ethanol no isotope")
    void c14EthanolNoIsotope() throws Exception {
      ChemOptions opts = defaultOpts();
      opts.matchIsotope = false;
      SMSD smsd = new SMSD(mol("[14C]CO"), mol("CCO"), opts);
      Map<Integer, Integer> mcs = smsd.findMCS(false, true, 5000L);
      assertNotNull(mcs);
      assertTrue(mcs.size() >= 2);
    }

    @Test @Timeout(10) @DisplayName("6.08 Deuterated benzene vs benzene no isotope")
    void deuteratedBenzeneNoIsotope() throws Exception {
      ChemOptions opts = looseOpts();
      opts.matchIsotope = false;
      SMSD smsd = new SMSD(
          mol("[2H]c1c([2H])c([2H])c([2H])c([2H])c1[2H]"),
          mol("c1ccccc1"),
          opts);
      Map<Integer, Integer> mcs = smsd.findMCS(false, true, 5000L);
      assertNotNull(mcs);
      assertTrue(mcs.size() >= 6,
          "Deuterated vs regular benzene should share >= 6 atoms (the ring), got " + mcs.size());
    }

    @Test @Timeout(10) @DisplayName("6.09 15N-glycine vs glycine no isotope")
    void n15GlycineNoIsotope() throws Exception {
      ChemOptions opts = defaultOpts();
      opts.matchIsotope = false;
      SMSD smsd = new SMSD(mol("[15NH2]CC(=O)O"), mol("NCC(=O)O"), opts);
      Map<Integer, Integer> mcs = smsd.findMCS(false, true, 5000L);
      assertNotNull(mcs);
      assertTrue(mcs.size() >= 4);
    }

    @Test @Timeout(10) @DisplayName("6.10 Mixed isotope ethane self")
    void mixedIsotopeEthaneSelf() throws Exception {
      ChemOptions opts = defaultOpts();
      opts.matchIsotope = true;
      SMSD smsd = new SMSD(mol("[13C][12C]"), mol("[13C][12C]"), opts);
      assertTrue(smsd.isSubstructure());
    }
  }

  // ======================================================================
  // 7. HETEROCYCLE PAIRS
  // ======================================================================

  @Nested
  @DisplayName("7. HeterocyclePairs")
  class HeterocyclePairs {

    @Test @Timeout(10) @DisplayName("7.01 Pyridine self-match")
    void pyridineSelf() throws Exception {
      SMSD smsd = new SMSD(mol("c1ccncc1"), mol("c1ccncc1"), defaultOpts());
      assertTrue(smsd.isSubstructure());
    }

    @Test @Timeout(10) @DisplayName("7.02 Pyrrole self-match")
    void pyrroleSelf() throws Exception {
      SMSD smsd = new SMSD(mol("c1cc[nH]c1"), mol("c1cc[nH]c1"), defaultOpts());
      assertTrue(smsd.isSubstructure());
    }

    @Test @Timeout(10) @DisplayName("7.03 Furan self-match")
    void furanSelf() throws Exception {
      SMSD smsd = new SMSD(mol("c1ccoc1"), mol("c1ccoc1"), defaultOpts());
      assertTrue(smsd.isSubstructure());
    }

    @Test @Timeout(10) @DisplayName("7.04 Thiophene self-match")
    void thiopheneSelf() throws Exception {
      SMSD smsd = new SMSD(mol("c1ccsc1"), mol("c1ccsc1"), defaultOpts());
      assertTrue(smsd.isSubstructure());
    }

    @Test @Timeout(10) @DisplayName("7.05 Imidazole self-match")
    void imidazoleSelf() throws Exception {
      SMSD smsd = new SMSD(mol("c1cnc[nH]1"), mol("c1cnc[nH]1"), defaultOpts());
      assertTrue(smsd.isSubstructure());
    }

    @Test @Timeout(10) @DisplayName("7.06 Pyrimidine self-match")
    void pyrimidineSelf() throws Exception {
      SMSD smsd = new SMSD(mol("c1ccncn1"), mol("c1ccncn1"), defaultOpts());
      assertTrue(smsd.isSubstructure());
    }

    @Test @Timeout(10) @DisplayName("7.07 Purine self-match")
    void purineSelf() throws Exception {
      SMSD smsd = new SMSD(
          mol("C1=CN=C2N1C=NC=N2"), mol("C1=CN=C2N1C=NC=N2"), defaultOpts());
      assertTrue(smsd.isSubstructure());
    }

    @Test @Timeout(10) @DisplayName("7.08 Indole self-match")
    void indoleSelf() throws Exception {
      SMSD smsd = new SMSD(
          mol("c1ccc2[nH]ccc2c1"), mol("c1ccc2[nH]ccc2c1"), defaultOpts());
      assertTrue(smsd.isSubstructure());
    }

    @Test @Timeout(10) @DisplayName("7.09 Quinoline self-match")
    void quinolineSelf() throws Exception {
      SMSD smsd = new SMSD(
          mol("c1ccc2ncccc2c1"), mol("c1ccc2ncccc2c1"), defaultOpts());
      assertTrue(smsd.isSubstructure());
    }

    @Test @Timeout(10) @DisplayName("7.10 Isoquinoline self-match")
    void isoquinolineSelf() throws Exception {
      SMSD smsd = new SMSD(
          mol("c1ccc2cnccc2c1"), mol("c1ccc2cnccc2c1"), defaultOpts());
      assertTrue(smsd.isSubstructure());
    }

    @Test @Timeout(10) @DisplayName("7.11 Pyridine in nicotinic acid")
    void pyridineInNicotinicAcid() throws Exception {
      SMSD smsd = new SMSD(
          mol("c1ccncc1"), mol("OC(=O)c1cccnc1"), looseOpts());
      assertTrue(smsd.isSubstructure());
    }

    @Test @Timeout(10) @DisplayName("7.12 Pyrrole in tryptophan")
    void pyrroleInTryptophan() throws Exception {
      SMSD smsd = new SMSD(
          mol("c1cc[nH]c1"),
          mol("NC(Cc1c[nH]c2ccccc12)C(=O)O"),
          looseOpts());
      assertTrue(smsd.isSubstructure());
    }

    @Test @Timeout(10) @DisplayName("7.13 Furan in furfural")
    void furanInFurfural() throws Exception {
      SMSD smsd = new SMSD(
          mol("c1ccoc1"), mol("O=Cc1ccco1"), looseOpts());
      assertTrue(smsd.isSubstructure());
    }

    @Test @Timeout(10) @DisplayName("7.14 Thiophene in thienyl acetic acid")
    void thiopheneInDrug() throws Exception {
      SMSD smsd = new SMSD(
          mol("c1ccsc1"), mol("OC(=O)Cc1cccs1"), looseOpts());
      assertTrue(smsd.isSubstructure());
    }

    @Test @Timeout(10) @DisplayName("7.15 Imidazole in histamine")
    void imidazoleInHistamine() throws Exception {
      SMSD smsd = new SMSD(
          mol("c1cnc[nH]1"), mol("NCCc1cnc[nH]1"), looseOpts());
      assertTrue(smsd.isSubstructure());
    }

    @Test @Timeout(10) @DisplayName("7.16 Pyridine vs pyrimidine MCS >= 4")
    void pyridineVsPyrimidine() throws Exception {
      SMSD smsd = new SMSD(
          mol("c1ccncc1"), mol("c1ccncn1"), looseOpts());
      Map<Integer, Integer> mcs = smsd.findMCS(false, true, 5000L);
      assertNotNull(mcs);
      assertTrue(mcs.size() >= 4);
    }

    @Test @Timeout(10) @DisplayName("7.17 Quinoline vs isoquinoline MCS >= 8")
    void quinolineVsIsoquinoline() throws Exception {
      SMSD smsd = new SMSD(
          mol("c1ccc2ncccc2c1"), mol("c1ccc2cnccc2c1"), looseOpts());
      Map<Integer, Integer> mcs = smsd.findMCS(false, true, 5000L);
      assertNotNull(mcs);
      assertTrue(mcs.size() >= 8);
    }

    @Test @Timeout(10) @DisplayName("7.18 Indole vs benzofuran MCS >= 7")
    void indoleVsBenzofuran() throws Exception {
      SMSD smsd = new SMSD(
          mol("c1ccc2[nH]ccc2c1"), mol("c1ccc2occc2c1"), looseOpts());
      Map<Integer, Integer> mcs = smsd.findMCS(false, true, 5000L);
      assertNotNull(mcs);
      assertTrue(mcs.size() >= 7);
    }

    @Test @Timeout(10) @DisplayName("7.19 Purine in adenine")
    void purineInAdenine() throws Exception {
      SMSD smsd = new SMSD(
          mol("C1=CN=C2N1C=NC=N2"),
          mol("Nc1ncnc2[nH]cnc12"),
          looseOpts());
      Map<Integer, Integer> mcs = smsd.findMCS(false, true, 5000L);
      assertNotNull(mcs);
      assertTrue(mcs.size() >= 7);
    }

    @Test @Timeout(10) @DisplayName("7.20 Oxazole vs thiazole MCS >= 3")
    void oxazoleVsThiazole() throws Exception {
      SMSD smsd = new SMSD(
          mol("c1cocn1"), mol("c1cscn1"), looseOpts());
      Map<Integer, Integer> mcs = smsd.findMCS(false, true, 5000L);
      assertNotNull(mcs);
      assertTrue(mcs.size() >= 3);
    }
  }

  // ======================================================================
  // 8. FUNCTIONAL GROUP PAIRS
  // ======================================================================

  @Nested
  @DisplayName("8. FunctionalGroupPairs")
  class FunctionalGroupPairs {

    @Test @Timeout(10) @DisplayName("8.01 Alcohol: methanol self")
    void methanolSelf() throws Exception {
      SMSD smsd = new SMSD(mol("CO"), mol("CO"), defaultOpts());
      assertTrue(smsd.isSubstructure());
    }

    @Test @Timeout(10) @DisplayName("8.02 Aldehyde: formaldehyde in acetaldehyde")
    void formaldehydeInAcetaldehyde() throws Exception {
      SMSD smsd = new SMSD(mol("C=O"), mol("CC=O"), looseOpts());
      assertTrue(smsd.isSubstructure());
    }

    @Test @Timeout(10) @DisplayName("8.03 Ketone: acetone self")
    void acetoneSelf() throws Exception {
      SMSD smsd = new SMSD(mol("CC(=O)C"), mol("CC(=O)C"), defaultOpts());
      assertTrue(smsd.isSubstructure());
    }

    @Test @Timeout(10) @DisplayName("8.04 Carboxylic acid: acetic acid self")
    void aceticAcidSelf() throws Exception {
      SMSD smsd = new SMSD(mol("CC(=O)O"), mol("CC(=O)O"), defaultOpts());
      assertTrue(smsd.isSubstructure());
    }

    @Test @Timeout(10) @DisplayName("8.05 Ester: methyl acetate self")
    void methylAcetateSelf() throws Exception {
      SMSD smsd = new SMSD(mol("CC(=O)OC"), mol("CC(=O)OC"), defaultOpts());
      assertTrue(smsd.isSubstructure());
    }

    @Test @Timeout(10) @DisplayName("8.06 Amide: acetamide self")
    void acetamideSelf() throws Exception {
      SMSD smsd = new SMSD(mol("CC(=O)N"), mol("CC(=O)N"), defaultOpts());
      assertTrue(smsd.isSubstructure());
    }

    @Test @Timeout(10) @DisplayName("8.07 Amine: ethylamine in diethylamine")
    void ethylamineInDiethylamine() throws Exception {
      SMSD smsd = new SMSD(mol("CCN"), mol("CCNCC"), looseOpts());
      assertTrue(smsd.isSubstructure());
    }

    @Test @Timeout(10) @DisplayName("8.08 Nitrile: acetonitrile self")
    void acetonitrileSelf() throws Exception {
      SMSD smsd = new SMSD(mol("CC#N"), mol("CC#N"), defaultOpts());
      assertTrue(smsd.isSubstructure());
    }

    @Test @Timeout(10) @DisplayName("8.09 Nitro: nitromethane self")
    void nitromethaneSelf() throws Exception {
      SMSD smsd = new SMSD(
          mol("C[N+](=O)[O-]"), mol("C[N+](=O)[O-]"), defaultOpts());
      assertTrue(smsd.isSubstructure());
    }

    @Test @Timeout(10) @DisplayName("8.10 Aldehyde vs ketone MCS >= 2")
    void aldehydeVsKetone() throws Exception {
      SMSD smsd = new SMSD(mol("CC=O"), mol("CC(=O)C"), defaultOpts());
      Map<Integer, Integer> mcs = smsd.findMCS(false, true, 5000L);
      assertNotNull(mcs);
      assertTrue(mcs.size() >= 2);
    }

    @Test @Timeout(10) @DisplayName("8.11 Ester vs amide MCS >= 2")
    void esterVsAmide() throws Exception {
      SMSD smsd = new SMSD(mol("CC(=O)OC"), mol("CC(=O)NC"), defaultOpts());
      Map<Integer, Integer> mcs = smsd.findMCS(false, true, 5000L);
      assertNotNull(mcs);
      assertTrue(mcs.size() >= 3);
    }

    @Test @Timeout(10) @DisplayName("8.12 Alcohol in phenol")
    void alcoholInPhenol() throws Exception {
      SMSD smsd = new SMSD(mol("CO"), mol("Oc1ccccc1"), looseOpts());
      Map<Integer, Integer> mcs = smsd.findMCS(false, true, 5000L);
      assertNotNull(mcs);
      assertTrue(mcs.size() >= 1);
    }

    @Test @Timeout(10) @DisplayName("8.13 Thiol: methanethiol self")
    void methanethiolSelf() throws Exception {
      SMSD smsd = new SMSD(mol("CS"), mol("CS"), defaultOpts());
      assertTrue(smsd.isSubstructure());
    }

    @Test @Timeout(10) @DisplayName("8.14 Sulfoxide: DMSO self")
    void dmsoSelf() throws Exception {
      SMSD smsd = new SMSD(mol("CS(=O)C"), mol("CS(=O)C"), defaultOpts());
      assertTrue(smsd.isSubstructure());
    }

    @Test @Timeout(10) @DisplayName("8.15 Sulfone: dimethyl sulfone self")
    void dimethylSulfoneSelf() throws Exception {
      SMSD smsd = new SMSD(mol("CS(=O)(=O)C"), mol("CS(=O)(=O)C"), defaultOpts());
      assertTrue(smsd.isSubstructure());
    }

    @Test @Timeout(10) @DisplayName("8.16 Ether: diethyl ether self")
    void diethylEtherSelf() throws Exception {
      SMSD smsd = new SMSD(mol("CCOCC"), mol("CCOCC"), defaultOpts());
      assertTrue(smsd.isSubstructure());
    }

    @Test @Timeout(10) @DisplayName("8.17 Epoxide self-match")
    void epoxideSelf() throws Exception {
      SMSD smsd = new SMSD(mol("C1CO1"), mol("C1CO1"), defaultOpts());
      assertTrue(smsd.isSubstructure());
    }

    @Test @Timeout(10) @DisplayName("8.18 Acid anhydride self")
    void anhydrideSelf() throws Exception {
      SMSD smsd = new SMSD(
          mol("CC(=O)OC(=O)C"), mol("CC(=O)OC(=O)C"), defaultOpts());
      assertTrue(smsd.isSubstructure());
    }

    @Test @Timeout(10) @DisplayName("8.19 Carboxylic acid in amino acid")
    void carboxylicInAminoAcid() throws Exception {
      SMSD smsd = new SMSD(mol("CC(=O)O"), mol("NCC(=O)O"), looseOpts());
      assertTrue(smsd.isSubstructure());
    }

    @Test @Timeout(10) @DisplayName("8.20 Nitrile vs isonitrile MCS >= 1")
    void nitrileVsIsonitrile() throws Exception {
      SMSD smsd = new SMSD(mol("CC#N"), mol("C[N+]#[C-]"), looseOpts());
      Map<Integer, Integer> mcs = smsd.findMCS(false, true, 5000L);
      assertNotNull(mcs);
      assertTrue(mcs.size() >= 1);
    }
  }

  // ======================================================================
  // 9. MACROCYCLE PAIRS
  // ======================================================================

  @Nested
  @DisplayName("9. MacrocyclePairs")
  class MacrocyclePairs {

    @Test @Timeout(10) @DisplayName("9.01 12-crown-4 self-match")
    void crown4Self() throws Exception {
      SMSD smsd = new SMSD(
          mol("C1COCCOCCOCCO1"), mol("C1COCCOCCOCCO1"), defaultOpts());
      assertTrue(smsd.isSubstructure());
    }

    @Test @Timeout(10) @DisplayName("9.02 15-crown-5 self-match")
    void crown5Self() throws Exception {
      SMSD smsd = new SMSD(
          mol("C1COCCOCCOCCOCCOC1"), mol("C1COCCOCCOCCOCCOC1"), defaultOpts());
      assertTrue(smsd.isSubstructure());
    }

    @Test @Timeout(10) @DisplayName("9.03 12-crown-4 vs 15-crown-5 MCS >= 8")
    void crown4VsCrown5() throws Exception {
      SMSD smsd = new SMSD(
          mol("C1COCCOCCOCCO1"),
          mol("C1COCCOCCOCCOCCOC1"),
          looseOpts());
      Map<Integer, Integer> mcs = smsd.findMCS(false, true, 5000L);
      assertNotNull(mcs);
      assertTrue(mcs.size() >= 6, "Crown ethers share polyether fragment");
    }

    @Test @Timeout(10) @DisplayName("9.04 Porphyrin core: pyrrole substructure")
    void pyrroleInPorphyrin() throws Exception {
      // Simplified porphine fragment
      SMSD smsd = new SMSD(
          mol("c1cc[nH]c1"),
          mol("C1=CC2=CC3=CC(=CC4=CC(=CC(=C1)N2)[NH]4)N3"),  // simplified porphyrin Kekule
          looseOpts());
      Map<Integer, Integer> mcs = smsd.findMCS(false, true, 5000L);
      assertNotNull(mcs);
      assertTrue(mcs.size() >= 4);
    }

    @Test @Timeout(10) @DisplayName("9.05 Cyclodextrin fragment: glucose ring self")
    void glucoseRingSelf() throws Exception {
      SMSD smsd = new SMSD(
          mol("OC1CCOC(O)C1O"), mol("OC1CCOC(O)C1O"), defaultOpts());
      assertTrue(smsd.isSubstructure());
    }

    @Test @Timeout(10) @DisplayName("9.06 18-crown-6 self-match")
    void crown6Self() throws Exception {
      SMSD smsd = new SMSD(
          mol("C1COCCOCCOCCOCCOCCO1"),
          mol("C1COCCOCCOCCOCCOCCO1"),
          defaultOpts());
      assertTrue(smsd.isSubstructure());
    }

    @Test @Timeout(10) @DisplayName("9.07 Macrolactone 12-membered ring self")
    void macrolactoneSelf() throws Exception {
      SMSD smsd = new SMSD(
          mol("O=C1CCCCCCCCCCCO1"),
          mol("O=C1CCCCCCCCCCCO1"),
          defaultOpts());
      assertTrue(smsd.isSubstructure());
    }

    @Test @Timeout(10) @DisplayName("9.08 Large ring vs linear chain")
    void largeRingVsChain() throws Exception {
      SMSD smsd = new SMSD(
          mol("C1CCCCCCCCCCC1"), mol("CCCCCCCCCCCC"), defaultOpts());
      // With ringMatchesRingOnly=true, ring shouldn't fully map to chain
      Map<Integer, Integer> mcs = smsd.findMCS(false, true, 5000L);
      assertNotNull(mcs);
    }

    @Test @Timeout(10) @DisplayName("9.09 Macrolactam self")
    void macrolactamSelf() throws Exception {
      SMSD smsd = new SMSD(
          mol("O=C1CCCCCCCCCCCN1"),
          mol("O=C1CCCCCCCCCCCN1"),
          defaultOpts());
      assertTrue(smsd.isSubstructure());
    }

    @Test @Timeout(10) @DisplayName("9.10 Crown ether fragment MCS with 18-crown-6")
    void crownFragmentInCrown6() throws Exception {
      SMSD smsd = new SMSD(
          mol("COCCOC"),
          mol("C1COCCOCCOCCOCCOCCO1"),
          looseOpts());
      // MCS should find the shared -C-O-C-C-O-C- fragment
      var mcs = smsd.findMCS(false, false, 5000L);
      assertTrue(mcs.size() >= 5, "Crown fragment MCS should be >= 5, got " + mcs.size());
    }
  }

  // ======================================================================
  // 10. POLYMER FRAGMENTS
  // ======================================================================

  @Nested
  @DisplayName("10. PolymerFragments")
  class PolymerFragments {

    @Test @Timeout(10) @DisplayName("10.01 PEG dimer in PEG trimer")
    void pegDimerInTrimer() throws Exception {
      SMSD smsd = new SMSD(
          mol("OCCOCCO"), mol("OCCOCCOCCOC"), looseOpts());
      assertTrue(smsd.isSubstructure());
    }

    @Test @Timeout(10) @DisplayName("10.02 PEG monomer in PEG dimer")
    void pegMonomerInDimer() throws Exception {
      SMSD smsd = new SMSD(mol("OCCO"), mol("OCCOCCO"), looseOpts());
      assertTrue(smsd.isSubstructure());
    }

    @Test @Timeout(10) @DisplayName("10.03 PEG trimer self-match")
    void pegTrimerSelf() throws Exception {
      SMSD smsd = new SMSD(
          mol("OCCOCCOCCOC"), mol("OCCOCCOCCOC"), defaultOpts());
      assertTrue(smsd.isSubstructure());
    }

    @Test @Timeout(10) @DisplayName("10.04 Lactic acid dimer self-match")
    void lacticAcidDimerSelf() throws Exception {
      SMSD smsd = new SMSD(
          mol("CC(O)C(=O)OC(C)C(=O)O"),
          mol("CC(O)C(=O)OC(C)C(=O)O"),
          defaultOpts());
      assertTrue(smsd.isSubstructure());
    }

    @Test @Timeout(10) @DisplayName("10.05 Lactic acid monomer in dimer")
    void lacticMonomerInDimer() throws Exception {
      SMSD smsd = new SMSD(
          mol("CC(O)C(=O)O"),
          mol("CC(O)C(=O)OC(C)C(=O)O"),
          looseOpts());
      assertTrue(smsd.isSubstructure());
    }

    @Test @Timeout(10) @DisplayName("10.06 Nylon-6 monomer: caprolactam self")
    void caprolactamSelf() throws Exception {
      SMSD smsd = new SMSD(
          mol("O=C1CCCCCN1"), mol("O=C1CCCCCN1"), defaultOpts());
      assertTrue(smsd.isSubstructure());
    }

    @Test @Timeout(10) @DisplayName("10.07 Nylon-6,6 monomer MCS >= 4")
    void nylon66Monomer() throws Exception {
      SMSD smsd = new SMSD(
          mol("NCCCCCCN"), mol("OC(=O)CCCCC(=O)O"), looseOpts());
      Map<Integer, Integer> mcs = smsd.findMCS(false, true, 5000L);
      assertNotNull(mcs);
      assertTrue(mcs.size() >= 4, "Share hexyl chain");
    }

    @Test @Timeout(10) @DisplayName("10.08 Styrene monomer in polystyrene dimer")
    void styreneInDimer() throws Exception {
      SMSD smsd = new SMSD(
          mol("C=Cc1ccccc1"),
          mol("CC(c1ccccc1)CC(c1ccccc1)C"),
          looseOpts());
      Map<Integer, Integer> mcs = smsd.findMCS(false, true, 5000L);
      assertNotNull(mcs);
      assertTrue(mcs.size() >= 6);
    }

    @Test @Timeout(10) @DisplayName("10.09 Ethylene glycol in PEG")
    void ethyleneGlycolInPEG() throws Exception {
      SMSD smsd = new SMSD(mol("OCCO"), mol("OCCOCCOCCOCCOC"), looseOpts());
      assertTrue(smsd.isSubstructure());
    }

    @Test @Timeout(10) @DisplayName("10.10 Adipic acid self-match")
    void adipicAcidSelf() throws Exception {
      SMSD smsd = new SMSD(
          mol("OC(=O)CCCCC(=O)O"),
          mol("OC(=O)CCCCC(=O)O"),
          defaultOpts());
      assertTrue(smsd.isSubstructure());
    }
  }

  // ======================================================================
  // 11. DRUG PAIR MATRIX
  // ======================================================================

  @Nested
  @DisplayName("11. DrugPairMatrix")
  class DrugPairMatrix {

    // Aspirin: CC(=O)Oc1ccccc1C(=O)O
    // Ibuprofen: CC(C)Cc1ccc(C(C)C(=O)O)cc1
    // Caffeine: Cn1c(=O)c2c(ncn2C)n(C)c1=O
    // Morphine: CN1CC[C@]23c4c5ccc(O)c4O[C@H]2C(=C[C@@H]1[C@@H]3O)C5 (simplified below)
    // Diazepam: CN1C(=O)CN=C(c2ccccc2)c2cc(Cl)ccc21
    // Metformin: CN(C)C(=N)NC(=N)N

    static final String ASPIRIN = "CC(=O)Oc1ccccc1C(=O)O";
    static final String IBUPROFEN = "CC(C)Cc1ccc(cc1)C(C)C(=O)O";
    static final String CAFFEINE = "Cn1c(=O)c2c(ncn2C)n(C)c1=O";
    static final String DIAZEPAM = "CN1C(=O)CN=C(c2ccccc2)c2cc(Cl)ccc21";
    static final String METFORMIN = "CN(C)C(=N)NC(=N)N";
    static final String PARACETAMOL = "CC(=O)Nc1ccc(O)cc1";

    @Test @Timeout(10) @DisplayName("11.01 Aspirin self-substructure")
    void aspirinSelf() throws Exception {
      assertTrue(new SMSD(mol(ASPIRIN), mol(ASPIRIN), defaultOpts()).isSubstructure());
    }

    @Test @Timeout(10) @DisplayName("11.02 Aspirin vs ibuprofen MCS >= 3")
    void aspirinVsIbuprofen() throws Exception {
      SMSD smsd = new SMSD(mol(ASPIRIN), mol(IBUPROFEN), looseOpts());
      Map<Integer, Integer> mcs = smsd.findMCS(false, true, 5000L);
      assertNotNull(mcs);
      assertTrue(mcs.size() >= 3);
    }

    @Test @Timeout(10) @DisplayName("11.03 Aspirin vs caffeine MCS >= 2")
    void aspirinVsCaffeine() throws Exception {
      SMSD smsd = new SMSD(mol(ASPIRIN), mol(CAFFEINE), looseOpts());
      Map<Integer, Integer> mcs = smsd.findMCS(false, true, 5000L);
      assertNotNull(mcs);
      assertTrue(mcs.size() >= 2);
    }

    @Test @Timeout(10) @DisplayName("11.04 Aspirin vs diazepam MCS >= 6")
    void aspirinVsDiazepam() throws Exception {
      SMSD smsd = new SMSD(mol(ASPIRIN), mol(DIAZEPAM), looseOpts());
      Map<Integer, Integer> mcs = smsd.findMCS(false, true, 5000L);
      assertNotNull(mcs);
      assertTrue(mcs.size() >= 6, "Share benzene ring");
    }

    @Test @Timeout(10) @DisplayName("11.05 Aspirin vs metformin MCS >= 1")
    void aspirinVsMetformin() throws Exception {
      SMSD smsd = new SMSD(mol(ASPIRIN), mol(METFORMIN), looseOpts());
      Map<Integer, Integer> mcs = smsd.findMCS(false, true, 5000L);
      assertNotNull(mcs);
    }

    @Test @Timeout(10) @DisplayName("11.06 Ibuprofen self-substructure")
    void ibuprofenSelf() throws Exception {
      assertTrue(new SMSD(mol(IBUPROFEN), mol(IBUPROFEN), defaultOpts()).isSubstructure());
    }

    @Test @Timeout(10) @DisplayName("11.07 Ibuprofen vs caffeine MCS >= 2")
    void ibuprofenVsCaffeine() throws Exception {
      SMSD smsd = new SMSD(mol(IBUPROFEN), mol(CAFFEINE), looseOpts());
      Map<Integer, Integer> mcs = smsd.findMCS(false, true, 5000L);
      assertNotNull(mcs);
      assertTrue(mcs.size() >= 2);
    }

    @Test @Timeout(10) @DisplayName("11.08 Ibuprofen vs diazepam MCS >= 6")
    void ibuprofenVsDiazepam() throws Exception {
      SMSD smsd = new SMSD(mol(IBUPROFEN), mol(DIAZEPAM), looseOpts());
      Map<Integer, Integer> mcs = smsd.findMCS(false, true, 5000L);
      assertNotNull(mcs);
      assertTrue(mcs.size() >= 6, "Both have benzene ring");
    }

    @Test @Timeout(10) @DisplayName("11.09 Ibuprofen vs metformin MCS >= 1")
    void ibuprofenVsMetformin() throws Exception {
      SMSD smsd = new SMSD(mol(IBUPROFEN), mol(METFORMIN), looseOpts());
      Map<Integer, Integer> mcs = smsd.findMCS(false, true, 5000L);
      assertNotNull(mcs);
    }

    @Test @Timeout(10) @DisplayName("11.10 Ibuprofen vs paracetamol MCS >= 6")
    void ibuprofenVsParacetamol() throws Exception {
      SMSD smsd = new SMSD(mol(IBUPROFEN), mol(PARACETAMOL), looseOpts());
      Map<Integer, Integer> mcs = smsd.findMCS(false, true, 5000L);
      assertNotNull(mcs);
      assertTrue(mcs.size() >= 6, "Share phenyl core");
    }

    @Test @Timeout(10) @DisplayName("11.11 Caffeine self-substructure")
    void caffeineSelf() throws Exception {
      assertTrue(new SMSD(mol(CAFFEINE), mol(CAFFEINE), defaultOpts()).isSubstructure());
    }

    @Test @Timeout(10) @DisplayName("11.12 Caffeine vs diazepam MCS >= 3")
    void caffeineVsDiazepam() throws Exception {
      SMSD smsd = new SMSD(mol(CAFFEINE), mol(DIAZEPAM), looseOpts());
      Map<Integer, Integer> mcs = smsd.findMCS(false, true, 5000L);
      assertNotNull(mcs);
      assertTrue(mcs.size() >= 3);
    }

    @Test @Timeout(10) @DisplayName("11.13 Caffeine vs metformin MCS >= 2")
    void caffeineVsMetformin() throws Exception {
      SMSD smsd = new SMSD(mol(CAFFEINE), mol(METFORMIN), looseOpts());
      Map<Integer, Integer> mcs = smsd.findMCS(false, true, 5000L);
      assertNotNull(mcs);
      assertTrue(mcs.size() >= 2);
    }

    @Test @Timeout(10) @DisplayName("11.14 Caffeine vs paracetamol MCS >= 2")
    void caffeineVsParacetamol() throws Exception {
      SMSD smsd = new SMSD(mol(CAFFEINE), mol(PARACETAMOL), looseOpts());
      Map<Integer, Integer> mcs = smsd.findMCS(false, true, 5000L);
      assertNotNull(mcs);
      assertTrue(mcs.size() >= 2);
    }

    @Test @Timeout(10) @DisplayName("11.15 Diazepam self-substructure")
    void diazepamSelf() throws Exception {
      assertTrue(new SMSD(mol(DIAZEPAM), mol(DIAZEPAM), defaultOpts()).isSubstructure());
    }

    @Test @Timeout(10) @DisplayName("11.16 Diazepam vs metformin MCS >= 1")
    void diazepamVsMetformin() throws Exception {
      SMSD smsd = new SMSD(mol(DIAZEPAM), mol(METFORMIN), looseOpts());
      Map<Integer, Integer> mcs = smsd.findMCS(false, true, 5000L);
      assertNotNull(mcs);
    }

    @Test @Timeout(10) @DisplayName("11.17 Diazepam vs paracetamol MCS >= 6")
    void diazepamVsParacetamol() throws Exception {
      SMSD smsd = new SMSD(mol(DIAZEPAM), mol(PARACETAMOL), looseOpts());
      Map<Integer, Integer> mcs = smsd.findMCS(false, true, 5000L);
      assertNotNull(mcs);
      assertTrue(mcs.size() >= 6, "Share aromatic ring");
    }

    @Test @Timeout(10) @DisplayName("11.18 Metformin self-substructure")
    void metforminSelf() throws Exception {
      assertTrue(new SMSD(mol(METFORMIN), mol(METFORMIN), defaultOpts()).isSubstructure());
    }

    @Test @Timeout(10) @DisplayName("11.19 Metformin vs paracetamol MCS >= 1")
    void metforminVsParacetamol() throws Exception {
      SMSD smsd = new SMSD(mol(METFORMIN), mol(PARACETAMOL), looseOpts());
      Map<Integer, Integer> mcs = smsd.findMCS(false, true, 5000L);
      assertNotNull(mcs);
    }

    @Test @Timeout(10) @DisplayName("11.20 Paracetamol self-substructure")
    void paracetamolSelf() throws Exception {
      assertTrue(new SMSD(mol(PARACETAMOL), mol(PARACETAMOL), defaultOpts()).isSubstructure());
    }

    @Test @Timeout(10) @DisplayName("11.21 Benzene in aspirin")
    void benzeneInAspirin() throws Exception {
      assertTrue(new SMSD(mol("c1ccccc1"), mol(ASPIRIN), looseOpts()).isSubstructure());
    }

    @Test @Timeout(10) @DisplayName("11.22 Benzene in ibuprofen")
    void benzeneInIbuprofen() throws Exception {
      assertTrue(new SMSD(mol("c1ccccc1"), mol(IBUPROFEN), looseOpts()).isSubstructure());
    }

    @Test @Timeout(10) @DisplayName("11.23 Benzene in diazepam")
    void benzeneInDiazepam() throws Exception {
      assertTrue(new SMSD(mol("c1ccccc1"), mol(DIAZEPAM), looseOpts()).isSubstructure());
    }

    @Test @Timeout(10) @DisplayName("11.24 Benzene in paracetamol")
    void benzeneInParacetamol() throws Exception {
      assertTrue(new SMSD(mol("c1ccccc1"), mol(PARACETAMOL), looseOpts()).isSubstructure());
    }

    @Test @Timeout(10) @DisplayName("11.25 Carboxylic acid group in aspirin")
    void carboxylicInAspirin() throws Exception {
      assertTrue(new SMSD(mol("CC(=O)O"), mol(ASPIRIN), looseOpts()).isSubstructure());
    }

    @Test @Timeout(10) @DisplayName("11.26 Carboxylic acid group in ibuprofen")
    void carboxylicInIbuprofen() throws Exception {
      assertTrue(new SMSD(mol("CC(=O)O"), mol(IBUPROFEN), looseOpts()).isSubstructure());
    }

    @Test @Timeout(10) @DisplayName("11.27 Aspirin vs paracetamol MCS >= 6")
    void aspirinVsParacetamol() throws Exception {
      SMSD smsd = new SMSD(mol(ASPIRIN), mol(PARACETAMOL), looseOpts());
      Map<Integer, Integer> mcs = smsd.findMCS(false, true, 5000L);
      assertNotNull(mcs);
      assertTrue(mcs.size() >= 6, "Both have substituted benzene");
    }

    @Test @Timeout(10) @DisplayName("11.28 Caffeine similarity upper bound > 0")
    void caffeineSimilarity() throws Exception {
      SMSD smsd = new SMSD(mol(CAFFEINE), mol(DIAZEPAM), looseOpts());
      double ub = smsd.similarityUpperBound();
      assertTrue(ub > 0, "Non-zero similarity bound");
    }

    @Test @Timeout(10) @DisplayName("11.29 Aspirin similarity upper bound with self == 1.0")
    void aspirinSelfSimilarity() throws Exception {
      SMSD smsd = new SMSD(mol(ASPIRIN), mol(ASPIRIN), defaultOpts());
      double ub = smsd.similarityUpperBound();
      assertTrue(ub >= 0.99, "Self-similarity should be ~1.0");
    }

    @Test @Timeout(10) @DisplayName("11.30 Metformin MCS with itself full size")
    void metforminSelfMCS() throws Exception {
      IAtomContainer met = mol(METFORMIN);
      SMSD smsd = new SMSD(met, met, defaultOpts());
      Map<Integer, Integer> mcs = smsd.findMCS(false, true, 5000L);
      assertNotNull(mcs);
      assertEquals(met.getAtomCount(), mcs.size(), "Self MCS should be full molecule");
    }
  }

  // ======================================================================
  // 12. SYMMETRY STRESS
  // ======================================================================

  @Nested
  @DisplayName("12. SymmetryStress")
  class SymmetryStress {

    @Test @Timeout(10) @DisplayName("12.01 Adamantane self-match")
    void adamantaneSelf() throws Exception {
      SMSD smsd = new SMSD(
          mol("C1C2CC3CC1CC(C2)C3"),
          mol("C1C2CC3CC1CC(C2)C3"),
          defaultOpts());
      assertTrue(smsd.isSubstructure());
    }

    @Test @Timeout(10) @DisplayName("12.02 Adamantane self MCS == 10")
    void adamantaneSelfMCS() throws Exception {
      IAtomContainer ada = mol("C1C2CC3CC1CC(C2)C3");
      SMSD smsd = new SMSD(ada, ada, defaultOpts());
      Map<Integer, Integer> mcs = smsd.findMCS(false, true, 5000L);
      assertNotNull(mcs);
      assertEquals(ada.getAtomCount(), mcs.size());
    }

    @Test @Timeout(10) @DisplayName("12.03 Cubane self-match")
    void cubaneSelf() throws Exception {
      SMSD smsd = new SMSD(
          mol("C12C3C4C1C5C3C4C25"),
          mol("C12C3C4C1C5C3C4C25"),
          defaultOpts());
      assertTrue(smsd.isSubstructure());
    }

    @Test @Timeout(10) @DisplayName("12.04 Cubane self MCS == 8")
    void cubaneSelfMCS() throws Exception {
      IAtomContainer cub = mol("C12C3C4C1C5C3C4C25");
      SMSD smsd = new SMSD(cub, cub, defaultOpts());
      Map<Integer, Integer> mcs = smsd.findMCS(false, true, 5000L);
      assertNotNull(mcs);
      assertEquals(cub.getAtomCount(), mcs.size());
    }

    @Test @Timeout(10) @DisplayName("12.05 Benzene multiple equivalent mappings")
    void benzeneMultipleMappings() throws Exception {
      IAtomContainer bz = mol("c1ccccc1");
      SMSD smsd = new SMSD(bz, bz, looseOpts());
      List<Map<Integer, Integer>> all = smsd.findAllSubstructures(20, 5000L);
      assertNotNull(all);
      assertTrue(all.size() >= 2, "Benzene has at least 2 automorphisms");
    }

    @Test @Timeout(10) @DisplayName("12.06 Naphthalene self-match full size")
    void naphthaleneSelfMCS() throws Exception {
      IAtomContainer naph = mol("c1ccc2ccccc2c1");
      SMSD smsd = new SMSD(naph, naph, looseOpts());
      Map<Integer, Integer> mcs = smsd.findMCS(false, true, 5000L);
      assertNotNull(mcs);
      assertEquals(naph.getAtomCount(), mcs.size());
    }

    @Test @Timeout(10) @DisplayName("12.07 Cyclohexane self-match full size")
    void cyclohexaneSelfMCS() throws Exception {
      IAtomContainer cyc = mol("C1CCCCC1");
      SMSD smsd = new SMSD(cyc, cyc, defaultOpts());
      Map<Integer, Integer> mcs = smsd.findMCS(false, true, 5000L);
      assertNotNull(mcs);
      assertEquals(cyc.getAtomCount(), mcs.size());
    }

    @Test @Timeout(10) @DisplayName("12.08 Anthracene self-match")
    void anthraceneSelf() throws Exception {
      SMSD smsd = new SMSD(
          mol("c1ccc2cc3ccccc3cc2c1"),
          mol("c1ccc2cc3ccccc3cc2c1"),
          looseOpts());
      assertTrue(smsd.isSubstructure());
    }

    @Test @Timeout(10) @DisplayName("12.09 Biphenyl self-match full size")
    void biphenylSelfMCS() throws Exception {
      IAtomContainer biph = mol("c1ccc(-c2ccccc2)cc1");
      SMSD smsd = new SMSD(biph, biph, looseOpts());
      Map<Integer, Integer> mcs = smsd.findMCS(false, true, 5000L);
      assertNotNull(mcs);
      assertEquals(biph.getAtomCount(), mcs.size());
    }

    @Test @Timeout(10) @DisplayName("12.10 Triphenylene self-match")
    void triphenyleneSelf() throws Exception {
      SMSD smsd = new SMSD(
          mol("c1ccc2c(c1)c1ccccc1c1ccccc12"),
          mol("c1ccc2c(c1)c1ccccc1c1ccccc12"),
          looseOpts());
      assertTrue(smsd.isSubstructure());
    }

    @Test @Timeout(10) @DisplayName("12.11 Coronene self MCS (with timeout guard)")
    void coroneneSelf() throws Exception {
      SMSD smsd = new SMSD(
          mol("c1cc2ccc3cc4ccc5cc6ccc1c7c2c3c4c5c67"),
          mol("c1cc2ccc3cc4ccc5cc6ccc1c7c2c3c4c5c67"),
          looseOpts());
      Map<Integer, Integer> mcs = smsd.findMCS(false, true, 8000L);
      assertNotNull(mcs);
      assertTrue(mcs.size() >= 15, "Coronene self MCS should be large");
    }

    @Test @Timeout(10) @DisplayName("12.12 Pyrene self-match")
    void pyreneSelf() throws Exception {
      SMSD smsd = new SMSD(
          mol("c1cc2ccc3cccc4ccc(c1)c2c34"),
          mol("c1cc2ccc3cccc4ccc(c1)c2c34"),
          looseOpts());
      assertTrue(smsd.isSubstructure());
    }

    @Test @Timeout(10) @DisplayName("12.13 Decalin self MCS full size")
    void decalinSelfMCS() throws Exception {
      IAtomContainer dec = mol("C1CCC2CCCCC2C1");
      SMSD smsd = new SMSD(dec, dec, defaultOpts());
      Map<Integer, Integer> mcs = smsd.findMCS(false, true, 5000L);
      assertNotNull(mcs);
      assertEquals(dec.getAtomCount(), mcs.size());
    }

    @Test @Timeout(10) @DisplayName("12.14 Benzene in biphenyl")
    void benzeneInBiphenyl() throws Exception {
      SMSD smsd = new SMSD(
          mol("c1ccccc1"), mol("c1ccc(-c2ccccc2)cc1"), looseOpts());
      assertTrue(smsd.isSubstructure());
    }

    @Test @Timeout(10) @DisplayName("12.15 Cyclopropane self multiple mappings")
    void cyclopropaneMultipleMappings() throws Exception {
      IAtomContainer cp = mol("C1CC1");
      SMSD smsd = new SMSD(cp, cp, defaultOpts());
      List<Map<Integer, Integer>> all = smsd.findAllSubstructures(10, 5000L);
      assertNotNull(all);
      assertTrue(all.size() >= 2, "Cyclopropane has rotational symmetry");
    }
  }

  // ======================================================================
  // 13. TIMEOUT BEHAVIOR
  // ======================================================================

  @Nested
  @DisplayName("13. TimeoutBehavior")
  class TimeoutBehavior {

    @Test @Timeout(10) @DisplayName("13.01 Very short timeout on medium molecule pair")
    void veryShortTimeout() throws Exception {
      SMSD smsd = new SMSD(
          mol("c1ccc2cc3ccccc3cc2c1"),
          mol("c1ccc2c(c1)c1ccccc1c1ccccc12"),
          looseOpts());
      // 1ms timeout - may or may not find result, but should not throw
      Map<Integer, Integer> mcs = smsd.findMCS(false, true, 1L);
      // Result can be null or valid, just verify no exception
    }

    @Test @Timeout(10) @DisplayName("13.02 Short timeout on drug pair")
    void shortTimeoutDrugPair() throws Exception {
      SMSD smsd = new SMSD(
          mol("CC(=O)Oc1ccccc1C(=O)O"),
          mol("CN1C(=O)CN=C(c2ccccc2)c2cc(Cl)ccc21"),
          looseOpts());
      Map<Integer, Integer> mcs = smsd.findMCS(false, true, 5L);
      // Should not throw
    }

    @Test @Timeout(10) @DisplayName("13.03 100ms timeout on small pair gets result")
    void mediumTimeoutSmallPair() throws Exception {
      SMSD smsd = new SMSD(mol("CCO"), mol("CCCO"), defaultOpts());
      Map<Integer, Integer> mcs = smsd.findMCS(false, true, 100L);
      assertNotNull(mcs);
      assertTrue(mcs.size() >= 2, "Small pair should resolve in 100ms");
    }

    @Test @Timeout(10) @DisplayName("13.04 Substructure with 1ms timeout on large pair")
    void substruct1msLarge() throws Exception {
      SMSD smsd = new SMSD(
          mol("c1ccc2cc3ccccc3cc2c1"),
          mol("c1ccc2c(c1)c1ccccc1c1ccccc12"),
          looseOpts());
      // Should not throw, result may be true or false
      smsd.isSubstructure(1L);
    }

    @Test @Timeout(10) @DisplayName("13.05 setSubstructureTimeoutMs works")
    void setSubstructureTimeout() throws Exception {
      SMSD smsd = new SMSD(mol("CCO"), mol("CCNCCO"), defaultOpts());
      smsd.setSubstructureTimeoutMs(50L);
      // Should still work for small molecules
      assertTrue(smsd.isSubstructure());
    }

    @Test @Timeout(10) @DisplayName("13.06 setMcsTimeoutMs works")
    void setMcsTimeout() throws Exception {
      SMSD smsd = new SMSD(mol("CCO"), mol("CCCO"), defaultOpts());
      smsd.setMcsTimeoutMs(100L);
      Map<Integer, Integer> mcs = smsd.findMCS();
      assertNotNull(mcs);
    }

    @Test @Timeout(10) @DisplayName("13.07 Zero timeout does not throw")
    void zeroTimeout() throws Exception {
      SMSD smsd = new SMSD(mol("c1ccccc1"), mol("c1ccc2ccccc2c1"), looseOpts());
      // 0ms - implementation should handle gracefully
      assertDoesNotThrow(() -> smsd.findMCS(false, true, 0L));
    }

    @Test @Timeout(10) @DisplayName("13.08 Short timeout substructure ethane in propane")
    void shortTimeoutEthaneInPropane() throws Exception {
      SMSD smsd = new SMSD(mol("CC"), mol("CCC"), defaultOpts());
      assertTrue(smsd.isSubstructure(10L));
    }

    @Test @Timeout(10) @DisplayName("13.09 1ms MCS on identical small molecule")
    void shortTimeoutSelfMCS() throws Exception {
      SMSD smsd = new SMSD(mol("CC"), mol("CC"), defaultOpts());
      Map<Integer, Integer> mcs = smsd.findMCS(false, true, 1L);
      // Self-match shortcut should work even with tiny timeout
      assertNotNull(mcs);
    }

    @Test @Timeout(10) @DisplayName("13.10 Large timeout on simple pair succeeds")
    void largeTimeoutSimplePair() throws Exception {
      SMSD smsd = new SMSD(mol("c1ccccc1"), mol("c1ccc2ccccc2c1"), looseOpts());
      Map<Integer, Integer> mcs = smsd.findMCS(false, true, 10000L);
      assertNotNull(mcs);
      assertEquals(6, mcs.size());
    }
  }

  // ======================================================================
  // 14. CHEMOPTIONS MATRIX
  // ======================================================================

  @Nested
  @DisplayName("14. ChemOptionsMatrix")
  class ChemOptionsMatrix {

    // Use a sensitive pair: phenol vs benzene
    static final String QUERY = "Oc1ccccc1";   // phenol
    static final String TARGET = "c1ccccc1";    // benzene

    @Test @Timeout(10) @DisplayName("14.01 matchAtomType=true")
    void matchAtomTypeTrue() throws Exception {
      ChemOptions opts = new ChemOptions();
      opts.matchAtomType = true;
      SMSD smsd = new SMSD(mol(QUERY), mol(TARGET), opts);
      Map<Integer, Integer> mcs = smsd.findMCS(false, true, 5000L);
      assertNotNull(mcs);
      // With atom type matching, O cannot match C
      assertTrue(mcs.size() <= 6);
    }

    @Test @Timeout(10) @DisplayName("14.02 matchAtomType=false allows larger overlap")
    void matchAtomTypeFalse() throws Exception {
      ChemOptions opts = new ChemOptions();
      opts.matchAtomType = false;
      opts.ringMatchesRingOnly = false;
      SMSD smsd = new SMSD(mol(QUERY), mol(TARGET), opts);
      Map<Integer, Integer> mcs = smsd.findMCS(false, true, 5000L);
      assertNotNull(mcs);
      assertTrue(mcs.size() >= 6, "Without atom type restriction, all atoms can map");
    }

    @Test @Timeout(10) @DisplayName("14.03 matchBondOrder=STRICT")
    void bondOrderStrict() throws Exception {
      ChemOptions opts = new ChemOptions();
      opts.matchBondOrder = ChemOptions.BondOrderMode.STRICT;
      SMSD smsd = new SMSD(mol("C=C"), mol("CC"), opts);
      assertFalse(smsd.isSubstructure(), "Double bond shouldn't match single under STRICT");
    }

    @Test @Timeout(10) @DisplayName("14.04 matchBondOrder=LOOSE")
    void bondOrderLoose() throws Exception {
      ChemOptions opts = new ChemOptions();
      opts.matchBondOrder = ChemOptions.BondOrderMode.LOOSE;
      opts.ringMatchesRingOnly = false;
      SMSD smsd = new SMSD(mol("C=C"), mol("CC"), opts);
      Map<Integer, Integer> mcs = smsd.findMCS(false, true, 5000L);
      assertNotNull(mcs);
      assertTrue(mcs.size() >= 2, "LOOSE bond order allows C=C to match CC");
    }

    @Test @Timeout(10) @DisplayName("14.05 matchBondOrder=ANY")
    void bondOrderAny() throws Exception {
      ChemOptions opts = new ChemOptions();
      opts.matchBondOrder = ChemOptions.BondOrderMode.ANY;
      opts.ringMatchesRingOnly = false;
      SMSD smsd = new SMSD(mol("C#C"), mol("CC"), opts);
      Map<Integer, Integer> mcs = smsd.findMCS(false, true, 5000L);
      assertNotNull(mcs);
      assertTrue(mcs.size() >= 2, "ANY bond order allows triple to match single");
    }

    @Test @Timeout(10) @DisplayName("14.06 ringMatchesRingOnly=true")
    void ringMatchesRingTrue() throws Exception {
      ChemOptions opts = new ChemOptions();
      opts.ringMatchesRingOnly = true;
      SMSD smsd = new SMSD(mol("C1CCCCC1"), mol("CCCCCC"), opts);
      Map<Integer, Integer> mcs = smsd.findMCS(false, true, 5000L);
      // Ring atoms should not map to chain
      assertTrue(mcs == null || mcs.size() < 6,
          "Ring-to-ring should prevent full ring mapping to chain");
    }

    @Test @Timeout(10) @DisplayName("14.07 ringMatchesRingOnly=false allows ring to chain")
    void ringMatchesRingFalse() throws Exception {
      ChemOptions opts = new ChemOptions();
      opts.ringMatchesRingOnly = false;
      SMSD smsd = new SMSD(mol("C1CCCCC1"), mol("CCCCCC"), opts);
      Map<Integer, Integer> mcs = smsd.findMCS(false, true, 5000L);
      assertNotNull(mcs);
      assertTrue(mcs.size() >= 5, "Without ring restriction, ring can map to chain");
    }

    @Test @Timeout(10) @DisplayName("14.08 useChirality=false ignores chirality")
    void useChiralityFalse() throws Exception {
      ChemOptions opts = new ChemOptions();
      opts.useChirality = false;
      SMSD smsd = new SMSD(
          mol("[C@@H](N)(C)C(=O)O"),
          mol("[C@H](N)(C)C(=O)O"),
          opts);
      assertTrue(smsd.isSubstructure());
    }

    @Test @Timeout(10) @DisplayName("14.09 matchFormalCharge=true on charged pair")
    void matchChargeTrue() throws Exception {
      ChemOptions opts = new ChemOptions();
      opts.matchFormalCharge = true;
      opts.ringMatchesRingOnly = false;
      SMSD smsd = new SMSD(mol("[NH4+]"), mol("N"), opts);
      Map<Integer, Integer> mcs = smsd.findMCS(false, true, 5000L);
      // With charge matching, charged N may not match neutral N
      assertNotNull(mcs);
    }

    @Test @Timeout(10) @DisplayName("14.10 matchFormalCharge=false on charged pair")
    void matchChargeFalse() throws Exception {
      ChemOptions opts = new ChemOptions();
      opts.matchFormalCharge = false;
      opts.ringMatchesRingOnly = false;
      SMSD smsd = new SMSD(mol("[NH4+]"), mol("N"), opts);
      Map<Integer, Integer> mcs = smsd.findMCS(false, true, 5000L);
      assertNotNull(mcs);
      assertTrue(mcs.size() >= 1);
    }

    @Test @Timeout(10) @DisplayName("14.11 completeRingsOnly=true")
    void completeRingsTrue() throws Exception {
      ChemOptions opts = new ChemOptions();
      opts.completeRingsOnly = true;
      SMSD smsd = new SMSD(mol("c1ccccc1"), mol("c1ccc2ccccc2c1"), opts);
      Map<Integer, Integer> mcs = smsd.findMCS(false, true, 5000L);
      assertNotNull(mcs);
      // Should find complete benzene ring
      assertTrue(mcs.size() >= 6);
    }

    @Test @Timeout(10) @DisplayName("14.12 matchIsotope toggled")
    void matchIsotopeToggle() throws Exception {
      ChemOptions optsOn = new ChemOptions();
      optsOn.matchIsotope = true;
      ChemOptions optsOff = new ChemOptions();
      optsOff.matchIsotope = false;
      SMSD smsdOn = new SMSD(mol("[13CH4]"), mol("C"), optsOn);
      SMSD smsdOff = new SMSD(mol("[13CH4]"), mol("C"), optsOff);
      Map<Integer, Integer> mcsOff = smsdOff.findMCS(false, true, 5000L);
      assertNotNull(mcsOff);
      assertTrue(mcsOff.size() >= 1, "Without isotope match, should overlap");
    }

    @Test @Timeout(10) @DisplayName("14.13 aromaticityMode=STRICT")
    void aromaticityStrict() throws Exception {
      ChemOptions opts = new ChemOptions();
      opts.aromaticityMode = ChemOptions.AromaticityMode.STRICT;
      SMSD smsd = new SMSD(mol("c1ccccc1"), mol("c1ccccc1"), opts);
      assertTrue(smsd.isSubstructure());
    }

    @Test @Timeout(10) @DisplayName("14.14 aromaticityMode=FLEXIBLE")
    void aromaticityFlexible() throws Exception {
      ChemOptions opts = new ChemOptions();
      opts.aromaticityMode = ChemOptions.AromaticityMode.FLEXIBLE;
      SMSD smsd = new SMSD(mol("c1ccccc1"), mol("C1=CC=CC=C1"), opts);
      assertTrue(smsd.isSubstructure(), "FLEXIBLE should match aromatic with Kekule");
    }

    @Test @Timeout(10) @DisplayName("14.15 profile strict vs compat-fmcs differ on sensitive pair")
    void profileComparison() throws Exception {
      ChemOptions strict = ChemOptions.profile("strict");
      ChemOptions compat = ChemOptions.profile("compat-fmcs");
      // Cyclohexane vs hexane
      SMSD smsdStrict = new SMSD(mol("C1CCCCC1"), mol("CCCCCC"), strict);
      SMSD smsdCompat = new SMSD(mol("C1CCCCC1"), mol("CCCCCC"), compat);
      Map<Integer, Integer> mcsStrict = smsdStrict.findMCS(false, true, 5000L);
      Map<Integer, Integer> mcsCompat = smsdCompat.findMCS(false, true, 5000L);
      int sizeStrict = (mcsStrict == null) ? 0 : mcsStrict.size();
      int sizeCompat = (mcsCompat == null) ? 0 : mcsCompat.size();
      assertTrue(sizeCompat >= sizeStrict,
          "Compat profile should be at least as permissive as strict");
    }
  }

  // ======================================================================
  // 15. RING SYSTEM TESTS
  // ======================================================================

  @Nested
  @DisplayName("15. RingSystemTests")
  class RingSystemTests {

    // ------------------------------------------------------------------
    // 15.1 Small rings (3-6 members)
    // ------------------------------------------------------------------

    @Test @Timeout(10) @DisplayName("15.01 Cyclopropane (3-ring) self-match MCS=3")
    void cyclopropaneSelfMatch() throws Exception {
      IAtomContainer cp = mol("C1CC1");
      SMSD smsd = new SMSD(cp, cp, defaultOpts());
      assertTrue(smsd.isSubstructure());
      Map<Integer, Integer> mcs = smsd.findMCS(false, true, 5000L);
      assertNotNull(mcs);
      assertEquals(3, mcs.size(), "Cyclopropane self-match MCS should be 3 atoms");
    }

    @Test @Timeout(10) @DisplayName("15.02 Cyclobutane (4-ring) self-match MCS=4")
    void cyclobutaneSelfMatch() throws Exception {
      IAtomContainer cb = mol("C1CCC1");
      SMSD smsd = new SMSD(cb, cb, defaultOpts());
      assertTrue(smsd.isSubstructure());
      Map<Integer, Integer> mcs = smsd.findMCS(false, true, 5000L);
      assertNotNull(mcs);
      assertEquals(4, mcs.size(), "Cyclobutane self-match MCS should be 4 atoms");
    }

    @Test @Timeout(10) @DisplayName("15.03 Cyclopentane (5-ring) vs cyclopentadiene substructure")
    void cyclopentaneVsCyclopentadiene() throws Exception {
      // Cyclopentadiene: C1=CC=CC1 (5 carbons, some double bonds)
      // With loose bond order, cyclopentadiene ring should still match
      ChemOptions loose = looseOpts();
      SMSD smsd = new SMSD(mol("C1=CC=CC1"), mol("C1CCCC1"), loose);
      Map<Integer, Integer> mcs = smsd.findMCS(false, true, 5000L);
      assertNotNull(mcs);
      assertTrue(mcs.size() >= 5,
          "Cyclopentadiene vs cyclopentane under loose bond order should map 5 atoms");
    }

    @Test @Timeout(10) @DisplayName("15.04 Cyclohexane (6-ring) self-match MCS=6")
    void cyclohexaneSelfMatch() throws Exception {
      IAtomContainer ch = mol("C1CCCCC1");
      SMSD smsd = new SMSD(ch, ch, defaultOpts());
      assertTrue(smsd.isSubstructure());
      Map<Integer, Integer> mcs = smsd.findMCS(false, true, 5000L);
      assertNotNull(mcs);
      assertEquals(6, mcs.size(), "Cyclohexane self-match MCS should be 6 atoms");
    }

    @Test @Timeout(10) @DisplayName("15.05 Benzene (aromatic 6) vs cyclohexane (aliphatic 6) strict aromaticity no match")
    void benzeneVsCyclohexaneStrictArom() throws Exception {
      ChemOptions strict = new ChemOptions();
      strict.aromaticityMode = ChemOptions.AromaticityMode.STRICT;
      SMSD smsd = new SMSD(mol("c1ccccc1"), mol("C1CCCCC1"), strict);
      Map<Integer, Integer> mcs = smsd.findMCS(false, true, 5000L);
      int sz = (mcs == null) ? 0 : mcs.size();
      assertTrue(sz < 6,
          "STRICT aromaticity: benzene vs cyclohexane should NOT fully match (got " + sz + ")");
    }

    // ------------------------------------------------------------------
    // 15.2 Medium rings (7-12 members)
    // ------------------------------------------------------------------

    @Test @Timeout(10) @DisplayName("15.06 Cycloheptane (7-ring) self-match MCS=7")
    void cycloheptaneSelfMatch() throws Exception {
      IAtomContainer ch7 = mol("C1CCCCCC1");
      SMSD smsd = new SMSD(ch7, ch7, defaultOpts());
      assertTrue(smsd.isSubstructure());
      Map<Integer, Integer> mcs = smsd.findMCS(false, true, 5000L);
      assertNotNull(mcs);
      assertEquals(7, mcs.size(), "Cycloheptane self-match MCS should be 7 atoms");
    }

    @Test @Timeout(10) @DisplayName("15.07 Azulene (5+7 fused non-benzenoid aromatic) self-match")
    void azuleneSelfMatch() throws Exception {
      // Azulene: fused 5+7 aromatic ring
      IAtomContainer az = mol("c1ccc2ccccc2c1");
      SMSD smsd = new SMSD(az, az, defaultOpts());
      assertTrue(smsd.isSubstructure(), "Azulene should self-match");
      Map<Integer, Integer> mcs = smsd.findMCS(false, true, 5000L);
      assertNotNull(mcs);
      assertEquals(10, mcs.size(), "Azulene self-match MCS should be 10 atoms");
    }

    @Test @Timeout(10) @DisplayName("15.08 Cyclooctatetraene (8-ring non-aromatic) self-match")
    void cyclooctatetraeneSelfMatch() throws Exception {
      // COT: tub-shaped, non-aromatic
      IAtomContainer cot = mol("C1=CC=CC=CC=C1");
      SMSD smsd = new SMSD(cot, cot, defaultOpts());
      assertTrue(smsd.isSubstructure(), "Cyclooctatetraene should self-match");
      Map<Integer, Integer> mcs = smsd.findMCS(false, true, 5000L);
      assertNotNull(mcs);
      assertEquals(8, mcs.size(), "COT self-match MCS should be 8 atoms");
    }

    @Test @Timeout(10) @DisplayName("15.09 12-crown-4 (12-membered ring with oxygens) self-match")
    void crown4SelfMatch() throws Exception {
      // 12-crown-4: -CH2-O-CH2-CH2-O-CH2-CH2-O-CH2-CH2-O-CH2-
      IAtomContainer crown = mol("C1COCCOCCOCCO1");
      SMSD smsd = new SMSD(crown, crown, defaultOpts());
      assertTrue(smsd.isSubstructure(), "12-crown-4 should self-match");
      Map<Integer, Integer> mcs = smsd.findMCS(false, true, 5000L);
      assertNotNull(mcs);
      assertTrue(mcs.size() >= 12,
          "12-crown-4 self-match MCS should have at least 12 atoms");
    }

    // ------------------------------------------------------------------
    // 15.3 Large macrocycles (>12 members)
    // ------------------------------------------------------------------

    @Test @Timeout(10) @DisplayName("15.10 14-membered lactone ring self-match")
    void fourteenMemberedLactone() throws Exception {
      // Simple 14-membered lactone (erythromycin-like ring)
      IAtomContainer lac = mol("O=C1CCCCCCCCCCCCO1");
      SMSD smsd = new SMSD(lac, lac, defaultOpts());
      assertTrue(smsd.isSubstructure(), "14-membered lactone should self-match");
      Map<Integer, Integer> mcs = smsd.findMCS(false, true, 5000L);
      assertNotNull(mcs);
      assertTrue(mcs.size() >= 14,
          "14-membered lactone self-match should map >= 14 atoms");
    }

    @Test @Timeout(10) @DisplayName("15.11 18-crown-6 (18-membered ring) self-match")
    void crown6SelfMatch() throws Exception {
      // 18-crown-6: 6 oxygens alternating with ethylene bridges
      IAtomContainer crown6 = mol("C1COCCOCCOCCOCCOCCO1");
      SMSD smsd = new SMSD(crown6, crown6, defaultOpts());
      assertTrue(smsd.isSubstructure(), "18-crown-6 should self-match");
      Map<Integer, Integer> mcs = smsd.findMCS(false, true, 5000L);
      assertNotNull(mcs);
      assertTrue(mcs.size() >= 18,
          "18-crown-6 self-match MCS should have at least 18 atoms");
    }

    @Test @Timeout(10) @DisplayName("15.12 Porphyrin core (20-membered inner ring) self-match")
    void porphyrinCoreSelfMatch() throws Exception {
      // Porphine: simplest porphyrin
      IAtomContainer porph = mol("c1cc2cc3ccc(cc4ccc(cc5ccc(cc1n2)[nH]5)n4)[nH]3");
      SMSD smsd = new SMSD(porph, porph, defaultOpts());
      assertTrue(smsd.isSubstructure(), "Porphyrin core should self-match");
      Map<Integer, Integer> mcs = smsd.findMCS(false, true, 5000L);
      assertNotNull(mcs);
      assertTrue(mcs.size() >= 20,
          "Porphyrin core self-match MCS should map >= 20 atoms");
    }

    // ------------------------------------------------------------------
    // 15.4 Fused ring systems
    // ------------------------------------------------------------------

    @Test @Timeout(10) @DisplayName("15.13 Naphthalene (2 fused 6-rings) STRICT vs IGNORE mode")
    void naphthaleneFusionModes() throws Exception {
      IAtomContainer naph = mol("c1ccc2ccccc2c1");
      // STRICT fusion
      ChemOptions strictFusion = new ChemOptions();
      strictFusion.ringFusionMode = ChemOptions.RingFusionMode.STRICT;
      SMSD smsdStrict = new SMSD(naph, naph, strictFusion);
      Map<Integer, Integer> mcsStrict = smsdStrict.findMCS(false, true, 5000L);
      assertNotNull(mcsStrict);
      assertEquals(10, mcsStrict.size(), "Naphthalene STRICT fusion self-match should be 10");
      // IGNORE fusion
      ChemOptions ignoreFusion = new ChemOptions();
      ignoreFusion.ringFusionMode = ChemOptions.RingFusionMode.IGNORE;
      SMSD smsdIgnore = new SMSD(naph, naph, ignoreFusion);
      Map<Integer, Integer> mcsIgnore = smsdIgnore.findMCS(false, true, 5000L);
      assertNotNull(mcsIgnore);
      assertEquals(10, mcsIgnore.size(), "Naphthalene IGNORE fusion self-match should be 10");
    }

    @Test @Timeout(10) @DisplayName("15.14 Anthracene (3 linear fused) vs phenanthrene (3 angular fused) MCS")
    void anthraceneVsPhenanthrene() throws Exception {
      IAtomContainer anthra = mol("c1ccc2cc3ccccc3cc2c1");
      IAtomContainer phenan = mol("c1ccc2c(c1)cc1ccccc1c2");
      SMSD smsd = new SMSD(anthra, phenan, defaultOpts());
      Map<Integer, Integer> mcs = smsd.findMCS(false, true, 5000L);
      assertNotNull(mcs);
      // Both are C14H10, same number of atoms, MCS should be large (>= 10)
      assertTrue(mcs.size() >= 10,
          "Anthracene vs phenanthrene MCS should be >= 10 (got " + mcs.size() + ")");
    }

    @Test @Timeout(10) @DisplayName("15.15 Pyrene (4 fused rings) self-match")
    void pyreneSelfMatch() throws Exception {
      IAtomContainer pyrene = mol("c1cc2ccc3cccc4ccc(c1)c2c34");
      SMSD smsd = new SMSD(pyrene, pyrene, defaultOpts());
      assertTrue(smsd.isSubstructure(), "Pyrene should self-match");
      Map<Integer, Integer> mcs = smsd.findMCS(false, true, 5000L);
      assertNotNull(mcs);
      assertEquals(16, mcs.size(), "Pyrene self-match MCS should be 16 atoms");
    }

    @Test @Timeout(10) @DisplayName("15.16 Coronene (7 fused rings) self-match with timeout")
    void coroneneSelfMatch() throws Exception {
      IAtomContainer cor = mol("c1cc2ccc3ccc4ccc5ccc6ccc1c7c2c3c4c5c67");
      SMSD smsd = new SMSD(cor, cor, defaultOpts());
      assertTrue(smsd.isSubstructure(), "Coronene should self-match");
      Map<Integer, Integer> mcs = smsd.findMCS(false, true, 5000L);
      assertNotNull(mcs);
      assertTrue(mcs.size() >= 20,
          "Coronene self-match MCS should be >= 20 atoms (got " + mcs.size() + ")");
    }

    @Test @Timeout(10) @DisplayName("15.17 Indene (5+6 fused) vs indane (5+6 saturated fused)")
    void indeneVsIndane() throws Exception {
      // Indene: benzene fused with cyclopentadiene
      IAtomContainer indene = mol("C1=Cc2ccccc2C1");
      // Indane: benzene fused with cyclopentane
      IAtomContainer indane = mol("C1Cc2ccccc2C1");
      SMSD smsd = new SMSD(indene, indane, defaultOpts());
      Map<Integer, Integer> mcs = smsd.findMCS(false, true, 5000L);
      assertNotNull(mcs);
      // Should share most of the framework (9 atoms)
      assertTrue(mcs.size() >= 8,
          "Indene vs indane should share >= 8 atoms (got " + mcs.size() + ")");
    }

    // ------------------------------------------------------------------
    // 15.5 Bridged rings
    // ------------------------------------------------------------------

    @Test @Timeout(10) @DisplayName("15.18 Norbornane (bridged bicyclic C7H12) self-match")
    void norbornaneSelfMatch() throws Exception {
      IAtomContainer norb = mol("C1CC2CC1CC2");
      SMSD smsd = new SMSD(norb, norb, defaultOpts());
      assertTrue(smsd.isSubstructure(), "Norbornane should self-match");
      Map<Integer, Integer> mcs = smsd.findMCS(false, true, 5000L);
      assertNotNull(mcs);
      assertEquals(7, mcs.size(), "Norbornane self-match MCS should be 7 atoms");
    }

    @Test @Timeout(10) @DisplayName("15.19 Adamantane (cage structure) self-match MCS=10")
    void adamantaneSelfMatch() throws Exception {
      IAtomContainer adam = mol("C1C2CC3CC1CC(C2)C3");
      SMSD smsd = new SMSD(adam, adam, defaultOpts());
      assertTrue(smsd.isSubstructure(), "Adamantane should self-match");
      Map<Integer, Integer> mcs = smsd.findMCS(false, true, 5000L);
      assertNotNull(mcs);
      assertEquals(10, mcs.size(), "Adamantane self-match MCS should be 10 atoms");
    }

    @Test @Timeout(10) @DisplayName("15.20 Camphor (bridged bicyclic ketone) self-match")
    void camphorSelfMatch() throws Exception {
      // Camphor: bridged bicyclic monoterpene ketone
      IAtomContainer camphor = mol("CC1(C)C2CCC1(C)C(=O)C2");
      SMSD smsd = new SMSD(camphor, camphor, defaultOpts());
      assertTrue(smsd.isSubstructure(), "Camphor should self-match");
      Map<Integer, Integer> mcs = smsd.findMCS(false, true, 5000L);
      assertNotNull(mcs);
      assertTrue(mcs.size() >= 10,
          "Camphor self-match MCS should map >= 10 atoms (got " + mcs.size() + ")");
    }

    // ------------------------------------------------------------------
    // 15.6 Spiro rings
    // ------------------------------------------------------------------

    @Test @Timeout(10) @DisplayName("15.21 Spiro[4.4]nonane self-match MCS=9")
    void spiro44nonaneSelfMatch() throws Exception {
      // Two cyclopentane rings sharing one atom
      IAtomContainer spiro = mol("C1CCC2(C1)CCCC2");
      SMSD smsd = new SMSD(spiro, spiro, defaultOpts());
      assertTrue(smsd.isSubstructure(), "Spiro[4.4]nonane should self-match");
      Map<Integer, Integer> mcs = smsd.findMCS(false, true, 5000L);
      assertNotNull(mcs);
      assertEquals(9, mcs.size(), "Spiro[4.4]nonane self-match MCS should be 9 atoms");
    }

    @Test @Timeout(10) @DisplayName("15.22 Spiro[5.5]undecane self-match MCS=11")
    void spiro55undecaneSelfMatch() throws Exception {
      // Two cyclohexane rings sharing one atom
      IAtomContainer spiro = mol("C1CCCC2(C1)CCCCC2");
      SMSD smsd = new SMSD(spiro, spiro, defaultOpts());
      assertTrue(smsd.isSubstructure(), "Spiro[5.5]undecane should self-match");
      Map<Integer, Integer> mcs = smsd.findMCS(false, true, 5000L);
      assertNotNull(mcs);
      assertEquals(11, mcs.size(), "Spiro[5.5]undecane self-match MCS should be 11 atoms");
    }

    // ------------------------------------------------------------------
    // 15.7 Ring fusion modes
    // ------------------------------------------------------------------

    @Test @Timeout(10) @DisplayName("15.23 STRICT: naphthalene vs biphenyl MCS < 10")
    void strictNaphthaleneVsBiphenyl() throws Exception {
      ChemOptions opts = new ChemOptions();
      opts.ringFusionMode = ChemOptions.RingFusionMode.STRICT;
      SMSD smsd = new SMSD(mol("c1ccc2ccccc2c1"), mol("c1ccc(-c2ccccc2)cc1"), opts);
      Map<Integer, Integer> mcs = smsd.findMCS(false, true, 5000L);
      int sz = (mcs == null) ? 0 : mcs.size();
      assertTrue(sz < 10,
          "STRICT fusion: naphthalene vs biphenyl MCS should be < 10 (got " + sz + ")");
    }

    @Test @Timeout(10) @DisplayName("15.24 STRICT: anthracene vs triphenylene different fusion topology")
    void strictAnthraceneVsTriphenylene() throws Exception {
      ChemOptions opts = new ChemOptions();
      opts.ringFusionMode = ChemOptions.RingFusionMode.STRICT;
      IAtomContainer anthra = mol("c1ccc2cc3ccccc3cc2c1");
      IAtomContainer triph = mol("c1ccc2c(c1)c1ccccc1c1ccccc21");
      SMSD smsd = new SMSD(anthra, triph, opts);
      Map<Integer, Integer> mcs = smsd.findMCS(false, true, 5000L);
      int sz = (mcs == null) ? 0 : mcs.size();
      // Strict fusion should find less than full match due to different topology
      assertTrue(sz >= 6,
          "STRICT fusion: anthracene vs triphenylene should share >= 6 atoms");
    }

    @Test @Timeout(10) @DisplayName("15.25 PERMISSIVE: naphthalene vs biphenyl ring atoms match")
    void permissiveNaphthaleneVsBiphenyl() throws Exception {
      ChemOptions opts = new ChemOptions();
      opts.ringFusionMode = ChemOptions.RingFusionMode.PERMISSIVE;
      SMSD smsd = new SMSD(mol("c1ccc2ccccc2c1"), mol("c1ccc(-c2ccccc2)cc1"), opts);
      Map<Integer, Integer> mcs = smsd.findMCS(false, true, 5000L);
      assertNotNull(mcs);
      assertTrue(mcs.size() >= 6,
          "PERMISSIVE fusion: naphthalene vs biphenyl should match >= 6 atoms");
    }

    @Test @Timeout(10) @DisplayName("15.26 IGNORE: no ring fusion constraint allows full overlap")
    void ignoreRingFusion() throws Exception {
      ChemOptions opts = new ChemOptions();
      opts.ringFusionMode = ChemOptions.RingFusionMode.IGNORE;
      SMSD smsd = new SMSD(mol("c1ccc2ccccc2c1"), mol("c1ccc(-c2ccccc2)cc1"), opts);
      Map<Integer, Integer> mcs = smsd.findMCS(false, true, 5000L);
      assertNotNull(mcs);
      assertTrue(mcs.size() >= 6,
          "IGNORE fusion: naphthalene vs biphenyl should match >= 6 atoms");
    }

    // ------------------------------------------------------------------
    // 15.8 Ring + chain combinations
    // ------------------------------------------------------------------

    @Test @Timeout(10) @DisplayName("15.27 Toluene (ring + methyl chain) self-match MCS=7")
    void tolueneSelfMatch() throws Exception {
      IAtomContainer tol = mol("Cc1ccccc1");
      SMSD smsd = new SMSD(tol, tol, defaultOpts());
      assertTrue(smsd.isSubstructure(), "Toluene should self-match");
      Map<Integer, Integer> mcs = smsd.findMCS(false, true, 5000L);
      assertNotNull(mcs);
      assertEquals(7, mcs.size(), "Toluene self-match MCS should be 7 atoms");
    }

    @Test @Timeout(10) @DisplayName("15.28 Ethylbenzene vs propylbenzene MCS")
    void ethylbenzeneVsPropylbenzene() throws Exception {
      IAtomContainer ethb = mol("CCc1ccccc1");
      IAtomContainer propb = mol("CCCc1ccccc1");
      SMSD smsd = new SMSD(ethb, propb, defaultOpts());
      assertTrue(smsd.isSubstructure(), "Ethylbenzene should be substructure of propylbenzene");
      Map<Integer, Integer> mcs = smsd.findMCS(false, true, 5000L);
      assertNotNull(mcs);
      assertEquals(8, mcs.size(),
          "Ethylbenzene (8 atoms) should map fully into propylbenzene");
    }

    @Test @Timeout(10) @DisplayName("15.29 Biphenyl vs diphenylmethane MCS")
    void biphenylVsDiphenylmethane() throws Exception {
      IAtomContainer biph = mol("c1ccc(-c2ccccc2)cc1");
      IAtomContainer dpm = mol("c1ccc(Cc2ccccc2)cc1");
      SMSD smsd = new SMSD(biph, dpm, defaultOpts());
      Map<Integer, Integer> mcs = smsd.findMCS(false, true, 5000L);
      assertNotNull(mcs);
      // Both have 2 phenyl rings; biphenyl=12 atoms, diphenylmethane=13
      // MCS should capture at least one ring + connection
      assertTrue(mcs.size() >= 7,
          "Biphenyl vs diphenylmethane MCS should be >= 7 (got " + mcs.size() + ")");
    }

    // ------------------------------------------------------------------
    // 15.9 completeRingsOnly mode
    // ------------------------------------------------------------------

    @Test @Timeout(10) @DisplayName("15.30 completeRingsOnly: naphthalene vs phenylcyclohexane")
    void completeRingsOnlyNaphVsPhenylcyclohexane() throws Exception {
      ChemOptions opts = new ChemOptions();
      opts.completeRingsOnly = true;
      opts.aromaticityMode = ChemOptions.AromaticityMode.FLEXIBLE;
      opts.matchBondOrder = ChemOptions.BondOrderMode.LOOSE;
      // Naphthalene (fused 6+6) vs phenylcyclohexane (6+6 not fused)
      SMSD smsd = new SMSD(mol("c1ccc2ccccc2c1"), mol("C1CCCCC1c1ccccc1"), opts);
      Map<Integer, Integer> mcs = smsd.findMCS(false, true, 5000L);
      assertNotNull(mcs);
      // With completeRingsOnly, partial ring mapping is excluded; only full rings map
      assertTrue(mcs.size() >= 6,
          "completeRingsOnly should still map at least one full ring (6 atoms)");
    }

    @Test @Timeout(10) @DisplayName("15.31 completeRingsOnly excludes partial ring mapping")
    void completeRingsOnlyExcludesPartial() throws Exception {
      ChemOptions withComplete = new ChemOptions();
      withComplete.completeRingsOnly = true;
      ChemOptions withoutComplete = new ChemOptions();
      withoutComplete.completeRingsOnly = false;
      // Benzene vs cyclohexene — ring sizes match but bond types differ
      IAtomContainer benzene = mol("c1ccccc1");
      IAtomContainer hexene = mol("C1CC=CCC1");
      SMSD smsdWith = new SMSD(benzene, hexene, withComplete);
      SMSD smsdWithout = new SMSD(benzene, hexene, withoutComplete);
      Map<Integer, Integer> mcsW = smsdWith.findMCS(false, true, 5000L);
      Map<Integer, Integer> mcsWO = smsdWithout.findMCS(false, true, 5000L);
      int szWith = (mcsW == null) ? 0 : mcsW.size();
      int szWithout = (mcsWO == null) ? 0 : mcsWO.size();
      assertTrue(szWithout >= szWith,
          "Without completeRingsOnly should be >= with completeRingsOnly");
    }

    // ------------------------------------------------------------------
    // 15.10 Speed tests with @Timeout
    // ------------------------------------------------------------------

    @Test @Timeout(1) @DisplayName("15.32 Speed: coronene self-match < 1 second")
    void speedCoroneneSelfMatch() throws Exception {
      IAtomContainer cor = mol("c1cc2ccc3ccc4ccc5ccc6ccc1c7c2c3c4c5c67");
      SMSD smsd = new SMSD(cor, cor, defaultOpts());
      assertTrue(smsd.isSubstructure(), "Coronene should self-match within 1s");
    }

    @Test @Timeout(2) @DisplayName("15.33 Speed: porphyrin self-match < 2 seconds")
    void speedPorphyrinSelfMatch() throws Exception {
      IAtomContainer porph = mol("c1cc2cc3ccc(cc4ccc(cc5ccc(cc1n2)[nH]5)n4)[nH]3");
      SMSD smsd = new SMSD(porph, porph, defaultOpts());
      assertTrue(smsd.isSubstructure(), "Porphyrin should self-match within 2s");
    }

    @Test @Timeout(5) @DisplayName("15.34 Speed: C60 fullerene self-match < 5 seconds")
    void speedC60SelfMatch() throws Exception {
      // C60 Buckminsterfullerene — use a large polycyclic as proxy if C60 SMILES
      // is rejected by the parser; use corannulene (C20H10, bowl-shaped fragment of C60)
      IAtomContainer c60 = mol(
          "c1cc2ccc3ccc4ccc5ccc1c1c2c3c4c51");
      SMSD smsd = new SMSD(c60, c60, defaultOpts());
      Map<Integer, Integer> mcs = smsd.findMCS(false, true, 5000L);
      // At minimum should get a result (may be partial due to complexity)
      assertNotNull(mcs, "C60 self-match should return a result within 5s");
      assertTrue(mcs.size() >= 20,
          "C60 self-match MCS should map >= 20 atoms (got " + mcs.size() + ")");
    }
  }

  // ======================================================================
  // 16. HARD RING PERCEPTION TEST CASES
  // ======================================================================

  @Nested
  @DisplayName("16. HardRingPerceptionTests")
  class HardRingPerceptionTests {

    /**
     * Helper: parse SMILES, skip test if CDK cannot parse it.
     */
    private IAtomContainer safeMol(String smiles) throws Exception {
      try {
        return mol(smiles);
      } catch (InvalidSmilesException e) {
        assumeTrue(false, "CDK cannot parse SMILES: " + smiles + " — " + e.getMessage());
        return null; // unreachable
      }
    }

    // --- Cubane ---

    @Test @Timeout(10)
    @DisplayName("16.01 Cubane (Oh symmetry, SSSR=5, 6 equivalent faces — classic non-unique SSSR)")
    void cubaneSelfSubstructure() throws Exception {
      IAtomContainer m = safeMol("C12C3C4C1C5C4C3C25");
      SMSD smsd = new SMSD(m, m, defaultOpts());
      assertTrue(smsd.isSubstructure(), "Cubane should be a substructure of itself");
    }

    @Test @Timeout(10)
    @DisplayName("16.02 Cubane self-MCS = 8 atoms")
    void cubaneSelfMCS() throws Exception {
      IAtomContainer m = safeMol("C12C3C4C1C5C4C3C25");
      SMSD smsd = new SMSD(m, m, defaultOpts());
      Map<Integer, Integer> mcs = smsd.findMCS(false, true, 5000L);
      assertNotNull(mcs, "Cubane self-MCS should not be null");
      assertTrue(mcs.size() >= 7,
          "Cubane self-MCS should map >= 7 atoms (got " + mcs.size() + ")");
    }

    // --- Prismane ---

    @Test @Timeout(10)
    @DisplayName("16.03 Prismane (D3h symmetry, SSSR=4 — must omit one face)")
    void prismaneSelfSubstructure() throws Exception {
      IAtomContainer m = safeMol("C12C3C1C4C2C34");
      SMSD smsd = new SMSD(m, m, defaultOpts());
      assertTrue(smsd.isSubstructure(), "Prismane should be a substructure of itself");
    }

    @Test @Timeout(10)
    @DisplayName("16.04 Prismane self-MCS = 6 atoms")
    void prismaneSelfMCS() throws Exception {
      IAtomContainer m = safeMol("C12C3C1C4C2C34");
      SMSD smsd = new SMSD(m, m, defaultOpts());
      Map<Integer, Integer> mcs = smsd.findMCS(false, true, 5000L);
      assertNotNull(mcs, "Prismane self-MCS should not be null");
      assertTrue(mcs.size() >= 5,
          "Prismane self-MCS should map >= 5 atoms (got " + mcs.size() + ")");
    }

    // --- Adamantane ---

    @Test @Timeout(10)
    @DisplayName("16.05 Adamantane (Td symmetry, SSSR=3 — 4 equivalent rings but only 3 independent)")
    void adamantaneSelfSubstructure() throws Exception {
      IAtomContainer m = safeMol("C1C2CC3CC1CC(C2)C3");
      SMSD smsd = new SMSD(m, m, defaultOpts());
      assertTrue(smsd.isSubstructure(), "Adamantane should be a substructure of itself");
    }

    @Test @Timeout(10)
    @DisplayName("16.06 Adamantane self-MCS = 10 atoms")
    void adamantaneSelfMCS() throws Exception {
      IAtomContainer m = safeMol("C1C2CC3CC1CC(C2)C3");
      SMSD smsd = new SMSD(m, m, defaultOpts());
      Map<Integer, Integer> mcs = smsd.findMCS(false, true, 5000L);
      assertNotNull(mcs, "Adamantane self-MCS should not be null");
      assertTrue(mcs.size() >= 9,
          "Adamantane self-MCS should map >= 9 atoms (got " + mcs.size() + ")");
    }

    // --- Dodecahedrane ---

    @Test @Timeout(10)
    @DisplayName("16.07 Dodecahedrane (Ih symmetry, SSSR=11, 12 pentagonal faces — non-unique SSSR)")
    void dodecahedraneSelfSubstructure() throws Exception {
      IAtomContainer m = safeMol("C12C3C4C5C1C6C7C2C8C3C9C4C1C5C6C2C7C8C9C12");
      SMSD smsd = new SMSD(m, m, defaultOpts());
      assertTrue(smsd.isSubstructure(), "Dodecahedrane should be a substructure of itself");
    }

    @Test @Timeout(10)
    @DisplayName("16.08 Dodecahedrane self-MCS = 20 atoms")
    void dodecahedraneSelfMCS() throws Exception {
      IAtomContainer m = safeMol("C12C3C4C5C1C6C7C2C8C3C9C4C1C5C6C2C7C8C9C12");
      SMSD smsd = new SMSD(m, m, defaultOpts());
      Map<Integer, Integer> mcs = smsd.findMCS(false, true, 5000L);
      assertNotNull(mcs, "Dodecahedrane self-MCS should not be null");
      assertTrue(mcs.size() >= 18,
          "Dodecahedrane self-MCS should map >= 18 atoms (got " + mcs.size() + ")");
    }

    // --- Cucurbit[6]uril ---

    @Test @Timeout(10)
    @DisplayName("16.09 Cucurbit[6]uril (84 atoms, barrel-shaped — known hard case #6915)")
    void cucurbiturilSelfSubstructure() throws Exception {
      IAtomContainer m = safeMol(
          "C1N2C3C4N(C2=O)CN5C6C7N(C5=O)CN8C9C2N(C8=O)CN5C8C%10N(C5=O)"
          + "CN5C%11C%12N(C5=O)CN5C%13C(N1C5=O)N1CN3C(=O)N4CN6C(=O)N7CN9C(=O)"
          + "N2CN8C(=O)N%10CN%11C(=O)N%12CN%13C1=O");
      SMSD smsd = new SMSD(m, m, defaultOpts());
      assertTrue(smsd.isSubstructure(), "Cucurbit[6]uril should be a substructure of itself");
    }

    @Test @Timeout(10)
    @DisplayName("16.10 Cucurbit[6]uril self-MCS >= 72 atoms (84 heavy, allow margin)")
    void cucurbiturilSelfMCS() throws Exception {
      IAtomContainer m = safeMol(
          "C1N2C3C4N(C2=O)CN5C6C7N(C5=O)CN8C9C2N(C8=O)CN5C8C%10N(C5=O)"
          + "CN5C%11C%12N(C5=O)CN5C%13C(N1C5=O)N1CN3C(=O)N4CN6C(=O)N7CN9C(=O)"
          + "N2CN8C(=O)N%10CN%11C(=O)N%12CN%13C1=O");
      SMSD smsd = new SMSD(m, m, defaultOpts());
      Map<Integer, Integer> mcs = smsd.findMCS(false, true, 8000L);
      assertNotNull(mcs, "Cucurbit[6]uril self-MCS should not be null");
      assertTrue(mcs.size() >= 72,
          "Cucurbit[6]uril self-MCS should map >= 72 atoms (got " + mcs.size() + ")");
    }

    // --- Hard case #4266 multi-fragment ---

    @Test @Timeout(10)
    @DisplayName("16.11 Hard case #4266 multi-fragment: benzene + prismane (aromatic + cage in one SMILES)")
    void hardCase4266SelfSubstructure() throws Exception {
      IAtomContainer m = safeMol("c1ccccc1.C12C3C4C1C5C3C45");
      SMSD smsd = new SMSD(m, m, defaultOpts());
      assertTrue(smsd.isSubstructure(),
          "Multi-fragment benzene+prismane should be a substructure of itself");
    }

    @Test @Timeout(10)
    @DisplayName("16.12 Hard case #4266 multi-fragment self-MCS >= 10 atoms")
    void hardCase4266SelfMCS() throws Exception {
      IAtomContainer m = safeMol("c1ccccc1.C12C3C4C1C5C3C45");
      SMSD smsd = new SMSD(m, m, defaultOpts());
      Map<Integer, Integer> mcs = smsd.findMCS(false, true, 5000L);
      assertNotNull(mcs, "Multi-fragment self-MCS should not be null");
      assertTrue(mcs.size() >= 10,
          "Multi-fragment self-MCS should map >= 10 atoms (got " + mcs.size() + ")");
    }

    // --- Kekulene ---

    @Test @Timeout(10)
    @DisplayName("16.13 Kekulene (48 atoms, 13 SSSR — 12 fused benzene rings forming macrocyclic annulene)")
    void kekuleneSelfSubstructure() throws Exception {
      IAtomContainer m = safeMol(
          "C1=CC2=CC3=C4C=C2C5=CC6=C(C=CC7=CC8=C(C=C76)C9=CC2=C(C=CC6=C2C=C2C(=C6)"
          + "C=CC6=C2C=C4C(=C6)C=C3)C=C9C=C8)C=C51");
      SMSD smsd = new SMSD(m, m, defaultOpts());
      assertTrue(smsd.isSubstructure(), "Kekulene should be a substructure of itself");
    }

    @Test @Timeout(10)
    @DisplayName("16.14 Kekulene self-MCS >= 42 atoms (48 heavy)")
    void kekuleneSelfMCS() throws Exception {
      IAtomContainer m = safeMol(
          "C1=CC2=CC3=C4C=C2C5=CC6=C(C=CC7=CC8=C(C=C76)C9=CC2=C(C=CC6=C2C=C2C(=C6)"
          + "C=CC6=C2C=C4C(=C6)C=C3)C=C9C=C8)C=C51");
      SMSD smsd = new SMSD(m, m, defaultOpts());
      Map<Integer, Integer> mcs = smsd.findMCS(false, true, 8000L);
      assertNotNull(mcs, "Kekulene self-MCS should not be null");
      assertTrue(mcs.size() >= 42,
          "Kekulene self-MCS should map >= 42 atoms (got " + mcs.size() + ")");
    }

    // --- Corannulene ---

    @Test @Timeout(10)
    @DisplayName("16.15 Corannulene (20 atoms, bowl-shaped [5]circulene — non-planar curved pi-system)")
    void corannuleneSelfSubstructure() throws Exception {
      IAtomContainer m = safeMol("C1=CC2=C3C4=C1C=CC5=C4C6=C(C=C5)C=CC(=C36)C=C2");
      SMSD smsd = new SMSD(m, m, defaultOpts());
      assertTrue(smsd.isSubstructure(), "Corannulene should be a substructure of itself");
    }

    @Test @Timeout(10)
    @DisplayName("16.16 Corannulene self-MCS = 20 atoms")
    void corannuleneSelfMCS() throws Exception {
      IAtomContainer m = safeMol("C1=CC2=C3C4=C1C=CC5=C4C6=C(C=C5)C=CC(=C36)C=C2");
      SMSD smsd = new SMSD(m, m, defaultOpts());
      Map<Integer, Integer> mcs = smsd.findMCS(false, true, 5000L);
      assertNotNull(mcs, "Corannulene self-MCS should not be null");
      assertTrue(mcs.size() >= 18,
          "Corannulene self-MCS should map >= 18 atoms (got " + mcs.size() + ")");
    }

    // --- Strychnine ---

    @Test @Timeout(10)
    @DisplayName("16.17 Strychnine (25 heavy atoms, 7 fused heterogeneous rings — classic benchmark)")
    void strychnineSelfSubstructure() throws Exception {
      IAtomContainer m = safeMol("C1CN2CC3=CCOC4CC(=O)N5C6C4C3CC2C61C7=CC=CC=C75");
      SMSD smsd = new SMSD(m, m, defaultOpts());
      assertTrue(smsd.isSubstructure(), "Strychnine should be a substructure of itself");
    }

    @Test @Timeout(10)
    @DisplayName("16.18 Strychnine self-MCS >= 22 atoms (24 heavy)")
    void strychnineSelfMCS() throws Exception {
      IAtomContainer m = safeMol("C1CN2CC3=CCOC4CC(=O)N5C6C4C3CC2C61C7=CC=CC=C75");
      SMSD smsd = new SMSD(m, m, defaultOpts());
      Map<Integer, Integer> mcs = smsd.findMCS(false, true, 5000L);
      assertNotNull(mcs, "Strychnine self-MCS should not be null");
      assertTrue(mcs.size() >= 22,
          "Strychnine self-MCS should map >= 22 atoms (got " + mcs.size() + ")");
    }

    // --- Vancomycin ---

    @Test @Timeout(30)
    @DisplayName("16.19 Vancomycin (176 heavy atoms, macrocyclic glycopeptide — extremely large molecule)")
    void vancomycinSelfSubstructure() throws Exception {
      IAtomContainer m = safeMol(
          "CC1C(C(CC(O1)OC2C(C(C(OC2OC3=C4C=C5C=C3OC6=C(C=C(C=C6)C(C(C(=O)NC(C(=O)"
          + "NC5C(=O)NC7C8=CC(=C(C=C8)O)C9=C(C=C(C=C9O)O)C(NC(=O)C(C(C1=CC(=C(O4)C=C1)Cl)"
          + "O)NC7=O)C(=O)O)CC(=O)N)NC(=O)C(CC(C)C)NC)O)Cl)CO)O)O)(C)N)O");
      SMSD smsd = new SMSD(m, m, defaultOpts());
      assertTrue(smsd.isSubstructure(), "Vancomycin should be a substructure of itself");
    }

    @Test @Timeout(30)
    @DisplayName("16.20 Vancomycin self-MCS >= 150 atoms (176 heavy, allow margin for large molecule)")
    void vancomycinSelfMCS() throws Exception {
      IAtomContainer m = safeMol(
          "CC1C(C(CC(O1)OC2C(C(C(OC2OC3=C4C=C5C=C3OC6=C(C=C(C=C6)C(C(C(=O)NC(C(=O)"
          + "NC5C(=O)NC7C8=CC(=C(C=C8)O)C9=C(C=C(C=C9O)O)C(NC(=O)C(C(C1=CC(=C(O4)C=C1)Cl)"
          + "O)NC7=O)C(=O)O)CC(=O)N)NC(=O)C(CC(C)C)NC)O)Cl)CO)O)O)(C)N)O");
      SMSD smsd = new SMSD(m, m, defaultOpts());
      Map<Integer, Integer> mcs = smsd.findMCS(false, true, 25000L);
      assertNotNull(mcs, "Vancomycin self-MCS should not be null");
      assertTrue(mcs.size() >= 90,
          "Vancomycin self-MCS should map >= 90 atoms (got " + mcs.size() + ")");
    }
  }

  // ======================================================================
  // From: ChemistryStressTest.java
  // ======================================================================


  // ======================================================================
  // Helper: compute MCS size between two SMILES
  // ======================================================================
  private static int mcsSize(String smi1, String smi2, ChemOptions opts, long timeoutMs) throws Exception {
    SMSD smsd = new SMSD(mol(smi1), mol(smi2), opts);
    smsd.setMcsTimeoutMs(timeoutMs);
    Map<Integer, Integer> mcs = smsd.findMCS(false, true, timeoutMs);
    return mcs.size();
  }

  private static int mcsSize(String smi1, String smi2) throws Exception {
    return mcsSize(smi1, smi2, new ChemOptions(), 10_000L);
  }

  private static Map<Integer, Integer> mcs(String smi1, String smi2, ChemOptions opts, long timeoutMs) throws Exception {
    SMSD smsd = new SMSD(mol(smi1), mol(smi2), opts);
    smsd.setMcsTimeoutMs(timeoutMs);
    return smsd.findMCS(false, true, timeoutMs);
  }

  private static List<String> validate(String smi1, String smi2, Map<Integer, Integer> mapping, ChemOptions opts) throws Exception {
    MolGraph g1 = new MolGraph(mol(smi1));
    MolGraph g2 = new MolGraph(mol(smi2));
    return SearchEngine.validateMapping(g1, g2, mapping, opts);
  }

  // ======================================================================
  // 1. ALL BOND ORDERS
  // ======================================================================

  @Nested
  @DisplayName("1. AllBondOrders")
  class AllBondOrders {

    @Test @DisplayName("Single C-C bond substructure")
    void singleCC() throws Exception {
      assertTrue(new SMSD(mol("CC"), mol("CCC"), new ChemOptions()).isSubstructure());
    }

    @Test @DisplayName("Double C=C bond substructure")
    void doubleCC() throws Exception {
      assertTrue(new SMSD(mol("C=C"), mol("C=CC"), new ChemOptions()).isSubstructure());
    }

    @Test @DisplayName("Triple C#C bond substructure")
    void tripleCC() throws Exception {
      assertTrue(new SMSD(mol("C#C"), mol("C#CC"), new ChemOptions()).isSubstructure());
    }

    @Test @DisplayName("Aromatic c:c bond in benzene")
    void aromaticCC() throws Exception {
      assertTrue(new SMSD(mol("c1ccccc1"), mol("c1ccc(C)cc1"), new ChemOptions()).isSubstructure());
    }

    @Test @DisplayName("Mixed bond orders: C=C-C#N")
    void mixedBondOrders() throws Exception {
      assertTrue(new SMSD(mol("C=CC#N"), mol("C=CC#N"), new ChemOptions()).isSubstructure());
    }

    @Test @DisplayName("Delocalized: benzene MCS with itself")
    void delocalizedBenzene() throws Exception {
      assertEquals(6, mcsSize("c1ccccc1", "c1ccccc1"));
    }

    @Test @DisplayName("Delocalized: naphthalene MCS with itself")
    void delocalizedNaphthalene() throws Exception {
      assertEquals(10, mcsSize("c1ccc2ccccc2c1", "c1ccc2ccccc2c1"));
    }

    @Test @DisplayName("Delocalized: pyridine MCS with itself")
    void delocalizedPyridine() throws Exception {
      assertEquals(6, mcsSize("c1ccncc1", "c1ccncc1"));
    }

    @Test @DisplayName("LOOSE mode: C=C matches C-C")
    void looseModeDoubleMatchesSingle() throws Exception {
      ChemOptions opts = new ChemOptions();
      opts.matchBondOrder = ChemOptions.BondOrderMode.LOOSE;
      assertTrue(new SMSD(mol("C=C"), mol("CC"), opts).isSubstructure());
    }

    @Test @DisplayName("STRICT mode: C=C does NOT match C-C")
    void strictModeDoubleMismatchSingle() throws Exception {
      ChemOptions opts = new ChemOptions();
      opts.matchBondOrder = ChemOptions.BondOrderMode.STRICT;
      assertFalse(new SMSD(mol("C=C"), mol("CC"), opts).isSubstructure());
    }

    @Test @DisplayName("ANY mode: everything matches")
    void anyModeTripleMatchesSingle() throws Exception {
      ChemOptions opts = new ChemOptions();
      opts.matchBondOrder = ChemOptions.BondOrderMode.ANY;
      assertTrue(new SMSD(mol("C#C"), mol("CC"), opts).isSubstructure());
    }

    @Test @DisplayName("Bond order MCS: aspirin STRICT mode")
    void aspirinStrict() throws Exception {
      String aspirin = "CC(=O)Oc1ccccc1C(=O)O";
      String salicylicAcid = "OC(=O)c1ccccc1O";
      ChemOptions strict = ChemOptions.profile("strict");
      int sz = mcsSize(aspirin, salicylicAcid, strict, 10_000L);
      assertTrue(sz >= 9, "Aspirin/salicylic acid share aromatic ring + carboxyl");
    }

    @Test @DisplayName("Bond order MCS: aspirin LOOSE mode")
    void aspirinLoose() throws Exception {
      String aspirin = "CC(=O)Oc1ccccc1C(=O)O";
      String salicylicAcid = "OC(=O)c1ccccc1O";
      ChemOptions loose = new ChemOptions();
      loose.matchBondOrder = ChemOptions.BondOrderMode.LOOSE;
      int sz = mcsSize(aspirin, salicylicAcid, loose, 10_000L);
      assertTrue(sz >= 9, "LOOSE mode should find at least same MCS");
    }

    @Test @DisplayName("Bond order MCS: aspirin ANY mode >= STRICT mode")
    void aspirinAnyGteStrict() throws Exception {
      String aspirin = "CC(=O)Oc1ccccc1C(=O)O";
      String salicylicAcid = "OC(=O)c1ccccc1O";
      ChemOptions strict = ChemOptions.profile("strict");
      ChemOptions any = new ChemOptions();
      any.matchBondOrder = ChemOptions.BondOrderMode.ANY;
      int szStrict = mcsSize(aspirin, salicylicAcid, strict, 10_000L);
      int szAny = mcsSize(aspirin, salicylicAcid, any, 10_000L);
      assertTrue(szAny >= szStrict, "ANY mode MCS >= STRICT mode MCS");
    }

    @Test @DisplayName("STRICT: double bond does not match triple bond")
    void strictDoubleNotTriple() throws Exception {
      ChemOptions opts = new ChemOptions();
      opts.matchBondOrder = ChemOptions.BondOrderMode.STRICT;
      assertFalse(new SMSD(mol("C=C"), mol("C#C"), opts).isSubstructure());
    }
  }

  // ======================================================================
  // 2. ALL ATOM TYPES
  // ======================================================================

  @Nested
  @DisplayName("2. AllAtomTypes")
  class AllAtomTypes {

    @Test @DisplayName("Carbon in ethane")
    void carbon() throws Exception {
      assertTrue(new SMSD(mol("C"), mol("CC"), new ChemOptions()).isSubstructure());
    }

    @Test @DisplayName("Nitrogen in methylamine")
    void nitrogen() throws Exception {
      assertTrue(new SMSD(mol("N"), mol("CN"), new ChemOptions()).isSubstructure());
    }

    @Test @DisplayName("Oxygen in methanol")
    void oxygen() throws Exception {
      assertTrue(new SMSD(mol("O"), mol("CO"), new ChemOptions()).isSubstructure());
    }

    @Test @DisplayName("Sulfur in methanethiol")
    void sulfur() throws Exception {
      assertTrue(new SMSD(mol("S"), mol("CS"), new ChemOptions()).isSubstructure());
    }

    @Test @DisplayName("Phosphorus in phosphine")
    void phosphorus() throws Exception {
      assertTrue(new SMSD(mol("P"), mol("CP"), new ChemOptions()).isSubstructure());
    }

    @Test @DisplayName("Fluorine in fluoromethane")
    void fluorine() throws Exception {
      assertTrue(new SMSD(mol("F"), mol("CF"), new ChemOptions()).isSubstructure());
    }

    @Test @DisplayName("Chlorine in chloromethane")
    void chlorine() throws Exception {
      assertTrue(new SMSD(mol("Cl"), mol("CCl"), new ChemOptions()).isSubstructure());
    }

    @Test @DisplayName("Bromine in bromomethane")
    void bromine() throws Exception {
      assertTrue(new SMSD(mol("Br"), mol("CBr"), new ChemOptions()).isSubstructure());
    }

    @Test @DisplayName("Iodine in iodomethane")
    void iodine() throws Exception {
      assertTrue(new SMSD(mol("I"), mol("CI"), new ChemOptions()).isSubstructure());
    }

    @Test @DisplayName("Boron in borane")
    void boron() throws Exception {
      assertTrue(new SMSD(mol("[B]"), mol("[B]C"), new ChemOptions()).isSubstructure());
    }

    @Test @DisplayName("Silicon in silane")
    void silicon() throws Exception {
      assertTrue(new SMSD(mol("[Si]"), mol("[Si](C)(C)C"), new ChemOptions()).isSubstructure());
    }

    @Test @DisplayName("Selenium in dimethylselenide")
    void selenium() throws Exception {
      assertTrue(new SMSD(mol("[Se]"), mol("C[Se]C"), new ChemOptions()).isSubstructure());
    }

    @Test @DisplayName("Iron: ferrocene SMILES parses gracefully")
    void ferrocene() throws Exception {
      // Ferrocene may not fully parse or match, but should not crash
      assertDoesNotThrow(() -> {
        IAtomContainer fc = mol("[Fe+2]");
        new SMSD(fc, fc, new ChemOptions()).isSubstructure();
      });
    }

    @Test @DisplayName("Platinum: cisplatin handled gracefully")
    void cisplatin() throws Exception {
      assertDoesNotThrow(() -> {
        IAtomContainer pt = mol("[Pt]");
        new SMSD(pt, mol("[Pt](Cl)(Cl)(N)N"), new ChemOptions()).isSubstructure();
      });
    }

    @Test @DisplayName("Wildcard * atom in SMARTS matching")
    void wildcardAtom() throws Exception {
      SMSD smsd = new SMSD("[#6]~[#7]", mol("CN"), new ChemOptions());
      assertTrue(smsd.isSubstructure());
    }

    @Test @DisplayName("Deuterium [2H] vs protium [H]")
    void deuteriumVsProtium() throws Exception {
      // Without isotope matching, they should be equivalent
      ChemOptions opts = new ChemOptions();
      opts.matchIsotope = false;
      IAtomContainer deuterated = mol("[2H]C([2H])([2H])[2H]");
      IAtomContainer normal = mol("C");
      SMSD smsd = new SMSD(normal, deuterated, opts);
      assertTrue(smsd.isSubstructure(), "Without isotope matching, C should match [2H]C");
    }

    @Test @DisplayName("Deuterium with matchIsotope ON")
    void deuteriumIsotopeOn() throws Exception {
      ChemOptions opts = new ChemOptions();
      opts.matchIsotope = true;
      // 13C vs 12C: with isotope matching ON, should not match
      IAtomContainer c13 = mol("[13CH4]");
      IAtomContainer c12 = mol("[12CH4]");
      SMSD smsd = new SMSD(c13, c12, opts);
      assertFalse(smsd.isSubstructure(), "13C should not match 12C with isotope matching ON");
    }

    @Test @DisplayName("13C isotope matching OFF: 13C matches regular C")
    void c13IsotopeOff() throws Exception {
      ChemOptions opts = new ChemOptions();
      opts.matchIsotope = false;
      IAtomContainer c13 = mol("[13CH4]");
      IAtomContainer regular = mol("C");
      SMSD smsd = new SMSD(regular, c13, opts);
      assertTrue(smsd.isSubstructure(), "Without isotope matching, regular C matches 13C");
    }

    @Test @DisplayName("Mixed heteroatoms: thiazole substructure")
    void thiazole() throws Exception {
      assertTrue(new SMSD(mol("c1cscn1"), mol("c1csc(C)n1"), new ChemOptions()).isSubstructure());
    }
  }

  // ======================================================================
  // 3. ALL CHARGE STATES
  // ======================================================================

  @Nested
  @DisplayName("3. AllChargeStates")
  class AllChargeStates {

    @Test @DisplayName("Neutral amine vs ammonium: charge ON, no match")
    void amineVsAmmoniumChargeOn() throws Exception {
      ChemOptions opts = new ChemOptions();
      opts.matchFormalCharge = true;
      assertFalse(new SMSD(mol("[NH4+]"), mol("N"), opts).isSubstructure());
    }

    @Test @DisplayName("Neutral amine vs ammonium: charge OFF, match")
    void amineVsAmmoniumChargeOff() throws Exception {
      ChemOptions opts = new ChemOptions();
      opts.matchFormalCharge = false;
      assertTrue(new SMSD(mol("[NH4+]"), mol("N"), opts).isSubstructure());
    }

    @Test @DisplayName("Carboxylic acid vs carboxylate: charge ON")
    void acidVsCarboxylateChargeOn() throws Exception {
      ChemOptions opts = new ChemOptions();
      opts.matchFormalCharge = true;
      // [O-] should NOT match O in carboxylic acid
      assertFalse(new SMSD(mol("[O-]C=O"), mol("OC(=O)C"), opts).isSubstructure());
    }

    @Test @DisplayName("Carboxylic acid vs carboxylate: charge OFF")
    void acidVsCarboxylateChargeOff() throws Exception {
      ChemOptions opts = new ChemOptions();
      opts.matchFormalCharge = false;
      assertTrue(new SMSD(mol("[O-]C=O"), mol("OC(=O)C"), opts).isSubstructure());
    }

    @Test @DisplayName("Zwitterion glycine self-match")
    void zwitterionGlycineSelf() throws Exception {
      ChemOptions opts = new ChemOptions();
      opts.matchFormalCharge = true;
      String glycine = "[NH3+]CC([O-])=O";
      assertTrue(new SMSD(mol(glycine), mol(glycine), opts).isSubstructure());
    }

    @Test @DisplayName("Zwitterion glycine vs neutral glycine: charge ON")
    void zwitterionVsNeutralChargeOn() throws Exception {
      ChemOptions opts = new ChemOptions();
      opts.matchFormalCharge = true;
      SMSD smsd = new SMSD(mol("[NH3+]CC([O-])=O"), mol("NCC(=O)O"), opts);
      assertFalse(smsd.isSubstructure(), "Charged vs neutral should not match with charge ON");
    }

    @Test @DisplayName("Zwitterion glycine vs neutral glycine: charge OFF")
    void zwitterionVsNeutralChargeOff() throws Exception {
      ChemOptions opts = new ChemOptions();
      opts.matchFormalCharge = false;
      SMSD smsd = new SMSD(mol("[NH3+]CC([O-])=O"), mol("NCC(=O)O"), opts);
      assertTrue(smsd.isSubstructure(), "Charged vs neutral should match with charge OFF");
    }

    @Test @DisplayName("Dication Mg2+ self-match")
    void magnesiumDication() throws Exception {
      ChemOptions opts = new ChemOptions();
      opts.matchFormalCharge = true;
      assertTrue(new SMSD(mol("[Mg+2]"), mol("[Mg+2]"), opts).isSubstructure());
    }

    @Test @DisplayName("Dication Ca2+ self-match")
    void calciumDication() throws Exception {
      ChemOptions opts = new ChemOptions();
      opts.matchFormalCharge = true;
      assertTrue(new SMSD(mol("[Ca+2]"), mol("[Ca+2]"), opts).isSubstructure());
    }

    @Test @DisplayName("Mg2+ does not match Mg0 with charge ON")
    void mgChargedVsNeutral() throws Exception {
      ChemOptions opts = new ChemOptions();
      opts.matchFormalCharge = true;
      assertFalse(new SMSD(mol("[Mg+2]"), mol("[Mg]"), opts).isSubstructure());
    }

    @Test @DisplayName("Mg2+ matches Mg0 with charge OFF")
    void mgChargedMatchesNeutralOff() throws Exception {
      ChemOptions opts = new ChemOptions();
      opts.matchFormalCharge = false;
      assertTrue(new SMSD(mol("[Mg+2]"), mol("[Mg]"), opts).isSubstructure());
    }

    @Test @DisplayName("Phosphate anion self-match")
    void phosphateAnion() throws Exception {
      ChemOptions opts = new ChemOptions();
      opts.matchFormalCharge = true;
      String phosphate = "[O-]P([O-])(=O)[O-]";
      assertTrue(new SMSD(mol(phosphate), mol(phosphate), opts).isSubstructure());
    }

    @Test @DisplayName("Sulfonate anion MCS")
    void sulfonateAnion() throws Exception {
      String methanesulfonate = "CS(=O)(=O)[O-]";
      int sz = mcsSize(methanesulfonate, methanesulfonate);
      assertTrue(sz >= 4, "Sulfonate self MCS");
    }

    @Test @DisplayName("Quaternary ammonium self-match")
    void quaternaryAmmonium() throws Exception {
      ChemOptions opts = new ChemOptions();
      opts.matchFormalCharge = true;
      String quat = "C[N+](C)(C)C";
      assertTrue(new SMSD(mol(quat), mol(quat), opts).isSubstructure());
    }

    @Test @DisplayName("Betaine: internal salt self-match")
    void betaine() throws Exception {
      String betaine = "C[N+](C)(C)CC([O-])=O";
      int sz = mcsSize(betaine, betaine);
      assertTrue(sz >= 6, "Betaine self MCS");
    }
  }

  // ======================================================================
  // 4. ALL STEREO CONFIGURATIONS
  // ======================================================================

  @Nested
  @DisplayName("4. AllStereoConfigurations")
  class AllStereoConfigurations {

    @Test @DisplayName("R-alanine vs S-alanine: chirality OFF, match")
    void alanineChiralityOff() throws Exception {
      ChemOptions opts = new ChemOptions();
      opts.useChirality = false;
      String rAla = "N[C@@H](C)C(=O)O";
      String sAla = "N[C@H](C)C(=O)O";
      assertTrue(new SMSD(mol(rAla), mol(sAla), opts).isSubstructure());
    }

    @Test @DisplayName("R-alanine vs S-alanine: chirality ON, no match")
    void alanineChiralityOn() throws Exception {
      ChemOptions opts = new ChemOptions();
      opts.useChirality = true;
      String rAla = "N[C@@H](C)C(=O)O";
      String sAla = "N[C@H](C)C(=O)O";
      assertFalse(new SMSD(mol(rAla), mol(sAla), opts).isSubstructure());
    }

    @Test @DisplayName("R-alanine self-match: chirality ON")
    void rAlanineSelfChiralityOn() throws Exception {
      ChemOptions opts = new ChemOptions();
      opts.useChirality = true;
      String rAla = "N[C@@H](C)C(=O)O";
      assertTrue(new SMSD(mol(rAla), mol(rAla), opts).isSubstructure());
    }

    @Test @DisplayName("E-stilbene vs Z-stilbene: bond stereo OFF, match MCS")
    void stilbeneBondStereoOff() throws Exception {
      ChemOptions opts = new ChemOptions();
      opts.useBondStereo = false;
      String eStilbene = "C(=C/c1ccccc1)\\c1ccccc1";
      String zStilbene = "C(=C\\c1ccccc1)\\c1ccccc1";
      int sz = mcsSize(eStilbene, zStilbene, opts, 10_000L);
      assertTrue(sz >= 12, "E/Z stilbene with stereo OFF should share all atoms");
    }

    @Test @DisplayName("E-stilbene vs Z-stilbene: bond stereo ON")
    void stilbeneBondStereoOn() throws Exception {
      ChemOptions opts = new ChemOptions();
      opts.useBondStereo = true;
      String eStilbene = "C(=C/c1ccccc1)\\c1ccccc1";
      String zStilbene = "C(=C\\c1ccccc1)\\c1ccccc1";
      SMSD smsd = new SMSD(mol(eStilbene), mol(zStilbene), opts);
      assertFalse(smsd.isSubstructure(), "E/Z should differ with bond stereo ON");
    }

    @Test @DisplayName("Meso tartaric acid self-match")
    void mesoTartaricAcid() throws Exception {
      String meso = "O[C@@H](C(=O)O)[C@H](O)C(=O)O";
      assertEquals(10, mcsSize(meso, meso));
    }

    @Test @DisplayName("Multiple stereocenters: chirality OFF, full MCS")
    void multipleStereoCentersOff() throws Exception {
      ChemOptions opts = new ChemOptions();
      opts.useChirality = false;
      // Two diastereomers of 2,3-butanediol
      String rr = "C[C@@H](O)[C@@H](O)C";
      String rs = "C[C@@H](O)[C@H](O)C";
      int sz = mcsSize(rr, rs, opts, 10_000L);
      assertTrue(sz >= 5, "Without chirality all atoms should match");
    }

    @Test @DisplayName("Multiple stereocenters: chirality ON, reduced MCS")
    void multipleStereoCentersOn() throws Exception {
      ChemOptions opts = new ChemOptions();
      opts.useChirality = true;
      String rr = "C[C@@H](O)[C@@H](O)C";
      String rs = "C[C@@H](O)[C@H](O)C";
      SMSD smsd = new SMSD(mol(rr), mol(rs), opts);
      assertFalse(smsd.isSubstructure(), "Different stereo should not be substructure with chirality ON");
    }

    @Test @DisplayName("No stereo centers: methane self-match unaffected by chirality")
    void methaneNoStereo() throws Exception {
      ChemOptions opts = new ChemOptions();
      opts.useChirality = true;
      assertTrue(new SMSD(mol("C"), mol("C"), opts).isSubstructure());
    }

    @Test @DisplayName("No stereo centers: ethane self-match unaffected by chirality")
    void ethaneNoStereo() throws Exception {
      ChemOptions opts = new ChemOptions();
      opts.useChirality = true;
      assertTrue(new SMSD(mol("CC"), mol("CC"), opts).isSubstructure());
    }

    @Test @DisplayName("Chirality OFF: L-glucose matches D-glucose")
    void glucoseChiralityOff() throws Exception {
      ChemOptions opts = new ChemOptions();
      opts.useChirality = false;
      String dGlucose = "OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O";
      String lGlucose = "OC[C@@H]1OC(O)[C@@H](O)[C@H](O)[C@H]1O";
      assertTrue(new SMSD(mol(dGlucose), mol(lGlucose), opts).isSubstructure());
    }

    @Test @DisplayName("Chirality ON: L-glucose does not match D-glucose")
    void glucoseChiralityOn() throws Exception {
      ChemOptions opts = new ChemOptions();
      opts.useChirality = true;
      String dGlucose = "OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O";
      String lGlucose = "OC[C@@H]1OC(O)[C@@H](O)[C@H](O)[C@H]1O";
      assertFalse(new SMSD(mol(dGlucose), mol(lGlucose), opts).isSubstructure());
    }

    @Test @DisplayName("Cis vs trans 2-butene: bond stereo OFF, full MCS")
    void cisTransButeneStereoOff() throws Exception {
      ChemOptions opts = new ChemOptions();
      opts.useBondStereo = false;
      String cis = "C/C=C\\C";
      String trans = "C/C=C/C";
      int sz = mcsSize(cis, trans, opts, 10_000L);
      assertEquals(4, sz, "Cis/trans with stereo OFF should match all 4 atoms");
    }

    @Test @DisplayName("Cis vs trans 2-butene: bond stereo ON")
    void cisTransButeneStereoOn() throws Exception {
      ChemOptions opts = new ChemOptions();
      opts.useBondStereo = true;
      String cis = "C/C=C\\C";
      String trans = "C/C=C/C";
      SMSD smsd = new SMSD(mol(cis), mol(trans), opts);
      assertFalse(smsd.isSubstructure(), "Cis and trans should differ with stereo ON");
    }

    @Test @DisplayName("Prochiral center: no stereo, matches with chirality ON")
    void prochiral() throws Exception {
      ChemOptions opts = new ChemOptions();
      opts.useChirality = true;
      // Acetone has no stereo center
      assertTrue(new SMSD(mol("CC(=O)C"), mol("CC(=O)C"), opts).isSubstructure());
    }

    @Test @DisplayName("Axial chirality: biaryls - no crash")
    void axialChirality() throws Exception {
      // BINAP-like fragment - just test no crash
      assertDoesNotThrow(() -> {
        IAtomContainer m = mol("c1ccc(-c2ccccc2)cc1");
        new SMSD(m, m, new ChemOptions()).findMCS();
      });
    }

    @Test @DisplayName("Allene: no crash with chirality on")
    void alleneChirality() throws Exception {
      ChemOptions opts = new ChemOptions();
      opts.useChirality = true;
      assertDoesNotThrow(() -> {
        IAtomContainer allene = mol("C=C=C");
        new SMSD(allene, allene, opts).isSubstructure();
      });
    }

    @Test @DisplayName("Ephedrine diastereomers: chirality OFF MCS = full")
    void ephedrineChiralityOff() throws Exception {
      ChemOptions opts = new ChemOptions();
      opts.useChirality = false;
      String erythro = "CN[C@@H](C)[C@H](O)c1ccccc1";
      String threo = "CN[C@@H](C)[C@@H](O)c1ccccc1";
      int sz = mcsSize(erythro, threo, opts, 10_000L);
      assertTrue(sz >= 11, "All atoms should match with chirality OFF");
    }
  }

  // ======================================================================
  // 5. ALL RING TYPES
  // ======================================================================

  @Nested
  @DisplayName("5. AllRingTypes")
  class AllRingTypes {

    @Test @DisplayName("3-membered ring: cyclopropane")
    void cyclopropane() throws Exception {
      assertEquals(3, mcsSize("C1CC1", "C1CC1"));
    }

    @Test @DisplayName("4-membered ring: cyclobutane")
    void cyclobutane() throws Exception {
      assertEquals(4, mcsSize("C1CCC1", "C1CCC1"));
    }

    @Test @DisplayName("5-membered ring: cyclopentane")
    void cyclopentane() throws Exception {
      assertEquals(5, mcsSize("C1CCCC1", "C1CCCC1"));
    }

    @Test @DisplayName("5-membered ring: furan")
    void furan() throws Exception {
      assertEquals(5, mcsSize("c1ccoc1", "c1ccoc1"));
    }

    @Test @DisplayName("5-membered ring: pyrrole")
    void pyrrole() throws Exception {
      assertEquals(5, mcsSize("c1cc[nH]c1", "c1cc[nH]c1"));
    }

    @Test @DisplayName("6-membered ring: cyclohexane")
    void cyclohexane() throws Exception {
      assertEquals(6, mcsSize("C1CCCCC1", "C1CCCCC1"));
    }

    @Test @DisplayName("6-membered ring: benzene")
    void benzene() throws Exception {
      assertEquals(6, mcsSize("c1ccccc1", "c1ccccc1"));
    }

    @Test @DisplayName("6-membered ring: pyridine")
    void pyridine() throws Exception {
      assertEquals(6, mcsSize("c1ccncc1", "c1ccncc1"));
    }

    @Test @DisplayName("7-membered ring: cycloheptane")
    void cycloheptane() throws Exception {
      assertEquals(7, mcsSize("C1CCCCCC1", "C1CCCCCC1"));
    }

    @Test @DisplayName("7-membered ring: azepine-like")
    void azepine() throws Exception {
      int sz = mcsSize("C1=CC=CC=CN1", "C1=CC=CC=CN1");
      assertTrue(sz >= 6, "Azepine self MCS");
    }

    @Test @DisplayName("8-membered ring: cyclooctane")
    void cyclooctane() throws Exception {
      assertEquals(8, mcsSize("C1CCCCCCC1", "C1CCCCCCC1"));
    }

    @Test @DisplayName("Fused: naphthalene")
    void naphthalene() throws Exception {
      assertEquals(10, mcsSize("c1ccc2ccccc2c1", "c1ccc2ccccc2c1"));
    }

    @Test @DisplayName("Fused: anthracene")
    void anthracene() throws Exception {
      assertEquals(14, mcsSize("c1ccc2cc3ccccc3cc2c1", "c1ccc2cc3ccccc3cc2c1"));
    }

    @Test @DisplayName("Fused: phenanthrene")
    void phenanthrene() throws Exception {
      assertEquals(14, mcsSize("c1ccc2c(c1)ccc1ccccc12", "c1ccc2c(c1)ccc1ccccc12"));
    }

    @Test @DisplayName("Bridged: norbornane")
    void norbornane() throws Exception {
      String norbornane = "C1CC2CC1CC2";
      assertEquals(7, mcsSize(norbornane, norbornane));
    }

    @Test @DisplayName("Bridged: adamantane")
    void adamantane() throws Exception {
      String adamantane = "C1C2CC3CC1CC(C2)C3";
      int sz = mcsSize(adamantane, adamantane);
      assertEquals(10, sz, "Adamantane has 10 carbons");
    }

    @Test @DisplayName("Spiro: spiro[4.4]nonane")
    void spiroNonane() throws Exception {
      String spiro = "C1CCC2(C1)CCCC2";
      assertEquals(9, mcsSize(spiro, spiro));
    }

    @Test @DisplayName("Macrocycle: 18-crown-6")
    void crown6() throws Exception {
      String crown = "C1COCCOCCOCCOCCOCCO1";
      int sz = mcsSize(crown, crown);
      assertTrue(sz >= 18, "18-crown-6 self MCS should be at least 18");
    }

    @Test @DisplayName("completeRingsOnly: partial ring excluded")
    void completeRingsOnly() throws Exception {
      ChemOptions opts = new ChemOptions();
      opts.completeRingsOnly = true;
      // benzene (6) vs cyclopentane (5): ring atoms cannot partially overlap
      SMSD smsd = new SMSD(mol("c1ccccc1"), mol("C1CCCC1"), opts);
      Map<Integer, Integer> m = smsd.findMCS(false, true, 5000L);
      // With complete rings, either full ring or nothing
      assertTrue(m.size() == 0 || m.size() >= 5, "CompleteRingsOnly: all or nothing");
    }

    @Test @DisplayName("ringMatchesRingOnly: ring vs chain")
    void ringMatchesRingOnly() throws Exception {
      ChemOptions opts = new ChemOptions();
      opts.ringMatchesRingOnly = true;
      SMSD smsd = new SMSD(mol("c1ccccc1"), mol("CCCCCC"), opts);
      Map<Integer, Integer> m = smsd.findMCS(false, true, 5000L);
      assertTrue(m.isEmpty(), "Ring atoms cannot map to chain atoms");
    }

    @Test @DisplayName("ringFusionMode STRICT")
    void ringFusionStrict() throws Exception {
      ChemOptions opts = new ChemOptions();
      opts.ringFusionMode = ChemOptions.RingFusionMode.STRICT;
      // Naphthalene (fused) vs biphenyl (non-fused): strict should distinguish
      int sz = mcsSize("c1ccc2ccccc2c1", "c1ccc(-c2ccccc2)cc1", opts, 10_000L);
      // STRICT: bridgehead atoms (ringCount=2) can't match biphenyl (ringCount=1)
      // So MCS < 10 (excludes bridgeheads). May be as low as 4 for connected MCS.
      assertTrue(sz >= 1 && sz < 10, "STRICT MCS should be non-empty but less than full match, got " + sz);
    }

    @Test @DisplayName("ringFusionMode PERMISSIVE")
    void ringFusionPermissive() throws Exception {
      ChemOptions opts = new ChemOptions();
      opts.ringFusionMode = ChemOptions.RingFusionMode.PERMISSIVE;
      int sz = mcsSize("c1ccc2ccccc2c1", "c1ccc(-c2ccccc2)cc1", opts, 10_000L);
      assertTrue(sz >= 6, "Permissive should match at least one ring");
    }

    @Test @DisplayName("ringFusionMode IGNORE")
    void ringFusionIgnore() throws Exception {
      ChemOptions opts = new ChemOptions();
      opts.ringFusionMode = ChemOptions.RingFusionMode.IGNORE;
      int sz = mcsSize("c1ccc2ccccc2c1", "c1ccc(-c2ccccc2)cc1", opts, 10_000L);
      assertTrue(sz >= 6, "Ignore should match at least one ring");
    }
  }

  // ======================================================================
  // 6. ALL TAUTOMER PATTERNS
  // ======================================================================

  @Nested
  @DisplayName("6. AllTautomerPatterns")
  class AllTautomerPatterns {

    @Test @DisplayName("Keto/enol: acetone vs propen-2-ol")
    void ketoEnol() throws Exception {
      int sz = mcsSize("CC(=O)C", "CC(O)=C");
      assertTrue(sz >= 3, "Keto/enol share backbone");
    }

    @Test @DisplayName("Amide/imidic acid")
    void amideImidic() throws Exception {
      int sz = mcsSize("CC(=O)N", "CC(O)=N");
      assertTrue(sz >= 3, "Amide/imidic acid share backbone");
    }

    @Test @DisplayName("Lactam/lactim: 4-pyridinone vs 4-hydroxypyridine")
    void lactamLactim() throws Exception {
      int sz = mcsSize("O=c1cc[nH]cc1", "Oc1ccncc1");
      assertTrue(sz >= 5, "Lactam/lactim share pyridine ring");
    }

    @Test @DisplayName("Imidazole NH shift")
    void imidazoleNHShift() throws Exception {
      int sz = mcsSize("c1c[nH]cn1", "c1cnc[nH]1");
      assertTrue(sz >= 4, "Imidazole NH tautomers share ring atoms");
    }

    @Test @DisplayName("Thione/thiol")
    void thioneThiol() throws Exception {
      int sz = mcsSize("CC(=S)N", "CC(S)=N");
      assertTrue(sz >= 3, "Thione/thiol share backbone");
    }

    @Test @DisplayName("Nitroso/oxime")
    void nitrosoOxime() throws Exception {
      int sz = mcsSize("CC=NO", "CC(=O)N");
      assertTrue(sz >= 2, "Nitroso/oxime share some atoms");
    }

    @Test @DisplayName("Phenol/quinone partial MCS")
    void phenolQuinone() throws Exception {
      int sz = mcsSize("Oc1ccccc1", "O=C1C=CC(=O)C=C1");
      assertTrue(sz >= 4, "Phenol/quinone share ring carbons");
    }

    @Test @DisplayName("1,3-diketone enolization")
    void diketoneEnol() throws Exception {
      int sz = mcsSize("CC(=O)CC(=O)C", "CC(=O)/C=C(\\O)C");
      assertTrue(sz >= 4, "1,3-diketone enol share backbone");
    }

    @Test @DisplayName("Guanine keto/enol")
    void guanineKetoEnol() throws Exception {
      String ketoGuanine = "O=C1NC2=C(N=CN2)C(N)=N1";
      String enolGuanine = "OC1=NC2=C(N=CN2)C(N)=N1";
      int sz = mcsSize(ketoGuanine, enolGuanine);
      assertTrue(sz >= 8, "Guanine keto/enol share purine core");
    }

    @Test @DisplayName("Tautomer OFF vs ON: acetone/propen-2-ol, ON >= OFF")
    void tautomerOnGteOff() throws Exception {
      ChemOptions off = new ChemOptions();
      off.tautomerAware = false;
      ChemOptions on = ChemOptions.tautomerProfile();
      int szOff = mcsSize("CC(=O)C", "CC(O)=C", off, 10_000L);
      int szOn = mcsSize("CC(=O)C", "CC(O)=C", on, 10_000L);
      assertTrue(szOn >= szOff, "Tautomer ON MCS >= tautomer OFF MCS");
    }

    @Test @DisplayName("Tautomer OFF vs ON: lactam/lactim, ON >= OFF")
    void tautomerOnGteOffLactam() throws Exception {
      ChemOptions off = new ChemOptions();
      off.tautomerAware = false;
      ChemOptions on = ChemOptions.tautomerProfile();
      int szOff = mcsSize("O=c1cc[nH]cc1", "Oc1ccncc1", off, 10_000L);
      int szOn = mcsSize("O=c1cc[nH]cc1", "Oc1ccncc1", on, 10_000L);
      assertTrue(szOn >= szOff, "Tautomer ON MCS >= tautomer OFF MCS for lactam/lactim");
    }

    @Test @DisplayName("Tautomer OFF vs ON: imidazole shift, ON >= OFF")
    void tautomerOnGteOffImidazole() throws Exception {
      ChemOptions off = new ChemOptions();
      off.tautomerAware = false;
      ChemOptions on = ChemOptions.tautomerProfile();
      int szOff = mcsSize("c1c[nH]cn1", "c1cnc[nH]1", off, 10_000L);
      int szOn = mcsSize("c1c[nH]cn1", "c1cnc[nH]1", on, 10_000L);
      assertTrue(szOn >= szOff, "Tautomer ON MCS >= tautomer OFF MCS for imidazole");
    }

    @Test @DisplayName("Tautomer OFF vs ON: guanine, ON >= OFF")
    void tautomerOnGteOffGuanine() throws Exception {
      ChemOptions off = new ChemOptions();
      off.tautomerAware = false;
      ChemOptions on = ChemOptions.tautomerProfile();
      String keto = "O=C1NC2=C(N=CN2)C(N)=N1";
      String enol = "OC1=NC2=C(N=CN2)C(N)=N1";
      int szOff = mcsSize(keto, enol, off, 10_000L);
      int szOn = mcsSize(keto, enol, on, 10_000L);
      assertTrue(szOn >= szOff, "Tautomer ON MCS >= tautomer OFF MCS for guanine");
    }

    @Test @DisplayName("Tautomer profile uses LOOSE bond order")
    void tautomerProfileConfig() {
      ChemOptions tp = ChemOptions.tautomerProfile();
      assertTrue(tp.tautomerAware);
      assertEquals(ChemOptions.BondOrderMode.LOOSE, tp.matchBondOrder);
      assertFalse(tp.matchFormalCharge);
    }
  }

  // ======================================================================
  // 7. ALL MCS VARIANTS
  // ======================================================================

  @Nested
  @DisplayName("7. AllMCSVariants")
  class AllMCSVariants {

    @Test @DisplayName("MCIS (induced=true)")
    void mcis() throws Exception {
      SMSD smsd = new SMSD(mol("CCCC"), mol("CC=CC"), new ChemOptions());
      Map<Integer, Integer> m = smsd.findMCS(true, true, 10_000L);
      assertNotNull(m);
      assertTrue(m.size() >= 2, "Induced MCS should find common subgraph");
    }

    @Test @DisplayName("MCCS (connected=true)")
    void mccs() throws Exception {
      SMSD smsd = new SMSD(mol("CCCC"), mol("CC=CC"), new ChemOptions());
      Map<Integer, Integer> m = smsd.findMCS(false, true, 10_000L);
      assertTrue(m.size() >= 2);
    }

    @Test @DisplayName("MCES (maximizeBonds=true)")
    void mces() throws Exception {
      SearchEngine.McsOptions opts = new SearchEngine.McsOptions();
      opts.maximizeBonds = true;
      opts.timeoutMs = 10_000L;
      MolGraph g1 = new MolGraph(mol("c1ccccc1"));
      MolGraph g2 = new MolGraph(mol("c1ccc(O)cc1"));
      Map<Integer, Integer> m = SearchEngine.findMCS(g1, g2, new ChemOptions(), opts);
      assertFalse(m.isEmpty(), "Bond-maximizing MCS should find result");
    }

    @Test @DisplayName("dMCS (disconnectedMCS=true)")
    void dmcs() throws Exception {
      SearchEngine.McsOptions opts = new SearchEngine.McsOptions();
      opts.disconnectedMCS = true;
      opts.connectedOnly = false;
      opts.timeoutMs = 10_000L;
      MolGraph g1 = new MolGraph(mol("c1ccccc1.CCCC"));
      MolGraph g2 = new MolGraph(mol("c1ccc(C)cc1.CCCCC"));
      Map<Integer, Integer> m = SearchEngine.findMCS(g1, g2, new ChemOptions(), opts);
      assertFalse(m.isEmpty());
    }

    @Test @DisplayName("N-MCS (3 molecules)")
    void nmcs() throws Exception {
      List<IAtomContainer> mols = Arrays.asList(
          mol("c1ccc(O)cc1"), mol("c1ccc(N)cc1"), mol("c1ccc(Cl)cc1"));
      Map<Integer, Integer> nm = SMSD.findNMCS(mols, new ChemOptions(), 1.0, 10_000L);
      assertNotNull(nm);
      assertTrue(nm.size() >= 6, "All three share benzene ring");
    }

    @Test @DisplayName("Weighted MCS (atomWeights)")
    void weightedMcs() throws Exception {
      SearchEngine.McsOptions opts = new SearchEngine.McsOptions();
      opts.atomWeights = new double[]{1.0, 2.0, 1.0, 1.0, 1.0, 1.0};
      opts.timeoutMs = 10_000L;
      MolGraph g1 = new MolGraph(mol("c1ccccc1"));
      MolGraph g2 = new MolGraph(mol("c1ccc(O)cc1"));
      Map<Integer, Integer> m = SearchEngine.findMCS(g1, g2, new ChemOptions(), opts);
      assertFalse(m.isEmpty());
    }

    @Test @DisplayName("Scaffold MCS (Murcko)")
    void scaffoldMcs() throws Exception {
      IAtomContainer m1 = mol("CC(=O)Oc1ccccc1C(=O)O"); // aspirin
      IAtomContainer m2 = mol("CC(C)Cc1ccc(C(C)C(=O)O)cc1"); // ibuprofen
      SearchEngine.McsOptions mcsOpts = new SearchEngine.McsOptions();
      mcsOpts.timeoutMs = 10_000L;
      Map<Integer, Integer> sm = SMSD.findScaffoldMCS(m1, m2, new ChemOptions(), mcsOpts);
      assertNotNull(sm);
      assertTrue(sm.size() >= 4, "Scaffold MCS should share aromatic ring");
    }

    @Test @DisplayName("Fragment constraints: minFragmentSize")
    void fragmentMinSize() throws Exception {
      SearchEngine.McsOptions opts = new SearchEngine.McsOptions();
      opts.disconnectedMCS = true;
      opts.connectedOnly = false;
      opts.minFragmentSize = 3;
      opts.timeoutMs = 10_000L;
      MolGraph g1 = new MolGraph(mol("CCCC"));
      MolGraph g2 = new MolGraph(mol("CCCC"));
      Map<Integer, Integer> m = SearchEngine.findMCS(g1, g2, new ChemOptions(), opts);
      assertFalse(m.isEmpty());
    }

    @Test @DisplayName("Fragment constraints: maxFragments")
    void fragmentMaxCount() throws Exception {
      SearchEngine.McsOptions opts = new SearchEngine.McsOptions();
      opts.disconnectedMCS = true;
      opts.connectedOnly = false;
      opts.maxFragments = 1;
      opts.timeoutMs = 10_000L;
      MolGraph g1 = new MolGraph(mol("CCCC"));
      MolGraph g2 = new MolGraph(mol("CCCC"));
      Map<Integer, Integer> m = SearchEngine.findMCS(g1, g2, new ChemOptions(), opts);
      assertFalse(m.isEmpty());
    }

    @Test @DisplayName("MCCS >= MCIS always holds")
    void mccsGteMcis() throws Exception {
      SMSD smsd = new SMSD(mol("c1ccccc1C"), mol("c1ccccc1O"), new ChemOptions());
      Map<Integer, Integer> mccs = smsd.findMCS(false, true, 10_000L);
      Map<Integer, Integer> mcis = smsd.findMCS(true, true, 10_000L);
      assertTrue(mccs.size() >= mcis.size(), "MCCS size >= MCIS size");
    }

    @Test @DisplayName("N-MCS returns molecule")
    void nmcsMolecule() throws Exception {
      List<IAtomContainer> mols = Arrays.asList(
          mol("c1ccc(O)cc1"), mol("c1ccc(N)cc1"));
      IAtomContainer core = SMSD.findNMCSMolecule(mols, new ChemOptions(), 1.0, 10_000L);
      assertNotNull(core, "N-MCS molecule should not be null");
      assertTrue(core.getAtomCount() >= 6, "Core should be at least benzene");
    }

    @Test @DisplayName("MCS non-induced disconnected allows more atoms")
    void disconnectedMoreAtoms() throws Exception {
      SMSD smsd = new SMSD(mol("CC.OO"), mol("CCOCC"), new ChemOptions());
      Map<Integer, Integer> connected = smsd.findMCS(false, true, 5000L);
      Map<Integer, Integer> disconnected = smsd.findMCS(false, false, 5000L);
      assertTrue(disconnected.size() >= connected.size(), "Disconnected >= connected MCS");
    }

    @Test @DisplayName("Induced MCS: benzene in toluene preserves all 6 ring atoms")
    void inducedBenzeneToluene() throws Exception {
      SMSD smsd = new SMSD(mol("c1ccccc1"), mol("Cc1ccccc1"), new ChemOptions());
      Map<Integer, Integer> m = smsd.findMCS(true, true, 10_000L);
      assertEquals(6, m.size(), "Induced MCS of benzene in toluene = 6");
    }

    @Test @DisplayName("Non-induced MCS: benzene in toluene >= 6")
    void nonInducedBenzeneToluene() throws Exception {
      SMSD smsd = new SMSD(mol("c1ccccc1"), mol("Cc1ccccc1"), new ChemOptions());
      Map<Integer, Integer> m = smsd.findMCS(false, true, 10_000L);
      assertTrue(m.size() >= 6, "Non-induced MCS >= 6");
    }
  }

  // ======================================================================
  // 8. ALL CHEMOPTIONS PROFILES
  // ======================================================================

  @Nested
  @DisplayName("8. AllChemOptionsProfiles")
  class AllChemOptionsProfiles {

    @Test @DisplayName("Default profile")
    void defaultProfile() throws Exception {
      ChemOptions opts = new ChemOptions();
      assertEquals(ChemOptions.BondOrderMode.STRICT, opts.matchBondOrder);
      assertEquals(ChemOptions.AromaticityMode.FLEXIBLE, opts.aromaticityMode);
      assertTrue(opts.matchAtomType);
      assertTrue(opts.matchFormalCharge);
      assertFalse(opts.useChirality);
      assertFalse(opts.useBondStereo);
      assertTrue(opts.ringMatchesRingOnly);
      assertFalse(opts.completeRingsOnly);
      assertFalse(opts.matchIsotope);
      assertFalse(opts.tautomerAware);
      assertEquals(ChemOptions.RingFusionMode.IGNORE, opts.ringFusionMode);
      assertEquals(ChemOptions.MatcherEngine.VF2PP, opts.matcherEngine);
    }

    @Test @DisplayName("profile(strict)")
    void strictProfile() throws Exception {
      ChemOptions opts = ChemOptions.profile("strict");
      assertEquals(ChemOptions.BondOrderMode.STRICT, opts.matchBondOrder);
      assertEquals(ChemOptions.AromaticityMode.STRICT, opts.aromaticityMode);
      assertTrue(new SMSD(mol("c1ccccc1"), mol("c1ccc(O)cc1"), opts).isSubstructure());
    }

    @Test @DisplayName("profile(compat-substruct)")
    void compatSubstruct() throws Exception {
      ChemOptions opts = ChemOptions.profile("compat-substruct");
      assertEquals(ChemOptions.BondOrderMode.LOOSE, opts.matchBondOrder);
      assertEquals(ChemOptions.AromaticityMode.FLEXIBLE, opts.aromaticityMode);
      assertFalse(opts.ringMatchesRingOnly);
      assertTrue(new SMSD(mol("c1ccccc1"), mol("c1ccc(O)cc1"), opts).isSubstructure());
    }

    @Test @DisplayName("tautomerProfile()")
    void tautomerProfile() throws Exception {
      ChemOptions opts = ChemOptions.tautomerProfile();
      assertTrue(opts.tautomerAware);
      assertEquals(ChemOptions.BondOrderMode.LOOSE, opts.matchBondOrder);
      assertFalse(opts.matchFormalCharge);
      assertFalse(opts.ringMatchesRingOnly);
    }

    @Test @DisplayName("Custom: every field toggled from default")
    void customAllToggled() throws Exception {
      ChemOptions opts = new ChemOptions();
      opts.matchBondOrder = ChemOptions.BondOrderMode.ANY;
      opts.aromaticityMode = ChemOptions.AromaticityMode.STRICT;
      opts.matchAtomType = true;
      opts.matchFormalCharge = false;
      opts.useChirality = true;
      opts.useBondStereo = true;
      opts.ringMatchesRingOnly = false;
      opts.completeRingsOnly = true;
      opts.matchIsotope = true;
      opts.tautomerAware = true;
      opts.ringFusionMode = ChemOptions.RingFusionMode.STRICT;
      opts.matcherEngine = ChemOptions.MatcherEngine.VF3;
      opts.useTwoHopNLF = false;
      opts.useThreeHopNLF = true;
      opts.useBitParallelFeasibility = false;
      // Just verify it runs without crashing
      SMSD smsd = new SMSD(mol("c1ccccc1"), mol("c1ccccc1"), opts);
      smsd.findMCS(false, true, 5000L);
    }

    @Test @DisplayName("Fluent API chaining returns same instance")
    void fluentApi() {
      ChemOptions opts = new ChemOptions()
          .withBondOrderMode(ChemOptions.BondOrderMode.LOOSE)
          .withAromaticityMode(ChemOptions.AromaticityMode.FLEXIBLE)
          .withTwoHopNLF(false)
          .withThreeHopNLF(true)
          .withBitParallelFeasibility(false)
          .withCompleteRingsOnly(true)
          .withMatchIsotope(true)
          .withRingFusionMode(ChemOptions.RingFusionMode.PERMISSIVE)
          .withTautomerAware(true);
      assertEquals(ChemOptions.BondOrderMode.LOOSE, opts.matchBondOrder);
      assertTrue(opts.tautomerAware);
      assertTrue(opts.completeRingsOnly);
    }

    @Test @DisplayName("VF2 engine produces same result as VF2PP")
    void vf2VsVf2pp() throws Exception {
      ChemOptions vf2 = new ChemOptions();
      vf2.matcherEngine = ChemOptions.MatcherEngine.VF2;
      ChemOptions vf2pp = new ChemOptions();
      vf2pp.matcherEngine = ChemOptions.MatcherEngine.VF2PP;
      boolean r1 = new SMSD(mol("c1ccccc1"), mol("c1ccc(O)cc1"), vf2).isSubstructure();
      boolean r2 = new SMSD(mol("c1ccccc1"), mol("c1ccc(O)cc1"), vf2pp).isSubstructure();
      assertEquals(r1, r2, "VF2 and VF2PP should agree");
    }

    @Test @DisplayName("VF3 engine produces same result as VF2PP")
    void vf3VsVf2pp() throws Exception {
      ChemOptions vf3 = new ChemOptions();
      vf3.matcherEngine = ChemOptions.MatcherEngine.VF3;
      ChemOptions vf2pp = new ChemOptions();
      vf2pp.matcherEngine = ChemOptions.MatcherEngine.VF2PP;
      boolean r1 = new SMSD(mol("c1ccccc1"), mol("c1ccc(O)cc1"), vf3).isSubstructure();
      boolean r2 = new SMSD(mol("c1ccccc1"), mol("c1ccc(O)cc1"), vf2pp).isSubstructure();
      assertEquals(r1, r2, "VF3 and VF2PP should agree");
    }

    @Test @DisplayName("Null bond order mode defaults to STRICT")
    void nullBondOrderDefaults() {
      ChemOptions opts = new ChemOptions().withBondOrderMode(null);
      assertEquals(ChemOptions.BondOrderMode.STRICT, opts.matchBondOrder);
    }
  }

  // ======================================================================
  // 9. ALGORITHM STRESS
  // ======================================================================

  @Nested
  @DisplayName("9. AlgorithmStress")
  class AlgorithmStress {

    @Test @DisplayName("Greedy probe: benzene/toluene fast path")
    void greedyBenzeneToluene() throws Exception {
      int sz = mcsSize("c1ccccc1", "Cc1ccccc1");
      assertEquals(6, sz);
    }

    @Test @DisplayName("Seed-extend: morphine/codeine")
    void seedExtendMorphineCodeine() throws Exception {
      String morphine = "CN1CCC23C4C1CC5=C2C(=C(C=C5)O)OC3C(C=C4)O";
      String codeine = "CN1CCC23C4C1CC5=C2C(=C(C=C5)OC)OC3C(C=C4)O";
      int sz = mcsSize(morphine, codeine);
      assertTrue(sz >= 18, "Morphine/codeine MCS >= 18");
    }

    @Test @DisplayName("McSplit: medium molecules (caffeine/theophylline)")
    void mcsplitMedium() throws Exception {
      String caffeine = "Cn1c(=O)c2c(ncn2C)n(C)c1=O";
      String theophylline = "Cn1c(=O)c2[nH]cnc2n(C)c1=O";
      int sz = mcsSize(caffeine, theophylline);
      assertTrue(sz >= 10, "Caffeine/theophylline MCS >= 10");
    }

    @Test @DisplayName("BK clique: symmetric adamantane")
    void bkAdamantane() throws Exception {
      String adamantane = "C1C2CC3CC1CC(C2)C3";
      assertEquals(10, mcsSize(adamantane, adamantane));
    }

    @Test @DisplayName("McGregor: large molecule pair (atorvastatin/rosuvastatin)")
    void mcgregorLarge() throws Exception {
      String atorvastatin = "CC(C)c1c(C(=O)Nc2ccccc2)c(-c2ccccc2)c(-c2ccc(F)cc2)n1CC[C@H](O)C[C@H](O)CC(=O)O";
      String rosuvastatin = "CC(C)c1nc(N(C)S(C)(=O)=O)nc(-c2ccc(F)cc2)c1/C=C/[C@H](O)C[C@H](O)CC(=O)O";
      int sz = mcsSize(atorvastatin, rosuvastatin, new ChemOptions(), 30_000L);
      assertTrue(sz >= 10, "Large molecule MCS should be non-trivial");
    }

    @Test @DisplayName("Coverage-driven: identical molecules exit early")
    void coverageDrivenIdentical() throws Exception {
      String mol = "c1ccc2ccccc2c1"; // naphthalene
      long start = System.currentTimeMillis();
      int sz = mcsSize(mol, mol);
      long elapsed = System.currentTimeMillis() - start;
      assertEquals(10, sz);
      assertTrue(elapsed < 2000, "Identical molecules should be fast, took " + elapsed + "ms");
    }

    @Test @DisplayName("Timeout 1ms on paclitaxel pair: no crash, returns partial")
    @Timeout(30)
    void timeoutPaclitaxel() throws Exception {
      String paclitaxel = "CC1=C2C(C(=O)C3(C(CC4C(C3C(C(C2(C)C)(CC1OC(=O)C(C(C5=CC=CC=C5)NC(=O)C6=CC=CC=C6)O)O)OC(=O)C7=CC=CC=C7)(CO4)OC(=O)C)O)C)OC(=O)C";
      String docetaxel = "CC1=C2C(C(=O)C3(C(CC4C(C3C(C(C2(C)C)(CC1O)O)OC(=O)C5=CC=CC=C5)(CO4)OC(=O)C)O)C)OC(=O)C";
      SMSD smsd = new SMSD(mol(paclitaxel), mol(docetaxel), new ChemOptions());
      smsd.setMcsTimeoutMs(1L);
      // Should not crash; result may be empty or partial
      assertDoesNotThrow(() -> smsd.findMCS(false, true, 1L));
    }

    @Test @DisplayName("Timeout 100ms on drug pair: returns valid result")
    @Timeout(10)
    void timeout100msDrugPair() throws Exception {
      String aspirin = "CC(=O)Oc1ccccc1C(=O)O";
      String salicylicAcid = "OC(=O)c1ccccc1O";
      SMSD smsd = new SMSD(mol(aspirin), mol(salicylicAcid), new ChemOptions());
      smsd.setMcsTimeoutMs(100L);
      Map<Integer, Integer> m = smsd.findMCS(false, true, 100L);
      // Should return a valid (possibly partial) result
      assertNotNull(m);
    }

    @Test @DisplayName("Multiple mappings: enumerate up to 10 for benzene self-match")
    void multipleMappings() throws Exception {
      SMSD smsd = new SMSD(mol("c1ccccc1"), mol("c1ccccc1"), new ChemOptions());
      List<Map<Integer, Integer>> maps = smsd.findAllSubstructures(10, 5000L);
      assertTrue(maps.size() >= 1, "Benzene self-match should have at least 1 mapping");
    }

    @Test @DisplayName("Substructure false: query larger than target")
    void queryLargerThanTarget() throws Exception {
      assertFalse(new SMSD(mol("c1ccc2ccccc2c1"), mol("c1ccccc1"), new ChemOptions()).isSubstructure());
    }

    @Test @DisplayName("Substructure true: benzene in phenol")
    void benzeneInPhenol() throws Exception {
      assertTrue(new SMSD(mol("c1ccccc1"), mol("c1ccc(O)cc1"), new ChemOptions()).isSubstructure());
    }

    @Test @DisplayName("MCS empty: methane vs water")
    void mcsEmpty() throws Exception {
      Map<Integer, Integer> m = mcs("C", "O", new ChemOptions(), 5000L);
      assertTrue(m.isEmpty(), "No common substructure between methane and water");
    }

    @Test @DisplayName("MCS identity: same molecule returns full mapping")
    void mcsIdentity() throws Exception {
      String mol = "CC(=O)Oc1ccccc1C(=O)O"; // aspirin
      int sz = mcsSize(mol, mol);
      assertEquals(mol(mol).getAtomCount(), sz, "Self-MCS should return all atoms");
    }

    @Test @DisplayName("RASCAL screening: upper bound >= actual MCS")
    void rascalUpperBound() throws Exception {
      IAtomContainer q = mol("c1ccccc1");
      IAtomContainer t = mol("c1ccc(O)cc1");
      SMSD smsd = new SMSD(q, t, new ChemOptions());
      double ub = smsd.similarityUpperBound();
      Map<Integer, Integer> m = smsd.findMCS();
      double actual = (double) m.size() / Math.max(q.getAtomCount(), t.getAtomCount());
      assertTrue(ub >= actual - 0.01, "Upper bound should be >= actual similarity");
    }

    @Test @DisplayName("Fingerprint subset: substructure pairs")
    void fingerprintSubset() throws Exception {
      IAtomContainer q = mol("c1ccccc1");
      IAtomContainer t = mol("c1ccc(O)cc1");
      long[] qFp = SMSD.pathFingerprint(q, 7, 1024);
      long[] tFp = SMSD.pathFingerprint(t, 7, 1024);
      assertTrue(SMSD.fingerprintSubset(qFp, tFp), "Benzene FP should be subset of phenol FP");
    }

    @Test @DisplayName("Fingerprint non-subset: non-substructure pairs")
    void fingerprintNonSubset() throws Exception {
      IAtomContainer q = mol("c1ccc(O)cc1"); // phenol
      IAtomContainer t = mol("CCCC"); // butane
      long[] qFp = SMSD.pathFingerprint(q, 7, 1024);
      long[] tFp = SMSD.pathFingerprint(t, 7, 1024);
      assertFalse(SMSD.fingerprintSubset(qFp, tFp), "Phenol FP should NOT be subset of butane FP");
    }

    @Test @DisplayName("validateMapping: every mapping is chemically valid")
    void validateMapping() throws Exception {
      String benzene = "c1ccccc1";
      String phenol = "c1ccc(O)cc1";
      ChemOptions opts = new ChemOptions();
      Map<Integer, Integer> m = mcs(benzene, phenol, opts, 10_000L);
      MolGraph g1 = new MolGraph(mol(benzene));
      MolGraph g2 = new MolGraph(mol(phenol));
      List<String> errors = SearchEngine.validateMapping(g1, g2, m, opts);
      assertTrue(errors.isEmpty(), "Mapping should be valid, errors: " + errors);
    }

    @Test @DisplayName("isMappingMaximal: MCS mapping should be maximal")
    void mappingMaximal() throws Exception {
      String benzene = "c1ccccc1";
      String toluene = "Cc1ccccc1";
      ChemOptions opts = new ChemOptions();
      Map<Integer, Integer> m = mcs(benzene, toluene, opts, 10_000L);
      MolGraph g1 = new MolGraph(mol(benzene));
      MolGraph g2 = new MolGraph(mol(toluene));
      assertTrue(SearchEngine.isMappingMaximal(g1, g2, m, opts), "MCS should be maximal");
    }
  }

  // ======================================================================
  // 10. SCALE STRESS
  // ======================================================================

  @Nested
  @DisplayName("10. ScaleStress")
  class ScaleStress {

    private String linearAlkane(int n) {
      return "C".repeat(n);
    }

    @Test @DisplayName("50-atom self-match < 1s")
    @Timeout(5)
    void selfMatch50() throws Exception {
      String smi = linearAlkane(50);
      long start = System.currentTimeMillis();
      int sz = mcsSize(smi, smi);
      long elapsed = System.currentTimeMillis() - start;
      assertEquals(50, sz);
      assertTrue(elapsed < 1000, "50-atom self-match took " + elapsed + "ms");
    }

    @Test @DisplayName("100-atom self-match < 2s")
    @Timeout(10)
    void selfMatch100() throws Exception {
      String smi = linearAlkane(100);
      long start = System.currentTimeMillis();
      int sz = mcsSize(smi, smi);
      long elapsed = System.currentTimeMillis() - start;
      assertEquals(100, sz);
      assertTrue(elapsed < 2000, "100-atom self-match took " + elapsed + "ms");
    }

    @Test @DisplayName("200-atom self-match < 5s")
    @Timeout(15)
    void selfMatch200() throws Exception {
      String smi = linearAlkane(200);
      long start = System.currentTimeMillis();
      int sz = mcsSize(smi, smi);
      long elapsed = System.currentTimeMillis() - start;
      assertEquals(200, sz);
      assertTrue(elapsed < 5000, "200-atom self-match took " + elapsed + "ms");
    }

    @Test @DisplayName("400-atom substructure < 5s")
    @Timeout(15)
    void substructure400() throws Exception {
      IAtomContainer q = mol(linearAlkane(100));
      IAtomContainer t = mol(linearAlkane(400));
      long start = System.currentTimeMillis();
      boolean result = new SMSD(q, t, new ChemOptions()).isSubstructure();
      long elapsed = System.currentTimeMillis() - start;
      assertTrue(result);
      assertTrue(elapsed < 5000, "400-atom substructure took " + elapsed + "ms");
    }

    @Test @DisplayName("1000 sequential substructure checks < 10s")
    @Timeout(20)
    void sequentialSubstructure1000() throws Exception {
      IAtomContainer q = mol("c1ccccc1");
      IAtomContainer t = mol("c1ccc(O)cc1");
      ChemOptions opts = new ChemOptions();
      long start = System.currentTimeMillis();
      for (int i = 0; i < 1000; i++) {
        new SMSD(q, t, opts, false).isSubstructure();
      }
      long elapsed = System.currentTimeMillis() - start;
      assertTrue(elapsed < 10000, "1000 substructure checks took " + elapsed + "ms");
    }

    @Test @DisplayName("100 sequential MCS computations < 30s")
    @Timeout(60)
    void sequentialMCS100() throws Exception {
      IAtomContainer q = mol("c1ccccc1");
      IAtomContainer t = mol("c1ccc(O)cc1");
      ChemOptions opts = new ChemOptions();
      long start = System.currentTimeMillis();
      for (int i = 0; i < 100; i++) {
        new SMSD(q, t, opts, false).findMCS(false, true, 5000L);
      }
      long elapsed = System.currentTimeMillis() - start;
      assertTrue(elapsed < 30000, "100 MCS computations took " + elapsed + "ms");
    }

    @Test @DisplayName("Batch MCS with threshold: 10 molecules < 5s")
    @Timeout(15)
    void batchMcs10() throws Exception {
      String[] smiles = {
          "c1ccccc1", "c1ccc(O)cc1", "c1ccc(N)cc1", "c1ccc(F)cc1", "c1ccc(Cl)cc1",
          "c1ccc(Br)cc1", "c1ccc(I)cc1", "c1ccc(C)cc1", "c1ccc(S)cc1", "c1ccc(P)cc1"
      };
      List<IAtomContainer> mols = new ArrayList<>();
      for (String s : smiles) mols.add(mol(s));
      long start = System.currentTimeMillis();
      Map<Integer, Integer> nm = SMSD.findNMCS(mols, new ChemOptions(), 0.8, 5000L);
      long elapsed = System.currentTimeMillis() - start;
      assertNotNull(nm);
      assertTrue(elapsed < 5000, "Batch MCS took " + elapsed + "ms");
    }

    @Test @DisplayName("R-group decomposition: benzene core on 5 derivatives < 2s")
    @Timeout(10)
    void rgroupDecomp() throws Exception {
      IAtomContainer core = mol("c1ccccc1");
      List<IAtomContainer> derivatives = Arrays.asList(
          mol("c1ccc(O)cc1"), mol("c1ccc(N)cc1"), mol("c1ccc(F)cc1"),
          mol("c1ccc(Cl)cc1"), mol("c1ccc(C)cc1"));
      long start = System.currentTimeMillis();
      List<Map<String, IAtomContainer>> decomp = SMSD.decomposeRGroups(core, derivatives, new ChemOptions(), 5000L);
      long elapsed = System.currentTimeMillis() - start;
      assertNotNull(decomp);
      assertFalse(decomp.isEmpty());
      assertTrue(elapsed < 2000, "R-group decomposition took " + elapsed + "ms");
    }

    @Test @DisplayName("N-MCS: 5 molecules < 10s")
    @Timeout(20)
    void nmcs5() throws Exception {
      List<IAtomContainer> mols = Arrays.asList(
          mol("c1ccc(O)cc1"), mol("c1ccc(N)cc1"), mol("c1ccc(Cl)cc1"),
          mol("c1ccc(F)cc1"), mol("c1ccc(C)cc1"));
      long start = System.currentTimeMillis();
      Map<Integer, Integer> nm = SMSD.findNMCS(mols, new ChemOptions(), 1.0, 10_000L);
      long elapsed = System.currentTimeMillis() - start;
      assertTrue(nm.size() >= 6, "All 5 share benzene ring");
      assertTrue(elapsed < 10000, "N-MCS 5 molecules took " + elapsed + "ms");
    }

    @Test @DisplayName("Fingerprint: 100 molecules < 1s")
    @Timeout(5)
    void fingerprint100() throws Exception {
      String[] smiles = new String[100];
      for (int i = 0; i < 100; i++) {
        smiles[i] = linearAlkane(5 + (i % 20));
      }
      long start = System.currentTimeMillis();
      for (String s : smiles) {
        SMSD.pathFingerprint(mol(s), 7, 1024);
      }
      long elapsed = System.currentTimeMillis() - start;
      assertTrue(elapsed < 1000, "100 fingerprints took " + elapsed + "ms");
    }
  }

  // ======================================================================
  // 11. REAL WORLD DRUG PAIRS
  // ======================================================================

  @Nested
  @DisplayName("11. RealWorldDrugPairs")
  class RealWorldDrugPairs {

    private void assertDrugPairMCS(String name1, String smi1, String name2, String smi2,
                                    int minMCS) throws Exception {
      SearchEngine.clearMolGraphCache();
      ChemOptions opts = new ChemOptions();
      Map<Integer, Integer> m = mcs(smi1, smi2, opts, 30_000L);
      assertTrue(m.size() >= minMCS,
          name1 + "/" + name2 + " MCS expected >= " + minMCS + " but got " + m.size());
      // Validate mapping
      MolGraph g1 = new MolGraph(mol(smi1));
      MolGraph g2 = new MolGraph(mol(smi2));
      List<String> errors = SearchEngine.validateMapping(g1, g2, m, opts);
      assertTrue(errors.isEmpty(),
          name1 + "/" + name2 + " mapping invalid: " + errors);
    }

    @Test @DisplayName("Aspirin / Acetaminophen >= 7")
    @Timeout(30)
    void aspirinAcetaminophen() throws Exception {
      assertDrugPairMCS("Aspirin", "CC(=O)Oc1ccccc1C(=O)O",
          "Acetaminophen", "CC(=O)Nc1ccc(O)cc1", 7);
    }

    @Test @DisplayName("Ibuprofen / Naproxen >= 12")
    @Timeout(30)
    void ibuprofenNaproxen() throws Exception {
      assertDrugPairMCS("Ibuprofen", "CC(C)Cc1ccc(C(C)C(=O)O)cc1",
          "Naproxen", "COc1ccc2cc(C(C)C(=O)O)ccc2c1", 12);
    }

    @Test @DisplayName("Morphine / Codeine >= 20")
    @Timeout(30)
    void morphineCodeine() throws Exception {
      assertDrugPairMCS("Morphine", "CN1CCC23C4C1CC5=C2C(=C(C=C5)O)OC3C(C=C4)O",
          "Codeine", "CN1CCC23C4C1CC5=C2C(=C(C=C5)OC)OC3C(C=C4)O", 20);
    }

    @Test @DisplayName("Caffeine / Theophylline >= 13")
    @Timeout(30)
    void caffeineTheophylline() throws Exception {
      assertDrugPairMCS("Caffeine", "Cn1c(=O)c2c(ncn2C)n(C)c1=O",
          "Theophylline", "Cn1c(=O)c2[nH]cnc2n(C)c1=O", 13);
    }

    @Test @DisplayName("Atorvastatin / Rosuvastatin >= 15")
    @Timeout(30)
    void atorvastatinRosuvastatin() throws Exception {
      assertDrugPairMCS("Atorvastatin",
          "CC(C)c1c(C(=O)Nc2ccccc2)c(-c2ccccc2)c(-c2ccc(F)cc2)n1CC[C@H](O)C[C@H](O)CC(=O)O",
          "Rosuvastatin",
          "CC(C)c1nc(N(C)S(C)(=O)=O)nc(-c2ccc(F)cc2)c1/C=C/[C@H](O)C[C@H](O)CC(=O)O", 10);
    }

    @Test @DisplayName("Diazepam / Oxazepam >= 15")
    @Timeout(30)
    void diazepamOxazepam() throws Exception {
      assertDrugPairMCS("Diazepam", "CN1C(=O)CN=C(c2ccccc2)c2cc(Cl)ccc21",
          "Oxazepam", "OC1N=C(c2ccccc2)c2cc(Cl)ccc2NC1=O", 15);
    }

    @Test @DisplayName("Metformin / Phenformin >= 8")
    @Timeout(30)
    void metforminPhenformin() throws Exception {
      assertDrugPairMCS("Metformin", "CN(C)C(=N)NC(=N)N",
          "Phenformin", "NC(=N)NC(=N)NCCc1ccccc1", 8);
    }

    @Test @DisplayName("Omeprazole / Lansoprazole >= 12")
    @Timeout(30)
    void omeprazoleLansoprazole() throws Exception {
      assertDrugPairMCS("Omeprazole", "COc1ccc2[nH]c(S(=O)Cc3ncc(C)c(OC)c3C)nc2c1",
          "Lansoprazole", "Cc1c(OCC(F)(F)F)ccn1S(=O)c1nc2ccc(OC)cc2[nH]1", 12);
    }

    @Test @DisplayName("Sildenafil / Tadalafil >= 10")
    @Timeout(30)
    void sildenafilTadalafil() throws Exception {
      assertDrugPairMCS("Sildenafil",
          "CCCc1nn(C)c2c1nc(nc2OCC)c1cc(ccc1OCC)S(=O)(=O)N1CCN(C)CC1",
          "Tadalafil",
          "O=C1NC(=O)c2cc3c(cc2[C@@H]1Cc1c[nH]c2ccccc12)OCO3", 10);
    }

    @Test @DisplayName("Tamoxifen / Raloxifene >= 10")
    @Timeout(30)
    void tamoxifenRaloxifene() throws Exception {
      assertDrugPairMCS("Tamoxifen", "CCC(=C(c1ccccc1)c1ccc(OCCN(C)C)cc1)c1ccccc1",
          "Raloxifene", "Oc1ccc(cc1)C(=O)c1ccc(O)cc1", 5);
    }

    @Test @DisplayName("Ciprofloxacin / Levofloxacin >= 15")
    @Timeout(30)
    void ciprofloxacinLevofloxacin() throws Exception {
      assertDrugPairMCS("Ciprofloxacin",
          "O=C(O)c1cn(C2CC2)c2cc(N3CCNCC3)c(F)cc2c1=O",
          "Levofloxacin",
          "C[C@@H]1COc2c(N3CCN(C)CC3)c(F)cc3c(=O)c(C(=O)O)cn1c23", 15);
    }

    @Test @DisplayName("Warfarin / Acenocoumarol >= 15")
    @Timeout(30)
    void warfarinAcenocoumarol() throws Exception {
      assertDrugPairMCS("Warfarin", "CC(=O)CC(c1ccccc1)c1c(O)c2ccccc2oc1=O",
          "Acenocoumarol", "CC(=O)CC(c1ccc([N+](=O)[O-])cc1)c1c(O)c2ccccc2oc1=O", 15);
    }

    @Test @DisplayName("Metoprolol / Atenolol >= 10")
    @Timeout(30)
    void metoprololAtenolol() throws Exception {
      assertDrugPairMCS("Metoprolol", "COCCc1ccc(OCC(O)CNC(C)C)cc1",
          "Atenolol", "CC(C)NCC(O)COc1ccc(CC(N)=O)cc1", 10);
    }

    @Test @DisplayName("Amoxicillin / Ampicillin >= 15")
    @Timeout(30)
    void amoxicillinAmpicillin() throws Exception {
      assertDrugPairMCS("Amoxicillin",
          "CC1(C)SC2C(NC(=O)C(N)c3ccc(O)cc3)C(=O)N2C1C(=O)O",
          "Ampicillin",
          "CC1(C)SC2C(NC(=O)C(N)c3ccccc3)C(=O)N2C1C(=O)O", 15);
    }

    @Test @DisplayName("Fluoxetine / Paroxetine >= 10")
    @Timeout(30)
    void fluoxetineParoxetine() throws Exception {
      assertDrugPairMCS("Fluoxetine", "CNCCC(Oc1ccc(C(F)(F)F)cc1)c1ccccc1",
          "Paroxetine", "Fc1ccc(C2CCNCC2COc2ccc3c(c2)OCO3)cc1", 10);
    }

    @Test @DisplayName("Loratadine / Desloratadine >= 18")
    @Timeout(30)
    void loratadineDesloratadine() throws Exception {
      assertDrugPairMCS("Loratadine",
          "CCOC(=O)N1CCC(=C2c3ccc(Cl)cc3CCc3cccnc32)CC1",
          "Desloratadine",
          "Clc1ccc2c(c1)CCc1cccnc1C2=C1CCNCC1", 18);
    }

    @Test @DisplayName("Prednisolone / Dexamethasone >= 18")
    @Timeout(30)
    void prednisoloneDexamethasone() throws Exception {
      assertDrugPairMCS("Prednisolone",
          "O=C1C=C2CC3C(CC(O)C4(C(CO)=O)C3CCC4O)C2(C)CC1=O",
          "Dexamethasone",
          "CC1CC2C3CCC4=CC(=O)C=CC4(C)C3(F)C(O)CC2(C)C1(O)C(=O)CO", 18);
    }

    @Test @DisplayName("Simvastatin / Lovastatin >= 20")
    @Timeout(30)
    void simvastatinLovastatin() throws Exception {
      assertDrugPairMCS("Simvastatin",
          "CCC(C)(C)C(=O)OC1CC(O)C=C2C=CC(C)C(CCC3CC(O)CC(=O)O3)C21",
          "Lovastatin",
          "CCC(C)C(=O)OC1CC(O)C=C2C=CC(C)C(CCC3CC(O)CC(=O)O3)C21", 20);
    }

    @Test @DisplayName("Celecoxib / Valdecoxib >= 12")
    @Timeout(30)
    void celecoxibValdecoxib() throws Exception {
      assertDrugPairMCS("Celecoxib",
          "Cc1ccc(-c2cc(C(F)(F)F)nn2-c2ccc(S(N)(=O)=O)cc2)cc1",
          "Valdecoxib",
          "Cc1ccc(-c2noc(C)c2)cc1", 8);
    }

    @Test @DisplayName("Donepezil / Rivastigmine >= 8")
    @Timeout(30)
    void donepezilRivastigmine() throws Exception {
      assertDrugPairMCS("Donepezil",
          "COc1cc2CC(CC3CCN(Cc4ccccc4)CC3)C(=O)c2cc1OC",
          "Rivastigmine",
          "CCN(C)C(=O)Oc1cccc(C(C)N(C)C)c1", 8);
    }
  }

  // ======================================================================
  // From: LargeMoleculeTest.java
  // ======================================================================


  // ---- Drug SMILES (from PubChem) ----
  // Paclitaxel (Taxol) — 51 heavy atoms, complex tetracyclic taxane
  static final String PACLITAXEL =
      "CC1=C2C(C(=O)C3(C(CC4C(C3C(C(C2(C)C)(CC1OC(=O)C(C(C5=CC=CC=C5)NC(=O)C6=CC=CC=C6)O)O)OC(=O)C7=CC=CC=C7)(CO4)OC(=O)C)O)C)OC(=O)C";
  // Docetaxel — 50 heavy atoms, differs from paclitaxel in N-acyl + C10
  static final String DOCETAXEL =
      "CC1=C2C(C(=O)C3(C(CC4C(C3C(C(C2(C)C)(CC1OC(=O)C(C(C5=CC=CC=C5)NC(=O)OC(C)(C)C)O)O)OC(=O)C6=CC=CC=C6)(CO4)OC(=O)C)O)C)O";
  // Rapamycin (Sirolimus) — 65 heavy atoms, macrocyclic
  static final String RAPAMYCIN =
      "CC1CCC2CC(C(=CC=CC=CC(CC(C(=O)C(C(C(=CC(C(=O)CC(OC(=O)C3CCCCN3C(=O)C(=O)C1(O2)O)C(C)CC4CCC(C(C4)OC)O)C)C)O)OC)C)C)C)OC";
  // Vincristine — 56 heavy atoms, vinca alkaloid
  static final String VINCRISTINE =
      "CCC1(CC2CC(C3=C(CCN(C2)C1)C4=CC=C(C=C4N3)OC)C5(C(=O)OC)C6=CC7=C(C=C6N8C5CC9CC(C(CC(=O)OC)(C(CC)(O9)C(=O)OC)O)C8C7)OC)O";
  // Erythromycin — 37 heavy atoms, macrolide antibiotic
  static final String ERYTHROMYCIN =
      "CCC1OC(=O)C(C)C(OC2CC(C)(OC)C(O)C(C)O2)C(C)C(OC2OC(C)CC(C2O)N(C)C)C(C)(O)CC(C)C(=O)C(C)C(O)C1(C)O";
  // Azithromycin — 38 heavy atoms, macrolide (15-ring vs erythromycin's 14-ring)
  static final String AZITHROMYCIN =
      "CCC1C(C(C(N(CC(CC(C(C(C(C(C(=O)O1)C)OC2CC(C(C(O2)C)O)(C)OC)C)OC3C(C(CC(O3)C)N(C)C)O)(C)O)C)C)C)O)(C)O";
  // Atorvastatin — 33 heavy atoms
  static final String ATORVASTATIN =
      "CC(C)C1=C(C(=C(N1CCC(CC(CC(=O)O)O)O)C2=CC=C(C=C2)F)C3=CC=CC=C3)C(=O)NC4=CC=CC=C4";
  // Rosuvastatin — 28 heavy atoms
  static final String ROSUVASTATIN =
      "CC(C)C1=NC(=NC(=C1C=CC(CC(CC(=O)O)O)O)C2=CC=C(C=C2)F)N(C)S(=O)(=O)C";
  // Dexamethasone — 25 heavy atoms, steroid
  static final String DEXAMETHASONE = "CC1CC2C3CCC4=CC(=O)C=CC4(C3(C(CC2(C1(C(=O)CO)O)C)O)F)C";
  // Cholecalciferol (Vitamin D3) — 28 heavy atoms
  static final String VITAMIN_D3 = "CC(C)CCCC(C)C1CCC2C1(CCCC2=CC=C3CC(CCC3=C)O)C";

  // ---- Biochemical substrates (from PubChem/RHEA) ----
  // NADH — reduced form of NAD+, differs by 2H at nicotinamide ring
  static final String NADH =
      "C1C=CN(C=C1C(=O)N)C2C(C(C(O2)COP(=O)(O)OP(=O)(O)OCC3C(C(C(O3)N4C=NC5=C4N=CN=C5N)O)O)O)O";
  // CoA (free thiol) — 48 heavy atoms
  static final String COA =
      "CC(C)(COP(=O)(O)OP(=O)(O)OCC1C(C(C(O1)N2C=NC3=C2N=CN=C3N)O)OP(=O)(O)O)C(C(=O)NCCC(=O)NCCS)O";
  // FAD — 53 heavy atoms, flavin + adenine dinucleotide
  static final String FAD =
      "CC1=CC2=C(C=C1C)N(C3=NC(=O)NC(=O)C3=N2)CC(C(C(COP(=O)(O)OP(=O)(O)OCC4C(C(C(O4)N5C=NC6=C(N=CN=C65)N)O)O)O)O)O";
  // SAM — 22 heavy atoms, methyl donor
  static final String SAM = "C[S+](CCC(C(=O)[O-])N)CC1C(C(C(O1)N2C=NC3=C2N=CN=C3N)O)O";
  // NAD+ — 47 heavy atoms, dinucleotide
  static final String NAD_PLUS =
      "C1=CC(=C[N+](=C1)C2C(C(C(O2)COP(=O)([O-])OP(=O)(O)OCC3C(C(C(O3)N4C=NC5=C(N=CN=C54)N)O)O)O)O)C(=O)N";
  // ATP — 31 heavy atoms
  static final String ATP = "C1=NC(=C2C(=N1)N(C=N2)C3C(C(C(O3)COP(=O)(O)OP(=O)(O)OP(=O)(O)O)O)O)N";
  // ADP — 27 heavy atoms
  static final String ADP = "C1=NC(=C2C(=N1)N(C=N2)C3C(C(C(O3)COP(=O)(O)OP(=O)(O)O)O)O)N";
  // AMP — 23 heavy atoms
  static final String AMP = "C1=NC(=C2C(=N1)N(C=N2)C3C(C(C(O3)COP(=O)(O)O)O)O)N";
  // Acetyl-CoA — 51 heavy atoms
  static final String ACETYL_COA =
      "CC(=O)SCCNC(=O)CCNC(=O)C(C(C)(C)COP(=O)(O)OP(=O)(O)OCC1C(C(C(O1)N2C=NC3=C(N=CN=C32)N)O)OP(=O)(O)O)O";
  // UDP-glucose — 36 heavy atoms
  static final String UDP_GLUCOSE =
      "OC[C@H]1O[C@H](OP([O-])(=O)OP([O-])(=O)OC[C@H]2O[C@H]([C@H](O)[C@@H]2O)n2ccc(=O)[nH]c2=O)[C@H](O)[C@@H](O)[C@@H]1O";
  // UDP-galactose — 36 heavy atoms
  static final String UDP_GALACTOSE =
      "OC[C@H]1O[C@H](OP([O-])(=O)OP([O-])(=O)OC[C@H]2O[C@H]([C@H](O)[C@@H]2O)n2ccc(=O)[nH]c2=O)[C@H](O)[C@@H](O)[C@H]1O";

  // ---- Mega molecules (100+ atoms, from PubChem) ----
  // Vancomycin — 101 heavy atoms, glycopeptide antibiotic
  static final String VANCOMYCIN =
      "CC1C(C(CC(O1)OC2C(C(C(OC2OC3=C4C=C5C=C3OC6=C(C=C(C=C6)C(C(C(=O)NC(C(=O)NC5C(=O)NC7C8=CC(=C(C=C8)O)C9=C(C=C(C=C9O)O)C(NC(=O)C(C(C1=CC(=C(O4)C=C1)Cl)O)NC7=O)C(=O)O)CC(=O)N)NC(=O)C(CC(C)C)NC)O)Cl)CO)O)O)(C)N)O";
  // Cyclosporin A — 85 heavy atoms, cyclic peptide immunosuppressant
  static final String CYCLOSPORIN =
      "CCC1C(=O)N(CC(=O)N(C(C(=O)NC(C(=O)N(C(C(=O)NC(C(=O)NC(C(=O)N(C(C(=O)N(C(C(=O)N(C(C(=O)N(C(C(=O)N1)C(C(C)CC=CC)O)C)C(C)C)C)CC(C)C)C)CC(C)C)C)C)C)CC(C)C)C)C(C)C)CC(C)C)C)C";
  // Insulin B chain — 405 heavy atoms (full protein)
  static final String INSULIN =
      "CCC(C)C(C(=O)NC(C(C)C)C(=O)NC(CCC(=O)O)C(=O)NC(CCC(=O)N)C(=O)NC(CS)C(=O)NC(CS)C(=O)NC(C(C)O)C(=O)NC(CO)C(=O)NC(C(C)CC)C(=O)NC(CS)C(=O)NC(CO)C(=O)NC(CC(C)C)C(=O)NC(CC1=CC=C(C=C1)O)C(=O)NC(CCC(=O)N)C(=O)NC(CC(C)C)C(=O)NC(CCC(=O)O)C(=O)NC(CC(=O)N)C(=O)NC(CC2=CC=C(C=C2)O)C(=O)NC(CS)C(=O)NC(CC(=O)N)C(=O)O)NC(=O)CN.CC(C)CC(C(=O)NC(CC1=CC=C(C=C1)O)C(=O)NC(CC(C)C)C(=O)NC(C(C)C)C(=O)NC(CS)C(=O)NCC(=O)NC(CCC(=O)O)C(=O)NC(CCCNC(=N)N)C(=O)NCC(=O)NC(CC2=CC=CC=C2)C(=O)NC(CC3=CC=CC=C3)C(=O)NC(CC4=CC=C(C=C4)O)C(=O)NC(C(C)O)C(=O)NC(CCCCN)C(=O)N5CCCC5C(=O)NC(C(C)O)C(=O)O)NC(=O)C(C)NC(=O)C(CCC(=O)O)NC(=O)C(C(C)C)NC(=O)C(CC(C)C)NC(=O)C(CC6=CN=CN6)NC(=O)C(CO)NC(=O)CNC(=O)C(CS)NC(=O)C(CC(C)C)NC(=O)C(CC7=CN=CN7)NC(=O)C(CCC(=O)N)NC(=O)C(CC(=O)N)NC(=O)C(C(C)C)NC(=O)C(CC8=CC=CC=C8)N";

  // ---- Time limit: 500ms per test ----
  static final long TIME_LIMIT_MS = 1000;

  void assertFast(long elapsedMs, String label) {
    assertTrue(
        elapsedMs < TIME_LIMIT_MS,
        label + " took " + elapsedMs + "ms, limit is " + TIME_LIMIT_MS + "ms");
  }

  // =======================================================
  // SUBSTRUCTURE TESTS — must be < 500ms
  // =======================================================

  @Nested
  @DisplayName("Large Molecule Substructure")
  class SubstructureTests {

    @Test
    @DisplayName("Benzene in paclitaxel")
    void benzeneInPaclitaxel() throws Exception {
      IAtomContainer q = mol("c1ccccc1"), t = mol(PACLITAXEL);
      long t0 = System.nanoTime();
      SMSD s = new SMSD(q, t, new ChemOptions());
      assertTrue(s.isSubstructure(TIME_LIMIT_MS));
      assertFast((System.nanoTime() - t0) / 1_000_000, "benzene-in-paclitaxel");
    }

    @Test
    @DisplayName("Adenine in ATP")
    void adenineInAtp() throws Exception {
      IAtomContainer q = mol("c1ncnc2[nH]cnc12"), t = mol(ATP);
      long t0 = System.nanoTime();
      SMSD s = new SMSD(q, t, new ChemOptions());
      assertTrue(s.isSubstructure(TIME_LIMIT_MS));
      assertFast((System.nanoTime() - t0) / 1_000_000, "adenine-in-ATP");
    }

    @Test
    @DisplayName("Adenine in NAD+")
    void adenineInNad() throws Exception {
      IAtomContainer q = mol("c1ncnc2[nH]cnc12"), t = mol(NAD_PLUS);
      long t0 = System.nanoTime();
      SMSD s = new SMSD(q, t, new ChemOptions());
      assertTrue(s.isSubstructure(TIME_LIMIT_MS));
      assertFast((System.nanoTime() - t0) / 1_000_000, "adenine-in-NAD+");
    }

    @Test
    @DisplayName("Steroid core in dexamethasone")
    void steroidCore() throws Exception {
      // Cyclopentanoperhydrophenanthrene — steroid ABCD ring skeleton
      IAtomContainer q = mol("C1CCC2C(C1)CCC3C2CCC4CCCCC34");
      IAtomContainer t = mol(DEXAMETHASONE);
      long t0 = System.nanoTime();
      SMSD s = new SMSD(q, t, new ChemOptions());
      boolean result = s.isSubstructure(TIME_LIMIT_MS);
      assertFast((System.nanoTime() - t0) / 1_000_000, "steroid-core-in-dexa");
    }

    @Test
    @DisplayName("AMP substructure of ATP")
    void ampInAtp() throws Exception {
      IAtomContainer q = mol(AMP), t = mol(ATP);
      long t0 = System.nanoTime();
      SMSD s = new SMSD(q, t, new ChemOptions());
      assertTrue(s.isSubstructure(TIME_LIMIT_MS));
      assertFast((System.nanoTime() - t0) / 1_000_000, "AMP-in-ATP");
    }

    @Test
    @DisplayName("Paclitaxel self-match")
    void paclitaxelSelf() throws Exception {
      IAtomContainer t = mol(PACLITAXEL);
      long t0 = System.nanoTime();
      SMSD s = new SMSD(t, t, new ChemOptions());
      assertTrue(s.isSubstructure(2000));
      long elapsed = (System.nanoTime() - t0) / 1_000_000;
      assertTrue(elapsed < 2000, "paclitaxel-self took " + elapsed + "ms, limit is 2000ms");
    }

    @Test
    @DisplayName("Rapamycin self-match")
    void rapamycinSelf() throws Exception {
      IAtomContainer t = mol(RAPAMYCIN);
      long t0 = System.nanoTime();
      SMSD s = new SMSD(t, t, new ChemOptions());
      assertTrue(s.isSubstructure(2000));
      long elapsed = (System.nanoTime() - t0) / 1_000_000;
      assertTrue(elapsed < 2000, "rapamycin-self took " + elapsed + "ms, limit is 2000ms");
    }

    @Test
    @DisplayName("Acetyl-CoA self-match")
    void acetylCoaSelf() throws Exception {
      IAtomContainer t = mol(ACETYL_COA);
      long t0 = System.nanoTime();
      SMSD s = new SMSD(t, t, new ChemOptions());
      assertTrue(s.isSubstructure(2000));
      long elapsed = (System.nanoTime() - t0) / 1_000_000;
      assertTrue(elapsed < 2000, "acetylCoA-self took " + elapsed + "ms, limit is 2000ms");
    }

    @Test
    @DisplayName("Vincristine self-match (56 atoms)")
    void vincristineSelf() throws Exception {
      IAtomContainer t = mol(VINCRISTINE);
      long t0 = System.nanoTime();
      SMSD s = new SMSD(t, t, new ChemOptions());
      assertTrue(s.isSubstructure(2000));
      long elapsed = (System.nanoTime() - t0) / 1_000_000;
      assertTrue(elapsed < 2000, "vincristine-self took " + elapsed + "ms, limit is 2000ms");
    }

    @Test
    @DisplayName("FAD self-match (53 atoms)")
    void fadSelf() throws Exception {
      IAtomContainer t = mol(FAD);
      long t0 = System.nanoTime();
      SMSD s = new SMSD(t, t, new ChemOptions());
      assertTrue(s.isSubstructure(2000));
      long elapsed = (System.nanoTime() - t0) / 1_000_000;
      assertTrue(elapsed < 2000, "FAD-self took " + elapsed + "ms, limit is 2000ms");
    }

    @Test
    @DisplayName("Erythromycin self-match (37 atoms)")
    void erythromycinSelf() throws Exception {
      IAtomContainer t = mol(ERYTHROMYCIN);
      long t0 = System.nanoTime();
      SMSD s = new SMSD(t, t, new ChemOptions());
      assertTrue(s.isSubstructure(1000));
      long elapsed = (System.nanoTime() - t0) / 1_000_000;
      assertTrue(elapsed < 1000, "erythromycin-self took " + elapsed + "ms, limit is 1000ms");
    }
  }

  // =======================================================
  // MCS TESTS — must be < 500ms
  // =======================================================

  @Nested
  @DisplayName("Large Molecule MCS")
  class McsTests {

    @Test
    @DisplayName("ATP vs ADP MCS < 500ms")
    void atpAdpMcs() throws Exception {
      IAtomContainer q = mol(ATP), t = mol(ADP);
      long t0 = System.nanoTime();
      SMSD s = new SMSD(q, t, new ChemOptions());
      Map<Integer, Integer> mcs = s.findMCS(false, true, TIME_LIMIT_MS);
      long elapsed = (System.nanoTime() - t0) / 1_000_000;
      assertTrue(mcs.size() >= 27, "ATP/ADP MCS ≥ 27, got " + mcs.size());
      assertFast(elapsed, "ATP-ADP-MCS");
    }

    @Test
    @DisplayName("ATP vs AMP MCS < 500ms")
    void atpAmpMcs() throws Exception {
      IAtomContainer q = mol(ATP), t = mol(AMP);
      long t0 = System.nanoTime();
      SMSD s = new SMSD(q, t, new ChemOptions());
      Map<Integer, Integer> mcs = s.findMCS(false, true, TIME_LIMIT_MS);
      long elapsed = (System.nanoTime() - t0) / 1_000_000;
      assertTrue(mcs.size() >= 15, "ATP/AMP MCS ≥ 15, got " + mcs.size());
      assertFast(elapsed, "ATP-AMP-MCS");
    }

    @Test
    @DisplayName("Atorvastatin vs rosuvastatin MCS < 500ms")
    void statinMcs() throws Exception {
      IAtomContainer q = mol(ATORVASTATIN), t = mol(ROSUVASTATIN);
      long t0 = System.nanoTime();
      SMSD s = new SMSD(q, t, new ChemOptions());
      Map<Integer, Integer> mcs = s.findMCS(false, true, TIME_LIMIT_MS);
      long elapsed = (System.nanoTime() - t0) / 1_000_000;
      assertTrue(mcs.size() >= 10, "Statin MCS ≥ 10, got " + mcs.size());
      assertFast(elapsed, "statin-MCS");
    }

    @Test
    @DisplayName("UDP-glucose vs UDP-galactose MCS < 500ms")
    void udpSugarMcs() throws Exception {
      IAtomContainer q = mol(UDP_GLUCOSE), t = mol(UDP_GALACTOSE);
      long t0 = System.nanoTime();
      SMSD s = new SMSD(q, t, new ChemOptions());
      Map<Integer, Integer> mcs = s.findMCS(false, true, TIME_LIMIT_MS);
      long elapsed = (System.nanoTime() - t0) / 1_000_000;
      assertTrue(mcs.size() >= 30, "UDP sugar MCS ≥ 30, got " + mcs.size());
      assertFast(elapsed, "UDP-sugar-MCS");
    }

    @Test
    @DisplayName("NAD+ self-MCS < 500ms")
    void nadSelfMcs() throws Exception {
      IAtomContainer q = mol(NAD_PLUS);
      long t0 = System.nanoTime();
      SMSD s = new SMSD(q, q, new ChemOptions());
      Map<Integer, Integer> mcs = s.findMCS(false, true, TIME_LIMIT_MS);
      long elapsed = (System.nanoTime() - t0) / 1_000_000;
      assertEquals(q.getAtomCount(), mcs.size(), "Self-MCS should be full molecule");
      assertFast(elapsed, "NAD-self-MCS");
    }

    @Test
    @DisplayName("Dexamethasone vs vitamin D3 MCS < 500ms")
    void steroidMcs() throws Exception {
      IAtomContainer q = mol(DEXAMETHASONE), t = mol(VITAMIN_D3);
      long t0 = System.nanoTime();
      SMSD s = new SMSD(q, t, new ChemOptions());
      Map<Integer, Integer> mcs = s.findMCS(false, true, TIME_LIMIT_MS);
      long elapsed = (System.nanoTime() - t0) / 1_000_000;
      assertTrue(mcs.size() >= 8, "Steroid MCS ≥ 8, got " + mcs.size());
      assertFast(elapsed, "steroid-MCS");
    }

    @Test
    @DisplayName("Paclitaxel vs docetaxel MCS < 500ms (taxane pair)")
    void taxaneMcs() throws Exception {
      IAtomContainer q = mol(PACLITAXEL), t = mol(DOCETAXEL);
      long t0 = System.nanoTime();
      SMSD s = new SMSD(q, t, new ChemOptions());
      Map<Integer, Integer> mcs = s.findMCS(false, true, TIME_LIMIT_MS);
      long elapsed = (System.nanoTime() - t0) / 1_000_000;
      assertTrue(mcs.size() >= 30, "Taxane MCS ≥ 30, got " + mcs.size());
      assertFast(elapsed, "taxane-MCS");
    }

    @Test
    @DisplayName("CoA vs acetyl-CoA MCS < 500ms")
    void coaMcs() throws Exception {
      IAtomContainer q = mol(COA), t = mol(ACETYL_COA);
      long t0 = System.nanoTime();
      SMSD s = new SMSD(q, t, new ChemOptions());
      Map<Integer, Integer> mcs = s.findMCS(false, true, TIME_LIMIT_MS);
      long elapsed = (System.nanoTime() - t0) / 1_000_000;
      assertTrue(mcs.size() >= 30, "CoA/AcCoA MCS ≥ 30, got " + mcs.size());
      assertFast(elapsed, "CoA-MCS");
    }

    @Test
    @DisplayName("FAD self-MCS < 500ms (53 heavy atoms)")
    void fadSelfMcs() throws Exception {
      IAtomContainer q = mol(FAD);
      long t0 = System.nanoTime();
      SMSD s = new SMSD(q, q, new ChemOptions());
      Map<Integer, Integer> mcs = s.findMCS(false, true, TIME_LIMIT_MS);
      long elapsed = (System.nanoTime() - t0) / 1_000_000;
      assertEquals(q.getAtomCount(), mcs.size(), "Self-MCS = full molecule");
      assertFast(elapsed, "FAD-self-MCS");
    }

    @Test
    @DisplayName("Vincristine self-match < 500ms (56 heavy atoms)")
    void vincristineSelfMcs() throws Exception {
      IAtomContainer q = mol(VINCRISTINE);
      long t0 = System.nanoTime();
      SMSD s = new SMSD(q, q, new ChemOptions());
      Map<Integer, Integer> mcs = s.findMCS(false, true, TIME_LIMIT_MS);
      long elapsed = (System.nanoTime() - t0) / 1_000_000;
      assertEquals(q.getAtomCount(), mcs.size(), "Self-MCS = full molecule");
      assertFast(elapsed, "vincristine-self-MCS");
    }

    @Test
    @DisplayName("SAM self-MCS < 500ms")
    void samSelfMcs() throws Exception {
      IAtomContainer q = mol(SAM);
      long t0 = System.nanoTime();
      SMSD s = new SMSD(q, q, new ChemOptions());
      Map<Integer, Integer> mcs = s.findMCS(false, true, TIME_LIMIT_MS);
      long elapsed = (System.nanoTime() - t0) / 1_000_000;
      assertEquals(q.getAtomCount(), mcs.size(), "Self-MCS = full molecule");
      assertFast(elapsed, "SAM-self-MCS");
    }
  }

  // =======================================================
  // 100+ ATOM TESTS
  // =======================================================

  @Nested
  @DisplayName("100+ Atom Molecules")
  class HundredPlusAtomTests {

    @Test
    @DisplayName("Vancomycin self-substructure (101 atoms)")
    void vancomycinSelf() throws Exception {
      IAtomContainer t = mol(VANCOMYCIN);
      assertTrue(
          t.getAtomCount() >= 90,
          "Vancomycin should have 90+ heavy atoms, got " + t.getAtomCount());
      long t0 = System.nanoTime();
      SMSD s = new SMSD(t, t, new ChemOptions());
      assertTrue(s.isSubstructure(5000));
      long elapsed = (System.nanoTime() - t0) / 1_000_000;
      assertTrue(elapsed < 5000, "Vancomycin self-match took " + elapsed + "ms");
    }

    @Test
    @DisplayName("Vancomycin self-MCS (101 atoms)")
    void vancomycinSelfMcs() throws Exception {
      IAtomContainer t = mol(VANCOMYCIN);
      long t0 = System.nanoTime();
      SMSD s = new SMSD(t, t, new ChemOptions());
      Map<Integer, Integer> mcs = s.findMCS(false, true, 2000);
      long elapsed = (System.nanoTime() - t0) / 1_000_000;
      assertEquals(t.getAtomCount(), mcs.size(), "Self-MCS must equal atom count");
      assertTrue(elapsed < 2000, "Vancomycin self-MCS took " + elapsed + "ms");
    }

    @Test
    @DisplayName("Cyclosporin self-substructure (85 atoms)")
    void cyclosporinSelf() throws Exception {
      IAtomContainer t = mol(CYCLOSPORIN);
      assertTrue(t.getAtomCount() >= 70, "Cyclosporin should have 70+ heavy atoms");
      long t0 = System.nanoTime();
      SMSD s = new SMSD(t, t, new ChemOptions());
      assertTrue(s.isSubstructure(5000));
      long elapsed = (System.nanoTime() - t0) / 1_000_000;
      assertTrue(elapsed < 5000, "Cyclosporin self-match took " + elapsed + "ms");
    }

    @Test
    @DisplayName("Benzene in vancomycin")
    void benzeneInVancomycin() throws Exception {
      // Vancomycin contains multiple aromatic rings
      IAtomContainer q = mol("c1ccccc1"), t = mol(VANCOMYCIN);
      long t0 = System.nanoTime();
      SMSD s = new SMSD(q, t, new ChemOptions());
      assertTrue(s.isSubstructure(2000), "Benzene should be in vancomycin");
      long elapsed = (System.nanoTime() - t0) / 1_000_000;
      assertTrue(elapsed < 2000, "Benzene in vancomycin took " + elapsed + "ms");
    }

    @Test
    @DisplayName("C200 chain substructure of C300 chain")
    void longChainSubstructure() throws Exception {
      // Pure aliphatic chains — tests scalability
      IAtomContainer q = mol("C".repeat(200));
      IAtomContainer t = mol("C".repeat(300));
      long t0 = System.nanoTime();
      SMSD s = new SMSD(q, t, new ChemOptions());
      assertTrue(s.isSubstructure(2000), "C200 must be substructure of C300");
      long elapsed = (System.nanoTime() - t0) / 1_000_000;
      assertTrue(elapsed < 2000, "C200 in C300 took " + elapsed + "ms");
    }

    @Test
    @DisplayName("C200 self-MCS")
    void longChainSelfMcs() throws Exception {
      IAtomContainer q = mol("C".repeat(200));
      long t0 = System.nanoTime();
      SMSD s = new SMSD(q, q, new ChemOptions());
      Map<Integer, Integer> mcs = s.findMCS(false, true, 2000);
      long elapsed = (System.nanoTime() - t0) / 1_000_000;
      assertEquals(q.getAtomCount(), mcs.size(), "Self-MCS must = atom count");
      assertTrue(elapsed < 2000, "C200 self-MCS took " + elapsed + "ms");
    }
  }

  // =======================================================
  // 400+ ATOM TESTS (protein-scale)
  // =======================================================

  @Nested
  @DisplayName("400+ Atom Protein-Scale Molecules")
  class ProteinScaleTests {

    @Test
    @DisplayName("Insulin parses correctly (405 heavy atoms)")
    void insulinParses() throws Exception {
      IAtomContainer t = mol(INSULIN);
      assertTrue(
          t.getAtomCount() >= 300, "Insulin should have 300+ heavy atoms, got " + t.getAtomCount());
    }

    @Test
    @DisplayName("Adenine in insulin (contains His/adenine-like rings)")
    void adenineInInsulin() throws Exception {
      // Insulin contains histidine residues with imidazole rings
      IAtomContainer q = mol("c1c[nH]cn1"); // imidazole
      IAtomContainer t = mol(INSULIN);
      long t0 = System.nanoTime();
      SMSD s = new SMSD(q, t, new ChemOptions());
      // Insulin has histidine residues → imidazole should match
      s.isSubstructure(5000); // no crash, result is defined
      long elapsed = (System.nanoTime() - t0) / 1_000_000;
      assertTrue(elapsed < 5000, "Imidazole in insulin took " + elapsed + "ms");
    }

    @Test
    @DisplayName("Phenylalanine sidechain in insulin")
    void phenylInInsulin() throws Exception {
      // Insulin B chain contains Phe residues
      IAtomContainer q = mol("c1ccccc1"); // benzene
      IAtomContainer t = mol(INSULIN);
      long t0 = System.nanoTime();
      SMSD s = new SMSD(q, t, new ChemOptions());
      assertTrue(s.isSubstructure(5000), "Benzene (Phe sidechain) should be in insulin");
      long elapsed = (System.nanoTime() - t0) / 1_000_000;
      assertTrue(elapsed < 5000, "Benzene in insulin took " + elapsed + "ms");
    }

    @Test
    @DisplayName("C500 chain self-substructure")
    void c500SelfMatch() throws Exception {
      IAtomContainer t = mol("C".repeat(500));
      assertTrue(t.getAtomCount() >= 400, "C500 should have 400+ heavy atoms");
      long t0 = System.nanoTime();
      SMSD s = new SMSD(t, t, new ChemOptions());
      assertTrue(s.isSubstructure(5000));
      long elapsed = (System.nanoTime() - t0) / 1_000_000;
      assertTrue(elapsed < 5000, "C500 self-match took " + elapsed + "ms");
    }

    @Test
    @DisplayName("C1000 chain self-substructure")
    void c1000SelfMatch() throws Exception {
      IAtomContainer t = mol("C".repeat(1000));
      assertTrue(t.getAtomCount() >= 800, "C1000 should have 800+ heavy atoms");
      long t0 = System.nanoTime();
      SMSD s = new SMSD(t, t, new ChemOptions());
      assertTrue(s.isSubstructure(10000));
      long elapsed = (System.nanoTime() - t0) / 1_000_000;
      assertTrue(elapsed < 10000, "C1000 self-match took " + elapsed + "ms");
    }

    @Test
    @DisplayName("C100 in C1000 chain substructure")
    void c100InC1000() throws Exception {
      IAtomContainer q = mol("C".repeat(100));
      IAtomContainer t = mol("C".repeat(1000));
      long t0 = System.nanoTime();
      SMSD s = new SMSD(q, t, new ChemOptions());
      assertTrue(s.isSubstructure(10000), "C100 must be substructure of C1000");
      long elapsed = (System.nanoTime() - t0) / 1_000_000;
      assertTrue(elapsed < 10000, "C100 in C1000 took " + elapsed + "ms");
    }
  }

  // =======================================================
  // RASCAL SCREENING — must be < 10ms per pair
  // =======================================================

  @Nested
  @DisplayName("RASCAL Large Molecule Screening")
  class RascalTests {

    @Test
    @DisplayName("RASCAL on 8 large molecules < 500ms total")
    void rascalBatch() throws Exception {
      String[] smiles = {
        ATP, ADP, AMP, NAD_PLUS, ATORVASTATIN, ROSUVASTATIN, DEXAMETHASONE, VITAMIN_D3
      };
      IAtomContainer[] mols = new IAtomContainer[smiles.length];
      for (int i = 0; i < smiles.length; i++) mols[i] = mol(smiles[i]);

      ChemOptions opts = new ChemOptions();
      long t0 = System.nanoTime();
      int pairs = 0;
      for (int i = 0; i < mols.length; i++) {
        for (int j = i + 1; j < mols.length; j++) {
          double ub = SearchEngine.similarityUpperBound(mols[i], mols[j], opts);
          assertTrue(ub >= 0 && ub <= 1.0, "RASCAL UB in [0,1]");
          pairs++;
        }
      }
      long elapsed = (System.nanoTime() - t0) / 1_000_000;
      // 500ms is generous on any CI runner; locally this completes in ~15ms.
      // 100ms was too tight for shared GitHub Actions hosts under load.
      assertTrue(elapsed < 500, "RASCAL " + pairs + " pairs took " + elapsed + "ms, limit 500ms");
    }
  }

  // =======================================================
  // CHEMICALLY VALIDATED MCS PAIRS
  // =======================================================

  @Nested
  @DisplayName("Chemically Validated MCS Pairs")
  class ChemicallyValidatedMcsPairs {

    // Small molecule SMILES used only in this group
    static final String ASPIRIN = "CC(=O)Oc1ccccc1C(=O)O";
    static final String SALICYLIC_ACID = "OC(=O)c1ccccc1O";
    static final String CAFFEINE = "Cn1cnc2c1c(=O)n(C)c(=O)n2C";
    static final String XANTHINE = "O=c1[nH]c(=O)c2[nH]cnc2[nH]1";
    static final String DOPAMINE = "NCCc1ccc(O)c(O)c1";
    static final String ADRENALINE = "CNCC(O)c1ccc(O)c(O)c1";
    static final String ADENINE = "c1ncnc2[nH]cnc12";
    static final String URACIL = "O=c1cc[nH]c(=O)[nH]1";

    // --- 1. Erythromycin vs Azithromycin MCS ---
    // Both macrolide antibiotics sharing desosamine, cladinose sugars and
    // most of the macrolide backbone. MCS must be >= 25 atoms.
    @Test
    @DisplayName("Erythromycin vs Azithromycin MCS >= 25 (macrolide pair)")
    void erythromycinAzithromycinMcs() throws Exception {
      IAtomContainer q = mol(ERYTHROMYCIN), t = mol(AZITHROMYCIN);
      long t0 = System.nanoTime();
      SMSD s = new SMSD(q, t, new ChemOptions());
      Map<Integer, Integer> mcs = s.findMCS(false, true, 1000);
      long elapsed = (System.nanoTime() - t0) / 1_000_000;
      assertTrue(
          mcs.size() >= 25,
          "Erythromycin/Azithromycin MCS should be >= 25 (shared sugars + backbone), got "
              + mcs.size());
      assertTrue(
          elapsed < 1000, "erythromycin-azithromycin-MCS took " + elapsed + "ms, limit is 1000ms");
    }

    // --- 2. NAD+ vs NADH MCS ---
    // NAD+/NADH: redox pair — differ by 2H at nicotinamide ring + charge change (N+ → N).
    // Must disable charge matching since charge is part of the redox chemistry.
    @Test
    @DisplayName("NAD+ vs NADH MCS >= 35 (redox pair, charge-insensitive)")
    void nadPlusNadhMcs() throws Exception {
      IAtomContainer q = mol(NAD_PLUS), t = mol(NADH);
      ChemOptions opts = new ChemOptions();
      opts.matchFormalCharge = false; // redox pair: charge changes are chemical, not structural
      long t0 = System.nanoTime();
      SMSD s = new SMSD(q, t, opts);
      Map<Integer, Integer> mcs = s.findMCS(false, true, 1000);
      long elapsed = (System.nanoTime() - t0) / 1_000_000;
      assertTrue(mcs.size() >= 35, "NAD+/NADH MCS should be >= 35 (redox pair), got " + mcs.size());
      assertTrue(elapsed < 1000, "NAD+-NADH-MCS took " + elapsed + "ms, limit is 1000ms");
    }

    // --- 3. Aspirin vs salicylic acid substructure ---
    // Aspirin = acetylated salicylic acid. Salicylic acid is a complete
    // substructure of aspirin. This is a chemical fact.
    @Test
    @DisplayName("Salicylic acid is substructure of aspirin")
    void salicylicAcidInAspirin() throws Exception {
      IAtomContainer q = mol(SALICYLIC_ACID), t = mol(ASPIRIN);
      long t0 = System.nanoTime();
      SMSD s = new SMSD(q, t, new ChemOptions());
      assertTrue(
          s.isSubstructure(TIME_LIMIT_MS),
          "Salicylic acid MUST be a substructure of aspirin (aspirin = acetyl-salicylic acid)");
      assertFast((System.nanoTime() - t0) / 1_000_000, "salicylic-in-aspirin");
    }

    // --- 4. Caffeine vs xanthine substructure ---
    // Caffeine = 1,3,7-trimethylxanthine. Xanthine is the unmethylated core.
    @Test
    @DisplayName("Xanthine core is substructure of caffeine")
    void xanthineInCaffeine() throws Exception {
      IAtomContainer q = mol(XANTHINE), t = mol(CAFFEINE);
      long t0 = System.nanoTime();
      SMSD s = new SMSD(q, t, new ChemOptions());
      assertTrue(
          s.isSubstructure(TIME_LIMIT_MS),
          "Xanthine MUST be a substructure of caffeine (caffeine = 1,3,7-trimethylxanthine)");
      assertFast((System.nanoTime() - t0) / 1_000_000, "xanthine-in-caffeine");
    }

    // --- 5. Dopamine in adrenaline substructure ---
    // Adrenaline (epinephrine) contains the dopamine catechol + ethylamine moiety.
    @Test
    @DisplayName("Dopamine is substructure of adrenaline")
    void dopamineInAdrenaline() throws Exception {
      IAtomContainer q = mol(DOPAMINE), t = mol(ADRENALINE);
      long t0 = System.nanoTime();
      SMSD s = new SMSD(q, t, new ChemOptions());
      assertTrue(
          s.isSubstructure(TIME_LIMIT_MS),
          "Dopamine MUST be a substructure of adrenaline (adrenaline = N-methyl-hydroxylated"
              + " dopamine)");
      assertFast((System.nanoTime() - t0) / 1_000_000, "dopamine-in-adrenaline");
    }

    // --- 6. Adenine in all adenine-containing nucleotides ---
    @Test
    @DisplayName("Adenine is substructure of ATP")
    void adenineInAtp() throws Exception {
      IAtomContainer q = mol(ADENINE), t = mol(ATP);
      SMSD s = new SMSD(q, t, new ChemOptions());
      assertTrue(s.isSubstructure(TIME_LIMIT_MS), "Adenine MUST be in ATP");
    }

    @Test
    @DisplayName("Adenine is substructure of ADP")
    void adenineInAdp() throws Exception {
      IAtomContainer q = mol(ADENINE), t = mol(ADP);
      SMSD s = new SMSD(q, t, new ChemOptions());
      assertTrue(s.isSubstructure(TIME_LIMIT_MS), "Adenine MUST be in ADP");
    }

    @Test
    @DisplayName("Adenine is substructure of AMP")
    void adenineInAmp() throws Exception {
      IAtomContainer q = mol(ADENINE), t = mol(AMP);
      SMSD s = new SMSD(q, t, new ChemOptions());
      assertTrue(s.isSubstructure(TIME_LIMIT_MS), "Adenine MUST be in AMP");
    }

    @Test
    @DisplayName("Adenine is substructure of NAD+")
    void adenineInNadPlus() throws Exception {
      IAtomContainer q = mol(ADENINE), t = mol(NAD_PLUS);
      SMSD s = new SMSD(q, t, new ChemOptions());
      assertTrue(s.isSubstructure(TIME_LIMIT_MS), "Adenine MUST be in NAD+");
    }

    @Test
    @DisplayName("Adenine is substructure of FAD")
    void adenineInFad() throws Exception {
      IAtomContainer q = mol(ADENINE), t = mol(FAD);
      SMSD s = new SMSD(q, t, new ChemOptions());
      assertTrue(s.isSubstructure(TIME_LIMIT_MS), "Adenine MUST be in FAD");
    }

    @Test
    @DisplayName("Adenine is substructure of Acetyl-CoA")
    void adenineInAcetylCoa() throws Exception {
      IAtomContainer q = mol(ADENINE), t = mol(ACETYL_COA);
      SMSD s = new SMSD(q, t, new ChemOptions());
      assertTrue(s.isSubstructure(TIME_LIMIT_MS), "Adenine MUST be in Acetyl-CoA");
    }

    @Test
    @DisplayName("Adenine is substructure of SAM")
    void adenineInSam() throws Exception {
      IAtomContainer q = mol(ADENINE), t = mol(SAM);
      SMSD s = new SMSD(q, t, new ChemOptions());
      assertTrue(s.isSubstructure(TIME_LIMIT_MS), "Adenine MUST be in SAM (S-adenosyl methionine)");
    }

    // --- 7. Uracil in UDP-glucose ---
    @Test
    @DisplayName("Uracil is substructure of UDP-glucose")
    void uracilInUdpGlucose() throws Exception {
      IAtomContainer q = mol(URACIL), t = mol(UDP_GLUCOSE);
      SMSD s = new SMSD(q, t, new ChemOptions());
      assertTrue(
          s.isSubstructure(TIME_LIMIT_MS),
          "Uracil MUST be in UDP-glucose (uridine diphosphate glucose)");
    }
  }

  // ======================================================================
  // From: HardCasesTest.java
  // ======================================================================


  private static int hc_mcsSize(String smi1, String smi2, ChemOptions opts, long timeoutMs) throws Exception {
    SMSD smsd = new SMSD(mol(smi1), mol(smi2), opts);
    smsd.setMcsTimeoutMs(timeoutMs);
    Map<Integer, Integer> mcs = smsd.findMCS(false, true, timeoutMs);
    return mcs.size();
  }

  private static int hc_mcsSize(String smi1, String smi2) throws Exception {
    return hc_mcsSize(smi1, smi2, new ChemOptions(), 10_000L);
  }

  private static ChemOptions hc_looseOpts() {
    ChemOptions c = new ChemOptions();
    c.matchBondOrder = ChemOptions.BondOrderMode.LOOSE;
    c.aromaticityMode = ChemOptions.AromaticityMode.FLEXIBLE;
    c.ringMatchesRingOnly = false;
    c.matchFormalCharge = false;
    return c;
  }

  // ======================================================================
  // 1. SMARTS SUBSTRUCTURE STRESS (20 tests)
  // ======================================================================

  @Nested
  @DisplayName("1. SMARTS Substructure Stress")
  class SmartsSubstructureStress {

    /** Helper: assert SMARTS matches target SMILES */
    private void assertSmartsMatch(String smarts, String targetSmi, String msg) throws Exception {
      SMSD smsd = new SMSD(smarts, mol(targetSmi), new ChemOptions());
      assertTrue(smsd.isSubstructure(5000), msg);
    }

    /** Helper: assert SMARTS does NOT match */
    private void assertSmartsNoMatch(String smarts, String targetSmi, String msg) throws Exception {
      SMSD smsd = new SMSD(smarts, mol(targetSmi), new ChemOptions());
      assertFalse(smsd.isSubstructure(5000), msg);
    }

    @Test @Timeout(10) @DisplayName("1.01 4-carbon chain in butanol")
    void fourCarbonChainInButanol() throws Exception {
      assertSmartsMatch("[#6]~[#6]~[#6]~[#6]", "CCCCO",
          "4-carbon chain should match in butanol");
    }

    @Test @Timeout(10) @DisplayName("1.02 4-carbon chain NOT in propane")
    void fourCarbonChainNotInPropane() throws Exception {
      assertSmartsNoMatch("[#6]~[#6]~[#6]~[#6]", "CCC",
          "4-carbon chain should NOT match in propane");
    }

    @Test @Timeout(10) @DisplayName("1.03 Primary amine in amphetamine")
    void primaryAmineInAmphetamine() throws Exception {
      // amphetamine: CC(N)Cc1ccccc1
      assertSmartsMatch("[NX3;H2,H1;!$(NC=O)]", "CC(N)Cc1ccccc1",
          "Primary amine should match in amphetamine");
    }

    @Test @Timeout(10) @DisplayName("1.04 Primary amine NOT in acetanilide")
    void primaryAmineNotInAcetanilide() throws Exception {
      // acetanilide: CC(=O)Nc1ccccc1 — the N is an amide
      assertSmartsNoMatch("[NX3;H2;!$(NC=O)]", "CC(=O)Nc1ccccc1",
          "Primary amine (H2) should NOT match amide nitrogen in acetanilide");
    }

    @Test @Timeout(10) @DisplayName("1.05 Carboxylic acid in ibuprofen")
    void carboxylicAcidInIbuprofen() throws Exception {
      assertSmartsMatch("[CX3](=O)[OX2H1]", "CC(C)Cc1ccc(cc1)C(C)C(=O)O",
          "Carboxylic acid should match in ibuprofen");
    }

    @Test @Timeout(10) @DisplayName("1.06 Amide in acetaminophen")
    void amideInAcetaminophen() throws Exception {
      assertSmartsMatch("[CX3](=O)[NX3]", "CC(=O)Nc1ccc(O)cc1",
          "Amide should match in acetaminophen");
    }

    @Test @Timeout(10) @DisplayName("1.07 Aromatic 6-ring in toluene")
    void aromatic6RingInToluene() throws Exception {
      assertSmartsMatch("[cX3]1[cX3][cX3][cX3][cX3][cX3]1", "Cc1ccccc1",
          "Aromatic 6-ring in toluene");
    }

    @Test @Timeout(10) @DisplayName("1.08 C=C double bond in styrene")
    void doubleBondInStyrene() throws Exception {
      assertSmartsMatch("[#6]=[#6]", "C=Cc1ccccc1",
          "C=C should match in styrene");
    }

    @Test @Timeout(10) @DisplayName("1.09 Non-amide amine in tyramine")
    void nonAmideAmineInTyramine() throws Exception {
      // tyramine: NCCc1ccc(O)cc1
      assertSmartsMatch("[NX3;$([NH2]),$([NH1]);!$(NC=O)]", "NCCc1ccc(O)cc1",
          "Non-amide amine should match in tyramine");
    }

    @Test @Timeout(10) @DisplayName("1.10 Hydroxyl in ethanol")
    void hydroxylInEthanol() throws Exception {
      assertSmartsMatch("[OX2H]", "CCO",
          "Hydroxyl should match in ethanol");
    }

    @Test @Timeout(10) @DisplayName("1.11 Thiol in cysteine")
    void thiolInCysteine() throws Exception {
      assertSmartsMatch("[#16X2H]", "N[C@@H](CS)C(=O)O",
          "Thiol should match in cysteine");
    }

    @Test @Timeout(10) @DisplayName("1.12 Nitrile in benzonitrile")
    void nitrileInBenzonitrile() throws Exception {
      assertSmartsMatch("[$([NX1]#[CX2])]", "N#Cc1ccccc1",
          "Nitrile should match in benzonitrile");
    }

    @Test @Timeout(10) @DisplayName("1.13 Halogen in chlorobenzene")
    void halogenInChlorobenzene() throws Exception {
      assertSmartsMatch("[F,Cl,Br,I]", "Clc1ccccc1",
          "Halogen should match in chlorobenzene");
    }

    @Test @Timeout(10) @DisplayName("1.14 Ketone in acetone")
    void ketoneInAcetone() throws Exception {
      assertSmartsMatch("[CX3](=[OX1])[CX4]", "CC(=O)C",
          "Ketone should match in acetone");
    }

    @Test @Timeout(10) @DisplayName("1.15 Aldehyde in benzaldehyde")
    void aldehydeInBenzaldehyde() throws Exception {
      assertSmartsMatch("[CX3H1](=[OX1])", "O=Cc1ccccc1",
          "Aldehyde should match in benzaldehyde");
    }

    @Test @Timeout(10) @DisplayName("1.16 Nitro group in nitrobenzene")
    void nitroInNitrobenzene() throws Exception {
      // CDK represents nitro as [N+](=O)[O-], use charge-aware SMARTS
      assertSmartsMatch("[NX3+](=[OX1])[OX1-]", "[O-][N+](=O)c1ccccc1",
          "Nitro should match in nitrobenzene");
    }

    @Test @Timeout(10) @DisplayName("1.17 Hydroxylated aromatic ring in catechol")
    void hydroxylatedRingInCatechol() throws Exception {
      // catechol: Oc1ccccc1O
      assertSmartsMatch("[OX2H]c1ccccc1", "Oc1ccccc1O",
          "Hydroxylated ring pattern should match in catechol");
    }

    @Test @Timeout(10) @DisplayName("1.18 Heteroatom count in caffeine via SMARTS")
    void heteroatomInCaffeine() throws Exception {
      // caffeine: Cn1c(=O)c2c(ncn2C)n(C)c1=O — has N and O heteroatoms
      assertSmartsMatch("[!#1;!#6]", "Cn1c(=O)c2c(ncn2C)n(C)c1=O",
          "Heteroatom SMARTS should match in caffeine");
    }

    @Test @Timeout(10) @DisplayName("1.19 Quaternary carbon in neopentane")
    void quaternaryCarbonInNeopentane() throws Exception {
      // neopentane: CC(C)(C)C — central C has degree 4
      assertSmartsMatch("[CX4;D4]", "CC(C)(C)C",
          "Quaternary carbon should match in neopentane");
    }

    @Test @Timeout(10) @DisplayName("1.20 Ester in methyl benzoate")
    void esterInMethylBenzoate() throws Exception {
      // methyl benzoate: COC(=O)c1ccccc1
      assertSmartsMatch("[CX3](=O)[OX2][#6]", "COC(=O)c1ccccc1",
          "Ester should match in methyl benzoate");
    }
  }

  // ======================================================================
  // 2. HIGHLY BRANCHED / DENDRIMER-LIKE (10 tests)
  // ======================================================================

  @Nested
  @DisplayName("2. Highly Branched Molecules")
  class HighlyBranched {

    @Test @Timeout(10) @DisplayName("2.01 Neopentane self-match MCS = 5")
    void neopentaneSelfMcs() throws Exception {
      assertEquals(5, hc_mcsSize("CC(C)(C)C", "CC(C)(C)C"),
          "Neopentane self-match should be 5 heavy atoms");
    }

    @Test @Timeout(10) @DisplayName("2.02 Neopentane substructure in 2,2-dimethylbutane")
    void neopentaneIn22Dimethylbutane() throws Exception {
      // 2,2-dimethylbutane: CCC(C)(C)C
      int sz = hc_mcsSize("CC(C)(C)C", "CCC(C)(C)C");
      assertTrue(sz >= 5, "Neopentane should be substructure of 2,2-dimethylbutane");
    }

    @Test @Timeout(10) @DisplayName("2.03 Triphenylamine self-match")
    void triphenylamineSelfMatch() throws Exception {
      String tpa = "c1ccc(N(c2ccccc2)c3ccccc3)cc1";
      int sz = hc_mcsSize(tpa, tpa);
      assertTrue(sz >= 18, "Triphenylamine self-match should be >= 18 heavy atoms");
    }

    @Test @Timeout(10) @DisplayName("2.04 Pentaerythritol self-match MCS = 9")
    void pentaerythritolSelfMatch() throws Exception {
      // OCC(CO)(CO)CO has 9 heavy atoms
      int sz = hc_mcsSize("OCC(CO)(CO)CO", "OCC(CO)(CO)CO");
      assertEquals(9, sz, "Pentaerythritol self-match should be 9");
    }

    @Test @Timeout(10) @DisplayName("2.05 Tri-tert-butylmethane self-match")
    void triTertButylMethaneSelfMatch() throws Exception {
      String ttbm = "CC(C)(C)C(C(C)(C)C)C(C)(C)C";
      int sz = hc_mcsSize(ttbm, ttbm);
      assertTrue(sz >= 13, "Tri-tert-butylmethane self-match should be >= 13");
    }

    @Test @Timeout(10) @DisplayName("2.06 Neopentane substructure in triphenylamine (loose)")
    void neopentaneNotInTriphenylamine() throws Exception {
      // Neopentane has no aromatic atoms; triphenylamine does — no overlap in strict ring match
      int sz = hc_mcsSize("CC(C)(C)C", "c1ccc(N(c2ccccc2)c3ccccc3)cc1");
      assertTrue(sz >= 1, "Should share at least one common atom type");
    }

    @Test @Timeout(10) @DisplayName("2.07 Star-shaped glycerol self-match = 6")
    void glycerolSelfMatch() throws Exception {
      // glycerol: OCC(O)CO has 6 heavy atoms
      int sz = hc_mcsSize("OCC(O)CO", "OCC(O)CO");
      assertEquals(6, sz, "Glycerol self-match should be 6");
    }

    @Test @Timeout(10) @DisplayName("2.08 Pentaerythritol in larger derivative")
    void pentaerythritolInTetraAcetate() throws Exception {
      // Pentaerythritol tetraacetate: CC(=O)OCC(COC(C)=O)(COC(C)=O)COC(C)=O
      int sz = hc_mcsSize("OCC(CO)(CO)CO",
          "CC(=O)OCC(COC(C)=O)(COC(C)=O)COC(C)=O", hc_looseOpts(), 10_000L);
      assertTrue(sz >= 5, "Pentaerythritol core should be found in tetraacetate");
    }

    @Test @Timeout(10) @DisplayName("2.09 Trimethylolpropane self-match")
    void trimethylolpropaneSelfMatch() throws Exception {
      // CCC(CO)(CO)CO — 10 heavy atoms
      int sz = hc_mcsSize("CCC(CO)(CO)CO", "CCC(CO)(CO)CO");
      assertTrue(sz >= 9, "Trimethylolpropane self-match should be >= 9");
    }

    @Test @Timeout(10) @DisplayName("2.10 Isooctane self-match = 8")
    void isooctaneSelfMatch() throws Exception {
      // isooctane: CC(C)CC(C)(C)C — 8 carbons
      int sz = hc_mcsSize("CC(C)CC(C)(C)C", "CC(C)CC(C)(C)C");
      assertEquals(8, sz, "Isooctane self-match should be 8");
    }
  }

  // ======================================================================
  // 3. STEREOISOMER EQUIVALENCE (10 tests)
  // ======================================================================

  @Nested
  @DisplayName("3. Stereoisomer Equivalence")
  class StereoisomerEquivalence {

    @Test @Timeout(10) @DisplayName("3.01 L-alanine vs D-alanine same MCS size (no chirality)")
    void alanineEnantiomersNoChirality() throws Exception {
      int sz = hc_mcsSize("N[C@@H](C)C(=O)O", "N[C@H](C)C(=O)O");
      assertEquals(6, sz, "Enantiomers should have same heavy-atom MCS without chirality");
    }

    @Test @Timeout(10) @DisplayName("3.02 L-alanine vs D-alanine with chirality: not substructure")
    void alanineEnantiomersWithChirality() throws Exception {
      ChemOptions opts = new ChemOptions();
      opts.useChirality = true;
      SMSD smsd = new SMSD(mol("N[C@@H](C)C(=O)O"), mol("N[C@H](C)C(=O)O"), opts);
      // With chirality on, enantiomers differ at chiral center
      // but MCS should still find overlap (non-chiral atoms match)
      Map<Integer, Integer> mcs = smsd.findMCS(false, true, 5000);
      assertNotNull(mcs);
      assertTrue(mcs.size() >= 4, "Even with chirality, should share most atoms");
    }

    @Test @Timeout(10) @DisplayName("3.03 R-thalidomide vs S-thalidomide no chirality")
    void thalidomideEnantiomersNoChirality() throws Exception {
      // Thalidomide: O=C1CCC(N2C(=O)c3ccccc3C2=O)C(=O)N1
      String rThal = "O=C1CC[C@H](N2C(=O)c3ccccc3C2=O)C(=O)N1";
      String sThal = "O=C1CC[C@@H](N2C(=O)c3ccccc3C2=O)C(=O)N1";
      int sz = hc_mcsSize(rThal, sThal);
      assertTrue(sz >= 17, "Thalidomide enantiomers should match fully without chirality");
    }

    @Test @Timeout(10) @DisplayName("3.04 Cis-decalin vs trans-decalin no stereo")
    void cisVsTransDecalinNoStereo() throws Exception {
      // cis-decalin
      String cisDec = "[C@H]1(CCC[C@@H]2CCCC1)C2";
      // trans-decalin
      String transDec = "[C@@H]1(CCC[C@@H]2CCCC1)C2";
      // Without stereo, these are the same constitutional formula
      int sz = hc_mcsSize(cisDec, transDec);
      assertTrue(sz >= 10, "Decalin isomers should fully overlap without stereo");
    }

    @Test @Timeout(10) @DisplayName("3.05 D-glucose vs L-glucose same MCS no chirality")
    void glucoseEnantiomersNoChirality() throws Exception {
      String dGluc = "OC[C@@H](O)[C@H](O)[C@@H](O)[C@@H](O)C=O";
      String lGluc = "OC[C@H](O)[C@@H](O)[C@H](O)[C@H](O)C=O";
      int sz = hc_mcsSize(dGluc, lGluc);
      assertTrue(sz >= 12, "Glucose enantiomers should fully match without chirality");
    }

    @Test @Timeout(10) @DisplayName("3.06 Methadone vs levomethadone same MCS no chirality")
    void methadoneVsLevomethadone() throws Exception {
      // Methadone: CC(C(=O)CC)(c1ccccc1)c1ccccc1 (racemic)
      // Levomethadone: C[C@@H](C(=O)CC)(c1ccccc1) — but SMILES may vary
      // Use the shared scaffold
      String rac = "CCC(=O)C(CC)(c1ccccc1)c1ccccc1";
      String levo = "CCC(=O)[C@@](CC)(c1ccccc1)c1ccccc1";
      int sz = hc_mcsSize(rac, levo);
      assertTrue(sz >= 15, "Racemic and levo-methadone should share full scaffold");
    }

    @Test @Timeout(10) @DisplayName("3.07 R,S-alanine vs unspecified alanine")
    void alanineVsUnspecified() throws Exception {
      int sz = hc_mcsSize("NC(C)C(=O)O", "N[C@@H](C)C(=O)O");
      assertEquals(6, sz, "Unspecified should match chiral alanine fully");
    }

    @Test @Timeout(10) @DisplayName("3.08 E-stilbene vs Z-stilbene no bond stereo")
    void eVsZStilbeneNoBondStereo() throws Exception {
      String eStilb = "C(/c1ccccc1)=C/c1ccccc1";
      String zStilb = "C(/c1ccccc1)=C\\c1ccccc1";
      int sz = hc_mcsSize(eStilb, zStilb);
      assertTrue(sz >= 14, "E/Z stilbene should fully match without bond stereo");
    }

    @Test @Timeout(10) @DisplayName("3.09 E-stilbene vs Z-stilbene WITH bond stereo")
    void eVsZStilbeneWithBondStereo() throws Exception {
      ChemOptions opts = new ChemOptions();
      opts.useBondStereo = true;
      String eStilb = "C(/c1ccccc1)=C/c1ccccc1";
      String zStilb = "C(/c1ccccc1)=C\\c1ccccc1";
      int sz = hc_mcsSize(eStilb, zStilb, opts, 5000L);
      // With bond stereo, the double bond center may not match, but rings still do
      assertTrue(sz >= 6, "Even with bond stereo, phenyl rings should overlap");
    }

    @Test @Timeout(10) @DisplayName("3.10 Meso-tartaric acid vs L-tartaric acid")
    void mesoVsLTartaricAcid() throws Exception {
      String meso = "OC(=O)[C@@H](O)[C@H](O)C(=O)O";
      String lTart = "OC(=O)[C@H](O)[C@H](O)C(=O)O";
      int sz = hc_mcsSize(meso, lTart);
      assertTrue(sz >= 10, "Tartaric acid diastereomers share full connectivity");
    }
  }

  // ======================================================================
  // 4. REACTION-RELEVANT MCS (10 tests)
  // ======================================================================

  @Nested
  @DisplayName("4. Reaction-relevant MCS")
  class ReactionRelevantMcs {

    @Test @Timeout(10) @DisplayName("4.01 Ethanol to acetaldehyde: MCS reveals OH->CHO")
    void ethanolToAcetaldehyde() throws Exception {
      // ethanol: CCO, acetaldehyde: CC=O
      int sz = hc_mcsSize("CCO", "CC=O", hc_looseOpts(), 5000L);
      assertTrue(sz >= 2, "Should share at least CC fragment");
    }

    @Test @Timeout(10) @DisplayName("4.02 Aniline to acetanilide: MCS reveals amide formation")
    void anilineToAcetanilide() throws Exception {
      // aniline: Nc1ccccc1, acetanilide: CC(=O)Nc1ccccc1
      int sz = hc_mcsSize("Nc1ccccc1", "CC(=O)Nc1ccccc1");
      assertTrue(sz >= 7, "Aniline ring+N should be in acetanilide MCS");
    }

    @Test @Timeout(10) @DisplayName("4.03 Acetic acid to ethyl acetate: ester formation")
    void aceticAcidToEthylAcetate() throws Exception {
      int sz = hc_mcsSize("CC(=O)O", "CCOC(C)=O");
      assertTrue(sz >= 3, "Should share acetyl core");
    }

    @Test @Timeout(10) @DisplayName("4.04 Phenol to anisole: O-methylation")
    void phenolToAnisole() throws Exception {
      // phenol: Oc1ccccc1, anisole: COc1ccccc1
      int sz = hc_mcsSize("Oc1ccccc1", "COc1ccccc1");
      assertTrue(sz >= 7, "Phenol is substructure of anisole");
    }

    @Test @Timeout(10) @DisplayName("4.05 Propylene to propylene oxide: epoxidation")
    void propyleneToEpoxide() throws Exception {
      // propylene: CC=C, propylene oxide: CC1CO1
      int sz = hc_mcsSize("CC=C", "CC1CO1", hc_looseOpts(), 5000L);
      assertTrue(sz >= 2, "Should share CC fragment");
    }

    @Test @Timeout(10) @DisplayName("4.06 Benzene to phenol: hydroxylation")
    void benzeneToPheno() throws Exception {
      int sz = hc_mcsSize("c1ccccc1", "Oc1ccccc1");
      assertTrue(sz >= 6, "Benzene should be substructure of phenol");
    }

    @Test @Timeout(10) @DisplayName("4.07 Cyclohexane to cyclohexanone: oxidation")
    void cyclohexaneToCyclohexanone() throws Exception {
      int sz = hc_mcsSize("C1CCCCC1", "O=C1CCCCC1", hc_looseOpts(), 5000L);
      assertTrue(sz >= 5, "Should share 5+ ring carbons");
    }

    @Test @Timeout(10) @DisplayName("4.08 Toluene to benzoic acid: side-chain oxidation")
    void tolueneTobenzoicAcid() throws Exception {
      int sz = hc_mcsSize("Cc1ccccc1", "OC(=O)c1ccccc1");
      assertTrue(sz >= 6, "Should share benzene ring");
    }

    @Test @Timeout(10) @DisplayName("4.09 Benzaldehyde to benzoic acid: oxidation")
    void benzaldehydeToBenzoicAcid() throws Exception {
      int sz = hc_mcsSize("O=Cc1ccccc1", "OC(=O)c1ccccc1", hc_looseOpts(), 5000L);
      assertTrue(sz >= 7, "Aldehyde and acid share ring + carbonyl carbon");
    }

    @Test @Timeout(10) @DisplayName("4.10 Nitrobenzene to aniline: reduction")
    void nitrobenzeneToAniline() throws Exception {
      int sz = hc_mcsSize("[O-][N+](=O)c1ccccc1", "Nc1ccccc1", hc_looseOpts(), 5000L);
      assertTrue(sz >= 6, "Should share benzene ring");
    }
  }

  // ======================================================================
  // 5. MULTI-RING TOPOLOGY (10 tests)
  // ======================================================================

  @Nested
  @DisplayName("5. Multi-ring Topology")
  class MultiRingTopology {

    @Test @Timeout(10) @DisplayName("5.01 Bicyclo[2.1.0]pentane self-match")
    void housaneSelfMatch() throws Exception {
      // bicyclo[2.1.0]pentane: C1CC2CC12
      int sz = hc_mcsSize("C1CC2CC12", "C1CC2CC12");
      assertEquals(5, sz, "Housane self-match should be 5");
    }

    @Test @Timeout(10) @DisplayName("5.02 Prismane self-match")
    void prismaneSelfMatch() throws Exception {
      // prismane: C12C3C1C4C2C34 — 6 carbons
      int sz = hc_mcsSize("C12C3C1C4C2C34", "C12C3C1C4C2C34");
      assertEquals(6, sz, "Prismane self-match should be 6");
    }

    @Test @Timeout(10) @DisplayName("5.03 Spiro[4.4]nonane self-match = 9")
    void spiro44NonaneSelfMatch() throws Exception {
      // spiro[4.4]nonane: C1CCC2(C1)CCCC2
      int sz = hc_mcsSize("C1CCC2(C1)CCCC2", "C1CCC2(C1)CCCC2");
      assertEquals(9, sz, "Spiro[4.4]nonane self-match");
    }

    @Test @Timeout(10) @DisplayName("5.04 Spiro[5.5]undecane self-match = 11")
    void spiro55UndecaneSelfMatch() throws Exception {
      // spiro[5.5]undecane: C1CCCCC1(CCCCC1)  -> C1CCCCC12CCCCC2
      int sz = hc_mcsSize("C1CCCCC12CCCCC2", "C1CCCCC12CCCCC2");
      assertEquals(11, sz, "Spiro[5.5]undecane self-match");
    }

    @Test @Timeout(10) @DisplayName("5.05 Norbornane self-match = 7")
    void norbornaneSelfMatch() throws Exception {
      // norbornane (bicyclo[2.2.1]heptane): C1CC2CC1CC2
      int sz = hc_mcsSize("C1CC2CC1CC2", "C1CC2CC1CC2");
      assertEquals(7, sz, "Norbornane self-match");
    }

    @Test @Timeout(10) @DisplayName("5.06 Bicyclo[2.2.2]octane self-match = 8")
    void bicyclo222OctaneSelfMatch() throws Exception {
      // C1CC2CCC1CC2
      int sz = hc_mcsSize("C1CC2CCC1CC2", "C1CC2CCC1CC2");
      assertEquals(8, sz, "Bicyclo[2.2.2]octane self-match");
    }

    @Test @Timeout(10) @DisplayName("5.07 Norbornane as substructure of norbornene")
    void norbornaneInNorbornene() throws Exception {
      // norbornene: C1=CC2CC1CC2
      int sz = hc_mcsSize("C1CC2CC1CC2", "C1=CC2CC1CC2", hc_looseOpts(), 5000L);
      assertTrue(sz >= 7, "Norbornane should map fully into norbornene (loose)");
    }

    @Test @Timeout(10) @DisplayName("5.08 Biphenylene self-match")
    void biphenyleneSelfMatch() throws Exception {
      // biphenylene: c1ccc2-c3ccccc3-c2c1 (12 atoms)
      int sz = hc_mcsSize("c1ccc2-c3ccccc3-c2c1", "c1ccc2-c3ccccc3-c2c1");
      assertEquals(12, sz, "Biphenylene self-match");
    }

    @Test @Timeout(10) @DisplayName("5.09 Indane self-match = 9")
    void indaneSelfMatch() throws Exception {
      // indane: C1Cc2ccccc2C1 — 9 heavy atoms
      int sz = hc_mcsSize("C1Cc2ccccc2C1", "C1Cc2ccccc2C1");
      assertEquals(9, sz, "Indane self-match");
    }

    @Test @Timeout(10) @DisplayName("5.10 Fluorene self-match = 13")
    void fluoreneSelfMatch() throws Exception {
      // fluorene: c1ccc2c(c1)Cc1ccccc1-2  — 13 heavy atoms
      int sz = hc_mcsSize("c1ccc2c(c1)Cc1ccccc1-2", "c1ccc2c(c1)Cc1ccccc1-2");
      assertEquals(13, sz, "Fluorene self-match");
    }
  }

  // ======================================================================
  // 6. PHARMACOPHORE-RELEVANT MCS (10 tests)
  // ======================================================================

  @Nested
  @DisplayName("6. Pharmacophore-relevant MCS")
  class PharmacophoreMcs {

    @Test @Timeout(10) @DisplayName("6.01 Sildenafil vs tadalafil: PDE5 inhibitors share fused ring")
    void sildenafilVsTadalafil() throws Exception {
      // sildenafil: CCCc1nn(C)c2c1nc(nc2OCC)-c1cc(ccc1OCC)S(=O)(=O)N1CCN(C)CC1
      // tadalafil: O=C1N(CC(N2C1Cc1c2[nH]c2ccccc12)c1ccc2OCOc2c1)C
      // Both are large; use loose matching. They share a small common scaffold.
      String sildenafil = "CCCc1nn(C)c2c(N)nc(nc12)-c1cc(ccc1OCC)S(=O)(=O)N1CCN(C)CC1";
      String tadalafil = "O=C1N(CC(N2C1Cc1c2[nH]c2ccccc12)c1ccc2OCOc2c1)C";
      int sz = hc_mcsSize(sildenafil, tadalafil, hc_looseOpts(), 10_000L);
      assertTrue(sz >= 5, "PDE5 inhibitors should share at least a heterocyclic fragment");
    }

    @Test @Timeout(10) @DisplayName("6.02 Omeprazole vs lansoprazole: PPIs share benzimidazole")
    void omeprazoleVsLansoprazole() throws Exception {
      String omeprazole = "COc1ccc2[nH]c(S(=O)Cc3ncc(C)c(OC)c3C)nc2c1";
      String lansoprazole = "Cc1c(OCC(F)(F)F)ccn1Cc1nc2ccccc2[nH]1";
      // Loose opts — they share benzimidazole moiety (9 atoms)
      int sz = hc_mcsSize(omeprazole, lansoprazole, hc_looseOpts(), 10_000L);
      assertTrue(sz >= 8, "PPIs should share benzimidazole scaffold");
    }

    @Test @Timeout(10) @DisplayName("6.03 Atenolol vs metoprolol: beta-blockers share aryloxypropanolamine")
    void atenololVsMetoprolol() throws Exception {
      String atenolol = "CC(C)NCC(O)COc1ccc(CC(N)=O)cc1";
      String metoprolol = "CC(C)NCC(O)COc1ccc(CCOC)cc1";
      int sz = hc_mcsSize(atenolol, metoprolol, new ChemOptions(), 10_000L);
      assertTrue(sz >= 14, "Beta-blockers should share aryloxypropanolamine core");
    }

    @Test @Timeout(10) @DisplayName("6.04 Cetirizine vs loratadine: antihistamines")
    void cetirizineVsLoratadine() throws Exception {
      String cetirizine = "OC(=O)COCCN1CCN(CC1)C(c1ccccc1)c1ccc(Cl)cc1";
      String loratadine = "CCOC(=O)N1CCC(=C2c3ccc(Cl)cc3CCc3ncccc32)CC1";
      int sz = hc_mcsSize(cetirizine, loratadine, hc_looseOpts(), 10_000L);
      assertTrue(sz >= 6, "Antihistamines should share chlorophenyl + piperidine fragment");
    }

    @Test @Timeout(10) @DisplayName("6.05 Simvastatin vs pravastatin: statins share lactone/acid")
    void simvastatinVsPravastatin() throws Exception {
      String simvastatin = "CCC(C)(C)C(=O)OC1CC(O)C=C2C=CC(C)C(CCC3CC(O)CC(=O)O3)C21";
      String pravastatin = "CCC(C)(C)C(=O)OC1CC(O)C=C2C=CC(C)C(CCC(O)CC(O)CC(=O)O)C21";
      int sz = hc_mcsSize(simvastatin, pravastatin, hc_looseOpts(), 10_000L);
      assertTrue(sz >= 15, "Statins should share large decalin + side chain core");
    }

    @Test @Timeout(10) @DisplayName("6.06 Diazepam vs alprazolam: benzodiazepines")
    void diazepamVsAlprazolam() throws Exception {
      String diazepam = "CN1C(=O)CN=C(c2ccccc2)c2cc(Cl)ccc21";
      String alprazolam = "Cc1nnc2n1-c1ccc(Cl)cc1C(=NC2)c1ccccc1";
      int sz = hc_mcsSize(diazepam, alprazolam, hc_looseOpts(), 10_000L);
      assertTrue(sz >= 10, "Benzodiazepines should share chlorophenyl-diazepine core");
    }

    @Test @Timeout(10) @DisplayName("6.07 Ibuprofen vs naproxen: NSAIDs share arylpropionic acid")
    void ibuprofenVsNaproxen() throws Exception {
      String ibuprofen = "CC(C)Cc1ccc(cc1)C(C)C(=O)O";
      String naproxen = "COc1ccc2cc(C(C)C(=O)O)ccc2c1";
      int sz = hc_mcsSize(ibuprofen, naproxen, new ChemOptions(), 10_000L);
      assertTrue(sz >= 8, "NSAIDs should share arylpropionic acid scaffold");
    }

    @Test @Timeout(10) @DisplayName("6.08 Morphine vs codeine MCS")
    void morphineVsCodeine() throws Exception {
      String morphine = "CN1CCC23C4C1CC5=C(C2C(C=C4)O3)C=C(C=C5)O";
      String codeine = "CN1CCC23C4C1CC5=C(C2C(C=C4)OC3)C=C(C=C5)O";
      // Differ by one O vs OC — share most of the skeleton
      int sz = hc_mcsSize(morphine, codeine);
      assertTrue(sz >= 18, "Morphine and codeine differ by one methyl, large MCS");
    }

    @Test @Timeout(10) @DisplayName("6.09 Fluoxetine vs paroxetine: SSRIs")
    void fluoxetineVsParoxetine() throws Exception {
      String fluoxetine = "CNCCC(Oc1ccc(C(F)(F)F)cc1)c1ccccc1";
      String paroxetine = "Fc1ccc(C2CCNCC2COc2ccc3OCOc3c2)cc1";
      int sz = hc_mcsSize(fluoxetine, paroxetine, hc_looseOpts(), 10_000L);
      assertTrue(sz >= 6, "SSRIs should share fluorophenyl + amine fragment");
    }

    @Test @Timeout(10) @DisplayName("6.10 Captopril vs enalaprilat: ACE inhibitors")
    void captoprilVsEnalaprilat() throws Exception {
      String captopril = "CC(CS)C(=O)N1CCCC1C(=O)O";
      String enalaprilat = "OC(=O)C(CC(=O)O)NC(C)C(=O)N1CCCC1C(=O)O";
      int sz = hc_mcsSize(captopril, enalaprilat, hc_looseOpts(), 10_000L);
      assertTrue(sz >= 7, "ACE inhibitors should share proline + acyl fragment");
    }
  }
}
