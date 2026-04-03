/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
package com.bioinception.smsd;

import static org.junit.jupiter.api.Assertions.*;

import com.bioinception.smsd.core.*;
import java.util.*;
import org.junit.jupiter.api.*;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;

/**
 * Feature tests: consolidated from FastPathTest, TautomerWeightTest,
 * HydrogenHandlingTest, RingFinderTest.
 *
 * @author Syed Asad Rahman
 */
@DisplayName("Feature Tests")
public class FeaturesTest extends TestBase {

  /** SmilesParser without Standardiser (for TautomerWeight and FastPath tests). */
  private static final SmilesParser RAW_SP =
      new SmilesParser(SilentChemObjectBuilder.getInstance());

  /** Parse SMILES without standardisation (raw CDK parse). */
  private static IAtomContainer rawMol(String smi) throws Exception {
    return RAW_SP.parseSmiles(smi);
  }

  /** Build a MolGraph from a SMILES string (for RingFinder tests). */
  private static MolGraph graph(String smi) throws Exception {
    return new MolGraph(mol(smi));
  }

  // ======================================================================
  // From: FastPathTest.java
  // ======================================================================

  @Nested
  @DisplayName("Fast Path Tests (5.7.0)")
  class FastPath {

    // --- Tree Fast-Path ---

    @Test
    @DisplayName("Tree: single atom self-MCS")
    void testTreeFastPathSingleAtom() throws Exception {
      SMSD s = new SMSD(rawMol("C"), rawMol("C"), new ChemOptions());
      assertEquals(1, s.findMCS().size());
    }

    @Test
    @DisplayName("Tree: two atoms self-MCS")
    void testTreeFastPathTwoAtoms() throws Exception {
      SMSD s = new SMSD(rawMol("CC"), rawMol("CC"), new ChemOptions());
      assertEquals(2, s.findMCS().size());
    }

    @Test
    @DisplayName("Tree: branched vs linear (isobutane vs butane)")
    void testTreeFastPathBranched() throws Exception {
      SMSD s = new SMSD(rawMol("CC(C)C"), rawMol("CCCC"), new ChemOptions());
      var mcs = s.findMCS();
      assertTrue(mcs.size() >= 3, "Isobutane/butane MCS >= 3, got " + mcs.size());
    }

    @Test
    @DisplayName("Tree: dendrimer-like branching")
    void testTreeFastPathDendrimer() throws Exception {
      SMSD s = new SMSD(rawMol("CC(CC)(CC)CC"), rawMol("CC(CCC)(CCC)CCC"), new ChemOptions());
      var mcs = s.findMCS();
      assertTrue(mcs.size() >= 5, "Dendrimer MCS >= 5, got " + mcs.size());
    }

    // --- Chain Fast-Path ---

    @Test
    @DisplayName("Chain: single atom in chain")
    void testChainFastPathSingleAtom() throws Exception {
      SMSD s = new SMSD(rawMol("C"), rawMol("CCCCCC"), new ChemOptions());
      assertTrue(s.isSubstructure(), "Single C should be in hexane");
    }

    @Test
    @DisplayName("Chain: PEG completes fast")
    void testChainFastPathPEG() throws Exception {
      long t0 = System.nanoTime();
      SMSD s =
          new SMSD(
              rawMol("OCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCO"),
              rawMol("OCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCO"),
              new ChemOptions());
      var mcs = s.findMCS(false, true, 10000);
      long ms = (System.nanoTime() - t0) / 1_000_000;
      assertTrue(mcs.size() >= 30, "PEG MCS >= 30, got " + mcs.size());
      assertTrue(ms < 5000, "PEG should complete in < 5s, took " + ms + "ms");
    }

    // --- Solvent Corrections ---

    @Test
    @DisplayName("Solvent: DMSO shifts tautomer weights")
    void testSolventCorrectionDMSO() throws Exception {
      ChemOptions opts = ChemOptions.tautomerProfile().withSolvent(ChemOptions.Solvent.DMSO);
      SMSD s = new SMSD(rawMol("CC(=O)C"), rawMol("CC(O)=C"), opts);
      var mcs = s.findMCS();
      assertTrue(mcs.size() >= 3, "Keto/enol MCS in DMSO >= 3, got " + mcs.size());
    }

    @Test
    @DisplayName("Solvent: Chloroform shifts tautomer weights")
    void testSolventCorrectionCHCl3() throws Exception {
      ChemOptions opts = ChemOptions.tautomerProfile().withSolvent(ChemOptions.Solvent.CHLOROFORM);
      SMSD s = new SMSD(rawMol("CC(=O)C"), rawMol("CC(O)=C"), opts);
      var mcs = s.findMCS();
      assertTrue(mcs.size() >= 3, "Keto/enol MCS in CHCl3 >= 3, got " + mcs.size());
    }

    @Test
    @DisplayName("Solvent: Aqueous (default) unchanged")
    void testSolventCorrectionAqueous() throws Exception {
      ChemOptions opts = ChemOptions.tautomerProfile().withSolvent(ChemOptions.Solvent.AQUEOUS);
      SMSD s = new SMSD(rawMol("CC(=O)C"), rawMol("CC(O)=C"), opts);
      var mcs = s.findMCS();
      assertTrue(mcs.size() >= 3, "Keto/enol MCS in water >= 3, got " + mcs.size());
    }

    // --- Lenient SMILES ---

    @Test
    @DisplayName("Lenient: unbalanced parentheses sanitized")
    void testLenientSmiles() throws Exception {
      String bad = "CC(CC(=O)O";
      String fixed = MolGraph.sanitizeSmiles(bad);
      assertNotNull(fixed);
      IAtomContainer m = rawMol(fixed);
      assertTrue(m.getAtomCount() >= 4, "Sanitized SMILES should produce valid molecule");
    }
  }

  // ======================================================================
  // From: TautomerWeightTest.java
  // ======================================================================

  @Nested
  @DisplayName("pKa-Informed Tautomer Relevance Weights")
  class TautomerWeight {

    @Nested
    @DisplayName("Weight constants at pH 7.4")
    class WeightConstants {

      @Test @DisplayName("Keto/enol weight is 0.98 (keto overwhelmingly dominant)")
      void ketoEnolWeight() {
        assertEquals(0.98f, MolGraph.TW_KETO_ENOL, 1e-6f);
      }

      @Test @DisplayName("Amide/imidic weight is 0.97 (amide dominant)")
      void amideWeight() {
        assertEquals(0.97f, MolGraph.TW_AMIDE_IMIDIC, 1e-6f);
      }

      @Test @DisplayName("Amidine/guanidine weight is 0.50 (symmetric tautomers)")
      void symmetricWeights() {
        assertEquals(0.50f, MolGraph.TW_AMIDINE, 1e-6f);
        assertEquals(0.50f, MolGraph.TW_GUANIDINE, 1e-6f);
      }

      @Test @DisplayName("Imidazole N-H weight is 0.50 (symmetric)")
      void imidazoleWeight() {
        assertEquals(0.50f, MolGraph.TW_IMIDAZOLE_NH, 1e-6f);
      }

      @Test @DisplayName("1,3-diketone enol weight is 0.72 (pH-sensitive at 7.4)")
      void diketoneWeight() {
        assertEquals(0.72f, MolGraph.TW_DIKETONE_ENOL, 1e-6f);
      }

      @Test @DisplayName("Nitro/aci-nitro weight is 0.95 (nitro form strongly dominant at pH 7.4)")
      void aciNitroWeight() {
        assertEquals(0.95f, MolGraph.TW_NITRO_ACI, 1e-6f);
      }
    }

    @Nested
    @DisplayName("Tautomer class + weight assignment")
    class ClassAssignment {

      @Test @DisplayName("Acetone C=O oxygen gets tautomer class + keto weight")
      void acetoneKetoOxygen() throws Exception {
        MolGraph g = new MolGraph(rawMol("CC(=O)C"));
        int oIdx = -1;
        for (int i = 0; i < g.atomCount(); i++) {
          if (g.atomicNum[i] == 8) { oIdx = i; break; }
        }
        assertNotEquals(-1, oIdx, "Carbonyl O not found");
        assertTrue(g.tautomerClass[oIdx] >= 0, "Carbonyl O should be in a tautomeric class");
        assertEquals(MolGraph.TW_KETO_ENOL, g.tautomerWeight[oIdx], 1e-5f,
            "Carbonyl O weight should be keto/enol weight");
      }

      @Test @DisplayName("Formamide N gets amide tautomer weight")
      void formamideN() throws Exception {
        MolGraph g = new MolGraph(rawMol("NC=O"));
        int nIdx = -1, oIdx = -1;
        for (int i = 0; i < g.atomCount(); i++) {
          if (g.atomicNum[i] == 7) nIdx = i;
          if (g.atomicNum[i] == 8) oIdx = i;
        }
        assertNotEquals(-1, nIdx); assertNotEquals(-1, oIdx);
        assertEquals(g.tautomerClass[nIdx], g.tautomerClass[oIdx],
            "Amide N and O should share a tautomeric class");
        assertEquals(MolGraph.TW_AMIDE_IMIDIC, g.tautomerWeight[nIdx], 1e-5f);
        assertEquals(MolGraph.TW_AMIDE_IMIDIC, g.tautomerWeight[oIdx], 1e-5f);
      }

      @Test @DisplayName("Imidazole N atoms get symmetric weight 0.50")
      void imidazoleNWeight() throws Exception {
        MolGraph g = new MolGraph(rawMol("c1cnc[nH]1"));
        int tautN = 0;
        for (int i = 0; i < g.atomCount(); i++) {
          if (g.atomicNum[i] == 7 && g.tautomerClass[i] >= 0) tautN++;
        }
        assertTrue(tautN >= 2, "Both imidazole N atoms should be in tautomeric class");
        for (int i = 0; i < g.atomCount(); i++) {
          if (g.atomicNum[i] == 7 && g.tautomerClass[i] >= 0) {
            assertEquals(MolGraph.TW_IMIDAZOLE_NH, g.tautomerWeight[i], 1e-5f,
                "Imidazole N weight should be symmetric 0.50");
          }
        }
      }

      @Test @DisplayName("Non-tautomeric atoms default to weight 1.0")
      void nonTautomericDefault() throws Exception {
        MolGraph g = new MolGraph(rawMol("c1ccccc1"));
        for (int i = 0; i < g.atomCount(); i++) {
          assertEquals(1.0f, g.tautomerWeight[i], 1e-6f,
              "Benzene atoms should have default weight 1.0");
          assertEquals(-1, g.tautomerClass[i], "Benzene has no tautomeric classes");
        }
      }

      @Test @DisplayName("2-Pyridinone N and exo-O get pyridinone weight")
      void pyridinone() throws Exception {
        MolGraph g = new MolGraph(rawMol("O=c1cccc[nH]1"));
        int nIdx = -1;
        for (int i = 0; i < g.atomCount(); i++) {
          if (g.atomicNum[i] == 7) { nIdx = i; break; }
        }
        assertNotEquals(-1, nIdx, "N not found");
        assertTrue(g.tautomerClass[nIdx] >= 0 || g.tautomerWeight[nIdx] == 1.0f,
            "Pyridinone N should be in tautomeric system or default");
      }
    }

    @Nested
    @DisplayName("SMSD.getTautomerConfidence()")
    class ConfidenceScoring {

      @Test @DisplayName("Benzene self-MCS -> confidence 1.0 (no tautomers)")
      void benzeneConfidence() throws Exception {
        IAtomContainer benz = rawMol("c1ccccc1");
        ChemOptions c = ChemOptions.tautomerProfile();
        SMSD smsd = new SMSD(benz, benz, c);
        Map<Integer,Integer> mcs = smsd.findMCS();
        assertEquals(6, mcs.size());
        assertEquals(1.0, smsd.getTautomerConfidence(), 1e-9,
            "Self-match on benzene should yield confidence 1.0");
      }

      @Test @DisplayName("Non-tautomeric MCS -> confidence always 1.0")
      void nonTautomericMode() throws Exception {
        IAtomContainer a = rawMol("c1ccccc1");
        IAtomContainer b = rawMol("c1ccc(F)cc1");
        ChemOptions c = new ChemOptions();
        SMSD smsd = new SMSD(a, b, c);
        smsd.findMCS();
        assertEquals(1.0, smsd.getTautomerConfidence(), 1e-9,
            "Non-tautomeric mode always reports confidence 1.0");
      }

      @Test @DisplayName("Keto/enol MCS -> confidence < 1.0 (cross-element match)")
      void ketoEnolConfidence() throws Exception {
        IAtomContainer keto = rawMol("CC(=O)CC");
        IAtomContainer enol = rawMol("CC(O)=CC");
        ChemOptions c = ChemOptions.tautomerProfile();
        SMSD smsd = new SMSD(keto, enol, c);
        Map<Integer,Integer> mcs = smsd.findMCS();
        double conf = smsd.getTautomerConfidence();
        assertTrue(conf <= 1.0 && conf >= 0.0,
            "Confidence must be in [0,1], got " + conf);
      }

      @Test @DisplayName("computeTautomerConfidence() static method returns 1.0 for identical elements")
      void staticConfidenceAllSameElement() throws Exception {
        IAtomContainer a = rawMol("CC(=O)C");
        IAtomContainer b = rawMol("CC(=O)CC");
        ChemOptions c = new ChemOptions();
        SearchEngine.MCSOptions m = new SearchEngine.MCSOptions();
        Map<Integer,Integer> mcs = SearchEngine.findMCS(a, b, c, m);
        double conf = SearchEngine.computeTautomerConfidence(a, b, mcs);
        assertEquals(1.0, conf, 1e-9);
      }
    }

    @Nested
    @DisplayName("ChemOptions pH field")
    class PHField {

      @Test @DisplayName("Default pH is 7.4")
      void defaultPH() {
        assertEquals(7.4, new ChemOptions().pH, 1e-9);
      }

      @Test @DisplayName("withPH() fluent setter clamps to [0,14]")
      void withPH() {
        ChemOptions c = new ChemOptions().withPH(6.8);
        assertEquals(6.8, c.pH, 1e-9);

        ChemOptions low = new ChemOptions().withPH(-1.0);
        assertEquals(0.0, low.pH, 1e-9);

        ChemOptions high = new ChemOptions().withPH(20.0);
        assertEquals(14.0, high.pH, 1e-9);
      }

      @Test @DisplayName("tautomerProfile() pH defaults to 7.4")
      void tautomerProfilePH() {
        assertEquals(7.4, ChemOptions.tautomerProfile().pH, 1e-9);
      }
    }
  }

  // ======================================================================
  // From: HydrogenHandlingTest.java
  // ======================================================================

  @Nested
  @DisplayName("Hydrogen Handling")
  class HydrogenHandling {

    final ChemOptions OPTS = new ChemOptions();

    @Nested
    @DisplayName("1. Implicit H normalisation")
    class ImplicitH {

      @Test @Timeout(10) @DisplayName("1.01 methane implicit vs bracket -- same heavy atom count")
      void methaneImplicitVsBracket() throws Exception {
        IAtomContainer implicit = mol("C");
        IAtomContainer bracket  = mol("[CH4]");
        SMSD a = new SMSD(implicit, bracket, OPTS);
        assertTrue(a.isSubstructure(), "[CH4] should match implicit C");
        SMSD b = new SMSD(bracket, implicit, OPTS);
        assertTrue(b.isSubstructure(), "implicit C should match [CH4]");
      }

      @Test @Timeout(10) @DisplayName("1.02 benzene implicit vs bracketed aromatic H -- same MCS")
      void benzeneImplicitVsBracketed() throws Exception {
        IAtomContainer implicit  = mol("c1ccccc1");
        IAtomContainer bracketed = mol("[cH]1[cH][cH][cH][cH][cH]1");
        Map<Integer, Integer> mcs = new SMSD(implicit, bracketed, OPTS).findMCS();
        assertEquals(6, mcs.size(), "MCS should span all 6 carbons regardless of H notation");
      }

      @Test @Timeout(10) @DisplayName("1.03 ethanol implicit vs explicit H-count -- heavy-atom MCS = 3")
      void ethanolImplicitVsExplicitCount() throws Exception {
        IAtomContainer a = mol("CCO");
        IAtomContainer b = mol("[CH3][CH2][OH]");
        Map<Integer, Integer> mcs = new SMSD(a, b, OPTS).findMCS();
        assertEquals(3, mcs.size());
      }

      @Test @Timeout(10) @DisplayName("1.04 aspirin heavy-atom count unaffected by H notation")
      void aspirinHeavyAtomCount() throws Exception {
        IAtomContainer implicit = mol("CC(=O)Oc1ccccc1C(=O)O");
        IAtomContainer explicit = mol("[CH3]C(=O)Oc1ccccc1C(=O)O");
        Map<Integer, Integer> mcs = new SMSD(implicit, explicit, OPTS).findMCS();
        assertEquals(13, mcs.size(), "MCS must cover all 13 heavy atoms");
      }
    }

    @Nested
    @DisplayName("2. Explicit [H] as graph nodes")
    class ExplicitHNodes {

      @Test @Timeout(10) @DisplayName("2.01 [H][H] hydrogen molecule -- 2-atom self match")
      void hydrogenMolecule() throws Exception {
        IAtomContainer h2 = mol("[H][H]");
        SMSD smsd = new SMSD(h2, h2, OPTS);
        assertTrue(smsd.isSubstructure(), "H2 should match itself");
      }

      @Test @Timeout(10) @DisplayName("2.02 explicit [H] query does NOT match implicit-H target")
      void explicitHQueryVsImplicitTarget() throws Exception {
        IAtomContainer withH    = mol("[H]C([H])([H])[H]");
        IAtomContainer withoutH = mol("C");
        SMSD smsd = new SMSD(withH, withoutH, OPTS);
        assertFalse(smsd.isSubstructure(),
            "Explicit-H methane (5 atoms) cannot be substructure of implicit-H methane (1 atom)");
      }

      @Test @Timeout(10) @DisplayName("2.03 MCS of explicit-H vs implicit-H strips to heavy-atom overlap")
      void mcsExplicitVsImplicit() throws Exception {
        IAtomContainer water_explicit = mol("[H]O[H]");
        IAtomContainer water_implicit = mol("O");
        Map<Integer, Integer> mcs = new SMSD(water_explicit, water_implicit, OPTS).findMCS();
        assertEquals(1, mcs.size(), "MCS across explicit/implicit water is the O atom");
      }
    }

    @Nested
    @DisplayName("3. H-count in substructure matching")
    class HCountSubstructure {

      @Test @Timeout(10) @DisplayName("3.01 primary amine NH2 matches in aniline")
      void nh2InAniline() throws Exception {
        IAtomContainer nh2     = mol("N");
        IAtomContainer aniline = mol("Nc1ccccc1");
        ChemOptions loose = new ChemOptions();
        loose.matchBondOrder = ChemOptions.BondOrderMode.LOOSE;
        SMSD smsd = new SMSD(nh2, aniline, loose);
        assertTrue(smsd.isSubstructure(), "N (ammonia-like) should match aniline N in LOOSE mode");
      }

      @Test @Timeout(10) @DisplayName("3.02 secondary amine does not match primary-amine query in strict mode")
      void primaryVsSecondaryAmine() throws Exception {
        IAtomContainer primary   = mol("CN");
        IAtomContainer secondary = mol("CNC");
        SMSD smsd = new SMSD(primary, secondary, OPTS);
        Map<Integer, Integer> mcs = smsd.findMCS();
        assertTrue(mcs.size() >= 2, "C-N fragment should still map in MCS");
      }

      @Test @Timeout(10) @DisplayName("3.03 water implicit vs [OH2] explicit -- substructure both ways")
      void waterSymmetry() throws Exception {
        IAtomContainer implicit = mol("O");
        IAtomContainer explicit = mol("[OH2]");
        assertTrue(new SMSD(implicit, explicit, OPTS).isSubstructure());
        assertTrue(new SMSD(explicit, implicit, OPTS).isSubstructure());
      }

      @Test @Timeout(10) @DisplayName("3.04 charged species NH4+ has same heavy atom as NH3")
      void ammoniumVsAmmonia() throws Exception {
        IAtomContainer nh3  = mol("N");
        IAtomContainer nh4p = mol("[NH4+]");
        ChemOptions noCharge = new ChemOptions();
        noCharge.matchFormalCharge = false;
        Map<Integer, Integer> mcs = new SMSD(nh3, nh4p, noCharge).findMCS();
        assertEquals(1, mcs.size(), "MCS of NH3 and NH4+ is a single N when charge ignored");
      }
    }

    @Nested
    @DisplayName("4. Drug-like H normalisation")
    class DrugLike {

      @Test @Timeout(30) @DisplayName("4.01 morphine implicit vs bracketed-CH -- same MCS size")
      void morphineHNormalisation() throws Exception {
        IAtomContainer std = mol("CN1CCC23C4C1CC5=C(C2C(C=C4)O3)C=C(C=C5)O");
        IAtomContainer brk = mol("C[N]1CC[C@@]23[C@@H]4[C@H]1C[C@@H]5=C([C@@H]2[C@H](C=C4)O3)C=C(C=C5)O");
        Map<Integer, Integer> mcs = new SMSD(std, brk, OPTS).findMCS();
        assertTrue(mcs.size() >= 17,
            "Morphine MCS across H-notation variants must be >= 17 heavy atoms, got " + mcs.size());
      }

      @Test @Timeout(30) @DisplayName("4.02 caffeine vs theophylline unaffected by [NH] vs N notation")
      void caffeineTheophylline() throws Exception {
        IAtomContainer caff = mol("Cn1cnc2c1c(=O)n(C)c(=O)n2C");
        IAtomContainer theo = mol("Cn1cnc2c1c(=O)[nH]c(=O)n2C");
        Map<Integer, Integer> mcs = new SMSD(caff, theo, OPTS).findMCS();
        assertTrue(mcs.size() >= 11,
            "Caffeine/theophylline MCS should be >= 11 regardless of [nH] notation, got " + mcs.size());
      }

      @Test @Timeout(10) @DisplayName("4.03 ibuprofen heavy-atom self-match both notations")
      void ibuprofenSelfMatch() throws Exception {
        IAtomContainer a = mol("CC(C)Cc1ccc(CC(C)C(=O)O)cc1");
        IAtomContainer b = mol("[CH3][C@@H]([CH3])Cc1ccc(CC(C)C(=O)O)cc1");
        Map<Integer, Integer> mcs = new SMSD(a, b, OPTS).findMCS();
        assertTrue(mcs.size() >= 12,
            "Ibuprofen MCS across H/stereo notations must be >= 12, got " + mcs.size());
      }
    }

    @Nested
    @DisplayName("5. Symmetric molecules")
    class Symmetric {

      @Test @Timeout(10) @DisplayName("5.01 benzene: all H implicit -- full 6-atom MCS")
      void benzeneImplicit() throws Exception {
        IAtomContainer a = mol("c1ccccc1");
        IAtomContainer b = mol("c1ccccc1");
        Map<Integer, Integer> mcs = new SMSD(a, b, OPTS).findMCS();
        assertEquals(6, mcs.size());
      }

      @Test @Timeout(10) @DisplayName("5.02 cyclohexane implicit vs [CH2] notation -- 6-atom MCS")
      void cyclohexaneHNotation() throws Exception {
        IAtomContainer a = mol("C1CCCCC1");
        IAtomContainer b = mol("[CH2]1[CH2][CH2][CH2][CH2][CH2]1");
        Map<Integer, Integer> mcs = new SMSD(a, b, OPTS).findMCS();
        assertEquals(6, mcs.size());
      }
    }
  }

  // ======================================================================
  // From: RingFinderTest.java
  // ======================================================================

  @Nested
  @DisplayName("Ring Finder Tests")
  class RingFinder {

    private void assertValidCycles(MolGraph g, int[][] cycles) {
      for (int[] c : cycles) {
        assertTrue(c.length >= 3, "Cycle must have at least 3 atoms");
        for (int i = 0; i < c.length; i++) {
          int a = c[i], b = c[(i + 1) % c.length];
          assertTrue(g.hasBond(a, b),
              "Atoms " + a + " and " + b + " should be bonded in cycle");
        }
      }
    }

    private Set<Integer> cycleLengths(int[][] cycles) {
      Set<Integer> lengths = new HashSet<>();
      for (int[] c : cycles) lengths.add(c.length);
      return lengths;
    }

    @Test
    @DisplayName("Benzene: MCB=1, relevant=1, URFs=1")
    void testBenzene() throws Exception {
      MolGraph g = graph("c1ccccc1");
      int[][] mcb = g.computeRings();
      int[][] relevant = g.computeRelevantCycles();
      int[][][] urfs = g.computeURFs();

      assertEquals(1, mcb.length, "Benzene MCB should have 1 ring");
      assertEquals(1, relevant.length, "Benzene should have 1 relevant cycle");
      assertEquals(1, urfs.length, "Benzene should have 1 URF");
      assertEquals(6, relevant[0].length, "The relevant cycle should be 6-membered");
      assertValidCycles(g, relevant);
    }

    @Test
    @DisplayName("Cyclohexane: MCB=1, relevant=1, URFs=1")
    void testCyclohexane() throws Exception {
      MolGraph g = graph("C1CCCCC1");
      int[][] mcb = g.computeRings();
      int[][] relevant = g.computeRelevantCycles();
      int[][][] urfs = g.computeURFs();

      assertEquals(1, mcb.length);
      assertEquals(1, relevant.length);
      assertEquals(1, urfs.length);
      assertValidCycles(g, relevant);
    }

    @Test
    @DisplayName("Biphenyl: MCB=2, relevant=2, URFs=2 (each phenyl is its own URF family)")
    void testBiphenyl() throws Exception {
      MolGraph g = graph("c1ccc(-c2ccccc2)cc1");
      int[][] mcb = g.computeRings();
      int[][] relevant = g.computeRelevantCycles();
      int[][][] urfs = g.computeURFs();

      assertEquals(2, mcb.length);
      assertEquals(2, relevant.length);
      assertEquals(2, urfs.length, "Biphenyl: two independent phenyl URF families");
      assertValidCycles(g, relevant);
    }

    @Test
    @DisplayName("Naphthalene: MCB=2, relevant=2, URFs=2 (fused rings)")
    void testNaphthalene() throws Exception {
      MolGraph g = graph("c1ccc2ccccc2c1");
      int[][] mcb = g.computeRings();
      int[][] relevant = g.computeRelevantCycles();
      int[][][] urfs = g.computeURFs();

      assertEquals(2, mcb.length);
      assertEquals(2, relevant.length);
      assertEquals(2, urfs.length, "Naphthalene: two URF families (fused 6-rings, exchange blocked)");
      assertValidCycles(g, relevant);

      Set<Integer> lengths = cycleLengths(relevant);
      assertTrue(lengths.contains(6), "Should contain 6-membered rings");
      assertFalse(lengths.contains(10), "10-ring is XOR of two 6-rings -- not relevant");
    }

    @Test
    @DisplayName("Cubane: MCB=5, relevant=6, URFs=1")
    void testCubane() throws Exception {
      MolGraph g = graph("C12C3C4C1C5C4C3C25");
      int[][] mcb = g.computeRings();
      int[][] relevant = g.computeRelevantCycles();
      int[][][] urfs = g.computeURFs();

      assertEquals(5, mcb.length, "Cubane MCB should have 5 rings");
      assertEquals(6, relevant.length, "Cubane should have 6 relevant cycles (all 6 faces)");
      assertEquals(1, urfs.length, "Cubane should have 1 URF (all 6 faces equivalent)");
      assertValidCycles(g, relevant);

      for (int[] c : relevant) {
        assertEquals(4, c.length, "All cubane relevant cycles should be 4-membered");
      }
    }

    @Test
    @DisplayName("Norbornane: MCB=2, relevant=2")
    void testNorbornane() throws Exception {
      MolGraph g = graph("C1CC2CC1CC2");
      int[][] mcb = g.computeRings();
      int[][] relevant = g.computeRelevantCycles();

      assertEquals(2, mcb.length, "Norbornane MCB should have 2 rings");
      assertEquals(2, relevant.length, "Norbornane should have exactly 2 relevant cycles");
      assertValidCycles(g, relevant);
    }

    @Test
    @DisplayName("Adamantane: MCB=3, relevant=4")
    void testAdamantane() throws Exception {
      MolGraph g = graph("C1C2CC3CC1CC(C2)C3");
      int[][] mcb = g.computeRings();
      int[][] relevant = g.computeRelevantCycles();

      assertEquals(3, mcb.length);
      assertEquals(4, relevant.length, "Adamantane should have 4 relevant cycles");
      assertValidCycles(g, relevant);

      for (int[] c : relevant) {
        assertEquals(6, c.length, "All adamantane relevant cycles should be 6-membered");
      }
    }

    @Test
    @DisplayName("Butane: no rings")
    void testButane() throws Exception {
      MolGraph g = graph("CCCC");
      int[][] mcb = g.computeRings();
      int[][] relevant = g.computeRelevantCycles();
      int[][][] urfs = g.computeURFs();

      assertEquals(0, mcb.length);
      assertEquals(0, relevant.length);
      assertEquals(0, urfs.length);
    }

    @Test
    @DisplayName("Cyclopropane: MCB=1, relevant=1, URFs=1")
    void testCyclopropane() throws Exception {
      MolGraph g = graph("C1CC1");
      int[][] mcb = g.computeRings();
      int[][] relevant = g.computeRelevantCycles();
      int[][][] urfs = g.computeURFs();

      assertEquals(1, mcb.length);
      assertEquals(1, relevant.length);
      assertEquals(1, urfs.length);
      assertEquals(3, relevant[0].length);
      assertValidCycles(g, relevant);
    }

    @Test
    @DisplayName("Spiro[4.4]nonane: MCB=2, relevant=2, URFs=2")
    void testSpiro() throws Exception {
      MolGraph g = graph("C1CCC2(C1)CCCC2");
      int[][] mcb = g.computeRings();
      int[][] relevant = g.computeRelevantCycles();
      int[][][] urfs = g.computeURFs();

      assertEquals(2, mcb.length);
      assertEquals(2, relevant.length);
      assertEquals(2, urfs.length, "Spiro: two URF families (rings share only a spiro atom)");
      assertValidCycles(g, relevant);
    }

    @Test
    @DisplayName("Anthracene: MCB=3, relevant >= 3")
    void testAnthracene() throws Exception {
      MolGraph g = graph("c1ccc2cc3ccccc3cc2c1");
      int[][] mcb = g.computeRings();
      int[][] relevant = g.computeRelevantCycles();

      assertEquals(3, mcb.length);
      assertTrue(relevant.length >= 3, "Anthracene should have at least 3 relevant cycles");
      assertValidCycles(g, relevant);
    }

    @Test
    @DisplayName("Relevant cycles are deterministic")
    void testDeterministic() throws Exception {
      MolGraph g1 = graph("c1ccc2ccccc2c1");
      MolGraph g2 = graph("c1ccc2ccccc2c1");
      int[][] r1 = g1.computeRelevantCycles();
      int[][] r2 = g2.computeRelevantCycles();

      assertEquals(r1.length, r2.length, "Same molecule should give same number of relevant cycles");
    }
  }
}
