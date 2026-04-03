/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
package com.bioinception.smsd;

import com.bioinception.smsd.core.*;
import com.bioinception.smsd.core.ChemOptions;
import com.bioinception.smsd.core.SMSD;
import com.bioinception.smsd.core.Standardiser;
import java.util.*;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.concurrent.Callable;
import java.util.concurrent.CyclicBarrier;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;
import java.util.stream.Stream;
import org.junit.jupiter.api.*;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Nested;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.Arguments;
import org.junit.jupiter.params.provider.MethodSource;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.aromaticity.ElectronDonation;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.ringsearch.RingSearch;
import org.openscience.cdk.silent.Atom;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import static org.junit.jupiter.api.Assertions.*;
import static org.junit.jupiter.params.provider.Arguments.arguments;

/**
 * Core SMSD Tests: consolidated from SMSDTest.java, SMSDRobustTest.java, SMSDCasesTest.java, SMSDComprehensiveTest.java.
 *
 * @author Syed Asad Rahman
 */
@DisplayName("Core SMSD Tests")
public class CoreTest {

  // ======================================================================
  // From: SMSDTest.java
  // ======================================================================


  private static final SmilesParser SP = new SmilesParser(DefaultChemObjectBuilder.getInstance());
  private static final String GOLDEN_843_EDUCT =
      "O=C(C#CC)N(C1=CC=C(OC(C)(C)C)C=C1)C(C(=O)NC(C)(C)C)C=2C=CC=CC2";
  private static final String GOLDEN_843_PRODUCT =
      "O=C1C=CC2(C=C1)C(=CC(=O)N2C(C=3C=CC=CC3)C(=O)NC(C)(C)C)C";

  private static boolean hasAny() {
    try {
      ChemOptions.BondOrderMode.valueOf("ANY");
      return true;
    } catch (IllegalArgumentException e) {
      return false;
    }
  }

  private static ChemOptions opts(
      boolean atomTypeStrict, boolean bondOrderStrict, boolean ringToRingOnly, boolean stereo) {
    ChemOptions c = new ChemOptions();
    c.matchAtomType = atomTypeStrict;
    c.matchBondOrder =
        bondOrderStrict
            ? ChemOptions.BondOrderMode.STRICT
            : (hasAny()
                ? ChemOptions.BondOrderMode.valueOf("ANY")
                : ChemOptions.BondOrderMode.LOOSE);
    // Ring parity policy: matched atoms/bonds must agree on ring membership.
    c.ringMatchesRingOnly = ringToRingOnly;
    c.useChirality = stereo;
    return c;
  }

  private static IAtomContainer mol(String smiles) {
    try {
      IAtomContainer m = SP.parseSmiles(smiles);
      // Perceive types, rings, and aromaticity so isInRing()/isAromatic() are reliable
      AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(m);
      Cycles.markRingAtomsAndBonds(m);
      @SuppressWarnings(
          "deprecation") // CDK 2.11 transitional; use Aromaticity.Model.Daylight in 2.12+
      boolean applied = new Aromaticity(ElectronDonation.daylight(), Cycles.all()).apply(m);
      return m;
    } catch (Exception e) {
      throw new AssertionError("Failed to parse SMILES: " + smiles, e);
    }
  }

  // Java 11: materialise the bonds iterable (no Iterable#toList())
  private static List<IBond> bonds(IAtomContainer m) {
    List<IBond> out = new ArrayList<>();
    for (IBond b : m.bonds()) out.add(b);
    return out;
  }

  private static int mappedBonds(IAtomContainer A, IAtomContainer B, Map<Integer, Integer> map) {
    int cnt = 0;
    for (IBond qb : bonds(A)) {
      int qi = A.indexOf(qb.getAtom(0));
      int qj = A.indexOf(qb.getAtom(1));
      Integer ti = map.get(qi), tj = map.get(qj);
      if (ti == null || tj == null) continue;
      IBond tb = B.getBond(B.getAtom(ti), B.getAtom(tj));
      if (tb != null) cnt++;
    }
    return cnt;
  }

  private static int mappedRingBonds(
      IAtomContainer A, IAtomContainer B, Map<Integer, Integer> map) {
    int cnt = 0;
    for (IBond qb : bonds(A)) {
      int qi = A.indexOf(qb.getAtom(0));
      int qj = A.indexOf(qb.getAtom(1));
      Integer ti = map.get(qi), tj = map.get(qj);
      if (ti == null || tj == null) continue;
      IBond tb = B.getBond(B.getAtom(ti), B.getAtom(tj));
      if (tb != null && tb.isInRing()) cnt++;
    }
    return cnt;
  }

  private static int countComponentsOnQuery(IAtomContainer q, Collection<Integer> mapped) {
    if (mapped.isEmpty()) return 0;
    List<Integer> idx = new ArrayList<>(mapped);
    Map<Integer, Integer> pos = new HashMap<>();
    for (int i = 0; i < idx.size(); i++) pos.put(idx.get(i), i);
    int n = idx.size(), comp = 0;
    boolean[] seen = new boolean[n];
    for (int i = 0; i < n; i++) {
      if (seen[i]) continue;
      comp++;
      Deque<Integer> dq = new ArrayDeque<>();
      dq.add(i);
      seen[i] = true;
      while (!dq.isEmpty()) {
        int u = dq.poll();
        int qu = idx.get(u);
        for (IBond b : bonds(q)) {
          int a = q.indexOf(b.getAtom(0));
          int c = q.indexOf(b.getAtom(1));
          if (a == qu || c == qu) {
            int vq = (a == qu ? c : a);
            Integer vi = pos.get(vq);
            if (vi != null && !seen[vi]) {
              seen[vi] = true;
              dq.add(vi);
            }
          }
        }
      }
    }
    return comp;
  }

  private static void assertNonEmpty(Map<Integer, Integer> map, String msg) {
    assertNotNull(map, msg + " (missing map)");
    assertTrue(!map.isEmpty(), msg + " (expected non-empty mapping)");
  }

  // ─────────────────────────────────────────────────────────────────────────────
  // A — quick sanity tests
  // ─────────────────────────────────────────────────────────────────────────────

  @Test
  @DisplayName("A.1  Ethanol is a substructure (quick sanity)")
  void a1_substructure_exists() {
    SMSD smsd = new SMSD(mol("CCO"), mol("CCNCCO"), new ChemOptions());
    assertTrue(smsd.isSubstructure(1200L));
  }

  @Test
  @DisplayName("A.2  Benzene is not a substructure of hexane")
  void a2_substructure_not_exists() {
    SMSD smsd = new SMSD(mol("c1ccccc1"), mol("CCCCCC"), opts(true, true, false, false));
    assertFalse(smsd.isSubstructure(1200L));
  }

  @Test
  @DisplayName("A.3  MCS(bz, naph) has 6 atoms")
  void a3_mcs_size_benzene_naphthalene() {
    SMSD smsd = new SMSD(mol("c1ccccc1"), mol("c1ccc2ccccc2c1"), opts(true, false, false, false));
    Map<Integer, Integer> m = smsd.findMCS(false, true, 2000L);
    assertNonEmpty(m, "benzene vs naphthalene");
    assertEquals(6, m.size(), "MCS(bz, naph) should be 6");
  }

  @Test
  @DisplayName("A.4  MolGraph builder accepts dense bond matrices")
  void a4_molgraph_builder_accepts_dense_bond_matrices() {
    MolGraph g =
        new MolGraph.Builder()
            .atomCount(2)
            .atomicNumbers(new int[] {6, 6})
            .ringFlags(new boolean[] {false, false})
            .aromaticFlags(new boolean[] {false, false})
            .neighbors(new int[][] {{1}, {0}})
            .bondOrders(new int[][] {{0, 2}, {2, 0}})
            .doubleBondStereo(new int[][] {{0, 1}, {1, 0}})
            .build();

    assertTrue(g.hasBond(0, 1));
    assertEquals("C=C", g.toCanonicalSmiles());
  }

  @Test
  @DisplayName("A.5  MolGraph builder rejects malformed adjacency and IDs")
  void a5_molgraph_builder_rejects_malformed_graphs() {
    assertThrows(
        IllegalArgumentException.class,
        () ->
            new MolGraph.Builder()
                .atomCount(2)
                .atomicNumbers(new int[] {6, 6})
                .neighbors(new int[][] {{1}, {}})
                .build());

    assertThrows(
        IllegalArgumentException.class,
        () ->
            new MolGraph.Builder()
                .atomCount(2)
                .atomicNumbers(new int[] {6, 6})
                .neighbors(new int[][] {{1, 1}, {0}})
                .build());

    assertThrows(
        IllegalArgumentException.class,
        () ->
            new MolGraph.Builder()
                .atomCount(2)
                .atomicNumbers(new int[] {6, 6})
                .neighbors(new int[][] {{1}, {0}})
                .bondOrders(new int[][] {{1}, {2}})
                .build());

    assertThrows(
        IllegalArgumentException.class,
        () ->
            new MolGraph.Builder()
                .atomCount(2)
                .atomicNumbers(new int[] {6, 6})
                .neighbors(new int[][] {{1}, {0}})
                .atomIds(new int[] {5, 5})
                .build());
  }

  @Test
  @DisplayName("A.6  Canonical hash distinguishes bond order and charge chemistry")
  void a6_molgraph_canonical_hash_tracks_bond_and_charge_chemistry() {
    MolGraph ethane = new MolGraph(mol("CC"));
    MolGraph ethene = new MolGraph(mol("C=C"));
    MolGraph methylamine = new MolGraph(mol("CN"));
    MolGraph methylammonium = new MolGraph(mol("C[NH3+]"));

    assertNotEquals(ethane.getCanonicalHash(), ethene.getCanonicalHash());
    assertNotEquals(methylamine.getCanonicalHash(), methylammonium.getCanonicalHash());
  }

  // ─────────────────────────────────────────────────────────────────────────────
  // B — named use-cases (easy to relate; small budgets)
  // ─────────────────────────────────────────────────────────────────────────────

  @Test
  @DisplayName("B.1  Benzene vs toluene: aromatic core survives")
  void b1_benzene_vs_toluene() {
    IAtomContainer a = mol("c1ccccc1"); // benzene
    IAtomContainer b = mol("Cc1ccccc1"); // toluene
    SMSD smsd = new SMSD(a, b, opts(true, true, false, false));
    Map<Integer, Integer> m = smsd.findMCS(false, true, 2000L);
    assertNonEmpty(m, "Phenyl core expected");
    assertTrue(m.size() >= 6, "At least the 6 aromatic carbons should match");
    assertTrue(countComponentsOnQuery(a, m.keySet()) <= 1, "Phenyl core should be connected");
  }

  @Test
  @DisplayName("B.2  Water vs benzene: empty under strict atom typing")
  void b2_water_vs_benzene_strict_atoms() {
    SMSD smsd = new SMSD(mol("O"), mol("c1ccccc1"), opts(true, true, false, false));
    Map<Integer, Integer> m = smsd.findMCS(false, true, 1500L);
    assertTrue(
        m == null || m.isEmpty(), "O vs benzene should yield empty core under strict typing");
  }

  @Test
  @DisplayName("B.3  Cyclohexene vs ethylbenzene: ring difference often collapses to a chain")
  void b3_ring_vs_ring_collapses_chain() {
    IAtomContainer a = mol("C1CCCC=C1CC"); // cyclohexene–ethyl
    IAtomContainer b = mol("CCc1ccccc1"); // ethylbenzene
    SMSD smsd = new SMSD(a, b, opts(true, true, false, false));
    Map<Integer, Integer> m = smsd.findMCS(false, true, 2500L);
    assertNonEmpty(m, "Expect some common acyclic chain");
    assertTrue(
        mappedRingBonds(a, b, m) == 0 || m.size() >= 4,
        "Likely a chain MCS; if ring bonds appear, ensure size ≥ 4");
  }

  @Test
  @DisplayName("B.4  Ring→ring policy blocks ring↔chain matches")
  void b4_ring_to_ring_blocks_chain() {
    IAtomContainer ring6 = mol("C1CCCCC1");
    IAtomContainer chain = mol("CCCCCC");
    SMSD smsd = new SMSD(ring6, chain, opts(true, true, true, false)); // ringOnly=true
    Map<Integer, Integer> m = smsd.findMCS(false, true, 1500L);
    assertTrue(
        m == null || m.isEmpty(), "Ring→ring policy should prevent mapping a ring to a chain");
  }

  @Test
  @DisplayName("B.5  Acetonitrile vs iminoethane: relaxed bond order grows the core")
  void b5_strict_vs_relaxed_bond_order() {
    IAtomContainer acn = mol("CC#N");
    IAtomContainer imi = mol("CC=N");
    SMSD sStrict = new SMSD(acn, imi, opts(true, true, false, false));
    SMSD sRelax = new SMSD(acn, imi, opts(true, false, false, false));
    Map<Integer, Integer> mStrict = sStrict.findMCS(false, true, 2000L);
    Map<Integer, Integer> mRelax = sRelax.findMCS(true, true, 2000L); // MCIS-like
    assertNonEmpty(mStrict, "Strict should find a non-empty core");
    assertNonEmpty(mRelax, "Relaxed should match at least the strict core");
    assertTrue(
        mappedBonds(acn, imi, mRelax) >= mappedBonds(acn, imi, mStrict),
        "Relaxing bond order (and MCIS-like) should not reduce matched bonds");
  }

  @Test
  @DisplayName("B.6  Benzyl bromide → benzyl alcohol: aromatic backbone conserved")
  void b6_reaction_backbone_conserved() {
    IAtomContainer rx = mol("c1ccc(cc1)CBr");
    IAtomContainer pr = mol("c1ccc(cc1)CO");
    SMSD smsd = new SMSD(rx, pr, opts(true, true, false, false));
    Map<Integer, Integer> m = smsd.findMCS(false, true, 2000L);
    assertNonEmpty(m, "Backbone conserved across reaction");
    assertTrue(m.size() >= 6, "Expect ring (6) mapped");
  }

  @Test
  @DisplayName("B.7  Aromatic ethers: MCS large enough to seed alignment")
  void b7_aromatic_ether_seed() {
    IAtomContainer a = mol("COc1ccccc1O");
    IAtomContainer b = mol("COc1ccc(OC)cc1");
    SMSD smsd = new SMSD(a, b, opts(true, true, false, false));
    Map<Integer, Integer> m = smsd.findMCS(false, true, 2000L);
    assertNonEmpty(m, "Aromatic ether MCS should exist");
    assertTrue(m.size() >= 7, "Seed large enough for 3D overlay");
  }

  @Test
  @DisplayName("B.8  Two families cluster tighter than cross-pairs")
  void b8_clustering_hint() {
    String[] arom = {"c1ccccc1", "Cc1ccccc1", "OCCc1ccccc1"};
    String[] alco = {"CCO", "CCCO", "CCCCO"};
    double thr = 0.60;
    Map<String, Integer> nAtoms = new HashMap<>();
    String[] all = concat(arom, alco);
    for (String s : all) nAtoms.put(s, mol(s).getAtomCount());

    int intraA = 0, intraB = 0, cross = 0;
    for (int i = 0; i < all.length; i++) {
      for (int j = i + 1; j < all.length; j++) {
        IAtomContainer x = mol(all[i]), y = mol(all[j]);
        ChemOptions chemOpts = opts(true, false, true, false);
        chemOpts.aromaticityMode = ChemOptions.AromaticityMode.STRICT;
        SMSD smsd = new SMSD(x, y, chemOpts);
        Map<Integer, Integer> m = smsd.findMCS(false, true, 1500L);
        boolean ok = (m != null && !m.isEmpty());
        double sim =
            ok ? ((double) m.size() / Math.min(nAtoms.get(all[i]), nAtoms.get(all[j]))) : 0.0;
        boolean edge = ok && sim >= thr;

        boolean xInA = contains(arom, all[i]);
        boolean yInA = contains(arom, all[j]);
        if (xInA && yInA) {
          if (edge) intraA++;
        } else if (!xInA && !yInA) {
          if (edge) intraB++;
        } else {
          if (edge) cross++;
        }
      }
    }
    assertTrue(intraA >= cross, "Aromatics should cluster at least as tight as cross-pairs");
    assertTrue(intraB >= cross, "Aliphatic alcohols should cluster at least as tight as cross-pairs");
  }

  // ─────────────────────────────────────────────────────────────────────────────
  // C — targeted extras (induced/MCCS, allowed ring→ring, timeout, stereo, etc.)
  // ─────────────────────────────────────────────────────────────────────────────

  @Test
  @DisplayName("C.1  MCIS (induced) is ≥ MCCS on a simple aliphatic pair")
  void c1_induced_vs_mccs() {
    IAtomContainer a = mol("CCCC"); // butane
    IAtomContainer b = mol("CC=CC"); // but-2-ene
    SMSD mccs = new SMSD(a, b, opts(true, true, false, false));
    SMSD mcis = new SMSD(a, b, opts(true, false, false, false)); // relaxed order for MCIS-like
    Map<Integer, Integer> m1 = mccs.findMCS(false, true, 1500L); // MCCS
    Map<Integer, Integer> m2 = mcis.findMCS(true, true, 1500L); // MCIS
    int s1 = (m1 == null ? 0 : m1.size());
    int s2 = (m2 == null ? 0 : m2.size());
    assertTrue(s2 >= s1, "MCIS should not be smaller than MCCS on simple pairs");
  }

  @Test
  @DisplayName("C.2  Ring→ring policy allows ring↔ring mapping (benzene vs kekulé benzene)")
  void c2_ring_to_ring_allows_ring_to_ring() {
    IAtomContainer ringA = mol("c1ccccc1");
    IAtomContainer ringB = mol("C1=CC=CC=C1");
    SMSD smsd = new SMSD(ringA, ringB, opts(true, true, true, false)); // ringOnly=true
    Map<Integer, Integer> m = smsd.findMCS(false, true, 1500L);
    assertTrue(m != null && !m.isEmpty(), "Ring→ring policy should allow ring↔ring matching");
  }

  @Test
  @DisplayName("C.3  Timeout budget is honoured (tiny budget)")
  void c3_timeout_honoured() {
    SMSD smsd =
        new SMSD(
            mol("c1ccc2ccccc2c1"), mol("c1cccc2ccc3cccc3c2c1"), opts(true, true, false, false));
    Map<Integer, Integer> m = smsd.findMCS(false, true, 1L); // 1 ms
    assertNotNull(m, "Should return (possibly empty) mapping without hanging");
  }

  @Test
  @Timeout(value = 15, unit = TimeUnit.SECONDS, threadMode = Timeout.ThreadMode.SEPARATE_THREAD)
  @DisplayName("C.3b  GOLDEN_843 MAX-style pair completes inside timeout budget")
  void c3b_golden843_max_style_pair_completes_within_timeout() {
    IAtomContainer educt = mol(GOLDEN_843_EDUCT);
    IAtomContainer product = mol(GOLDEN_843_PRODUCT);

    // Mirrors the benchmark flags: atomType=false, bondMatch=true,
    // ringMatch=false, ringSizeMatch=false (MAX algorithm).
    ChemOptions chem = new ChemOptions();
    chem.matchAtomType = false;
    chem.matchBondOrder = ChemOptions.BondOrderMode.STRICT;
    chem.ringMatchesRingOnly = false;
    chem.ringFusionMode = ChemOptions.RingFusionMode.IGNORE;

    SMSD smsd = new SMSD(educt, product, chem);

    long start = System.nanoTime();
    Map<Integer, Integer> mapping = smsd.findMCS(false, true, 10_000L);
    long elapsedMs = (System.nanoTime() - start) / 1_000_000L;

    assertNotNull(mapping, "GOLDEN_843 should return a mapping instead of hanging");
    assertTrue(
        elapsedMs < 12_000L,
        "GOLDEN_843 exceeded the 10s budget by too much: " + elapsedMs + "ms");
    assertTrue(
        mapping.size() >= 12,
        "GOLDEN_843 should keep a substantial MCS under MAX-style flags, got "
            + mapping.size());
  }

  @Test
  @Timeout(value = 30, unit = TimeUnit.SECONDS, threadMode = Timeout.ThreadMode.SEPARATE_THREAD)
  @DisplayName("C.3c  GOLDEN_843 shared molecules stay stable under parallel reuse")
  void c3c_golden843_shared_molecules_parallel_reuse_is_stable() throws Exception {
    final IAtomContainer educt = mol(GOLDEN_843_EDUCT);
    final IAtomContainer product = mol(GOLDEN_843_PRODUCT);

    ChemOptions chem = new ChemOptions();
    chem.matchAtomType = false;
    chem.matchBondOrder = ChemOptions.BondOrderMode.STRICT;
    chem.ringMatchesRingOnly = false;
    chem.ringFusionMode = ChemOptions.RingFusionMode.IGNORE;

    SearchEngine.MCSOptions baselineOpts = new SearchEngine.MCSOptions();
    baselineOpts.induced = false;
    baselineOpts.connectedOnly = true;
    baselineOpts.timeoutMs = 15_000L;

    SearchEngine.clearMolGraphCache();
    Map<Integer, Integer> baseline = SearchEngine.findMCS(educt, product, chem, baselineOpts);
    assertNotNull(baseline, "Baseline GOLDEN_843 call should return a mapping");
    assertTrue(
        baseline.size() >= 12,
        "Baseline GOLDEN_843 should keep a substantial MCS, got " + baseline.size());

    int expectedSize = baseline.size();
    int expectedMappedBonds = mappedBonds(educt, product, baseline);

    final int threads = 8;
    final int rounds = 4;
    for (int round = 0; round < rounds; round++) {
      SearchEngine.clearMolGraphCache();

      CyclicBarrier barrier = new CyclicBarrier(threads);
      ExecutorService pool = Executors.newFixedThreadPool(threads);
      try {
        List<Callable<Map<Integer, Integer>>> tasks = new ArrayList<>();
        for (int i = 0; i < threads; i++) {
          tasks.add(
              () -> {
                barrier.await();
                SearchEngine.MCSOptions opts = new SearchEngine.MCSOptions();
                opts.induced = false;
                opts.connectedOnly = true;
                opts.timeoutMs = 15_000L;
                return SearchEngine.findMCS(educt, product, chem, opts);
              });
        }

        List<Future<Map<Integer, Integer>>> futures = pool.invokeAll(tasks, 20, TimeUnit.SECONDS);
        for (int i = 0; i < futures.size(); i++) {
          Future<Map<Integer, Integer>> future = futures.get(i);
          assertFalse(future.isCancelled(), "Round " + round + " task " + i + " timed out");

          Map<Integer, Integer> mapping = future.get();
          assertNotNull(mapping, "Round " + round + " task " + i + " returned null");
          assertEquals(
              expectedSize,
              mapping.size(),
              "Round " + round + " task " + i + " returned an unstable MCS size");
          assertEquals(
              expectedMappedBonds,
              mappedBonds(educt, product, mapping),
              "Round " + round + " task " + i + " returned an unstable mapped-bond count");
        }
      } finally {
        pool.shutdownNow();
        assertTrue(
            pool.awaitTermination(5, TimeUnit.SECONDS),
            "Worker pool did not terminate cleanly");
      }
    }
  }

  @Test
  @DisplayName("C.4  Stereo sensitivity: trans vs cis reduces match when stereo is ON")
  void c4_stereo_sensitivity_trans_vs_cis() {
    IAtomContainer trans = mol("C/C=C\\C");
    IAtomContainer cis = mol("C/C=C/C");
    SMSD sOn = new SMSD(trans, cis, opts(true, true, false, true)); // stereo ON
    SMSD sOff = new SMSD(trans, cis, opts(true, true, false, false)); // stereo OFF
    Map<Integer, Integer> mOn = sOn.findMCS(false, true, 2000L);
    Map<Integer, Integer> mOff = sOff.findMCS(false, true, 2000L);
    int aOn = (mOn == null ? 0 : mOn.size());
    int aOff = (mOff == null ? 0 : mOff.size());
    assertTrue(aOff >= aOn, "With stereo OFF, mapping should be ≥ stereo ON");
  }

  @Test
  @DisplayName("C.5  Connected vs disconnected: dMCS may yield ≥ fragments")
  void c5_connected_vs_disconnected() {
    IAtomContainer a = mol("CCOC(C)=O");
    IAtomContainer b = mol("CCOC(C)C=O");
    SMSD cMCS = new SMSD(a, b, opts(true, true, false, false));
    SMSD dMCS = new SMSD(a, b, opts(true, true, false, false));
    Map<Integer, Integer> mc = cMCS.findMCS(false, true, 2000L); // connectedOnly
    Map<Integer, Integer> md = dMCS.findMCS(false, false, 2000L); // disconnected
    int fc = countComponentsOnQuery(a, (mc == null ? Collections.emptySet() : mc.keySet()));
    int fd = countComponentsOnQuery(a, (md == null ? Collections.emptySet() : md.keySet()));
    assertNotNull(md, "dMCS returns a mapping");
    assertTrue(fd >= fc, "dMCS may have ≥ fragments than cMCS");
  }

  @Test
  @DisplayName("C.6  Aromaticity flexible mode tolerates aromatic vs kekulé forms")
  void c6_aromaticity_flexible() {
    IAtomContainer arom = mol("c1ccccc1");
    IAtomContainer kekul = mol("C1=CC=CC=C1");
    ChemOptions flex = opts(true, true, false, false);
    flex.aromaticityMode = ChemOptions.AromaticityMode.FLEXIBLE;
    SMSD smsd = new SMSD(arom, kekul, flex);
    Map<Integer, Integer> m = smsd.findMCS(false, true, 1500L);
    assertTrue(m != null && m.size() >= 6, "Flexible mode should match the full benzene core");
  }

  @Test
  @DisplayName("C.7  Mapping integrity: values are unique")
  void c7_mapping_integrity() {
    SMSD smsd = new SMSD(mol("c1ccccc1"), mol("c1ccc2ccccc2c1"), opts(true, false, false, false));
    Map<Integer, Integer> m = smsd.findMCS(false, true, 2000L);
    assertNotNull(m);
    assertEquals(m.size(), new HashSet<>(m.values()).size(), "Target indices must be unique");
  }

  @Test
  @DisplayName("C.8  Bond-count monotonicity check on a second pair")
  void c8_bond_count_monotonicity_recap() {
    IAtomContainer a = mol("CC=CC"); // 2-butene
    IAtomContainer b = mol("CCCC"); // butane
    SMSD strict = new SMSD(a, b, opts(true, true, false, false));
    SMSD relax = new SMSD(a, b, opts(true, false, false, false));
    Map<Integer, Integer> ms = strict.findMCS(false, true, 2000L);
    Map<Integer, Integer> mr = relax.findMCS(true, true, 2000L);
    int bs = mappedBonds(a, b, (ms == null ? Collections.emptyMap() : ms));
    int br = mappedBonds(a, b, (mr == null ? Collections.emptyMap() : mr));
    assertTrue(br >= bs, "Relaxed (MCIS-like) should not reduce matched bonds vs strict");
  }

  // tiny utilities
  private static String[] concat(String[] a, String[] b) {
    String[] z = Arrays.copyOf(a, a.length + b.length);
    System.arraycopy(b, 0, z, a.length, b.length);
    return z;
  }

  private static boolean contains(String[] arr, String x) {
    for (String s : arr) if (Objects.equals(s, x)) return true;
    return false;
  }

  // ======================================================================
  // From: SMSDRobustTest.java
  // ======================================================================


  // -----------------------------
  // Case generators
  // -----------------------------

  /** Alkane chains: C×k is always a substructure of C×(k+2). */
  static Stream<Arguments> chainCases() {
    List<Arguments> out = new ArrayList<>();
    for (int k = 2; k <= 21; k++) {
      out.add(arguments("chain_" + k, "C".repeat(k), "C".repeat(k + 2), true));
    }
    return out.stream();
  }

  /** Aromatic: benzene in naphthalene (true); aromatic vs Kekulé benzene (true). */
  static Stream<Arguments> aromCases() {
    return Stream.of(
        arguments("arom_0", "c1ccccc1", "c1ccc2ccccc2c1", true),
        arguments("arom_1", "c1ccccc1", "C1=CC=CC=C1", true));
  }

  /** Stereo: no stereo flag enabled → topological match → both true. */
  static Stream<Arguments> stereoCases() {
    return Stream.of(
        arguments("stereo_0", "C/C=C\\C", "C/C=C\\C", true),
        arguments("stereo_1", "C/C=C\\C", "C/C=C/C", true));
  }

  /** Disconnected fragments: two separate C atoms found in ethane. */
  static Stream<Arguments> disconnCases() {
    return Stream.of(arguments("disc_0", "C.C", "CC", true));
  }

  /**
   * Functional group tests: fx_0: amide NC(=O) is in NCC(=O)NCCC → true fx_1: [O-]C=O vs OC=O →
   * false (matchFormalCharge=true by default: [O-] ≠ O)
   */
  static Stream<Arguments> functionalCases() {
    return Stream.of(
        arguments("fx_0", "NC(=O)", "NCC(=O)NCCC", true),
        arguments("fx_1", "[O-]C=O", "OC=O", false));
  }

  /**
   * Extra diverse cases. Expected results are chemically validated: - extra_g: CF3 (fluoroform) is
   * NOT a substructure of CC(C)(F)F (only 2 fluorines)
   */
  static Stream<Arguments> extraCases() {
    return Stream.of(
        arguments("extra_a", "c1ccccc1CC", "c1ccccc1CCCC", true),
        arguments("extra_b", "CCN(CC)CC", "CCN(CC)CCO", true),
        // Cyclohexane is NOT a substructure of cycloheptane with ringMatchesRingOnly=true:
        // a 6-membered ring cannot map into a 7-membered ring as an induced subgraph
        arguments("extra_c", "C1CCCCC1", "C1CCCCCC1", false),
        arguments("extra_d", "[13CH3]C", "CC", true),
        arguments("extra_e", "c1ncccc1", "c1ncccc1Cl", true),
        arguments("extra_f", "N1CCOCC1", "O=C(N1CCOCC1)C", true),
        arguments("extra_g", "FC(F)F", "CC(C)(F)F", false),
        arguments("extra_h", "O=S(=O)N", "O=S(=O)NCC", true),
        arguments("extra_i", "P(=O)(O)O", "OP(=O)(O)OCC", true),
        arguments("extra_j", "Clc1ccccc1", "Clc1ccc(Cl)cc1", true));
  }

  /** All SMILES-based substructure cases with expected boolean. */
  static Stream<Arguments> allSmilesCases() {
    return Stream.of(
            chainCases(),
            aromCases(),
            stereoCases(),
            disconnCases(),
            functionalCases(),
            extraCases())
        .flatMap(s -> s);
  }

  /** SMARTS query cases: must use the SMARTS constructor for correct matching. */
  static Stream<Arguments> smartsCases() {
    return Stream.of(
        // rec_0: recursive SMARTS for carboxyl carbon matches acetic acid
        arguments("rec_0", "[C;$(C(=O)O)]", "CC(=O)O", true),
        // rec_1: named predicate $isAmideN matches amide nitrogen in target
        arguments("rec_1", "[N;$isAmideN]", "CC(=O)NCC", true));
  }

  /** MCS subset: first 30 SMILES cases. */
  static Stream<Arguments> mcsSubset() {
    return allSmilesCases().limit(30);
  }

  // Removed: extraCasesRepeated() was running 90 duplicate tests (9 copies of 10 cases).

  // -----------------------------
  // Tests
  // -----------------------------

  @ParameterizedTest(name = "{0}")
  @MethodSource("allSmilesCases")
  @DisplayName("Substructure SMILES cases with chemical assertions")
  void test_substructure_smiles(String name, String q, String t, boolean expected)
      throws Exception {
    SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
    IAtomContainer queryMol = sp.parseSmiles(q);
    IAtomContainer targetMol = sp.parseSmiles(t);
    SMSD smsd = new SMSD(queryMol, targetMol, new ChemOptions());
    boolean ok = smsd.isSubstructure();
    assertEquals(expected, ok, name + ": '" + q + "' vs '" + t + "' expected " + expected);
  }

  // Removed: test_substructure_extra_repeated ran 90 duplicate tests identical to extraCases.

  @ParameterizedTest(name = "{0}")
  @MethodSource("smartsCases")
  @DisplayName("SMARTS query substructure (uses SMARTS constructor)")
  void test_substructure_smarts(String name, String smarts, String t, boolean expected)
      throws Exception {
    SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
    IAtomContainer targetMol = sp.parseSmiles(t);
    targetMol = Standardiser.standardise(targetMol, Standardiser.TautomerMode.NONE);
    // Expand named predicates (e.g., $isAmideN) before matching
    Standardiser.PredicateRegistry reg = new Standardiser.PredicateRegistry();
    String expanded = Standardiser.expandNamedPredicates(smarts, reg);
    SMSD smsd = new SMSD(expanded, targetMol, new ChemOptions());
    boolean ok = smsd.isSubstructure();
    assertEquals(
        expected,
        ok,
        name
            + ": SMARTS '"
            + smarts
            + "' (expanded: '"
            + expanded
            + "') vs '"
            + t
            + "' expected "
            + expected);
  }

  @ParameterizedTest(name = "{0}")
  @MethodSource("mcsSubset")
  @DisplayName("MCS on subset (induced, connected)")
  void test_mcs_runs_mcisd(String name, String q, String t, boolean ignoredExpected)
      throws Exception {
    SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
    IAtomContainer queryMol = sp.parseSmiles(q);
    IAtomContainer targetMol = sp.parseSmiles(t);
    SMSD smsd = new SMSD(queryMol, targetMol, new ChemOptions());
    Map<Integer, Integer> mr = smsd.findMCS(true, true, 2500L);
    assertNotNull(mr, name + ": MCS should not be null");
    // Query and target always share at least some atoms
    // For cases where both molecules share atom types, MCS should be non-empty
    // For cases where they don't (e.g., fx_1: [O-]C=O vs OC=O with charge matching),
    // MCS may be empty. Just verify it returned a valid map.
    assertNotNull(mr, name + ": MCS result should not be null");
  }

  // Removed: test_cli_json_shape_comment was a trivially-passing placeholder (assertTrue(true)).

  // ======================================================================
  // From: SMSDCasesTest.java
  // ======================================================================


  // Alias for local use (original code used instance field 'sp')
  private final org.openscience.cdk.smiles.SmilesParser sp = SP;

  @Nested
  @DisplayName("Core Substructure Tests")
  class SubstructureTests {

    @Test
    @DisplayName("should find a simple alkane substructure")
    void findSimpleAlkane() throws Exception {
      IAtomContainer q = sp.parseSmiles("CCC"); // Propane
      IAtomContainer t = sp.parseSmiles("CCCCCC"); // Hexane
      SMSD smsd =
          new SMSD(
              Standardiser.standardise(q, Standardiser.TautomerMode.NONE),
              Standardiser.standardise(t, Standardiser.TautomerMode.NONE),
              ChemOptions.profile("compat-substruct"));
      assertTrue(smsd.isSubstructure());
    }

    @Test
    @DisplayName("should not find a longer chain in a shorter one")
    void notFindLongerChain() throws Exception {
      IAtomContainer q = sp.parseSmiles("CCCCCCCC"); // Octane
      IAtomContainer t = sp.parseSmiles("CCCCCC"); // Hexane
      SMSD smsd =
          new SMSD(
              Standardiser.standardise(q, Standardiser.TautomerMode.NONE),
              Standardiser.standardise(t, Standardiser.TautomerMode.NONE),
              ChemOptions.profile("compat-substruct"));
      assertFalse(smsd.isSubstructure());
    }

    @Test
    @DisplayName("should find a heterocycle substructure")
    void findHeterocycle() throws Exception {
      IAtomContainer q = sp.parseSmiles("c1ncccc1"); // Pyridine
      IAtomContainer t = sp.parseSmiles("Clc1ncccc1"); // 2-Chloropyridine
      SMSD smsd =
          new SMSD(
              Standardiser.standardise(q, Standardiser.TautomerMode.NONE),
              Standardiser.standardise(t, Standardiser.TautomerMode.NONE),
              ChemOptions.profile("compat-substruct"));
      assertTrue(smsd.isSubstructure());
    }

    @Test
    @DisplayName("should match formal charges when required")
    void matchFormalCharge() throws Exception {
      IAtomContainer q = sp.parseSmiles("[O-]C=O"); // Carboxylate
      IAtomContainer t = sp.parseSmiles("CC(=O)[O-]"); // Acetate
      ChemOptions options = ChemOptions.profile("compat-substruct");
      options.matchFormalCharge = true;
      SMSD smsd =
          new SMSD(
              Standardiser.standardise(q, Standardiser.TautomerMode.NONE),
              Standardiser.standardise(t, Standardiser.TautomerMode.NONE),
              options);
      assertTrue(smsd.isSubstructure());
    }

    @Test
    @DisplayName("should ignore formal charges when not required")
    void ignoreFormalCharge() throws Exception {
      IAtomContainer q = sp.parseSmiles("OC=O"); // Carboxylic acid
      IAtomContainer t = sp.parseSmiles("CC(=O)[O-]"); // Acetate
      ChemOptions options = ChemOptions.profile("compat-substruct");
      options.matchFormalCharge = false;
      SMSD smsd =
          new SMSD(
              Standardiser.standardise(q, Standardiser.TautomerMode.NONE),
              Standardiser.standardise(t, Standardiser.TautomerMode.NONE),
              options);
      assertTrue(smsd.isSubstructure());
    }
  }

  @Nested
  @DisplayName("Advanced SMARTS Feature Tests")
  class SmartsFeatureTests {

    @Test
    @DisplayName("should handle recursive SMARTS for carboxyl carbon")
    void recursiveSmartsCarboxyl() throws Exception {

      String qSmarts = "[C;$(C(=O)O)]";
      IAtomContainer t = sp.parseSmiles("CC(=O)O");
      SMSD smsd =
          new SMSD(
              qSmarts,
              Standardiser.standardise(t, Standardiser.TautomerMode.NONE),
              ChemOptions.profile("compat-substruct"));
      assertTrue(smsd.isSubstructure());
    }

    @Test
    @DisplayName("should handle recursive SMARTS for amide nitrogen")
    void recursiveSmartsAmide() throws Exception {

      String qSmarts = "[N;$(NC=O)]";
      IAtomContainer t = sp.parseSmiles("CC(=O)NCC");
      SMSD smsd =
          new SMSD(
              qSmarts,
              Standardiser.standardise(t, Standardiser.TautomerMode.NONE),
              ChemOptions.profile("compat-substruct"));
      assertTrue(smsd.isSubstructure());
    }

    @Test
    @DisplayName("should match disconnected query fragments in a target")
    void disconnectedFragments() throws Exception {
      IAtomContainer q = sp.parseSmiles("C.C");
      IAtomContainer t = sp.parseSmiles("CCC");
      SMSD smsd =
          new SMSD(
              Standardiser.standardise(q, Standardiser.TautomerMode.NONE),
              Standardiser.standardise(t, Standardiser.TautomerMode.NONE),
              ChemOptions.profile("compat-substruct"));
      assertTrue(smsd.isSubstructure());
    }
  }

  @Nested
  @DisplayName("Stereo and Aromaticity Tests")
  class StereoAndAromaticityTests {

    @Test
    @DisplayName("should match identical stereoisomers when stereo is ignored")
    void stereoMatch() throws Exception {
      IAtomContainer q = sp.parseSmiles("C/C=C\\C");
      IAtomContainer t = sp.parseSmiles("C/C=C\\C");
      ChemOptions options = ChemOptions.profile("strict");
      options.useBondStereo = false; // Disabling stereo check to allow topological pass
      SMSD smsd =
          new SMSD(
              Standardiser.standardise(q, Standardiser.TautomerMode.NONE),
              Standardiser.standardise(t, Standardiser.TautomerMode.NONE),
              options);
      assertTrue(smsd.isSubstructure());
    }

    @Test
    @DisplayName("should not match different stereoisomers with strict settings")
    void stereoMismatch() throws Exception {
      IAtomContainer q = sp.parseSmiles("C/C=C\\C");
      IAtomContainer t = sp.parseSmiles("C/C=C/C");
      ChemOptions options = ChemOptions.profile("strict");
      options.useBondStereo = true;
      SMSD smsd =
          new SMSD(
              Standardiser.standardise(q, Standardiser.TautomerMode.NONE),
              Standardiser.standardise(t, Standardiser.TautomerMode.NONE),
              options);
      assertFalse(smsd.isSubstructure());
    }

    @Test
    @DisplayName("should match aromatic and Kekule forms with flexible settings")
    void aromaticFlexibleMatch() throws Exception {
      IAtomContainer q = sp.parseSmiles("c1ccccc1");
      IAtomContainer t = sp.parseSmiles("C1=CC=CC=C1");
      SMSD smsd =
          new SMSD(
              Standardiser.standardise(q, Standardiser.TautomerMode.NONE),
              Standardiser.standardise(t, Standardiser.TautomerMode.NONE),
              ChemOptions.profile("compat-substruct"));
      assertTrue(smsd.isSubstructure());
    }

    @Test
    @DisplayName("should not match aromatic and non-aromatic rings with strict settings")
    void aromaticStrictMismatch() throws Exception {
      IAtomContainer q = sp.parseSmiles("c1ccccc1");
      IAtomContainer t = sp.parseSmiles("C1CCCCC1"); // Cyclohexane
      SMSD smsd =
          new SMSD(
              Standardiser.standardise(q, Standardiser.TautomerMode.NONE),
              Standardiser.standardise(t, Standardiser.TautomerMode.NONE),
              ChemOptions.profile("strict"));
      assertFalse(smsd.isSubstructure());
    }
  }

  @Nested
  @DisplayName("Maximum Common Substructure (MCS) Tests")
  class MCSTests {

    @Test
    @DisplayName("should find benzene as the MCS of benzene and naphthalene")
    void mcsBenzeneNaphthalene() throws Exception {
      IAtomContainer a = sp.parseSmiles("c1ccccc1");
      IAtomContainer b = sp.parseSmiles("c1ccc2ccccc2c1");
      SMSD smsd =
          new SMSD(
              Standardiser.standardise(a, Standardiser.TautomerMode.NONE),
              Standardiser.standardise(b, Standardiser.TautomerMode.NONE),
              ChemOptions.profile("compat-fmcs"));
      Map<Integer, Integer> m = smsd.findMCS(true, true, 2500L);
      assertNotNull(m);
      assertEquals(6, m.size());
    }

    @Test
    @DisplayName("should find a common core between two different structures")
    void mcsCommonCore() throws Exception {
      IAtomContainer a = sp.parseSmiles("CC(C)c1ccc(C(=O)O)cc1"); // Ibuprofen core
      IAtomContainer b = sp.parseSmiles("CC(C(=O)O)c1ccc2ccccc2c1"); // Naproxen core
      SMSD smsd =
          new SMSD(
              Standardiser.standardise(a, Standardiser.TautomerMode.NONE),
              Standardiser.standardise(b, Standardiser.TautomerMode.NONE),
              ChemOptions.profile("compat-fmcs"));
      Map<Integer, Integer> m = smsd.findMCS(true, true, 2500L);
      assertNotNull(m);
      assertTrue(m.size() >= 4);
    }
  }

  @Nested
  @DisplayName("Drug and Real-World Molecule Tests")
  class DrugMoleculeTests {

    @Test
    @DisplayName("should find the salicylic acid core in Aspirin")
    void testSalicylicAcidInAspirin() throws Exception {
      IAtomContainer query = sp.parseSmiles("c1(C(=O)O)ccccc1O"); // Salicylic acid
      IAtomContainer target = sp.parseSmiles("CC(=O)Oc1ccccc1C(=O)O"); // Aspirin
      SMSD smsd =
          new SMSD(
              Standardiser.standardise(query, Standardiser.TautomerMode.NONE),
              Standardiser.standardise(target, Standardiser.TautomerMode.NONE),
              ChemOptions.profile("compat-substruct"));
      assertTrue(smsd.isSubstructure(), "Salicylic acid is a substructure of Aspirin");
    }

    @Test
    @DisplayName("should find the catechol pharmacophore in Dopamine")
    void testCatecholInDopamine() throws Exception {
      IAtomContainer query = sp.parseSmiles("c1ccc(O)c(O)c1"); // Catechol
      IAtomContainer target = sp.parseSmiles("C1=CC(=C(C=C1CCN)O)O"); // Dopamine
      SMSD smsd =
          new SMSD(
              Standardiser.standardise(query, Standardiser.TautomerMode.NONE),
              Standardiser.standardise(target, Standardiser.TautomerMode.NONE),
              ChemOptions.profile("compat-substruct"));
      assertTrue(smsd.isSubstructure(), "Catechol is a substructure of Dopamine");
    }

    @Test
    @DisplayName("should confirm a complex fragment is not present")
    void testIndoleNotInAspirin() throws Exception {
      IAtomContainer query = sp.parseSmiles("c1ccc2[nH]ccc2c1"); // Indole
      IAtomContainer target = sp.parseSmiles("CC(=O)Oc1ccccc1C(=O)O"); // Aspirin
      SMSD smsd =
          new SMSD(
              Standardiser.standardise(query, Standardiser.TautomerMode.NONE),
              Standardiser.standardise(target, Standardiser.TautomerMode.NONE),
              ChemOptions.profile("compat-substruct"));
      assertFalse(smsd.isSubstructure(), "Indole is not a substructure of Aspirin");
    }

    @Test
    @DisplayName("should find a large MCS between Morphine and Codeine")
    void testMcsMorphineVsCodeine() throws Exception {
      // Morphine (21 atoms) and Codeine (22 atoms) differ at the ether bridge:
      // Morphine has a direct C-O bridge; Codeine has C-O-CH3.
      // Induced MCS = 19 (bridge topology mismatch), non-induced = 20.
      IAtomContainer morphine = sp.parseSmiles("CN1CCC23C4C1CC5=C(C2C(C=C4)O3)C=C(C=C5)O");
      IAtomContainer codeine = sp.parseSmiles("CN1CCC23C4C1CC5=C(C2C(C=C4)OC3)C=C(C=C5)O");

      SMSD smsd =
          new SMSD(
              Standardiser.standardise(morphine, Standardiser.TautomerMode.NONE),
              Standardiser.standardise(codeine, Standardiser.TautomerMode.NONE),
              new ChemOptions());

      Map<Integer, Integer> mcsInduced = smsd.findMCS(true, true, 30000L);
      assertNotNull(mcsInduced);
      assertTrue(
          mcsInduced.size() >= 19,
          "Induced MCS of Morphine/Codeine should be >= 19 (got " + mcsInduced.size() + ")");

      Map<Integer, Integer> mcsNonInduced = smsd.findMCS(false, true, 30000L);
      assertNotNull(mcsNonInduced);
      assertTrue(
          mcsNonInduced.size() >= 19,
          "Non-induced MCS of Morphine/Codeine should be >= 19 (got " + mcsNonInduced.size() + ")");
    }
  }

  // ======================================================================
  // From: SMSDComprehensiveTest.java
  // ======================================================================


  // ======================================================================
  // 1. HETEROCYCLES (O, S, mixed)
  // ======================================================================

  @Nested
  @DisplayName("Heterocycle Substructure Tests")
  class HeterocycleTests {
    @Test
    void furanInBenzofuran() throws Exception {
      assertTrue(
          new SMSD(mol("c1ccoc1"), mol("c1ccc2occc2c1"), new ChemOptions()).isSubstructure());
    }

    @Test
    void thiopheneInBenzothiophene() throws Exception {
      assertTrue(
          new SMSD(mol("c1ccsc1"), mol("c1ccc2sccc2c1"), new ChemOptions()).isSubstructure());
    }

    @Test
    void imidazoleInHistidine() throws Exception {
      // Imidazole ring in histidine side chain
      assertTrue(
          new SMSD(mol("c1cnc[nH]1"), mol("NC(Cc1cnc[nH]1)C(=O)O"), new ChemOptions())
              .isSubstructure());
    }

    @Test
    void pyrimidineRingInCytosine() throws Exception {
      // Cytosine contains a pyrimidine ring — use the 1,3-diazine core
      SMSD smsd = new SMSD(mol("C=1C=NC=NC1"), mol("Nc1ccnc(=O)[nH]1"), new ChemOptions());
      Map<Integer, Integer> mcs = smsd.findMCS(false, true, 5000L);
      assertTrue(mcs.size() >= 4, "Should share pyrimidine atoms");
    }

    @Test
    void morpholineInDrug() throws Exception {
      assertTrue(
          new SMSD(mol("C1COCCN1"), mol("O=C(N1CCOCC1)c1ccccc1"), new ChemOptions())
              .isSubstructure());
    }
  }

  // ======================================================================
  // 2. BRIDGED, SPIRO, COMPLEX RING SYSTEMS
  // ======================================================================

  @Nested
  @DisplayName("Complex Ring Systems")
  class ComplexRingTests {
    @Test
    void norbornaneSubstructure() throws Exception {
      // Bicyclo[2.2.1]heptane — bridged ring
      SMSD smsd = new SMSD(mol("C1CC2CC1C2"), mol("C1CC2CC1CC2C"), new ChemOptions());
      // Both are bridged, should find common subgraph
      Map<Integer, Integer> mcs = smsd.findMCS(false, true, 5000L);
      assertNotNull(mcs);
      assertTrue(mcs.size() >= 4, "Bridged ring MCS should find common atoms");
    }

    @Test
    void spiroCompoundMCS() throws Exception {
      // Spiro[4.4]nonane
      SMSD smsd = new SMSD(mol("C1CCC2(C1)CCCC2"), mol("C1CCC2(C1)CCCCC2"), new ChemOptions());
      Map<Integer, Integer> mcs = smsd.findMCS(false, true, 5000L);
      assertTrue(mcs.size() >= 5, "Spiro MCS should share at least one ring");
    }
    // Removed: adamantaneSelfMatch — duplicated by AdversarialTest.adamantaneSelfMatch.
  }

  // ======================================================================
  // 3. FUNCTIONAL GROUPS — systematic substructure
  // ======================================================================

  @Nested
  @DisplayName("Functional Group Substructure")
  class FunctionalGroupTests {
    @Test
    void primaryAmineInAlanine() throws Exception {
      assertTrue(new SMSD(mol("CN"), mol("CC(N)C(=O)O"), new ChemOptions()).isSubstructure());
    }

    @Test
    void aldehydeInBenzaldehyde() throws Exception {
      assertTrue(new SMSD(mol("C=O"), mol("O=Cc1ccccc1"), new ChemOptions()).isSubstructure());
    }

    @Test
    void nitrileInAcetonitrile() throws Exception {
      assertTrue(new SMSD(mol("C#N"), mol("CC#N"), new ChemOptions()).isSubstructure());
    }

    @Test
    void nitroInNitrobenzene() throws Exception {
      assertTrue(
          new SMSD(mol("[N+](=O)[O-]"), mol("c1ccc([N+](=O)[O-])cc1"), new ChemOptions())
              .isSubstructure());
    }

    @Test
    void sulfonylInDMSO() throws Exception {
      assertTrue(new SMSD(mol("S=O"), mol("CS(=O)C"), new ChemOptions()).isSubstructure());
    }

    @Test
    void thioetherInMethionine() throws Exception {
      assertTrue(new SMSD(mol("CSC"), mol("CSCCC(N)C(=O)O"), new ChemOptions()).isSubstructure());
    }

    @Test
    void phenolNotInCyclohexanol() throws Exception {
      // Aromatic OH should not match aliphatic OH with strict aromaticity
      SMSD smsd = new SMSD(mol("Oc1ccccc1"), mol("OC1CCCCC1"), ChemOptions.profile("strict"));
      assertFalse(smsd.isSubstructure());
    }

    @Test
    void epoxideSubstructure() throws Exception {
      assertTrue(new SMSD(mol("C1OC1"), mol("C1OC1CC"), new ChemOptions()).isSubstructure());
    }
  }

  // ======================================================================
  // 4. REACTION / TRANSFORMATION MAPPING
  // ======================================================================

  @Nested
  @DisplayName("Reaction Transformation MCS")
  class ReactionTests {
    @Test
    void alcoholToAldehyde() throws Exception {
      // Oxidation: ethanol → acetaldehyde — share C-C backbone
      SMSD smsd = new SMSD(mol("CCO"), mol("CC=O"), new ChemOptions());
      Map<Integer, Integer> mcs = smsd.findMCS(false, true, 5000L);
      assertTrue(mcs.size() >= 2, "Oxidation preserves C-C backbone");
    }

    @Test
    void acidToAmide() throws Exception {
      // Amide bond formation: acetic acid → acetamide
      SMSD smsd = new SMSD(mol("CC(=O)O"), mol("CC(=O)N"), new ChemOptions());
      Map<Integer, Integer> mcs = smsd.findMCS(false, true, 5000L);
      assertTrue(mcs.size() >= 3, "Acid→amide preserves acyl group");
    }

    @Test
    void halogenationMCS() throws Exception {
      // Chlorination: benzene → chlorobenzene
      SMSD smsd = new SMSD(mol("c1ccccc1"), mol("Clc1ccccc1"), new ChemOptions());
      Map<Integer, Integer> mcs = smsd.findMCS(false, true, 5000L);
      assertEquals(6, mcs.size(), "Benzene ring preserved in chlorination");
    }

    @Test
    void demethylation() throws Exception {
      // Codeine → morphine (loss of methyl from ether)
      SMSD smsd =
          new SMSD(
              mol("CN1CCC23C4C1CC5=C(C2C(C=C4)OC3)C=C(C=C5)O"),
              mol("CN1CCC23C4C1CC5=C(C2C(C=C4)O3)C=C(C=C5)O"),
              new ChemOptions());
      Map<Integer, Integer> mcs = smsd.findMCS(false, true, 30000L);
      assertTrue(mcs.size() >= 19, "Demethylation preserves core scaffold");
    }

    @Test
    void ketoEnolTautomer() throws Exception {
      // Keto-enol: acetone ↔ propen-2-ol
      SMSD smsd = new SMSD(mol("CC(=O)C"), mol("CC(O)=C"), new ChemOptions());
      Map<Integer, Integer> mcs = smsd.findMCS(false, true, 5000L);
      assertTrue(mcs.size() >= 3, "Tautomers share backbone");
    }

    @Test
    void ringOpeningMCS() throws Exception {
      // Ring opening: cyclopropane vs propane — need ring constraint OFF
      ChemOptions opts = new ChemOptions();
      opts.ringMatchesRingOnly = false;
      SMSD smsd = new SMSD(mol("C1CC1"), mol("CCC"), opts);
      Map<Integer, Integer> mcs = smsd.findMCS(false, true, 5000L);
      assertTrue(mcs.size() >= 2, "Ring opening preserves C-C bonds when ring constraint off");
    }

    @Test
    void esterHydrolysis() throws Exception {
      // Ethyl acetate → acetic acid + ethanol (shared acyl)
      SMSD smsd = new SMSD(mol("CC(=O)OCC"), mol("CC(=O)O"), new ChemOptions());
      Map<Integer, Integer> mcs = smsd.findMCS(false, true, 5000L);
      assertTrue(mcs.size() >= 3, "Ester hydrolysis preserves acyl group");
    }
  }

  // ======================================================================
  // 5. CHARGED SPECIES
  // ======================================================================

  @Nested
  @DisplayName("Charged Species")
  class ChargedTests {
    @Test
    void ammoniumNotMatchesAmine() throws Exception {
      // [NH4+] should NOT match NH3 with charge matching ON
      ChemOptions opts = new ChemOptions();
      opts.matchFormalCharge = true;
      assertFalse(new SMSD(mol("[NH4+]"), mol("N"), opts).isSubstructure());
    }

    @Test
    void ammoniumMatchesAmineChargeOff() throws Exception {
      ChemOptions opts = new ChemOptions();
      opts.matchFormalCharge = false;
      assertTrue(new SMSD(mol("[NH4+]"), mol("N"), opts).isSubstructure());
    }

    @Test
    void zwitterionSubstructure() throws Exception {
      // Glycine zwitterion
      ChemOptions opts = new ChemOptions();
      opts.matchFormalCharge = true;
      assertTrue(new SMSD(mol("[NH3+]CC([O-])=O"), mol("[NH3+]CC([O-])=O"), opts).isSubstructure());
    }
  }

  // ======================================================================
  // 6. CHEMOPTIONS CONFIGURATIONS
  // ======================================================================

  @Nested
  @DisplayName("ChemOptions Configuration Tests")
  class ChemOptionsTests {
    @Test
    void vf2MatcherEngine() throws Exception {
      ChemOptions opts = new ChemOptions();
      opts.matcherEngine = ChemOptions.MatcherEngine.VF2;
      assertTrue(new SMSD(mol("c1ccccc1"), mol("c1ccc(O)cc1"), opts).isSubstructure());
    }

    @Test
    void vf2ppMatcherEngine() throws Exception {
      ChemOptions opts = new ChemOptions();
      opts.matcherEngine = ChemOptions.MatcherEngine.VF2PP;
      assertTrue(new SMSD(mol("c1ccccc1"), mol("c1ccc(O)cc1"), opts).isSubstructure());
    }

    @Test
    void vf3MatcherEngine() throws Exception {
      ChemOptions opts = new ChemOptions();
      opts.matcherEngine = ChemOptions.MatcherEngine.VF3;
      assertTrue(new SMSD(mol("c1ccccc1"), mol("c1ccc(O)cc1"), opts).isSubstructure());
    }

    @Test
    void bondOrderAny() throws Exception {
      ChemOptions opts = new ChemOptions();
      opts.matchBondOrder = ChemOptions.BondOrderMode.ANY;
      // Single bond query should match double bond target
      assertTrue(new SMSD(mol("CC"), mol("C=C"), opts).isSubstructure());
    }

    @Test
    void pruningTogglesTwoHop() throws Exception {
      ChemOptions opts = new ChemOptions();
      opts.useTwoHopNLF = false;
      assertTrue(new SMSD(mol("c1ccccc1"), mol("c1ccc(O)cc1"), opts).isSubstructure());
    }

    @Test
    void pruningTogglesThreeHop() throws Exception {
      ChemOptions opts = new ChemOptions();
      opts.useThreeHopNLF = false;
      assertTrue(new SMSD(mol("c1ccccc1"), mol("c1ccc(O)cc1"), opts).isSubstructure());
    }

    @Test
    void pruningTogglesBitParallel() throws Exception {
      ChemOptions opts = new ChemOptions();
      opts.useBitParallelFeasibility = false;
      assertTrue(new SMSD(mol("c1ccccc1"), mol("c1ccc(O)cc1"), opts).isSubstructure());
    }

    @Test
    void ringMatchDisabled() throws Exception {
      ChemOptions opts = new ChemOptions();
      opts.ringMatchesRingOnly = false;
      // With ring constraint OFF, cyclopentane (5) should map into cyclohexane (6)
      // because ring atoms are allowed to match non-ring atoms
      SMSD smsd = new SMSD(mol("C1CCCC1"), mol("CCCCCC"), opts);
      Map<Integer, Integer> mcs = smsd.findMCS(false, true, 5000L);
      assertTrue(mcs.size() >= 4, "Ring match disabled allows ring-to-chain mapping in MCS");
    }

    @Test
    void fluentApi() throws Exception {
      ChemOptions opts =
          new ChemOptions()
              .withBondOrderMode(ChemOptions.BondOrderMode.LOOSE)
              .withAromaticityMode(ChemOptions.AromaticityMode.FLEXIBLE)
              .withTwoHopNLF(true)
              .withThreeHopNLF(false)
              .withBitParallelFeasibility(true);
      assertTrue(new SMSD(mol("c1ccccc1"), mol("c1ccc(O)cc1"), opts).isSubstructure());
    }
  }

  // ======================================================================
  // 7. EDGE CASES
  // ======================================================================

  @Nested
  @DisplayName("Edge Cases")
  class EdgeCaseTests {
    @Test
    void singleAtomQuery() throws Exception {
      assertTrue(new SMSD(mol("C"), mol("CCCC"), new ChemOptions()).isSubstructure());
    }

    @Test
    void singleAtomNoMatch() throws Exception {
      // Nitrogen not in all-carbon chain
      assertFalse(new SMSD(mol("N"), mol("CCCC"), new ChemOptions()).isSubstructure());
    }

    @Test
    void selfMatchSubstructure() throws Exception {
      IAtomContainer m = mol("c1ccc(O)cc1");
      assertTrue(new SMSD(m, m, new ChemOptions()).isSubstructure());
    }

    @Test
    void selfMatchMCS() throws Exception {
      IAtomContainer m = mol("c1ccc(O)cc1");
      Map<Integer, Integer> mcs = new SMSD(m, m, new ChemOptions()).findMCS(true, true, 5000L);
      assertEquals(m.getAtomCount(), mcs.size(), "Self-MCS should match all atoms");
    }

    @Test
    void noCommonAtomsMCS() throws Exception {
      // Methane (C) vs water (O) — no shared atoms with strict typing
      SMSD smsd = new SMSD(mol("C"), mol("O"), new ChemOptions());
      Map<Integer, Integer> mcs = smsd.findMCS(false, true, 2000L);
      assertTrue(mcs.isEmpty(), "No common atoms between methane and water");
    }

    @Test
    void disconnectedQuery() throws Exception {
      assertTrue(new SMSD(mol("C.N"), mol("CCN"), new ChemOptions()).isSubstructure());
    }

    @Test
    void disconnectedMCSInduced() throws Exception {
      SMSD smsd = new SMSD(mol("CC.OO"), mol("CCOCC"), new ChemOptions());
      Map<Integer, Integer> mcs = smsd.findMCS(true, false, 5000L);
      assertNotNull(mcs);
    }

    @Test
    void veryShortTimeout() throws Exception {
      SMSD smsd = new SMSD(mol("c1ccccc1"), mol("c1ccc2ccccc2c1"), new ChemOptions());
      // Should not hang; may or may not find result
      smsd.findMCS(true, true, 1L);
    }
  }

  // ======================================================================
  // 8. SMSD API METHODS
  // ======================================================================

  @Nested
  @DisplayName("SMSD API Coverage")
  class ApiTests {
    @Test
    void findAllSubstructuresCount() throws Exception {
      SMSD smsd = new SMSD(mol("C"), mol("CCCCCC"), new ChemOptions());
      List<Map<Integer, Integer>> maps = smsd.findAllSubstructures(100, 5000L);
      assertEquals(6, maps.size(), "Single C should match all 6 positions in hexane");
    }

    @Test
    void findAllSubstructuresMaxLimit() throws Exception {
      SMSD smsd = new SMSD(mol("C"), mol("CCCCCC"), new ChemOptions());
      List<Map<Integer, Integer>> maps = smsd.findAllSubstructures(3, 5000L);
      assertEquals(3, maps.size(), "Should respect maxSolutions limit");
    }

    @Test
    void isSubstructureWithStatsReturnsData() throws Exception {
      SMSD smsd = new SMSD(mol("c1ccccc1"), mol("c1ccc2ccccc2c1"), new ChemOptions());
      SearchEngine.SubstructureResult res = smsd.isSubstructureWithStats(5000L);
      assertTrue(res.exists());
      assertTrue(res.stats().nodesVisited() > 0, "Should report nodes visited");
      assertTrue(res.stats().timeMillis() >= 0, "Should report elapsed time");
    }

    @Test
    void findAllSubstructuresWithStatsReturnsData() throws Exception {
      SMSD smsd = new SMSD(mol("c1ccccc1"), mol("c1ccc(O)cc1"), new ChemOptions());
      SearchEngine.SubstructureResult res = smsd.findAllSubstructuresWithStats(100, 5000L);
      assertTrue(res.exists());
      assertFalse(res.mappings().isEmpty());
      assertTrue(res.stats().solutions() > 0);
    }

    @Test
    void timeoutSetters() throws Exception {
      SMSD smsd = new SMSD(mol("c1ccccc1"), mol("c1ccc2ccccc2c1"), new ChemOptions());
      smsd.setSubstructureTimeoutMs(100L);
      smsd.setMcsTimeoutMs(100L);
      // Should work with custom timeouts
      assertTrue(smsd.isSubstructure());
      assertNotNull(smsd.findMCS());
    }

    @Test
    void smartsConstructor() throws Exception {
      IAtomContainer target = mol("c1ccc(O)cc1");
      SMSD smsd = new SMSD("[OH]", target, new ChemOptions());
      assertTrue(smsd.isSubstructure(), "SMARTS [OH] should match phenol");
      assertThrows(UnsupportedOperationException.class, smsd::getQuery, "SMARTS query should throw on getQuery()");
    }

    @Test
    void smartsUnsupportedForEnumeration() throws Exception {
      IAtomContainer target = mol("CCO");
      SMSD smsd = new SMSD("[OH]", target, new ChemOptions());
      assertThrows(UnsupportedOperationException.class, () -> smsd.findAllSubstructures(10, 1000L));
    }

    @Test
    void smilesConvenienceConstructor() throws Exception {
      SMSD smsd = new SMSD("c1ccccc1", "c1ccc(O)cc1", new ChemOptions());
      assertTrue(smsd.isSubstructure());
      assertNotNull(smsd.getQuery());
      assertNotNull(smsd.getTarget());
    }
  }

  // ======================================================================
  // 9. NAMED SMARTS PREDICATES
  // ======================================================================

  @Nested
  @DisplayName("Named SMARTS Predicates")
  class PredicateTests {
    private boolean matchPredicate(String sig, String targetSmi) throws Exception {
      Standardiser.PredicateRegistry reg = new Standardiser.PredicateRegistry();
      String expanded = Standardiser.expandNamedPredicates(sig, reg);
      IAtomContainer t = mol(targetSmi);
      return !Standardiser.matchAll(expanded, t).isEmpty();
    }

    @Test
    void isKetone() throws Exception {
      assertTrue(matchPredicate("[C;$isKetone]", "CC(=O)C"));
    }

    @Test
    void isKetoneNotAldehyde() throws Exception {
      assertFalse(matchPredicate("[C;$isKetone]", "CC=O"));
    }

    @Test
    void isCarboxylC() throws Exception {
      assertTrue(matchPredicate("[C;$isCarboxylC]", "CC(=O)O"));
    }

    @Test
    void isEster() throws Exception {
      assertTrue(matchPredicate("[C;$isEster]", "CC(=O)OC"));
    }

    @Test
    void isHalogen() throws Exception {
      assertTrue(matchPredicate("[*;$isHalogen]", "CCCl"));
    }

    @Test
    void isHalogenBromine() throws Exception {
      assertTrue(matchPredicate("[*;$isHalogen]", "CCBr"));
    }
  }

  // ======================================================================
  // 10. MCS SCENARIOS
  // ======================================================================

  @Nested
  @DisplayName("MCS Scenarios")
  class MCSScenarioTests {
    @Test
    void inducedDisconnected() throws Exception {
      SMSD smsd = new SMSD(mol("CC.CC"), mol("CCCCCC"), new ChemOptions());
      Map<Integer, Integer> mcs = smsd.findMCS(true, false, 5000L);
      assertNotNull(mcs);
      assertTrue(mcs.size() >= 2);
    }

    @Test
    void mccsVsMcis() throws Exception {
      // MCCS (non-induced) should be >= MCIS (induced)
      SMSD smsd = new SMSD(mol("CCCC"), mol("CC=CC"), new ChemOptions());
      Map<Integer, Integer> mccs = smsd.findMCS(false, true, 5000L);
      Map<Integer, Integer> mcis = smsd.findMCS(true, true, 5000L);
      assertTrue(mccs.size() >= mcis.size(), "MCCS >= MCIS");
    }

    @Test
    void ringVsChainMCS() throws Exception {
      IAtomContainer ring = mol("C1CCCCC1");
      IAtomContainer chain = mol("CCCCC");

      ChemOptions strict = new ChemOptions();
      strict.ringMatchesRingOnly = true;
      Map<Integer, Integer> strictRingToChain = new SMSD(ring, chain, strict).findMCS(false, true, 5000L);
      Map<Integer, Integer> strictChainToRing = new SMSD(chain, ring, strict).findMCS(false, true, 5000L);
      assertTrue(strictRingToChain.isEmpty(), "Ring-constrained MCS should reject ring->chain");
      assertTrue(strictChainToRing.isEmpty(), "Ring-constrained MCS should reject chain->ring");

      ChemOptions relaxed = new ChemOptions();
      relaxed.ringMatchesRingOnly = false;
      Map<Integer, Integer> relaxedRingToChain = new SMSD(ring, chain, relaxed).findMCS(false, true, 5000L);
      Map<Integer, Integer> relaxedChainToRing = new SMSD(chain, ring, relaxed).findMCS(false, true, 5000L);
      assertEquals(5, relaxedRingToChain.size(), "Relaxed ring matching should recover the 5-atom path MCS");
      assertEquals(5, relaxedChainToRing.size(), "Relaxed ring matching should recover the 5-atom path MCS");
    }

    @Test
    void chainInRingSubstructureHonorsRingConstraint() throws Exception {
      IAtomContainer ring = mol("C1CCCCC1");
      IAtomContainer chain = mol("CCCCC");

      ChemOptions strict = new ChemOptions();
      strict.ringMatchesRingOnly = true;
      assertFalse(new SMSD(chain, ring, strict).isSubstructure(5000L),
          "Ring-constrained substructure should reject chain->ring matches");

      ChemOptions relaxed = new ChemOptions();
      relaxed.ringMatchesRingOnly = false;
      assertTrue(new SMSD(chain, ring, relaxed).isSubstructure(5000L),
          "Relaxed ring matching should allow chain->ring substructure matches");
    }

    @Test
    void largeMolMCS() throws Exception {
      // Aspirin vs ibuprofen — share aromatic ring
      SMSD smsd =
          new SMSD(
              mol("CC(=O)Oc1ccccc1C(=O)O"), // aspirin
              mol("CC(C)Cc1ccc(C(C)C(=O)O)cc1"), // ibuprofen
              new ChemOptions());
      Map<Integer, Integer> mcs = smsd.findMCS(false, true, 10000L);
      assertTrue(mcs.size() >= 6, "Should share at least a benzene ring");
    }
  }

  // ======================================================================
  // 11. DRUG MOLECULES — real-world pairs
  // ======================================================================

  @Nested
  @DisplayName("Drug Molecule Pairs")
  class DrugPairTests {
    @Test
    void caffeineContainsPurine() throws Exception {
      // Purine core in caffeine
      assertTrue(
          new SMSD(
                  mol("c1ncnc2[nH]cnc12"), // purine
                  mol("Cn1c(=O)c2c(ncn2C)n(C)c1=O"), // caffeine
                  ChemOptions.profile("compat-substruct"))
              .isSubstructure());
    }

    @Test
    void diazepamVsClonazepam() throws Exception {
      // Benzodiazepine core shared
      SMSD smsd =
          new SMSD(
              mol("c1ccc2c(c1)C(=NCC(=O)N2)c1ccccc1Cl"), // diazepam
              mol("c1cc2c(cc1[N+](=O)[O-])C(=NCC(=O)N2)c1ccccc1Cl"), // clonazepam
              new ChemOptions());
      Map<Integer, Integer> mcs = smsd.findMCS(false, true, 15000L);
      assertTrue(mcs.size() >= 15, "Benzodiazepine core should be shared");
    }

    @Test
    void ibuprofenVsNaproxen() throws Exception {
      SMSD smsd =
          new SMSD(
              mol("CC(C)Cc1ccc(C(C)C(=O)O)cc1"), // ibuprofen
              mol("COc1ccc2cc(C(C)C(=O)O)ccc2c1"), // naproxen
              new ChemOptions());
      Map<Integer, Integer> mcs = smsd.findMCS(false, true, 10000L);
      assertTrue(mcs.size() >= 8, "Share propionic acid + aromatic");
    }
  }

  // ======================================================================
  // 12. ROBUSTNESS & EDGE CASE TESTS
  // ======================================================================

  @Nested
  @DisplayName("Robustness & Edge Case Tests")
  class RobustnessEdgeCaseTests {

    @Test
    @DisplayName("Null query molecule throws NullPointerException")
    void nullQueryMolecule() {
      IAtomContainer target = assertDoesNotThrow(() -> mol("C"));
      ChemOptions opts = new ChemOptions();
      assertThrows(NullPointerException.class, () -> new SMSD((IAtomContainer) null, target, opts));
    }

    @Test
    @DisplayName("Null target molecule throws NullPointerException")
    void nullTargetMolecule() {
      IAtomContainer query = assertDoesNotThrow(() -> mol("C"));
      ChemOptions opts = new ChemOptions();
      assertThrows(NullPointerException.class, () -> new SMSD(query, (IAtomContainer) null, opts));
    }

    @Test
    @DisplayName("Empty SMILES for both query and target does not throw")
    void emptySmilesBothMolecules() throws Exception {
      IAtomContainer q = SP.parseSmiles("");
      IAtomContainer t = SP.parseSmiles("");
      ChemOptions opts = new ChemOptions();
      SMSD smsd = new SMSD(q, t, opts);
      // Should not throw
      smsd.isSubstructure();
    }

    @Test
    @DisplayName("Single atom C in single atom C is substructure")
    void singleAtomCinC() throws Exception {
      assertTrue(new SMSD(mol("C"), mol("C"), new ChemOptions()).isSubstructure());
    }

    @Test
    @DisplayName("Single atom C in single atom N is not substructure")
    void singleAtomCinN() throws Exception {
      assertFalse(new SMSD(mol("C"), mol("N"), new ChemOptions()).isSubstructure());
    }

    @Test
    @DisplayName("Query larger than target: naphthalene not in benzene")
    void queryLargerThanTarget() throws Exception {
      assertFalse(
          new SMSD(mol("c1ccc2ccccc2c1"), mol("c1ccccc1"), new ChemOptions()).isSubstructure());
    }

    @Test
    @DisplayName("Timeout 1ms on large pair does not hang")
    void timeoutDoesNotHang() throws Exception {
      IAtomContainer q = mol("c1ccc2ccccc2c1"); // naphthalene
      IAtomContainer t = mol("c1ccc2cc3ccccc3cc2c1"); // anthracene
      ChemOptions opts = new ChemOptions();
      SMSD smsd = new SMSD(q, t, opts);
      smsd.setSubstructureTimeoutMs(1L);
      long start = System.currentTimeMillis();
      // May or may not find the result, but must not hang
      smsd.isSubstructure();
      long elapsed = System.currentTimeMillis() - start;
      assertTrue(elapsed < 5000, "Should complete within 5 seconds, took " + elapsed + "ms");
    }

    @Test
    @DisplayName("Molecules with explicit hydrogens: [H]O[H] water")
    void explicitHydrogensWater() throws Exception {
      IAtomContainer water = SP.parseSmiles("[H]O[H]");
      IAtomContainer target = SP.parseSmiles("[H]O[H]");
      ChemOptions opts = new ChemOptions();
      SMSD smsd = new SMSD(water, target, opts);
      // Should parse and check without error
      boolean result = smsd.isSubstructure();
      assertTrue(result, "Water should be substructure of itself");
    }

    @Test
    @DisplayName("3+ disconnected fragments: C.C.C in CCCC handled gracefully")
    void disconnectedFragments() throws Exception {
      IAtomContainer q = mol("C.C.C");
      IAtomContainer t = mol("CCCC");
      ChemOptions opts = new ChemOptions();
      SMSD smsd = new SMSD(q, t, opts);
      // Should handle gracefully without exception
      smsd.isSubstructure();
    }

    @Test
    @DisplayName("Very long chain: C*200 contains C*100 substructure, completes < 5s")
    void veryLongChainSubstructure() throws Exception {
      String longChain200 = "C".repeat(200);
      String longChain100 = "C".repeat(100);
      IAtomContainer q = SP.parseSmiles(longChain100);
      IAtomContainer t = SP.parseSmiles(longChain200);
      ChemOptions opts = new ChemOptions();
      SMSD smsd = new SMSD(q, t, opts);
      long start = System.currentTimeMillis();
      boolean result = smsd.isSubstructure();
      long elapsed = System.currentTimeMillis() - start;
      assertTrue(result, "C*100 should be substructure of C*200");
      assertTrue(elapsed < 5000, "Should complete within 5 seconds, took " + elapsed + "ms");
    }
  }

  // ======================================================================
  // 13. STANDARDISER & TAUTOMERMODE TESTS
  // ======================================================================

  @Nested
  @DisplayName("Standardiser & TautomerMode Tests")
  class StandardiserTautomerTests {

    @Test
    @DisplayName("Standardiser.standardise(null, NONE) throws NullPointerException")
    void standardiseNullThrows() {
      assertThrows(
          NullPointerException.class,
          () -> Standardiser.standardise(null, Standardiser.TautomerMode.NONE));
    }

    @Test
    @DisplayName("TautomerMode.NONE: standardise benzene preserves 6 atoms")
    void tautomerNoneBenzene() throws Exception {
      IAtomContainer benzene = SP.parseSmiles("c1ccccc1");
      IAtomContainer result = Standardiser.standardise(benzene, Standardiser.TautomerMode.NONE);
      assertEquals(6, result.getAtomCount(), "Benzene should have 6 atoms after standardisation");
    }

    @Test
    @DisplayName("Standardiser handles molecule without implicit H counts")
    void moleculeWithoutImplicitH() throws Exception {
      IAtomContainer mol = SilentChemObjectBuilder.getInstance().newAtomContainer();
      Atom c1 = new Atom("C");
      Atom c2 = new Atom("C");
      // Deliberately do NOT set implicit H count
      mol.addAtom(c1);
      mol.addAtom(c2);
      mol.addBond(0, 1, IBond.Order.SINGLE);
      // Should not throw
      IAtomContainer result = Standardiser.standardise(mol, Standardiser.TautomerMode.NONE);
      assertNotNull(result);
    }

    @Test
    @DisplayName("Standardiser preserves ring info for naphthalene")
    void preservesRingInfo() throws Exception {
      IAtomContainer naphthalene = SP.parseSmiles("c1ccc2ccccc2c1");
      IAtomContainer result = Standardiser.standardise(naphthalene, Standardiser.TautomerMode.NONE);
      RingSearch rs = new RingSearch(result);
      assertTrue(
          rs.numRings() >= 2, "Naphthalene should have at least 2 rings after standardisation");
    }

    @Test
    @DisplayName("Salt stripping: largest fragment of [Na+].[Cl-].c1ccccc1 has 6 carbons")
    void saltStripping() throws Exception {
      IAtomContainer saltMix = SP.parseSmiles("[Na+].[Cl-].c1ccccc1");
      IAtomContainer result = Standardiser.standardise(saltMix, Standardiser.TautomerMode.NONE);
      // Count carbon atoms in the result
      long carbonCount = 0;
      for (int i = 0; i < result.getAtomCount(); i++) {
        if ("C".equals(result.getAtom(i).getSymbol())) {
          carbonCount++;
        }
      }
      assertEquals(
          6, carbonCount, "After salt stripping, largest fragment should have 6 carbon atoms");
    }
  }

  // ======================================================================
  // Profiled MCS (Phase 0 Instrumentation)
  // ======================================================================

  @Nested
  @DisplayName("MCS Stage Profiling")
  class MCSStageProfilingTest {

    @Test
    @DisplayName("findMCSProfiledFromSmiles returns valid timers")
    void profiledFromSmiles_validTimers() throws Exception {
      var pr = SearchEngine.findMCSProfiledFromSmiles(
          "c1ccccc1", "c1ccc(O)cc1", new ChemOptions(), new SearchEngine.MCSOptions());

      assertNotNull(pr, "profiled result must not be null");
      assertNotNull(pr.result(), "MCSResult must not be null");
      assertNotNull(pr.timers(), "timers must not be null");

      var t = pr.timers();
      assertTrue(t.totalUs() > 0, "totalUs must be positive");
      assertTrue(t.orientationUs() >= 0, "orientationUs must be non-negative");
      assertTrue(t.seedsUs() >= 0, "seedsUs must be non-negative");
      assertTrue(t.mcSplitUs() >= 0, "mcSplitUs must be non-negative");
      assertTrue(t.bkUs() >= 0, "bkUs must be non-negative");
      assertTrue(t.mcGregorUs() >= 0, "mcGregorUs must be non-negative");
      assertTrue(t.repairUs() >= 0, "repairUs must be non-negative");

      // The mapping should be non-empty for benzene vs phenol
      assertTrue(pr.result().size() >= 6, "benzene/phenol MCS should be >= 6 atoms");
    }

    @Test
    @DisplayName("findMCSProfiled (MolGraph API) returns valid timers")
    void profiledMolGraph_validTimers() throws Exception {
      var sp = new SmilesParser(SilentChemObjectBuilder.getInstance());
      var mol1 = sp.parseSmiles("c1ccccc1");
      var mol2 = sp.parseSmiles("c1ccc(O)cc1");

      var pr = SearchEngine.findMCSProfiled(mol1, mol2, new ChemOptions(), new SearchEngine.MCSOptions());

      assertNotNull(pr);
      var t = pr.timers();
      assertTrue(t.totalUs() > 0, "totalUs must be positive");
      assertTrue(t.orientationUs() >= 0, "orientationUs must be non-negative");
      assertTrue(t.seedsUs() >= 0, "seedsUs must be non-negative");
      assertTrue(t.mcSplitUs() >= 0, "mcSplitUs must be non-negative");
      assertTrue(t.bkUs() >= 0, "bkUs must be non-negative");
      assertTrue(t.mcGregorUs() >= 0, "mcGregorUs must be non-negative");
      assertTrue(t.repairUs() >= 0, "repairUs must be non-negative");
    }

    @Test
    @DisplayName("bestAfter* fields are non-negative and consistent")
    void profiledBestAfterFields() throws Exception {
      var pr = SearchEngine.findMCSProfiledFromSmiles(
          "CC(=O)Oc1ccccc1C(O)=O", "CC(=O)Oc1ccccc1",
          new ChemOptions(), new SearchEngine.MCSOptions());

      var t = pr.timers();
      assertTrue(t.bestAfterGreedy() >= 0, "bestAfterGreedy must be non-negative");
      assertTrue(t.bestAfterSeed() >= 0, "bestAfterSeed must be non-negative");
      assertTrue(t.bestAfterBK() >= 0, "bestAfterBK must be non-negative");
      assertTrue(t.bestAfterMcGregor() >= 0, "bestAfterMcGregor must be non-negative");
      // The final result size should be >= the best seen at any intermediate stage
      int finalSize = pr.result().size();
      assertTrue(finalSize >= t.bestAfterGreedy(), "final >= bestAfterGreedy");
      assertTrue(finalSize >= t.bestAfterSeed(), "final >= bestAfterSeed");
    }
  }
}
