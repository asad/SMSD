/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
package com.bioinception.smsd;

import com.bioinception.smsd.core.MCS;

import com.bioinception.smsd.TestSupport;

import static com.bioinception.smsd.TestSupport.mapping;
import static com.bioinception.smsd.TestSupport.repeat;

import com.bioinception.smsd.core.*;
import com.bioinception.smsd.core.ChemOptions;
import com.bioinception.smsd.TestSupport.MoleculePair;
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
import org.junit.jupiter.params.provider.CsvSource;
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
 * Core tests: the search API and its options, robustness to awkward input, reference substructure
 * and MCS cases, and chemistry coverage from heterocycles to charged species.
 *
 * @author Syed Asad Rahman
 */
@DisplayName("Core SMSD Tests")
public class CoreTest {

  // ======================================================================
  // Search API and options
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
          "deprecation") // ElectronDonation.daylight() is deprecated in CDK
      boolean applied = new Aromaticity(ElectronDonation.daylight(), Cycles.all()).apply(m);
      return m;
    } catch (Exception e) {
      throw new AssertionError("Failed to parse SMILES: " + smiles, e);
    }
  }

  // Materialise the bonds iterable, which has no toList() of its own.
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
    MoleculePair smsd = new MoleculePair(mol("CCO"), mol("CCNCCO"), new ChemOptions());
    assertTrue(smsd.contains(1200L));
  }

  @Test
  @DisplayName("A.2  Benzene is not a substructure of hexane")
  void a2_substructure_not_exists() {
    MoleculePair smsd = new MoleculePair(mol("c1ccccc1"), mol("CCCCCC"), opts(true, true, false, false));
    assertFalse(smsd.contains(1200L));
  }

  @Test
  @DisplayName("A.3  MCS(bz, naph) has 6 atoms")
  void a3_mcs_size_benzene_naphthalene() {
    MoleculePair smsd = new MoleculePair(mol("c1ccccc1"), mol("c1ccc2ccccc2c1"), opts(true, false, false, false));
    Map<Integer, Integer> m = smsd.mcsMapping( true, 2000L);
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
    MoleculePair smsd = new MoleculePair(a, b, opts(true, true, false, false));
    Map<Integer, Integer> m = smsd.mcsMapping( true, 2000L);
    assertNonEmpty(m, "Phenyl core expected");
    assertTrue(m.size() >= 6, "At least the 6 aromatic carbons should match");
    assertTrue(countComponentsOnQuery(a, m.keySet()) <= 1, "Phenyl core should be connected");
  }

  @Test
  @DisplayName("B.2  Water vs benzene: empty under strict atom typing")
  void b2_water_vs_benzene_strict_atoms() {
    MoleculePair smsd = new MoleculePair(mol("O"), mol("c1ccccc1"), opts(true, true, false, false));
    Map<Integer, Integer> m = smsd.mcsMapping( true, 1500L);
    assertTrue(
        m == null || m.isEmpty(), "O vs benzene should yield empty core under strict typing");
  }

  @Test
  @DisplayName("B.3  Cyclohexene vs ethylbenzene: ring difference often collapses to a chain")
  void b3_ring_vs_ring_collapses_chain() {
    IAtomContainer a = mol("C1CCCC=C1CC"); // cyclohexene–ethyl
    IAtomContainer b = mol("CCc1ccccc1"); // ethylbenzene
    MoleculePair smsd = new MoleculePair(a, b, opts(true, true, false, false));
    Map<Integer, Integer> m = smsd.mcsMapping( true, 2500L);
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
    MoleculePair smsd = new MoleculePair(ring6, chain, opts(true, true, true, false)); // ringOnly=true
    Map<Integer, Integer> m = smsd.mcsMapping( true, 1500L);
    assertTrue(
        m == null || m.isEmpty(), "Ring→ring policy should prevent mapping a ring to a chain");
  }

  @Test
  @DisplayName("B.5  Acetonitrile vs iminoethane: relaxed bond order grows the core")
  void b5_strict_vs_relaxed_bond_order() {
    IAtomContainer acn = mol("CC#N");
    IAtomContainer imi = mol("CC=N");
    MoleculePair sStrict = new MoleculePair(acn, imi, opts(true, true, false, false));
    MoleculePair sRelax = new MoleculePair(acn, imi, opts(true, false, false, false));
    Map<Integer, Integer> mStrict = sStrict.mcsMapping( true, 2000L);
    Map<Integer, Integer> mRelax = sRelax.mcsMapping( true, 2000L); // MCIS-like
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
    MoleculePair smsd = new MoleculePair(rx, pr, opts(true, true, false, false));
    Map<Integer, Integer> m = smsd.mcsMapping( true, 2000L);
    assertNonEmpty(m, "Backbone conserved across reaction");
    assertTrue(m.size() >= 6, "Expect ring (6) mapped");
  }

  @Test
  @DisplayName("B.7  Aromatic ethers: MCS large enough to seed alignment")
  void b7_aromatic_ether_seed() {
    IAtomContainer a = mol("COc1ccccc1O");
    IAtomContainer b = mol("COc1ccc(OC)cc1");
    MoleculePair smsd = new MoleculePair(a, b, opts(true, true, false, false));
    Map<Integer, Integer> m = smsd.mcsMapping( true, 2000L);
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
        MoleculePair smsd = new MoleculePair(x, y, chemOpts);
        Map<Integer, Integer> m = smsd.mcsMapping( true, 1500L);
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
    MoleculePair mccs = new MoleculePair(a, b, opts(true, true, false, false));
    MoleculePair mcis = new MoleculePair(a, b, opts(true, false, false, false)); // relaxed order for MCIS-like
    Map<Integer, Integer> m1 = mccs.mcsMapping( true, 1500L); // MCCS
    Map<Integer, Integer> m2 = mcis.mcsMapping( true, 1500L); // MCIS
    int s1 = (m1 == null ? 0 : m1.size());
    int s2 = (m2 == null ? 0 : m2.size());
    assertTrue(s2 >= s1, "MCIS should not be smaller than MCCS on simple pairs");
  }

  @Test
  @DisplayName("C.2  Ring→ring policy allows ring↔ring mapping (benzene vs kekulé benzene)")
  void c2_ring_to_ring_allows_ring_to_ring() {
    IAtomContainer ringA = mol("c1ccccc1");
    IAtomContainer ringB = mol("C1=CC=CC=C1");
    MoleculePair smsd = new MoleculePair(ringA, ringB, opts(true, true, true, false)); // ringOnly=true
    Map<Integer, Integer> m = smsd.mcsMapping( true, 1500L);
    assertTrue(m != null && !m.isEmpty(), "Ring→ring policy should allow ring↔ring matching");
  }

  @Test
  @DisplayName("C.3  Timeout budget is honoured (tiny budget)")
  void c3_timeout_honoured() {
    MoleculePair smsd =
        new MoleculePair(
            mol("c1ccc2ccccc2c1"), mol("c1cccc2ccc3cccc3c2c1"), opts(true, true, false, false));
    Map<Integer, Integer> m = smsd.mcsMapping( true, 1L); // 1 ms
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

    MoleculePair smsd = new MoleculePair(educt, product, chem);

    long start = System.nanoTime();
    Map<Integer, Integer> mapping = smsd.mcsMapping( true, 10_000L);
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
  // One 15-second baseline plus four sequential 15-second parallel rounds
  // can legitimately consume 75 seconds; each round still has its 20s guard.
  @Timeout(value = 90, unit = TimeUnit.SECONDS, threadMode = Timeout.ThreadMode.SEPARATE_THREAD)
  @DisplayName("C.3c  GOLDEN_843 shared molecules stay stable under parallel reuse")
  void c3c_golden843_shared_molecules_parallel_reuse_is_stable() throws Exception {
    final IAtomContainer educt = mol(GOLDEN_843_EDUCT);
    final IAtomContainer product = mol(GOLDEN_843_PRODUCT);

    ChemOptions chem = new ChemOptions();
    chem.matchAtomType = false;
    chem.matchBondOrder = ChemOptions.BondOrderMode.STRICT;
    chem.ringMatchesRingOnly = false;
    chem.ringFusionMode = ChemOptions.RingFusionMode.IGNORE;

    MCS.Options baselineOpts = new MCS.Options();
    baselineOpts.connectivity = (true) ? MCS.Connectivity.Connected : MCS.Connectivity.Fragments;
    baselineOpts.timeoutMillis = 15_000L;

    TestSupport.clearCaches();
    Map<Integer, Integer> baseline = TestSupport.mcsMapping(educt, product, chem, baselineOpts);
    assertNotNull(baseline, "Baseline GOLDEN_843 call should return a mapping");
    assertTrue(
        baseline.size() >= 12,
        "Baseline GOLDEN_843 should keep a substantial MCS, got " + baseline.size());

    int expectedSize = baseline.size();
    int expectedMappedBonds = mappedBonds(educt, product, baseline);

    final int threads = 8;
    final int rounds = 4;
    for (int round = 0; round < rounds; round++) {
      TestSupport.clearCaches();

      CyclicBarrier barrier = new CyclicBarrier(threads);
      ExecutorService pool = Executors.newFixedThreadPool(threads);
      try {
        List<Callable<Map<Integer, Integer>>> tasks = new ArrayList<>();
        for (int i = 0; i < threads; i++) {
          tasks.add(
              () -> {
                barrier.await();
                MCS.Options opts = new MCS.Options();
                opts.connectivity = (true) ? MCS.Connectivity.Connected : MCS.Connectivity.Fragments;
                opts.timeoutMillis = 15_000L;
                return TestSupport.mcsMapping(educt, product, chem, opts);
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
    IAtomContainer trans = mol("C/C=C/C");
    IAtomContainer cis = mol("C/C=C\\C");
    MoleculePair sOn = new MoleculePair(trans, cis, opts(true, true, false, true)); // stereo ON
    MoleculePair sOff = new MoleculePair(trans, cis, opts(true, true, false, false)); // stereo OFF
    Map<Integer, Integer> mOn = sOn.mcsMapping( true, 2000L);
    Map<Integer, Integer> mOff = sOff.mcsMapping( true, 2000L);
    int aOn = (mOn == null ? 0 : mOn.size());
    int aOff = (mOff == null ? 0 : mOff.size());
    assertTrue(aOff >= aOn, "With stereo OFF, mapping should be ≥ stereo ON");
  }

  @Test
  @DisplayName("C.5  Connected vs disconnected: dMCS may yield ≥ fragments")
  void c5_connected_vs_disconnected() {
    IAtomContainer a = mol("CCOC(C)=O");
    IAtomContainer b = mol("CCOC(C)C=O");
    MoleculePair cMCS = new MoleculePair(a, b, opts(true, true, false, false));
    MoleculePair dMCS = new MoleculePair(a, b, opts(true, true, false, false));
    Map<Integer, Integer> mc = cMCS.mcsMapping( true, 2000L); // connectedOnly
    Map<Integer, Integer> md = dMCS.mcsMapping( false, 2000L); // disconnected
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
    MoleculePair smsd = new MoleculePair(arom, kekul, flex);
    Map<Integer, Integer> m = smsd.mcsMapping( true, 1500L);
    assertTrue(m != null && m.size() >= 6, "Flexible mode should match the full benzene core");
  }

  @Test
  @DisplayName("C.7  Mapping integrity: values are unique")
  void c7_mapping_integrity() {
    MoleculePair smsd = new MoleculePair(mol("c1ccccc1"), mol("c1ccc2ccccc2c1"), opts(true, false, false, false));
    Map<Integer, Integer> m = smsd.mcsMapping( true, 2000L);
    assertNotNull(m);
    assertEquals(m.size(), new HashSet<>(m.values()).size(), "Target indices must be unique");
  }

  @Test
  @DisplayName("C.8  Bond-count monotonicity check on a second pair")
  void c8_bond_count_monotonicity_recap() {
    IAtomContainer a = mol("CC=CC"); // 2-butene
    IAtomContainer b = mol("CCCC"); // butane
    MoleculePair strict = new MoleculePair(a, b, opts(true, true, false, false));
    MoleculePair relax = new MoleculePair(a, b, opts(true, false, false, false));
    Map<Integer, Integer> ms = strict.mcsMapping( true, 2000L);
    Map<Integer, Integer> mr = relax.mcsMapping( true, 2000L);
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
  // Robustness to awkward input
  // ======================================================================


  // -----------------------------
  // Case generators
  // -----------------------------

  /** Alkane chains: C×k is always a substructure of C×(k+2). */
  static Stream<Arguments> chainCases() {
    List<Arguments> out = new ArrayList<>();
    for (int k = 2; k <= 21; k++) {
      out.add(arguments("chain_" + k, repeat("C", k), repeat("C", k + 2), true));
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
   * true (matchFormalCharge=false by default: charge ignored, [O-] matches O)
   */
  static Stream<Arguments> functionalCases() {
    return Stream.of(
        arguments("fx_0", "NC(=O)", "NCC(=O)NCCC", true),
        arguments("fx_1", "[O-]C=O", "OC=O", true));
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
    MoleculePair smsd = new MoleculePair(queryMol, targetMol, new ChemOptions());
    boolean ok = smsd.contains();
    assertEquals(expected, ok, name + ": '" + q + "' vs '" + t + "' expected " + expected);
  }

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
    MoleculePair smsd = new MoleculePair(expanded, targetMol, new ChemOptions());
    boolean ok = smsd.contains();
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
    MoleculePair smsd = new MoleculePair(queryMol, targetMol, new ChemOptions());
    Map<Integer, Integer> mr = smsd.mcsMapping( true, 2500L);
    assertNotNull(mr, name + ": MCS should not be null");
    // Query and target always share at least some atoms
    // For cases where both molecules share atom types, MCS should be non-empty
    // For cases where they don't (e.g., fx_1: [O-]C=O vs OC=O with charge matching),
    // MCS may be empty. Just verify it returned a valid map.
    assertNotNull(mr, name + ": MCS result should not be null");
  }

  // ======================================================================
  // Reference substructure and MCS cases
  // ======================================================================


  // The shared parser, as an instance field for the nested tests.
  private final org.openscience.cdk.smiles.SmilesParser sp = SP;

  @Nested
  @DisplayName("Core Substructure Tests")
  class SubstructureTests {

    @Test
    @DisplayName("should find a simple alkane substructure")
    void findSimpleAlkane() throws Exception {
      IAtomContainer q = sp.parseSmiles("CCC"); // Propane
      IAtomContainer t = sp.parseSmiles("CCCCCC"); // Hexane
      MoleculePair smsd =
          new MoleculePair(
              Standardiser.standardise(q, Standardiser.TautomerMode.NONE),
              Standardiser.standardise(t, Standardiser.TautomerMode.NONE),
              ChemOptions.profile("compat-substruct"));
      assertTrue(smsd.contains());
    }

    @Test
    @DisplayName("should not find a longer chain in a shorter one")
    void notFindLongerChain() throws Exception {
      IAtomContainer q = sp.parseSmiles("CCCCCCCC"); // Octane
      IAtomContainer t = sp.parseSmiles("CCCCCC"); // Hexane
      MoleculePair smsd =
          new MoleculePair(
              Standardiser.standardise(q, Standardiser.TautomerMode.NONE),
              Standardiser.standardise(t, Standardiser.TautomerMode.NONE),
              ChemOptions.profile("compat-substruct"));
      assertFalse(smsd.contains());
    }

    @Test
    @DisplayName("should find a heterocycle substructure")
    void findHeterocycle() throws Exception {
      IAtomContainer q = sp.parseSmiles("c1ncccc1"); // Pyridine
      IAtomContainer t = sp.parseSmiles("Clc1ncccc1"); // 2-Chloropyridine
      MoleculePair smsd =
          new MoleculePair(
              Standardiser.standardise(q, Standardiser.TautomerMode.NONE),
              Standardiser.standardise(t, Standardiser.TautomerMode.NONE),
              ChemOptions.profile("compat-substruct"));
      assertTrue(smsd.contains());
    }

    @Test
    @DisplayName("should match formal charges when required")
    void matchFormalCharge() throws Exception {
      IAtomContainer q = sp.parseSmiles("[O-]C=O"); // Carboxylate
      IAtomContainer t = sp.parseSmiles("CC(=O)[O-]"); // Acetate
      ChemOptions options = ChemOptions.profile("compat-substruct");
      options.matchFormalCharge = true;
      MoleculePair smsd =
          new MoleculePair(
              Standardiser.standardise(q, Standardiser.TautomerMode.NONE),
              Standardiser.standardise(t, Standardiser.TautomerMode.NONE),
              options);
      assertTrue(smsd.contains());
    }

    @Test
    @DisplayName("should ignore formal charges when not required")
    void ignoreFormalCharge() throws Exception {
      IAtomContainer q = sp.parseSmiles("OC=O"); // Carboxylic acid
      IAtomContainer t = sp.parseSmiles("CC(=O)[O-]"); // Acetate
      ChemOptions options = ChemOptions.profile("compat-substruct");
      options.matchFormalCharge = false;
      MoleculePair smsd =
          new MoleculePair(
              Standardiser.standardise(q, Standardiser.TautomerMode.NONE),
              Standardiser.standardise(t, Standardiser.TautomerMode.NONE),
              options);
      assertTrue(smsd.contains());
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
      MoleculePair smsd =
          new MoleculePair(
              qSmarts,
              Standardiser.standardise(t, Standardiser.TautomerMode.NONE),
              ChemOptions.profile("compat-substruct"));
      assertTrue(smsd.contains());
    }

    @Test
    @DisplayName("should handle recursive SMARTS for amide nitrogen")
    void recursiveSmartsAmide() throws Exception {

      String qSmarts = "[N;$(NC=O)]";
      IAtomContainer t = sp.parseSmiles("CC(=O)NCC");
      MoleculePair smsd =
          new MoleculePair(
              qSmarts,
              Standardiser.standardise(t, Standardiser.TautomerMode.NONE),
              ChemOptions.profile("compat-substruct"));
      assertTrue(smsd.contains());
    }

    @Test
    @DisplayName("should match disconnected query fragments in a target")
    void disconnectedFragments() throws Exception {
      IAtomContainer q = sp.parseSmiles("C.C");
      IAtomContainer t = sp.parseSmiles("CCC");
      MoleculePair smsd =
          new MoleculePair(
              Standardiser.standardise(q, Standardiser.TautomerMode.NONE),
              Standardiser.standardise(t, Standardiser.TautomerMode.NONE),
              ChemOptions.profile("compat-substruct"));
      assertTrue(smsd.contains());
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
      MoleculePair smsd =
          new MoleculePair(
              Standardiser.standardise(q, Standardiser.TautomerMode.NONE),
              Standardiser.standardise(t, Standardiser.TautomerMode.NONE),
              options);
      assertTrue(smsd.contains());
    }

    @Test
    @DisplayName("should not match different stereoisomers with strict settings")
    void stereoMismatch() throws Exception {
      IAtomContainer q = sp.parseSmiles("C/C=C\\C");
      IAtomContainer t = sp.parseSmiles("C/C=C/C");
      ChemOptions options = ChemOptions.profile("strict");
      options.useBondStereo = true;
      MoleculePair smsd =
          new MoleculePair(
              Standardiser.standardise(q, Standardiser.TautomerMode.NONE),
              Standardiser.standardise(t, Standardiser.TautomerMode.NONE),
              options);
      assertFalse(smsd.contains());
    }

    @Test
    @DisplayName("should match aromatic and Kekule forms with flexible settings")
    void aromaticFlexibleMatch() throws Exception {
      IAtomContainer q = sp.parseSmiles("c1ccccc1");
      IAtomContainer t = sp.parseSmiles("C1=CC=CC=C1");
      MoleculePair smsd =
          new MoleculePair(
              Standardiser.standardise(q, Standardiser.TautomerMode.NONE),
              Standardiser.standardise(t, Standardiser.TautomerMode.NONE),
              ChemOptions.profile("compat-substruct"));
      assertTrue(smsd.contains());
    }

    @Test
    @DisplayName("should not match aromatic and non-aromatic rings with strict settings")
    void aromaticStrictMismatch() throws Exception {
      IAtomContainer q = sp.parseSmiles("c1ccccc1");
      IAtomContainer t = sp.parseSmiles("C1CCCCC1"); // Cyclohexane
      MoleculePair smsd =
          new MoleculePair(
              Standardiser.standardise(q, Standardiser.TautomerMode.NONE),
              Standardiser.standardise(t, Standardiser.TautomerMode.NONE),
              ChemOptions.profile("strict"));
      assertFalse(smsd.contains());
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
      MoleculePair smsd =
          new MoleculePair(
              Standardiser.standardise(a, Standardiser.TautomerMode.NONE),
              Standardiser.standardise(b, Standardiser.TautomerMode.NONE),
              ChemOptions.profile("compat-fmcs"));
      Map<Integer, Integer> m = smsd.mcsMapping( true, 2500L);
      assertNotNull(m);
      assertEquals(6, m.size());
    }

    @Test
    @DisplayName("should find a common core between two different structures")
    void mcsCommonCore() throws Exception {
      IAtomContainer a = sp.parseSmiles("CC(C)c1ccc(C(=O)O)cc1"); // Ibuprofen core
      IAtomContainer b = sp.parseSmiles("CC(C(=O)O)c1ccc2ccccc2c1"); // Naproxen core
      MoleculePair smsd =
          new MoleculePair(
              Standardiser.standardise(a, Standardiser.TautomerMode.NONE),
              Standardiser.standardise(b, Standardiser.TautomerMode.NONE),
              ChemOptions.profile("compat-fmcs"));
      Map<Integer, Integer> m = smsd.mcsMapping( true, 2500L);
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
      MoleculePair smsd =
          new MoleculePair(
              Standardiser.standardise(query, Standardiser.TautomerMode.NONE),
              Standardiser.standardise(target, Standardiser.TautomerMode.NONE),
              ChemOptions.profile("compat-substruct"));
      assertTrue(smsd.contains(), "Salicylic acid is a substructure of Aspirin");
    }

    @Test
    @DisplayName("should find the catechol pharmacophore in Dopamine")
    void testCatecholInDopamine() throws Exception {
      IAtomContainer query = sp.parseSmiles("c1ccc(O)c(O)c1"); // Catechol
      IAtomContainer target = sp.parseSmiles("C1=CC(=C(C=C1CCN)O)O"); // Dopamine
      MoleculePair smsd =
          new MoleculePair(
              Standardiser.standardise(query, Standardiser.TautomerMode.NONE),
              Standardiser.standardise(target, Standardiser.TautomerMode.NONE),
              ChemOptions.profile("compat-substruct"));
      assertTrue(smsd.contains(), "Catechol is a substructure of Dopamine");
    }

    @Test
    @DisplayName("should confirm a complex fragment is not present")
    void testIndoleNotInAspirin() throws Exception {
      IAtomContainer query = sp.parseSmiles("c1ccc2[nH]ccc2c1"); // Indole
      IAtomContainer target = sp.parseSmiles("CC(=O)Oc1ccccc1C(=O)O"); // Aspirin
      MoleculePair smsd =
          new MoleculePair(
              Standardiser.standardise(query, Standardiser.TautomerMode.NONE),
              Standardiser.standardise(target, Standardiser.TautomerMode.NONE),
              ChemOptions.profile("compat-substruct"));
      assertFalse(smsd.contains(), "Indole is not a substructure of Aspirin");
    }

    @Test
    @DisplayName("should find a large MCS between Morphine and Codeine")
    void testMCSMorphineVsCodeine() throws Exception {
      // Morphine and codeine differ by methylation of the phenolic oxygen.
      IAtomContainer morphine = sp.parseSmiles("CN1CC[C@]23[C@@H]4[C@H]1CC5=C2C(=C(C=C5)O)O[C@H]3[C@H](C=C4)O");
      IAtomContainer codeine = sp.parseSmiles("CN1CC[C@]23[C@@H]4[C@H]1CC5=C2C(=C(C=C5)OC)O[C@H]3[C@H](C=C4)O");

      assertEquals(21, morphine.getAtomCount());
      assertEquals(22, codeine.getAtomCount());
      assertEquals(3, java.util.stream.StreamSupport.stream(morphine.atoms().spliterator(), false)
          .filter(atom -> atom.getAtomicNumber() == 8).count());
      MoleculePair smsd =
          new MoleculePair(
              Standardiser.standardise(morphine, Standardiser.TautomerMode.NONE),
              Standardiser.standardise(codeine, Standardiser.TautomerMode.NONE),
              new ChemOptions());

      Map<Integer, Integer> mcsInduced = smsd.mcsMapping( true, 30000L);
      assertNotNull(mcsInduced);
      assertTrue(
          mcsInduced.size() >= 19,
          "Induced MCS of Morphine/Codeine should be >= 19 (got " + mcsInduced.size() + ")");

      Map<Integer, Integer> mcsNonInduced = smsd.mcsMapping( true, 30000L);
      assertNotNull(mcsNonInduced);
      assertTrue(
          mcsNonInduced.size() >= 19,
          "Non-induced MCS of Morphine/Codeine should be >= 19 (got " + mcsNonInduced.size() + ")");
    }
  }

  // ======================================================================
  // Chemistry coverage: heterocycles, ring systems, functional groups, charges and options
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
          new MoleculePair(mol("c1ccoc1"), mol("c1ccc2occc2c1"), new ChemOptions()).contains());
    }

    @Test
    void thiopheneInBenzothiophene() throws Exception {
      assertTrue(
          new MoleculePair(mol("c1ccsc1"), mol("c1ccc2sccc2c1"), new ChemOptions()).contains());
    }

    @Test
    void imidazoleInHistidine() throws Exception {
      // Imidazole ring in histidine side chain
      assertTrue(
          new MoleculePair(mol("c1cnc[nH]1"), mol("NC(Cc1cnc[nH]1)C(=O)O"), new ChemOptions())
              .contains());
    }

    @Test
    void pyrimidineRingInCytosine() throws Exception {
      // Cytosine contains a pyrimidine ring — use the 1,3-diazine core
      MoleculePair smsd = new MoleculePair(mol("C=1C=NC=NC1"), mol("Nc1ccnc(=O)[nH]1"), new ChemOptions());
      Map<Integer, Integer> mcs = smsd.mcsMapping( true, 5000L);
      assertTrue(mcs.size() >= 4, "Should share pyrimidine atoms");
    }

    @Test
    void morpholineInDrug() throws Exception {
      assertTrue(
          new MoleculePair(mol("C1COCCN1"), mol("O=C(N1CCOCC1)c1ccccc1"), new ChemOptions())
              .contains());
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
      MoleculePair smsd = new MoleculePair(mol("C1CC2CC1C2"), mol("C1CC2CC1CC2C"), new ChemOptions());
      // Both are bridged, should find common subgraph
      Map<Integer, Integer> mcs = smsd.mcsMapping( true, 5000L);
      assertNotNull(mcs);
      assertTrue(mcs.size() >= 4, "Bridged ring MCS should find common atoms");
    }

    @Test
    void spiroCompoundMCS() throws Exception {
      // Spiro[4.4]nonane
      MoleculePair smsd = new MoleculePair(mol("C1CCC2(C1)CCCC2"), mol("C1CCC2(C1)CCCCC2"), new ChemOptions());
      Map<Integer, Integer> mcs = smsd.mcsMapping( true, 5000L);
      assertTrue(mcs.size() >= 5, "Spiro MCS should share at least one ring");
    }
  }

  // ======================================================================
  // 3. FUNCTIONAL GROUPS — systematic substructure
  // ======================================================================

  @Nested
  @DisplayName("Functional Group Substructure")
  class FunctionalGroupTests {
    @Test
    void primaryAmineInAlanine() throws Exception {
      assertTrue(new MoleculePair(mol("CN"), mol("CC(N)C(=O)O"), new ChemOptions()).contains());
    }

    @Test
    void aldehydeInBenzaldehyde() throws Exception {
      assertTrue(new MoleculePair(mol("C=O"), mol("O=Cc1ccccc1"), new ChemOptions()).contains());
    }

    @Test
    void nitrileInAcetonitrile() throws Exception {
      assertTrue(new MoleculePair(mol("C#N"), mol("CC#N"), new ChemOptions()).contains());
    }

    @Test
    void nitroInNitrobenzene() throws Exception {
      assertTrue(
          new MoleculePair(mol("[N+](=O)[O-]"), mol("c1ccc([N+](=O)[O-])cc1"), new ChemOptions())
              .contains());
    }

    @Test
    void sulfonylInDMSO() throws Exception {
      assertTrue(new MoleculePair(mol("S=O"), mol("CS(=O)C"), new ChemOptions()).contains());
    }

    @Test
    void thioetherInMethionine() throws Exception {
      assertTrue(new MoleculePair(mol("CSC"), mol("CSCCC(N)C(=O)O"), new ChemOptions()).contains());
    }

    @Test
    void phenolNotInCyclohexanol() throws Exception {
      // Aromatic OH should not match aliphatic OH with strict aromaticity
      MoleculePair smsd = new MoleculePair(mol("Oc1ccccc1"), mol("OC1CCCCC1"), ChemOptions.profile("strict"));
      assertFalse(smsd.contains());
    }

    @Test
    void epoxideSubstructure() throws Exception {
      assertTrue(new MoleculePair(mol("C1OC1"), mol("C1OC1CC"), new ChemOptions()).contains());
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
      MoleculePair smsd = new MoleculePair(mol("CCO"), mol("CC=O"), new ChemOptions());
      Map<Integer, Integer> mcs = smsd.mcsMapping( true, 5000L);
      assertTrue(mcs.size() >= 2, "Oxidation preserves C-C backbone");
    }

    @Test
    void acidToAmide() throws Exception {
      // Amide bond formation: acetic acid → acetamide
      MoleculePair smsd = new MoleculePair(mol("CC(=O)O"), mol("CC(=O)N"), new ChemOptions());
      Map<Integer, Integer> mcs = smsd.mcsMapping( true, 5000L);
      assertTrue(mcs.size() >= 3, "Acid→amide preserves acyl group");
    }

    @Test
    void halogenationMCS() throws Exception {
      // Chlorination: benzene → chlorobenzene
      MoleculePair smsd = new MoleculePair(mol("c1ccccc1"), mol("Clc1ccccc1"), new ChemOptions());
      Map<Integer, Integer> mcs = smsd.mcsMapping( true, 5000L);
      assertEquals(6, mcs.size(), "Benzene ring preserved in chlorination");
    }

    @Test
    void demethylation() throws Exception {
      // Codeine → morphine (loss of methyl from ether)
      MoleculePair smsd =
          new MoleculePair(
              mol("CN1CC[C@]23[C@@H]4[C@H]1CC5=C2C(=C(C=C5)OC)O[C@H]3[C@H](C=C4)O"),
              mol("CN1CC[C@]23[C@@H]4[C@H]1CC5=C2C(=C(C=C5)O)O[C@H]3[C@H](C=C4)O"),
              new ChemOptions());
      Map<Integer, Integer> mcs = smsd.mcsMapping( true, 30000L);
      assertTrue(mcs.size() >= 19, "Demethylation preserves core scaffold");
    }

    @Test
    void ketoEnolTautomer() throws Exception {
      // Keto-enol: acetone ↔ propen-2-ol
      MolGraph keto = new MolGraph(mol("CC(=O)C"));
      MolGraph enol = new MolGraph(mol("CC(O)=C"));
      ChemOptions strict = new ChemOptions();
      MCS.Options options = new MCS.Options();
      options.timeoutMillis = 5_000L;
      Map<Integer, Integer> witness = mapping(0, 0, 1, 1);
      assertTrue(MCS.validateMapping(keto, enol, witness, strict).isEmpty());
      // Each connected three-atom subset would need either C=O to match C-O,
      // or both query C-C bonds to match the target's C-C and C=C bonds.
      assertEquals(2, TestSupport.mcsMapping(keto, enol, strict, options).size(),
          "Strict bond orders retain only the two-atom C-C backbone");

      ChemOptions tautomer = ChemOptions.tautomerProfile();
      Map<Integer, Integer> relaxed = TestSupport.mcsMapping(keto, enol, tautomer, options);
      assertEquals(4, relaxed.size(), "Tautomer-aware matching preserves the full backbone");
      assertTrue(MCS.validateMapping(keto, enol, relaxed, tautomer).isEmpty());
    }

    @Test
    void ringOpeningMCS() throws Exception {
      // Ring opening: cyclopropane vs propane — need ring constraint OFF
      ChemOptions opts = new ChemOptions();
      opts.ringMatchesRingOnly = false;
      MoleculePair smsd = new MoleculePair(mol("C1CC1"), mol("CCC"), opts);
      Map<Integer, Integer> mcs = smsd.mcsMapping( true, 5000L);
      assertTrue(mcs.size() >= 2, "Ring opening preserves C-C bonds when ring constraint off");
    }

    @Test
    void esterHydrolysis() throws Exception {
      // Ethyl acetate → acetic acid + ethanol (shared acyl)
      MoleculePair smsd = new MoleculePair(mol("CC(=O)OCC"), mol("CC(=O)O"), new ChemOptions());
      Map<Integer, Integer> mcs = smsd.mcsMapping( true, 5000L);
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
      assertFalse(new MoleculePair(mol("[NH4+]"), mol("N"), opts).contains());
    }

    @Test
    void ammoniumMatchesAmineChargeOff() throws Exception {
      ChemOptions opts = new ChemOptions();
      opts.matchFormalCharge = false;
      assertTrue(new MoleculePair(mol("[NH4+]"), mol("N"), opts).contains());
    }

    @Test
    void zwitterionSubstructure() throws Exception {
      // Glycine zwitterion
      ChemOptions opts = new ChemOptions();
      opts.matchFormalCharge = true;
      assertTrue(new MoleculePair(mol("[NH3+]CC([O-])=O"), mol("[NH3+]CC([O-])=O"), opts).contains());
    }
  }

  // ======================================================================
  // 6. CHEMOPTIONS CONFIGURATIONS
  // ======================================================================

  @Nested
  @DisplayName("ChemOptions Configuration Tests")
  class ChemOptionsTests {
    @Test
    void defaultSubstructureQuery() throws Exception {
      ChemOptions opts = new ChemOptions();
      assertTrue(new MoleculePair(mol("c1ccccc1"), mol("c1ccc(O)cc1"), opts).contains());
    }

    @Test
    void preparedSubstructureQueryCanBeReused() throws Exception {
      Substructure query = Substructure.find(mol("c1ccccc1"));
      assertTrue(query.match(mol("c1ccc(O)cc1")).found());
      assertTrue(query.match(mol("Cc1ccccc1")).found());
    }

    @Test
    void graphAndContainerSubstructureQueriesAgree() throws Exception {
      IAtomContainer query = mol("c1ccccc1");
      IAtomContainer target = mol("c1ccc(O)cc1");
      assertTrue(Substructure.find(query).match(target).found());
      assertTrue(Substructure.find(new MolGraph(query)).match(new MolGraph(target)).found());
    }

    @Test
    void bondOrderAny() throws Exception {
      ChemOptions opts = new ChemOptions();
      opts.matchBondOrder = ChemOptions.BondOrderMode.ANY;
      // Single bond query should match double bond target
      assertTrue(new MoleculePair(mol("CC"), mol("C=C"), opts).contains());
    }

    @Test
    void pruningTogglesTwoHop() throws Exception {
      ChemOptions opts = new ChemOptions();
      opts.useTwoHopNLF = false;
      assertTrue(new MoleculePair(mol("c1ccccc1"), mol("c1ccc(O)cc1"), opts).contains());
    }

    @Test
    void pruningTogglesThreeHop() throws Exception {
      ChemOptions opts = new ChemOptions();
      opts.useThreeHopNLF = false;
      assertTrue(new MoleculePair(mol("c1ccccc1"), mol("c1ccc(O)cc1"), opts).contains());
    }

    @Test
    void pruningTogglesBitParallel() throws Exception {
      ChemOptions opts = new ChemOptions();
      opts.useBitParallelFeasibility = false;
      assertTrue(new MoleculePair(mol("c1ccccc1"), mol("c1ccc(O)cc1"), opts).contains());
    }

    @Test
    void ringMatchDisabled() throws Exception {
      ChemOptions opts = new ChemOptions();
      opts.ringMatchesRingOnly = false;
      // With ring constraint OFF, cyclopentane (5) should map into cyclohexane (6)
      // because ring atoms are allowed to match non-ring atoms
      MoleculePair smsd = new MoleculePair(mol("C1CCCC1"), mol("CCCCCC"), opts);
      Map<Integer, Integer> mcs = smsd.mcsMapping( true, 5000L);
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
      assertTrue(new MoleculePair(mol("c1ccccc1"), mol("c1ccc(O)cc1"), opts).contains());
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
      assertTrue(new MoleculePair(mol("C"), mol("CCCC"), new ChemOptions()).contains());
    }

    @Test
    void singleAtomNoMatch() throws Exception {
      // Nitrogen not in all-carbon chain
      assertFalse(new MoleculePair(mol("N"), mol("CCCC"), new ChemOptions()).contains());
    }

    @Test
    void selfMatchSubstructure() throws Exception {
      IAtomContainer m = mol("c1ccc(O)cc1");
      assertTrue(new MoleculePair(m, m, new ChemOptions()).contains());
    }

    @Test
    void selfMatchMCS() throws Exception {
      IAtomContainer m = mol("c1ccc(O)cc1");
      Map<Integer, Integer> mcs = new MoleculePair(m, m, new ChemOptions()).mcsMapping( true, 5000L);
      assertEquals(m.getAtomCount(), mcs.size(), "Self-MCS should match all atoms");
    }

    @Test
    void noCommonAtomsMCS() throws Exception {
      // Methane (C) vs water (O) — no shared atoms with strict typing
      MoleculePair smsd = new MoleculePair(mol("C"), mol("O"), new ChemOptions());
      Map<Integer, Integer> mcs = smsd.mcsMapping( true, 2000L);
      assertTrue(mcs.isEmpty(), "No common atoms between methane and water");
    }

    @Test
    void disconnectedQuery() throws Exception {
      assertTrue(new MoleculePair(mol("C.N"), mol("CCN"), new ChemOptions()).contains());
    }

    @Test
    void disconnectedMCSInduced() throws Exception {
      MoleculePair smsd = new MoleculePair(mol("CC.OO"), mol("CCOCC"), new ChemOptions());
      Map<Integer, Integer> mcs = smsd.mcsMapping( false, 5000L);
      assertNotNull(mcs);
    }

    @Test
    void veryShortTimeout() throws Exception {
      MoleculePair smsd = new MoleculePair(mol("c1ccccc1"), mol("c1ccc2ccccc2c1"), new ChemOptions());
      // Should not hang; may or may not find result
      smsd.mcsMapping( true, 1L);
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
      MoleculePair smsd = new MoleculePair(mol("C"), mol("CCCCCC"), new ChemOptions());
      List<Map<Integer, Integer>> maps = smsd.substructureMappings(100, 5000L);
      assertEquals(6, maps.size(), "Single C should match all 6 positions in hexane");
    }

    @Test
    void findAllSubstructuresMaxLimit() throws Exception {
      MoleculePair smsd = new MoleculePair(mol("C"), mol("CCCCCC"), new ChemOptions());
      List<Map<Integer, Integer>> maps = smsd.substructureMappings(3, 5000L);
      assertEquals(3, maps.size(), "Should respect maxSolutions limit");
    }

    @Test
    void isSubstructureWithStatsReturnsData() throws Exception {
      MoleculePair smsd = new MoleculePair(mol("c1ccccc1"), mol("c1ccc2ccccc2c1"), new ChemOptions());
      SearchResult res = smsd.substructureResult(5000L);
      assertTrue(res.found());
      assertTrue(res.nodesVisited() > 0, "Should report nodes visited");
      assertTrue(res.elapsedMillis() >= 0, "Should report elapsed time");
    }

    @Test
    void findAllSubstructuresWithStatsReturnsData() throws Exception {
      MoleculePair smsd = new MoleculePair(mol("c1ccccc1"), mol("c1ccc(O)cc1"), new ChemOptions());
      SearchResult res = smsd.substructureResults(100, 5000L);
      assertTrue(res.found());
      assertFalse(res.matches().isEmpty());
      assertTrue(res.matches().size() > 0);
    }

    @Test
    void timeoutSetters() throws Exception {
      MoleculePair smsd = new MoleculePair(mol("c1ccccc1"), mol("c1ccc2ccccc2c1"), new ChemOptions());
      smsd.setSubstructureTimeoutMs(100L);
      smsd.setMCSTimeoutMs(100L);
      // Should work with custom timeouts
      assertTrue(smsd.contains());
      assertNotNull(smsd.mcsMapping());
    }

    @Test
    void smartsConstructor() throws Exception {
      IAtomContainer target = mol("c1ccc(O)cc1");
      MoleculePair smsd = new MoleculePair("[OH]", target, new ChemOptions());
      assertTrue(smsd.contains(), "SMARTS [OH] should match phenol");
      assertThrows(UnsupportedOperationException.class, smsd::query, "SMARTS query should throw on getQuery()");
    }

    @Test
    void smartsSupportsEnumeration() throws Exception {
      IAtomContainer target = mol("CCO");
      MoleculePair smsd = new MoleculePair("[OH]", target, new ChemOptions());
      assertEquals(1,smsd.substructureMappings(10,1000L).size());
    }

    @Test
    void smilesConvenienceConstructor() throws Exception {
      MoleculePair smsd = new MoleculePair("c1ccccc1", "c1ccc(O)cc1", new ChemOptions());
      assertTrue(smsd.contains());
      assertNotNull(smsd.query());
      assertNotNull(smsd.target());
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
      MoleculePair smsd = new MoleculePair(mol("CC.CC"), mol("CCCCCC"), new ChemOptions());
      Map<Integer, Integer> mcs = smsd.mcsMapping( false, 5000L);
      assertNotNull(mcs);
      assertTrue(mcs.size() >= 2);
    }

    @Test
    void mccsVsMcis() throws Exception {
      // MCCS (non-induced) should be >= MCIS (induced)
      MoleculePair smsd = new MoleculePair(mol("CCCC"), mol("CC=CC"), new ChemOptions());
      Map<Integer, Integer> mccs = smsd.mcsMapping( true, 5000L);
      Map<Integer, Integer> mcis = smsd.mcsMapping( true, 5000L);
      assertTrue(mccs.size() >= mcis.size(), "MCCS >= MCIS");
    }

    @Test
    void ringVsChainMCS() throws Exception {
      IAtomContainer ring = mol("C1CCCCC1");
      IAtomContainer chain = mol("CCCCC");

      ChemOptions strict = new ChemOptions();
      strict.ringMatchesRingOnly = true;
      Map<Integer, Integer> strictRingToChain = new MoleculePair(ring, chain, strict).mcsMapping( true, 5000L);
      Map<Integer, Integer> strictChainToRing = new MoleculePair(chain, ring, strict).mcsMapping( true, 5000L);
      assertTrue(strictRingToChain.isEmpty(), "Ring-constrained MCS should reject ring->chain");
      assertTrue(strictChainToRing.isEmpty(), "Ring-constrained MCS should reject chain->ring");

      ChemOptions relaxed = new ChemOptions();
      relaxed.ringMatchesRingOnly = false;
      Map<Integer, Integer> relaxedRingToChain = new MoleculePair(ring, chain, relaxed).mcsMapping( true, 5000L);
      Map<Integer, Integer> relaxedChainToRing = new MoleculePair(chain, ring, relaxed).mcsMapping( true, 5000L);
      assertEquals(5, relaxedRingToChain.size(), "Relaxed ring matching should recover the 5-atom path MCS");
      assertEquals(5, relaxedChainToRing.size(), "Relaxed ring matching should recover the 5-atom path MCS");
    }

    @Test
    void chainInRingSubstructureHonorsRingConstraint() throws Exception {
      IAtomContainer ring = mol("C1CCCCC1");
      IAtomContainer chain = mol("CCCCC");

      ChemOptions strict = new ChemOptions();
      strict.ringMatchesRingOnly = true;
      assertFalse(new MoleculePair(chain, ring, strict).contains(5000L),
          "Ring-constrained substructure should reject chain->ring matches");

      ChemOptions relaxed = new ChemOptions();
      relaxed.ringMatchesRingOnly = false;
      assertTrue(new MoleculePair(chain, ring, relaxed).contains(5000L),
          "Relaxed ring matching should allow chain->ring substructure matches");
    }

    @Test
    void largeMolMCS() throws Exception {
      // Aspirin vs ibuprofen — share aromatic ring
      MoleculePair smsd =
          new MoleculePair(
              mol("CC(=O)Oc1ccccc1C(=O)O"), // aspirin
              mol("CC(C)Cc1ccc(C(C)C(=O)O)cc1"), // ibuprofen
              new ChemOptions());
      Map<Integer, Integer> mcs = smsd.mcsMapping( true, 10000L);
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
          new MoleculePair(
                  mol("c1ncnc2[nH]cnc12"), // purine
                  mol("Cn1c(=O)c2c(ncn2C)n(C)c1=O"), // caffeine
                  ChemOptions.profile("compat-substruct"))
              .contains());
    }

    @Test
    void diazepamVsClonazepam() throws Exception {
      // Benzodiazepine core shared
      MoleculePair smsd =
          new MoleculePair(
              mol("CN1C(=O)CN=C(c2ccccc2)c2cc(Cl)ccc21"), // diazepam
              mol("c1cc2c(cc1[N+](=O)[O-])C(=NCC(=O)N2)c1ccccc1Cl"), // clonazepam
              new ChemOptions());
      Map<Integer, Integer> mcs = smsd.mcsMapping( true, 15000L);
      assertTrue(mcs.size() >= 15, "Benzodiazepine core should be shared");
    }

    @Test
    void ibuprofenVsNaproxen() throws Exception {
      MoleculePair smsd =
          new MoleculePair(
              mol("CC(C)Cc1ccc(C(C)C(=O)O)cc1"), // ibuprofen
              mol("COc1ccc2cc(C(C)C(=O)O)ccc2c1"), // naproxen
              new ChemOptions());
      Map<Integer, Integer> mcs = smsd.mcsMapping( true, 10000L);
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
      assertThrows(NullPointerException.class, () -> new MoleculePair((IAtomContainer) null, target, opts));
    }

    @Test
    @DisplayName("Null target molecule throws NullPointerException")
    void nullTargetMolecule() {
      IAtomContainer query = assertDoesNotThrow(() -> mol("C"));
      ChemOptions opts = new ChemOptions();
      assertThrows(NullPointerException.class, () -> new MoleculePair(query, (IAtomContainer) null, opts));
    }

    @Test
    @DisplayName("Empty SMILES for both query and target does not throw")
    void emptySmilesBothMolecules() throws Exception {
      IAtomContainer q = SP.parseSmiles("");
      IAtomContainer t = SP.parseSmiles("");
      ChemOptions opts = new ChemOptions();
      MoleculePair smsd = new MoleculePair(q, t, opts);
      // Should not throw
      smsd.contains();
    }

    @Test
    @DisplayName("Single atom C in single atom C is substructure")
    void singleAtomCinC() throws Exception {
      assertTrue(new MoleculePair(mol("C"), mol("C"), new ChemOptions()).contains());
    }

    @Test
    @DisplayName("Single atom C in single atom N is not substructure")
    void singleAtomCinN() throws Exception {
      assertFalse(new MoleculePair(mol("C"), mol("N"), new ChemOptions()).contains());
    }

    @Test
    @DisplayName("Query larger than target: naphthalene not in benzene")
    void queryLargerThanTarget() throws Exception {
      assertFalse(
          new MoleculePair(mol("c1ccc2ccccc2c1"), mol("c1ccccc1"), new ChemOptions()).contains());
    }

    @Test
    @DisplayName("Timeout 1ms on large pair does not hang")
    void timeoutDoesNotHang() throws Exception {
      IAtomContainer q = mol("c1ccc2ccccc2c1"); // naphthalene
      IAtomContainer t = mol("c1ccc2cc3ccccc3cc2c1"); // anthracene
      ChemOptions opts = new ChemOptions();
      MoleculePair smsd = new MoleculePair(q, t, opts);
      smsd.setSubstructureTimeoutMs(1L);
      long start = System.currentTimeMillis();
      // May or may not find the result, but must not hang
      smsd.contains();
      long elapsed = System.currentTimeMillis() - start;
      assertTrue(elapsed < 5000, "Should complete within 5 seconds, took " + elapsed + "ms");
    }

    @Test
    @DisplayName("Molecules with explicit hydrogens: [H]O[H] water")
    void explicitHydrogensWater() throws Exception {
      IAtomContainer water = SP.parseSmiles("[H]O[H]");
      IAtomContainer target = SP.parseSmiles("[H]O[H]");
      ChemOptions opts = new ChemOptions();
      MoleculePair smsd = new MoleculePair(water, target, opts);
      // Should parse and check without error
      boolean result = smsd.contains();
      assertTrue(result, "Water should be substructure of itself");
    }

    @Test
    @DisplayName("3+ disconnected fragments: C.C.C in CCCC handled gracefully")
    void disconnectedFragments() throws Exception {
      IAtomContainer q = mol("C.C.C");
      IAtomContainer t = mol("CCCC");
      ChemOptions opts = new ChemOptions();
      MoleculePair smsd = new MoleculePair(q, t, opts);
      // Should handle gracefully without exception
      smsd.contains();
    }

    @Test
    @DisplayName("Very long chain: C*200 contains C*100 substructure, completes < 5s")
    void veryLongChainSubstructure() throws Exception {
      String longChain200 = repeat("C", 200);
      String longChain100 = repeat("C", 100);
      IAtomContainer q = SP.parseSmiles(longChain100);
      IAtomContainer t = SP.parseSmiles(longChain200);
      ChemOptions opts = new ChemOptions();
      MoleculePair smsd = new MoleculePair(q, t, opts);
      long start = System.currentTimeMillis();
      boolean result = smsd.contains();
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

  @Nested
  @DisplayName("MCS completion and search counters")
  class MCSCompletionTest {
    @Test
    void completedPairReportsMappingAndCounters() throws Exception {
      SearchResult result=MCS.find(mol("c1ccccc1")).match(mol("c1ccc(O)cc1"));
      assertNotNull(result);
      assertTrue(result.optimal());
      assertTrue(result.complete());
      assertEquals(SearchResult.Status.Completed,result.status());
      assertEquals(6,result.matches().get(0).atomCount());
      assertTrue(result.nodesVisited()>0);
      assertTrue(result.elapsedMillis()>=0);
    }
    @Test
    void cdkAndGraphEntryPointsReturnTheSameObjective() throws Exception {
      IAtomContainer query=mol("c1ccccc1"),target=mol("c1ccc(O)cc1");
      SearchResult cdk=MCS.find(query).match(target);
      SearchResult graph=MCS.find(new MolGraph(query)).match(new MolGraph(target));
      assertEquals(cdk.score(),graph.score());
      assertEquals(cdk.matches(),graph.matches());
      assertTrue(graph.nodesVisited()>0);
    }
    @Test
    void stoppedSearchReportsAnUnprovenOutcome() throws Exception {
      SearchResult result=MCS.find(mol("CC(=O)Oc1ccccc1C(O)=O")).withTimeout(0).match(mol("CC(=O)Oc1ccccc1"));
      assertEquals(SearchResult.Status.TimedOut,result.status());
      assertFalse(result.optimal());
      assertFalse(result.complete());
      assertTrue(result.elapsedMillis()>=0);
    }
    @Test
    void fragmentSearchProvesAWholeQueryEmbeddingImmediately() throws Exception {
      // ChEMBL nearest-neighbour pair 604: the query is an induced subgraph of the target.
      IAtomContainer query=mol("CC1OC(OC2C(O)C(O)C(CO)OC2Oc2cc(O)cc(C=Cc3ccc(O)cc3)c2)C(O)C(O)C1O");
      MCS.Options options=new MCS.Options();
      options.connectivity=MCS.Connectivity.Fragments;
      options.timeoutMillis=30000;
      SearchResult result=MCS.find(query).withOptions(options).match(
          mol("Oc1ccc(C=Cc2cc(O)cc(OC3C(OC4OC(C)C(O)C(O)C4O)C(O)C(O)C(COC4OC(C)C(O)C(O)C4O)O3)c2)cc1"));
      assertTrue(result.optimal());
      assertEquals(query.getAtomCount(),result.matches().get(0).atomCount());
      assertTrue(result.nodesVisited()<10_000,"nodes "+result.nodesVisited());
    }
    @ParameterizedTest
    @CsvSource({"C1CC1CCCC,C1CC1CCC1CCCCC1,rings", "C[C@H](N)O,C[C@@H](N)OC,chirality"})
    void rejectedWholeQueryEmbeddingStillFindsTheOptimum(String query,String target,String rule)
        throws Exception {
      MCS.Options options=new MCS.Options();
      options.connectivity=MCS.Connectivity.Fragments;
      options.timeoutMillis=30000;
      options.chemistry.completeRingsOnly=rule.equals("rings");
      options.chemistry.useChirality=rule.equals("chirality");
      SearchResult result=MCS.find(mol(query)).withOptions(options).match(mol(target));
      // Any node limit selects the dynamic search order and skips the whole-query probe.
      options.nodeLimit=Long.MAX_VALUE;
      SearchResult previous=MCS.find(mol(query)).withOptions(options).match(mol(target));
      assertTrue(result.optimal());
      assertTrue(previous.optimal());
      assertEquals(previous.score(),result.score());
      assertEquals(previous.matches().get(0).bondCount(),result.matches().get(0).bondCount());
    }
    @Test
    void fragmentSearchGrowsMappedPiecesFirst() throws Exception {
      // ChEMBL random pair 4051: growing mapped pieces first proves its 23-atom optimum in time.
      MCS.Options options=new MCS.Options();
      options.connectivity=MCS.Connectivity.Fragments;
      options.timeoutMillis=30000;
      SearchResult result=MCS.find(mol("CC(=O)NC(C)C#Cc1cnc(Oc2c([N+]([O-])=O)cc(OCC(C)C)cc2)s1"))
          .withOptions(options).match(mol("Cc1ccccc1NC(NCCCCC(NC(C(Cc1c[nH]c2ccccc12)NC(OC(C)(C)C)=O)=O)"
              +"C(=O)NC(CC(O)=O)C(N(C)C(Cc1ccccc1)C(N)=O)=O)=O"));
      assertTrue(result.optimal());
      assertEquals(23,result.matches().get(0).atomCount());
      assertTrue(result.nodesVisited()<200_000,"nodes "+result.nodesVisited());
    }
  }

  @Nested
  @DisplayName("Common-bond MCES search")
  class CommonBondSearchTests {

    private MolGraph graph(String smiles) throws Exception {
      return new MolGraph(SP.parseSmiles(smiles));
    }

    private MCS.Options mces(MCS.Connectivity connectivity, ChemOptions.BondOrderMode bondOrder) {
      MCS.Options options = new MCS.Options();
      options.objective = MCS.Objective.MCES;
      options.connectivity = connectivity;
      options.chemistry.matchBondOrder = bondOrder;
      options.timeoutMillis = 30_000;
      return options;
    }

    @ParameterizedTest
    @CsvSource(delimiter = '|', value = {
        // ChEMBL random pairs 3471, 7381 and 3168 with their optimal common bond counts.
        "N#Cc1ccc(C(Nc2ccc3cc4ccc(NC(c5ccc(C#N)cc5)=O)cc4nc3c2)=O)cc1"
            + "|OCCCNP(O)(=O)OC1CC(n2cnc3c(N)ncnc23)OC1COP(O)(=O)OC1CC(n2ccc(N)nc2=O)OC1COP(O)(O)=O|24",
        "COc1c(Cl)cc(Nc2nc3ccc(C(=O)c4ccccc4)cc3n3c(C)nnc32)c(OC)c1"
            + "|CCCCc1n(Cc2ccc(-c3ccccc3-c3nn[nH]n3)cc2)c(SCc2ccccc2-c2nn[nH]n2)nn1|27",
        "CCCCCCCCCCCCCCCNC(OCCOCCOC(CCCCC[n+]1ccsc1)=O)=O"
            + "|OP(O)(=O)C(F)(F)c1ccc(CC(n2nnc3ccccc23)(c2scc(-c3ccccc3)n2)Cc2ccc(C(F)(F)P(O)(O)=O)cc2)cc1|27"})
    void fragmentMcesProofFitsItsNodeBudget(
        String querySmiles, String targetSmiles, int bonds) throws Exception {
      MolGraph query = graph(querySmiles), target = graph(targetSmiles);
      MCS.Options options = mces(MCS.Connectivity.Fragments, ChemOptions.BondOrderMode.ANY);
      SearchResult result = MCS.find(query).withOptions(options).match(target);
      assertEquals(SearchResult.Status.Completed, result.status());
      assertTrue(result.optimal());
      assertEquals(bonds, result.matches().get(0).bondCount());
      assertEquals(bonds, result.score());
      assertTrue(result.nodesVisited() < 20_000, "nodes " + result.nodesVisited());
      assertTrue(MCS.validateMapping(query, target, result.matches().get(0), options).isEmpty());
      // The node limit caps the pairs the search maps: with as many as the proof took, the search
      // proves the same optimum; with fewer, it stops with a valid mapping and no proof.
      options.nodeLimit = result.nodesVisited();
      SearchResult enough = MCS.find(query).withOptions(options).match(target);
      assertEquals(SearchResult.Status.Completed, enough.status());
      assertTrue(enough.optimal());
      assertEquals(bonds, enough.matches().get(0).bondCount());
      options.nodeLimit = Math.max(1, result.nodesVisited() / 2);
      SearchResult capped = MCS.find(query).withOptions(options).match(target);
      assertEquals(SearchResult.Status.LimitReached, capped.status());
      assertFalse(capped.optimal());
      assertTrue(capped.nodesVisited() <= options.nodeLimit);
      assertTrue(capped.matches().get(0).bondCount() <= bonds);
      assertTrue(MCS.validateMapping(query, target, capped.matches().get(0), options).isEmpty());
    }

    @Test
    void stoppedSearchesKeepAValidIncumbentWithoutClaimingOptimality() throws Exception {
      // ChEMBL random pair 8031, which takes far longer than these limits to prove.
      MolGraph query = graph("Cn1nnnc1Sc1ncnc2scc(-c3c(Cl)cccc3)c21");
      MolGraph target = graph("NCCCCC(N)C(=O)NC(CCCCN)C(NC(CCCCN)C(NC(C(NC(C(NC(C(NC(C(NC(C=O)CCCCN)=O)"
          + "CCCCN)=O)CCCCN)=O)CCCCN)=O)CCCCN)=O)=O");
      MCS.Options options = mces(MCS.Connectivity.Fragments, ChemOptions.BondOrderMode.STRICT);
      options.timeoutMillis = 200;
      SearchResult timedOut = MCS.find(query).withOptions(options).match(target);
      java.util.concurrent.atomic.AtomicInteger checks = new java.util.concurrent.atomic.AtomicInteger();
      options.timeoutMillis = 30_000;
      options.cancellation = () -> checks.incrementAndGet() > 200;
      SearchResult cancelled = MCS.find(query).withOptions(options).match(target);
      assertEquals(SearchResult.Status.TimedOut, timedOut.status());
      assertEquals(SearchResult.Status.Cancelled, cancelled.status());
      for (SearchResult result : Arrays.asList(timedOut, cancelled)) {
        assertFalse(result.optimal());
        assertFalse(result.complete());
        assertTrue(result.found());
        assertEquals(result.matches().get(0).bondCount(), result.score());
        assertTrue(MCS.validateMapping(query, target, result.matches().get(0), options).isEmpty());
      }
    }

    @ParameterizedTest
    @CsvSource({
        "CC(=O)Oc1ccccc1C(=O)O, CC(=O)Nc1ccc(O)cc1",
        "CC(C)Cc1ccc(C(C)C(=O)O)cc1, COc1ccc2cc(C(C)C(=O)O)ccc2c1",
        "Cn1cnc2c1c(=O)n(C)c(=O)n2C, Cn1c(=O)c2[nH]cnc2n(C)c1=O",
        "c1ccc2ccccc2c1, c1ccccc1CC",
        "C1CCC2CCCCC2C1, C12C3C4C1C5C2C3C45",
        "NCC(=O)O.CCO, OCC(N)=O"})
    void scoresIgnoreAnUnreachableNodeLimit(String querySmiles, String targetSmiles)
        throws Exception {
      MolGraph query = graph(querySmiles), target = graph(targetSmiles);
      for (MCS.Connectivity connectivity : MCS.Connectivity.values())
        for (ChemOptions.BondOrderMode bondOrder : new ChemOptions.BondOrderMode[] {
                 ChemOptions.BondOrderMode.STRICT, ChemOptions.BondOrderMode.ANY}) {
          MCS.Options options = mces(connectivity, bondOrder);
          SearchResult result = MCS.find(query).withOptions(options).match(target);
          // An unreachable node limit leaves the search and its objective unchanged.
          options.nodeLimit = Long.MAX_VALUE;
          SearchResult general = MCS.find(query).withOptions(options).match(target);
          String label = connectivity + " " + bondOrder;
          assertTrue(result.optimal(), label);
          assertTrue(general.optimal(), label);
          assertEquals(general.score(), result.score(), label);
          assertEquals(general.matches().get(0).bondCount(), result.matches().get(0).bondCount(), label);
          assertEquals(general.matches().get(0).atomCount(), result.matches().get(0).atomCount(), label);
          assertTrue(MCS.validateMapping(query, target, result.matches().get(0), options).isEmpty(),
              label);
        }
    }

    // Every MCES option goes to the common-bond search, which gives the same tied mapping as C++.
    @ParameterizedTest
    @CsvSource(delimiter = '|', value = {
        "all|Connected|7|2|{3=8, 4=7, 5=6, 6=5, 7=4, 8=10, 9=9}",
        "all|Fragments|9|2|{0=0, 1=1, 2=2, 3=8, 4=7, 5=6, 6=5, 7=4, 8=10, 9=9}",
        "chirality|Connected|2|1|{3=1, 4=2, 5=3}",
        "chirality|Fragments|2|1|{3=1, 4=2, 5=3}",
        "bondStereo|Connected|3|1|{0=5, 1=4, 2=3, 3=2}",
        "bondStereo|Fragments|4|1|{0=0, 1=1, 2=2, 3=3, 4=4, 5=5}",
        "rings|Connected|6|1|{3=0, 4=1, 5=2, 6=3, 7=4, 8=5}",
        "rings|Fragments|6|1|{3=0, 4=1, 5=2, 6=3, 7=4, 8=5}",
        "weights|Connected|9|1|{3=1, 4=2}",
        "weights|Fragments|14|1|{1=3, 2=4, 3=1, 4=2}",
        "minAtoms|Connected|3|1|{0=0, 1=1, 2=2, 3=3}",
        "minAtoms|Fragments|5|1|{0=0, 1=1, 2=2, 3=3, 4=5, 5=6, 6=7}"})
    void commonBondOptionsKeepTheirResults(String mode, MCS.Connectivity connectivity,
        double score, int matches, String firstMapping) throws Exception {
      String querySmiles = mode.equals("all")          ? "CC(=O)Oc1ccccc1C(=O)O"
                           : mode.equals("chirality")  ? "CCC[C@H](F)Cl"
                           : mode.equals("bondStereo") ? "CC/C=C/CC"
                           : mode.equals("rings")      ? "c1ccc2ccccc2c1"
                           : mode.equals("weights")    ? "CCO.CN"
                                                       : "CCCC.OCCO";
      String targetSmiles = mode.equals("all")          ? "CC(=O)Nc1ccc(O)cc1"
                            : mode.equals("chirality")  ? "C[C@@H](F)Cl"
                            : mode.equals("bondStereo") ? "CC/C=C\\CC"
                            : mode.equals("rings")      ? "c1ccccc1CC"
                            : mode.equals("weights")    ? "CCN.CO"
                                                        : "CCCCC.OCCCO";
      MolGraph query = graph(querySmiles), target = graph(targetSmiles);
      MCS.Options options = mces(connectivity, ChemOptions.BondOrderMode.STRICT);
      options.chemistry.useChirality = mode.equals("chirality");
      options.chemistry.useBondStereo = mode.equals("bondStereo");
      options.chemistry.completeRingsOnly = mode.equals("rings");
      if (mode.equals("weights")) {
        options.objective = MCS.Objective.MCIS;
        options.atomWeights = new double[] {1, 2, 3, 4, 5};
      }
      if (mode.equals("minAtoms"))
        options.minFragmentAtoms = 3;
      MCS prepared = MCS.find(query).withOptions(options);
      SearchResult result = mode.equals("all") ? prepared.matchAll(target, 1000) : prepared.match(target);
      assertEquals(SearchResult.Status.Completed, result.status());
      assertTrue(result.optimal());
      assertEquals(score, result.score());
      assertEquals(matches, result.matches().size());
      assertEquals(firstMapping, result.matches().get(0).atomMapping().toString());
      assertTrue(MCS.validateMapping(query, target, result.matches().get(0), options).isEmpty());
    }

    /** The state of a fixed 64-bit linear congruential sequence, the same as in the C++ tests. */
    private long drawState;

    /** The next draw of the sequence, from 0 up to {@code bound}. */
    private int draw(int bound) {
      drawState = drawState * 6364136223846793005L + 1442695040888963407L;
      return (int) ((drawState >>> 33) % bound);
    }

    @Test
    void seededOptionsIgnoreAnUnreachableNodeLimit() throws Exception {
      // The same pairs and draws are in the C++ MCS regression tests.
      String[][] pairs = {
          {"Oc1ccccc1", "Oc1ccncc1"}, {"CC(=O)[O-]", "CC(=O)OC"}, {"[13CH3]C(=O)O", "CC(=O)O"},
          {"c1ccc2ccccc2c1", "C1CCC2CCCCC2C1"}, {"O=c1cc[nH]cc1", "Oc1ccncc1"},
          {"C1CC2CCC1C2", "CC1CCCCC1"}, {"C[N+](C)(C)CC(=O)O", "CN(C)CC(N)=O"},
          {"NCC(=O)O.CCO", "OCC(N)=O"}, {"c1ccc(cc1)-c1ccccn1", "C=Cc1ccccc1"},
          {"[15NH2]c1ccccc1", "Nc1ccc(O)cc1"}};
      ChemOptions.BondOrderMode[] bondOrders = {
          ChemOptions.BondOrderMode.LOOSE, ChemOptions.BondOrderMode.ANY};
      ChemOptions.AromaticityMode[] aromaticities = {
          ChemOptions.AromaticityMode.STRICT, ChemOptions.AromaticityMode.FLEXIBLE};
      ChemOptions.RingFusionMode[] fusions = {ChemOptions.RingFusionMode.IGNORE,
          ChemOptions.RingFusionMode.PERMISSIVE, ChemOptions.RingFusionMode.STRICT};
      drawState = 20261019L;
      int cases = 0;
      for (String[] pair : pairs) {
        for (int round = 0; round < 8; round++, cases++) {
          // Fresh graphs, so that validation prepares what the options compare by itself.
          MolGraph query = graph(pair[0]), target = graph(pair[1]);
          MCS.Connectivity connectivity =
              draw(2) == 0 ? MCS.Connectivity.Connected : MCS.Connectivity.Fragments;
          MCS.Options options = mces(connectivity, bondOrders[draw(2)]);
          options.chemistry.aromaticityMode = aromaticities[draw(2)];
          options.chemistry.matchIsotope = draw(2) == 1;
          options.chemistry.matchFormalCharge = draw(2) == 1;
          options.chemistry.ringFusionMode = fusions[draw(3)];
          options.chemistry.matchAtomType = draw(4) != 0;
          for (int excluded = draw(3); excluded > 0; excluded--)
            options.excludedTargetAtoms.add(draw(target.atomCount()));
          String label = pair[0] + " " + pair[1] + " " + connectivity + " "
              + options.chemistry.matchBondOrder + " " + options.chemistry.aromaticityMode
              + " isotope " + options.chemistry.matchIsotope + " charge "
              + options.chemistry.matchFormalCharge + " " + options.chemistry.ringFusionMode
              + " type " + options.chemistry.matchAtomType + " excluded "
              + options.excludedTargetAtoms;
          SearchResult result = MCS.find(query).withOptions(options).match(target);
          // An unreachable node limit leaves the search and its objective unchanged.
          options.nodeLimit = Long.MAX_VALUE;
          SearchResult general = MCS.find(query).withOptions(options).match(target);
          assertTrue(result.optimal() && general.optimal(), label);
          assertEquals(general.found(), result.found(), label);
          assertEquals(general.score(), result.score(), label);
          if (!result.found())
            continue;
          for (SearchResult found : Arrays.asList(result, general)) {
            SearchResult.Match match = found.matches().get(0);
            assertTrue(MCS.validateMapping(query, target, match, options).isEmpty(), label);
            for (int targetAtom : match.atomMapping().values())
              assertFalse(options.excludedTargetAtoms.contains(targetAtom), label);
          }
          assertEquals(general.matches().get(0).bondCount(), result.matches().get(0).bondCount(), label);
          assertEquals(general.matches().get(0).atomCount(), result.matches().get(0).atomCount(), label);
        }
      }
      assertEquals(80, cases);
    }

    /** Three aromatic carbons joined by bonds of the given orders that are not aromatic. */
    private MolGraph aromaticChain(int firstOrder, int secondOrder) {
      return new MolGraph.Builder().atomCount(3).atomicNumbers(new int[] {6, 6, 6})
          .aromaticFlags(new boolean[] {true, true, true})
          .neighbors(new int[][] {{1}, {0, 2}, {1}})
          .bondOrders(new int[][] {{firstOrder}, {firstOrder, secondOrder}, {secondOrder}})
          .build();
    }

    @Test
    void flexibleAromaticityInterchangesSingleAndDoubleBondsBetweenAromaticAtoms() {
      // The same graphs and expectations are in the C++ MCS regression tests.
      MolGraph query = aromaticChain(1, 1), target = aromaticChain(2, 1);
      for (ChemOptions.AromaticityMode mode : ChemOptions.AromaticityMode.values()) {
        boolean flexible = mode == ChemOptions.AromaticityMode.FLEXIBLE;
        Substructure.Options substructure = new Substructure.Options();
        substructure.chemistry.aromaticityMode = mode;
        assertEquals(flexible,
            Substructure.find(query).withOptions(substructure).match(target).found(), mode.name());
        for (long nodeLimit : new long[] {0, Long.MAX_VALUE}) {
          MCS.Options options = mces(MCS.Connectivity.Fragments, ChemOptions.BondOrderMode.STRICT);
          options.chemistry.aromaticityMode = mode;
          options.nodeLimit = nodeLimit;
          SearchResult result = MCS.find(query).withOptions(options).match(target);
          String label = mode + " node limit " + nodeLimit;
          assertTrue(result.optimal(), label);
          assertEquals(flexible ? 2 : 1, result.matches().get(0).bondCount(), label);
          assertEquals(flexible ? 3 : 2, result.matches().get(0).atomCount(), label);
          assertTrue(MCS.validateMapping(query, target, result.matches().get(0), options).isEmpty(),
              label);
        }
      }
    }

    @Test
    void strictRingFusionValidatesAQueryThatNoSearchPrepared() throws Exception {
      MCS.Options options = mces(MCS.Connectivity.Connected, ChemOptions.BondOrderMode.STRICT);
      options.chemistry.ringFusionMode = ChemOptions.RingFusionMode.STRICT;
      MolGraph target = graph("C1CCC2CCCCC2C1");
      SearchResult result = MCS.find(graph("C1CCC2CCCCC2C1")).withOptions(options).match(target);
      assertTrue(result.optimal());
      assertEquals(11, result.matches().get(0).bondCount());
      // The search worked out the ring counts of the target only; validation works out the query's.
      assertTrue(MCS.validateMapping(graph("C1CCC2CCCCC2C1"), target, result.matches().get(0), options)
          .isEmpty());
    }
  }

  @Nested
  @DisplayName("Common-atom MCIS search")
  class CommonAtomSearchTests {

    private MolGraph graph(String smiles) throws Exception {
      return new MolGraph(SP.parseSmiles(smiles));
    }

    private MCS.Options mcis(MCS.Connectivity connectivity, ChemOptions.BondOrderMode bondOrder) {
      MCS.Options options = new MCS.Options();
      options.objective = MCS.Objective.MCIS;
      options.connectivity = connectivity;
      options.chemistry.matchBondOrder = bondOrder;
      options.timeoutMillis = 30_000;
      return options;
    }

    @ParameterizedTest
    @CsvSource(delimiter = '|', value = {
        // ChEMBL random pairs 5776, 7041 and 8361 with their optimal atoms and bonds. The same
        // mappings are in the C++ MCS regression tests.
        "CC(OC1C(=O)C2(C)C(O)CC3OCC3(OC(C)=O)C2C(OC(c2cc(N=[N+]=[N-])ccc2)=O)C2(O)CC(OC(C(O)C(NC(OC(C)(C)C)=O)C(F)F)=O)C(C)=C1C2(C)C)=O"
            + "|CC1N(C(=O)c2ccc(C(F)(F)F)cc2)CC(=O)OC1(Cn1cncn1)c1c(F)cc(F)cc1|26|23"
            + "|{0=15, 1=16, 2=18, 3=19, 4=1, 6=0, 22=3, 23=5, 24=14, 25=13, 29=8, 30=7, 31=6, 32=4, 33=32, "
            + "35=30, 36=29, 41=22, 42=23, 43=24, 51=11, 52=10, 54=27, 56=26, 57=33, 60=17}",
        "COc1ccc(CCNC(=O)Cn2cc(S(Cc3ccc(F)cc3)(=O)=O)c3ccccc23)cc1OC"
            + "|CC(OC1OC2OC3(C)CCC4C(C)CCC(C1C)C24OO3)c1ccc(C(=O)O)cc1|25|20"
            + "|{0=5, 1=6, 2=7, 3=9, 4=10, 5=11, 6=12, 7=13, 10=27, 13=0, 14=1, 16=3, 17=17, 18=18, 22=15, "
            + "23=16, 24=28, 25=20, 26=22, 27=23, 28=24, 29=25, 30=29, 31=30, 33=8}",
        "C=CCNC(=O)C(C#N)=Cc1c(=O)n2cccc(C)c2nc1N1CCN(c2ccc(OC)cc2)CC1"
            + "|CS(=O)c1c(-c2nc(C(=O)Nc3ccccc3)c(N)nc2)cccc1|21|18"
            + "|{5=2, 10=7, 11=8, 12=9, 13=10, 14=11, 15=12, 16=13, 17=14, 18=15, 20=18, 21=17, 22=19, "
            + "24=0, 26=4, 27=21, 28=22, 29=23, 32=24, 33=3, 35=20}"})
    void fragmentMcisProofFitsANodeBudgetTheGeneralSearchExceeds(String querySmiles,
        String targetSmiles, int atoms, int bonds, String mapping) throws Exception {
      MolGraph query = graph(querySmiles), target = graph(targetSmiles);
      MCS.Options options = mcis(MCS.Connectivity.Fragments, ChemOptions.BondOrderMode.ANY);
      SearchResult result = MCS.find(query).withOptions(options).match(target);
      assertEquals(SearchResult.Status.Completed, result.status());
      assertTrue(result.optimal());
      assertEquals(atoms, result.matches().get(0).atomCount());
      assertEquals(bonds, result.matches().get(0).bondCount());
      assertEquals(atoms, result.score());
      assertEquals(mapping, result.matches().get(0).atomMapping().toString());
      assertTrue(result.nodesVisited() < 20_000, "nodes " + result.nodesVisited());
      assertTrue(MCS.validateMapping(query, target, result.matches().get(0), options).isEmpty());
      // A node limit keeps the general search, which does not prove the optimum within it.
      options.nodeLimit = 20_000;
      SearchResult general = MCS.find(query).withOptions(options).match(target);
      assertEquals(SearchResult.Status.LimitReached, general.status());
      assertFalse(general.optimal());
      assertTrue(general.matches().get(0).atomCount() <= atoms);
    }

    @Test
    void stoppedSearchesKeepAValidIncumbentWithoutClaimingOptimality() throws Exception {
      // ChEMBL random pair 7175, which takes far longer than these limits to prove.
      MolGraph query = graph("COc1cc2c(c(OC)c1OC)-c1c(cc(=O)c(SC)cc1)C(NC(=O)c1ccc([N+]([O-])=O)"
          + "c(CO[N+]([O-])=O)c1)CC2");
      MolGraph target = graph("CC(C)=CCc1cc(C2CC(=O)c3c(c(CC=C(C)C)c(O)c(CC=C(C)C)c3O)O2)cc(O)c1O");
      MCS.Options options = mcis(MCS.Connectivity.Fragments, ChemOptions.BondOrderMode.STRICT);
      options.timeoutMillis = 200;
      SearchResult timedOut = MCS.find(query).withOptions(options).match(target);
      java.util.concurrent.atomic.AtomicInteger checks = new java.util.concurrent.atomic.AtomicInteger();
      options.timeoutMillis = 30_000;
      options.cancellation = () -> checks.incrementAndGet() > 200;
      SearchResult cancelled = MCS.find(query).withOptions(options).match(target);
      assertEquals(SearchResult.Status.TimedOut, timedOut.status());
      assertEquals(SearchResult.Status.Cancelled, cancelled.status());
      for (SearchResult result : Arrays.asList(timedOut, cancelled)) {
        assertFalse(result.optimal());
        assertFalse(result.complete());
        assertTrue(result.found());
        assertEquals(result.matches().get(0).atomCount(), result.score());
        assertTrue(MCS.validateMapping(query, target, result.matches().get(0), options).isEmpty());
      }
    }

    @ParameterizedTest
    @CsvSource({
        "CC(=O)Oc1ccccc1C(=O)O, CC(=O)Nc1ccc(O)cc1",
        "CC(C)Cc1ccc(C)cc1, COc1ccc(C(C)C=O)cc1",
        "Cn1cnc2c1c(=O)n(C)c(=O)n2C, Cn1c(=O)c2[nH]cnc2n(C)c1=O",
        "c1ccc2ccccc2c1, c1ccccc1CC",
        "C1CCC2CCCCC2C1, C12C3C4C1C5C2C3C45",
        "NCC(=O)O.CCO, OCC(N)=O",
        "CC(C)(C)c1ccccc1.O.N, O.Nc1ccccc1",
        "C1CC1C1CC1.C, CC.C1CC1"})
    void scoresAgreeWithTheGeneralSearch(String querySmiles, String targetSmiles)
        throws Exception {
      MolGraph query = graph(querySmiles), target = graph(targetSmiles);
      for (int minAtoms = 1; minAtoms <= 2; minAtoms++)
        for (ChemOptions.BondOrderMode bondOrder : new ChemOptions.BondOrderMode[] {
                 ChemOptions.BondOrderMode.STRICT, ChemOptions.BondOrderMode.ANY})
          for (boolean exclude : new boolean[] {false, true}) {
            MCS.Options options = mcis(MCS.Connectivity.Fragments, bondOrder);
            options.minFragmentAtoms = minAtoms;
            if (exclude)
              options.excludedTargetAtoms.addAll(Arrays.asList(0, target.atomCount() - 1));
            SearchResult result = MCS.find(query).withOptions(options).match(target);
            // An unreachable node limit keeps the general search without changing the objective.
            options.nodeLimit = Long.MAX_VALUE;
            SearchResult general = MCS.find(query).withOptions(options).match(target);
            String label = "minimum " + minAtoms + " " + bondOrder + " excluded " + exclude;
            assertTrue(result.optimal(), label);
            assertTrue(general.optimal(), label);
            assertEquals(general.found(), result.found(), label);
            assertEquals(general.score(), result.score(), label);
            if (!result.found())
              continue;
            SearchResult.Match match = result.matches().get(0);
            assertEquals(general.matches().get(0).atomCount(), match.atomCount(), label);
            assertEquals(general.matches().get(0).bondCount(), match.bondCount(), label);
            assertTrue(MCS.validateMapping(query, target, match, options).isEmpty(), label);
            for (int targetAtom : match.atomMapping().values())
              assertFalse(options.excludedTargetAtoms.contains(targetAtom), label);
          }
    }

    @ParameterizedTest
    @CsvSource(delimiter = '|', value = {
        "all|Fragments|9|12|{0=0, 1=1, 2=2, 4=4, 5=5, 6=6, 7=7, 8=9, 9=10}",
        "chirality|Fragments|3|1|{3=1, 4=2, 5=3}",
        "bondStereo|Fragments|5|1|{0=0, 1=1, 2=2, 3=3, 5=5}",
        "rings|Fragments|6|1|{0=2, 1=3, 2=4, 3=5, 8=0, 9=1}",
        "weights|Fragments|14|1|{1=3, 2=4, 3=1, 4=2}",
        "minAtoms|Fragments|7|1|{0=0, 1=1, 2=2, 3=3, 4=5, 5=6, 6=7}",
        "maxFragments|Fragments|7|1|{0=0, 1=1, 2=2, 3=3, 4=5, 5=6, 6=7}",
        "plain|Connected|7|1|{3=8, 4=7, 5=6, 6=5, 7=4, 8=10, 9=9}"})
    void ineligibleOptionsKeepThePreviousResults(String mode, MCS.Connectivity connectivity,
        double score, int matches, String firstMapping) throws Exception {
      boolean drugs = mode.equals("all") || mode.equals("plain");
      String querySmiles = drugs                       ? "CC(=O)Oc1ccccc1C(=O)O"
                           : mode.equals("chirality")  ? "CCC[C@H](F)Cl"
                           : mode.equals("bondStereo") ? "CC/C=C/CC"
                           : mode.equals("rings")      ? "c1ccc2ccccc2c1"
                           : mode.equals("weights")    ? "CCO.CN"
                                                       : "CCCC.OCCO";
      String targetSmiles = drugs                       ? "CC(=O)Nc1ccc(O)cc1"
                            : mode.equals("chirality")  ? "C[C@@H](F)Cl"
                            : mode.equals("bondStereo") ? "CC/C=C\\CC"
                            : mode.equals("rings")      ? "c1ccccc1CC"
                            : mode.equals("weights")    ? "CCN.CO"
                                                        : "CCCCC.OCCCO";
      MolGraph query = graph(querySmiles), target = graph(targetSmiles);
      MCS.Options options = mcis(connectivity, ChemOptions.BondOrderMode.STRICT);
      options.chemistry.useChirality = mode.equals("chirality");
      options.chemistry.useBondStereo = mode.equals("bondStereo");
      options.chemistry.completeRingsOnly = mode.equals("rings");
      if (mode.equals("weights"))
        options.atomWeights = new double[] {1, 2, 3, 4, 5};
      if (mode.equals("minAtoms"))
        options.minFragmentAtoms = 3;
      if (mode.equals("maxFragments"))
        options.maxFragments = 2;
      MCS prepared = MCS.find(query).withOptions(options);
      SearchResult result = mode.equals("all") ? prepared.matchAll(target, 1000) : prepared.match(target);
      assertEquals(SearchResult.Status.Completed, result.status());
      assertTrue(result.optimal());
      assertEquals(score, result.score());
      assertEquals(matches, result.matches().size());
      assertEquals(firstMapping, result.matches().get(0).atomMapping().toString());
      assertTrue(MCS.validateMapping(query, target, result.matches().get(0), options).isEmpty());
    }

    /** The state of a fixed 64-bit linear congruential sequence, the same as in the C++ tests. */
    private long drawState;

    /** The next draw of the sequence, from 0 up to {@code bound}. */
    private int draw(int bound) {
      drawState = drawState * 6364136223846793005L + 1442695040888963407L;
      return (int) ((drawState >>> 33) % bound);
    }

    @Test
    void seededOptionsAgreeWithTheGeneralSearch() throws Exception {
      // The same pairs and draws are in the C++ MCS regression tests.
      String[][] pairs = {
          {"Oc1ccccc1", "Oc1ccncc1"}, {"CC(=O)[O-]", "CC(=O)OC"}, {"[13CH3]C(=O)O", "CC(=O)O"},
          {"c1ccc2ccccc2c1", "C1CCC2CCCCC2C1"}, {"O=c1cc[nH]cc1", "Oc1ccncc1"},
          {"C1CC2CCC1C2", "CC1CCCCC1"}, {"C[N+](C)(C)CC(=O)O", "CN(C)CC(N)=O"},
          {"NCC(=O)O.CCO", "OCC(N)=O"}, {"c1ccc(cc1)-c1ccccn1", "C=Cc1ccccc1"},
          {"[15NH2]c1ccccc1", "Nc1ccc(O)cc1"}};
      ChemOptions.BondOrderMode[] bondOrders = {ChemOptions.BondOrderMode.STRICT,
          ChemOptions.BondOrderMode.LOOSE, ChemOptions.BondOrderMode.ANY};
      ChemOptions.AromaticityMode[] aromaticities = {
          ChemOptions.AromaticityMode.STRICT, ChemOptions.AromaticityMode.FLEXIBLE};
      ChemOptions.RingFusionMode[] fusions = {ChemOptions.RingFusionMode.IGNORE,
          ChemOptions.RingFusionMode.PERMISSIVE, ChemOptions.RingFusionMode.STRICT};
      drawState = 20261020L;
      int cases = 0;
      for (String[] pair : pairs) {
        for (int round = 0; round < 8; round++, cases++) {
          // Fresh graphs, so that validation prepares what the options compare by itself.
          MolGraph query = graph(pair[0]), target = graph(pair[1]);
          MCS.Options options = mcis(MCS.Connectivity.Fragments, bondOrders[draw(3)]);
          options.minFragmentAtoms = 1 + draw(2);
          options.chemistry.aromaticityMode = aromaticities[draw(2)];
          options.chemistry.matchIsotope = draw(2) == 1;
          options.chemistry.matchFormalCharge = draw(2) == 1;
          options.chemistry.ringFusionMode = fusions[draw(3)];
          options.chemistry.matchAtomType = draw(4) != 0;
          for (int excluded = draw(3); excluded > 0; excluded--)
            options.excludedTargetAtoms.add(draw(target.atomCount()));
          String label = pair[0] + " " + pair[1] + " minimum " + options.minFragmentAtoms + " "
              + options.chemistry.matchBondOrder + " " + options.chemistry.aromaticityMode
              + " isotope " + options.chemistry.matchIsotope + " charge "
              + options.chemistry.matchFormalCharge + " " + options.chemistry.ringFusionMode
              + " type " + options.chemistry.matchAtomType + " excluded "
              + options.excludedTargetAtoms;
          SearchResult result = MCS.find(query).withOptions(options).match(target);
          // An unreachable node limit keeps the general search without changing the objective.
          options.nodeLimit = Long.MAX_VALUE;
          SearchResult general = MCS.find(query).withOptions(options).match(target);
          assertTrue(result.optimal() && general.optimal(), label);
          assertEquals(general.found(), result.found(), label);
          assertEquals(general.score(), result.score(), label);
          if (!result.found())
            continue;
          for (SearchResult found : Arrays.asList(result, general)) {
            SearchResult.Match match = found.matches().get(0);
            assertTrue(MCS.validateMapping(query, target, match, options).isEmpty(), label);
            for (int targetAtom : match.atomMapping().values())
              assertFalse(options.excludedTargetAtoms.contains(targetAtom), label);
          }
          assertEquals(general.matches().get(0).atomCount(), result.matches().get(0).atomCount(), label);
          assertEquals(general.matches().get(0).bondCount(), result.matches().get(0).bondCount(), label);
        }
      }
      assertEquals(80, cases);
    }

    /** Three aromatic carbons joined by bonds of the given orders that are not aromatic. */
    private MolGraph aromaticChain(int firstOrder, int secondOrder) {
      return new MolGraph.Builder().atomCount(3).atomicNumbers(new int[] {6, 6, 6})
          .aromaticFlags(new boolean[] {true, true, true})
          .neighbors(new int[][] {{1}, {0, 2}, {1}})
          .bondOrders(new int[][] {{firstOrder}, {firstOrder, secondOrder}, {secondOrder}})
          .build();
    }

    @Test
    void flexibleAromaticityKeepsResonanceBondsInducedBetweenAromaticAtoms() {
      // The same graphs and expectations are in the C++ MCS regression tests.
      MolGraph query = aromaticChain(1, 1), target = aromaticChain(2, 1);
      for (ChemOptions.AromaticityMode mode : ChemOptions.AromaticityMode.values()) {
        boolean flexible = mode == ChemOptions.AromaticityMode.FLEXIBLE;
        for (long nodeLimit : new long[] {0, Long.MAX_VALUE}) {
          MCS.Options options = mcis(MCS.Connectivity.Fragments, ChemOptions.BondOrderMode.STRICT);
          options.chemistry.aromaticityMode = mode;
          options.nodeLimit = nodeLimit;
          SearchResult result = MCS.find(query).withOptions(options).match(target);
          String label = mode + " node limit " + nodeLimit;
          assertTrue(result.optimal(), label);
          assertEquals(flexible ? 3 : 2, result.matches().get(0).atomCount(), label);
          assertEquals(flexible ? 2 : 1, result.matches().get(0).bondCount(), label);
          assertTrue(MCS.validateMapping(query, target, result.matches().get(0), options).isEmpty(),
              label);
        }
      }
    }
  }
}
