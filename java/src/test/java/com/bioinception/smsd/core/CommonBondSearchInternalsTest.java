/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
package com.bioinception.smsd.core;

import static com.bioinception.smsd.TestSupport.rawGraph;
import static org.junit.jupiter.api.Assertions.*;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import org.junit.jupiter.api.Test;

/** Internals of the common-bond search that the public searches reach only in long runs. */
class CommonBondSearchInternalsTest {
  /** Carbons in a cycle, each bonded to the next {@code step} atoms around it. */
  private static MolGraph circulant(int atoms, int step) {
    int[][] neighbors = new int[atoms][2 * step], orders = new int[atoms][2 * step];
    int[] listed = new int[atoms];
    for (int atom = 0; atom < atoms; atom++)
      for (int offset = 1; offset <= step; offset++) {
        int other = (atom + offset) % atoms;
        neighbors[atom][listed[atom]++] = other;
        neighbors[other][listed[other]++] = atom;
      }
    for (int[] row : orders)
      Arrays.fill(row, 1);
    int[] carbons = new int[atoms];
    Arrays.fill(carbons, 6);
    return new MolGraph.Builder().atomCount(atoms).atomicNumbers(carbons).neighbors(neighbors)
        .bondOrders(orders).build();
  }

  /** A carbon bonded to {@code leaves} atoms: an oxygen first when {@code oxygen}, and carbons otherwise. */
  private static MolGraph star(int leaves, boolean oxygen) {
    int[] elements = new int[leaves + 1];
    Arrays.fill(elements, 6);
    if (oxygen)
      elements[1] = 8;
    int[][] neighbors = new int[leaves + 1][], orders = new int[leaves + 1][];
    neighbors[0] = new int[leaves];
    orders[0] = new int[leaves];
    for (int leaf = 1; leaf <= leaves; leaf++) {
      neighbors[0][leaf - 1] = leaf;
      orders[0][leaf - 1] = 1;
      neighbors[leaf] = new int[] {0};
      orders[leaf] = new int[] {1};
    }
    return new MolGraph.Builder().atomCount(leaves + 1).atomicNumbers(elements).neighbors(neighbors)
        .bondOrders(orders).build();
  }

  /** What a common-atom or common-bond search found when its clock stopped, or at its end. */
  private static final class Polled {
    boolean finished;
    long polls;
    SearchResult.Match match;
  }

  /** The search of {@code options}, which the clock stops after {@code limit} polls; a negative limit runs it to its end. */
  private static Polled polled(MolGraph query, MolGraph target, MCS.Options options, long limit) {
    Polled result = new Polled();
    CommonBondGraphs.Clock clock = new CommonBondGraphs.Clock(() -> ++result.polls > limit && limit >= 0);
    CommonBondGraphs graph = null;
    CommonAtomSearch atoms = null;
    CommonBondSearch bonds = null;
    int[] found = null;
    try {
      graph = new CommonBondGraphs(query, target, options.chemistry, Collections.emptySet(), clock);
      if (options.objective == MCS.Objective.MCIS) {
        atoms = new CommonAtomSearch(graph, clock, options.minFragmentAtoms);
        found = atoms.run();
      } else {
        bonds = new CommonBondSearch(graph, clock, options.connectivity == MCS.Connectivity.Connected);
        List<int[]> all = bonds.run();
        found = all.isEmpty() ? null : all.get(0);
      }
      result.finished = true;
    } catch (CommonBondGraphs.Stop stop) {
      if (atoms != null)
        found = atoms.found();
      if (bonds != null)
        found = bonds.found().isEmpty() ? null : bonds.found().get(0);
    }
    if (found != null)
      result.match = options.objective == MCS.Objective.MCIS
          ? CommonAtomSearch.match(graph, found, false) : CommonBondSearch.match(graph, found);
    return result;
  }

  private static MCS.Options options(MCS.Objective objective, MCS.Connectivity connectivity) {
    MCS.Options options = new MCS.Options();
    options.objective = objective;
    options.connectivity = connectivity;
    options.timeoutMillis = 60_000;
    return options;
  }

  @Test
  void bondPairsBeyondTheLimitAreDeclined() {
    // 8192 target bonds take 128 words a query bond, so 32768 query bonds fill the limit exactly,
    // as in the native core.
    MolGraph target = circulant(8192, 1);
    assertEquals(CommonBondGraphs.MAX_BOND_PAIR_WORDS, 32768L * 128);
    assertTrue(CommonBondGraphs.supports(circulant(32768, 1), target));
    assertFalse(CommonBondGraphs.supports(circulant(32769, 1), target));
    assertTrue(CommonBondGraphs.supports(target, circulant(32768, 1)));
    // A dense pair beyond the limit is declined at once: no match, with the limit reached.
    MolGraph query = circulant(600, 40), dense = circulant(640, 40);
    assertFalse(CommonBondGraphs.supports(query, dense));
    for (MCS.Connectivity connectivity : MCS.Connectivity.values()) {
      MCS.Options options = options(MCS.Objective.MCES, connectivity);
      options.timeoutMillis = 200;
      SearchResult result = MCS.find(query).withOptions(options).match(dense);
      assertEquals(SearchResult.Status.LimitReached, result.status(), connectivity.toString());
      assertFalse(result.optimal(), connectivity.toString());
      assertTrue(result.matches().isEmpty(), connectivity.toString());
      assertTrue(result.elapsedMillis() < 5000, connectivity + " took " + result.elapsedMillis() + " ms");
    }
  }

  @Test
  void starSeedsGrowQuickly() {
    // A seed passes the taken targets next to a hub once, and looks for pairs only at the
    // frontier atoms that might beat the best pair so far, so seeding a star is not cubic in its
    // leaves.
    for (MCS.Objective objective : MCS.Objective.values()) {
      boolean bonds = objective == MCS.Objective.MCES;
      MolGraph query = star(2000, bonds), target = star(2000, false);
      MCS.Options options = options(objective, MCS.Connectivity.Fragments);
      Polled result = polled(query, target, options, -1);
      assertTrue(result.finished && result.polls < 100_000, objective + " polls " + result.polls);
      assertEquals(bonds ? 2000 : 2001, result.match.atomCount(), objective.toString());
      assertEquals(bonds ? 1999 : 2000, result.match.bondCount(), objective.toString());
      assertTrue(MCS.validateMapping(query, target, result.match, options).isEmpty(), objective.toString());
    }
  }

  @Test
  void searchesStoppedWhileSeedingKeepTheirSeed() {
    // Stopped a few polls into seeding, which takes most of these searches, each search returns
    // the seed it was growing: a valid mapping rather than nothing.
    MolGraph leaves = star(400, false), oxygen = star(400, true);
    for (MCS.Connectivity connectivity : MCS.Connectivity.values())
      for (MCS.Objective objective : MCS.Objective.values()) {
        boolean bonds = objective == MCS.Objective.MCES;
        if (!bonds && connectivity == MCS.Connectivity.Connected)
          continue;
        MolGraph query = bonds ? oxygen : leaves;
        MCS.Options options = options(objective, connectivity);
        String label = objective + " " + connectivity;
        Polled whole = polled(query, leaves, options, -1);
        assertTrue(whole.finished && whole.polls > 800, label + " polls " + whole.polls);
        Polled stopped = polled(query, leaves, options, 800);
        assertFalse(stopped.finished, label);
        assertNotNull(stopped.match, label);
        assertTrue(stopped.match.atomCount() > 0 && stopped.match.bondCount() > 0, label);
        assertTrue(MCS.validateMapping(query, leaves, stopped.match, options).isEmpty(), label);
      }
  }

  /** Oxygens bonded to nothing, with mass numbers 1, 2, ... when {@code labelled}. */
  private static MolGraph isolated(int atoms, boolean labelled) {
    int[] elements = new int[atoms], masses = new int[atoms];
    Arrays.fill(elements, 8);
    for (int atom = 0; atom < atoms; atom++)
      masses[atom] = labelled ? atom + 1 : 0;
    return new MolGraph.Builder().atomCount(atoms).atomicNumbers(elements).massNumbers(masses)
        .neighbors(new int[atoms][0]).build();
  }

  @Test
  void isolatedAtomSeedsGrowIncrementally() {
    // Each seed round of isolated atoms keeps one more atom. A round takes its roots from runs of
    // target cells and passes a target that can no longer root a seed only once, so seeding four
    // thousand atoms takes about forty thousand polls, as in the native core; scanning every pair
    // of atoms in each round would take about fifty million.
    MolGraph graph = isolated(4000, false);
    MCS.Options options = options(MCS.Objective.MCIS, MCS.Connectivity.Fragments);
    Polled whole = polled(graph, graph, options, -1);
    assertTrue(whole.finished && whole.polls < 200_000, "polls " + whole.polls);
    assertEquals(4000, whole.match.atomCount());
    assertTrue(MCS.validateMapping(graph, graph, whole.match, options).isEmpty());
    // Stopped part way through seeding, the search returns the atoms kept so far.
    Polled stopped = polled(graph, graph, options, 36_000);
    assertFalse(stopped.finished);
    assertTrue(stopped.match != null && stopped.match.atomCount() >= 2000,
        "atoms " + (stopped.match == null ? 0 : stopped.match.atomCount()));
    assertTrue(MCS.validateMapping(graph, graph, stopped.match, options).isEmpty());
  }

  @Test
  void kindsAddNoPairTable() {
    // Labelled atoms each have a kind of their own. Whether two kinds match is read from the
    // pair rows, so the graphs allocate no table of every pair of kinds, which would take 144 MB
    // here.
    java.lang.management.ThreadMXBean threads = java.lang.management.ManagementFactory.getThreadMXBean();
    org.junit.jupiter.api.Assumptions.assumeTrue(threads instanceof com.sun.management.ThreadMXBean);
    com.sun.management.ThreadMXBean allocations = (com.sun.management.ThreadMXBean) threads;
    org.junit.jupiter.api.Assumptions.assumeTrue(allocations.isThreadAllocatedMemorySupported()
        && allocations.isThreadAllocatedMemoryEnabled());
    MolGraph graph = isolated(12_000, true);
    ChemOptions chemistry = new ChemOptions();
    chemistry.matchIsotope = true;
    CommonBondGraphs.Clock clock = new CommonBondGraphs.Clock(() -> false);
    long thread = Thread.currentThread().getId();
    long before = allocations.getThreadAllocatedBytes(thread);
    CommonBondGraphs prepared = new CommonBondGraphs(graph, graph, chemistry, Collections.emptySet(), clock);
    long allocated = allocations.getThreadAllocatedBytes(thread) - before;
    // The pair rows and their transpose take 18 MB each.
    assertTrue(allocated < 80L << 20, "allocated " + (allocated >> 20) + " MiB");
    assertEquals(12_000, prepared.query.atomCount);
  }

  @Test
  void scoutingFromTheFirstNodesKeepsTheGeneralOptimum() throws Exception {
    String[][] pairs = {
        {"CC(=O)Oc1ccccc1C(=O)O", "CC(=O)Nc1ccc(O)cc1"},
        {"CC(C)Cc1ccc(C(C)C(=O)O)cc1", "COc1ccc2cc(C(C)C(=O)O)ccc2c1"},
        {"Cn1cnc2c1c(=O)n(C)c(=O)n2C", "Cn1c(=O)c2[nH]cnc2n(C)c1=O"},
        {"c1ccc2ccccc2c1", "c1ccccc1CC"},
        {"C1CCC2CCCCC2C1", "C12C3C4C1C5C2C3C45"},
        {"NCC(=O)O.CCO", "OCC(N)=O"},
        {"OC(=O)c1ccccc1O", "Nc1ccc(cc1)S(N)(=O)=O"},
        {"CCN(CC)CCOC(=O)c1ccc(N)cc1", "CN1CCC(CC1)OC(=O)c1ccccc1"}};
    CommonBondSearch.ScoutLimits limits = new CommonBondSearch.ScoutLimits(1, 2, 32, 4);
    long scoutings = 0, refuted = 0;
    for (String[] pair : pairs) {
      MolGraph query = rawGraph(pair[0]), target = rawGraph(pair[1]);
      for (MCS.Connectivity connectivity : MCS.Connectivity.values())
        for (ChemOptions.BondOrderMode bondOrder : new ChemOptions.BondOrderMode[] {
                 ChemOptions.BondOrderMode.STRICT, ChemOptions.BondOrderMode.ANY}) {
          MCS.Options options = new MCS.Options();
          options.objective = MCS.Objective.MCES;
          options.connectivity = connectivity;
          options.chemistry.matchBondOrder = bondOrder;
          options.timeoutMillis = 60_000;
          // An unreachable node limit leaves the search and its objective unchanged.
          options.nodeLimit = Long.MAX_VALUE;
          SearchResult general = MCS.find(query).withOptions(options).match(target);
          CommonBondGraphs.Clock clock = new CommonBondGraphs.Clock(() -> false);
          CommonBondGraphs graph = new CommonBondGraphs(
              query, target, options.chemistry, Collections.emptySet(), clock);
          CommonBondSearch search = new CommonBondSearch(
              graph, clock, connectivity == MCS.Connectivity.Connected, limits);
          List<int[]> found = search.run();
          String label = pair[0] + " " + pair[1] + " " + connectivity + " " + bondOrder;
          assertTrue(general.optimal(), label);
          assertFalse(found.isEmpty(), label);
          SearchResult.Match match = CommonBondSearch.match(graph, found.get(0));
          assertEquals(general.matches().get(0).bondCount(), match.bondCount(), label);
          assertEquals(general.matches().get(0).atomCount(), match.atomCount(), label);
          scoutings += search.scoutings();
          refuted += search.refutedRoots();
        }
    }
    assertTrue(scoutings > 0, "scoutings " + scoutings);
    assertTrue(refuted > 0, "refuted root pairs " + refuted);
  }

  /** A molecule from SMILES, standardised with no tautomers, as the benchmark runners read them. */
  private static MolGraph standardised(String smiles) throws Exception {
    return new MolGraph(Standardiser.standardise(
        com.bioinception.smsd.TestSupport.rawMolecule(smiles),
        Standardiser.TautomerMode.NONE));
  }

  /**
   * Fragment pairs searched with the default limits, with the bonds, atoms and nodes of their
   * optimum pinned to the counts of the native search, so both searches walk the same tree. The
   * chain pair starts its last search from the seed grown the other way round, which holds the
   * optimum, whichever molecule is the query. The fused pair is proved once the search the other
   * way round is raced against its last search.
   */
  @Test
  void defaultLimitsMatchTheOtherLanguage() throws Exception {
    String chain = "CCCCCCCCCCCCCCSCC(COP(O)(OP(O)(OCC1OC(n2ccc(N)nc2=O)C(O)C1O)=O)=O)OC(=O)CCCCCCCCCCC";
    String amine = "CN(CCOC(c1ccccc1)c1ccccc1)CCN(C)CCc1ccc(F)cc1";
    String fused = "CCC1(O)N2CCC3c4c(cccc4)NC3C(c3c(OC)cc4c(c3)C35CCN6CC=CC(CC)"
        + "(C(OC(CC7CCCN7)=O)C(O)(C(OC)=O)C3N4C)C56)C(C(OC)=O)C(CC1)C2";
    String triazine = "Oc1nc(NCc2cc3[nH]ccc3cc2)nc(Nc2cc(Oc3ccccc3)ccc2)n1";
    String[][] pinned = {
        {chain, amine, "26", "30", "27734"},
        {amine, chain, "26", "30", "27734"},
        {fused, triazine, "29", "31", "677291"},
        {triazine, fused, "29", "31", "677291"}};
    MCS.Options options = options(MCS.Objective.MCES, MCS.Connectivity.Fragments);
    // Elements and any bond order, with flexible aromaticity.
    options.chemistry.matchAtomType = true;
    options.chemistry.matchBondOrder = ChemOptions.BondOrderMode.ANY;
    options.chemistry.aromaticityMode = ChemOptions.AromaticityMode.FLEXIBLE;
    options.chemistry.matchFormalCharge = options.chemistry.matchIsotope = false;
    options.chemistry.useChirality = options.chemistry.useBondStereo = false;
    options.chemistry.ringMatchesRingOnly = options.chemistry.completeRingsOnly = false;
    options.chemistry.tautomerAware = false;
    options.chemistry.ringFusionMode = ChemOptions.RingFusionMode.IGNORE;
    for (String[] pair : pinned) {
      SearchResult result =
          MCS.find(standardised(pair[0])).withOptions(options).match(standardised(pair[1]));
      assertTrue(result.optimal(), pair[0]);
      assertEquals(Integer.parseInt(pair[2]), result.matches().get(0).bondCount(), pair[0]);
      assertEquals(Integer.parseInt(pair[3]), result.matches().get(0).atomCount(), pair[0]);
      assertEquals(Long.parseLong(pair[4]), result.nodesVisited(), pair[0]);
    }
  }

  /**
   * The fused pair stopped by node limits before the race, within the first nodes of the search
   * the other way round, within its longer run and one node before it proves the optimum, with
   * the bonds and atoms of the best mapping held there pinned to the native search; both searches
   * together map exactly the limit, whichever molecule is the query.
   */
  @Test
  void nodeLimitsInsideTheRaceMatchTheOtherLanguage() throws Exception {
    String fused = "CCC1(O)N2CCC3c4c(cccc4)NC3C(c3c(OC)cc4c(c3)C35CCN6CC=CC(CC)"
        + "(C(OC(CC7CCCN7)=O)C(O)(C(OC)=O)C3N4C)C56)C(C(OC)=O)C(CC1)C2";
    String triazine = "Oc1nc(NCc2cc3[nH]ccc3cc2)nc(Nc2cc(Oc3ccccc3)ccc2)n1";
    long[][] capped = {{530000, 27, 31}, {545000, 28, 30}, {600000, 29, 30}, {677290, 29, 31}};
    MCS.Options options = options(MCS.Objective.MCES, MCS.Connectivity.Fragments);
    options.chemistry.matchAtomType = true;
    options.chemistry.matchBondOrder = ChemOptions.BondOrderMode.ANY;
    options.chemistry.aromaticityMode = ChemOptions.AromaticityMode.FLEXIBLE;
    options.chemistry.matchFormalCharge = options.chemistry.matchIsotope = false;
    options.chemistry.useChirality = options.chemistry.useBondStereo = false;
    options.chemistry.ringMatchesRingOnly = options.chemistry.completeRingsOnly = false;
    options.chemistry.tautomerAware = false;
    options.chemistry.ringFusionMode = ChemOptions.RingFusionMode.IGNORE;
    for (long[] cap : capped) {
      options.nodeLimit = cap[0];
      for (String[] pair : new String[][] {{fused, triazine}, {triazine, fused}}) {
        SearchResult result =
            MCS.find(standardised(pair[0])).withOptions(options).match(standardised(pair[1]));
        String label = cap[0] + " " + pair[0];
        assertEquals(SearchResult.Status.LimitReached, result.status(), label);
        assertEquals(cap[0], result.nodesVisited(), label);
        assertEquals(cap[1], result.matches().get(0).bondCount(), label);
        assertEquals(cap[2], result.matches().get(0).atomCount(), label);
      }
    }
  }

  /**
   * With the search the other way round raced against the last fragment search after 64 nodes,
   * for 16 nodes and up to 64 while it beats it, fragment searches prove the optimum they prove
   * without a race, in the nodes that the native search counts for them.
   */
  @Test
  void racesAfterAFewNodesKeepTheOptimum() throws Exception {
    String[][] pairs = {
        {"COc1ccc(C(c2n(CC3OCCC3)nnn2)N2CCN(c3ccccc3)CC2)cc1", "COC(=O)C1=CC23C4N(CCC44c5ccccc5NC14CC2)CCC3", "6830", "7115"},
        {"C=CCN1CCC23CC(=O)CCC2C1Cc1cccc(OC(C)=O)c13", "CN1CCN(C(=O)C2=CC(c3ccc(C(F)(F)F)cc3)CC(OCCCCO)O2)CC1", "47001", "22312"},
        {"Nc1ncnc2n(C3OC(CO)C(F)C3O)cnc12", "Cc1ccccc1Cn1c2c(scc2)c(=O)n(CCCC(NCc2occc2)=O)c1=O", "1179", "39195"},
        {"COc1cc2c(cc1)C(=O)CC(c1ccc(C)cc1)O2", "NS(c1cc(NC(COC(C=Cc2nc3ccccc3s2)=O)=O)ccc1)(=O)=O", "5185", "15555"}};
    CommonBondSearch.ScoutLimits raced = new CommonBondSearch.ScoutLimits(1L << 15, 4, 1L << 15, 1L << 10, 64, 16, 64);
    CommonBondSearch.ScoutLimits never =
        new CommonBondSearch.ScoutLimits(1L << 15, 4, 1L << 15, 1L << 10, Long.MAX_VALUE, 1L << 14, 1L << 19);
    for (String[] pair : pairs) {
      MolGraph query = standardised(pair[0]), target = standardised(pair[1]);
      for (ChemOptions.BondOrderMode order : new ChemOptions.BondOrderMode[] {
               ChemOptions.BondOrderMode.STRICT, ChemOptions.BondOrderMode.ANY}) {
        MCS.Options options = options(MCS.Objective.MCES, MCS.Connectivity.Fragments);
        options.chemistry.matchBondOrder = order;
        SearchResult race = CommonBondSearch.search(query, target, options, 1, false, null,
            new SearchResult.Control(60_000, null), raced);
        SearchResult plain = CommonBondSearch.search(query, target, options, 1, false, null,
            new SearchResult.Control(60_000, null), never);
        String label = pair[0] + " " + order;
        assertTrue(race.optimal() && plain.optimal(), label);
        assertEquals(plain.matches().get(0).bondCount(), race.matches().get(0).bondCount(), label);
        assertEquals(plain.matches().get(0).atomCount(), race.matches().get(0).atomCount(), label);
        assertEquals(Long.parseLong(order == ChemOptions.BondOrderMode.STRICT ? pair[2] : pair[3]),
            race.nodesVisited(), label);
      }
    }
  }

  @Test
  void reachStampsClearEveryMarkWhenTheyWrap() throws Exception {
    MolGraph molecule = rawGraph("CCO");
    CommonBondGraphs.Clock clock = new CommonBondGraphs.Clock(() -> false);
    CommonBondGraphs graph = new CommonBondGraphs(
        molecule, molecule, new ChemOptions(), Collections.emptySet(), clock);
    CommonBondSearch.Labels labels = new CommonBondSearch.Labels(
        graph.query, new int[2 * graph.query.bondCount], 1, 1, 1, 1, new boolean[] {true});
    // After 2^32 - 1 reaches the stamp is -1, as the unsigned C++ stamp is 2^32 - 1.
    labels.stamp = -1;
    Arrays.fill(labels.reached, 1);
    labels.reached[0] = 0;
    assertEquals(1, labels.nextStamp());
    for (int mark : labels.reached) assertEquals(0, mark);
    assertEquals(2, labels.nextStamp());
  }
}
