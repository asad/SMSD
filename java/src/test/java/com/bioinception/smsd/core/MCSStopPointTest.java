/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
package com.bioinception.smsd.core;

import static com.bioinception.smsd.TestSupport.rawGraph;
import static org.junit.jupiter.api.Assertions.*;

import java.util.*;
import java.util.function.BooleanSupplier;
import java.util.function.Function;
import org.junit.jupiter.api.Tag;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.Arguments;
import org.junit.jupiter.params.provider.MethodSource;

/**
 * MCS searches stopped by their cancellation after 1, 2, 4, 8, ... polls. A search runs on the
 * calling thread and polls its cancellation at points that the molecules and options alone fix:
 * at each node and branch of the general search, and once every so much work in the common-atom
 * and common-bond searches, the race of the other orientation included. Every run therefore makes
 * the same polls up to the one that stops it, and a time limit far beyond the search never
 * decides where it ends. Each stopped search is cancelled and unproved, and keeps valid mappings
 * of one score, which it reports. Its best score and its nodes never fall as it runs longer, and
 * by its last poll it has its best mapping; the first run that is not stopped is the search
 * without cancellation.
 */
class MCSStopPointTest {
  /** A time limit that no search here comes near, so that only the cancellation stops it. */
  private static final long UNREACHED_TIMEOUT_MILLIS = 600_000;
  /** The fewest polls of a full search here, so that a sweep stops it at five points or more. */
  private static final long FEWEST_POLLS = 20;
  /** The most mappings that matchAll keeps here, its default limit. */
  private static final int MATCH_ALL_LIMIT = 1000;

  private static final String NAPHTHALENE = "c1ccc2ccccc2c1", BIPHENYL = "c1ccc(cc1)-c1ccccc1";
  private static final String ANTHRACENE = "c1ccc2cc3ccccc3cc2c1";
  private static final String PHENANTHRENE = "c1ccc2c(c1)ccc1ccccc12";
  private static final String CAFFEINE = "Cn1cnc2c1c(=O)n(C)c(=O)n2C";
  private static final String THEOBROMINE = "Cn1cnc2c1c(=O)[nH]c(=O)n2C";
  private static final String THEOPHYLLINE = "Cn1c(=O)c2[nH]cnc2n(C)c1=O";
  private static final String URIC_ACID = "O=c1[nH]c(=O)c2[nH]c(=O)[nH]c2[nH]1";
  /** Gly-Ala-Ser-Val and Ala-Gly-Ser-Leu. */
  private static final String GASV = "NCC(=O)NC(C)C(=O)NC(CO)C(=O)NC(C(C)C)C(=O)O";
  private static final String AGSL = "NC(C)C(=O)NCC(=O)NC(CO)C(=O)NC(CC(C)C)C(=O)O";

  /** A cancellation that answers false to its first {@code allowed} polls and true to the rest. */
  private static final class StopAfter implements BooleanSupplier {
    private final long allowed;
    private long polls;

    StopAfter(long allowed) {
      this.allowed = allowed;
    }

    @Override
    public boolean getAsBoolean() {
      return ++polls > allowed;
    }
  }

  /** A cancellation that never cancels, but interrupts the calling thread on poll {@code at}. */
  private static final class InterruptAt implements BooleanSupplier {
    private final long at;
    private long polls;

    InterruptAt(long at) {
      this.at = at;
    }

    @Override
    public boolean getAsBoolean() {
      if (++polls == at)
        Thread.currentThread().interrupt();
      return false;
    }
  }

  /**
   * The searches swept: for each objective and connectivity, pairs whose full search polls at
   * least {@link #FEWEST_POLLS} times and has its best mapping by its last poll. A single MCIS in
   * fragments goes to the common-atom search, any other MCIS to the general search, and every MCES
   * to the common-bond search, connected or in fragments; matchAll enumerates every best mapping.
   * The general search for one mapping polls no more once it has mapped the whole of the smaller
   * molecule, so pairs such as caffeine and theobromine, which find that mapping after their last
   * poll, are not swept.
   */
  static List<Arguments> searches() {
    return Arrays.asList(
        Arguments.of(NAPHTHALENE, BIPHENYL, MCS.Objective.MCIS, MCS.Connectivity.Connected, false),
        Arguments.of(CAFFEINE, URIC_ACID, MCS.Objective.MCIS, MCS.Connectivity.Connected, false),
        Arguments.of(GASV, AGSL, MCS.Objective.MCIS, MCS.Connectivity.Connected, false),
        Arguments.of(NAPHTHALENE, BIPHENYL, MCS.Objective.MCIS, MCS.Connectivity.Connected, true),
        Arguments.of(ANTHRACENE, PHENANTHRENE, MCS.Objective.MCIS, MCS.Connectivity.Fragments, false),
        Arguments.of(THEOBROMINE, THEOPHYLLINE, MCS.Objective.MCIS, MCS.Connectivity.Fragments, false),
        Arguments.of(GASV, AGSL, MCS.Objective.MCIS, MCS.Connectivity.Fragments, false),
        Arguments.of(NAPHTHALENE, BIPHENYL, MCS.Objective.MCIS, MCS.Connectivity.Fragments, true),
        Arguments.of(THEOBROMINE, THEOPHYLLINE, MCS.Objective.MCES, MCS.Connectivity.Connected, false),
        Arguments.of(GASV, AGSL, MCS.Objective.MCES, MCS.Connectivity.Connected, false),
        Arguments.of(NAPHTHALENE, BIPHENYL, MCS.Objective.MCES, MCS.Connectivity.Connected, true),
        Arguments.of(THEOBROMINE, THEOPHYLLINE, MCS.Objective.MCES, MCS.Connectivity.Fragments, false),
        Arguments.of(GASV, AGSL, MCS.Objective.MCES, MCS.Connectivity.Fragments, false),
        Arguments.of(NAPHTHALENE, BIPHENYL, MCS.Objective.MCES, MCS.Connectivity.Fragments, true));
  }

  private static MCS.Options options(MCS.Objective objective, MCS.Connectivity connectivity) {
    MCS.Options options = new MCS.Options();
    options.objective = objective;
    options.connectivity = connectivity;
    options.timeoutMillis = UNREACHED_TIMEOUT_MILLIS;
    return options;
  }

  /** {@link MCS#match}, or {@link MCS#matchAll} with its default limit, under the options given. */
  private static Function<MCS.Options, SearchResult> search(MolGraph query, MolGraph target, boolean all) {
    return options -> all ? MCS.find(query).withOptions(options).matchAll(target, MATCH_ALL_LIMIT)
                          : MCS.find(query).withOptions(options).match(target);
  }

  /** The most mappings that {@link #search} keeps. */
  private static int most(boolean all) {
    return all ? MATCH_ALL_LIMIT : 1;
  }

  /**
   * The score of a mapping as the objective ranks it, atoms then bonds for an MCIS and bonds then
   * atoms for an MCES.
   */
  private static long score(SearchResult.Match match, MCS.Objective objective) {
    boolean edges = objective == MCS.Objective.MCES;
    long first = edges ? match.bondCount() : match.atomCount();
    long second = edges ? match.atomCount() : match.bondCount();
    return first << 32 | second;
  }

  /** The best score of the matches, or 0 when there is none. */
  private static long bestScore(SearchResult result, MCS.Objective objective) {
    long best = 0;
    for (SearchResult.Match match : result.matches())
      best = Math.max(best, score(match, objective));
    return best;
  }

  /** A score of {@link #bestScore} as its two counts, the one ranked first first. */
  private static String describe(long score) {
    return (score >>> 32) + " then " + (score & 0xFFFFFFFFL);
  }

  /**
   * At most {@code most} mappings, each valid under the options and all of one score, which the
   * result reports as their atoms for an MCIS and their bonds for an MCES, or 0 without one.
   */
  private static void assertMappings(MolGraph query, MolGraph target, SearchResult result,
      MCS.Options options, int most, String label) {
    assertTrue(result.matches().size() <= most,
        label + ": " + result.matches().size() + " mappings");
    long best = bestScore(result, options.objective);
    for (SearchResult.Match match : result.matches()) {
      assertEquals(Collections.emptyList(), MCS.validateMapping(query, target, match, options), label);
      assertEquals(describe(best), describe(score(match, options.objective)), label);
    }
    assertEquals((double) (best >>> 32), result.score(), label);
  }

  /** The same matches in the same order, status, proof, score and nodes; the time may differ. */
  private static void assertSameSearch(SearchResult expected, SearchResult actual, String label) {
    assertEquals(expected.matches(), actual.matches(), label);
    assertEquals(expected.status(), actual.status(), label);
    assertEquals(expected.optimal(), actual.optimal(), label);
    assertEquals(expected.complete(), actual.complete(), label);
    assertEquals(expected.score(), actual.score(), label);
    assertEquals(expected.nodesVisited(), actual.nodesVisited(), label);
  }

  /**
   * The polls of the full search, counted by a cancellation that never cancels, whose result
   * must be that of the search without cancellation, {@code plain}.
   */
  private static long polls(
      Function<MCS.Options, SearchResult> search, MCS.Options options, SearchResult plain, String label) {
    StopAfter counting = new StopAfter(Long.MAX_VALUE);
    options.cancellation = counting;
    assertSameSearch(plain, search.apply(options), label);
    return counting.polls;
  }

  /**
   * Stops {@code search}, which keeps at most {@code most} mappings, after 1, 2, 4, ... polls of
   * its cancellation, or after each number of polls with {@code everyPoint}, and after all but its
   * last poll, then runs it to its end, and checks each result against the search without
   * cancellation.
   *
   * @return the polls of the full search
   */
  private static long sweep(String label, MolGraph query, MolGraph target, MCS.Options options,
      int most, Function<MCS.Options, SearchResult> search, boolean everyPoint) {
    options.cancellation = null;
    SearchResult plain = search.apply(options);
    assertEquals(SearchResult.Status.Completed, plain.status(), label);
    assertTrue(plain.optimal() && plain.complete() && plain.found(), label);
    assertMappings(query, target, plain, options, most, label);
    long best = bestScore(plain, options.objective);
    long polls = polls(search, options, plain, label);
    assertTrue(polls >= FEWEST_POLLS, label + " polls " + polls);
    List<Long> stops = new ArrayList<>();
    long allowed = 1;
    for (; allowed < polls; allowed = everyPoint ? allowed + 1 : 2 * allowed)
      stops.add(allowed);
    if (!stops.contains(polls - 1))
      stops.add(polls - 1);
    long score = 0, nodes = 0;
    for (long stop : stops) {
      StopAfter cancellation = new StopAfter(stop);
      options.cancellation = cancellation;
      SearchResult stopped = search.apply(options);
      String at = label + " stopped after " + stop + " of " + polls + " polls";
      // The poll that answers true stops the search, which polls no more.
      assertEquals(stop + 1, cancellation.polls, at);
      assertEquals(SearchResult.Status.Cancelled, stopped.status(), at);
      assertFalse(stopped.optimal(), at);
      assertFalse(stopped.complete(), at);
      assertMappings(query, target, stopped, options, most, at);
      long reached = bestScore(stopped, options.objective);
      assertTrue(reached >= score, at + ": the best score fell from " + describe(score) + " to "
              + describe(reached));
      score = reached;
      // A stopped search is the start of the full one, so its nodes never fall as it runs longer
      // and never pass those of the full search.
      assertTrue(stopped.nodesVisited() >= nodes && stopped.nodesVisited() <= plain.nodesVisited(),
          at + ": " + stopped.nodesVisited() + " nodes after " + nodes + ", of "
              + plain.nodesVisited());
      nodes = stopped.nodesVisited();
      // By its last poll every search here has found its best mapping, which it keeps.
      if (stop == polls - 1)
        assertEquals(describe(best), describe(reached), at);
    }
    StopAfter cancellation = new StopAfter(allowed);
    options.cancellation = cancellation;
    SearchResult full = search.apply(options);
    String end = label + " allowed " + allowed + " of " + polls + " polls";
    assertEquals(polls, cancellation.polls, end);
    assertSameSearch(plain, full, end);
    return polls;
  }

  private static String label(String query, String target, MCS.Objective objective,
      MCS.Connectivity connectivity, boolean all) {
    return (all ? "matchAll " : "match ") + objective + " " + connectivity + " " + query + " " + target;
  }

  /**
   * Stopped after 1, 2, 4, ... polls, and after all but the last, each search is cancelled and
   * unproved, with valid mappings of one score, which never falls and is the best by the last
   * poll; allowed every poll it needs, it is the search without cancellation.
   */
  @ParameterizedTest
  @MethodSource("searches")
  void stopPointsKeepValidMappingsThatNeverWorsen(String querySmiles, String targetSmiles,
      MCS.Objective objective, MCS.Connectivity connectivity, boolean all) throws Exception {
    MolGraph query = rawGraph(querySmiles), target = rawGraph(targetSmiles);
    sweep(label(querySmiles, targetSmiles, objective, connectivity, all), query, target,
        options(objective, connectivity), most(all), search(query, target, all), false);
  }

  /** The same checks with the search stopped after each number of polls in turn. */
  @ParameterizedTest
  @MethodSource("searches")
  @Tag("SlowTest")
  void everyStopPointKeepsValidMappingsThatNeverWorsen(String querySmiles, String targetSmiles,
      MCS.Objective objective, MCS.Connectivity connectivity, boolean all) throws Exception {
    MolGraph query = rawGraph(querySmiles), target = rawGraph(targetSmiles);
    sweep(label(querySmiles, targetSmiles, objective, connectivity, all), query, target,
        options(objective, connectivity), most(all), search(query, target, all), true);
  }

  /**
   * A fragment search for one mapping, of any bond order, that asks for the seed grown the other
   * way round once its levels of common bonds run out of budget, and then races the search the
   * other way round against its last search. Through match the race starts only after 2^19 pairs
   * of the last search; with the lower limits that the common-bond search takes from tests, it
   * starts after 64 pairs and maps 16, and up to 64 while it beats the last search, on the same
   * thread and clock.
   */
  private static void sweepAroundTheRace(boolean everyPoint) throws Exception {
    String first = "COc1cc2c(cc1)C(=O)CC(c1ccc(C)cc1)O2";
    String second = "NS(c1cc(NC(COC(C=Cc2nc3ccccc3s2)=O)=O)ccc1)(=O)=O";
    MolGraph query = rawGraph(first), target = rawGraph(second);
    MCS.Options options = options(MCS.Objective.MCES, MCS.Connectivity.Fragments);
    options.chemistry.matchBondOrder = ChemOptions.BondOrderMode.ANY;
    CommonBondSearch.ScoutLimits raced =
        new CommonBondSearch.ScoutLimits(1L << 15, 4, 1L << 15, 1L << 10, 64, 16, 64);
    sweep("race " + first + " " + second, query, target, options, 1,
        swept -> CommonBondSearch.search(query, target, swept, 1, false, null,
            new SearchResult.Control(swept.timeoutMillis, swept.cancellation), raced),
        everyPoint);
  }

  /** Stopped before the seed grown the other way round, and after the race. */
  @Test
  void stopPointsAroundTheRaceKeepValidMappingsThatNeverWorsen() throws Exception {
    sweepAroundTheRace(false);
  }

  /** Stopped at each poll in turn, those of the race among them. */
  @Test
  @Tag("SlowTest")
  void everyStopPointAroundTheRaceKeepsValidMappingsThatNeverWorsen() throws Exception {
    sweepAroundTheRace(true);
  }

  /**
   * Longer fragment searches through match: the chain pair starts its last search from the seed
   * grown the other way round, and the fused pair is proved by the search the other way round,
   * raced against its last search after 2^19 pairs, so that it ends inside the race.
   */
  @Test
  @Tag("SlowTest")
  void stopPointsOfLongFragmentSearchesKeepValidMappingsThatNeverWorsen() throws Exception {
    String[][] pairs = {
        {"CCCCCCCCCCCCCCSCC(COP(O)(OP(O)(OCC1OC(n2ccc(N)nc2=O)C(O)C1O)=O)=O)OC(=O)CCCCCCCCCCC",
            "CN(CCOC(c1ccccc1)c1ccccc1)CCN(C)CCc1ccc(F)cc1"},
        {"CCC1(O)N2CCC3c4c(cccc4)NC3C(c3c(OC)cc4c(c3)C35CCN6CC=CC(CC)"
                + "(C(OC(CC7CCCN7)=O)C(O)(C(OC)=O)C3N4C)C56)C(C(OC)=O)C(CC1)C2",
            "Oc1nc(NCc2cc3[nH]ccc3cc2)nc(Nc2cc(Oc3ccccc3)ccc2)n1"}};
    for (String[] pair : pairs) {
      MolGraph query = rawGraph(pair[0]), target = rawGraph(pair[1]);
      MCS.Options options = options(MCS.Objective.MCES, MCS.Connectivity.Fragments);
      options.chemistry.matchBondOrder = ChemOptions.BondOrderMode.ANY;
      options.chemistry.aromaticityMode = ChemOptions.AromaticityMode.FLEXIBLE;
      sweep("long " + pair[0], query, target, options, most(false), search(query, target, false),
          false);
    }
  }

  /**
   * A search called on an interrupted thread is cancelled before its first node, under each
   * objective and connectivity, and leaves the interrupt to its caller.
   */
  @Test
  void anInterruptedCallerGetsACancelledSearch() throws Exception {
    MolGraph query = rawGraph(NAPHTHALENE), target = rawGraph(BIPHENYL);
    for (MCS.Objective objective : MCS.Objective.values())
      for (MCS.Connectivity connectivity : MCS.Connectivity.values())
        for (boolean all : new boolean[] {false, true}) {
          String label = label(NAPHTHALENE, BIPHENYL, objective, connectivity, all);
          SearchResult result;
          Thread.currentThread().interrupt();
          try {
            result = search(query, target, all).apply(options(objective, connectivity));
            assertTrue(Thread.currentThread().isInterrupted(), label);
          } finally {
            Thread.interrupted();
          }
          assertEquals(SearchResult.Status.Cancelled, result.status(), label);
          assertFalse(result.optimal(), label);
          assertFalse(result.complete(), label);
          assertFalse(result.found(), label);
          assertEquals(0.0, result.score(), label);
          assertEquals(0, result.nodesVisited(), label);
        }
  }

  /**
   * An interrupt raised during a poll stops the search at the next poll, with the result of a
   * cancellation that answers true there, and stays for the caller.
   */
  @ParameterizedTest
  @MethodSource("searches")
  void anInterruptDuringTheSearchStopsItAtTheNextPoll(String querySmiles, String targetSmiles,
      MCS.Objective objective, MCS.Connectivity connectivity, boolean all) throws Exception {
    MolGraph query = rawGraph(querySmiles), target = rawGraph(targetSmiles);
    String label = label(querySmiles, targetSmiles, objective, connectivity, all);
    MCS.Options options = options(objective, connectivity);
    Function<MCS.Options, SearchResult> search = search(query, target, all);
    SearchResult plain = search.apply(options);
    long polls = polls(search, options, plain, label);
    for (long at : new long[] {1, polls / 2, polls - 1}) {
      String stop = label + " interrupted on poll " + at + " of " + polls;
      options.cancellation = new StopAfter(at);
      SearchResult cancelled = search.apply(options);
      options.cancellation = new InterruptAt(at);
      SearchResult interrupted;
      try {
        interrupted = search.apply(options);
        assertTrue(Thread.currentThread().isInterrupted(), stop);
      } finally {
        Thread.interrupted();
      }
      assertEquals(SearchResult.Status.Cancelled, interrupted.status(), stop);
      assertSameSearch(cancelled, interrupted, stop);
    }
  }
}
