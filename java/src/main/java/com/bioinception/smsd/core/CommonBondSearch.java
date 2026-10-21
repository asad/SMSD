/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
package com.bioinception.smsd.core;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.function.Supplier;

/**
 * The maximum common edge subgraph (MCES) search on the int graphs of a
 * {@link CommonBondGraphs}: the mapping from the query to the target with the most common bonds
 * and then the most atoms, or every distinct such mapping up to a limit. Each mapped atom has a
 * common bond, and every compatible bond between mapped atoms is common. The common bonds form
 * one fragment or, unless the search is connected, several fragments anywhere in either
 * molecule. Where the options also ask for a required mapping, fragment sizes or counts,
 * stereochemistry or whole rings, each mapping recorded is checked by
 * {@link CommonBondConditions}, under which it may select only some of its common bonds; the
 * common bonds still bound every selection, so the bounds below hold. Every pair mapped counts
 * towards the node limit.
 *
 * <p>When the whole query embeds in the target, an embedding is the answer ({@link Embedding}).
 * Otherwise the query atoms are ranked, fewest candidates first, greedy seeds give a score to
 * beat, and a depth-first branch and bound tries each query atom as a root, highest rank first,
 * and grows the mapping pair by pair on an explicit stack of {@link Level}s, one for each mapped
 * pair. After the first pair of a fragment, each pair maps a query atom next to the fragment onto
 * a free neighbour of the image of one of its mapped neighbours. A level holds its pair, what
 * taking it off restores, the candidates and targets it has left to try, and a bound from the
 * pairs its fragment reaches ({@link Reach}) and from the bond labels of the open atoms
 * ({@link EdgeBound}). A level that cannot beat the best score is taken off, and the pair it tried
 * is excluded below its later siblings, so each mapping is found once. Unless the search is
 * connected, the best score is first sought level by level of common bonds, down from an upper
 * bound, within a small budget; one search from the best mapping known then finishes the work,
 * now and then scouting the root pairs ahead for a better mapping ({@link #scout}). Every loop
 * ticks the {@link CommonBondGraphs.Clock}, which throws {@link CommonBondGraphs.Stop} when the
 * time limit is reached or the search is cancelled; the best mapping found by then is kept
 * ({@link #found}).
 */
final class CommonBondSearch {

  /** The most greedy seeds grown. */
  private static final int SEEDS = 64;
  /** The most free targets a seed looks at next to a hub. */
  private static final int HUB_TARGETS = 64;
  /** The levels of likeness of two atoms, from no radius alike to all of them. */
  private static final int LIKENESS_LEVELS = CommonBondGraphs.RADII + 1;
  /** The first capacity of a list of pairs that grows as needed. */
  private static final int INITIAL_PAIRS = 64;
  /** The fewest common bonds of a fragment. */
  private static final int MIN_BONDS = 1;
  /** The most pairs the levels of common bonds map between them, see {@link #collectByBonds}. */
  private static final long LEVEL_BUDGET = 1 << 13;
  /** The most counts in a histogram of the bond-label bound, one for each molecule. */
  private static final long MAX_HISTOGRAM = 1L << 24;
  /** The most atom pairs whose root bounds are kept, see {@link #rootRefuted}. */
  private static final long MAX_ROOT_BOUNDS = 1L << 20;
  /** An empty list of atoms, shared by every level that has nothing to try. */
  private static final int[] NONE = new int[0];

  /** The graphs of both molecules, with their compatible pairs and classes. */
  private final CommonBondGraphs graph;
  /** The query and target graphs of {@link #graph}. */
  private final CommonBondGraphs.Side query, target;
  /** Counts the work of every loop and ends the search with {@link CommonBondGraphs.Stop} when told to. */
  private final CommonBondGraphs.Clock clock;
  /** Whether the common bonds form a single fragment. */
  private final boolean connected;
  /** The numbers of query and target atoms, and of longs in a row of target atoms. */
  private final int queryCount, targetCount, words;
  /** The compatible pairs that are not excluded: a row for each query atom, with a bit for each target atom. */
  private final long[] allowed;
  /** The mapping both ways, -1 where unmapped. */
  private final int[] queryToTarget, targetToQuery;
  /** The mapped pairs in the order they were mapped, the first {@link #depth} of the array. */
  private final Level[] levels;
  /** The number of mapped pairs, and the common bonds between them. */
  private int depth, commonBonds;
  /** The pairs excluded so far, the first {@link #excludedCount}, put back as the levels that tried them go. */
  private int[] excluded = new int[INITIAL_PAIRS];
  /** The number of excluded pairs. */
  private int excludedCount;
  /** The atoms by rank, lowest first, and the rank of each atom, see {@link #rank}. */
  private int[] queryOrder, queryRank, targetOrder, targetRank;
  /** The cell of each atom: atoms that look alike share a cell. */
  private int[] queryCell, targetCell;
  /** The signatures of the atoms at each radius, which tell how alike two atoms are. */
  private long[][] querySignature, targetSignature;
  /** The next atom alike to each terminal atom, see {@link #twins}. */
  private int[] queryTwin, targetTwin;
  /**
   * The taken targets that lead the neighbours of each target atom, and the seed they were counted
   * for, see {@link #leadingTaken}.
   */
  private final int[] leadingTaken, leadingRound;
  /** The likest level of each query atom with any target, -1 until worked out, see {@link #likestLevel}. */
  private final int[] likestLevel;
  /** The best score so far, see {@link CommonBondGraphs#score}, and a score that no mapping beats. */
  private long bestScore, upperBound;
  /** The mappings held: the embeddings found, or those with the best score so far, see {@link #record}. */
  private List<int[]> mappings;
  /** The most mappings wanted. */
  private int maxMappings;
  /** The seed, held at the front of the mappings from the start until the search finds it again. */
  private int[] heldSeed;
  /**
   * The mapping a stopped search falls back on while the pass under way holds none: the seed,
   * or the best mapping a level of the search met below its threshold; null until the seed is
   * grown. See {@link #found}.
   */
  private int[] fallback;
  /** The score of {@link #fallback}, or one below every mapping. */
  private long fallbackScore;
  /** The pairs reached from a state, for its bound. */
  private final Reach reach;
  /** The bound from the bond labels of the open atoms. */
  private final EdgeBound edgeBound;
  /** The pairs mapped and the embedding atoms placed so far. */
  private long nodes;
  /** The number of {@link #nodes} at which the pass under way ends, see {@link #done}. */
  private long nodeLimit = Long.MAX_VALUE;
  /** The last bound worked out, before it was compared with the best score, see {@link #bound}. */
  private long lastBound;
  /**
   * The most common bonds below each root pair, by {@link #pair}, as the bound of its state
   * showed when it was last pruned, or -1. Null when the pairs are too many to keep.
   */
  private final int[] rootBonds;
  /** When the last search scouts and how far, see {@link #scout}. */
  private final ScoutLimits scoutLimits;
  /** The number of {@link #nodes} at which the search under way scouts, see {@link #scout}, and the gap to the next scouting. */
  private long scoutAt = Long.MAX_VALUE, scoutGap;
  /** The scoutings so far, and the root pairs skipped as refuted, see {@link #rootRefuted}. */
  private long scoutings, refutedRoots;
  /** The root pair of the search under way, as the rank of its query atom and the place of its target. */
  private int rootRank, rootPlace;
  /** What a mapping must satisfy beyond its common bonds, or null when nothing else is asked. */
  private final CommonBondConditions conditions;
  /** The most mappings the search holds, see {@link #run}. */
  private final int wanted;
  /** The most pairs mapped and embedding atoms placed before the search stops, see {@link #map}. */
  private final long nodeCap;
  /** Whether the search stopped at {@link #nodeCap}. */
  private boolean capReached;
  /** The lowest rank of a first root: with a required mapping, only that of its highest ranked atom. */
  private int lowestRoot;
  /** The target that a required mapping fixes for each query atom, -1 for none; null without one. */
  private final int[] requiredTarget;
  /** The selections of the mapping being recorded, and of one offered as the fallback. */
  private final List<int[]> selections = new ArrayList<>(), offered = new ArrayList<>();
  /** The mappings held, while several are wanted, see {@link #keep}. */
  private final Set<Record> held = new HashSet<>();
  /** The number of fragments of the mapping, the open one included, see {@link #nextRoot}. */
  private int fragments;
  /** Grows a seed of the same molecules the other way round, turned round, or null; see {@link #seedAlsoFrom}. */
  private Supplier<int[]> otherSeed;
  /** Scratch for sorting candidates and targets, used up before each returns, see {@link #candidates}. */
  private int[] layerScratch = new int[INITIAL_PAIRS];
  private long[] keyScratch = new long[INITIAL_PAIRS];

  /**
   * When the last search first scouts, see {@link #scout}, and how far each scouting goes.
   * Searches use {@link #DEFAULT}; tests lower the limits so that small searches also scout and
   * skip refuted root pairs.
   */
  static final class ScoutLimits {
    /** The limits of every search. */
    static final ScoutLimits DEFAULT = new ScoutLimits(1L << 15, 4, 1L << 15, 1L << 10);
    /** The pairs the last search maps before it first scouts, and the factor by which the gap to the next scouting grows. */
    final long after, growth;
    /** The most pairs mapped in one scouting, and below one root pair. */
    final long budget, root;

    ScoutLimits(long after, long growth, long budget, long root) {
      this.after = after;
      this.growth = growth;
      this.budget = budget;
      this.root = root;
    }
  }

  /**
   * A search of {@code graph}. With {@code connected} the common bonds form one fragment,
   * otherwise any number of fragments of at least one bond.
   */
  CommonBondSearch(CommonBondGraphs graph, CommonBondGraphs.Clock clock, boolean connected) {
    this(graph, clock, connected, ScoutLimits.DEFAULT);
  }

  /** A search of {@code graph} that scouts as {@code limits} say; only tests pass other limits than the default. */
  CommonBondSearch(CommonBondGraphs graph, CommonBondGraphs.Clock clock, boolean connected,
      ScoutLimits limits) {
    this(graph, clock, connected, limits, null, null, 1, Long.MAX_VALUE);
  }

  /**
   * A search of {@code graph} for at most {@code wanted} mappings that meet {@code conditions},
   * unless they are null, that keep the target of each query atom in {@code requiredTarget}, -1
   * for none, unless it is null, and that stops once it would map more than {@code nodeCap}
   * pairs.
   */
  CommonBondSearch(CommonBondGraphs graph, CommonBondGraphs.Clock clock, boolean connected,
      ScoutLimits limits, CommonBondConditions conditions, int[] requiredTarget, int wanted,
      long nodeCap) {
    this.graph = graph;
    this.clock = clock;
    this.connected = connected;
    this.conditions = conditions;
    this.requiredTarget = requiredTarget;
    this.wanted = wanted;
    this.nodeCap = nodeCap;
    scoutLimits = limits;
    scoutGap = limits.after;
    query = graph.query;
    target = graph.target;
    queryCount = query.atomCount;
    targetCount = target.atomCount;
    words = CommonBondGraphs.words(targetCount);
    allowed = graph.compatiblePairs();
    if (requiredTarget != null)
      keepRequiredPairs();
    queryToTarget = new int[queryCount];
    targetToQuery = new int[targetCount];
    Arrays.fill(queryToTarget, -1);
    Arrays.fill(targetToQuery, -1);
    levels = new Level[Math.min(queryCount, targetCount)];
    reach = new Reach();
    edgeBound = new EdgeBound();
    leadingTaken = new int[targetCount];
    leadingRound = new int[targetCount];
    likestLevel = new int[queryCount];
    Arrays.fill(likestLevel, -1);
    if ((long) queryCount * targetCount <= MAX_ROOT_BOUNDS) {
      rootBonds = new int[queryCount * targetCount];
      Arrays.fill(rootBonds, -1);
    } else {
      rootBonds = null;
    }
  }

  /**
   * Searches {@code query} in {@code target} for the best MCES under {@code options}: one
   * mapping, or with {@code all} every distinct best mapping up to {@code limit}, each
   * extending {@code required} unless it is null. The matches are sorted by atom mapping and
   * then by bond mapping. A search that holds more than {@code limit} best mappings keeps the
   * first {@code limit} it found and reports the limit reached. A search stopped by the time
   * limit, the cancellation or the node limit keeps what it found by then and is not optimal.
   * Graphs this search cannot take, as they are not well formed or too large (see
   * {@link CommonBondGraphs#supports} and {@link #fits}), give no match, with the limit reached.
   */
  static SearchResult search(MolGraph query, MolGraph target, MCS.Options options, int limit,
      boolean all, SearchResult.Match required, SearchResult.Control control) {
    boolean connected =
        options.connectivity == MCS.Connectivity.Connected || options.maxFragments == 1;
    if (control.stopped())
      return new SearchResult(Collections.emptyList(), control.status(), false, false, 0.0, 0,
          control.elapsedMillis());
    CommonBondGraphs.Clock clock = new CommonBondGraphs.Clock(control::stopped);
    CommonBondGraphs graph = null;
    CommonBondSearch search = null;
    List<int[]> found = Collections.emptyList();
    boolean declined = !CommonBondGraphs.supports(query, target);
    try {
      clock.poll();
      if (!declined) {
        graph = new CommonBondGraphs(query, target, pairedChemistry(options.chemistry),
            options.excludedTargetAtoms, clock);
        declined = !fits(graph);
      }
      if (!declined) {
        CommonBondConditions conditions =
            CommonBondConditions.of(query, target, graph, options, connected, required, clock);
        int wanted = all ? (limit == Integer.MAX_VALUE ? limit : limit + 1) : 1;
        int[] requiredTarget = requiredTargets(required, query.n);
        search = new CommonBondSearch(graph, clock, connected, ScoutLimits.DEFAULT, conditions,
            requiredTarget, wanted, options.nodeLimit == 0 ? Long.MAX_VALUE : options.nodeLimit);
        // Without a required mapping or excluded atoms, the search the other way round is of the
        // same mappings, so its seed, turned round, is a mapping of this search. A fragment search
        // grows it only once its levels run out of budget, see seedAlsoFrom().
        if (!connected && requiredTarget == null && options.excludedTargetAtoms.isEmpty()
            && CommonBondGraphs.supports(target, query)) {
          ChemOptions chemistry = pairedChemistry(options.chemistry);
          search.seedAlsoFrom(() -> {
            CommonBondGraphs reversed = new CommonBondGraphs(target, query, chemistry,
                options.excludedTargetAtoms, clock);
            if (!fits(reversed))
              return null;
            int[] seed = new CommonBondSearch(reversed, clock, connected).greedySeed();
            if (seed == null)
              return null;
            int[] turned = new int[query.n];
            Arrays.fill(turned, -1);
            for (int targetAtom = 0; targetAtom < target.n; targetAtom++)
              if (seed[targetAtom] >= 0)
                turned[seed[targetAtom]] = targetAtom;
            return turned;
          });
        }
        found = search.run();
      }
    } catch (CommonBondGraphs.Stop stop) {
      // The time ran out, the search was cancelled or it reached its node limit. What the search
      // found by then is the result.
      if (search != null)
        found = search.found();
    }
    List<SearchResult.Match> matches = new ArrayList<>();
    for (int[] record : found) {
      if (matches.size() == limit)
        break;
      matches.add(match(graph, record));
    }
    if (declined || search != null && search.capReached)
      control.limitReached();
    boolean proved = control.status() == SearchResult.Status.Completed;
    boolean overflow = all && found.size() > limit;
    if (proved && overflow)
      control.limitReached();
    matches.sort(CANONICAL);
    double score = matches.isEmpty() ? 0.0 : matches.get(0).score();
    return new SearchResult(matches, control.status(), proved, proved && !overflow, score,
        search == null ? 0 : search.nodes, control.elapsedMillis());
  }

  /**
   * The chemistry the atom and bond pairs are compared under. With complete rings a mapped ring
   * atom lies on a whole ring of both molecules, and so does a selected ring bond, so ring atoms
   * and bonds pair only with ring atoms and bonds, and chain ones with chain ones.
   */
  private static ChemOptions pairedChemistry(ChemOptions chemistry) {
    if (!chemistry.completeRingsOnly || chemistry.ringMatchesRingOnly)
      return chemistry;
    ChemOptions paired = ChemOptions.copyOf(chemistry);
    paired.ringMatchesRingOnly = true;
    return paired;
  }

  /** The target that {@code required} fixes for each of {@code atoms} query atoms, -1 for none, or null. */
  private static int[] requiredTargets(SearchResult.Match required, int atoms) {
    if (required == null || required.atomMapping().isEmpty())
      return null;
    int[] targets = new int[atoms];
    Arrays.fill(targets, -1);
    required.atomMapping().forEach((queryAtom, targetAtom) -> targets[queryAtom] = targetAtom);
    return targets;
  }

  /** Orders matches by atom mapping, then by bond mapping, each compared pair by pair. */
  static final Comparator<SearchResult.Match> CANONICAL = (first, second) -> {
    int order = compareRows(atomRows(first), atomRows(second));
    return order != 0 ? order : compareRows(bondRows(first), bondRows(second));
  };

  /** The atom pairs of {@code match} in query order. */
  private static int[][] atomRows(SearchResult.Match match) {
    int[][] rows = new int[match.atomCount()][];
    int index = 0;
    for (Map.Entry<Integer, Integer> pair : match.atomMapping().entrySet())
      rows[index++] = new int[] {pair.getKey(), pair.getValue()};
    return rows;
  }

  /** The bonds of {@code match}, each as its query atoms and their images, in query order. */
  private static int[][] bondRows(SearchResult.Match match) {
    int[][] rows = new int[match.bondCount()][];
    int index = 0;
    for (SearchResult.Bond bond : match.bondMapping().keySet())
      rows[index++] = new int[] {bond.first(), bond.second(), match.atomMapping().get(bond.first()),
          match.atomMapping().get(bond.second())};
    Arrays.sort(rows, CommonBondSearch::compareRow);
    return rows;
  }

  /** Compares two lists of rows element by element, a shorter prefix first. */
  private static int compareRows(int[][] first, int[][] second) {
    for (int index = 0; index < first.length && index < second.length; index++) {
      int order = compareRow(first[index], second[index]);
      if (order != 0)
        return order;
    }
    return Integer.compare(first.length, second.length);
  }

  /** Compares two rows element by element, a shorter prefix first. */
  private static int compareRow(int[] first, int[] second) {
    for (int index = 0; index < first.length && index < second.length; index++)
      if (first[index] != second[index])
        return Integer.compare(first[index], second[index]);
    return Integer.compare(first.length, second.length);
  }

  /**
   * Keeps for each query atom of the required mapping only the pair with its target, and takes
   * these targets from every other query atom.
   */
  private void keepRequiredPairs() {
    for (int queryAtom = 0; queryAtom < queryCount; queryAtom++) {
      int targetAtom = requiredTarget[queryAtom];
      if (targetAtom < 0)
        continue;
      boolean compatible = has(allowed, queryAtom, targetAtom);
      Arrays.fill(allowed, queryAtom * words, (queryAtom + 1) * words, 0L);
      if (compatible)
        set(allowed, queryAtom, targetAtom);
      for (int other = 0; other < queryCount; other++)
        if (other != queryAtom)
          clear(allowed, other, targetAtom);
    }
  }

  /**
   * Whether the histograms of the bond-label bound of {@code graph} stay within
   * {@link #MAX_HISTOGRAM} counts: a row of the most bonds at an atom, plus one, for each end
   * of each label. A search that needs more is declined.
   */
  static boolean fits(CommonBondGraphs graph) {
    return 2L * graph.labels * histogramWidth(graph) <= MAX_HISTOGRAM;
  }

  /** The width of a histogram row: the most bonds at any atom of either molecule, plus one. */
  private static int histogramWidth(CommonBondGraphs graph) {
    return Math.max(graph.query.maxDegree(), graph.target.maxDegree()) + 1;
  }

  /**
   * The match of {@code record}: its atoms, and its selected common bonds in the order of the
   * query bonds, each with its image, and the fragments these bonds form. A record is a mapping,
   * followed, when it leaves common bonds out, by their number and the query bonds left out in
   * ascending order, see {@link CommonBondConditions#select}. The score is the number of
   * selected bonds.
   */
  static SearchResult.Match match(CommonBondGraphs graph, int[] record) {
    CommonBondGraphs.Side query = graph.query, target = graph.target;
    Map<Integer, Integer> atoms = new TreeMap<>();
    for (int queryAtom = 0; queryAtom < query.atomCount; queryAtom++)
      if (record[queryAtom] >= 0)
        atoms.put(queryAtom, record[queryAtom]);
    Map<SearchResult.Bond, SearchResult.Bond> bonds = new LinkedHashMap<>();
    int skipped = query.atomCount + 1, skippedEnd = skipped + leftOut(record, query.atomCount);
    for (int queryBond = 0; queryBond < query.bondCount; queryBond++) {
      int first = record[query.bondFirst[queryBond]], second = record[query.bondSecond[queryBond]];
      if (first < 0 || second < 0)
        continue;
      if (skipped < skippedEnd && record[skipped] == queryBond) {
        skipped++;
        continue;
      }
      int targetBond = target.bondBetween(first, second);
      if (targetBond >= 0 && graph.bondsMatch(queryBond, targetBond))
        bonds.put(new SearchResult.Bond(query.bondFirst[queryBond], query.bondSecond[queryBond]),
            new SearchResult.Bond(first, second));
    }
    return new SearchResult.Match(atoms, bonds, MatchGraph.fragments(atoms, bonds.keySet()), bonds.size());
  }

  /** The number of common bonds that {@code record}, of a query of {@code atoms} atoms, leaves out. */
  static int leftOut(int[] record, int atoms) {
    return record.length > atoms ? record[atoms] : 0;
  }

  /**
   * The maximum mappings, at most as many as wanted, as records (see {@link #match}), or none
   * when no pair is compatible or either molecule has no bond. When the whole query embeds, the
   * embeddings are the result. Otherwise the atoms are ranked, a greedy seed gives the score to
   * start from, and the branch and bound finds the mappings. A search that holds as many as
   * wanted may have more of the same score.
   */
  List<int[]> run() {
    clock.poll();
    if (!anyCompatible() || Math.min(query.bondCount, target.bondCount) < MIN_BONDS)
      return new ArrayList<>();
    // Embeddings that differ by an exchange of twins meet the conditions alike, so under
    // conditions only one of each kind is tried until one meets them. Then, if more are wanted,
    // every embedding is tried.
    Embedding embedding = embedding(conditions != null && conditions.keepsTwins());
    List<int[]> embeddings;
    try {
      embeddings = embedding != null ? embedding.all(wanted) : new ArrayList<>();
      if (embedding != null && embedding.prunesTwins && wanted > 1 && !embeddings.isEmpty())
        embeddings = embedding(false).all(wanted);
    } catch (CommonBondGraphs.Stop stop) {
      // At the node limit before an embedding is found, the seeds, which take no nodes, still
      // give a mapping.
      if (capReached && mappings.isEmpty()) {
        rank();
        fallbackScore = belowEveryMapping();
        seed(bestPossible());
      }
      throw stop;
    }
    if (!embeddings.isEmpty())
      return embeddings;
    rank();
    long bound = bestPossible();
    fallbackScore = belowEveryMapping();
    seed(bound);
    return collectByBonds(bound);
  }

  /** The scoutings so far. */
  long scoutings() {
    return scoutings;
  }

  /** The root pairs skipped so far as refuted. */
  long refutedRoots() {
    return refutedRoots;
  }

  /**
   * The best mappings found so far, for a search the clock stopped: the mappings held or, while
   * the pass under way holds none, the one it falls back on. Empty when nothing was found. Each
   * is complete, as a mapping is held only once it is recorded whole.
   */
  List<int[]> found() {
    if (mappings != null && !mappings.isEmpty())
      return mappings;
    if (fallback != null)
      return Collections.singletonList(fallback);
    return Collections.emptyList();
  }

  /**
   * The best greedy seed, as {@link #run} grows it before the branch and bound, or null when
   * there is none. Another search of the same molecules may start from it, see
   * {@link #seedAlsoFrom}.
   */
  int[] greedySeed() {
    clock.poll();
    if (!anyCompatible() || Math.min(query.bondCount, target.bondCount) < MIN_BONDS)
      return null;
    rank();
    fallbackScore = belowEveryMapping();
    seed(bestPossible());
    return fallback;
  }

  /**
   * Lets the search ask {@code grow} for a second seed: a mapping of the same molecules grown by
   * a search the other way round and turned round, or null. The greedy seed grows from the
   * query, so the two seeds often differ, and the one grown from the other molecule may hold
   * more common bonds. A fragment search asks only once its levels of common bonds run out of
   * budget, as the last search is then long and starts from the better seed, see
   * {@link #collectByBonds}. The seed costs no nodes.
   */
  void seedAlsoFrom(Supplier<int[]> grow) {
    otherSeed = grow;
  }

  /**
   * A maximum mapping below {@code bound}. Unless the search is connected, the best score is
   * sought first level by level of common bonds, down from the bound to the score of the best
   * mapping known: each level asks for one mapping with at least that many bonds, and without
   * one the bound comes down. Twins can then be skipped all along, see {@link #twinExcluded}. A
   * level far above the best score prunes much and costs little, so this soon finds an optimum
   * close to the bound. Each level down repeats the work of the one above it, though, so the
   * levels together map at most {@link #LEVEL_BUDGET} pairs. A level cut short refutes nothing.
   * A connected search has no levels: its bound, see {@link EdgeBound#groupedBound}, prunes as
   * much near the best mapping known, so the levels would mostly repeat work. One search from
   * the best mapping known then finishes the work, its threshold rising as better mappings are
   * found. When the levels run out of budget, that search starts from the seed grown the other
   * way round if that is better, see {@link #seedAlsoFrom}: the search is long, and the higher
   * it starts the more it prunes.
   */
  private List<int[]> collectByBonds(long bound) {
    long budget = nodes + LEVEL_BUDGET;
    for (int triedBonds = bonds(bound); !connected && triedBonds > bonds(fallbackScore);
         triedBonds--) {
      nodeLimit = budget;
      List<int[]> mapping = collect(null, CommonBondGraphs.score(0, triedBonds) - 1, bound, 1);
      nodeLimit = Long.MAX_VALUE;
      if (nodes >= budget) {
        // A mapping found in a level cut short has the bonds of the level, but perhaps not the most atoms.
        if (!mapping.isEmpty() && score(mapping.get(0)) > fallbackScore) {
          fallback = mapping.get(0);
          fallbackScore = score(fallback);
        }
        if (otherSeed != null) {
          int[] other = otherSeed.get();
          if (other != null)
            offerFallback(other, score(other));
        }
        break;
      }
      if (!mapping.isEmpty()) {
        if (wanted == 1)
          return mapping;
        // The level holds a best mapping; one search from it finds the others of its score.
        fallback = mapping.get(0);
        fallbackScore = bound = score(fallback);
        break;
      }
      // No mapping has that many bonds. With one bond fewer there are at most as many atoms, or
      // one more for each further fragment.
      bound = CommonBondGraphs.score(Math.min(atoms(bound), mostAtoms(triedBonds - 1)), triedBonds - 1);
    }
    scoutAt = nodes + scoutLimits.after;
    return collect(fallback, fallbackScore, bound, wanted);
  }

  /** The most atoms of a mapping of {@code bonds} common bonds: one fragment of b bonds has at most b + 1 atoms. */
  private int mostAtoms(int bonds) {
    if (bonds <= 0)
      return 0;
    return connected ? bonds + 1 : bonds + bonds / MIN_BONDS;
  }

  /**
   * A score that no mapping beats, when the query does not fit whole. Every cap on the common
   * bonds or on the atoms holds for every mapping, so the two are capped apart: the least of
   * the scores would keep the atoms of the one with the fewest bonds, and lose a lower cap on
   * the atoms.
   */
  private long bestPossible() {
    int bondCap = Math.min(query.bondCount, target.bondCount);
    int atomCap = Math.min(queryCount, targetCount);
    // Without an embedding no mapping of a connected query has every query bond. A query in
    // parts may still have every bond.
    if (queryConnected())
      bondCap = Math.min(bondCap, query.bondCount - 1);
    long label = graph.labelBound(connected);
    bondCap = Math.min(bondCap, CommonBondGraphs.bonds(label));
    atomCap = Math.min(atomCap, CommonBondGraphs.atoms(label));
    // A mapping with at most so many bonds has at most so many atoms.
    return CommonBondGraphs.score(Math.min(atomCap, mostAtoms(bondCap)), bondCap);
  }

  /** A score below that of every mapping, as every mapping has a common bond. */
  private long belowEveryMapping() {
    return CommonBondGraphs.score(queryCount + 1, MIN_BONDS - 1);
  }

  /** Whether any pair of atoms is compatible. */
  private boolean anyCompatible() {
    for (long word : allowed)
      if (word != 0)
        return true;
    return false;
  }

  /**
   * Ranks the atoms and works out their signatures. Query atoms with the fewest candidates, then
   * the most bonds, rank highest; ties are broken by colour refinement, so that the order is the
   * same for the same graph. The search tries the roots from the highest rank down, and below a
   * root maps only atoms of lower rank. Target atoms are ranked by their cell only, in input
   * order within it.
   */
  void rank() {
    querySignature = query.signatures(clock);
    targetSignature = target.signatures(clock);
    int[] queryChoices = new int[queryCount], targetChoices = new int[targetCount];
    for (int queryAtom = 0; queryAtom < queryCount; queryAtom++) {
      for (int targetAtom = CommonBondGraphs.nextBit(allowed, queryAtom * words, words, 0); targetAtom >= 0;
           targetAtom = CommonBondGraphs.nextBit(allowed, queryAtom * words, words, targetAtom + 1)) {
        clock.tick();
        queryChoices[queryAtom]++;
        targetChoices[targetAtom]++;
      }
    }
    int[] colours = query.invariants(queryChoices);
    queryCell = CommonBondGraphs.refine(query.start, query.neighbour, colours, false, clock);
    int[] canon = CommonBondGraphs.refine(query.start, query.neighbour, colours, true, clock);
    // Ties among target atoms stay in input order. Breaking them one atom at a time costs a
    // round per atom on a symmetric target, and the order only has to be repeatable.
    targetCell = CommonBondGraphs.refine(
        target.start, target.neighbour, target.invariants(targetChoices), false, clock);
    int[] next = new int[targetCount];
    targetRank = new int[targetCount];
    for (int targetAtom = 0; targetAtom < targetCount; targetAtom++)
      next[targetCell[targetAtom]] = targetCell[targetAtom];
    for (int targetAtom = 0; targetAtom < targetCount; targetAtom++)
      targetRank[targetAtom] = next[targetCell[targetAtom]]++;
    targetOrder = CommonBondGraphs.inverse(targetRank);
    // The most choices come first, then the fewest bonds, then the lowest canonical rank. The
    // atoms of a required mapping come last of all, so that its highest ranked atom roots every
    // mapping that holds it. Each key packs these, low to high, and as the canonical ranks differ
    // so do the keys; the canonical rank at the bottom gives the atom back.
    int canonBits = bitsFor(queryCount - 1), degreeBits = bitsFor(query.maxDegree());
    int choiceBits = bitsFor(targetCount);
    long[] keys = new long[queryCount];
    for (int queryAtom = 0; queryAtom < queryCount; queryAtom++) {
      long key = requiredTarget != null && requiredTarget[queryAtom] >= 0 ? 1 : 0;
      key = key << choiceBits | targetCount - queryChoices[queryAtom];
      key = key << degreeBits | query.degree(queryAtom);
      keys[queryAtom] = key << canonBits | canon[queryAtom];
    }
    Arrays.sort(keys);
    int[] atomOfCanon = CommonBondGraphs.inverse(canon);
    long canonMask = (1L << canonBits) - 1;
    queryOrder = new int[queryCount];
    queryRank = new int[queryCount];
    for (int rank = 0; rank < queryCount; rank++) {
      queryOrder[rank] = atomOfCanon[(int) (keys[rank] & canonMask)];
      queryRank[queryOrder[rank]] = rank;
    }
    lowestRoot = requiredTarget != null ? queryCount - 1 : 0;
  }

  /** The number of bits that hold {@code value}, which is not negative. */
  private static int bitsFor(int value) {
    return 32 - Integer.numberOfLeadingZeros(value);
  }

  /**
   * The best greedy seed, if any, as the {@link #fallback}. A seed is a connected mapping grown
   * from a root pair, from each of at most {@link #SEEDS} of them or until one reaches
   * {@code bound}. Unless the search is connected, it is a set of fragments of at least one common
   * bond, grown one after another on the free atoms. The fallback follows each better seed, so a
   * search stopped while seeding still has the best seed grown by then, or the seed being grown if
   * that is better. It needs {@link #rank} first.
   */
  void seed(long bound) {
    int[] keptBonds = new int[1];
    try {
      growSeeds(bound, keptBonds);
    } catch (CommonBondGraphs.Stop stop) {
      // Each pair added to the seed being grown brought a common bond, so it is a fragment once it
      // has one. Its conditions are not checked once the search is to stop.
      long current = score(queryToTarget);
      if (bonds(current) - keptBonds[0] >= MIN_BONDS && current > fallbackScore
          && conditions == null) {
        fallback = queryToTarget.clone();
        fallbackScore = current;
      }
      throw stop;
    }
    Arrays.fill(queryToTarget, -1);
    Arrays.fill(targetToQuery, -1);
  }

  /**
   * Makes {@code mapping}, of score {@code current}, the {@link #fallback} if it beats it. Under
   * conditions its best selection that meets them, if any, is the one that has to beat it.
   */
  private void offerFallback(int[] mapping, long current) {
    if (current <= fallbackScore)
      return;
    if (conditions == null) {
      fallback = mapping.clone();
      fallbackScore = current;
      return;
    }
    offered.clear();
    long valid = conditions.select(mapping, fallbackScore + 1, 1, offered);
    if (valid >= 0) {
      fallback = offered.get(0);
      fallbackScore = valid;
    }
  }

  /**
   * The rounds of {@link #seed}, with the common bonds of the fragments kept from earlier rounds in
   * {@code keptBonds[0]}.
   */
  private void growSeeds(long bound, int[] keptBonds) {
    long[] taken = new long[queryCount * words];
    int[] roots = new int[LIKENESS_LEVELS * SEEDS], rootCount = new int[LIKENESS_LEVELS];
    int[] frontier = new int[queryCount], merged = new int[queryCount], mapped = new int[queryCount];
    int[] added = new int[queryCount];
    long[] keys = new long[queryCount];
    long bestSeed = belowEveryMapping();
    int[] bestMapping = null;
    for (int round = 0, keptAtoms = 0;;) {
      long before = bestSeed;
      seedRoots(roots, rootCount, taken);
      for (int level = CommonBondGraphs.RADII, tried = 0; level >= 0 && bestSeed != bound;
           level--) {
        for (int rootIndex = 0;
             rootIndex < rootCount[level] && tried < SEEDS && bestSeed != bound;
             rootIndex++, tried++) {
          long grown = growSeed(roots[level * SEEDS + rootIndex], frontier, keys, merged, mapped,
              added, ++round);
          int seedAtoms = atoms(grown), seedBonds = bonds(grown);
          if (seedBonds >= MIN_BONDS
              && CommonBondGraphs.score(keptAtoms + seedAtoms, keptBonds[0] + seedBonds) > bestSeed) {
            bestSeed = CommonBondGraphs.score(keptAtoms + seedAtoms, keptBonds[0] + seedBonds);
            bestMapping = queryToTarget.clone();
            offerFallback(bestMapping, score(bestMapping));
          }
          for (int mappedIndex = 0; mappedIndex < seedAtoms; mappedIndex++) {
            targetToQuery[queryToTarget[mapped[mappedIndex]]] = -1;
            queryToTarget[mapped[mappedIndex]] = -1;
          }
        }
      }
      if (connected || bestSeed == before || bestSeed == bound)
        break;
      // The fragments found so far stay mapped while the next one grows. It shares no common
      // bond with them, since each grew until no pair next to it was left.
      keptAtoms = atoms(bestSeed);
      keptBonds[0] = bonds(bestSeed);
      System.arraycopy(bestMapping, 0, queryToTarget, 0, queryCount);
      for (int queryAtom = 0; queryAtom < queryCount; queryAtom++)
        if (queryToTarget[queryAtom] >= 0)
          targetToQuery[queryToTarget[queryAtom]] = queryAtom;
    }
  }

  /**
   * Picks the seed roots by how alike their atoms are. The {@code rootCount[level]} roots of
   * each level are at {@code roots[level * SEEDS]}. Alike cells grow alike seeds, so there is
   * one root for each pair of cells, marked in {@code taken}.
   */
  private void seedRoots(int[] roots, int[] rootCount, long[] taken) {
    Arrays.fill(taken, 0);
    Arrays.fill(rootCount, 0);
    int likest = CommonBondGraphs.RADII;
    // The most alike roots are tried first. Once there are enough of them, the rest are not needed.
    for (int rank = queryCount - 1; rank >= 0 && rootCount[likest] < SEEDS; rank--) {
      int queryAtom = queryOrder[rank];
      for (int targetPosition = 0; targetPosition < targetCount && rootCount[likest] < SEEDS
           && queryToTarget[queryAtom] < 0;
           targetPosition++) {
        int targetAtom = targetOrder[targetPosition];
        clock.tick();
        if (targetToQuery[targetAtom] >= 0 || !has(allowed, queryAtom, targetAtom))
          continue;
        int level = likeness(queryAtom, targetAtom);
        if (rootCount[level] < SEEDS && !has(taken, queryCell[queryAtom], targetCell[targetAtom])) {
          set(taken, queryCell[queryAtom], targetCell[targetAtom]);
          roots[level * SEEDS + rootCount[level]++] = pair(queryAtom, targetAtom);
        }
      }
    }
  }

  /**
   * Grows a seed from the root pair, each time by the best frontier pair, see
   * {@link #bestSeedPair}. The frontier is the free query atoms next to the seed, each listed once
   * by marking it with {@code round} in {@code added}, highest rank first; {@code keys} and
   * {@code merged} are scratch for it. The best pair is the same in any order, and in this one
   * {@link #bestSeedPair} soon passes over the atoms that cannot beat the pair it holds. The seed
   * is left mapped, with its query atoms in {@code mapped} in order. The result is its score.
   */
  private long growSeed(int root, int[] frontier, long[] keys, int[] merged, int[] mapped,
      int[] added, int round) {
    int count = 0, frontierSize = 0, seedBonds = 0;
    for (int pair = root; pair >= 0;) {
      int queryAtom = queryAtom(pair), targetAtom = targetAtom(pair);
      queryToTarget[queryAtom] = targetAtom;
      targetToQuery[targetAtom] = queryAtom;
      mapped[count++] = queryAtom;
      int listed = frontierSize;
      for (int queryPlace = query.start[queryAtom]; queryPlace < query.start[queryAtom + 1];
           queryPlace++) {
        int queryNeighbour = query.neighbour[queryPlace];
        if (queryToTarget[queryNeighbour] < 0 && added[queryNeighbour] != round) {
          added[queryNeighbour] = round;
          frontier[frontierSize++] = queryNeighbour;
        }
      }
      CommonBondGraphs.mergeByRank(frontier, listed, frontierSize, queryRank, keys, merged, clock);
      long best = -1;
      int kept = 0;
      for (int frontierIndex = 0; frontierIndex < frontierSize; frontierIndex++) {
        int frontierAtom = frontier[frontierIndex];
        // The atom just mapped leaves the frontier.
        if (queryToTarget[frontierAtom] >= 0)
          continue;
        frontier[kept++] = frontierAtom;
        best = bestSeedPair(frontierAtom, best, round);
      }
      frontierSize = kept;
      pair = best < 0 ? -1 : (int) best;
      seedBonds += best < 0 ? 0 : (int) (best >>> 32);
    }
    return CommonBondGraphs.score(count, seedBonds);
  }

  /**
   * The best pair of {@code frontierAtom} with a free target next to the image of one of its
   * mapped neighbours, through a compatible bond, unless the pair {@code held} beats it: the
   * most alike pair, then the one with the most common bonds. Ties go to the higher ranked query
   * atom, then to the lower ranked target atom. The pair is returned in the low int, with its
   * common bonds in the high int, or -1 if there is none; {@code held} is read the same way.
   * Next to a hub, the first few free targets are as good for a seed as all of them. The seed of
   * {@code round} only takes targets, so the taken targets that lead the list of an image are
   * passed once for the seed, see {@link #leadingTaken}, and not again for each frontier atom.
   */
  private long bestSeedPair(int frontierAtom, long held, int round) {
    int pair = held < 0 ? -1 : (int) held, links = held < 0 ? 0 : (int) (held >>> 32);
    int likest = pair < 0 ? -1 : likeness(queryAtom(pair), targetAtom(pair));
    clock.tick(1 + query.degree(frontierAtom));
    if (pair >= 0) {
      // No pair of the atom is liker than its likest level, nor has more common bonds than the atom
      // has mapped neighbours, so its pairs are not looked at when none can beat the held one.
      // Its likest level is worked out only when it decides.
      int mappedNeighbours = 0;
      for (int queryPlace = query.start[frontierAtom]; queryPlace < query.start[frontierAtom + 1];
           queryPlace++)
        if (queryToTarget[query.neighbour[queryPlace]] >= 0)
          mappedNeighbours++;
      boolean fewerBonds = mappedNeighbours < links
          || mappedNeighbours == links && queryRank[frontierAtom] < queryRank[queryAtom(pair)];
      if (fewerBonds ? likest == CommonBondGraphs.RADII || likestLevel(frontierAtom) <= likest
                     : likest > 0 && likestLevel(frontierAtom) < likest)
        return held;
    }
    for (int queryPlace = query.start[frontierAtom]; queryPlace < query.start[frontierAtom + 1];
         queryPlace++) {
      int neighbourImage = queryToTarget[query.neighbour[queryPlace]];
      if (neighbourImage < 0)
        continue;
      for (int targetPlace = target.start[neighbourImage] + leadingTaken(neighbourImage, round),
               freeCount = 0;
           targetPlace < target.start[neighbourImage + 1] && freeCount < HUB_TARGETS;
           targetPlace++) {
        int candidateTarget = target.neighbour[targetPlace];
        clock.tick();
        if (targetToQuery[candidateTarget] >= 0 || !has(allowed, frontierAtom, candidateTarget)
            || !graph.bondsMatch(query.neighbourBond[queryPlace], target.neighbourBond[targetPlace]))
          continue;
        freeCount++;
        int level = likeness(frontierAtom, candidateTarget);
        if (level < likest)
          continue;
        int newBonds = common(frontierAtom, candidateTarget);
        if (level > likest || newBonds > links
            || newBonds == links
                && (queryRank[frontierAtom] > queryRank[queryAtom(pair)]
                    || frontierAtom == queryAtom(pair)
                        && targetRank[candidateTarget] < targetRank[targetAtom(pair)])) {
          pair = pair(frontierAtom, candidateTarget);
          likest = level;
          links = newBonds;
        }
      }
    }
    return pair < 0 ? -1 : (long) links << 32 | pair;
  }

  /**
   * The likest level of {@code queryAtom} with any target it is compatible with, worked out when
   * first needed.
   */
  private int likestLevel(int queryAtom) {
    if (likestLevel[queryAtom] < 0) {
      int level = 0;
      long[] compatible = edgeBound.compatible;
      for (int targetAtom = CommonBondGraphs.nextBit(compatible, queryAtom * words, words, 0);
           targetAtom >= 0 && level < CommonBondGraphs.RADII;
           targetAtom = CommonBondGraphs.nextBit(compatible, queryAtom * words, words, targetAtom + 1)) {
        clock.tick();
        level = Math.max(level, likeness(queryAtom, targetAtom));
      }
      likestLevel[queryAtom] = level;
    }
    return likestLevel[queryAtom];
  }

  /**
   * The number of taken targets that lead the list of neighbours of {@code targetAtom}, for the
   * seed of {@code round}. A seed only takes targets while it grows, so the count only rises, and
   * it is counted afresh for each seed.
   */
  private int leadingTaken(int targetAtom, int round) {
    if (leadingRound[targetAtom] != round) {
      leadingRound[targetAtom] = round;
      leadingTaken[targetAtom] = 0;
    }
    int taken = leadingTaken[targetAtom];
    while (taken < target.degree(targetAtom)
        && targetToQuery[target.neighbour[target.start[targetAtom] + taken]] >= 0) {
      clock.tick();
      taken++;
    }
    leadingTaken[targetAtom] = taken;
    return taken;
  }

  /** The number of radii, from 0, at which {@code queryAtom} and {@code targetAtom} look alike. */
  private int likeness(int queryAtom, int targetAtom) {
    int level = 0;
    while (level < CommonBondGraphs.RADII
        && querySignature[level][queryAtom] == targetSignature[level][targetAtom])
      level++;
    return level;
  }

  /**
   * The mappings with the best score, at most {@code maxMappings}. Each is found once, as a
   * mapping is reached only from its highest ranked query atom and under the first sibling pair
   * it contains; scouting may meet a mapping before its turn, so the mappings held are kept
   * distinct. The search starts from the score {@code start}, which the {@code seed} holds
   * unless it is null. It ends early once it holds enough mappings that reach
   * {@code upperBound}. It needs {@link #rank} first.
   */
  List<int[]> collect(int[] seed, long start, long upperBound, int maxMappings) {
    this.bestScore = start;
    this.upperBound = upperBound;
    this.maxMappings = maxMappings;
    heldSeed = seed;
    mappings = new ArrayList<>();
    held.clear();
    if (seed != null) {
      mappings.add(seed);
      held.add(new Record(seed));
    }
    findTwins();
    roots(queryCount - 1, targetCount - 1, Long.MAX_VALUE, Long.MAX_VALUE);
    return mappings;
  }

  /**
   * Searches below each root pair in turn, from the target at place {@code firstPlace} of the
   * query atom of rank {@code firstRank}, until the pass can end or {@code nodeEnd} nodes are
   * reached, and each root pair for at most {@code rootNodes} nodes. Below a root only lower
   * ranked query atoms are mapped, so the roots of low rank have less to reach.
   */
  private void roots(int firstRank, int firstPlace, long nodeEnd, long rootNodes) {
    long limit = nodeLimit;
    for (int rank = firstRank;
         rank >= lowestRoot && !done() && nodes < nodeEnd && canImprove(rootBound(rank)); rank--) {
      int queryAtom = queryOrder[rank];
      for (int targetPosition = rank == firstRank ? firstPlace : targetCount - 1;
           targetPosition >= 0 && !done() && nodes < nodeEnd; targetPosition--) {
        int targetAtom = targetOrder[targetPosition];
        clock.tick();
        if (!has(allowed, queryAtom, targetAtom))
          continue;
        if (rootRefuted(queryAtom, targetAtom)) {
          refutedRoots++;
          continue;
        }
        if (rootNodes == Long.MAX_VALUE) {
          rootRank = rank;
          rootPlace = targetPosition;
        } else {
          nodeLimit = Math.min(limit, Math.min(nodeEnd, nodes + rootNodes));
        }
        if (push(queryAtom, targetAtom, 0, upperBound))
          search();
        else if (rootBonds != null)
          rootBonds[pair(queryAtom, targetAtom)] = bonds(lastBound);
        nodeLimit = limit;
      }
    }
  }

  /**
   * Scouts once the search under way has mapped many pairs, to find a better mapping soon when
   * the first root pairs hold none: the state of the search is put aside, the root pairs after
   * the one it is below are searched for at most {@link ScoutLimits#root} nodes each and
   * {@link ScoutLimits#budget} in all, and the state is put back. A better mapping found raises
   * the best score, which only prunes more of what the search has left. Each gap to the next
   * scouting is {@link ScoutLimits#growth} times the one before, so scouting takes a small share
   * of a long search.
   */
  private void scout() {
    scoutAt = Long.MAX_VALUE;
    scoutings++;
    // The levels of the search are put aside, and scouting makes its own.
    Level[] saved = Arrays.copyOf(levels, depth);
    Arrays.fill(levels, 0, depth, null);
    int[] savedQueryToTarget = queryToTarget.clone(), savedTargetToQuery = targetToQuery.clone();
    int[] savedExcluded = Arrays.copyOf(excluded, excludedCount);
    long[] savedAllowed = allowed.clone();
    int savedDepth = depth, savedCommonBonds = commonBonds, savedExcludedCount = excludedCount;
    int savedFragments = fragments;
    depth = commonBonds = excludedCount = fragments = 0;
    Arrays.fill(queryToTarget, -1);
    Arrays.fill(targetToQuery, -1);
    System.arraycopy(edgeBound.compatible, 0, allowed, 0, allowed.length);
    int rank = rootRank, place = rootPlace;
    roots(place > 0 ? rank : rank - 1, place > 0 ? place - 1 : targetCount - 1,
        nodes + scoutLimits.budget, scoutLimits.root);
    rootRank = rank;
    rootPlace = place;
    for (int i = 0; i < savedDepth; i++) levels[i] = saved[i];
    System.arraycopy(savedQueryToTarget, 0, queryToTarget, 0, queryCount);
    System.arraycopy(savedTargetToQuery, 0, targetToQuery, 0, targetCount);
    System.arraycopy(savedExcluded, 0, excluded, 0, savedExcludedCount);
    System.arraycopy(savedAllowed, 0, allowed, 0, allowed.length);
    excludedCount = savedExcludedCount;
    depth = savedDepth;
    commonBonds = savedCommonBonds;
    fragments = savedFragments;
    scoutGap *= scoutLimits.growth;
    scoutAt = nodes + scoutGap;
  }

  /**
   * Whether a root pair can be skipped, as the bound of its state, when it was last pruned,
   * shows that nothing below it beats the best score. That bound holds for every later pass:
   * the pairs are never excluded at the top, and the upper bound only falls.
   */
  private boolean rootRefuted(int queryAtom, int targetAtom) {
    if (rootBonds == null)
      return false;
    int bonds = rootBonds[pair(queryAtom, targetAtom)];
    return bonds >= 0 && !canImprove(capped(Math.min(queryCount, targetCount), bonds));
  }

  /**
   * Searches depth first from the root pair until it is taken off again. The top level tries
   * the targets of each of its candidates in turn, and then, unless the search is connected,
   * the roots of new fragments. It is taken off once nothing is left or its bound is too low.
   * The levels are an explicit stack, so large molecules do not run out of stack.
   */
  private void search() {
    while (depth > 0) {
      clock.tick();
      if (nodes >= scoutAt)
        scout();
      Level level = levels[depth - 1];
      if (done() || !canImprove(level.bound)) {
        pop();
      } else if (level.nextTarget < level.targetsCount) {
        int queryAtom = level.triedAtom, targetAtom = level.targets[level.nextTarget++];
        if (twinExcluded(queryAtom, targetAtom)
            || !push(queryAtom, targetAtom, level.fragmentStart, level.bound))
          exclude(queryAtom, targetAtom);
      } else if (level.nextCandidate < level.candidateCount) {
        // The targets tried so far are now excluded, which can lower the bound.
        if (level.nextCandidate > 0 && (level.bound = bound(level.fragmentStart, level.bound)) < 0) {
          pop();
        } else {
          level.triedAtom = level.candidates[level.nextCandidate++];
          targets(level.triedAtom, level);
          level.nextTarget = 0;
        }
      } else if (!nextRoot(level)) {
        pop();
      }
    }
  }

  /**
   * Maps {@code queryAtom} onto {@code targetAtom} as a new level of the fragment that starts
   * at level {@code fragmentStart}, with the bound of the state it extends until it has its own.
   * False, with the pair left unmapped, if the bound prunes it.
   */
  private boolean push(int queryAtom, int targetAtom, int fragmentStart, long bound) {
    Level level = map(queryAtom, targetAtom, fragmentStart, bound);
    record(fragmentStart);
    return isLeaf() || branch(level, fragmentStart, bound);
  }

  /**
   * Maps {@code queryAtom} onto {@code targetAtom} as the level at {@link #depth}, of the
   * fragment that starts at level {@code fragmentStart} and with {@code bound} until the level
   * has its own. The level records what taking the pair off restores, the common bonds and the
   * excluded pairs before it, and for the bound the query and target bonds with no end mapped
   * once it is. It has nothing to try until {@link #branch} gives it candidates.
   */
  private Level map(int queryAtom, int targetAtom, int fragmentStart, long bound) {
    countNode();
    if (levels[depth] == null)
      levels[depth] = new Level();
    Level level = levels[depth];
    level.queryAtom = queryAtom;
    level.targetAtom = targetAtom;
    level.startsFragment = fragmentStart == depth;
    if (level.startsFragment)
      fragments++;
    level.commonBondsBefore = commonBonds;
    level.excludedBefore = excludedCount;
    level.freeQueryBonds = depth > 0 ? levels[depth - 1].freeQueryBonds : query.bondCount;
    level.freeTargetBonds = depth > 0 ? levels[depth - 1].freeTargetBonds : target.bondCount;
    for (int queryPlace = query.start[queryAtom]; queryPlace < query.start[queryAtom + 1];
         queryPlace++)
      if (queryToTarget[query.neighbour[queryPlace]] < 0)
        level.freeQueryBonds--;
    for (int targetPlace = target.start[targetAtom]; targetPlace < target.start[targetAtom + 1];
         targetPlace++)
      if (targetToQuery[target.neighbour[targetPlace]] < 0)
        level.freeTargetBonds--;
    level.fragmentStart = fragmentStart;
    level.bound = bound;
    level.candidateCount = level.targetsCount = 0;
    level.nextCandidate = level.nextTarget = level.rootRank = 0;
    commonBonds += common(queryAtom, targetAtom);
    queryToTarget[queryAtom] = targetAtom;
    targetToQuery[targetAtom] = queryAtom;
    depth++;
    return level;
  }

  /**
   * Whether the mapping is a leaf, with nothing left to try: the search holds what it wants,
   * every query atom up to the rank of the root is mapped, or every target atom is, or the
   * mapping has the common bonds of the upper bound.
   */
  private boolean isLeaf() {
    // A mapping that drops bonds for its conditions may select fewer than its common bonds, and
    // a larger mapping may then select more, up to the upper bound.
    return done() || depth == queryRank[levels[0].queryAtom] + 1 || depth == targetCount
        || commonBonds == bonds(upperBound) && (conditions == null || !conditions.dropsBonds());
  }

  /**
   * Gives the top {@code level}, whose fragment starts at level {@code fragmentStart} and whose
   * bound is {@code bound} so far, its own bound and the candidates whose targets it tries.
   * Unless the search is connected, a fragment with a bond also offers the roots of new
   * fragments once its candidates are done; these rank below its own root. False, with the pair
   * taken off, if the bound prunes the level.
   */
  private boolean branch(Level level, int fragmentStart, long bound) {
    level.bound = bound(fragmentStart, bound);
    if (level.bound < 0) {
      takeOff();
      return false;
    }
    candidates(level);
    clock.tick(level.candidateCount);
    if (!connected && fragmentBonds(fragmentStart) >= MIN_BONDS)
      level.rootRank = queryRank[levels[fragmentStart].queryAtom];
    return true;
  }

  /** The common bonds of the open fragment, which starts at level {@code fragmentStart}. */
  private int fragmentBonds(int fragmentStart) {
    return commonBonds - levels[fragmentStart].commonBondsBefore;
  }

  /** Takes the top level off, and excludes its pair for its later siblings and everything below them. */
  private void pop() {
    Level level = levels[depth - 1];
    takeOff();
    if (depth > 0)
      exclude(level.queryAtom, level.targetAtom);
  }

  /**
   * Counts a pair mapped or an embedding atom placed, and stops the search with
   * {@link CommonBondGraphs.Stop} instead once {@link #nodeCap} are counted.
   */
  private void countNode() {
    if (nodes >= nodeCap) {
      capReached = true;
      throw new CommonBondGraphs.Stop();
    }
    nodes++;
  }

  /** Unmaps the pair of the top level, and allows again the pairs excluded since it was mapped. */
  private void takeOff() {
    Level level = levels[--depth];
    if (level.startsFragment)
      fragments--;
    while (excludedCount > level.excludedBefore) {
      int pair = excluded[--excludedCount];
      set(allowed, queryAtom(pair), targetAtom(pair));
    }
    queryToTarget[level.queryAtom] = -1;
    targetToQuery[level.targetAtom] = -1;
    commonBonds = level.commonBondsBefore;
  }

  /** Excludes the pair of {@code queryAtom} and {@code targetAtom}, which must be allowed. */
  private void exclude(int queryAtom, int targetAtom) {
    clear(allowed, queryAtom, targetAtom);
    if (excludedCount == excluded.length)
      excluded = Arrays.copyOf(excluded, 2 * excluded.length);
    excluded[excludedCount++] = pair(queryAtom, targetAtom);
  }

  /**
   * Offers the next root of a new fragment once the candidates of the level are done. This
   * closes the fragment, as every pair that extends it is excluded. A root is the next free
   * query atom ranked below the root of the closed fragment, so a mapping is reached only from
   * the highest ranked atom of each fragment. False when no root is left, when the bound of
   * those left is too low, or when the search is connected. That bound falls with their rank.
   */
  private boolean nextRoot(Level level) {
    if (level.rootRank == 0)
      return false;
    // The fragment is closed once it is set to start at the next level, where the roots offered
    // start fragments of their own. Its atoms then lie in no larger fragment, and the fragments
    // only add up, so a fragment too small or one too many leaves nothing to find below.
    if (level.fragmentStart < depth) {
      if (hasPendantPair(level.fragmentStart) || conditions != null
          && (depth - level.fragmentStart < conditions.minFragmentAtoms()
              || fragments >= conditions.maxFragments()))
        return false;
      level.fragmentStart = depth;
    }
    while (--level.rootRank >= 0) {
      int queryAtom = queryOrder[level.rootRank];
      if (queryToTarget[queryAtom] >= 0)
        continue;
      level.bound = rootBound(level.rootRank);
      if (!canImprove(level.bound))
        return false;
      level.triedAtom = queryAtom;
      freeTargets(queryAtom, level);
      level.nextTarget = 0;
      if (level.targetsCount > 0)
        return true;
    }
    return false;
  }

  /**
   * Keeps the mapping if it beats the best score, or ties with it while fewer than wanted are
   * held, and otherwise keeps it as the fallback if it beats that. The open fragment, from
   * level {@code fragmentStart}, must have a common bond.
   */
  private void record(int fragmentStart) {
    long current = CommonBondGraphs.score(depth, commonBonds);
    boolean full = mappings.size() == maxMappings;
    if (current < bestScore || current == bestScore && full) {
      // Below the threshold of a level, a mapping better than the fallback becomes the fallback.
      if (fragmentBonds(fragmentStart) >= MIN_BONDS)
        offerFallback(queryToTarget, current);
      return;
    }
    // A fragment needs a common bond. The open one may still get it further down.
    if (fragmentBonds(fragmentStart) < MIN_BONDS)
      return;
    if (conditions == null) {
      keep(current, queryToTarget.clone());
      return;
    }
    // The best selections that meet the conditions score at most as much. One better than the
    // fallback serves as that while the best score is out of reach.
    selections.clear();
    long least = Math.min(fallbackScore + 1, full ? bestScore + 1 : bestScore);
    long valid = conditions.select(queryToTarget, least, maxMappings, selections);
    if (valid < 0)
      return;
    if (valid < bestScore || valid == bestScore && full) {
      if (valid > fallbackScore) {
        fallback = selections.get(0);
        fallbackScore = valid;
      }
      return;
    }
    for (int[] selection : selections)
      if (!keep(valid, selection))
        break;
  }

  /**
   * Holds {@code record}, of score {@code score}, which is no lower than the best: a better
   * one replaces those held. False, holding nothing, once as many as wanted are held.
   */
  private boolean keep(long score, int[] record) {
    if (score > bestScore) {
      bestScore = score;
      mappings.clear();
      held.clear();
      heldSeed = null;
    } else if (heldSeed != null && Arrays.equals(record, heldSeed)) {
      // The seed is held from the start, so it is not added again when found.
      heldSeed = null;
      return true;
    }
    if (mappings.size() == maxMappings)
      return false;
    // A single mapping is found once; several may meet one met before while scouting.
    if (maxMappings > 1 && !held.add(new Record(record)))
      return true;
    mappings.add(record);
    return true;
  }

  /** A record held, compared by its contents. */
  private static final class Record {
    final int[] values;
    final int hash;

    Record(int[] values) {
      this.values = values;
      hash = Arrays.hashCode(values);
    }

    @Override
    public boolean equals(Object other) {
      return other instanceof Record && Arrays.equals(values, ((Record) other).values);
    }

    @Override
    public int hashCode() {
      return hash;
    }
  }

  /**
   * Whether the pass can end, as it holds as many mappings as wanted with a score that no
   * mapping beats, or as it has mapped as many pairs as it may.
   */
  private boolean done() {
    return mappings.size() == maxMappings && bestScore == upperBound || nodes >= nodeLimit;
  }

  /**
   * The bound of the state whose open fragment starts at level {@code fragmentStart}, or -1 if
   * the state is pruned. A mapping below adds to the open fragment only open pairs that it
   * reaches through compatible bonds, so the bound counts these. The first layer of them gives
   * the candidates.
   */
  private long bound(int fragmentStart, long current) {
    reach.start();
    int open = 0;
    for (int i = fragmentStart; i < depth; i++) {
      int queryLinks = reach.queryLinks, targetLinks = reach.targetLinks;
      reach.link(levels[i].queryAtom, levels[i].targetAtom);
      open += Math.min(reach.queryLinks - queryLinks, reach.targetLinks - targetLinks);
    }
    reach.firstLayer = reach.count;
    Level top = levels[depth - 1];
    // The atoms of the upper bound cap only a mapping with as many bonds. A mapping with fewer
    // bonds may have more atoms.
    int mostAtoms = Math.min(queryCount, targetCount);
    // Centres that already disagree on three ligands keep disagreeing below, so a mapping there
    // leaves bonds out for them, see CommonBondConditions#persistentLoss.
    int loss = conditions != null ? conditions.persistentLoss(queryToTarget) : 0;
    int mostBonds = commonBonds + open + Math.min(top.freeQueryBonds, top.freeTargetBonds) - loss;
    long result = Math.min(current, capped(mostAtoms, mostBonds));
    // A new atom lies within as many bonds of the mapped atoms as a better mapping adds. One
    // that drops bonds for its conditions may keep fewer of the common bonds so far, so it may
    // add up to all its bonds.
    int newBonds = conditions != null && conditions.dropsBonds() ? bonds(result)
                                                                 : bonds(result) - commonBonds;
    if (canImprove(result))
      result = Math.min(result, edgeBound.bound(queryRank[levels[fragmentStart].queryAtom], open,
          newBonds));
    reach.finish();
    lastBound = result;
    return canImprove(result) ? result : -1;
  }

  /**
   * The bound of the roots offered up to {@code rankLimit}. The fragment is closed, so no bond
   * links it. Centres that already disagree on three ligands lower it, as in {@link #bound}.
   */
  private long rootBound(int rankLimit) {
    reach.start();
    long bound = edgeBound.bound(rankLimit, 0, 0);
    int loss = depth > 0 && conditions != null ? conditions.persistentLoss(queryToTarget) : 0;
    return loss == 0 ? bound : capped(atoms(bound), bonds(bound) - loss);
  }

  /** The atoms of a score of this search. */
  private static int atoms(long score) {
    return CommonBondGraphs.atoms(score);
  }

  /** The common bonds of a score of this search. */
  private static int bonds(long score) {
    return CommonBondGraphs.bonds(score);
  }

  /**
   * The score of {@code mapping}, a record (see {@link #match}): its atoms, and the compatible
   * bonds between them that the target has between their images, less those it leaves out.
   */
  long score(int[] mapping) {
    int atomCount = 0, bondCount = -leftOut(mapping, queryCount);
    for (int queryAtom = 0; queryAtom < queryCount; queryAtom++) {
      if (mapping[queryAtom] >= 0) {
        atomCount++;
        bondCount += commonBondsAt(mapping, queryAtom);
      }
    }
    return CommonBondGraphs.score(atomCount, bondCount);
  }

  /**
   * The compatible bonds from {@code queryAtom} to the lower mapped atoms of {@code mapping}
   * that the target has between their images, so that each common bond counts once.
   */
  private int commonBondsAt(int[] mapping, int queryAtom) {
    int targetAtom = mapping[queryAtom], count = 0;
    for (int queryPlace = query.start[queryAtom]; queryPlace < query.start[queryAtom + 1];
         queryPlace++) {
      int other = query.neighbour[queryPlace];
      if (other >= queryAtom || mapping[other] < 0)
        continue;
      int targetBond = target.bondBetween(targetAtom, mapping[other]);
      if (targetBond >= 0 && graph.bondsMatch(query.neighbourBond[queryPlace], targetBond))
        count++;
    }
    return count;
  }

  /**
   * A bound, capped by the upper bound: a mapping with the common bonds of the upper bound has
   * no more atoms than it.
   */
  private long capped(int atomCount, int bondCount) {
    int mostAtoms = atoms(upperBound), mostBonds = bonds(upperBound);
    if (bondCount >= mostBonds)
      return CommonBondGraphs.score(Math.min(atomCount, mostAtoms), mostBonds);
    return CommonBondGraphs.score(atomCount, bondCount);
  }

  /** Whether a state with this bound can hold a better mapping, or a tie while fewer than wanted are held. */
  private boolean canImprove(long bound) {
    return bound > bestScore || bound == bestScore && mappings.size() < maxMappings;
  }

  /**
   * Gives {@code level} as candidates the query atoms of the first layer of reached pairs: not
   * terminal first, then fewer pairs, then higher rank. The reached pairs keep their order, as
   * the bound reads them in it.
   */
  private void candidates(Level level) {
    int layerCount = reach.firstLayer;
    int[] layer = layerScratch = grown(layerScratch, layerCount);
    System.arraycopy(reach.pairs, 0, layer, 0, layerCount);
    Arrays.sort(layer, 0, layerCount);
    long[] keys = keyScratch = grown(keyScratch, layerCount);
    int count = 0;
    for (int from = 0, to; from < layerCount; from = to) {
      int queryAtom = queryAtom(layer[from]);
      to = from + 1;
      while (to < layerCount && queryAtom(layer[to]) == queryAtom) to++;
      keys[count++] = candidateKey(queryAtom, to - from);
    }
    Arrays.sort(keys, 0, count);
    int[] sorted = level.candidates = grown(level.candidates, count);
    for (int i = 0; i < count; i++) sorted[i] = candidateOf(keys[i]);
    level.candidateCount = count;
  }

  /** {@code array}, or a longer one with nothing kept when it holds fewer than {@code need}. */
  private static int[] grown(int[] array, int need) {
    return array.length >= need ? array : new int[Math.max(need, 2 * array.length)];
  }

  /** {@code array}, or a longer one with nothing kept when it holds fewer than {@code need}. */
  private static long[] grown(long[] array, int need) {
    return array.length >= need ? array : new long[Math.max(need, 2 * array.length)];
  }

  /** The sort key of candidate {@code queryAtom}: whether it is terminal, its number of pairs, then its rank in reverse. */
  private long candidateKey(int queryAtom, int pairs) {
    return (query.degree(queryAtom) == 1 ? 1L << 62 : 0) | (long) pairs << 31
        | queryCount - 1 - queryRank[queryAtom];
  }

  /** The candidate of a key from {@link #candidateKey}. */
  private int candidateOf(long key) {
    return queryOrder[queryCount - 1 - (int) (key & Integer.MAX_VALUE)];
  }

  /**
   * Gives {@code level} as targets those of query atom {@code queryAtom}. These are the free
   * neighbours of the images of its mapped neighbours, joined by compatible bonds, in pairs not
   * excluded. The most alike come first, then by descending rank. The order decides only which
   * mappings are found first.
   */
  private void targets(int queryAtom, Level level) {
    int capacity = 0;
    for (int queryPlace = query.start[queryAtom]; queryPlace < query.start[queryAtom + 1];
         queryPlace++)
      if (queryToTarget[query.neighbour[queryPlace]] >= 0)
        capacity += target.degree(queryToTarget[query.neighbour[queryPlace]]);
    long[] keys = keyScratch = grown(keyScratch, capacity);
    int count = 0;
    for (int queryPlace = query.start[queryAtom]; queryPlace < query.start[queryAtom + 1];
         queryPlace++) {
      int neighbourImage = queryToTarget[query.neighbour[queryPlace]];
      if (neighbourImage < 0)
        continue;
      clock.tick(1 + target.degree(neighbourImage));
      for (int targetPlace = target.start[neighbourImage];
           targetPlace < target.start[neighbourImage + 1]; targetPlace++) {
        int targetAtom = target.neighbour[targetPlace];
        if (targetToQuery[targetAtom] < 0 && has(allowed, queryAtom, targetAtom)
            && graph.bondsMatch(query.neighbourBond[queryPlace], target.neighbourBond[targetPlace]))
          keys[count++] = targetKey(queryAtom, targetAtom);
      }
    }
    // A target next to two mapped neighbours is listed twice, with the same key.
    Arrays.sort(keys, 0, count);
    int[] targets = level.targets = grown(level.targets, count);
    int unique = 0;
    for (int i = 0; i < count; i++)
      if (i == 0 || keys[i] != keys[i - 1])
        targets[unique++] = targetOf(keys[i]);
    level.targetsCount = unique;
  }

  /** The sort key of target {@code targetAtom} of query atom {@code queryAtom}: how unlike the two are, then the rank in reverse. */
  private long targetKey(int queryAtom, int targetAtom) {
    return (long) (CommonBondGraphs.RADII - likeness(queryAtom, targetAtom)) << 32
        | targetCount - 1 - targetRank[targetAtom];
  }

  /** The target of a key from {@link #targetKey}. */
  private int targetOf(long key) {
    return targetOrder[targetCount - 1 - (int) key];
  }

  /** Gives {@code level} as targets the free targets of root {@code queryAtom}, by descending rank. */
  private void freeTargets(int queryAtom, Level level) {
    int[] targets = level.targets = grown(level.targets, targetCount);
    int count = 0;
    clock.tick(targetCount);
    for (int targetPosition = targetCount - 1; targetPosition >= 0; targetPosition--)
      if (targetToQuery[targetOrder[targetPosition]] < 0
          && has(allowed, queryAtom, targetOrder[targetPosition]))
        targets[count++] = targetOrder[targetPosition];
    level.targetsCount = count;
  }

  /** The compatible bonds from {@code queryAtom} to mapped atoms that the target has from {@code targetAtom} to their images. */
  private int common(int queryAtom, int targetAtom) {
    int count = 0;
    for (int queryPlace = query.start[queryAtom]; queryPlace < query.start[queryAtom + 1];
         queryPlace++) {
      int neighbourImage = queryToTarget[query.neighbour[queryPlace]];
      if (neighbourImage < 0)
        continue;
      clock.tick(1 + target.degree(targetAtom));
      int targetBond = target.bondBetween(targetAtom, neighbourImage);
      if (targetBond >= 0 && graph.bondsMatch(query.neighbourBond[queryPlace], targetBond))
        count++;
    }
    return count;
  }

  /**
   * Finds the twins once, before any pair is excluded. A query that fits whole needs them only
   * while its embeddings are tried under conditions.
   */
  private void findTwins() {
    if (queryTwin == null) {
      queryTwin = twins(query, true);
      targetTwin = twins(target, false);
    }
  }

  /**
   * The twins of one molecule, each linked to the next in a ring, -1 for none. Twins are
   * terminal atoms with the same neighbour that match the same atoms by the same bonds and,
   * under conditions, carry with that neighbour no stereo tag that the conditions read, see
   * {@link CommonBondConditions#untagged}. Swapping two takes every mapping onto one of the same
   * score.
   */
  private int[] twins(CommonBondGraphs.Side side, boolean onQuery) {
    int[] next = new int[side.atomCount];
    Arrays.fill(next, -1);
    for (int parent = 0; parent < side.atomCount; parent++) {
      if (!untagged(onQuery, parent))
        continue;
      for (int place = side.start[parent]; place < side.start[parent + 1]; place++) {
        int first = side.neighbour[place], last = first;
        if (next[first] >= 0 || side.degree(first) != 1 || !untagged(onQuery, first))
          continue;
        for (int laterPlace = place + 1; laterPlace < side.start[parent + 1]; laterPlace++) {
          int other = side.neighbour[laterPlace];
          if (next[other] < 0 && side.degree(other) == 1 && untagged(onQuery, other)
              && alike(onQuery, first, side.neighbourBond[place], other,
                  side.neighbourBond[laterPlace])) {
            next[last] = other;
            last = other;
          }
        }
        if (last != first)
          next[last] = first;
      }
    }
    return next;
  }

  /** Whether {@code atom} of one molecule carries no stereo tag that the conditions read, if there are any. */
  private boolean untagged(boolean onQuery, int atom) {
    return conditions == null || conditions.untagged(onQuery, atom);
  }

  /** Whether {@code atom} and {@code other}, with their bonds {@code bond} and {@code otherBond}, match alike. */
  private boolean alike(boolean onQuery, int atom, int bond, int other, int otherBond) {
    clock.tick(words + queryCount + query.bondCount + target.bondCount);
    if (onQuery) {
      for (int word = 0; word < words; word++)
        if (allowed[atom * words + word] != allowed[other * words + word])
          return false;
      for (int targetBond = 0; targetBond < target.bondCount; targetBond++)
        if (graph.bondsMatch(bond, targetBond) != graph.bondsMatch(otherBond, targetBond))
          return false;
    } else {
      for (int queryAtom = 0; queryAtom < queryCount; queryAtom++)
        if (has(allowed, queryAtom, atom) != has(allowed, queryAtom, other))
          return false;
      for (int queryBond = 0; queryBond < query.bondCount; queryBond++)
        if (graph.bondsMatch(queryBond, bond) != graph.bondsMatch(queryBond, otherBond))
          return false;
    }
    return true;
  }

  /**
   * Whether the pair of {@code queryAtom} and {@code targetAtom} can be skipped, as a free twin
   * of one has its pair with the other excluded. Swapping the twins turns every mapping below
   * this pair into one of the same score below the excluded pair, which was tried before. So
   * this holds only while no mapping of that score could be recorded. That is when as many are
   * held as wanted, and also when nothing but the seed is held, see {@link #collectByBonds}.
   */
  private boolean twinExcluded(int queryAtom, int targetAtom) {
    if (conditions != null && !conditions.keepsTwins())
      return false;
    boolean full = mappings.size() == maxMappings;
    boolean onlySeed = mappings.size() == (heldSeed != null ? 1 : 0);
    if (!full && !onlySeed)
      return false;
    for (int twin = queryTwin[queryAtom]; twin >= 0 && twin != queryAtom; twin = queryTwin[twin])
      if (queryToTarget[twin] < 0 && !has(allowed, twin, targetAtom))
        return true;
    for (int twin = targetTwin[targetAtom]; twin >= 0 && twin != targetAtom;
         twin = targetTwin[twin])
      if (targetToQuery[twin] < 0 && !has(allowed, queryAtom, twin))
        return true;
    return false;
  }

  /**
   * Whether the open fragment, from level {@code fragmentStart}, leaves out a pendant pair. That
   * is a free terminal query atom next to a mapped atom whose image has a free terminal target
   * atom it may map onto by a compatible bond. No mapping below the roots the level offers is
   * then a maximum one, as adding that pair beats it.
   */
  private boolean hasPendantPair(int fragmentStart) {
    if (conditions != null && !conditions.keepsPendantPairs())
      return false;
    for (int level = fragmentStart; level < depth; level++) {
      int queryAtom = levels[level].queryAtom, targetAtom = levels[level].targetAtom;
      clock.tick(query.degree(queryAtom) * target.degree(targetAtom));
      for (int queryPlace = query.start[queryAtom]; queryPlace < query.start[queryAtom + 1];
           queryPlace++) {
        int queryNeighbour = query.neighbour[queryPlace];
        if (queryToTarget[queryNeighbour] >= 0 || query.degree(queryNeighbour) != 1)
          continue;
        for (int targetPlace = target.start[targetAtom];
             targetPlace < target.start[targetAtom + 1]; targetPlace++) {
          int targetNeighbour = target.neighbour[targetPlace];
          if (targetToQuery[targetNeighbour] < 0 && target.degree(targetNeighbour) == 1
              && has(edgeBound.compatible, queryNeighbour, targetNeighbour)
              && graph.bondsMatch(query.neighbourBond[queryPlace], target.neighbourBond[targetPlace]))
            return true;
        }
      }
    }
    return false;
  }

  /**
   * The search for the whole query, if it may fit. It may when the query is connected and no
   * larger than the target, and each query atom can have a target atom of its own with at least
   * as many bonds. Null if not. With {@code prunesTwins} it tries only embeddings that place the
   * twins of each molecule in order, see {@link Embedding#twinsOutOfOrder}.
   */
  Embedding embedding(boolean prunesTwins) {
    if (queryCount > targetCount || !queryConnected())
      return null;
    // A row for each number of bonds holds the target atoms with at least that many, up to the
    // most a query atom has. The last row is empty.
    int rows = Math.min(query.maxDegree(), target.maxDegree()) + 1;
    long[] atLeast = new long[(rows + 1) * words];
    for (int targetAtom = 0; targetAtom < targetCount; targetAtom++)
      for (int fewestBonds = Math.min(target.degree(targetAtom), rows); fewestBonds >= 0;
           fewestBonds--)
        set(atLeast, fewestBonds, targetAtom);
    long[] candidateRows = new long[queryCount * words];
    for (int queryAtom = 0; queryAtom < queryCount; queryAtom++) {
      clock.tick(words);
      int row = Math.min(query.degree(queryAtom), rows) * words;
      for (int i = 0; i < words; i++)
        candidateRows[queryAtom * words + i] = allowed[queryAtom * words + i] & atLeast[row + i];
    }
    return hasDistinctTargets(candidateRows, queryCount, words, clock)
        ? new Embedding(candidateRows, prunesTwins)
        : null;
  }

  /** Whether every query atom is reached from the first by bonds. */
  private boolean queryConnected() {
    if (queryCount == 0)
      return false;
    int[] queue = new int[queryCount];
    boolean[] reached = new boolean[queryCount];
    int count = 1;
    reached[0] = true;
    for (int head = 0; head < count; head++)
      for (int place = query.start[queue[head]]; place < query.start[queue[head] + 1]; place++)
        if (!reached[query.neighbour[place]]) {
          reached[query.neighbour[place]] = true;
          queue[count++] = query.neighbour[place];
        }
    return count == queryCount;
  }

  /**
   * The pair of {@code queryAtom} and {@code targetAtom} as one int,
   * {@code queryAtom * targetCount + targetAtom}, which fits as at most 2<sup>30</sup> pairs are
   * searched. Pairs sort by query atom, then by target atom.
   */
  private int pair(int queryAtom, int targetAtom) {
    return queryAtom * targetCount + targetAtom;
  }

  /** The query atom of a packed pair. */
  private int queryAtom(int pair) {
    return pair / targetCount;
  }

  /** The target atom of a packed pair. */
  private int targetAtom(int pair) {
    return pair % targetCount;
  }

  /**
   * A mapped pair of the branch and bound, one for each depth, with what it restores when taken
   * off and what it has left to try. The levels are made once and reused, see {@link #map}.
   */
  private static final class Level {

    /** The mapped pair. */
    int queryAtom, targetAtom;
    /** The common bonds and the excluded pairs before the pair was mapped, put back when it is taken off. */
    int commonBondsBefore, excludedBefore;
    /** The query and target bonds with no end mapped, once the pair is mapped. */
    int freeQueryBonds, freeTargetBonds;
    /** The level that starts the open fragment. Closing the fragment sets this to the next level. */
    int fragmentStart;
    /** A score that no mapping below the level beats, see {@link CommonBondSearch#bound}. */
    long bound;
    /**
     * The candidates whose targets the level tries, the first {@link #candidateCount} in order,
     * see {@link CommonBondSearch#candidates}. The array is kept for the next pair the level maps.
     */
    int[] candidates = NONE;
    /** The number of candidates, and the place in {@link #candidates} of the next to try. */
    int candidateCount, nextCandidate;
    /** The candidate being tried. */
    int triedAtom;
    /**
     * The targets of {@link #triedAtom}, the first {@link #targetsCount} in order, see
     * {@link CommonBondSearch#targets}. The array is kept as {@link #candidates} is.
     */
    int[] targets = NONE;
    /** The number of targets, and the place in {@link #targets} of the next to try. */
    int targetsCount, nextTarget;
    /** Once the fragment is closed, the roots of new fragments offered are ranked below this. */
    int rootRank;
    /** Whether the pair is the root of a fragment. */
    boolean startsFragment;
  }

  /**
   * The open pairs reached from the open fragment of a state through compatible bonds, see
   * {@link #bound}. The bonds that link them are counted once on each side.
   */
  private final class Reach {

    /** The reached pairs as bits, rows by query atom. */
    private final long[] seen = new long[queryCount * words];
    /** The reached pairs in the order they were reached, the first {@link #count} of the array. */
    private int[] pairs = new int[INITIAL_PAIRS];
    /** The number of reached pairs. */
    private int count;
    /** The number of pairs reached from the fragment itself, which come first and give the candidates. */
    private int firstLayer;
    /** The stamp of each query bond counted as a link, see {@link #stamp}. */
    private final int[] queryBondStamp = new int[query.bondCount];
    /** The stamp of each target bond counted as a link. */
    private final int[] targetBondStamp = new int[target.bondCount];
    /** The numbers of query and target bonds counted as links. */
    private int queryLinks, targetLinks;
    /** The stamp of this use, so that no marks have to be cleared. */
    private int stamp;

    /** Starts a use with nothing reached. Every mark is cleared when the stamp would overflow. */
    void start() {
      if (++stamp == Integer.MAX_VALUE) {
        Arrays.fill(queryBondStamp, 0);
        Arrays.fill(targetBondStamp, 0);
        stamp = 1;
      }
      count = queryLinks = targetLinks = 0;
    }

    /** Reaches the open pairs next to the pair of {@code queryAtom} and {@code targetAtom}. */
    void link(int queryAtom, int targetAtom) {
      int from = target.start[targetAtom], to = target.start[targetAtom + 1];
      int rootRank = queryRank[levels[0].queryAtom];
      for (int queryPlace = query.start[queryAtom]; queryPlace < query.start[queryAtom + 1];
           queryPlace++) {
        int queryNeighbour = query.neighbour[queryPlace], queryBond = query.neighbourBond[queryPlace];
        if (queryToTarget[queryNeighbour] >= 0 || queryRank[queryNeighbour] > rootRank)
          continue;
        clock.tick(1 + to - from);
        for (int targetPlace = from; targetPlace < to; targetPlace++) {
          int targetNeighbour = target.neighbour[targetPlace],
              targetBond = target.neighbourBond[targetPlace];
          if (targetToQuery[targetNeighbour] >= 0 || !has(allowed, queryNeighbour, targetNeighbour)
              || !graph.bondsMatch(queryBond, targetBond))
            continue;
          if (queryBondStamp[queryBond] != stamp) {
            queryBondStamp[queryBond] = stamp;
            queryLinks++;
          }
          if (targetBondStamp[targetBond] != stamp) {
            targetBondStamp[targetBond] = stamp;
            targetLinks++;
          }
          if (!has(seen, queryNeighbour, targetNeighbour)) {
            set(seen, queryNeighbour, targetNeighbour);
            if (count == pairs.length)
              pairs = Arrays.copyOf(pairs, 2 * pairs.length);
            pairs[count++] = pair(queryNeighbour, targetNeighbour);
          }
        }
      }
    }

    /** Ends this use by clearing the bits of the reached pairs. */
    void finish() {
      for (int i = 0; i < count; i++) clear(seen, queryAtom(pairs[i]), targetAtom(pairs[i]));
    }
  }

  /**
   * The bound from the bond labels. An open atom is one ranked up to a limit that has a pair
   * not excluded, and a mapping below the state maps open atoms onto the targets of these pairs.
   * A common bond and its image share a label. An open atom and its image get no more new bonds
   * of a label than either has open, so pairing the open atoms of each label bounds its new
   * bonds. Where a label holds fine labels that do not all match, its new bonds are also at most
   * the most bonds of matching fine labels that the open atoms can pair off, and the new atoms of
   * an atom class at most the most open atoms of matching kinds, see
   * {@link CommonBondGraphs.Groups}. In a connected search the open atoms are also cut to those
   * near the mapped atoms, and fall into groups that pair off apart, see {@link #groupedBound}.
   */
  private final class EdgeBound {

    /** The compatible pairs before any is excluded, see {@link #hasPendantPair}. */
    private final long[] compatible = allowed.clone();
    /** The labels of the query and of the target, with the open atoms of the state being bounded. */
    private final Labels queryLabels, targetLabels;
    /** The class at each slot. A label has a slot for each end, the second -1 if both ends share a class. */
    private final int[] slotSymbol;
    /** The most bonds of one slot at an atom, plus one, which is the width of a histogram row. */
    private final int histogramWidth;
    /** Scratch of {@link #countAtoms}: the open bonds of each slot at one atom, and the slots of that atom. */
    private final int[] slotCount, slotsAt;
    /** Scratch of {@link #bound}: the target atoms not used, and the targets of the open pairs. */
    private final long[] unusedTargets = new long[words], reachableTargets = new long[words];
    /** The groups of atom kinds and of fine labels. */
    private final CommonBondGraphs.Groups atomGroups = graph.atomGroups, labelGroups = graph.labelGroups;
    /**
     * The labels and the classes counted for a bound, the first {@link #labelsCounted} and
     * {@link #symbolsCounted}, each listed once while its flag is set, see {@link #pairCounted}.
     */
    private final int[] countedLabels = new int[graph.labels],
        countedSymbols = new int[graph.atomClasses];
    private final boolean[] labelListed = new boolean[graph.labels],
        symbolListed = new boolean[graph.atomClasses];
    private int labelsCounted, symbolsCounted;
    /**
     * The union-find of the open atoms of both molecules, the target atoms after the query
     * atoms, the group of each root while its mark is the stamp, the start of each group, the
     * molecules each group has atoms of, and the atoms by group, see {@link #groupedBound}.
     */
    private final int[] partParent = new int[queryCount + targetCount],
        groupMark = new int[queryCount + targetCount], groupOf = new int[queryCount + targetCount],
        groupStart = new int[queryCount + targetCount + 1],
        groupSides = new int[queryCount + targetCount], groupItems = new int[queryCount + targetCount];
    private int markStamp = 1;

    /** Sets up the slots of the labels of {@link CommonBondGraphs#labels}. */
    EdgeBound() {
      slotSymbol = new int[2 * graph.labels];
      for (int label = 0; label < graph.labels; label++) {
        slotSymbol[2 * label] = graph.labelLow[label];
        slotSymbol[2 * label + 1] =
            graph.labelLow[label] == graph.labelHigh[label] ? -1 : graph.labelHigh[label];
      }
      histogramWidth = histogramWidth(graph);
      int histogramSize = slotSymbol.length * histogramWidth, symbols = graph.atomClasses;
      queryLabels = new Labels(query, slots(query), histogramSize, symbols,
          atomGroups.queryGroup.length, labelGroups.queryGroup.length);
      targetLabels = new Labels(target, slots(target), histogramSize, symbols,
          atomGroups.targetGroup.length, labelGroups.targetGroup.length);
      slotCount = new int[slotSymbol.length];
      slotsAt = new int[histogramWidth];
    }

    /**
     * The bound of the state, where the open atoms are those ranked up to {@code rankLimit}. The
     * {@code open} new bonds at the mapped atoms are already counted. A new atom has a new common
     * bond and the class of its image. A connected extension of a mapping, of at most
     * {@code newBonds} new bonds, is bounded by groups of the open atoms near it instead, see
     * {@link #groupedBound}. A fragment of b new bonds has at most b + 1 atoms, and a connected
     * extension of a mapping at most b new atoms.
     */
    long bound(int rankLimit, int open, int newBonds) {
      Arrays.fill(unusedTargets, -1L);
      for (int i = 0; i < depth; i++) clear(unusedTargets, 0, levels[i].targetAtom);
      Arrays.fill(reachableTargets, 0);
      if (connected && depth > 0)
        return groupedBound(rankLimit, open, newBonds);
      for (int rank = 0; rank <= rankLimit; rank++) {
        int queryAtom = queryOrder[rank];
        if (queryToTarget[queryAtom] < 0 && markTargets(queryAtom))
          queryLabels.markOpen(queryAtom);
      }
      for (int targetAtom = CommonBondGraphs.nextBit(reachableTargets, 0, words, 0); targetAtom >= 0;
           targetAtom = CommonBondGraphs.nextBit(reachableTargets, 0, words, targetAtom + 1))
        targetLabels.markOpen(targetAtom);
      countAtoms(queryLabels, query, true, queryLabels.openAtoms, 0, queryLabels.openCount, 0,
          false);
      countAtoms(targetLabels, target, false, targetLabels.openAtoms, 0, targetLabels.openCount, 0,
          false);
      long counted = pairAll();
      int added = open + (int) (counted >>> 32), newAtoms = (int) counted;
      queryLabels.clearOpen();
      targetLabels.clearOpen();
      int mostNew = !connected ? added + added / MIN_BONDS : mostAtoms(added);
      return capped(depth + Math.min(newAtoms, mostNew), commonBonds + added);
    }

    /**
     * Adds the unused targets of the pairs of the free {@code queryAtom} to the reachable
     * targets, and tells whether it has any.
     */
    private boolean markTargets(int queryAtom) {
      long any = 0;
      clock.tick(words);
      for (int word = 0; word < words; word++) {
        long reachable = allowed[queryAtom * words + word] & unusedTargets[word];
        reachableTargets[word] |= reachable;
        any |= reachable;
      }
      return any != 0;
    }

    /**
     * The bound of a connected extension of the mapping by at most {@code newBonds} new bonds. A
     * new atom is joined to a mapped atom by a path of new atoms and common bonds. The first pair
     * of the path is one of the pairs reached from the mapped pairs, see {@link Reach#link}, and
     * the path goes on through open atoms and bonds that have a label, no more than
     * {@code newBonds} bonds in all. The query atoms so reached are open if they are ranked up to
     * {@code rankLimit} and have a pair not excluded, and the target atoms so reached if they are
     * the targets of such pairs of the open query atoms. The path lies in one part of the open
     * query atoms, linked by bonds that have a label, and its image in one part of the open target
     * atoms, and its first pair joins the two. So the parts fall into groups joined by the pairs
     * reached, and each group pairs off its own bonds and atoms; a group of one molecule only gets
     * nothing. Pairing off within each group pairs off no more than over all of them.
     */
    private long groupedBound(int rankLimit, int open, int newBonds) {
      int queryOpen = spread(queryLabels, query, true, rankLimit, newBonds);
      spread(targetLabels, target, false, rankLimit, newBonds);
      int work = reach.count;
      for (int i = 0; i < reach.count && newBonds > 0; i++)
        UnionFind.join(partParent, queryAtom(reach.pairs[i]),
            queryCount + targetAtom(reach.pairs[i]));
      // The open atoms of each group by a counting sort on its root, the query atoms first.
      int groups = 0, openCount = queryOpen + targetLabels.openCount;
      work += openCount;
      clock.tick(work);
      for (int index = 0; index < openCount; index++) {
        int root = UnionFind.root(partParent, openItem(index, queryOpen));
        if (groupMark[root] != markStamp) {
          groupMark[root] = markStamp;
          groupOf[root] = groups;
          groupStart[++groups] = 0;
          groupSides[groups - 1] = 0;
        }
        groupStart[groupOf[root] + 1]++;
        groupSides[groupOf[root]] |= index < queryOpen ? 1 : 2;
      }
      groupStart[0] = 0;
      for (int group = 0; group < groups; group++) groupStart[group + 1] += groupStart[group];
      for (int index = 0; index < openCount; index++) {
        int item = openItem(index, queryOpen);
        groupItems[groupStart[groupOf[UnionFind.root(partParent, item)]]++] = item;
      }
      int added = open, newAtoms = 0;
      for (int group = 0, from = 0; group < groups; group++) {
        int to = groupStart[group];
        if (groupSides[group] == 3) {
          // The query atoms of the group come first.
          int split = from;
          while (groupItems[split] < queryCount) split++;
          countAtoms(queryLabels, query, true, groupItems, from, split - from, 0, true);
          countAtoms(targetLabels, target, false, groupItems, split, to - split, queryCount, true);
          long counted = pairCounted();
          added += (int) (counted >>> 32);
          newAtoms += (int) counted;
        }
        from = to;
      }
      if (++markStamp == Integer.MAX_VALUE) {
        Arrays.fill(groupMark, 0);
        markStamp = 1;
      }
      queryLabels.clearOpen();
      targetLabels.clearOpen();
      return capped(depth + Math.min(newAtoms, added), commonBonds + added);
    }

    /**
     * The union-find item of the open atom at {@code index}: the query atoms first, and then the
     * target atoms after {@link #queryCount}.
     */
    private int openItem(int index, int queryOpen) {
      return index < queryOpen ? queryLabels.openAtoms[index]
                               : queryCount + targetLabels.openAtoms[index - queryOpen];
    }

    /**
     * Opens the atoms of {@code side} in layers: first those of the pairs reached from the mapped
     * pairs, and then their neighbours through open atoms and bonds that have a label, up to
     * {@code newBonds} layers. Atoms linked by such bonds are joined in the union-find of the
     * groups. Each open query atom adds the targets of its pairs to the reachable targets, which
     * are the target atoms that may open. Gives the number of atoms opened.
     */
    private int spread(Labels labels, CommonBondGraphs.Side side, boolean onQuery, int rankLimit,
        int newBonds) {
      int stamp = labels.nextStamp(), offset = onQuery ? 0 : queryCount, work = 0;
      for (int i = 0; i < depth; i++)
        labels.reached[onQuery ? levels[i].queryAtom : levels[i].targetAtom] = stamp;
      // The first layer holds the atoms of the pairs reached from the mapped pairs.
      for (int i = 0; i < reach.count && newBonds > 0; i++, work++)
        reach(labels, onQuery, onQuery ? queryAtom(reach.pairs[i]) : targetAtom(reach.pairs[i]),
            rankLimit, stamp, offset);
      for (int index = 0, layerEnd = labels.openCount, distance = 1; index < labels.openCount;
           index++) {
        if (index == layerEnd) {
          distance++;
          layerEnd = labels.openCount;
        }
        int atom = labels.openAtoms[index];
        work += 1 + side.degree(atom);
        for (int place = side.start[atom]; place < side.start[atom + 1]; place++) {
          int other = side.neighbour[place];
          if (labels.slot[place] < 0)
            continue;
          // An atom first met from the last layer would lie a bond too far.
          if (labels.reached[other] != stamp && distance < newBonds)
            reach(labels, onQuery, other, rankLimit, stamp, offset);
          if (labels.isOpen[other])
            UnionFind.join(partParent, offset + atom, offset + other);
        }
      }
      clock.tick(work);
      return labels.openCount;
    }

    /**
     * Marks {@code atom}, met in {@link #spread}, as reached, and opens it if it may be mapped: a
     * free query atom ranked up to {@code rankLimit} with a pair not excluded, or a reachable
     * target atom.
     */
    private void reach(Labels labels, boolean onQuery, int atom, int rankLimit, int stamp,
        int offset) {
      if (labels.reached[atom] == stamp)
        return;
      labels.reached[atom] = stamp;
      boolean mayOpen = onQuery
          ? queryToTarget[atom] < 0 && queryRank[atom] <= rankLimit && markTargets(atom)
          : has(reachableTargets, 0, atom);
      if (!mayOpen)
        return;
      labels.markOpen(atom);
      partParent[offset + atom] = offset + atom;
    }

    /**
     * Counts the {@code count} open atoms of {@code side} at {@code atoms[from]} onwards, each
     * less {@code offset}, by slot and by their number of open bonds in it. One with an open bond
     * or a link is also counted by class and, where a class holds kinds that do not all match, by
     * kind. Each bond to a higher open atom of a label whose fine labels do not all match is
     * counted by fine label. With {@code listed}, the labels and the classes counted are listed
     * for {@link #pairCounted}.
     */
    private void countAtoms(Labels labels, CommonBondGraphs.Side side, boolean onQuery, int[] atoms,
        int from, int count, int offset, boolean listed) {
      int[] linkStamp = onQuery ? reach.queryBondStamp : reach.targetBondStamp;
      int work = count;
      for (int index = from; index < from + count; index++) {
        int atom = atoms[index] - offset, slots = 0;
        boolean link = false;
        work += side.degree(atom);
        for (int place = side.start[atom]; place < side.start[atom + 1]; place++) {
          int slot = labels.slot[place], other = side.neighbour[place];
          if (slot < 0 || !labels.isOpen[other]) {
            link |= linkStamp[side.neighbourBond[place]] == reach.stamp;
            continue;
          }
          if (slotCount[slot]++ == 0)
            slotsAt[slots++] = slot;
          if (other > atom && !labelGroups.complete[slot / 2])
            labels.countFine(side.fineLabel[side.neighbourBond[place]]);
        }
        for (int slotIndex = 0; slotIndex < slots; slotIndex++) {
          int slot = slotsAt[slotIndex];
          labels.histogram[slot * histogramWidth + slotCount[slot]]++;
          slotCount[slot] = 0;
          if (listed && !labelListed[slot / 2]) {
            labelListed[slot / 2] = true;
            countedLabels[labelsCounted++] = slot / 2;
          }
        }
        if (slots > 0 || link) {
          int symbol = labels.symbol[atom];
          if (labels.bySymbol[symbol]++ == 0 && listed && !symbolListed[symbol]) {
            symbolListed[symbol] = true;
            countedSymbols[symbolsCounted++] = symbol;
          }
          if (!atomGroups.complete[symbol])
            labels.countKind(side.atomKind[atom]);
        }
      }
      clock.tick(work);
    }

    /**
     * The new bonds, in the high int, and the new atoms that the atoms counted can pair off: the
     * bonds of each label and the atoms of each class. Clears the counts by slot and by class.
     */
    private long pairAll() {
      int bonds = 0, atoms = 0;
      clock.tick(slotSymbol.length * histogramWidth + queryLabels.bySymbol.length);
      for (int slot = 0; slot < slotSymbol.length; slot += 2) {
        int label = slot / 2;
        int pairs = slotSymbol[slot + 1] < 0 ? slotBonds(slot) / 2
                                             : Math.min(slotBonds(slot), slotBonds(slot + 1));
        if (!labelGroups.complete[label] && pairs > 0)
          pairs = Math.min(pairs, labelGroups.bound(label, queryLabels.fineCount, targetLabels.fineCount));
        bonds += pairs;
      }
      for (int symbol = 0; symbol < queryLabels.bySymbol.length; symbol++) {
        int pairs = Math.min(queryLabels.bySymbol[symbol], targetLabels.bySymbol[symbol]);
        if (!atomGroups.complete[symbol] && pairs > 0)
          pairs = atomGroups.bound(symbol, queryLabels.kindCount, targetLabels.kindCount);
        atoms += pairs;
        queryLabels.bySymbol[symbol] = targetLabels.bySymbol[symbol] = 0;
      }
      return (long) bonds << 32 | atoms;
    }

    /**
     * The new bonds, in the high int, and the new atoms that the atoms listed by
     * {@link #countAtoms} can pair off: the bonds of each label listed and the atoms of each
     * class listed. Clears the counts.
     */
    private long pairCounted() {
      int bonds = 0, atoms = 0;
      clock.tick(labelsCounted * histogramWidth + symbolsCounted);
      for (int index = 0; index < labelsCounted; index++) {
        int label = countedLabels[index], slot = 2 * label;
        int pairs = slotSymbol[slot + 1] < 0 ? slotBonds(slot) / 2
                                             : Math.min(slotBonds(slot), slotBonds(slot + 1));
        if (!labelGroups.complete[label] && pairs > 0)
          pairs = Math.min(pairs, labelGroups.bound(label, queryLabels.fineCount, targetLabels.fineCount));
        bonds += pairs;
        labelListed[label] = false;
      }
      for (int index = 0; index < symbolsCounted; index++) {
        int symbol = countedSymbols[index];
        int pairs = Math.min(queryLabels.bySymbol[symbol], targetLabels.bySymbol[symbol]);
        if (!atomGroups.complete[symbol] && pairs > 0)
          pairs = atomGroups.bound(symbol, queryLabels.kindCount, targetLabels.kindCount);
        atoms += pairs;
        queryLabels.bySymbol[symbol] = targetLabels.bySymbol[symbol] = 0;
        symbolListed[symbol] = false;
      }
      labelsCounted = symbolsCounted = 0;
      queryLabels.clearCounts();
      targetLabels.clearCounts();
      return (long) bonds << 32 | atoms;
    }

    /** The most new bonds of a slot, pairing the open atoms with the most bonds of it on either side. This clears the counts of the slot. */
    private int slotBonds(int slot) {
      int sum = 0;
      for (int openBonds = histogramWidth - 1, queryAtoms = 0, targetAtoms = 0; openBonds > 0;
           openBonds--) {
        queryAtoms += queryLabels.histogram[slot * histogramWidth + openBonds];
        targetAtoms += targetLabels.histogram[slot * histogramWidth + openBonds];
        queryLabels.histogram[slot * histogramWidth + openBonds] = 0;
        targetLabels.histogram[slot * histogramWidth + openBonds] = 0;
        sum += Math.min(queryAtoms, targetAtoms);
      }
      return sum;
    }

    /** The slot at each place, or -1 for none: twice the label of its bond, plus one at the second class. */
    private int[] slots(CommonBondGraphs.Side side) {
      int[] slot = new int[side.neighbour.length];
      for (int atom = 0; atom < side.atomCount; atom++) {
        for (int place = side.start[atom]; place < side.start[atom + 1]; place++) {
          int label = side.label[side.neighbourBond[place]];
          slot[place] = label < 0 ? -1
              : 2 * label + (side.atomClass[atom] == slotSymbol[2 * label] ? 0 : 1);
        }
      }
      return slot;
    }
  }

  /**
   * The labels of one molecule in the bond-label bound, and the open atoms of a state with their
   * counts. Package-private so that tests can check the stamp wrap.
   */
  static final class Labels {

    /** The class of each atom, and the slot at each place of the neighbour lists, see {@link EdgeBound}. */
    final int[] symbol, slot;
    /** The open atoms by slot and by the number of open bonds of that slot. */
    final int[] histogram;
    /** The open atoms by class that may get a new bond. */
    final int[] bySymbol;
    /** The open atoms by kind, and the bonds between open atoms by fine label, where these are counted. */
    final int[] kindCount, fineCount;
    /** The kinds and fine labels counted, the first {@link #kinds} and {@link #fines}, to clear. */
    private final int[] countedKinds, countedFines;
    /** The open atoms, the first {@link #openCount} of the array. */
    final int[] openAtoms;
    /** Whether each atom is open. */
    final boolean[] isOpen;
    /** The stamp of the atoms reached from the mapped atoms, see {@link EdgeBound#spread}. */
    final int[] reached;
    /** The number of open atoms, the stamp of the last reach, and the numbers of kinds and fine labels counted. */
    int openCount, stamp, kinds, fines;

    /**
     * The labels of {@code side}, with the {@code slot} at each place, a histogram of
     * {@code histogramSize}, a count for each of the {@code symbols}, and counts for its
     * {@code atomKinds} kinds and {@code fineLabels} fine labels.
     */
    Labels(CommonBondGraphs.Side side, int[] slot, int histogramSize, int symbols, int atomKinds,
        int fineLabels) {
      symbol = side.atomClass;
      this.slot = slot;
      histogram = new int[histogramSize];
      bySymbol = new int[symbols];
      kindCount = new int[atomKinds];
      fineCount = new int[fineLabels];
      countedKinds = new int[atomKinds];
      countedFines = new int[fineLabels];
      openAtoms = new int[symbol.length];
      isOpen = new boolean[symbol.length];
      reached = new int[symbol.length];
    }

    /** Marks {@code atom}, which is not yet open, as open. */
    void markOpen(int atom) {
      isOpen[atom] = true;
      openAtoms[openCount++] = atom;
    }

    /** Counts one open atom of {@code kind}. */
    void countKind(int kind) {
      if (kindCount[kind]++ == 0)
        countedKinds[kinds++] = kind;
    }

    /** Counts one bond between open atoms of fine label {@code fine}. */
    void countFine(int fine) {
      if (fineCount[fine]++ == 0)
        countedFines[fines++] = fine;
    }

    /**
     * Starts a new reach, clearing every mark when the stamp wraps to zero, as the unsigned
     * 32-bit stamp of the C++ search does after the same number of reaches.
     */
    int nextStamp() {
      if (++stamp == 0) {
        Arrays.fill(reached, 0);
        stamp = 1;
      }
      return stamp;
    }

    /** Clears the counts by kind and fine label. */
    void clearCounts() {
      while (kinds > 0) kindCount[countedKinds[--kinds]] = 0;
      while (fines > 0) fineCount[countedFines[--fines]] = 0;
    }

    /** Closes every open atom and clears the counts by kind and fine label. */
    void clearOpen() {
      for (int openIndex = 0; openIndex < openCount; openIndex++) isOpen[openAtoms[openIndex]] = false;
      openCount = 0;
      clearCounts();
    }
  }

  /**
   * The embeddings of a connected query: maps of every query atom that take each query bond onto
   * a compatible target bond, with other target bonds allowed. They are found by backtracking
   * with forward checking. Placing an atom narrows the candidates of each unplaced neighbour to
   * the free neighbours of its target, and these neighbours must be able to have distinct
   * targets. The next atom placed is the one with the fewest free candidates. There is no
   * recursion, as the depth is the number of query atoms.
   */
  final class Embedding {

    /** The target atom each query atom is placed on, and the query atom placed on each target atom, or -1. */
    private final int[] placedOn, placedAtom;
    /**
     * The candidates of each query atom in the order they are tried, the target atoms it may
     * still have. They are set for the first atom at the start, and for any other once a
     * neighbour is placed. Until then they are null.
     */
    private final int[][] candidates;
    /** The number of candidates of each unplaced query atom that are free. */
    private final int[] freeCount;
    /** The unplaced query atoms that have candidates, the first {@link #frontierSize}, in any order. */
    private final int[] frontier;
    /** The number of atoms on the frontier. */
    private int frontierSize;
    /** Room for the candidates of a neighbour while they are found, before they are copied. */
    private final int[] scratch;
    /** Whether only embeddings with their twins in order are tried, and the first that counts ends the search. */
    final boolean prunesTwins;

    /** Starts from the query atom with the fewest targets it may have, then most bonds, then lowest index. */
    private Embedding(long[] candidateRows, boolean prunesTwins) {
      this.prunesTwins = prunesTwins;
      if (prunesTwins)
        findTwins();
      placedOn = new int[queryCount];
      placedAtom = new int[targetCount];
      Arrays.fill(placedOn, -1);
      Arrays.fill(placedAtom, -1);
      candidates = new int[queryCount][];
      freeCount = new int[queryCount];
      frontier = new int[queryCount];
      scratch = new int[target.maxDegree()];
      int root = 0, fewest = Integer.MAX_VALUE;
      for (int queryAtom = 0; queryAtom < queryCount; queryAtom++) {
        int count = 0;
        for (int i = 0; i < words; i++) count += Long.bitCount(candidateRows[queryAtom * words + i]);
        if (count < fewest || count == fewest && query.degree(queryAtom) > query.degree(root)) {
          root = queryAtom;
          fewest = count;
        }
      }
      int[] targets = new int[fewest];
      int row = root * words;
      int filled = 0;
      for (int targetAtom = CommonBondGraphs.nextBit(candidateRows, row, words, 0); targetAtom >= 0;
           targetAtom = CommonBondGraphs.nextBit(candidateRows, row, words, targetAtom + 1))
        targets[filled++] = targetAtom;
      candidates[root] = targets;
      freeCount[root] = fewest;
      frontier[frontierSize++] = root;
    }

    /**
     * The embeddings, at most {@code maxMappings}, held as the mappings as they are found. Under
     * conditions only an embedding that meets them with every query bond counts, as only it has
     * the score of an embedding; if none does, no mapping meets them with every query bond. When
     * it prunes twins, the first embedding that counts ends it. Called once only.
     */
    List<int[]> all(int maxMappings) {
      mappings = new ArrayList<>();
      long whole = CommonBondGraphs.score(queryCount, query.bondCount);
      Placement[] placements = new Placement[queryCount];
      placements[0] = new Placement().enter(frontier[0], 0);
      for (int placed = 1; placed > 0;) {
        clock.tick();
        Placement top = placements[placed - 1];
        if (placedOn[top.atom] >= 0)
          unplace(top);
        int targetAtom = nextTarget(top);
        if (targetAtom < 0) {
          placed--;
        } else if (place(top, targetAtom)) {
          int next = choose();
          if (next >= 0) {
            if (placements[placed] == null)
              placements[placed] = new Placement();
            placements[placed++].enter(frontier[next], next);
          } else if (placed == queryCount) {
            if (conditions == null)
              mappings.add(placedOn.clone());
            else
              conditions.select(placedOn, whole, maxMappings - mappings.size(), mappings);
            if (mappings.size() == maxMappings || prunesTwins && !mappings.isEmpty())
              return mappings;
          }
        }
      }
      return mappings;
    }

    /** The next free candidate of the atom of {@code placement}, or -1 once it has tried them all. */
    private int nextTarget(Placement placement) {
      int[] targets = candidates[placement.atom];
      while (placement.nextCandidate < targets.length) {
        int targetAtom = targets[placement.nextCandidate++];
        if (placedAtom[targetAtom] < 0 && !(prunesTwins && twinsOutOfOrder(placement.atom, targetAtom)))
          return targetAtom;
      }
      return -1;
    }

    /**
     * Whether placing {@code queryAtom} on {@code targetAtom} puts twins out of order: a placed
     * twin of the query atom lies on the other side of the target atom, by index, than of the
     * query atom, or a twin of the target atom in use has its query atom on the other side of
     * the query atom. Exchanging twins of either molecule takes an embedding onto one that meets
     * the conditions alike, and the exchanges take every embedding onto one with all its twins
     * in order: the least of them, each read as the targets of the query atoms in turn.
     */
    private boolean twinsOutOfOrder(int queryAtom, int targetAtom) {
      for (int twin = queryTwin[queryAtom]; twin >= 0 && twin != queryAtom; twin = queryTwin[twin]) {
        clock.tick();
        if (placedOn[twin] >= 0 && (twin < queryAtom) != (placedOn[twin] < targetAtom))
          return true;
      }
      for (int twin = targetTwin[targetAtom]; twin >= 0 && twin != targetAtom; twin = targetTwin[twin]) {
        clock.tick();
        if (placedAtom[twin] >= 0 && (twin < targetAtom) != (placedAtom[twin] < queryAtom))
          return true;
      }
      return false;
    }

    /** Places the atom of {@code placement} on {@code targetAtom} and records the changes in {@code placement}, to be undone on failure. */
    private boolean place(Placement placement, int targetAtom) {
      countNode();
      int queryAtom = placement.atom;
      frontier[placement.frontierPlace] = frontier[--frontierSize];
      placedOn[queryAtom] = targetAtom;
      placedAtom[targetAtom] = queryAtom;
      take(placement, targetAtom);
      int from = target.start[targetAtom], degree = target.degree(targetAtom), open = 0, narrowed = 0;
      clock.tick(1 + query.degree(queryAtom));
      for (int queryPlace = query.start[queryAtom]; queryPlace < query.start[queryAtom + 1];
           queryPlace++)
        if (placedOn[query.neighbour[queryPlace]] < 0)
          open++;
      // The candidates of each neighbour as places among the target's neighbours, for the check at the end.
      long[] rows = open > 1 && degree <= Long.SIZE ? new long[open] : null;
      for (int queryPlace = query.start[queryAtom]; queryPlace < query.start[queryAtom + 1];
           queryPlace++) {
        int queryNeighbour = query.neighbour[queryPlace], count = 0;
        if (placedOn[queryNeighbour] >= 0)
          continue;
        // A neighbour with candidates has another placed neighbour, so it closes a ring.
        boolean closes = candidates[queryNeighbour] != null;
        clock.tick(1 + (closes ? degree * (1 + query.degree(queryNeighbour)) : degree));
        long bits = 0;
        for (int offset = 0; offset < degree; offset++) {
          int candidateTarget = target.neighbour[from + offset];
          if (placedAtom[candidateTarget] < 0
              && graph.bondsMatch(query.neighbourBond[queryPlace], target.neighbourBond[from + offset])
              && mayMap(queryNeighbour, candidateTarget)
              && (!closes || bondsFit(queryNeighbour, candidateTarget, queryAtom))) {
            scratch[count++] = candidateTarget;
            bits |= 1L << offset;
          }
        }
        if (count == 0)
          return false;
        if (!closes)
          frontier[frontierSize++] = queryNeighbour;
        placement.remember(queryNeighbour, candidates[queryNeighbour], freeCount[queryNeighbour]);
        candidates[queryNeighbour] = Arrays.copyOf(scratch, count);
        freeCount[queryNeighbour] = count;
        if (rows != null)
          rows[narrowed++] = bits;
      }
      return rows == null || hasDistinctTargets(rows, narrowed, 1, clock);
    }

    /** Takes the atom of {@code placement} off its target, and puts back everything that placing it changed. */
    private void unplace(Placement placement) {
      int queryAtom = placement.atom;
      while (placement.changedCount > 0) {
        int changed = --placement.changedCount, changedAtom = placement.changedAtoms[changed];
        // A neighbour that had no candidates came on the frontier last.
        if (placement.oldCandidates[changed] == null)
          frontierSize--;
        candidates[changedAtom] = placement.oldCandidates[changed];
        freeCount[changedAtom] = placement.oldFreeCounts[changed];
      }
      placedAtom[placedOn[queryAtom]] = -1;
      placedOn[queryAtom] = -1;
      frontier[frontierSize++] = frontier[placement.frontierPlace];
      frontier[placement.frontierPlace] = queryAtom;
    }

    /**
     * The frontier place of the atom to place next, with the fewest free candidates, then the
     * most bonds, then the lowest index. It is -1 if the frontier is empty or an atom on it has
     * no free candidate left.
     */
    private int choose() {
      clock.tick(frontierSize);
      int pick = -1;
      for (int frontierIndex = 0; frontierIndex < frontierSize; frontierIndex++) {
        int candidate = frontier[frontierIndex], picked = pick < 0 ? -1 : frontier[pick];
        if (freeCount[candidate] == 0)
          return -1;
        if (picked < 0 || freeCount[candidate] < freeCount[picked]
            || freeCount[candidate] == freeCount[picked]
                && (query.degree(candidate) > query.degree(picked)
                    || query.degree(candidate) == query.degree(picked) && candidate < picked))
          pick = frontierIndex;
      }
      return pick;
    }

    /**
     * Takes target atom {@code targetAtom} from the free candidates of the unplaced atoms that
     * have it, remembering their counts in {@code placement}. Such an atom may map onto
     * {@code targetAtom} with compatible bonds to the targets of all its placed neighbours, and
     * is counted from the first of these.
     */
    private void take(Placement placement, int targetAtom) {
      int work = 1 + target.degree(targetAtom);
      for (int targetPlace = target.start[targetAtom]; targetPlace < target.start[targetAtom + 1];
           targetPlace++) {
        int placedNeighbour = placedAtom[target.neighbour[targetPlace]];
        if (placedNeighbour < 0)
          continue;
        for (int queryPlace = query.start[placedNeighbour];
             queryPlace < query.start[placedNeighbour + 1]; queryPlace++) {
          int queryNeighbour = query.neighbour[queryPlace];
          work++;
          if (placedOn[queryNeighbour] >= 0
              || !graph.bondsMatch(query.neighbourBond[queryPlace], target.neighbourBond[targetPlace])
              || !mayMap(queryNeighbour, targetAtom))
            continue;
          work += query.degree(queryNeighbour);
          if (firstPlacedNeighbour(queryNeighbour) == placedNeighbour
              && bondsFit(queryNeighbour, targetAtom, placedNeighbour)) {
            placement.remember(queryNeighbour, candidates[queryNeighbour], freeCount[queryNeighbour]);
            freeCount[queryNeighbour]--;
          }
        }
      }
      clock.tick(work);
    }

    /** Whether {@code queryAtom} may map onto {@code targetAtom}: they match, and the target has as many bonds. */
    private boolean mayMap(int queryAtom, int targetAtom) {
      return has(allowed, queryAtom, targetAtom) && target.degree(targetAtom) >= query.degree(queryAtom);
    }

    /** Whether {@code targetAtom} has compatible bonds to the images of the placed neighbours of {@code queryAtom} other than {@code skipped}. */
    private boolean bondsFit(int queryAtom, int targetAtom, int skipped) {
      for (int queryPlace = query.start[queryAtom]; queryPlace < query.start[queryAtom + 1];
           queryPlace++) {
        int queryNeighbour = query.neighbour[queryPlace], neighbourImage = placedOn[queryNeighbour];
        if (neighbourImage < 0 || queryNeighbour == skipped)
          continue;
        // Search the shorter list of neighbours, as one may be a hub.
        int targetBond = target.degree(targetAtom) <= target.degree(neighbourImage)
            ? target.bondBetween(targetAtom, neighbourImage)
            : target.bondBetween(neighbourImage, targetAtom);
        if (targetBond < 0 || !graph.bondsMatch(query.neighbourBond[queryPlace], targetBond))
          return false;
      }
      return true;
    }

    /** The first placed neighbour of query atom {@code queryAtom}, or -1 if none is placed. */
    private int firstPlacedNeighbour(int queryAtom) {
      for (int i = query.start[queryAtom]; i < query.start[queryAtom + 1]; i++)
        if (placedOn[query.neighbour[i]] >= 0)
          return query.neighbour[i];
      return -1;
    }
  }

  /**
   * A placed atom of the whole-query search, one for each depth. It holds the place the atom
   * had on the frontier, the next of its candidates to try, and the candidates and counts that
   * placing it changed as they were before, so that taking it off puts them back.
   */
  private static final class Placement {

    /** The placed atom, the place it had on the frontier, and the next of its candidates to try. */
    int atom, frontierPlace, nextCandidate;
    /** The number of changes recorded. */
    int changedCount;
    /** The atoms whose candidates or counts changed, the first {@link #changedCount}, and their old counts. */
    int[] changedAtoms = new int[4], oldFreeCounts = new int[4];
    /** The old candidates of each changed atom, null where it had none. */
    int[][] oldCandidates = new int[4][];

    /** Starts placing {@code atom}, from place {@code frontierPlace} of the frontier, at its first candidate. */
    Placement enter(int atom, int frontierPlace) {
      this.atom = atom;
      this.frontierPlace = frontierPlace;
      nextCandidate = 0;
      return this;
    }

    /** Remembers that {@code atom} had these candidates, and this count of free ones. */
    void remember(int atom, int[] candidates, int freeCount) {
      if (changedCount == changedAtoms.length) {
        changedAtoms = Arrays.copyOf(changedAtoms, 2 * changedCount);
        oldFreeCounts = Arrays.copyOf(oldFreeCounts, 2 * changedCount);
        oldCandidates = Arrays.copyOf(oldCandidates, 2 * changedCount);
      }
      changedAtoms[changedCount] = atom;
      oldCandidates[changedCount] = candidates;
      oldFreeCounts[changedCount++] = freeCount;
    }
  }

  /**
   * Whether each of the first {@code count} rows can have a column of its own. Each row is
   * {@code words} longs of {@code rows}, and its columns are their bits. A row takes a free
   * column, or else follows an augmenting path. There is no recursion, as there may be a row for
   * each query atom.
   */
  private static boolean hasDistinctTargets(long[] rows, int count, int words,
      CommonBondGraphs.Clock clock) {
    int columns = words << 6;
    int[] owner = new int[columns], visited = new int[columns], path = new int[count],
          resume = new int[count];
    long[] freeColumns = new long[words];
    Arrays.fill(owner, -1);
    Arrays.fill(freeColumns, -1L);
    for (int newRow = 0; newRow < count; newRow++) {
      clock.tick(words);
      int column = -1, step = 0;
      for (int word = 0; word < words && column < 0; word++) {
        long freeBits = rows[newRow * words + word] & freeColumns[word];
        if (freeBits != 0)
          column = (word << 6) + Long.numberOfTrailingZeros(freeBits);
      }
      path[0] = newRow;
      resume[0] = 0;
      while (column < 0) {
        clock.tick();
        int row = path[step];
        column = CommonBondGraphs.nextBit(rows, row * words, words, resume[step]);
        while (column >= 0 && visited[column] == newRow + 1) {
          clock.tick();
          column = CommonBondGraphs.nextBit(rows, row * words, words, column + 1);
        }
        if (column < 0) {
          if (--step < 0)
            return false;
          continue;
        }
        resume[step] = column + 1;
        visited[column] = newRow + 1;
        if (owner[column] >= 0) {
          path[++step] = owner[column];
          resume[step] = 0;
          column = -1;
        }
      }
      // Each row on the path takes the column it reached. The last row takes a free column.
      freeColumns[column >>> 6] &= ~(1L << column);
      owner[column] = path[step];
      for (int earlier = step - 1; earlier >= 0; earlier--) owner[resume[earlier] - 1] = path[earlier];
    }
    return true;
  }

  /** Whether bit {@code bit} of row {@code row} is set, in rows of {@link #words} longs as in {@link #allowed}. */
  private boolean has(long[] rows, int row, int bit) {
    return (rows[row * words + (bit >>> 6)] & 1L << bit) != 0;
  }

  /** Sets bit {@code bit} of row {@code row}. */
  private void set(long[] rows, int row, int bit) {
    rows[row * words + (bit >>> 6)] |= 1L << bit;
  }

  /** Clears bit {@code bit} of row {@code row}. */
  private void clear(long[] rows, int row, int bit) {
    rows[row * words + (bit >>> 6)] &= ~(1L << bit);
  }
}
