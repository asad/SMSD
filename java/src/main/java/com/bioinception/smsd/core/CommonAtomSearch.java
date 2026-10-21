/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
package com.bioinception.smsd.core;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;

/**
 * The maximum common induced subgraph (MCIS) search in fragments on the int graphs of a
 * {@link CommonBondGraphs}: the mapping from the query to the target with the most atoms and then
 * the most bonds, in which two mapped query atoms are bonded exactly when their images are, and
 * then by compatible bonds. The common subgraph is any number of fragments of at least
 * {@link #minAtoms} atoms anywhere in either molecule.
 *
 * <p>The unmapped atoms that may still be mapped are kept in {@link Blocks}: a block holds the
 * query and target atoms of one atom class with the same bonds, by bond class, to the mapped
 * atoms. An atom can be mapped only onto an atom of its own block, so the atoms a mapping can still
 * gain are at most those each block can pair off, which bounds a state. Mapping a pair splits each
 * block by the bonds to the new pair, and an atom whose block has nothing left on the other side
 * can no longer be mapped. The blocks and their bound are the label classes of C. McCreesh, P.
 * Prosser and J. Trimble, "A Partitioning Algorithm for Maximum Common Subgraph Problems", Proc.
 * IJCAI 2017, 712-719, with atom classes as vertex labels and bond classes as edge labels. A
 * depth-first branch and bound grows the mapping one fragment at a time
 * on an explicit stack of {@link Level}s, one for each mapped pair. A level maps a query atom of
 * the smallest block onto each target of the block in turn, and then leaves that atom unmapped for
 * the rest of the level, so each mapping is found once. While a fragment is open, the atom is one
 * bonded to it; once none is left, the fragment is closed and the atom is one bonded to no mapped
 * atom, the root of a new fragment. A pair that some best extension of a state maps is the only
 * one tried ({@link #forcedPair}). Greedy seeds give a score to beat. The best score is first
 * sought level by level of common atoms, down from an upper bound, while the levels grow in cost;
 * one search from the best mapping known then finishes the work. Every loop ticks the
 * {@link CommonBondGraphs.Clock}, which throws {@link CommonBondGraphs.Stop} when the time limit is
 * reached or the search is cancelled; the best mapping found by then is kept ({@link #found}).
 */
final class CommonAtomSearch {

  /** The most greedy seeds grown in each round. */
  private static final int SEEDS = 64;
  /** The most free targets a seed looks at next to a hub. */
  private static final int HUB_TARGETS = 64;
  /** The levels of likeness of two atoms, from no radius alike to all of them. */
  private static final int LIKENESS_LEVELS = CommonBondGraphs.RADII + 1;
  /** The most pairs of a linked block that {@link #blockComplete} checks one by one. */
  private static final int COMPLETE_PAIRS = 64;
  /** The first capacity of the block pool, which grows as needed. */
  private static final int INITIAL_BLOCKS = 64;
  /** The least growth in cost from one level of common atoms to the next that keeps the levels going, see {@link #collectByAtoms}. */
  private static final int LEVEL_GROWTH = 2;
  /** An empty list of atoms, shared by every level that has nothing to try. */
  private static final int[] NONE = new int[0];

  /** The graphs of both molecules, with their compatible pairs and classes. */
  private final CommonBondGraphs graph;
  /** The query and target graphs of {@link #graph}. */
  private final CommonBondGraphs.Side query, target;
  /** Counts the work of every loop and ends the search with {@link CommonBondGraphs.Stop} when told to. */
  private final CommonBondGraphs.Clock clock;
  /** The fewest atoms of a fragment, one or two. */
  private final int minAtoms;
  /** The numbers of query and target atoms, and of longs in a row of target atoms. */
  private final int queryCount, targetCount, words;
  /** The compatible atom pairs: a row for each query atom, with a bit for each target atom. */
  private final long[] compatible;
  /** Whether every query bond matches every target bond of its class, so that blocks need no bond check. */
  private final boolean bondClassesMatch;
  /** The mapping both ways, -1 where unmapped. */
  private final int[] queryToTarget, targetToQuery;
  /** The levels: the first holds nothing mapped, and level d the d-th pair mapped. */
  private final Level[] levels;
  /** The number of mapped pairs, which is the index of the top level, and the bonds between them. */
  private int depth, commonBonds;
  /** The serial of the last level set up, see {@link Level#serial}; a long, so that it never wraps. */
  private long serial;
  /** The serial of the level at which each query atom was last tried, see {@link #twinTried}. */
  private final long[] triedSerial;
  /** The atoms of the blocks: each block holds a run of each array, its live atoms first. */
  private final int[] queryAtoms, targetAtoms;
  /** The place of each atom in {@link #queryAtoms} and {@link #targetAtoms}. */
  private final int[] queryPlace, targetPlace;
  /** The blocks of the levels on the stack. */
  private final Blocks blocks = new Blocks();
  /** The class of the bond from the pair being mapped to each atom next to it, valid where the stamp is current. */
  private final int[] queryLinkClass, targetLinkClass;
  private final long[] queryLinkStamp, targetLinkStamp;
  /** The stamp of the pair being mapped, see {@link #queryLinkStamp}. Stamps are longs, so they never wrap. */
  private long linkStamp;
  /** Scratch of {@link #refine}: the live atoms bonded to the new pair on each side, and their blocks. */
  private final int[] linkedQuery, linkedQueryBlock, linkedTarget, linkedTargetBlock;
  /** Scratch of {@link #liveBonds}: the atoms of live blocks by stamp, and the live block of each. */
  private final long[] queryLive, targetLive;
  private final int[] queryBlock, targetBlock;
  /** The stamp of the last count of {@link #liveBonds}. */
  private long liveStamp;
  /** Scratch of {@link #liveBonds}: the bonds by block and bond class on each side, and the keys counted. */
  private int[] queryPairCount = NONE, targetPairCount = NONE, touchedPairs = NONE;
  /** Scratch of {@link #blockPairs}: the atoms of each kind in a block. */
  private final int[] queryKindCount, targetKindCount;
  /** Scratch of {@link #exactPairs}: the targets of each query atom, and the query atom of each target. */
  private final long[] matchRows = new long[Long.SIZE];
  private final int[] matchOwner = new int[Long.SIZE];
  /** Scratch of {@link #exactPairs}: the targets an augmenting path has passed. */
  private long matchSeen;
  /** The forced pair found by {@link #forcedPair}, and its block. */
  private int forcedBlock, forcedQuery, forcedTarget;
  /** The atoms by rank, lowest first, and the rank of each atom, see {@link #rank}. */
  private int[] queryOrder, queryRank, targetOrder, targetRank;
  /** The cell of each atom: atoms that look alike share a cell. */
  private int[] queryCell, targetCell;
  /**
   * The runs of target atoms of one cell in target order, from {@code targetRunStart[run]} up to
   * the next start, and in each the first place whose atom may still be a seed root.
   */
  private int[] targetRunStart, targetRunFree;
  /** The {@link #seedRoots} call in which each query cell had no run left to take a root from. */
  private long[] queryCellClosed;
  private long rootCall;
  /** The signatures of the atoms at each radius, which tell how alike two atoms are. */
  private long[][] querySignature, targetSignature;
  /** The next atom alike to each terminal atom, see {@link #twins}, and the first atom of each ring of target twins. */
  private int[] queryTwin, targetTwin, targetTwinFirst;
  /** Scratch of {@link #targets}: the stamp of each ring of twins seen, and the sort keys of the targets. */
  private final long[] twinSeen, targetKeys;
  /** The stamp of the targets being listed, see {@link #twinSeen}. */
  private long twinStamp;
  /**
   * The taken targets that lead the neighbours of each target atom, and the seed they were counted
   * for, see {@link #leadingTaken}.
   */
  private final int[] leadingTaken, leadingRound;
  /** The likest level of each query atom with any target, -1 until worked out, see {@link #likestLevel}. */
  private final int[] likestLevel;
  /** The best score so far, see {@link #score}, and a score that no mapping beats. */
  private long bestScore, upperBound;
  /** The mapping with the best score so far, or null while none beats the threshold of the pass. */
  private int[] best;
  /**
   * The mapping a stopped search falls back on while the pass under way holds none: the seed, or
   * the best mapping a level of the search met below its threshold; null until one is found. See
   * {@link #found}.
   */
  private int[] fallback;
  /** The score of {@link #fallback}, or below every mapping. */
  private long fallbackScore;
  /** The pairs mapped. */
  private long nodes;

  /** A search of {@code graph} for fragments of at least {@code minAtoms} atoms, one or two. */
  CommonAtomSearch(CommonBondGraphs graph, CommonBondGraphs.Clock clock, int minAtoms) {
    this.graph = graph;
    this.clock = clock;
    this.minAtoms = Math.max(1, minAtoms);
    query = graph.query;
    target = graph.target;
    queryCount = query.atomCount;
    targetCount = target.atomCount;
    words = CommonBondGraphs.words(targetCount);
    compatible = graph.compatiblePairs();
    bondClassesMatch = bondClassesMatch();
    queryToTarget = new int[queryCount];
    targetToQuery = new int[targetCount];
    Arrays.fill(queryToTarget, -1);
    Arrays.fill(targetToQuery, -1);
    levels = new Level[Math.min(queryCount, targetCount) + 1];
    triedSerial = new long[queryCount];
    queryAtoms = new int[queryCount];
    targetAtoms = new int[targetCount];
    queryPlace = new int[queryCount];
    targetPlace = new int[targetCount];
    queryLinkClass = new int[queryCount];
    targetLinkClass = new int[targetCount];
    queryLinkStamp = new long[queryCount];
    targetLinkStamp = new long[targetCount];
    linkedQuery = new int[query.maxDegree()];
    linkedQueryBlock = new int[linkedQuery.length];
    linkedTarget = new int[target.maxDegree()];
    linkedTargetBlock = new int[linkedTarget.length];
    queryLive = new long[queryCount];
    targetLive = new long[targetCount];
    queryBlock = new int[queryCount];
    targetBlock = new int[targetCount];
    queryKindCount = new int[graph.atomGroups.queryGroup.length];
    targetKindCount = new int[graph.atomGroups.targetGroup.length];
    twinSeen = new long[targetCount];
    targetKeys = new long[targetCount];
    leadingTaken = new int[targetCount];
    leadingRound = new int[targetCount];
    likestLevel = new int[queryCount];
    Arrays.fill(likestLevel, -1);
  }

  /**
   * Searches {@code query} in {@code target} for the single best MCIS in fragments of at least
   * {@code minFragmentAtoms} atoms, one or two, under {@code chemistry}, with the atoms of the
   * target in {@code excludedTargetAtoms} left out, and converts the mapping as the vertex search
   * of {@link MCS} reports it. Null when the graphs cannot be searched here, see
   * {@link CommonBondGraphs#supports}.
   */
  static SearchResult search(MolGraph query, MolGraph target, ChemOptions chemistry,
      Set<Integer> excludedTargetAtoms, int minFragmentAtoms, SearchResult.Control control) {
    if (!CommonBondGraphs.supports(query, target))
      return null;
    CommonBondGraphs.Clock clock = new CommonBondGraphs.Clock(control::stopped);
    CommonBondGraphs graph = null;
    CommonAtomSearch search = null;
    int[] found = null;
    // Each query atom is either mapped or left unmapped, so the molecule with fewer atoms is
    // searched as the query, unless target atoms are left out. The predicates are symmetric.
    boolean reversed = excludedTargetAtoms.isEmpty() && query.n > target.n;
    try {
      clock.poll();
      graph = reversed ? new CommonBondGraphs(target, query, chemistry, excludedTargetAtoms, clock)
                       : new CommonBondGraphs(query, target, chemistry, excludedTargetAtoms, clock);
      search = new CommonAtomSearch(graph, clock, minFragmentAtoms);
      found = search.run();
    } catch (CommonBondGraphs.Stop stop) {
      // The time ran out or the search was cancelled. What the search found by then is the result.
      if (search != null)
        found = search.found();
    }
    List<SearchResult.Match> matches = new ArrayList<>();
    if (found != null)
      matches.add(match(graph, found, reversed));
    boolean proved = control.status() == SearchResult.Status.Completed;
    double score = matches.isEmpty() ? 0.0 : matches.get(0).score();
    return new SearchResult(matches, control.status(), proved, proved, score,
        search == null ? 0 : search.nodes, control.elapsedMillis());
  }

  /**
   * The match of {@code found}, a mapping from the query of {@code graph} to its target, or from
   * the target to the query when {@code reversed}: its atoms, and the bonds between them in the
   * order of the query bonds, each with its image, and the fragments these bonds form. The score
   * is the number of atoms.
   */
  static SearchResult.Match match(CommonBondGraphs graph, int[] found, boolean reversed) {
    CommonBondGraphs.Side query = reversed ? graph.target : graph.query;
    int[] mapping = found;
    if (reversed) {
      mapping = new int[query.atomCount];
      Arrays.fill(mapping, -1);
      for (int atom = 0; atom < found.length; atom++)
        if (found[atom] >= 0)
          mapping[found[atom]] = atom;
    }
    Map<Integer, Integer> atoms = new TreeMap<>();
    for (int queryAtom = 0; queryAtom < query.atomCount; queryAtom++)
      if (mapping[queryAtom] >= 0)
        atoms.put(queryAtom, mapping[queryAtom]);
    Map<SearchResult.Bond, SearchResult.Bond> bonds = new LinkedHashMap<>();
    for (int queryBond = 0; queryBond < query.bondCount; queryBond++) {
      int first = mapping[query.bondFirst[queryBond]], second = mapping[query.bondSecond[queryBond]];
      if (first >= 0 && second >= 0)
        bonds.put(new SearchResult.Bond(query.bondFirst[queryBond], query.bondSecond[queryBond]),
            new SearchResult.Bond(first, second));
    }
    return new SearchResult.Match(atoms, bonds, MatchGraph.fragments(atoms, bonds.keySet()), atoms.size());
  }

  /**
   * A maximum mapping, or null when there is none, as no pair is compatible or, with fragments of
   * two atoms, no bond can be common. The atoms are ranked, a greedy seed gives the score to start
   * from, and the branch and bound finds the mapping.
   */
  int[] run() {
    clock.poll();
    if (!anyCompatible())
      return null;
    rank();
    long bound = bestPossible();
    fallbackScore = belowEveryMapping();
    seed(bound);
    return collectByAtoms(bound);
  }

  /**
   * The best mapping found so far, for a search the clock stopped: the best of the pass under way
   * or, while it holds none, the one it falls back on. Null when nothing was found. Each is
   * complete, as a mapping is kept only once it is recorded whole.
   */
  int[] found() {
    return best != null ? best : fallback;
  }

  /**
   * A maximum mapping below {@code bound}. The best score is sought first level by level of
   * common atoms, down from the bound to the score of the best mapping known: each level asks for
   * one mapping with at least that many atoms, and without one the bound comes down. A level far
   * above the best score prunes much and costs little, so this soon finds an optimum close to the
   * bound. Each level down repeats the work of the one above it, though, which pays only while
   * each level costs at least {@link #LEVEL_GROWTH} times the one above it: the levels then cost
   * little more than the last. Once a level costs less, the next would cost about as much, and one
   * search from the best mapping known finishes the work, its threshold rising as better mappings
   * are found.
   */
  private int[] collectByAtoms(long bound) {
    for (int triedAtoms = atoms(bound), previous = 0; triedAtoms > atoms(fallbackScore);
         triedAtoms--) {
      long before = nodes;
      int[] mapping = collect(null, score(triedAtoms, 0) - 1, bound);
      if (mapping != null)
        return mapping;
      // No mapping has that many atoms.
      bound = score(triedAtoms - 1, bonds(bound));
      long cost = nodes - before;
      if (previous > 0 && cost < (long) LEVEL_GROWTH * previous)
        break;
      previous = (int) Math.min(Integer.MAX_VALUE, Math.max(1, cost));
    }
    return collect(fallback, fallbackScore, bound);
  }

  /**
   * A score that no mapping beats. Every cap on the atoms or on the bonds holds for every mapping,
   * so the two are capped apart; every bond between mapped atoms is a common bond.
   */
  private long bestPossible() {
    long label = graph.labelBound(false);
    int atomCap = Math.min(Math.min(queryCount, targetCount), CommonBondGraphs.atoms(label));
    int bondCap = Math.min(Math.min(query.bondCount, target.bondCount), CommonBondGraphs.bonds(label));
    // A mapping with at most so many atoms has at most so many bonds between them.
    long pairs = (long) atomCap * (atomCap - 1) / 2;
    return score(atomCap, (int) Math.min(bondCap, pairs));
  }

  /** A score below that of every mapping, as every mapping has an atom. */
  private static long belowEveryMapping() {
    return score(0, 0);
  }

  /** Whether any pair of atoms is compatible. */
  private boolean anyCompatible() {
    for (long word : compatible)
      if (word != 0)
        return true;
    return false;
  }

  /**
   * Whether every query bond matches every target bond of its class. The bond predicate gives
   * bonds of one kind the same pairs, so a bond of each kind stands for the others.
   */
  private boolean bondClassesMatch() {
    int queryKinds = 0, targetKinds = 0;
    for (int kind : query.bondKind) queryKinds = Math.max(queryKinds, kind + 1);
    for (int kind : target.bondKind) targetKinds = Math.max(targetKinds, kind + 1);
    int[] queryExample = new int[queryKinds], targetExample = new int[targetKinds];
    Arrays.fill(queryExample, -1);
    Arrays.fill(targetExample, -1);
    for (int bond = query.bondCount - 1; bond >= 0; bond--) queryExample[query.bondKind[bond]] = bond;
    for (int bond = target.bondCount - 1; bond >= 0; bond--)
      targetExample[target.bondKind[bond]] = bond;
    for (int queryBond : queryExample) {
      clock.tick(1 + targetKinds);
      for (int targetBond : targetExample)
        if (queryBond >= 0 && targetBond >= 0 && query.bondClass[queryBond] >= 0
            && query.bondClass[queryBond] == target.bondClass[targetBond]
            && !graph.bondsMatch(queryBond, targetBond))
          return false;
    }
    return true;
  }

  /**
   * Ranks the atoms and works out their signatures. Query atoms with the fewest candidates, then
   * the most bonds, rank highest; ties are broken by colour refinement, so that the order is the
   * same for the same graph. Target atoms are ranked by their cell only, in input order within it.
   * The ranks order the seeds and break ties between candidates.
   */
  private void rank() {
    querySignature = query.signatures(clock);
    targetSignature = target.signatures(clock);
    int[] queryChoices = new int[queryCount], targetChoices = new int[targetCount];
    for (int queryAtom = 0; queryAtom < queryCount; queryAtom++) {
      for (int targetAtom = CommonBondGraphs.nextBit(compatible, queryAtom * words, words, 0); targetAtom >= 0;
           targetAtom = CommonBondGraphs.nextBit(compatible, queryAtom * words, words, targetAtom + 1)) {
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
    int runs = 0;
    int[] starts = new int[targetCount + 1];
    for (int position = 0; position < targetCount; position++)
      if (position == 0 || targetCell[targetOrder[position]] != targetCell[targetOrder[position - 1]])
        starts[runs++] = position;
    starts[runs++] = targetCount;
    targetRunStart = Arrays.copyOf(starts, runs);
    Integer[] byRank = new Integer[queryCount];
    for (int queryAtom = 0; queryAtom < queryCount; queryAtom++) byRank[queryAtom] = queryAtom;
    // The most choices come first, then the fewest bonds, then the lowest canonical rank.
    Comparator<Integer> byChoices = Comparator.comparingInt(queryAtom -> -queryChoices[queryAtom]);
    Arrays.sort(byRank,
        byChoices.thenComparingInt(query::degree).thenComparingInt(queryAtom -> canon[queryAtom]));
    queryOrder = new int[queryCount];
    queryRank = new int[queryCount];
    for (int rank = 0; rank < queryCount; rank++) {
      queryOrder[rank] = byRank[rank];
      queryRank[queryOrder[rank]] = rank;
    }
  }

  /**
   * The best greedy seed, if any, as the {@link #fallback}: a set of fragments grown one after
   * another, each on the atoms bonded to no mapped atom on either side. A fragment is grown from a
   * root pair, from each of at most {@link #SEEDS} of them or until one reaches {@code bound}, and
   * the best is kept. The fallback follows each better seed, so a search stopped while seeding
   * still has the best seed grown by then, or the seed being grown if that is better. It needs
   * {@link #rank} first.
   */
  private void seed(long bound) {
    int[] keptAtoms = new int[1];
    try {
      growSeeds(bound, keptAtoms);
    } catch (CommonBondGraphs.Stop stop) {
      // The fragments mapped when the clock stopped are induced, as each pair added kept them so.
      int mappedAtoms = 0;
      for (int image : queryToTarget)
        if (image >= 0)
          mappedAtoms++;
      if (mappedAtoms - keptAtoms[0] >= minAtoms && score(queryToTarget) > fallbackScore) {
        fallback = queryToTarget.clone();
        fallbackScore = score(fallback);
      }
      throw stop;
    }
    Arrays.fill(queryToTarget, -1);
    Arrays.fill(targetToQuery, -1);
  }

  /**
   * The rounds of {@link #seed}, with the atoms of the fragments kept from earlier rounds in
   * {@code keptAtoms[0]}.
   */
  private void growSeeds(long bound, int[] keptAtoms) {
    long[] taken = new long[queryCount * words];
    int[] roots = new int[LIKENESS_LEVELS * SEEDS], rootCount = new int[LIKENESS_LEVELS];
    targetRunFree = Arrays.copyOf(targetRunStart, targetRunStart.length - 1);
    queryCellClosed = new long[queryCount];
    int[] frontier = new int[queryCount], merged = new int[queryCount], mapped = new int[queryCount];
    int[] added = new int[queryCount];
    long[] keys = new long[queryCount];
    long bestSeed = belowEveryMapping();
    for (int round = 0, keptBonds = 0;;) {
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
          if (seedAtoms >= minAtoms
              && score(keptAtoms[0] + seedAtoms, keptBonds + seedBonds) > bestSeed) {
            bestSeed = score(keptAtoms[0] + seedAtoms, keptBonds + seedBonds);
            fallback = queryToTarget.clone();
            fallbackScore = score(fallback);
          }
          for (int mappedIndex = 0; mappedIndex < seedAtoms; mappedIndex++) {
            targetToQuery[queryToTarget[mapped[mappedIndex]]] = -1;
            queryToTarget[mapped[mappedIndex]] = -1;
          }
        }
      }
      if (bestSeed == before || bestSeed == bound)
        break;
      // The fragments found so far stay mapped while the next one grows from atoms bonded to none
      // of them. Each pair added keeps the mapping induced, so the fragments stay valid.
      keptAtoms[0] = atoms(bestSeed);
      keptBonds = bonds(bestSeed);
      System.arraycopy(fallback, 0, queryToTarget, 0, queryCount);
      for (int queryAtom = 0; queryAtom < queryCount; queryAtom++)
        if (queryToTarget[queryAtom] >= 0)
          targetToQuery[queryToTarget[queryAtom]] = queryAtom;
    }
  }

  /**
   * Picks the seed roots by how alike their atoms are, among the free atoms bonded to no mapped
   * atom. The {@code rootCount[level]} roots of each level are at {@code roots[level * SEEDS]}.
   * Alike cells grow alike seeds, so there is one root for each pair of cells, marked in
   * {@code taken}. Each query atom takes, from each run of target cells in target order, the
   * first free compatible target bonded to no mapped atom. Atoms of one cell look alike, so a run
   * whose level is full is passed whole, and a query cell with no run left is passed for the rest
   * of the call. The kept fragments only grow between calls, so a target that cannot be a root
   * stays so and the start of each run moves past it once.
   */
  private void seedRoots(int[] roots, int[] rootCount, long[] taken) {
    // Only the cell pairs of the last call's roots are taken.
    for (int level = 0; level < LIKENESS_LEVELS; level++)
      for (int index = 0; index < rootCount[level]; index++) {
        int packed = roots[level * SEEDS + index];
        int cell = queryCell[queryAtom(packed)], bit = targetCell[targetAtom(packed)];
        taken[cell * words + (bit >>> 6)] &= ~(1L << bit);
      }
    Arrays.fill(rootCount, 0);
    rootCall++;
    int likest = CommonBondGraphs.RADII, runs = targetRunStart.length - 1;
    // The most alike roots are tried first. Once there are enough of them, the rest are not needed.
    for (int rank = queryCount - 1; rank >= 0 && rootCount[likest] < SEEDS; rank--) {
      int queryAtom = queryOrder[rank], cell = queryCell[queryAtom];
      clock.tick();
      if (queryCellClosed[cell] == rootCall || queryToTarget[queryAtom] >= 0
          || mappedNeighbours(query, queryToTarget, queryAtom) > 0)
        continue;
      boolean closed = true;
      for (int run = 0; run < runs && rootCount[likest] < SEEDS; run++) {
        int first = targetRunFree[run], end = targetRunStart[run + 1];
        clock.tick();
        while (first < end && !rootTarget(targetOrder[first])) first++;
        targetRunFree[run] = first;
        if (first == end || has(taken, cell, targetCell[targetOrder[first]]))
          continue;
        int level = likeness(queryAtom, targetOrder[first]);
        if (rootCount[level] >= SEEDS)
          continue;
        boolean rooted = false;
        for (int position = first; position < end && !rooted; position++) {
          int targetAtom = targetOrder[position];
          clock.tick();
          if (!has(compatible, queryAtom, targetAtom) || !rootTarget(targetAtom))
            continue;
          set(taken, cell, targetCell[targetAtom]);
          roots[level * SEEDS + rootCount[level]++] = pair(queryAtom, targetAtom);
          rooted = true;
        }
        closed = closed && rooted;
      }
      if (closed)
        queryCellClosed[cell] = rootCall;
    }
  }

  /** Whether a target atom is free and bonded to no mapped atom, so that it may root a seed. */
  private boolean rootTarget(int targetAtom) {
    return targetToQuery[targetAtom] < 0 && mappedNeighbours(target, targetToQuery, targetAtom) == 0;
  }

  /** The number of mapped atoms bonded to {@code atom} of {@code side}, with {@code image} its mapping. */
  private int mappedNeighbours(CommonBondGraphs.Side side, int[] image, int atom) {
    int count = 0;
    clock.tick(1 + side.degree(atom));
    for (int place = side.start[atom]; place < side.start[atom + 1]; place++)
      if (image[side.neighbour[place]] >= 0)
        count++;
    return count;
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
    for (long pair = root; pair >= 0;) {
      int queryAtom = queryAtom((int) pair), targetAtom = targetAtom((int) pair);
      queryToTarget[queryAtom] = targetAtom;
      targetToQuery[targetAtom] = queryAtom;
      mapped[count++] = queryAtom;
      seedBonds += (int) (pair >>> 32);
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
      long held = -1;
      int kept = 0;
      for (int frontierIndex = 0; frontierIndex < frontierSize; frontierIndex++) {
        int frontierAtom = frontier[frontierIndex];
        // The atom just mapped leaves the frontier.
        if (queryToTarget[frontierAtom] >= 0)
          continue;
        frontier[kept++] = frontierAtom;
        held = bestSeedPair(frontierAtom, held, round);
      }
      frontierSize = kept;
      pair = held;
    }
    return score(count, seedBonds);
  }

  /**
   * The best pair of {@code frontierAtom} with a free target next to the image of one of its
   * mapped neighbours, unless the pair {@code held} beats it. The pair must keep the mapping
   * induced, see {@link #inducedLinks}. The most alike pair wins, then the one with the most new
   * bonds; ties go to the higher ranked query atom, then to the lower ranked target atom. The pair
   * is returned in the low int, with its new bonds in the high int, or -1 if there is none;
   * {@code held} is read the same way. Next to a hub, the first few free targets are as good for a
   * seed as all of them. The seed of {@code round} only takes targets, so the taken targets that
   * lead the list of an image are passed once for the seed, see {@link #leadingTaken}, and not
   * again for each frontier atom.
   */
  private long bestSeedPair(int frontierAtom, long held, int round) {
    int pair = held < 0 ? -1 : (int) held, links = held < 0 ? 0 : (int) (held >>> 32);
    int likest = pair < 0 ? -1 : likeness(queryAtom(pair), targetAtom(pair));
    clock.tick(1 + query.degree(frontierAtom));
    if (pair >= 0) {
      // No pair of the atom is liker than its likest level, nor has more new bonds than the atom
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
        if (targetToQuery[candidateTarget] >= 0 || !has(compatible, frontierAtom, candidateTarget))
          continue;
        freeCount++;
        int level = likeness(frontierAtom, candidateTarget);
        if (level < likest)
          continue;
        int newBonds = inducedLinks(frontierAtom, candidateTarget);
        if (newBonds < 0)
          continue;
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

  /**
   * The bonds that mapping {@code queryAtom} onto {@code targetAtom} adds, or -1 if the mapping
   * would not stay induced: the target atom must be bonded to the image of each mapped neighbour of
   * the query atom by a compatible bond, and to no other mapped atom. As the images of distinct
   * neighbours are distinct, the second holds when both atoms have as many mapped neighbours.
   */
  private int inducedLinks(int queryAtom, int targetAtom) {
    int queryLinks = 0;
    clock.tick(1 + query.degree(queryAtom) * (1 + target.degree(targetAtom)));
    for (int queryPlace = query.start[queryAtom]; queryPlace < query.start[queryAtom + 1];
         queryPlace++) {
      int neighbourImage = queryToTarget[query.neighbour[queryPlace]];
      if (neighbourImage < 0)
        continue;
      queryLinks++;
      int targetBond = target.bondBetween(targetAtom, neighbourImage);
      if (targetBond < 0 || !graph.bondsMatch(query.neighbourBond[queryPlace], targetBond))
        return -1;
    }
    return queryLinks == mappedNeighbours(target, targetToQuery, targetAtom) ? queryLinks : -1;
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
   * The best mapping that beats the score {@code start}, which the {@code seed} holds unless it is
   * null, or the seed when none does. The search ends early once it holds a mapping that reaches
   * {@code upperBound}. It needs {@link #rank} first.
   */
  private int[] collect(int[] seed, long start, long upperBound) {
    bestScore = start;
    this.upperBound = upperBound;
    best = seed;
    // The twins are found once, before any search.
    if (queryTwin == null) {
      queryTwin = twins(query, true);
      targetTwin = twins(target, false);
      targetTwinFirst = ringFirsts(targetTwin);
    }
    Level base = levels[0] != null ? levels[0] : (levels[0] = new Level());
    depth = 0;
    commonBonds = 0;
    base.queryAtom = base.targetAtom = -1;
    base.fragmentStart = 1;
    base.freeQueryBonds = query.bondCount;
    base.freeTargetBonds = target.bondCount;
    base.serial = ++serial;
    base.triedAtom = -1;
    base.forced = false;
    base.targetCount = base.nextTarget = 0;
    base.blocksFrom = blocks.count = 0;
    baseBlocks();
    base.blocksTo = blocks.count;
    base.bound = bound(base);
    search();
    return best;
  }

  /**
   * Sets up the blocks of the base level: for each atom class, its query and target atoms in
   * index order. The atoms that match nothing come after every block.
   */
  private void baseBlocks() {
    int classes = graph.atomClasses;
    int[] queryStart = new int[classes + 1], targetStart = new int[classes + 1];
    for (int queryAtom = 0; queryAtom < queryCount; queryAtom++)
      if (query.atomClass[queryAtom] >= 0)
        queryStart[query.atomClass[queryAtom] + 1]++;
    for (int targetAtom = 0; targetAtom < targetCount; targetAtom++)
      if (target.atomClass[targetAtom] >= 0)
        targetStart[target.atomClass[targetAtom] + 1]++;
    for (int atomClass = 0; atomClass < classes; atomClass++) {
      queryStart[atomClass + 1] += queryStart[atomClass];
      targetStart[atomClass + 1] += targetStart[atomClass];
    }
    int[] queryNext = Arrays.copyOf(queryStart, classes), targetNext = Arrays.copyOf(targetStart, classes);
    int queryOutside = queryStart[classes], targetOutside = targetStart[classes];
    clock.tick(queryCount + targetCount);
    for (int queryAtom = 0; queryAtom < queryCount; queryAtom++) {
      int atomClass = query.atomClass[queryAtom];
      int place = atomClass >= 0 ? queryNext[atomClass]++ : queryOutside++;
      queryAtoms[place] = queryAtom;
      queryPlace[queryAtom] = place;
    }
    for (int targetAtom = 0; targetAtom < targetCount; targetAtom++) {
      int atomClass = target.atomClass[targetAtom];
      int place = atomClass >= 0 ? targetNext[atomClass]++ : targetOutside++;
      targetAtoms[place] = targetAtom;
      targetPlace[targetAtom] = place;
    }
    for (int atomClass = 0; atomClass < classes; atomClass++) {
      int queryAtomsOfClass = queryStart[atomClass + 1] - queryStart[atomClass];
      int targetAtomsOfClass = targetStart[atomClass + 1] - targetStart[atomClass];
      if (queryAtomsOfClass > 0 && targetAtomsOfClass > 0)
        blocks.add(queryStart[atomClass], queryAtomsOfClass, targetStart[atomClass],
            targetAtomsOfClass, 0, false, atomClass);
    }
  }

  /**
   * Searches depth first from the base level until it is taken off. A level tries the targets of
   * each of its candidates in turn, see {@link #nextCandidate}. It is taken off once nothing is
   * left or its bound is too low. The levels are an explicit stack, so large molecules do not run
   * out of stack.
   */
  private void search() {
    while (depth >= 0) {
      clock.tick();
      Level level = levels[depth];
      if (done() || !canImprove(level.bound)) {
        pop();
      } else if (level.nextTarget < level.targetCount) {
        push(level.triedAtom, level.targets[level.nextTarget++], level.triedRoot);
      } else if (!nextCandidate(level)) {
        pop();
      }
    }
  }

  /**
   * Maps {@code queryAtom} onto {@code targetAtom} as a new level, which roots a new fragment when
   * {@code root} is set and otherwise joins the open one. The mapping is recorded, and the level
   * gets its blocks and its bound unless it is a leaf.
   */
  private void push(int queryAtom, int targetAtom, boolean root) {
    nodes++;
    Level parent = levels[depth];
    if (levels[depth + 1] == null)
      levels[depth + 1] = new Level();
    Level level = levels[depth + 1];
    level.queryAtom = queryAtom;
    level.targetAtom = targetAtom;
    level.fragmentStart = root ? depth + 1 : parent.fragmentStart;
    level.commonBondsBefore = commonBonds;
    level.freeQueryBonds = parent.freeQueryBonds;
    level.freeTargetBonds = parent.freeTargetBonds;
    clock.tick(1 + query.degree(queryAtom) + target.degree(targetAtom));
    for (int queryPlace = query.start[queryAtom]; queryPlace < query.start[queryAtom + 1];
         queryPlace++) {
      if (queryToTarget[query.neighbour[queryPlace]] < 0)
        level.freeQueryBonds--;
      else
        commonBonds++;
    }
    for (int targetPlace = target.start[targetAtom]; targetPlace < target.start[targetAtom + 1];
         targetPlace++)
      if (targetToQuery[target.neighbour[targetPlace]] < 0)
        level.freeTargetBonds--;
    queryToTarget[queryAtom] = targetAtom;
    targetToQuery[targetAtom] = queryAtom;
    depth++;
    level.serial = ++serial;
    level.triedAtom = -1;
    level.forced = false;
    level.targetCount = level.nextTarget = 0;
    level.blocksFrom = blocks.count;
    record();
    if (isLeaf()) {
      level.blocksTo = blocks.count;
      level.bound = -1;
      return;
    }
    refine(parent, queryAtom, targetAtom, root);
    level.blocksTo = blocks.count;
    level.bound = bound(level);
  }

  /** Takes the top level off, putting back the blocks and the mapping of its parent; the base level ends the search. */
  private void pop() {
    Level level = levels[depth];
    blocks.count = level.blocksFrom;
    if (depth > 0) {
      queryToTarget[level.queryAtom] = -1;
      targetToQuery[level.targetAtom] = -1;
      commonBonds = level.commonBondsBefore;
    }
    depth--;
  }

  /**
   * Whether the mapping is a leaf, with nothing left to add: the pass can end, or the mapping has
   * the atoms of the upper bound or every atom of either molecule.
   */
  private boolean isLeaf() {
    return done() || depth == atoms(upperBound) || depth == queryCount || depth == targetCount;
  }

  /**
   * Gives {@code level} its next candidate and the targets to try it on, or false when nothing is
   * left. The candidate tried so far is left unmapped below the level from now on, as every
   * mapping with it was tried, which can lower the bound. A forced pair, see {@link #forcedPair},
   * is the only one the level then tries. While the open fragment has query atoms bonded to it,
   * the next candidate is one of them, see {@link #candidate}. Once there are none the fragment is
   * closed and, unless it is too small, the next candidate is a query atom bonded to no mapped
   * atom, the root of a new fragment. A candidate with a twin tried at this level is skipped, see
   * {@link #twinTried}.
   */
  private boolean nextCandidate(Level level) {
    for (;;) {
      clock.tick();
      if (level.triedAtom >= 0) {
        // A forced pair holds a best extension of the state, so nothing else is needed.
        if (level.forced)
          return false;
        leaveUnmapped(level.triedBlock, level.triedAtom);
        triedSerial[level.triedAtom] = level.serial;
        level.triedAtom = -1;
        level.bound = bound(level);
        if (!canImprove(level.bound))
          return false;
      }
      if (forcedPair(level)) {
        level.triedAtom = forcedQuery;
        level.triedBlock = forcedBlock;
        level.triedRoot = false;
        level.forced = true;
        level.nextTarget = 0;
        if (level.targets.length == 0)
          level.targets = new int[1];
        level.targets[0] = forcedTarget;
        level.targetCount = 1;
        return true;
      }
      long chosen = depth > 0 ? candidate(level, true) : -1;
      boolean root = chosen < 0;
      if (root) {
        if (depth > 0 && depth - level.fragmentStart + 1 < minAtoms)
          return false;
        chosen = candidate(level, false);
        if (chosen < 0)
          return false;
      }
      int block = (int) (chosen >>> 32), atom = (int) chosen;
      level.triedAtom = atom;
      level.triedBlock = block;
      level.triedRoot = root;
      level.nextTarget = 0;
      if (twinTried(atom, level.serial)) {
        level.targetCount = 0;
        continue;
      }
      targets(level, block, atom);
      if (level.targetCount > 0)
        return true;
    }
  }

  /**
   * Whether the state of {@code level} has a forced pair, which some best extension of the state
   * maps, so that no other pair needs to be tried: a query and a target atom of one block, neither
   * bonded to an atom of a live block, so neither to an atom that may still be mapped, where
   * every pair of the block matches, see {@link #blockComplete}. Swapping such an atom with the one
   * an extension maps onto its partner, or adding the pair when neither is mapped, keeps the
   * extension induced without losing atoms or bonds. A pair bonded to no mapped atom forms a
   * fragment of one atom, so it is forced only when such fragments are allowed. The pair is left
   * in {@link #forcedBlock}, {@link #forcedQuery} and {@link #forcedTarget}: of the first block
   * that has one, the query atom of highest rank and the first such target of the run.
   */
  private boolean forcedPair(Level level) {
    for (int block = level.blocksFrom; block < level.blocksTo; block++) {
      clock.tick();
      int queryAtomCount = blocks.queryCount[block], targetAtomCount = blocks.targetCount[block];
      if (queryAtomCount == 0 || targetAtomCount == 0 || blocks.links[block] == 0 && minAtoms > 1
          || !blockComplete(block))
        continue;
      int chosen = -1;
      for (int place = blocks.queryStart[block]; place < blocks.queryStart[block] + queryAtomCount;
           place++) {
        int atom = queryAtoms[place];
        clock.tick();
        if ((chosen < 0 || queryRank[atom] > queryRank[chosen])
            && !hasLiveNeighbour(level, query, queryPlace, atom, true)
            && !twinTried(atom, level.serial))
          chosen = atom;
      }
      if (chosen < 0)
        continue;
      for (int place = blocks.targetStart[block];
           place < blocks.targetStart[block] + targetAtomCount; place++) {
        int atom = targetAtoms[place];
        clock.tick();
        if (!hasLiveNeighbour(level, target, targetPlace, atom, false)) {
          forcedBlock = block;
          forcedQuery = chosen;
          forcedTarget = atom;
          return true;
        }
      }
    }
    return false;
  }

  /**
   * Whether {@code atom} of {@code side}, with the atoms at {@code place}, is bonded to an atom of
   * a live block of {@code level}, one with atoms on both sides.
   */
  private boolean hasLiveNeighbour(Level level, CommonBondGraphs.Side side, int[] place, int atom,
      boolean onQuery) {
    for (int neighbourPlace = side.start[atom]; neighbourPlace < side.start[atom + 1];
         neighbourPlace++) {
      int block = blockOf(level, place[side.neighbour[neighbourPlace]], onQuery);
      if (block >= 0 && blocks.queryCount[block] > 0 && blocks.targetCount[block] > 0)
        return true;
    }
    return false;
  }

  /**
   * The next candidate from the open blocks of {@code level}, which joins the open fragment, or
   * from the others, which roots a new fragment, as its block in the high int and the atom in the
   * low int, or -1 if there is none. The block is the one with the fewest atoms on its larger side,
   * then the earliest; the atom is the one of its query atoms with the most bonds, then of highest
   * rank.
   */
  private long candidate(Level level, boolean open) {
    int chosenBlock = -1, chosenSize = Integer.MAX_VALUE;
    for (int block = level.blocksFrom; block < level.blocksTo; block++) {
      int queryAtomCount = blocks.queryCount[block], targetAtomCount = blocks.targetCount[block];
      if (blocks.open[block] != open || queryAtomCount == 0 || targetAtomCount == 0)
        continue;
      int size = Math.max(queryAtomCount, targetAtomCount);
      if (size < chosenSize) {
        chosenBlock = block;
        chosenSize = size;
      }
    }
    clock.tick(level.blocksTo - level.blocksFrom);
    if (chosenBlock < 0)
      return -1;
    int start = blocks.queryStart[chosenBlock], chosen = -1;
    clock.tick(blocks.queryCount[chosenBlock]);
    for (int place = start; place < start + blocks.queryCount[chosenBlock]; place++) {
      int atom = queryAtoms[place];
      if (chosen < 0 || query.degree(atom) > query.degree(chosen)
          || query.degree(atom) == query.degree(chosen) && queryRank[atom] > queryRank[chosen])
        chosen = atom;
    }
    return (long) chosenBlock << 32 | chosen;
  }

  /**
   * Whether a twin of {@code queryAtom} was tried at the level of serial {@code levelSerial}.
   * Swapping the two takes every mapping below the atom onto one of the same score below the
   * twin, which was searched before. The twin is then unmapped, and the atom is left so too.
   */
  private boolean twinTried(int queryAtom, long levelSerial) {
    for (int twin = queryTwin[queryAtom]; twin >= 0 && twin != queryAtom; twin = queryTwin[twin]) {
      clock.tick();
      if (triedSerial[twin] == levelSerial)
        return true;
    }
    return false;
  }

  /**
   * Lists in {@code level} the targets of {@code queryAtom}, which is in {@code block}: the target
   * atoms of the block it is compatible with, by compatible bonds to the images of its mapped
   * neighbours. The most alike come first, then by descending rank. Of target twins, only the first
   * is kept, as swapping two free twins takes every mapping onto one of the same score.
   */
  private void targets(Level level, int block, int queryAtom) {
    int start = blocks.targetStart[block], count = blocks.targetCount[block];
    if (level.targets.length < count)
      level.targets = new int[Math.max(count, 2 * level.targets.length)];
    int listed = 0;
    clock.tick(1 + count);
    for (int place = start; place < start + count; place++) {
      int targetAtom = targetAtoms[place];
      if (has(compatible, queryAtom, targetAtom)
          && (bondClassesMatch || bondsFit(queryAtom, targetAtom)))
        targetKeys[listed++] = targetKey(queryAtom, targetAtom);
    }
    Arrays.sort(targetKeys, 0, listed);
    long stamp = ++twinStamp;
    int kept = 0;
    for (int i = 0; i < listed; i++) {
      int targetAtom = targetOf(targetKeys[i]), first = targetTwinFirst[targetAtom];
      if (first >= 0) {
        if (twinSeen[first] == stamp)
          continue;
        twinSeen[first] = stamp;
      }
      level.targets[kept++] = targetAtom;
    }
    level.targetCount = kept;
  }

  /** Whether the bonds from {@code queryAtom} to its mapped neighbours match those from {@code targetAtom} to their images. */
  private boolean bondsFit(int queryAtom, int targetAtom) {
    clock.tick(1 + query.degree(queryAtom) * (1 + target.degree(targetAtom)));
    for (int queryPlace = query.start[queryAtom]; queryPlace < query.start[queryAtom + 1];
         queryPlace++) {
      int neighbourImage = queryToTarget[query.neighbour[queryPlace]];
      if (neighbourImage < 0)
        continue;
      int targetBond = target.bondBetween(targetAtom, neighbourImage);
      if (targetBond < 0 || !graph.bondsMatch(query.neighbourBond[queryPlace], targetBond))
        return false;
    }
    return true;
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

  /**
   * Leaves {@code queryAtom} of {@code block} unmapped below the level that holds the block: it
   * moves to the end of the live query atoms of the block, which no longer count it.
   */
  private void leaveUnmapped(int block, int queryAtom) {
    blocks.exactPairs[block] = blocks.complete[block] = -1;
    int last = blocks.queryStart[block] + --blocks.queryCount[block];
    swapQuery(queryPlace[queryAtom], last);
  }

  /**
   * Gives the new level the blocks of its parent split by the bonds to the new pair: the query
   * atoms of a block bonded to {@code queryAtom} and the target atoms bonded to {@code targetAtom}
   * by bonds of one class form a new block, open as it is bonded to the open fragment, and the
   * atoms bonded to neither stay together. The pair leaves its blocks, and a block with nothing
   * left on one side is dropped, as its atoms can no longer be mapped. A new fragment closes the
   * one before, so the open blocks of the parent, whose query atoms were all tried, are dropped
   * and no other block stays open. Only the blocks of the pair and of its neighbours change, so
   * the others are copied as they are; the blocks of a level are in order of their runs on both
   * sides, which finds the block of an atom by bisection.
   */
  private void refine(Level parent, int queryAtom, int targetAtom, boolean root) {
    long stamp = ++linkStamp;
    int queryLinked = gatherLinks(parent, query, queryAtom, true, stamp);
    int targetLinked = gatherLinks(parent, target, targetAtom, false, stamp);
    int queryPairBlock = blockOf(parent, queryPlace[queryAtom], true);
    int targetPairBlock = blockOf(parent, targetPlace[targetAtom], false);
    clock.tick(1 + parent.blocksTo - parent.blocksFrom);
    for (int block = parent.blocksFrom, queryNext = 0, targetNext = 0; block < parent.blocksTo;
         block++) {
      int queryStart = blocks.queryStart[block], queryAtomCount = blocks.queryCount[block];
      int targetStart = blocks.targetStart[block], targetAtomCount = blocks.targetCount[block];
      int links = blocks.links[block], atomClass = blocks.atomClass[block];
      boolean open = blocks.open[block];
      int queryFirstLinked = queryNext, targetFirstLinked = targetNext;
      while (queryNext < queryLinked && linkedQueryBlock[queryNext] == block) queryNext++;
      while (targetNext < targetLinked && linkedTargetBlock[targetNext] == block) targetNext++;
      if (queryAtomCount == 0 || targetAtomCount == 0 || root && open)
        continue;
      if (block != queryPairBlock && block != targetPairBlock && queryNext == queryFirstLinked
          && targetNext == targetFirstLinked) {
        blocks.add(queryStart, queryAtomCount, targetStart, targetAtomCount, links, open && !root,
            atomClass);
        continue;
      }
      if (block == queryPairBlock)
        swapQuery(queryPlace[queryAtom], queryStart + --queryAtomCount);
      if (block == targetPairBlock)
        swapTarget(targetPlace[targetAtom], targetStart + --targetAtomCount);
      if (queryAtomCount == 0 || targetAtomCount == 0)
        continue;
      // The linked atoms move to the end of the live runs, in order of bond class.
      int queryEnd = queryStart + queryAtomCount, targetEnd = targetStart + targetAtomCount;
      int queryUnlinked = queryEnd, targetUnlinked = targetEnd;
      for (int i = queryFirstLinked; i < queryNext; i++)
        swapQuery(queryPlace[linkedQuery[i]], --queryUnlinked);
      for (int i = targetFirstLinked; i < targetNext; i++)
        swapTarget(targetPlace[linkedTarget[i]], --targetUnlinked);
      sortByLink(queryAtoms, queryPlace, queryUnlinked, queryEnd, queryLinkClass);
      sortByLink(targetAtoms, targetPlace, targetUnlinked, targetEnd, targetLinkClass);
      if (queryUnlinked > queryStart && targetUnlinked > targetStart)
        blocks.add(queryStart, queryUnlinked - queryStart, targetStart, targetUnlinked - targetStart,
            links, open && !root, atomClass);
      // The linked atoms of both sides are in order of bond class; equal classes form a block.
      int queryFrom = queryUnlinked, targetFrom = targetUnlinked;
      while (queryFrom < queryEnd && targetFrom < targetEnd) {
        int queryClass = queryLinkClass[queryAtoms[queryFrom]];
        int targetClass = targetLinkClass[targetAtoms[targetFrom]];
        int queryTo = queryFrom, targetTo = targetFrom;
        if (queryClass <= targetClass)
          while (queryTo < queryEnd && queryLinkClass[queryAtoms[queryTo]] == queryClass) queryTo++;
        if (targetClass <= queryClass)
          while (targetTo < targetEnd && targetLinkClass[targetAtoms[targetTo]] == targetClass)
            targetTo++;
        if (queryClass == targetClass && queryClass >= 0)
          blocks.add(queryFrom, queryTo - queryFrom, targetFrom, targetTo - targetFrom, links + 1,
              true, atomClass);
        queryFrom = queryTo;
        targetFrom = targetTo;
      }
    }
  }

  /**
   * Lists the live atoms of the blocks of {@code parent} that are bonded to {@code atom} of
   * {@code side}, in {@link #linkedQuery} or {@link #linkedTarget}, with their blocks, in order of
   * block, and marks each with {@code stamp} and the class of its bond to the atom. Gives their
   * number.
   */
  private int gatherLinks(Level parent, CommonBondGraphs.Side side, int atom, boolean onQuery,
      long stamp) {
    int[] linked = onQuery ? linkedQuery : linkedTarget;
    int[] linkedBlock = onQuery ? linkedQueryBlock : linkedTargetBlock;
    int[] place = onQuery ? queryPlace : targetPlace;
    long[] linkStampOf = onQuery ? queryLinkStamp : targetLinkStamp;
    int[] linkClass = onQuery ? queryLinkClass : targetLinkClass;
    int count = 0;
    clock.tick(1 + side.degree(atom));
    for (int neighbourPlace = side.start[atom]; neighbourPlace < side.start[atom + 1];
         neighbourPlace++) {
      int other = side.neighbour[neighbourPlace], block = blockOf(parent, place[other], onQuery);
      if (block < 0)
        continue;
      linkStampOf[other] = stamp;
      linkClass[other] = side.bondClass[side.neighbourBond[neighbourPlace]];
      // An insertion by block, which keeps the order of the neighbours within a block.
      int i = count++;
      while (i > 0 && linkedBlock[i - 1] > block) {
        linked[i] = linked[i - 1];
        linkedBlock[i] = linkedBlock[i - 1];
        i--;
      }
      linked[i] = other;
      linkedBlock[i] = block;
    }
    return count;
  }

  /**
   * The block of {@code parent} whose live run of query atoms, or else of target atoms, holds
   * place {@code place}, or -1 if none does. The runs are in order, so this is a bisection.
   */
  private int blockOf(Level parent, int place, boolean onQuery) {
    int[] start = onQuery ? blocks.queryStart : blocks.targetStart;
    int[] count = onQuery ? blocks.queryCount : blocks.targetCount;
    int low = parent.blocksFrom, high = parent.blocksTo - 1, found = -1;
    while (low <= high) {
      int middle = (low + high) >>> 1;
      if (start[middle] <= place) {
        found = middle;
        low = middle + 1;
      } else {
        high = middle - 1;
      }
    }
    clock.tick();
    return found >= 0 && place < start[found] + count[found] ? found : -1;
  }

  /**
   * Sorts the atoms from place {@code from} up to {@code to} of {@code atoms} by the class of their
   * bond to the new pair, keeping their order otherwise. There are few, so an insertion sort will do.
   */
  private void sortByLink(int[] atoms, int[] place, int from, int to, int[] linkClass) {
    clock.tick((to - from) * (to - from));
    for (int i = from + 1; i < to; i++) {
      int atom = atoms[i], key = linkClass[atom], j = i;
      while (j > from && linkClass[atoms[j - 1]] > key) {
        atoms[j] = atoms[j - 1];
        place[atoms[j]] = j;
        j--;
      }
      atoms[j] = atom;
      place[atom] = j;
    }
  }

  /** Swaps the query atoms at places {@code first} and {@code second}. */
  private void swapQuery(int first, int second) {
    int atom = queryAtoms[first];
    queryAtoms[first] = queryAtoms[second];
    queryAtoms[second] = atom;
    queryPlace[queryAtoms[first]] = first;
    queryPlace[atom] = second;
  }

  /** Swaps the target atoms at places {@code first} and {@code second}. */
  private void swapTarget(int first, int second) {
    int atom = targetAtoms[first];
    targetAtoms[first] = targetAtoms[second];
    targetAtoms[second] = atom;
    targetPlace[targetAtoms[first]] = first;
    targetPlace[atom] = second;
  }

  /**
   * Keeps the mapping if it beats the best score, and otherwise as the fallback if it beats that.
   * The open fragment must have at least {@link #minAtoms} atoms; the closed ones had them.
   */
  private void record() {
    if (depth - levels[depth].fragmentStart + 1 < minAtoms)
      return;
    long current = score(depth, commonBonds);
    if (current > bestScore) {
      bestScore = current;
      best = queryToTarget.clone();
    } else if (current > fallbackScore) {
      fallback = queryToTarget.clone();
      fallbackScore = current;
    }
  }

  /** Whether the pass can end, as it holds a mapping with a score that no mapping beats. */
  private boolean done() {
    return best != null && bestScore >= upperBound;
  }

  /**
   * The bound of the state of {@code level}. Each block pairs off at most as many atoms as it has
   * on its smaller side, see {@link #blockPairs}. Each atom it adds brings a bond to each mapped
   * atom its block is linked to, and the bonds between the added atoms are at most the bonds with
   * no end mapped on either side, or more closely those of {@link #liveBonds}, which cost more and
   * are counted only when the bonds may decide.
   */
  private long bound(Level level) {
    int newAtoms = 0;
    long newLinks = 0;
    for (int block = level.blocksFrom; block < level.blocksTo; block++) {
      if (blocks.queryCount[block] == 0 || blocks.targetCount[block] == 0)
        continue;
      int pairs = blockPairs(block);
      newAtoms += pairs;
      newLinks += (long) pairs * blocks.links[block];
    }
    clock.tick(1 + level.blocksTo - level.blocksFrom);
    int atoms = depth + newAtoms;
    long loose = capped(atoms, (int) Math.min(Integer.MAX_VALUE,
        commonBonds + newLinks + Math.min(level.freeQueryBonds, level.freeTargetBonds)));
    // Bonds decide only between mappings of as many atoms as the best one.
    if (!canImprove(loose) || atoms(loose) > atoms(bestScore))
      return loose;
    return capped(atoms, (int) Math.min(Integer.MAX_VALUE, commonBonds + newLinks + liveBonds(level)));
  }

  /**
   * The most bonds between atoms the state of {@code level} may still map. A bond between two
   * such atoms has its image between atoms of the same two blocks, by a bond of the same class, so
   * for each two live blocks and bond class the fewer of these bonds on either side count.
   */
  private int liveBonds(Level level) {
    long stamp = ++liveStamp;
    int live = 0;
    for (int block = level.blocksFrom; block < level.blocksTo; block++) {
      int queryAtomCount = blocks.queryCount[block], targetAtomCount = blocks.targetCount[block];
      if (queryAtomCount == 0 || targetAtomCount == 0)
        continue;
      clock.tick(queryAtomCount + targetAtomCount);
      for (int place = blocks.queryStart[block]; place < blocks.queryStart[block] + queryAtomCount;
           place++) {
        queryLive[queryAtoms[place]] = stamp;
        queryBlock[queryAtoms[place]] = live;
      }
      for (int place = blocks.targetStart[block];
           place < blocks.targetStart[block] + targetAtomCount; place++) {
        targetLive[targetAtoms[place]] = stamp;
        targetBlock[targetAtoms[place]] = live;
      }
      live++;
    }
    int classes = Math.max(1, graph.bondClasses);
    if (queryPairCount.length < live * classes) {
      queryPairCount = new int[2 * live * classes];
      targetPairCount = new int[2 * live * classes];
      touchedPairs = new int[2 * live * classes];
    }
    int total = 0;
    for (int block = level.blocksFrom, index = 0; block < level.blocksTo; block++) {
      int queryAtomCount = blocks.queryCount[block], targetAtomCount = blocks.targetCount[block];
      if (queryAtomCount == 0 || targetAtomCount == 0)
        continue;
      int touched = countBonds(query, queryAtoms, blocks.queryStart[block], queryAtomCount,
          queryLive, queryBlock, stamp, index, classes, queryPairCount, targetPairCount, 0);
      touched = countBonds(target, targetAtoms, blocks.targetStart[block], targetAtomCount,
          targetLive, targetBlock, stamp, index, classes, targetPairCount, queryPairCount, touched);
      for (int i = 0; i < touched; i++) {
        int key = touchedPairs[i];
        total += Math.min(queryPairCount[key], targetPairCount[key]);
        queryPairCount[key] = targetPairCount[key] = 0;
      }
      index++;
    }
    return total;
  }

  /**
   * Counts in {@code counts}, by the live block of the other end and the bond class, the bonds of
   * {@code side} from the {@code atomCount} atoms from {@code start} of {@code atoms}, all of live
   * block {@code index}, to live atoms of this or a later block, each bond once. A key counted for
   * the first time on both sides is added to {@link #touchedPairs} after the first
   * {@code touched}, and the new number of keys there is returned.
   */
  private int countBonds(CommonBondGraphs.Side side, int[] atoms, int start, int atomCount,
      long[] live, int[] blockOfAtom, long stamp, int index, int classes, int[] counts,
      int[] otherCounts, int touched) {
    int work = 0;
    for (int place = start; place < start + atomCount; place++) {
      int atom = atoms[place];
      work += 1 + side.degree(atom);
      for (int other = side.start[atom]; other < side.start[atom + 1]; other++) {
        int neighbour = side.neighbour[other], bondClass = side.bondClass[side.neighbourBond[other]];
        if (live[neighbour] != stamp || bondClass < 0 || blockOfAtom[neighbour] < index
            || blockOfAtom[neighbour] == index && neighbour < atom)
          continue;
        int key = blockOfAtom[neighbour] * classes + bondClass;
        if (counts[key]++ == 0 && otherCounts[key] == 0)
          touchedPairs[touched++] = key;
      }
    }
    clock.tick(work);
    return touched;
  }

  /**
   * The most atoms {@code block} pairs off. A block linked to mapped atoms, where bonds of one
   * class need not match, gives its exact pairs, see {@link #exactPairs}. Otherwise it is the fewer
   * of its live atoms on either side or, where its class holds kinds that do not all match, the
   * b-matching of their kinds, see {@link CommonBondGraphs.Groups}.
   */
  private int blockPairs(int block) {
    int queryStart = blocks.queryStart[block], queryAtomCount = blocks.queryCount[block];
    int targetStart = blocks.targetStart[block], targetAtomCount = blocks.targetCount[block];
    if (!bondClassesMatch && blocks.links[block] > 0 && targetAtomCount <= Long.SIZE
        && queryAtomCount <= Long.SIZE) {
      if (blocks.exactPairs[block] < 0)
        blocks.exactPairs[block] = exactPairs(block);
      return blocks.exactPairs[block];
    }
    if (graph.atomGroups.complete[blocks.atomClass[block]])
      return Math.min(queryAtomCount, targetAtomCount);
    clock.tick(queryAtomCount + targetAtomCount);
    for (int place = queryStart; place < queryStart + queryAtomCount; place++)
      queryKindCount[query.atomKind[queryAtoms[place]]]++;
    for (int place = targetStart; place < targetStart + targetAtomCount; place++)
      targetKindCount[target.atomKind[targetAtoms[place]]]++;
    int pairs = graph.atomGroups.bound(blocks.atomClass[block], queryKindCount, targetKindCount);
    for (int place = queryStart; place < queryStart + queryAtomCount; place++)
      queryKindCount[query.atomKind[queryAtoms[place]]] = 0;
    for (int place = targetStart; place < targetStart + targetAtomCount; place++)
      targetKindCount[target.atomKind[targetAtoms[place]]] = 0;
    return pairs;
  }

  /**
   * The most pairs of {@code block}, a block linked to mapped atoms with at most 64 atoms on each
   * side, that match exactly: the atoms and the bonds to the mapped atoms. Bonds of one class need
   * not match, so the pairs are a maximum bipartite matching. Such blocks are small, as their atoms
   * are bonded to a mapped atom.
   */
  private int exactPairs(int block) {
    int queryStart = blocks.queryStart[block], queryAtomCount = blocks.queryCount[block];
    int targetStart = blocks.targetStart[block], targetAtomCount = blocks.targetCount[block];
    for (int row = 0; row < queryAtomCount; row++) {
      int queryAtom = queryAtoms[queryStart + row];
      long columns = 0;
      for (int column = 0; column < targetAtomCount; column++) {
        int targetAtom = targetAtoms[targetStart + column];
        if (has(compatible, queryAtom, targetAtom) && bondsFit(queryAtom, targetAtom))
          columns |= 1L << column;
      }
      matchRows[row] = columns;
    }
    Arrays.fill(matchOwner, 0, targetAtomCount, -1);
    int pairs = 0;
    for (int row = 0; row < queryAtomCount; row++) {
      matchSeen = 0;
      if (augment(row))
        pairs++;
    }
    return pairs;
  }

  /**
   * Whether row {@code row} of {@link #matchRows} gets a column, by an augmenting path through
   * the columns not yet seen. A path has at most one step for each row, at most 64.
   */
  private boolean augment(int row) {
    long columns = matchRows[row] & ~matchSeen;
    while (columns != 0) {
      clock.tick();
      int column = Long.numberOfTrailingZeros(columns);
      columns &= columns - 1;
      matchSeen |= 1L << column;
      if (matchOwner[column] < 0 || augment(matchOwner[column])) {
        matchOwner[column] = row;
        return true;
      }
    }
    return false;
  }

  /**
   * Whether every atom pair of {@code block} matches exactly, the atoms and the bonds to the mapped
   * atoms, which holds when bonds of one class always match and otherwise is checked for blocks
   * of at most {@link #COMPLETE_PAIRS} pairs.
   */
  private boolean blockComplete(int block) {
    if (!graph.atomGroups.complete[blocks.atomClass[block]])
      return false;
    if (bondClassesMatch || blocks.links[block] == 0)
      return true;
    if (blocks.complete[block] >= 0)
      return blocks.complete[block] == 1;
    blocks.complete[block] = 0;
    int queryStart = blocks.queryStart[block], queryAtomCount = blocks.queryCount[block];
    int targetStart = blocks.targetStart[block], targetAtomCount = blocks.targetCount[block];
    if ((long) queryAtomCount * targetAtomCount > COMPLETE_PAIRS)
      return false;
    for (int queryPlace = queryStart; queryPlace < queryStart + queryAtomCount; queryPlace++)
      for (int targetPlace = targetStart; targetPlace < targetStart + targetAtomCount; targetPlace++)
        if (!bondsFit(queryAtoms[queryPlace], targetAtoms[targetPlace]))
          return false;
    blocks.complete[block] = 1;
    return true;
  }

  /** The score of a mapping of {@code atomCount} atoms with {@code bondCount} bonds between them, which compares by atoms first. */
  private static long score(int atomCount, int bondCount) {
    return (long) atomCount << 32 | bondCount;
  }

  /** The atoms of a score of this search. */
  private static int atoms(long score) {
    return (int) (score >>> 32);
  }

  /** The bonds of a score of this search. */
  private static int bonds(long score) {
    return (int) score;
  }

  /** The score of {@code mapping}: its atoms, and the bonds between them. */
  private long score(int[] mapping) {
    int atomCount = 0, bondCount = 0;
    for (int queryAtom = 0; queryAtom < queryCount; queryAtom++) {
      if (mapping[queryAtom] < 0)
        continue;
      atomCount++;
      for (int place = query.start[queryAtom]; place < query.start[queryAtom + 1]; place++)
        if (query.neighbour[place] < queryAtom && mapping[query.neighbour[place]] >= 0)
          bondCount++;
    }
    return score(atomCount, bondCount);
  }

  /** A bound capped by the upper bound, which caps the atoms and the bonds of every mapping apart. */
  private long capped(int atomCount, int bondCount) {
    return score(Math.min(atomCount, atoms(upperBound)), Math.min(bondCount, bonds(upperBound)));
  }

  /** Whether a state with this bound can hold a better mapping. */
  private boolean canImprove(long bound) {
    return bound > bestScore;
  }

  /**
   * The twins of one molecule, each linked to the next in a ring, -1 for none. Twins are
   * terminal atoms with the same neighbour that match the same atoms by the same bonds.
   * Swapping two takes every mapping onto one of the same score.
   */
  private int[] twins(CommonBondGraphs.Side side, boolean onQuery) {
    int[] next = new int[side.atomCount];
    Arrays.fill(next, -1);
    for (int parent = 0; parent < side.atomCount; parent++) {
      for (int place = side.start[parent]; place < side.start[parent + 1]; place++) {
        int first = side.neighbour[place], last = first;
        if (next[first] >= 0 || side.degree(first) != 1)
          continue;
        for (int laterPlace = place + 1; laterPlace < side.start[parent + 1]; laterPlace++) {
          int other = side.neighbour[laterPlace];
          if (next[other] < 0 && side.degree(other) == 1
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

  /** The first atom of the ring of twins of each atom, its lowest index, or -1 for an atom with no twin. */
  private static int[] ringFirsts(int[] twin) {
    int[] first = new int[twin.length];
    Arrays.fill(first, -1);
    for (int atom = 0; atom < twin.length; atom++) {
      if (twin[atom] < 0 || first[atom] >= 0)
        continue;
      // The ring is met from its lowest atom first.
      for (int member = atom; first[member] < 0; member = twin[member]) first[member] = atom;
    }
    return first;
  }

  /** Whether {@code atom} and {@code other}, with their bonds {@code bond} and {@code otherBond}, match alike. */
  private boolean alike(boolean onQuery, int atom, int bond, int other, int otherBond) {
    clock.tick(words + queryCount + query.bondCount + target.bondCount);
    if (onQuery) {
      for (int word = 0; word < words; word++)
        if (compatible[atom * words + word] != compatible[other * words + word])
          return false;
      for (int targetBond = 0; targetBond < target.bondCount; targetBond++)
        if (graph.bondsMatch(bond, targetBond) != graph.bondsMatch(otherBond, targetBond))
          return false;
    } else {
      for (int queryAtom = 0; queryAtom < queryCount; queryAtom++)
        if (has(compatible, queryAtom, atom) != has(compatible, queryAtom, other))
          return false;
      for (int queryBond = 0; queryBond < query.bondCount; queryBond++)
        if (graph.bondsMatch(queryBond, bond) != graph.bondsMatch(queryBond, otherBond))
          return false;
    }
    return true;
  }

  /**
   * The pair of {@code queryAtom} and {@code targetAtom} as one int,
   * {@code queryAtom * targetCount + targetAtom}, which fits as at most 2<sup>30</sup> pairs are
   * searched.
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
   * A mapped pair of the branch and bound, one for each depth, or the base level with nothing
   * mapped, with what it restores when taken off and what it has left to try. The levels are made
   * once and reused, see {@link #push}.
   */
  private static final class Level {

    /** The mapped pair, -1 at the base level. */
    int queryAtom, targetAtom;
    /** The level of the root of the open fragment, or one at the base level, where no fragment is open. */
    int fragmentStart;
    /** The bonds between mapped atoms before the pair was mapped, put back when it is taken off. */
    int commonBondsBefore;
    /** The query and target bonds with no end mapped, once the pair is mapped. */
    int freeQueryBonds, freeTargetBonds;
    /** The blocks of the level, from {@code blocksFrom} up to {@code blocksTo}, see {@link Blocks}. */
    int blocksFrom, blocksTo;
    /** A score that no mapping below the level beats, see {@link CommonAtomSearch#bound}. */
    long bound;
    /** A number of its own, which tells the candidates tried at this level, see {@link CommonAtomSearch#twinTried}. */
    long serial;
    /** The candidate being tried, or -1, and its block. */
    int triedAtom, triedBlock;
    /** Whether the candidate roots a new fragment, and whether its pair is forced, see {@link CommonAtomSearch#forcedPair}. */
    boolean triedRoot, forced;
    /** The targets of {@link #triedAtom}, the first {@link #targetCount}, see {@link CommonAtomSearch#targets}. */
    int[] targets = NONE;
    /** The number of targets, and the place of the next target to try. */
    int targetCount, nextTarget;
  }

  /**
   * The pool of blocks of every level on the stack, those of each level after those of its
   * parent. A block holds runs of the query and of the target atoms: the live atoms of each side
   * from its start, and as many as its count. Its links are the mapped atoms its atoms are bonded
   * to, it is open when these include an atom of the open fragment, and all its atoms have one atom
   * class. A level that leaves an atom unmapped moves it after the live atoms of its block, so the
   * runs of its parent keep the same atoms. The blocks of a level are in the order of their runs on
   * either side.
   */
  private static final class Blocks {

    /** The start and the number of live atoms of each side of each block. */
    int[] queryStart = new int[INITIAL_BLOCKS], queryCount = new int[INITIAL_BLOCKS],
          targetStart = new int[INITIAL_BLOCKS], targetCount = new int[INITIAL_BLOCKS];
    /** The number of mapped atoms the atoms of each block are bonded to, and their atom class. */
    int[] links = new int[INITIAL_BLOCKS], atomClass = new int[INITIAL_BLOCKS];
    /** Whether each block is bonded to the open fragment. */
    boolean[] open = new boolean[INITIAL_BLOCKS];
    /**
     * What is known of each block, which depends only on its atoms: its exact pairs, see
     * {@link CommonAtomSearch#exactPairs}, and whether all its pairs match, 1 or 0, see
     * {@link CommonAtomSearch#blockComplete}; -1 while unknown.
     */
    int[] exactPairs = new int[INITIAL_BLOCKS], complete = new int[INITIAL_BLOCKS];
    /** The number of blocks in use. */
    int count;

    /** Adds a block. */
    void add(int queryFrom, int queryAtoms, int targetFrom, int targetAtoms, int linked,
        boolean isOpen, int ofClass) {
      if (count == queryStart.length) {
        int capacity = 2 * count;
        queryStart = Arrays.copyOf(queryStart, capacity);
        queryCount = Arrays.copyOf(queryCount, capacity);
        targetStart = Arrays.copyOf(targetStart, capacity);
        targetCount = Arrays.copyOf(targetCount, capacity);
        links = Arrays.copyOf(links, capacity);
        atomClass = Arrays.copyOf(atomClass, capacity);
        open = Arrays.copyOf(open, capacity);
        exactPairs = Arrays.copyOf(exactPairs, capacity);
        complete = Arrays.copyOf(complete, capacity);
      }
      queryStart[count] = queryFrom;
      queryCount[count] = queryAtoms;
      targetStart[count] = targetFrom;
      targetCount[count] = targetAtoms;
      links[count] = linked;
      open[count] = isOpen;
      atomClass[count] = ofClass;
      exactPairs[count] = complete[count] = -1;
      count++;
    }
  }

  /** Whether bit {@code bit} of row {@code row} is set, in rows of {@link #words} longs as in {@link #compatible}. */
  private boolean has(long[] rows, int row, int bit) {
    return (rows[row * words + (bit >>> 6)] & 1L << bit) != 0;
  }

  /** Sets bit {@code bit} of row {@code row}. */
  private void set(long[] rows, int row, int bit) {
    rows[row * words + (bit >>> 6)] |= 1L << bit;
  }
}
