/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
package com.bioinception.smsd.core;

import java.util.*;

/**
 * The general MCIS search, with per-invocation state; prepared queries do not retain mutable
 * search data. {@link MCS} sends here every MCIS that the common-atom search does not take, and
 * every MCES goes to {@link CommonBondSearch}.
 */
final class VertexSearch {
  private final MolGraph query, target;
  private final MCS.Options options;
  private final SearchResult.Control control;
  private final boolean all, connected;
  private final int limit;
  private final int[] queryToTarget, targetToQuery;
  private final int[] domainSize;
  private final long[][] domainBits;
  private final long[] usedTargetBits, availableTargetBits;
  private final boolean[] availableQuery, availableTarget, reachableQuery, reachableTarget;
  private final int[] queryQueue, targetQueue;
  private final int[] mappedQueryAtoms;
  private int mappedCount;
  private final double[] weights;
  private final List<SearchResult.Bond> queryEdges, targetEdges;
  private final boolean[] excludedAtoms;
  private final List<SearchResult.Match> best = new ArrayList<>();
  private final Set<String> bestKeys = new HashSet<>();
  private double bestWeight;
  private long nodes;
  private int bestAtoms, bestBonds, initialAtoms, initialBonds;
  private boolean overflow, upperBoundReached, allowSkipping = true;
  private long probeNodeLimit;
  private int[][] candidateOrder;
  private int[] candidateSupport;
  private int frameDepth;
  private RingIndex queryRings, targetRings;
  /** The stamp of the last complete-ring test. Package-private so that tests can check the wrap. */
  int ringStamp;
  // Unweighted, uncapped single-result fragment searches take the native order; connected,
  // weighted, matchAll and node-capped searches take the dynamic smallest-domain order.
  private final boolean nativeOrder;
  /** Whether each new pair is tested against the centres and double bonds it touches. */
  private final boolean centrePairs, doubleBondPairs;
  /** The tetrahedral ligands of the tagged atoms, or null; see {@link #centreBroken}. */
  private final int[][] queryLigands, targetLigands;
  private final int[] permutation = new int[4];

  VertexSearch(MolGraph query, MolGraph target, MCS.Options options, int limit, boolean all) {
    this(query, target, options, limit, all, null);
  }

  VertexSearch(MolGraph query, MolGraph target, MCS.Options options, int limit, boolean all,
      SearchResult.Control sharedControl) {
    this.query = query;
    this.target = target;
    this.options = options;
    connected = options.connectivity == MCS.Connectivity.Connected || options.maxFragments == 1;
    nativeOrder = !connected && !all && options.atomWeights == null && options.nodeLimit == 0;
    this.limit = limit;
    this.all = all;
    control = sharedControl != null
        ? sharedControl
        : new SearchResult.Control(
              SearchResult.Control.resolveTimeout(options.timeoutMillis, query.n, target.n, true),
              options.cancellation);
    queryToTarget = new int[query.n];
    Arrays.fill(queryToTarget, -1);
    targetToQuery = new int[target.n];
    Arrays.fill(targetToQuery, -1);
    excludedAtoms = new boolean[query.n];
    domainSize = new int[query.n];
    domainBits = new long[query.n][target.words];
    usedTargetBits = new long[target.words];
    availableQuery = new boolean[query.n];
    availableTarget = new boolean[target.n];
    reachableQuery = new boolean[query.n];
    reachableTarget = new boolean[target.n];
    queryQueue = new int[query.n];
    targetQueue = new int[target.n];
    mappedQueryAtoms = new int[query.n];
    weights = new double[query.n];
    for (int atom = 0; atom < query.n; atom++)
      weights[atom] = options.atomWeights == null ? 1.0 : options.atomWeights[atom];
    centrePairs = options.chemistry.useChirality && tagged(query) && tagged(target);
    queryLigands = centrePairs ? ligands(query) : null;
    targetLigands = centrePairs ? ligands(target) : null;
    doubleBondPairs = options.chemistry.useBondStereo && configured(query) && configured(target);
    if (control.stopped()) {
      queryEdges = Collections.emptyList();
      targetEdges = Collections.emptyList();
      availableTargetBits = null;
      return;
    }
    if (options.chemistry.ringFusionMode == ChemOptions.RingFusionMode.STRICT) {
      query.ensureRingCounts();
      target.ensureRingCounts();
    }
    if (options.chemistry.completeRingsOnly && !query.ringsTruncated() && !target.ringsTruncated()) {
      queryRings = new RingIndex(query);
      targetRings = new RingIndex(target);
    }
    boolean[] excludedTargets = new boolean[target.n];
    for (Integer targetAtom : options.excludedTargetAtoms)
      if (targetAtom != null && targetAtom >= 0 && targetAtom < target.n)
        excludedTargets[targetAtom] = true;
    for (int queryAtom = 0; queryAtom < query.n && !control.stopped(); queryAtom++) {
      long[] domain = domainBits[queryAtom];
      int count = 0;
      for (int targetAtom = 0; targetAtom < target.n; targetAtom++)
        if (!excludedTargets[targetAtom]
            && (queryRings == null
                || queryRings.ringAtom(query, queryAtom) == targetRings.ringAtom(target, targetAtom))
            && MolGraph.ChemOps.atomsCompatFast(
                query, queryAtom, target, targetAtom, options.chemistry)) {
          domain[targetAtom >>> 6] |= 1L << (targetAtom & 63);
          count++;
        }
      domainSize[queryAtom] = count;
    }
    queryEdges = edges(query);
    targetEdges = edges(target);
    availableTargetBits = new long[target.words];
  }

  private static boolean tagged(MolGraph graph) {
    if (graph.tetraChirality != null)
      for (int tag : graph.tetraChirality)
        if (tag != 0)
          return true;
    return false;
  }

  private static int[][] ligands(MolGraph graph) {
    int[][] result = new int[graph.n][];
    for (int atom = 0; atom < graph.n; atom++)
      if (graph.tetraChirality[atom] != 0)
        result[atom] = graph.tetrahedralLigands(atom);
    return result;
  }

  private static boolean configured(MolGraph graph) {
    for (int atom = 0; atom < graph.n; atom++)
      for (int neighbor : graph.neighbors[atom])
        if (neighbor > atom && graph.dbStereo(atom, neighbor) != 0)
          return true;
    return false;
  }

  SearchResult search() {
    if (!control.stopped()) {
      // As in the native search, first try to embed the whole query; an induced embedding of
      // every query atom is optimal and gives the fragment search an incumbent.
      if (nativeOrder && initialAtoms == 0 && query.n > 0 && query.n <= target.n) {
        allowSkipping = false;
        probeNodeLimit = nodes + 2048;
        atoms(0, 0);
        allowSkipping = true;
        probeNodeLimit = 0;
      }
      atoms(initialAtoms, initialBonds);
    }
    boolean proved = control.status() == SearchResult.Status.Completed;
    if (proved && overflow && all)
      control.limitReached();
    double score = best.isEmpty() ? 0.0 : best.get(0).score();
    return new SearchResult(best, control.status(), proved, proved && !(overflow && all), score,
        nodes, control.elapsedMillis());
  }

  /** Start the search from the atom mapping of {@code match}, in place of any earlier seed. */
  void seed(SearchResult.Match match) {
    while (mappedCount > 0) {
      int queryAtom = mappedQueryAtoms[--mappedCount], targetAtom = queryToTarget[queryAtom];
      queryToTarget[queryAtom] = -1;
      targetToQuery[targetAtom] = -1;
      usedTargetBits[targetAtom >>> 6] &= ~(1L << (targetAtom & 63));
    }
    match.atomMapping().forEach((queryAtom, targetAtom) -> {
      queryToTarget[queryAtom] = targetAtom;
      targetToQuery[targetAtom] = queryAtom;
      mappedQueryAtoms[mappedCount++] = queryAtom;
      usedTargetBits[targetAtom >>> 6] |= 1L << (targetAtom & 63);
    });
    initialAtoms = match.atomCount();
    initialBonds = match.bondCount();
  }

  private boolean visit() {
    if (upperBoundReached || control.stopped())
      return false;
    if (probeNodeLimit != 0 && nodes >= probeNodeLimit)
      return false;
    if (options.nodeLimit != 0 && nodes >= options.nodeLimit) {
      control.limitReached();
      return false;
    }
    nodes++;
    return true;
  }

  private static List<SearchResult.Bond> edges(MolGraph graph) {
    List<SearchResult.Bond> result = new ArrayList<>();
    for (int atom = 0; atom < graph.n; atom++)
      for (int neighbor : graph.neighbors[atom])
        if (atom < neighbor)
          result.add(new SearchResult.Bond(atom, neighbor));
    return result;
  }

  private long feasibleTargets(int queryAtom, int word) {
    long candidates = domainBits[queryAtom][word] & ~usedTargetBits[word];
    boolean hasMappedNeighbor = false;
    for (int index = 0; index < mappedCount && candidates != 0; index++) {
      int mappedQuery = mappedQueryAtoms[index];
      long adjacent = target.adjLong[queryToTarget[mappedQuery]][word];
      if (query.bondOrder(queryAtom, mappedQuery) != 0) {
        candidates &= adjacent;
        hasMappedNeighbor = true;
      } else {
        candidates &= ~adjacent;
      }
    }
    if (!hasMappedNeighbor)
      return candidates;
    long remaining = candidates;
    while (remaining != 0) {
      int offset = Long.numberOfTrailingZeros(remaining);
      int targetAtom = (word << 6) + offset;
      remaining &= remaining - 1;
      for (int neighbor : query.neighbors[queryAtom]) {
        int mappedTarget = queryToTarget[neighbor];
        if (mappedTarget >= 0 && !MolGraph.ChemOps.bondsCompatible(
                query, queryAtom, neighbor, target, targetAtom, mappedTarget, options.chemistry)) {
          candidates &= ~(1L << offset);
          break;
        }
      }
    }
    return candidates;
  }

  private void atoms(int mappedAtoms, int mappedBonds) {
    if (!visit())
      return;
    int remaining = 0;
    Arrays.fill(availableTargetBits, 0L);
    for (int queryAtom = 0; queryAtom < query.n; queryAtom++)
      availableQuery[queryAtom] = queryToTarget[queryAtom] >= 0;

    int chosen = -1, smallestDomain = Integer.MAX_VALUE, mostNeighbors = -1;
    for (int queryAtom = 0; queryAtom < query.n; queryAtom++) {
      if (excludedAtoms[queryAtom] || queryToTarget[queryAtom] >= 0)
        continue;
      int count = 0;
      for (int word = 0; word < target.words; word++) {
        long candidates = feasibleTargets(queryAtom, word);
        availableTargetBits[word] |= candidates;
        count += Long.bitCount(candidates);
      }
      if (count > 0) {
        remaining++;
        availableQuery[queryAtom] = true;
      }
      int mappedNeighbors = 0;
      for (int neighbor : query.neighbors[queryAtom])
        if (queryToTarget[neighbor] >= 0)
          mappedNeighbors++;
      if (connected && mappedAtoms > 0
          && mappedNeighbors == 0)
        continue;
      if (nativeOrder) {
        // As in the native search: grow mapped pieces first, then take the smallest
        // static domain, then the highest degree.
        if (count == 0)
          continue;
        boolean frontier = mappedNeighbors > 0, chosenFrontier = mostNeighbors > 0;
        if (chosen < 0 || (frontier && !chosenFrontier)
            || (frontier == chosenFrontier
                && (domainSize[queryAtom] < domainSize[chosen]
                    || (domainSize[queryAtom] == domainSize[chosen]
                        && query.neighbors[queryAtom].length > query.neighbors[chosen].length)))) {
          chosen = queryAtom;
          mostNeighbors = mappedNeighbors;
        }
        continue;
      }
      boolean preferred = count < smallestDomain
          || (count == smallestDomain && mappedNeighbors > mostNeighbors);
      // Weighted fragment searches grow mapped pieces first, then break ties by degree.
      if (!connected && !all && options.nodeLimit == 0 && chosen >= 0) {
        boolean frontier = mappedNeighbors > 0;
        boolean chosenFrontier = mostNeighbors > 0;
        if (frontier != chosenFrontier)
          preferred = frontier;
        else if (count == smallestDomain && mappedNeighbors == mostNeighbors)
          preferred = query.degree[queryAtom] > query.degree[chosen];
      }
      if (chosen < 0 || preferred) {
        chosen = queryAtom;
        smallestDomain = count;
        mostNeighbors = mappedNeighbors;
      }
    }
    if (!allowSkipping && mappedAtoms + remaining < query.n)
      return;
    int targetCapacity = 0;
    for (long candidates : availableTargetBits)
      targetCapacity += Long.bitCount(candidates);
    int atomBound = mappedAtoms + Math.min(remaining, targetCapacity);
    for (int targetAtom = 0; targetAtom < target.n; targetAtom++)
      availableTarget[targetAtom] = targetToQuery[targetAtom] >= 0
          || (availableTargetBits[targetAtom >>> 6] & (1L << (targetAtom & 63))) != 0;
    // Feasible domains only shrink; every connected extension stays in these components.
    if (connected) {
      int queryRoot = mappedCount == 0 ? -1 : mappedQueryAtoms[0];
      int targetRoot = mappedCount == 0 ? -1 : queryToTarget[queryRoot];
      atomBound = Math.min(atomBound,
          connectedBound(query, availableQuery, reachableQuery, queryQueue, queryRoot));
      atomBound = Math.min(atomBound,
          connectedBound(target, availableTarget, reachableTarget, targetQueue, targetRoot));
    }
    int queryBondBound = 0, targetBondBound = 0;
    for (int edge = 0; edge < queryEdges.size(); edge++) {
      SearchResult.Bond bond = queryEdges.get(edge);
      if (availableQuery[bond.first()] && availableQuery[bond.second()]
          && (!connected
              || (reachableQuery[bond.first()] && reachableQuery[bond.second()])))
        queryBondBound++;
    }
    for (int edge = 0; edge < targetEdges.size(); edge++) {
      SearchResult.Bond bond = targetEdges.get(edge);
      if (availableTarget[bond.first()] && availableTarget[bond.second()]
          && (!connected
              || (reachableTarget[bond.first()] && reachableTarget[bond.second()])))
        targetBondBound++;
    }
    int bondBound = (int) Math.min(Math.min(queryBondBound, targetBondBound),
        (long) atomBound * (atomBound - 1) / 2);
    double weightBound = options.atomWeights == null ? 0.0 : weightBound();
    if (!Double.isNaN(weightBound)) {
      int boundOrder = compare(atomBound, bondBound, weightBound);
      if (boundOrder < 0 || (!all && boundOrder == 0 && !best.isEmpty()))
        return;
    }
    consider(mappedAtoms, mappedBonds,
        options.atomWeights == null ? mappedAtoms : orderedWeight());
    if (chosen < 0)
      return;
    int addedBonds = 0;
    for (int neighbor : query.neighbors[chosen])
      if (queryToTarget[neighbor] >= 0)
        addedBonds++;
    if (nativeOrder) {
      int depth = frameDepth, size = orderCandidates(chosen);
      int[] order = candidateOrder[depth];
      boolean stop = false;
      for (int position = 0; position < size && !stop; position++)
        stop = branch(chosen, order[position], mappedAtoms, mappedBonds + addedBonds);
      frameDepth--;
      if (stop)
        return;
    } else {
      for (int word = 0; word < target.words; word++) {
        long candidates = feasibleTargets(chosen, word);
        while (candidates != 0) {
          int targetAtom = (word << 6) + Long.numberOfTrailingZeros(candidates);
          candidates &= candidates - 1;
          if (branch(chosen, targetAtom, mappedAtoms, mappedBonds + addedBonds))
            return;
        }
      }
    }
    if (!allowSkipping)
      return;
    excludedAtoms[chosen] = true;
    atoms(mappedAtoms, mappedBonds);
    excludedAtoms[chosen] = false;
  }

  private boolean branch(int chosen, int targetAtom, int mappedAtoms, int mappedBonds) {
    int word = targetAtom >>> 6;
    queryToTarget[chosen] = targetAtom;
    targetToQuery[targetAtom] = chosen;
    usedTargetBits[word] |= 1L << (targetAtom & 63);
    mappedQueryAtoms[mappedCount++] = chosen;
    if (!stereoBroken(chosen))
      atoms(mappedAtoms + 1, mappedBonds);
    mappedCount--;
    usedTargetBits[word] &= ~(1L << (targetAtom & 63));
    targetToQuery[targetAtom] = -1;
    queryToTarget[chosen] = -1;
    return upperBoundReached || control.stopped()
        || (probeNodeLimit != 0 && nodes >= probeNodeLimit);
  }

  /**
   * Whether the stereo test of {@link #consider}, {@link Substructure#stereoMismatchAtom}, fails at
   * a centre or double bond that reads the newly mapped {@code atom}: the atom itself, its mapped
   * neighbours, and the double bonds at either. Mapped pairs stay in every extension, so such a
   * failure stays too, and the search backtracks at once.
   */
  private boolean stereoBroken(int atom) {
    if (centrePairs) {
      if (centreBroken(atom))
        return true;
      for (int neighbor : query.neighbors[atom])
        if (queryToTarget[neighbor] >= 0 && centreBroken(neighbor))
          return true;
    }
    if (doubleBondPairs)
      for (int neighbor : query.neighbors[atom]) {
        if (queryToTarget[neighbor] < 0)
          continue;
        if (doubleBondBroken(atom, neighbor))
          return true;
        for (int other : query.neighbors[neighbor])
          if (other != atom && queryToTarget[other] >= 0 && doubleBondBroken(neighbor, other))
            return true;
      }
    return false;
  }

  /**
   * The tetrahedral test at a mapped centre. Once three ligands are mapped the permutation is
   * fixed, as a fourth can only take the remaining target ligand, and a ligand whose image is not
   * a target ligand stays so.
   */
  private boolean centreBroken(int centre) {
    int targetAtom = queryToTarget[centre];
    int queryTag = query.tetraChirality[centre], targetTag = target.tetraChirality[targetAtom];
    int[] queryOrder = queryLigands[centre], targetOrder = targetLigands[targetAtom];
    if (queryTag == 0 || targetTag == 0 || queryOrder == null || targetOrder == null
        || queryOrder.length != 4 || targetOrder.length != 4)
      return false;
    int mapped = 0, usedSlots = 0;
    for (int slot = 0; slot < 4; slot++) {
      permutation[slot] = -1;
      int ligand = queryOrder[slot];
      if (ligand < 0 || queryToTarget[ligand] < 0)
        continue;
      int targetSlot = 0;
      while (targetSlot < 4 && targetOrder[targetSlot] != queryToTarget[ligand])
        targetSlot++;
      if (targetSlot == 4)
        return true;
      permutation[slot] = targetSlot;
      usedSlots |= 1 << targetSlot;
      mapped++;
    }
    if (mapped < 3)
      return false;
    for (int slot = 0; slot < 4; slot++)
      if (permutation[slot] < 0)
        permutation[slot] = Integer.numberOfTrailingZeros(~usedSlots);
    int inversions = 0;
    for (int first = 0; first < 4; first++)
      for (int second = first + 1; second < 4; second++)
        if (permutation[first] > permutation[second])
          inversions++;
    return (queryTag == targetTag) != ((inversions & 1) == 0);
  }

  /**
   * The double-bond test between two mapped atoms. It reads the first mapped substituent at each
   * end; since the mapping is induced, any other gives the same relation when an end has one
   * substituent in either graph, or at most two in both, and the bond is left to {@link
   * #consider} otherwise.
   */
  private boolean doubleBondBroken(int first, int second) {
    int queryConfiguration = query.dbStereo(first, second);
    if (queryConfiguration == 0)
      return false;
    int targetConfiguration = target.dbStereo(queryToTarget[first], queryToTarget[second]);
    if (targetConfiguration == 0)
      return false;
    for (int end = 0; end < 2; end++) {
      int queryEnd = end == 0 ? first : second, otherQueryEnd = end == 0 ? second : first;
      int targetEnd = queryToTarget[queryEnd], otherTargetEnd = queryToTarget[otherQueryEnd];
      int querySubstituents = query.neighbors[queryEnd].length - 1;
      int targetSubstituents = target.neighbors[targetEnd].length - 1;
      if (Math.min(querySubstituents, targetSubstituents) > 1
          && Math.max(querySubstituents, targetSubstituents) > 2)
        return false;
      int queryReference = query.doubleBondReference(queryEnd, otherQueryEnd);
      int targetReference = target.doubleBondReference(targetEnd, otherTargetEnd);
      int queryLigand = -1;
      for (int neighbor : query.neighbors[queryEnd])
        if (neighbor != otherQueryEnd && queryToTarget[neighbor] >= 0) {
          queryLigand = neighbor;
          break;
        }
      if (queryLigand < 0 || queryReference < 0 || targetReference < 0)
        return false;
      if ((queryLigand != queryReference) != (queryToTarget[queryLigand] != targetReference))
        queryConfiguration = 3 - queryConfiguration;
    }
    return queryConfiguration != targetConfiguration;
  }

  /**
   * Order the feasible targets of {@code chosen} in the row of the current depth: first those
   * that leave more unmapped query neighbours extendable, then those whose degree is closest to
   * the query atom's (after the native ordering, which also counts mapped neighbours). A row is
   * allocated when its depth is first reached and grows to the largest domain ordered there.
   *
   * @return the number of candidates
   */
  private int orderCandidates(int chosen) {
    if (candidateOrder == null) {
      candidateOrder = new int[query.n + 1][];
      candidateSupport = new int[target.n];
    }
    // The feasible targets are a subset of the domain.
    int[] order = candidateOrder[frameDepth];
    if (order == null || order.length < domainSize[chosen])
      candidateOrder[frameDepth] = order = new int[domainSize[chosen]];
    frameDepth++;
    int size = 0, degree = query.neighbors[chosen].length;
    for (int word = 0; word < target.words; word++) {
      long candidates = feasibleTargets(chosen, word);
      while (candidates != 0) {
        int candidate = (word << 6) + Long.numberOfTrailingZeros(candidates);
        candidates &= candidates - 1;
        int support = 0;
        for (int neighbor : query.neighbors[chosen]) {
          if (excludedAtoms[neighbor] || queryToTarget[neighbor] >= 0)
            continue;
          for (int adjacent : target.neighbors[candidate])
            if ((usedTargetBits[adjacent >>> 6] & (1L << (adjacent & 63))) == 0
                && (domainBits[neighbor][adjacent >>> 6] & (1L << (adjacent & 63))) != 0
                && MolGraph.ChemOps.bondsCompatible(
                    query, chosen, neighbor, target, candidate, adjacent, options.chemistry)) {
              support++;
              break;
            }
        }
        candidateSupport[candidate] = support;
        int position = size++;
        while (position > 0 && before(candidate, order[position - 1], degree)) {
          order[position] = order[position - 1];
          position--;
        }
        order[position] = candidate;
      }
    }
    return size;
  }

  private boolean before(int first, int second, int degree) {
    if (candidateSupport[first] != candidateSupport[second])
      return candidateSupport[first] > candidateSupport[second];
    return Math.abs(target.neighbors[first].length - degree)
        < Math.abs(target.neighbors[second].length - degree);
  }

  private static int connectedBound(
      MolGraph graph, boolean[] available, boolean[] reached, int[] queue, int root) {
    Arrays.fill(reached, false);
    int largest = 0;
    for (int start = root < 0 ? 0 : root; start < graph.n; start++) {
      if (!available[start] || reached[start])
        continue;
      int head = 0, tail = 0;
      queue[tail++] = start;
      reached[start] = true;
      while (head < tail) {
        int atom = queue[head++];
        for (int neighbor : graph.neighbors[atom])
          if (available[neighbor] && !reached[neighbor]) {
            reached[neighbor] = true;
            queue[tail++] = neighbor;
          }
      }
      if (root >= 0)
        return tail;
      largest = Math.max(largest, tail);
    }
    return largest;
  }

  private double orderedWeight() {
    double result = 0.0;
    for (int queryAtom = 0; queryAtom < query.n; queryAtom++)
      if (queryToTarget[queryAtom] >= 0)
        result += weights[queryAtom];
    return result;
  }

  private double weightBound() {
    double result = 0.0;
    for (int queryAtom = 0; queryAtom < query.n; queryAtom++) {
      double contribution = queryToTarget[queryAtom] >= 0 ? weights[queryAtom]
          : excludedAtoms[queryAtom]                      ? 0.0
                                                          : Math.max(0.0, weights[queryAtom]);
      result = Math.nextUp(result + contribution);
    }
    return result;
  }

  private int compare(int atoms, int bonds, double weight) {
    if (options.atomWeights != null) {
      int order = Double.compare(weight, bestWeight);
      if (order != 0)
        return order;
    }
    int order = Integer.compare(atoms, bestAtoms);
    return order != 0 ? order : Integer.compare(bonds, bestBonds);
  }

  private void consider(int atomCount, int bondCount, double weight) {
    int objectiveOrder = compare(atomCount, bondCount, weight);
    if (atomCount == 0 || objectiveOrder < 0 || (!all && objectiveOrder == 0 && !best.isEmpty()))
      return;
    // The stereo and ring rules read the mapping itself, so they reject before the maps are built.
    if ((options.chemistry.useChirality || options.chemistry.useBondStereo)
        && Substructure.stereoMismatchAtom(query, target, queryToTarget, options.chemistry)
            >= 0)
      return;
    if (options.chemistry.completeRingsOnly && !ringsComplete())
      return;
    Map<Integer, Integer> atoms = new TreeMap<>();
    for (int queryAtom = 0; queryAtom < query.n; queryAtom++)
      if (queryToTarget[queryAtom] >= 0)
        atoms.put(queryAtom, queryToTarget[queryAtom]);
    Map<SearchResult.Bond, SearchResult.Bond> bonds = new LinkedHashMap<>();
    for (int index = 0; index < queryEdges.size(); index++) {
      SearchResult.Bond bond = queryEdges.get(index);
      if (atoms.containsKey(bond.first()) && atoms.containsKey(bond.second()))
        bonds.put(bond,
            new SearchResult.Bond(queryToTarget[bond.first()], queryToTarget[bond.second()]));
    }
    List<Set<Integer>> components = MatchGraph.fragments(atoms, bonds.keySet());
    if (options.connectivity == MCS.Connectivity.Connected && components.size() != 1)
      return;
    if (components.size() > options.maxFragments)
      return;
    for (Set<Integer> component : components)
      if (component.size() < options.minFragmentAtoms)
        return;
    int comparison = compare(atomCount, bondCount, weight);
    if (comparison > 0) {
      best.clear();
      bestKeys.clear();
      overflow = false;
      bestAtoms = atomCount;
      bestBonds = bondCount;
      bestWeight = weight;
    }
    String key = atoms.toString() + bonds.keySet();
    if (!bestKeys.add(key))
      return;
    if (best.size() >= limit) {
      if (all)
        overflow = true;
      return;
    }
    double score = options.atomWeights == null ? atomCount : weight;
    best.add(new SearchResult.Match(atoms, bonds, components, score));
    if (!all && options.atomWeights == null && atomCount == Math.min(query.n, target.n)
        && bondCount == Math.min(queryEdges.size(), targetEdges.size()))
      upperBoundReached = true;
  }

  /**
   * {@link #completeRings} for both graphs over the current mapping, which selects the mapped
   * query atoms, every query bond between them, and their images in the target.
   */
  boolean ringsComplete() {
    if (queryRings == null) {
      queryRings = new RingIndex(query);
      targetRings = new RingIndex(target);
    }
    if (ringStamp == Integer.MAX_VALUE) {
      queryRings.clearStamps();
      targetRings.clearStamps();
      ringStamp = 0;
    }
    int stamp = ++ringStamp;
    // A query ring is selected whole when each of its bonds joins two mapped atoms.
    for (int ring = 0; ring < queryRings.cycles.length; ring++) {
      int[] cycle = queryRings.cycles[ring];
      boolean whole = true;
      for (int index = 0; index < cycle.length && whole; index++) {
        int first = cycle[index], second = cycle[(index + 1) % cycle.length];
        whole = queryToTarget[first] >= 0 && queryToTarget[second] >= 0
            && bonded(query, first, second);
      }
      if (whole)
        queryRings.cover(ring, stamp);
    }
    // A target ring is selected whole when each of its bonds is the image of a selected bond.
    for (int ring = 0; ring < targetRings.cycles.length; ring++) {
      int[] cycle = targetRings.cycles[ring];
      boolean whole = true;
      for (int index = 0; index < cycle.length && whole; index++) {
        int first = targetToQuery[cycle[index]];
        int second = targetToQuery[cycle[(index + 1) % cycle.length]];
        whole = first >= 0 && second >= 0 && bonded(query, first, second);
      }
      if (whole)
        targetRings.cover(ring, stamp);
    }
    for (int index = 0; index < mappedCount; index++) {
      int queryAtom = mappedQueryAtoms[index], targetAtom = queryToTarget[queryAtom];
      if (!queryRings.atomCovered(query, queryAtom, stamp)
          || !targetRings.atomCovered(target, targetAtom, stamp))
        return false;
      for (int neighbor : query.neighbors[queryAtom])
        if (neighbor > queryAtom && queryToTarget[neighbor] >= 0
            && (!queryRings.bondCovered(query, queryAtom, neighbor, stamp)
                || !targetRings.bondCovered(
                    target, targetAtom, queryToTarget[neighbor], stamp)))
          return false;
    }
    return true;
  }

  /** Whether the bond list of {@code graph} holds the bond between two atoms, as {@link #edges}. */
  private static boolean bonded(MolGraph graph, int first, int second) {
    int low = Math.min(first, second), high = Math.max(first, second);
    for (int neighbor : graph.neighbors[low])
      if (neighbor == high)
        return true;
    return false;
  }

  /**
   * The relevant cycles of one graph with their bonds numbered, for {@link #ringsComplete}. The
   * bonds are numbered by their lower atom, then their higher atom; a stamp marks the atoms and
   * bonds of the cycles selected whole in the current test.
   */
  static final class RingIndex {
    final int[][] cycles;
    /** The number of the bond from each cycle position to the next. */
    private final int[][] cycleBonds;
    private final boolean[] onCycle;
    /** The bonds whose lower atom is {@code atom} are numbered from {@code start[atom]}. */
    private final int[] start, higher;
    private final int[] atomStamp, bondStamp;

    RingIndex(MolGraph graph) {
      cycles = graph.computeRelevantCycles();
      onCycle = new boolean[graph.n];
      int positions = 0;
      for (int[] cycle : cycles)
        positions += cycle.length;
      long[] pairs = new long[positions];
      int count = 0;
      for (int[] cycle : cycles)
        for (int index = 0; index < cycle.length; index++) {
          int first = cycle[index], second = cycle[(index + 1) % cycle.length];
          onCycle[first] = true;
          pairs[count++] = pair(first, second);
        }
      Arrays.sort(pairs);
      int distinct = 0;
      for (int index = 0; index < count; index++)
        if (distinct == 0 || pairs[index] != pairs[distinct - 1])
          pairs[distinct++] = pairs[index];
      start = new int[graph.n + 1];
      higher = new int[distinct];
      for (int bond = 0; bond < distinct; bond++) {
        start[(int) (pairs[bond] >>> 32) + 1]++;
        higher[bond] = (int) pairs[bond];
      }
      for (int atom = 0; atom < graph.n; atom++)
        start[atom + 1] += start[atom];
      cycleBonds = new int[cycles.length][];
      for (int ring = 0; ring < cycles.length; ring++) {
        int[] cycle = cycles[ring];
        cycleBonds[ring] = new int[cycle.length];
        for (int index = 0; index < cycle.length; index++)
          cycleBonds[ring][index] = bond(cycle[index], cycle[(index + 1) % cycle.length]);
      }
      atomStamp = new int[graph.n];
      bondStamp = new int[distinct];
    }

    private static long pair(int first, int second) {
      return ((long) Math.min(first, second) << 32) | Math.max(first, second);
    }

    /** The number of the cycle bond between two atoms, or -1 when no cycle has that bond. */
    private int bond(int first, int second) {
      int low = Math.min(first, second), high = Math.max(first, second);
      for (int bond = start[low]; bond < start[low + 1]; bond++)
        if (higher[bond] == high)
          return bond;
      return -1;
    }

    void cover(int ring, int stamp) {
      for (int atom : cycles[ring])
        atomStamp[atom] = stamp;
      for (int bond : cycleBonds[ring])
        bondStamp[bond] = stamp;
    }

    /** Whether an atom is flagged or lies on a relevant cycle. */
    boolean ringAtom(MolGraph graph, int atom) {
      return graph.ring[atom] || onCycle[atom];
    }

    /** Whether a bond is flagged or lies on a relevant cycle. */
    boolean ringBond(MolGraph graph, int first, int second) {
      return graph.bondInRing(first, second) || bond(first, second) >= 0;
    }

    /** Whether a selected atom is off every ring, or on a ring selected whole. */
    boolean atomCovered(MolGraph graph, int atom, int stamp) {
      return !ringAtom(graph, atom) || atomStamp[atom] == stamp;
    }

    /** Whether a selected bond is off every ring, or on a ring selected whole. */
    boolean bondCovered(MolGraph graph, int first, int second, int stamp) {
      int bond = bond(first, second);
      return !(graph.bondInRing(Math.min(first, second), Math.max(first, second)) || bond >= 0)
          || (bond >= 0 && bondStamp[bond] == stamp);
    }

    void clearStamps() {
      Arrays.fill(atomStamp, 0);
      Arrays.fill(bondStamp, 0);
    }
  }

  /**
   * Whether every selected ring atom and bond of {@code graph} lies on a relevant ring that is
   * selected whole.
   */
  static boolean completeRings(
      MolGraph graph, Set<Integer> selectedAtoms, Set<SearchResult.Bond> selectedBonds) {
    Set<Integer> coveredAtoms = new HashSet<>();
    Set<SearchResult.Bond> coveredBonds = new HashSet<>();
    Set<Integer> ringAtoms = new HashSet<>();
    Set<SearchResult.Bond> ringBonds = new HashSet<>();
    for (int[] ring : graph.computeRelevantCycles()) {
      for (int index = 0; index < ring.length; index++) {
        ringAtoms.add(ring[index]);
        ringBonds.add(new SearchResult.Bond(ring[index], ring[(index + 1) % ring.length]));
      }
      boolean complete = true;
      for (int index = 0; index < ring.length; index++)
        if (!selectedBonds.contains(
                new SearchResult.Bond(ring[index], ring[(index + 1) % ring.length]))) {
          complete = false;
          break;
        }
      if (complete)
        for (int index = 0; index < ring.length; index++) {
          coveredAtoms.add(ring[index]);
          coveredBonds.add(new SearchResult.Bond(ring[index], ring[(index + 1) % ring.length]));
        }
    }
    for (int atom : selectedAtoms)
      if ((graph.ring[atom] || ringAtoms.contains(atom)) && !coveredAtoms.contains(atom))
        return false;
    for (SearchResult.Bond bond : selectedBonds)
      if ((graph.bondInRing(bond.first(), bond.second()) || ringBonds.contains(bond))
          && !coveredBonds.contains(bond))
        return false;
    return true;
  }
}
