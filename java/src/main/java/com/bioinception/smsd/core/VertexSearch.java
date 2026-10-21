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
  private final int[][] domains;
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
  // Unweighted, uncapped single-result fragment searches take the native order; connected,
  // weighted, matchAll and node-capped searches take the dynamic smallest-domain order.
  private final boolean nativeOrder;

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
    if (control.stopped()) {
      domains = new int[query.n][];
      queryEdges = Collections.emptyList();
      targetEdges = Collections.emptyList();
      availableTargetBits = null;
      return;
    }
    if (options.chemistry.ringFusionMode == ChemOptions.RingFusionMode.STRICT) {
      query.ensureRingCounts();
      target.ensureRingCounts();
    }
    domains = new int[query.n][];
    Arrays.fill(domains, new int[0]);
    for (int queryAtom = 0; queryAtom < query.n && !control.stopped(); queryAtom++) {
      int[] candidates = new int[target.n];
      int count = 0;
      for (int targetAtom = 0; targetAtom < target.n; targetAtom++)
        if (!options.excludedTargetAtoms.contains(targetAtom)
            && MolGraph.ChemOps.atomsCompatFast(
                query, queryAtom, target, targetAtom, options.chemistry))
          candidates[count++] = targetAtom;
      for (int index = 0; index < count; index++) {
        int targetAtom = candidates[index];
        domainBits[queryAtom][targetAtom >>> 6] |= 1L << (targetAtom & 63);
      }
      domains[queryAtom] = Arrays.copyOf(candidates, count);
    }
    queryEdges = edges(query);
    targetEdges = edges(target);
    availableTargetBits = new long[target.words];
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

  void seed(SearchResult.Match match) {
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
                && (domains[queryAtom].length < domains[chosen].length
                    || (domains[queryAtom].length == domains[chosen].length
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
      int[] order = orderedCandidates(chosen);
      boolean stop = false;
      for (int position = 0; position < order[target.n] && !stop; position++)
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
    atoms(mappedAtoms + 1, mappedBonds);
    mappedCount--;
    usedTargetBits[word] &= ~(1L << (targetAtom & 63));
    targetToQuery[targetAtom] = -1;
    queryToTarget[chosen] = -1;
    return upperBoundReached || control.stopped()
        || (probeNodeLimit != 0 && nodes >= probeNodeLimit);
  }

  /**
   * Try first the targets that leave more unmapped query neighbours extendable, then those
   * whose degree is closest to the query atom's (after the native ordering, which also counts
   * mapped neighbours). The candidate count is stored in the last slot.
   */
  private int[] orderedCandidates(int chosen) {
    if (candidateOrder == null) {
      candidateOrder = new int[query.n + 1][target.n + 1];
      candidateSupport = new int[target.n];
    }
    int[] order = candidateOrder[frameDepth++];
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
    order[target.n] = size;
    return order;
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
    if ((options.chemistry.useChirality || options.chemistry.useBondStereo)
        && Substructure.stereoMismatchAtom(query, target, queryToTarget, options.chemistry)
            >= 0)
      return;
    List<Set<Integer>> components = MatchGraph.fragments(atoms, bonds.keySet());
    if (options.connectivity == MCS.Connectivity.Connected && components.size() != 1)
      return;
    if (components.size() > options.maxFragments)
      return;
    for (Set<Integer> component : components)
      if (component.size() < options.minFragmentAtoms)
        return;
    if (options.chemistry.completeRingsOnly
        && (!completeRings(query, atoms.keySet(), bonds.keySet())
            || !completeRings(
                target, new HashSet<>(atoms.values()), new HashSet<>(bonds.values()))))
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
