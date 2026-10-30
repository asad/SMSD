/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
package com.bioinception.smsd.core;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * An atom bound for connected mappings from true bridges. Roots with only bridge bonds pair
 * their disjoint branches injectively; cyclic roots use part counts. Actual predicates are
 * used at roots and links, without assuming raw chemical label equality.
 */
final class BridgeBound {
  private static final int MAX_BRANCHES = 8, MAX_ATOMS = 512;
  private static final long MAX_PART_PAIRS = 250000, MAX_WORK = 50000000, MAX_COUNTS_PER_GRAPH = 1000000;

  private static final class Part {
    final int root, parent;
    int size;
    int[] branches;
    int[] atomCounts;

    Part(int root, int parent) {
      this.root = root;
      this.parent = parent;
    }
  }

  private static final class Parts {
    final List<Part> parts = new ArrayList<>();
    Integer[] order;
    boolean[] chainRoot;
    boolean acyclic = true, hasBridge;

    boolean build(MolGraph graph, int[] elements, SearchResult.Control control) {
      if (graph.n <= 0 || graph.n > MAX_ATOMS || graph.neighbors.length != graph.n)
        return false;
      int[] parent = new int[graph.n], queue = new int[graph.n];
      Arrays.fill(parent, -2);
      for (int atom = 0; atom < graph.n; atom++) {
        if (control.stopped() || graph.neighbors[atom].length > MAX_BRANCHES)
          return false;
        int[] neighbors = graph.neighbors[atom];
        for (int place = 0; place < neighbors.length; place++) {
          int neighbor = neighbors[place];
          if (neighbor < 0 || neighbor >= graph.n || neighbor == atom)
            return false;
          for (int earlier = 0; earlier < place; earlier++)
            if (neighbors[earlier] == neighbor)
              return false;
        }
        // Certify the topology actually read by the induced search, including bit adjacency.
        for (int other = 0; other < graph.n; other++) {
          boolean listed = contains(neighbors, other);
          if (graph.hasBond(atom, other) != listed || (graph.bondOrder(atom, other) != 0) != listed
              || graph.bondOrder(atom, other) != graph.bondOrder(other, atom))
            return false;
          if (listed && (!contains(graph.neighbors[other], atom)
              || graph.bondAromatic(atom, other) != graph.bondAromatic(other, atom)
              || graph.bondInRing(atom, other) != graph.bondInRing(other, atom)))
            return false;
        }
      }
      int[] entered = new int[graph.n], low = new int[graph.n], next = new int[graph.n];
      Arrays.fill(entered, -1);
      boolean[][] bridge = new boolean[graph.n][];
      for (int atom = 0; atom < graph.n; atom++)
        bridge[atom] = new boolean[graph.neighbors[atom].length];
      int clock = 0;
      for (int start = 0; start < graph.n; start++) {
        if (entered[start] >= 0)
          continue;
        int size = 1;
        queue[0] = start;
        parent[start] = -1;
        entered[start] = low[start] = clock++;
        while (size > 0) {
          if (control.stopped())
            return false;
          int atom = queue[size - 1];
          if (next[atom] < graph.neighbors[atom].length) {
            int neighbor = graph.neighbors[atom][next[atom]++];
            if (neighbor == parent[atom])
              continue;
            if (entered[neighbor] < 0) {
              parent[neighbor] = atom;
              entered[neighbor] = low[neighbor] = clock++;
              queue[size++] = neighbor;
            } else {
              low[atom] = Math.min(low[atom], entered[neighbor]);
            }
          } else {
            size--;
            int previous = parent[atom];
            if (previous < 0)
              continue;
            low[previous] = Math.min(low[previous], low[atom]);
            if (low[atom] > entered[previous]) {
              for (int place = 0; place < graph.neighbors[previous].length; place++)
                if (graph.neighbors[previous][place] == atom)
                  bridge[previous][place] = true;
              for (int place = 0; place < graph.neighbors[atom].length; place++)
                if (graph.neighbors[atom][place] == previous)
                  bridge[atom][place] = true;
              hasBridge = true;
            }
          }
        }
      }
      if (!hasBridge)
        for (int[] adjacent : graph.neighbors)
          if (adjacent.length != 0)
            return false;
      chainRoot = new boolean[graph.n];
      Arrays.fill(chainRoot, true);
      int[][] beyond = new int[graph.n][];
      for (int atom = 0; atom < graph.n; atom++)
        parts.add(new Part(atom, -1));
      for (int atom = 0; atom < graph.n; atom++) {
        beyond[atom] = new int[graph.neighbors[atom].length];
        Arrays.fill(beyond[atom], -1);
        for (int place = 0; place < graph.neighbors[atom].length; place++) {
          if (bridge[atom][place]) {
            beyond[atom][place] = parts.size();
            parts.add(new Part(graph.neighbors[atom][place], atom));
          } else {
            chainRoot[atom] = false;
            acyclic = false;
          }
        }
      }
      if ((long) elements.length * parts.size() > MAX_COUNTS_PER_GRAPH)
        return false;
      for (Part part : parts) {
        if (control.stopped())
          return false;
        int[] neighbors = graph.neighbors[part.root];
        int[] branches = new int[neighbors.length];
        int branchCount = 0;
        for (int place = 0; place < neighbors.length; place++)
          if (neighbors[place] != part.parent && beyond[part.root][place] >= 0)
            branches[branchCount++] = beyond[part.root][place];
        part.branches = Arrays.copyOf(branches, branchCount);
        part.atomCounts = new int[elements.length];
        Arrays.fill(parent, -2);
        parent[part.root] = -1;
        queue[0] = part.root;
        int size = 1;
        for (int head = 0; head < size; head++) {
          if ((head & 255) == 0 && control.stopped())
            return false;
          int atom = queue[head];
          if (elements.length != 0)
            part.atomCounts[Arrays.binarySearch(elements, graph.atomicNum[atom])]++;
          for (int neighbor : graph.neighbors[atom])
            if (!(atom == part.root && neighbor == part.parent) && parent[neighbor] == -2) {
              parent[neighbor] = atom;
              queue[size++] = neighbor;
            }
        }
        part.size = size;
      }
      order = new Integer[parts.size()];
      for (int part = 0; part < order.length; part++)
        order[part] = part;
      Arrays.sort(order, (first, second) -> Integer.compare(parts.get(first).size, parts.get(second).size));
      return true;
    }

    private static boolean contains(int[] values, int value) {
      for (int candidate : values)
        if (candidate == value)
          return true;
      return false;
    }
  }

  private final int targetParts;
  private final int[] pairBounds;
  private int upperBound;
  private boolean acyclic;

  private BridgeBound(int queryParts, int targetParts) {
    this.targetParts = targetParts;
    pairBounds = new int[queryParts * targetParts];
  }

  /** Decline graphs without useful bridges, high degree or resource limits; publish only a finished table. */
  static BridgeBound create(MolGraph query, MolGraph target, ChemOptions chemistry, SearchResult.Control control) {
    if (query.n <= 0 || target.n <= 0 || query.n > MAX_ATOMS || target.n > MAX_ATOMS)
      return null;
    int[] elements = new int[0];
    if (chemistry.matchAtomType) {
      elements = new int[query.n + target.n];
      System.arraycopy(query.atomicNum, 0, elements, 0, query.n);
      System.arraycopy(target.atomicNum, 0, elements, query.n, target.n);
      Arrays.sort(elements);
      int distinct = 0;
      for (int atom = 0; atom < elements.length; atom++)
        if (distinct == 0 || elements[atom] != elements[distinct - 1])
          elements[distinct++] = elements[atom];
      elements = Arrays.copyOf(elements, distinct);
    }
    Parts queryParts = new Parts(), targetParts = new Parts();
    if (!queryParts.build(query, elements, control) || !targetParts.build(target, elements, control)
        || (long) queryParts.parts.size() * targetParts.parts.size() > MAX_PART_PAIRS)
      return null;
    BridgeBound result = new BridgeBound(queryParts.parts.size(), targetParts.parts.size());
    result.acyclic = queryParts.acyclic && targetParts.acyclic;
    int[] best = new int[1 << MAX_BRANCHES];
    int[] weights = new int[MAX_BRANCHES];
    long work = 0;
    for (int queryIndex : queryParts.order) {
      for (int targetIndex : targetParts.order) {
        if (control.stopped())
          return null;
        Part queryPart = queryParts.parts.get(queryIndex), targetPart = targetParts.parts.get(targetIndex);
        if (!MolGraph.ChemOps.atomsCompatFast(query, queryPart.root, target, targetPart.root, chemistry))
          continue;
        if (elements.length + 1L > MAX_WORK - work)
          return null;
        work += elements.length + 1L;
        int bound = Math.min(queryPart.size, targetPart.size);
        if (elements.length != 0) {
          bound = 0;
          for (int element = 0; element < elements.length; element++)
            bound += Math.min(queryPart.atomCounts[element], targetPart.atomCounts[element]);
        }
        if (!queryParts.chainRoot[queryPart.root] || !targetParts.chainRoot[targetPart.root]) {
          result.pairBounds[queryIndex * result.targetParts + targetIndex] = bound;
          if (queryIndex < query.n && targetIndex < target.n)
            result.upperBound = Math.max(result.upperBound, bound);
          continue;
        }
        int sets = 1 << targetPart.branches.length;
        long units = (queryPart.branches.length + 1L) * sets * (targetPart.branches.length + 1);
        if (units > MAX_WORK - work)
          return null;
        work += units;
        Arrays.fill(best, 0, sets, -1);
        best[0] = 0;
        for (int queryBranch : queryPart.branches) {
          if (control.stopped())
            return null;
          int queryChild = queryParts.parts.get(queryBranch).root;
          Arrays.fill(weights, 0);
          for (int branch = 0; branch < targetPart.branches.length; branch++) {
            int targetBranch = targetPart.branches[branch];
            if (MolGraph.ChemOps.bondsCompatible(query, queryPart.root, queryChild,
                target, targetPart.root, targetParts.parts.get(targetBranch).root, chemistry))
              weights[branch] = result.pairBounds[queryBranch * result.targetParts + targetBranch];
          }
          // Supersets first: each query branch is paired at most once; skipping is retained.
          for (int used = sets - 1; used >= 0; used--) {
            if ((used & 255) == 255 && control.stopped())
              return null;
            if (best[used] < 0)
              continue;
            for (int branch = 0; branch < targetPart.branches.length; branch++) {
              if ((used & (1 << branch)) == 0 && weights[branch] != 0) {
                int with = used | (1 << branch);
                best[with] = Math.max(best[with], best[used] + weights[branch]);
              }
            }
          }
        }
        int branchBound = 1;
        for (int used = 0; used < sets; used++)
          branchBound = Math.max(branchBound, 1 + best[used]);
        bound = Math.min(bound, branchBound);
        result.pairBounds[queryIndex * result.targetParts + targetIndex] = bound;
        if (queryIndex < query.n && targetIndex < target.n)
          result.upperBound = Math.max(result.upperBound, bound);
      }
    }
    return result;
  }

  int upper() {
    return upperBound;
  }

  boolean acyclic() {
    return acyclic;
  }

  /** The full-component bound for a root pair, rather than a directed edge side. */
  int bound(int queryAtom, int targetAtom) {
    return pairBounds[queryAtom * targetParts + targetAtom];
  }
}
