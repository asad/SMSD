/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
package com.bioinception.smsd.core;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * Test only: the pairwise candidate search of release 8.1.0 as it stood before the search read
 * its candidates from one tree per atom, kept verbatim, with a replay that reports the charge of
 * every pair. The shared helpers it calls (ring systems, the bounded pipeline, the certificate and
 * the final selection) are unchanged.
 */
final class RingPerceptionReference {
  private RingPerceptionReference() {}

  static RingPerception.Rings perceive(int n, int[][] neighbors, boolean enumerateRelevant, long workLimit) {
    if (n == 0) return new RingPerception.Rings(new int[0][], new int[0][]);

    int cycleRank = RingPerception.cycleRank(n, neighbors);
    if (cycleRank <= 0) return new RingPerception.Rings(new int[0][], new int[0][]);

    // Candidates: for every vertex v and every edge (eU, eV),
    // find shortest path from v to eU and from v to eV (excluding edge eU-eV),
    // form a cycle from both paths + the edge. These candidates contain a minimum
    // cycle basis; relevant cycles need the separate all-shortest-path enumeration.

    // Build edge index
    Map<Long, Integer> edgeIdx = new HashMap<>();
    List<int[]> edgeList = new ArrayList<>();
    for (int i = 0; i < n; i++)
      for (int j : neighbors[i])
        if (j > i) {
          edgeIdx.put(MolGraph.bondKey(i, j), edgeList.size());
          edgeList.add(new int[]{i, j});
        }
    int totalEdges = edgeList.size();

    // Helper: BFS shortest path from start to target, avoiding edge (avoidU, avoidV)
    // Returns atom path, or empty list if unreachable
    List<boolean[]> candidateMasks = new ArrayList<>();
    List<int[]> candidateAtoms = new ArrayList<>();
    List<Integer> candidateSizes = new ArrayList<>();
    Set<BitSet> seenMasks = new HashSet<>();

    // Only pairs of an atom and a bond that can add a candidate are searched, in the same
    // order. A bridge leaves one of its atoms unreachable; from outside the ring system of a
    // bond both paths enter the system through the same atom; and in a system that is one
    // simple ring every pair repeats the ring that its lowest atom and first bond give.
    RingPerception.Systems systems = RingPerception.systems(n, neighbors);
    int[][] systemEdges = null;
    int[] lowestAtom = null;
    if (systems != null) {
      int count = systems.atoms.length;
      int[] sizes = new int[count];
      for (int[] edge : edgeList) {
        int system = systems.ofAtom[edge[0]];
        if (system >= 0 && systems.ofAtom[edge[1]] == system) sizes[system]++;
      }
      systemEdges = new int[count][];
      for (int system = 0; system < count; system++) systemEdges[system] = new int[sizes[system]];
      Arrays.fill(sizes, 0);
      for (int edge = 0; edge < totalEdges; edge++) {
        int system = systems.ofAtom[edgeList.get(edge)[0]];
        if (system >= 0 && systems.ofAtom[edgeList.get(edge)[1]] == system) systemEdges[system][sizes[system]++] = edge;
      }
      lowestAtom = new int[count];
      Arrays.fill(lowestAtom, -1);
      for (int atom = n - 1; atom >= 0; atom--)
        if (systems.ofAtom[atom] >= 0) lowestAtom[systems.ofAtom[atom]] = atom;
    }
    int[] allEdges = new int[totalEdges];
    for (int edge = 0; edge < totalEdges; edge++) allEdges[edge] = edge;
    PathSearch search = new PathSearch(n, neighbors, systems);

    for (int v = 0; v < n; v++) {
      int[] pairs = allEdges;
      if (systems != null) {
        int system = systems.ofAtom[v];
        pairs = system < 0 ? new int[0]
            : !systems.simpleRing[system] ? systemEdges[system]
            : v == lowestAtom[system] && systemEdges[system].length > 0 ? new int[] {systemEdges[system][0]}
            : new int[0];
      }
      for (int pair : pairs) {
        if (search.work > workLimit) {
          int[][][] bounded = RingPerception.boundedRings(n, neighbors, systems, cycleRank);
          if (!RingPerception.certifiedMinimumBasis(n, neighbors, bounded[0], cycleRank, workLimit))
            return new RingPerception.Rings(bounded[0], bounded[1], true);
          // The work bound stays recorded: fusion counts still come from the bounded candidates.
          return RingPerception.selected(n, neighbors, systems, cycleRank, bounded[0], bounded[1], enumerateRelevant, true);
        }
        int[] edge = edgeList.get(pair);
        int eU = edge[0], eV = edge[1];

        // Find both shortest paths while avoiding edge (eU, eV).
        List<Integer> path1 = new ArrayList<>(), path2 = new ArrayList<>();
        search.paths(v, eU, eV, eU, eV, path1, path2);
        if (path1.isEmpty()) continue;
        if (path2.isEmpty()) continue;
        search.work += path1.size() + path2.size();

        // Check paths are vertex-disjoint (only share vertex v)
        Set<Integer> p1Interior = new HashSet<>();
        for (int i = 1; i < path1.size(); i++) p1Interior.add(path1.get(i));
        boolean disjoint = true;
        for (int i = 1; i < path2.size(); i++) {
          if (p1Interior.contains(path2.get(i))) { disjoint = false; break; }
        }
        if (!disjoint) continue;

        // Build edge mask for this cycle
        boolean[] mask = new boolean[totalEdges];
        int edgeId = edgeIdx.get(MolGraph.bondKey(eU, eV));
        mask[edgeId] = true;

        boolean valid = true;
        for (int i = 0; i < path1.size() - 1 && valid; i++) {
          Integer idx = edgeIdx.get(MolGraph.bondKey(Math.min(path1.get(i), path1.get(i+1)),
                                            Math.max(path1.get(i), path1.get(i+1))));
          if (idx == null) valid = false; else mask[idx] = true;
        }
        for (int i = 0; i < path2.size() - 1 && valid; i++) {
          Integer idx = edgeIdx.get(MolGraph.bondKey(Math.min(path2.get(i), path2.get(i+1)),
                                            Math.max(path2.get(i), path2.get(i+1))));
          if (idx == null) valid = false; else mask[idx] = true;
        }
        if (!valid) continue;

        // Deduplicate by edge mask
        BitSet maskKey = new BitSet(totalEdges);
        for (int bit = 0; bit < totalEdges; bit++) if (mask[bit]) maskKey.set(bit);
        search.work += totalEdges;
        if (!seenMasks.add(maskKey)) continue;

        int ringSize = path1.size() + path2.size() - 1;
        if (ringSize < 3) continue;

        candidateMasks.add(mask);
        candidateSizes.add(ringSize);

        // Build atom list: path1(v→eU) + eV + reverse(path2 interior, eV→...→v excluded)
        // path1 = [v, ..., eU], path2 = [v, ..., eV]
        // Cycle: v → ... → eU → eV → ... → v (backwards through path2)
        List<Integer> atoms = new ArrayList<>(path1); // [v, ..., eU]
        // Append path2 in reverse, skipping the first element (v) since path1 starts at v
        for (int i = path2.size() - 1; i >= 1; i--) atoms.add(path2.get(i));
        candidateAtoms.add(atoms.stream().mapToInt(Integer::intValue).toArray());
      }
    }

    // Sort candidates by size (smallest first)
    Integer[] order = new Integer[candidateMasks.size()];
    for (int i = 0; i < order.length; i++) order[i] = i;
    Arrays.sort(order, Comparator.comparingInt(candidateSizes::get));

    // Interleaved ring basis and relevant-cycle selection, in two reductions:
    // Phase 1: Check relevance (independence from STRICTLY SHORTER cycles only)
    // Phase 2: Check SSSR inclusion (independence from ALL basis cycles)
    List<int[]> sssr = new ArrayList<>();
    List<int[]> legacyRelevant = new ArrayList<>();
    List<boolean[]> basisMatrix = new ArrayList<>();
    List<Integer> basisSizes = new ArrayList<>();

    // O(E) leading-bit lookup for row-echelon reduction
    int[] basisByLeadingBit = new int[totalEdges];
    Arrays.fill(basisByLeadingBit, -1);

    for (int idx : order) {
      boolean[] mask = candidateMasks.get(idx);
      int cSize = candidateSizes.get(idx);

      // Phase 1: Check relevance — reduce only against STRICTLY SHORTER basis cycles
      boolean[] phase1Mask = mask.clone();
      for (int b = 0; b < totalEdges; b++) {
        if (phase1Mask[b]) {
          int basisIdx = basisByLeadingBit[b];
          if (basisIdx != -1 && basisSizes.get(basisIdx) < cSize) {
            for (int k = b; k < totalEdges; k++) {
              phase1Mask[k] ^= basisMatrix.get(basisIdx)[k];
            }
          }
        }
      }

      boolean relevant = false;
      for (boolean bit : phase1Mask) if (bit) { relevant = true; break; }

      if (relevant) {
        legacyRelevant.add(candidateAtoms.get(idx));
        // Phase 2: Check SSSR inclusion — reduce against ALL basis cycles
        boolean[] phase2Mask = phase1Mask.clone();
        for (int b = 0; b < totalEdges; b++) {
          if (phase2Mask[b]) {
            int basisIdx = basisByLeadingBit[b];
            if (basisIdx != -1) {
              for (int k = b; k < totalEdges; k++) {
                phase2Mask[k] ^= basisMatrix.get(basisIdx)[k];
              }
            }
          }
        }

        boolean independent = false;
        int newLeadingBit = -1;
        for (int b = 0; b < totalEdges; b++) {
          if (phase2Mask[b]) { independent = true; newLeadingBit = b; break; }
        }

        if (independent && sssr.size() < cycleRank) {
          sssr.add(candidateAtoms.get(idx));
          basisMatrix.add(phase2Mask); // Add the FULLY REDUCED mask
          basisSizes.add(cSize);
          basisByLeadingBit[newLeadingBit] = basisMatrix.size() - 1;
        }
      }
    }

    return RingPerception.selected(n, neighbors, systems, cycleRank, sssr.toArray(new int[0][]),
                    legacyRelevant.toArray(new int[0][]), enumerateRelevant, false);
  }

  /** Every atom-bond pair of the pairwise search, with its paths and its charge. */
  static final class Replay {
    final List<int[]> pairs = new ArrayList<>();
    final List<int[]> firstPaths = new ArrayList<>(), secondPaths = new ArrayList<>();
    final List<Boolean> valid = new ArrayList<>();
    final List<Long> searchCharge = new ArrayList<>(), charge = new ArrayList<>();
    /** The work counted before the last pair, and the charge of the last pair. */
    long beforeLast, last;

    /** Whether the pairwise search falls back at this work limit. */
    boolean fallsBack(long workLimit) {
      return !pairs.isEmpty() && beforeLast > workLimit;
    }
  }

  /** The pairs of a graph with ring systems, in the order of the pairwise search. */
  static Replay replay(int n, int[][] neighbors) {
    Replay replay = new Replay();
    RingPerception.Systems systems = RingPerception.systems(n, neighbors);
    if (systems == null || RingPerception.cycleRank(n, neighbors) <= 0) return replay;
    List<int[]> edgeList = new ArrayList<>();
    for (int i = 0; i < n; i++)
      for (int j : neighbors[i])
        if (j > i) edgeList.add(new int[] {i, j});
    int totalEdges = edgeList.size();
    int count = systems.atoms.length;
    List<List<Integer>> systemEdges = new ArrayList<>();
    for (int s = 0; s < count; s++) systemEdges.add(new ArrayList<>());
    for (int edge = 0; edge < totalEdges; edge++) {
      int s = systems.ofAtom[edgeList.get(edge)[0]];
      if (s >= 0 && systems.ofAtom[edgeList.get(edge)[1]] == s) systemEdges.get(s).add(edge);
    }
    int[] lowestAtom = new int[count];
    Arrays.fill(lowestAtom, -1);
    for (int atom = n - 1; atom >= 0; atom--) if (systems.ofAtom[atom] >= 0) lowestAtom[systems.ofAtom[atom]] = atom;
    PathSearch search = new PathSearch(n, neighbors, systems);
    for (int v = 0; v < n; v++) {
      int s = systems.ofAtom[v];
      if (s < 0) continue;
      List<Integer> pairs = !systems.simpleRing[s] ? systemEdges.get(s)
          : v == lowestAtom[s] && !systemEdges.get(s).isEmpty() ? Collections.singletonList(systemEdges.get(s).get(0))
          : Collections.<Integer>emptyList();
      for (int pair : pairs) {
        int[] edge = edgeList.get(pair);
        List<Integer> path1 = new ArrayList<>(), path2 = new ArrayList<>();
        search.work = 0;
        search.paths(v, edge[0], edge[1], edge[0], edge[1], path1, path2);
        long cost = search.work;
        boolean disjoint = !path1.isEmpty() && !path2.isEmpty();
        if (disjoint) {
          cost += path1.size() + path2.size();
          Set<Integer> interior = new HashSet<>(path1.subList(1, path1.size()));
          for (int i = 1; i < path2.size() && disjoint; i++) disjoint = !interior.contains(path2.get(i));
          if (disjoint) cost += totalEdges;
        }
        replay.pairs.add(new int[] {v, pair});
        replay.firstPaths.add(path1.stream().mapToInt(Integer::intValue).toArray());
        replay.secondPaths.add(path2.stream().mapToInt(Integer::intValue).toArray());
        replay.valid.add(disjoint);
        replay.searchCharge.add(search.work);
        replay.charge.add(cost);
      }
    }
    int pairs = replay.pairs.size();
    for (int k = 0; k + 1 < pairs; k++) replay.beforeLast += replay.charge.get(k);
    if (pairs > 0) replay.last = replay.charge.get(pairs - 1);
    return replay;
  }

  /**
   * Breadth-first shortest paths avoiding one bond, with reused buffers and a running count of
   * the atoms and neighbours read. With ring systems known, the search stays in the ring
   * system of the target: its atoms are reached in the same order and with the same parents.
   */
  static final class PathSearch {
    final int[][] neighbors;
    final RingPerception.Systems systems;
    final int[] parent, visited, queue;
    int stamp;
    long work;

    PathSearch(int n, int[][] neighbors, RingPerception.Systems systems) {
      this.neighbors = neighbors;
      this.systems = systems;
      parent = new int[n];
      visited = new int[n];
      queue = new int[n];
    }

    /** Shortest paths to both targets, avoiding the bond (avoidU, avoidV). */
    void paths(int start, int firstTarget, int secondTarget, int avoidU, int avoidV,
               List<Integer> firstPath, List<Integer> secondPath) {
      boolean firstReached = start == firstTarget, secondReached = start == secondTarget;
      if (firstReached && secondReached) {
        firstPath.add(start);
        secondPath.add(start);
        return;
      }
      if (++stamp == 0) { Arrays.fill(visited, 0); stamp = 1; }
      int system = systems == null ? -1 : systems.ofAtom[firstTarget];
      visited[start] = stamp;
      parent[start] = -1;
      int head = 0, tail = 0;
      queue[tail++] = start;
      while (head < tail && !(firstReached && secondReached)) {
        int u = queue[head++];
        work += 1 + neighbors[u].length;
        for (int v : neighbors[u]) {
          if ((u == avoidU && v == avoidV) || (u == avoidV && v == avoidU)) continue;
          if (visited[v] == stamp || (systems != null && systems.ofAtom[v] != system)) continue;
          visited[v] = stamp;
          parent[v] = u;
          if (v == firstTarget) firstReached = true;
          if (v == secondTarget) secondReached = true;
          if (firstReached && secondReached) break;
          queue[tail++] = v;
        }
      }
      if (!firstReached) return;
      for (int c = firstTarget; c != -1; c = parent[c]) firstPath.add(c);
      Collections.reverse(firstPath);
      if (!secondReached) return;
      for (int c = secondTarget; c != -1; c = parent[c]) secondPath.add(c);
      Collections.reverse(secondPath);
    }
  }
}
