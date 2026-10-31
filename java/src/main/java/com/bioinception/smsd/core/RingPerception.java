/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
package com.bioinception.smsd.core;

import java.math.BigInteger;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.Collections;
import java.util.Comparator;
import java.util.Deque;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;

/**
 * Ring perception for {@link MolGraph}: the smallest set of smallest rings and the relevant
 * cycles, ring families, ring systems and ring orders for drawing. The searches mirror the
 * native ring basis and relevant-cycle enumeration, stay bounded on very large ring systems and
 * are iterative, so no search depth depends on the molecule size.
 */
final class RingPerception {
  /** Work allowance of one relevant-cycle enumeration, as in the native core. */
  static final long RELEVANT_CYCLE_WORK_LIMIT = 1L << 26;

  private RingPerception() {}

  /**
   * Ring systems found from the bridges in linear time. Atoms joined by bonds that lie on a
   * cycle share a system; an atom whose bonds are all bridges has none (-1). Every cycle lies
   * inside one system, a bond is a bridge exactly when its atoms are not in the same system,
   * and a system with as many bonds as atoms is one simple ring.
   */
  static final class Systems {
    final int[] ofAtom;
    final int[] atoms;
    final boolean[] simpleRing;

    Systems(int[] ofAtom, int[] atoms, boolean[] simpleRing) {
      this.ofAtom = ofAtom;
      this.atoms = atoms;
      this.simpleRing = simpleRing;
    }
  }

  /**
   * Marks the atoms and bonds of a container that lie on a ring, as CDK's
   * {@code Cycles.markRingAtomsAndBonds} does, with an explicit stack instead of one call per
   * atom of a path, so that long chains and macrocycles need no deep call stack. Returns the
   * number of ring-closing bonds found.
   */
  static int markRingAtomsAndBonds(IAtomContainer mol) {
    int n = mol.getAtomCount();
    int[] depth = new int[n], low = new int[n];
    IAtom[] stackAtom = new IAtom[n];
    IBond[] stackBond = new IBond[n];
    @SuppressWarnings("unchecked")
    Iterator<IBond>[] stackBonds = new Iterator[n];
    int closures = 0, remaining = n;
    for (IAtom start : mol.atoms()) {
      if (depth[start.getIndex()] != 0) continue;
      int top = 0;
      stackAtom[0] = start;
      stackBond[0] = null;
      stackBonds[0] = start.bonds().iterator();
      depth[start.getIndex()] = 1;
      low[start.getIndex()] = 2;
      remaining--;
      while (top >= 0) {
        IAtom atom = stackAtom[top];
        int index = atom.getIndex();
        if (stackBonds[top].hasNext()) {
          IBond bond = stackBonds[top].next();
          if (bond == stackBond[top]) continue;
          IAtom other = bond.getOther(atom);
          int reached = depth[other.getIndex()];
          if (reached == 0) {
            top++;
            stackAtom[top] = other;
            stackBond[top] = bond;
            stackBonds[top] = other.bonds().iterator();
            depth[other.getIndex()] = depth[index] + 1;
            low[other.getIndex()] = depth[index] + 2;
            remaining--;
          } else if (reached < depth[index]) {
            closures++;
            bond.setIsInRing(true);
            low[index] = Math.min(low[index], reached);
          }
          continue;
        }
        atom.setIsInRing(low[index] <= depth[index]);
        IBond parentBond = stackBond[top];
        stackBonds[top] = null;
        top--;
        if (top >= 0) {
          int parent = stackAtom[top].getIndex();
          parentBond.setIsInRing(low[index] <= depth[parent]);
          low[parent] = Math.min(low[parent], low[index]);
        }
      }
      if (remaining == 0) break;
    }
    return closures;
  }

  /** The ring systems, or null when an atom lists itself or a neighbour twice. */
  static Systems systems(int n, int[][] neighbors) {
    int[] listedBy = new int[n];
    Arrays.fill(listedBy, -1);
    for (int atom = 0; atom < n; atom++)
      for (int neighbor : neighbors[atom]) {
        if (neighbor == atom || listedBy[neighbor] == atom) return null;
        listedBy[neighbor] = atom;
      }
    // Iterative depth-first search: a tree bond is a bridge when no back bond from the
    // child's subtree reaches its parent or above.
    int[] order = new int[n], low = new int[n], treeParent = new int[n];
    boolean[] bridgeToParent = new boolean[n];
    Arrays.fill(order, -1);
    Arrays.fill(treeParent, -1);
    int[] stackAtom = new int[n], stackNext = new int[n];
    int counter = 0;
    for (int root = 0; root < n; root++) {
      if (order[root] >= 0) continue;
      order[root] = low[root] = counter++;
      int depth = 0;
      stackAtom[0] = root;
      stackNext[0] = 0;
      while (depth >= 0) {
        int atom = stackAtom[depth];
        if (stackNext[depth] < neighbors[atom].length) {
          int neighbor = neighbors[atom][stackNext[depth]++];
          if (neighbor == treeParent[atom]) continue;
          if (order[neighbor] >= 0) {
            low[atom] = Math.min(low[atom], order[neighbor]);
            continue;
          }
          treeParent[neighbor] = atom;
          order[neighbor] = low[neighbor] = counter++;
          depth++;
          stackAtom[depth] = neighbor;
          stackNext[depth] = 0;
          continue;
        }
        depth--;
        int parent = treeParent[atom];
        if (parent >= 0) {
          low[parent] = Math.min(low[parent], low[atom]);
          bridgeToParent[atom] = low[atom] > order[parent];
        }
      }
    }
    int[] ofAtom = new int[n];
    Arrays.fill(ofAtom, -1);
    List<Integer> atoms = new ArrayList<>();
    List<Boolean> simple = new ArrayList<>();
    int[] queue = new int[n];
    for (int start = 0; start < n; start++) {
      if (ofAtom[start] >= 0) continue;
      boolean cyclic = false;
      for (int neighbor : neighbors[start])
        if (!bridge(start, neighbor, treeParent, bridgeToParent)) { cyclic = true; break; }
      if (!cyclic) continue;
      int system = atoms.size();
      long count = 0, bondEnds = 0;
      int head = 0, tail = 0;
      queue[tail++] = start;
      ofAtom[start] = system;
      while (head < tail) {
        int atom = queue[head++];
        count++;
        for (int neighbor : neighbors[atom]) {
          if (bridge(atom, neighbor, treeParent, bridgeToParent)) continue;
          bondEnds++;
          if (ofAtom[neighbor] < 0) {
            ofAtom[neighbor] = system;
            queue[tail++] = neighbor;
          }
        }
      }
      atoms.add((int) count);
      simple.add(bondEnds / 2 == count);
    }
    int[] sizes = new int[atoms.size()];
    boolean[] simpleRing = new boolean[atoms.size()];
    for (int system = 0; system < sizes.length; system++) {
      sizes[system] = atoms.get(system);
      simpleRing[system] = simple.get(system);
    }
    return new Systems(ofAtom, sizes, simpleRing);
  }

  private static boolean bridge(int first, int second, int[] treeParent, boolean[] bridgeToParent) {
    return (treeParent[second] == first && bridgeToParent[second])
        || (treeParent[first] == second && bridgeToParent[first]);
  }

  /** Edge-set identity of a cycle. */
  private static final class EdgeSet {
    final long[] keys;
    final int hash;

    EdgeSet(int[] cycle) {
      keys = new long[cycle.length];
      for (int i = 0; i < cycle.length; i++)
        keys[i] = MolGraph.bondKey(cycle[i], cycle[(i + 1) % cycle.length]);
      Arrays.sort(keys);
      hash = Arrays.hashCode(keys);
    }

    @Override public boolean equals(Object other) {
      return other instanceof EdgeSet && Arrays.equals(keys, ((EdgeSet) other).keys);
    }

    @Override public int hashCode() { return hash; }
  }

  /**
   * The ring basis and relevant cycles of a graph whose exhaustive candidate search would
   * exceed its work bound: {@code {basis, relevant}}. Candidates are the relevant cycles of
   * a bounded enumeration and the rings of a greedy basis, selected by the same two-phase
   * rule as the exhaustive search.
   */
  static int[][][] boundedRings(int n, int[][] neighbors, Systems systems, int cycleRank) {
    int[][] basis = greedyBasis(n, neighbors, systems, cycleRank);
    List<int[]> candidates = relevantCycles(n, neighbors, systems, basis, cycleRank, new boolean[1], false);
    Set<EdgeSet> seen = new HashSet<>();
    for (int[] cycle : candidates) seen.add(new EdgeSet(cycle));
    for (int[] cycle : basis) if (seen.add(new EdgeSet(cycle))) candidates.add(cycle);
    return select(n, neighbors, candidates, cycleRank);
  }

  /** Index of every bond, by lower atom and then neighbour order. */
  private static Map<Long, Integer> edgeIndex(int n, int[][] neighbors) {
    Map<Long, Integer> index = new HashMap<>();
    for (int i = 0; i < n; i++)
      for (int j : neighbors[i])
        if (j > i) index.putIfAbsent(MolGraph.bondKey(i, j), index.size());
    return index;
  }

  /**
   * A cycle basis chosen greedily from the shortest cycle through each bond, smallest first,
   * keeping cycles that add a bond. The search for each bond stays inside its ring system; a
   * system that is one simple ring yields that ring.
   */
  private static int[][] greedyBasis(int n, int[][] neighbors, Systems systems, int cycleRank) {
    int[] parent = new int[n], reached = new int[n], queue = new int[n];
    int[] stamp = {0};
    List<int[]> candidates = new ArrayList<>();
    for (int start = 0; start < n; start++)
      for (int nb : neighbors[start]) {
        if (nb <= start) continue;
        int size = shortestCycle(start, nb, null, neighbors, systems, parent, reached, queue, stamp);
        if (size >= 3) candidates.add(new int[] {size, start, nb});
      }
    candidates.sort((a, b) -> Integer.compare(a[0], b[0]));
    List<int[]> basis = new ArrayList<>();
    Set<Long> allBasis = new HashSet<>();
    List<Integer> cycle = new ArrayList<>();
    for (int[] candidate : candidates) {
      if (basis.size() >= cycleRank) break;
      int system = systems == null ? -1 : systems.ofAtom[candidate[1]];
      if (system >= 0 && systems.simpleRing[system] && !basis.isEmpty()
          && allBasis.contains(MolGraph.bondKey(candidate[1], candidate[2]))) continue;
      shortestCycle(candidate[1], candidate[2], cycle, neighbors, systems, parent, reached, queue, stamp);
      boolean hasNewEdge = false;
      long[] edges = new long[cycle.size()];
      for (int i = 0; i < cycle.size(); i++) {
        edges[i] = MolGraph.bondKey(cycle.get(i), cycle.get((i + 1) % cycle.size()));
        if (!allBasis.contains(edges[i])) hasNewEdge = true;
      }
      if (hasNewEdge || basis.isEmpty()) {
        int[] ring = new int[cycle.size()];
        for (int i = 0; i < ring.length; i++) ring[i] = cycle.get(i);
        basis.add(ring);
        for (long edge : edges) allBasis.add(edge);
      }
    }
    return basis.toArray(new int[0][]);
  }

  /** Size of the shortest cycle through bond (start, nb), and its atoms from nb when path is given. */
  private static int shortestCycle(int start, int nb, List<Integer> path, int[][] neighbors, Systems systems,
                                   int[] parent, int[] reached, int[] queue, int[] stamp) {
    int system = systems == null ? -1 : systems.ofAtom[start];
    if (systems != null) {
      if (system < 0 || systems.ofAtom[nb] != system) return 0;
      if (path == null && systems.simpleRing[system]) return systems.atoms[system];
    }
    if (++stamp[0] == 0) { Arrays.fill(reached, 0); stamp[0] = 1; }
    reached[start] = stamp[0];
    parent[start] = -1;
    int head = 0, tail = 0;
    queue[tail++] = start;
    boolean found = false;
    while (head < tail && !found) {
      int u = queue[head++];
      for (int v : neighbors[u]) {
        if ((u == start && v == nb) || (v == start && u == nb)) continue;
        if (reached[v] == stamp[0]) continue;
        if (systems != null && systems.ofAtom[v] != system) continue;
        reached[v] = stamp[0];
        parent[v] = u;
        if (v == nb) { found = true; break; }
        queue[tail++] = v;
      }
    }
    if (!found) return 0;
    if (path != null) path.clear();
    int size = 0, cur = nb;
    while (cur != start && cur != -1) {
      if (path != null) path.add(cur);
      size++;
      cur = parent[cur];
    }
    if (cur != start) return 0;
    if (path != null) path.add(start);
    return size + 1;
  }

  /**
   * Relevant cycles: cycles of at most the longest basis length that are independent of
   * strictly shorter cycles, enumerated from shortest paths out of every root (at most 50 per
   * atom) inside its ring system, with the work bounded by RELEVANT_CYCLE_WORK_LIMIT.
   */
  private static List<int[]> relevantCycles(int n, int[][] neighbors, Systems systems, int[][] basis,
                                           int cycleRank, boolean[] truncated, boolean boundedResources) {
    truncated[0] = basis.length != cycleRank;
    List<int[]> allCycles = new ArrayList<>();
    if (basis.length == 0) return allCycles;
    int maxLen = boundedResources && basis.length != cycleRank ? n : 0;
    for (int[] cycle : basis) maxLen = Math.max(maxLen, cycle.length);
    long[] work = {0};
    Set<EdgeSet> unique = new HashSet<>();
    for (int[] cycle : basis) {
      if (boundedResources && !chargeRelevant(work, 128 + 16L * cycle.length, truncated))
        return new ArrayList<>(Arrays.asList(basis));
      if (unique.add(new EdgeSet(cycle))) allCycles.add(cycle);
    }

    int[] dist = new int[n], queue = new int[n];
    int[][] parents = new int[n][];
    int[] parentCount = new int[n];
    int[][][] pathsTo = new int[n][][];
    int[] pathCount = new int[n];
    boolean[] listed = new boolean[n];
    int[] mark = new int[n];
    int stamp = 0;
    int[] pathScratch = new int[n + 1];
    int[] stackNode = new int[n + 1], stackNext = new int[n + 1];
    for (int root = 0; root < n && work[0] < RELEVANT_CYCLE_WORK_LIMIT; root++) {
      int system = systems == null ? -1 : systems.ofAtom[root];
      if (systems != null && (system < 0 || systems.simpleRing[system])) continue;
      Arrays.fill(dist, -1);
      Arrays.fill(parentCount, 0);
      Arrays.fill(listed, false);
      work[0] += n;
      dist[root] = 0;
      int head = 0, tail = 0;
      queue[tail++] = root;
      while (head < tail && (!boundedResources || work[0] < RELEVANT_CYCLE_WORK_LIMIT)) {
        int u = queue[head++];
        if (dist[u] >= maxLen / 2) continue;
        work[0] += 1 + neighbors[u].length;
        for (int v : neighbors[u]) {
          if (systems != null && systems.ofAtom[v] != system) continue;
          if (dist[v] == -1) {
            dist[v] = dist[u] + 1;
            addParent(parents, parentCount, v, u);
            queue[tail++] = v;
          } else if (dist[v] == dist[u] + 1) {
            addParent(parents, parentCount, v, u);
          }
        }
      }
      // Odd-length cycles from edges; a tree edge whose farther atom has no other parent
      // closes no cycle.
      for (int u = 0; u < n && work[0] < RELEVANT_CYCLE_WORK_LIMIT; u++) {
        if (dist[u] == -1) continue;
        for (int v : neighbors[u]) {
          if (boundedResources && work[0] >= RELEVANT_CYCLE_WORK_LIMIT) break;
          if (v <= u || dist[v] == -1) continue;
          int cycleLen = dist[u] + dist[v] + 1;
          if (cycleLen < 3 || cycleLen > maxLen) continue;
          if (dist[u] != dist[v] && (dist[u] > dist[v]
              ? onlyParent(parents, parentCount, u, v) : onlyParent(parents, parentCount, v, u))) continue;
          for (int atom : new int[] {u, v})
            if (!listed[atom]) {
              pathsTo[atom] = shortestPaths(root, atom, dist, parents, parentCount, pathScratch, stackNode, stackNext, truncated, work, boundedResources);
              pathCount[atom] = pathsTo[atom].length;
              listed[atom] = true;
              work[0] += (long) pathCount[atom] * (dist[atom] + 1);
            }
          for (int pathU = 0; pathU < pathCount[u] && (!boundedResources || work[0] < RELEVANT_CYCLE_WORK_LIMIT); pathU++) {
            int[] pu = pathsTo[u][pathU];
            int onPath = ++stamp;
            for (int x : pu) mark[x] = onPath;
            work[0] += pu.length;
            for (int pathV = 0; pathV < pathCount[v] && (!boundedResources || work[0] < RELEVANT_CYCLE_WORK_LIMIT); pathV++) {
              int[] pv = pathsTo[v][pathV];
              boolean overlap = false;
              int read = 1;
              for (; read < pv.length; read++) if (mark[pv[read]] == onPath) { overlap = true; break; }
              work[0] += read;
              if (overlap) continue;
              int[] cycle = new int[pu.length + pv.length - 1];
              System.arraycopy(pu, 0, cycle, 0, pu.length);
              for (int k = pv.length - 1, at = pu.length; k >= 1; k--) cycle[at++] = pv[k];
              long cycleWork = remember(cycle, maxLen, unique, allCycles, work, truncated, boundedResources);
              work[0] += cycleWork;
            }
          }
        }
      }
      // Even-length cycles from two disjoint paths to the same atom.
      for (int u = 0; u < n && work[0] < RELEVANT_CYCLE_WORK_LIMIT; u++) {
        if (dist[u] == -1 || parentCount[u] < 2) continue;
        int cycleLen = 2 * dist[u];
        if (cycleLen < 3 || cycleLen > maxLen) continue;
        if (!listed[u]) {
          pathsTo[u] = shortestPaths(root, u, dist, parents, parentCount, pathScratch, stackNode, stackNext, truncated, work, boundedResources);
          pathCount[u] = pathsTo[u].length;
          listed[u] = true;
          work[0] += (long) pathCount[u] * (dist[u] + 1);
        }
        int[][] paths = pathsTo[u];
        for (int pi = 0; pi < pathCount[u] && (!boundedResources || work[0] < RELEVANT_CYCLE_WORK_LIMIT); pi++) {
          int[] p1 = paths[pi];
          int interior = ++stamp;
          for (int k = 1; k + 1 < p1.length; k++) mark[p1[k]] = interior;
          work[0] += p1.length;
          for (int pj = pi + 1; pj < pathCount[u] && (!boundedResources || work[0] < RELEVANT_CYCLE_WORK_LIMIT); pj++) {
            int[] p2 = paths[pj];
            boolean disjoint = true;
            int read = 1;
            for (; read + 1 < p2.length; read++) if (mark[p2[read]] == interior) { disjoint = false; break; }
            work[0] += read;
            if (!disjoint) continue;
            int[] cycle = new int[p1.length + p2.length - 2];
            System.arraycopy(p1, 0, cycle, 0, p1.length);
            for (int k = p2.length - 2, at = p1.length; k >= 1; k--) cycle[at++] = p2[k];
            long cycleWork = remember(cycle, maxLen, unique, allCycles, work, truncated, boundedResources);
            work[0] += cycleWork;
          }
        }
      }
    }
    if (work[0] >= RELEVANT_CYCLE_WORK_LIMIT) truncated[0] = true;
    int levels = 0;
    for (int count = allCycles.size(); count > 1; count = (count + 1) / 2) levels++;
    if (boundedResources && !chargeRelevant(work, 2L * allCycles.size() * levels, truncated))
      return new ArrayList<>(Arrays.asList(basis));
    allCycles.sort((a, b) -> Integer.compare(a.length, b.length));
    // Keep each cycle independent of strictly shorter cycles; rows keep their lowest bond as pivot.
    Map<Long, Integer> edges = edgeIndex(n, neighbors);
    int words = CommonBondGraphs.words(edges.size());
    List<long[]> rows = new ArrayList<>();
    int[] rowByPivot = new int[edges.size()];
    Arrays.fill(rowByPivot, -1);
    List<int[]> relevant = new ArrayList<>();
    int start = 0;
    while (start < allCycles.size()) {
      int end = start;
      while (end < allCycles.size() && allCycles.get(end).length == allCycles.get(start).length) end++;
      for (int k = start; k < end; k++) {
        int[] cycle = allCycles.get(k);
        if (boundedResources && !chargeRelevant(work, 24 + 8L * words + cycle.length, truncated))
          return new ArrayList<>(Arrays.asList(basis));
        if ((boundedResources ? reduceRelevant(bits(cycle, edges, words), rows, rowByPivot, work, truncated) : reduce(bits(cycle, edges, words), rows, rowByPivot, -1, null)) >= 0)
          relevant.add(cycle);
        if (boundedResources && work[0] >= RELEVANT_CYCLE_WORK_LIMIT) return new ArrayList<>(Arrays.asList(basis));
      }
      for (int k = start; k < end; k++) {
        int[] cycle = allCycles.get(k);
        if (boundedResources && !chargeRelevant(work, 24 + 8L * words + cycle.length, truncated))
          return new ArrayList<>(Arrays.asList(basis));
        long[] vector = bits(cycle, edges, words);
        int pivot = boundedResources ? reduceRelevant(vector, rows, rowByPivot, work, truncated) : reduce(vector, rows, rowByPivot, -1, null);
        if (boundedResources && work[0] >= RELEVANT_CYCLE_WORK_LIMIT) return new ArrayList<>(Arrays.asList(basis));
        if (pivot < 0) continue;
        rowByPivot[pivot] = rows.size();
        rows.add(vector);
      }
      start = end;
    }
    return relevant;
  }

  private static void addParent(int[][] parents, int[] parentCount, int atom, int parent) {
    int[] list = parents[atom];
    if (list == null) parents[atom] = list = new int[4];
    else if (parentCount[atom] == list.length) parents[atom] = list = Arrays.copyOf(list, list.length * 2);
    list[parentCount[atom]++] = parent;
  }

  private static boolean onlyParent(int[][] parents, int[] parentCount, int atom, int parent) {
    for (int k = 0; k < parentCount[atom]; k++) if (parents[atom][k] != parent) return false;
    return true;
  }

  private static boolean chargeRelevant(long[] work, long amount, boolean[] truncated) {
    if (work[0] >= RELEVANT_CYCLE_WORK_LIMIT || amount >= RELEVANT_CYCLE_WORK_LIMIT - work[0]) {
      work[0] = RELEVANT_CYCLE_WORK_LIMIT;
      truncated[0] = true;
      return false;
    }
    work[0] += amount;
    return true;
  }

  private static long remember(int[] cycle, int maxLen, Set<EdgeSet> unique, List<int[]> allCycles,
                               long[] work, boolean[] truncated, boolean boundedResources) {
    if (cycle.length < 3 || cycle.length > maxLen) return 0;
    EdgeSet key = new EdgeSet(cycle);
    if (!unique.contains(key)) {
      if (boundedResources && !chargeRelevant(work, 128 + 16L * cycle.length, truncated)) return 0;
      unique.add(key);
      allCycles.add(cycle);
    }
    return cycle.length;
  }

  private static int reduceRelevant(long[] vector, List<long[]> rows, int[] rowByPivot,
                                    long[] work, boolean[] truncated) {
    for (int word = 0; word < vector.length; word++) {
      if (!chargeRelevant(work, 1, truncated)) return -1;
      while (vector[word] != 0) {
        int pivot = (word << 6) + Long.numberOfTrailingZeros(vector[word]);
        int row = rowByPivot[pivot];
        if (row < 0) return pivot;
        long[] basisRow = rows.get(row);
        if (!chargeRelevant(work, vector.length - word, truncated)) return -1;
        for (int other = word; other < vector.length; other++) vector[other] ^= basisRow[other];
      }
    }
    return -1;
  }

  /**
   * Up to 50 shortest paths from root to target, following parents depth first in list order
   * on an explicit stack; each path runs from root to target.
   */
  private static int[][] shortestPaths(int root, int target, int[] dist, int[][] parents, int[] parentCount,
                                       int[] current, int[] stackNode, int[] stackNext, boolean[] truncated,
                                       long[] work, boolean boundedResources) {
    List<int[]> paths = new ArrayList<>();
    int length = 0;
    current[length++] = target;
    if (target == root) {
      if (boundedResources && !chargeRelevant(work, 28, truncated)) return new int[0][];
      paths.add(new int[] {root});
      return paths.toArray(new int[0][]);
    }
    int depth = 0;
    stackNode[0] = target;
    stackNext[0] = 0;
    while (depth >= 0 && paths.size() < 50) {
      int node = stackNode[depth];
      while (stackNext[depth] < parentCount[node] && dist[parents[node][stackNext[depth]]] != dist[node] - 1)
        stackNext[depth]++;
      if (stackNext[depth] == parentCount[node]) {
        depth--;
        length--;
        continue;
      }
      int next = parents[node][stackNext[depth]++];
      current[length++] = next;
      if (next == root) {
        if (boundedResources && !chargeRelevant(work, 24 + 4L * length, truncated)) break;
        int[] path = new int[length];
        for (int k = 0; k < length; k++) path[k] = current[length - 1 - k];
        paths.add(path);
        length--;
        continue;
      }
      depth++;
      stackNode[depth] = next;
      stackNext[depth] = 0;
    }
    if (paths.size() == 50)
      for (int frame = 0; frame <= depth && !truncated[0]; frame++) {
        int node = stackNode[frame];
        for (int parent = stackNext[frame]; parent < parentCount[node] && !truncated[0]; parent++)
          truncated[0] = dist[parents[node][parent]] == dist[node] - 1;
      }
    return paths.toArray(new int[0][]);
  }

  private static long[] bits(int[] cycle, Map<Long, Integer> edges, int words) {
    long[] vector = new long[words];
    for (int i = 0; i < cycle.length; i++) {
      Integer index = edges.get(MolGraph.bondKey(cycle[i], cycle[(i + 1) % cycle.length]));
      if (index != null) vector[index >>> 6] ^= 1L << (index & 63);
    }
    return vector;
  }

  /**
   * Reduces vector by the rows whose size is below size (all rows when size is -1); returns
   * the lowest bond left without such a row, or -1 when the vector is in their span.
   */
  private static int reduce(long[] vector, List<long[]> rows, int[] rowByPivot, int size, List<Integer> rowSizes) {
    for (int word = 0; word < vector.length; word++) {
      while (vector[word] != 0) {
        int pivot = (word << 6) + Long.numberOfTrailingZeros(vector[word]);
        int row = rowByPivot[pivot];
        if (row < 0 || (size >= 0 && rowSizes.get(row) >= size)) return pivot;
        long[] basisRow = rows.get(row);
        for (int other = word; other < vector.length; other++) vector[other] ^= basisRow[other];
      }
    }
    return -1;
  }

  /**
   * The two-phase selection over candidates in ascending size: a candidate independent of
   * strictly shorter selected cycles is relevant, and a relevant one independent of every
   * selected cycle joins the basis until it has cycleRank cycles.
   */
  private static int[][][] select(int n, int[][] neighbors, List<int[]> candidates, int cycleRank) {
    candidates.sort((a, b) -> Integer.compare(a.length, b.length));
    Map<Long, Integer> edges = edgeIndex(n, neighbors);
    int words = CommonBondGraphs.words(edges.size());
    List<long[]> rows = new ArrayList<>();
    List<Integer> rowSizes = new ArrayList<>();
    int[] rowByPivot = new int[edges.size()];
    Arrays.fill(rowByPivot, -1);
    List<int[]> basis = new ArrayList<>(), relevant = new ArrayList<>();
    for (int[] cycle : candidates) {
      long[] vector = bits(cycle, edges, words);
      if (reduce(vector, rows, rowByPivot, cycle.length, rowSizes) < 0) continue;
      int pivot = reduce(vector, rows, rowByPivot, -1, null);
      if (pivot >= 0 && basis.size() < cycleRank) {
        rowByPivot[pivot] = rows.size();
        rows.add(vector);
        rowSizes.add(cycle.length);
        basis.add(cycle);
      }
      relevant.add(cycle);
    }
    return new int[][][] {basis.toArray(new int[0][]), relevant.toArray(new int[0][])};
  }

  // ---- Rings of a graph, perceived once and kept by it ----

  /**
   * Work bound of the exhaustive candidate search of {@link #perceive}. The work is that of the
   * search that pairs every atom with every bond of its ring system by one breadth-first search
   * per pair, counted in atoms and neighbours read: one plus the degree of every atom such a
   * search takes from its queue, the atoms of both paths of every pair and the bonds of the graph
   * for every candidate cycle. The candidates themselves are read from one shortest-path tree per
   * atom ({@link RingCandidateSearch}), which finds the same count without spending it. Natural
   * products, large cyclodextrins and fused polycyclic sheets of several hundred atoms stay within
   * the bound. A graph that exceeds it, such as a very large ring system, takes its candidates from
   * a bounded relevant-cycle enumeration and a greedy ring basis, and the same selection rules
   * apply. A full basis found that way is exact only when {@link #certifiedMinimumBasis} proves it
   * minimum within the same bound.
   */
  static final long RING_SEARCH_WORK_LIMIT = 1L << 29;

  /**
   * The smallest set of smallest rings and the relevant cycles (the union of all minimum cycle
   * bases) of a graph, perceived together and never changed.
   */
  static final class Rings {
    final int[][] smallest, relevant, legacyRelevant, numericBasis;
    /** The basis is not proved to be a minimum cycle basis of full cycle rank. */
    final boolean basisTruncated;
    final boolean relevantTruncated;
    /** The original Horton work-bound flag; fusion counts are unproved while it is set. */
    final boolean truncated;
    /**
     * The basis the bounded fallback chose. It is the same array as smallest unless a basis the
     * fallback left unproved was completed; the relevant cycles and tautomer classes then keep
     * reading this one.
     */
    final int[][] baseSmallest;
    /** A basis the fallback left unproved was completed and proved minimum. */
    final boolean basisLifted;

    Rings(int[][] smallest, int[][] relevant) {
      this(smallest, relevant, false);
    }

    Rings(int[][] smallest, int[][] relevant, boolean truncated) {
      this(smallest, relevant, relevant, smallest, truncated, truncated);
    }

    Rings(int[][] smallest, int[][] relevant, int[][] legacyRelevant, int[][] numericBasis,
          boolean basisTruncated, boolean relevantTruncated) {
      this(smallest, relevant, legacyRelevant, numericBasis, basisTruncated, relevantTruncated, basisTruncated);
    }

    Rings(int[][] smallest, int[][] relevant, int[][] legacyRelevant, int[][] numericBasis,
          boolean basisTruncated, boolean relevantTruncated, boolean truncated) {
      this(smallest, relevant, legacyRelevant, numericBasis, basisTruncated, relevantTruncated, truncated,
           smallest, false);
    }

    Rings(int[][] smallest, int[][] relevant, int[][] legacyRelevant, int[][] numericBasis,
          boolean basisTruncated, boolean relevantTruncated, boolean truncated,
          int[][] baseSmallest, boolean basisLifted) {
      this.smallest = smallest;
      this.relevant = relevant != null ? relevant : new int[0][];
      this.legacyRelevant = legacyRelevant != null ? legacyRelevant : new int[0][];
      this.numericBasis = numericBasis;
      this.basisTruncated = basisTruncated;
      this.relevantTruncated = relevantTruncated;
      this.truncated = truncated;
      this.baseSmallest = baseSmallest;
      this.basisLifted = basisLifted;
    }
  }

  /**
   * Whether {@code cycles}, {@code cycleRank} of them, form a minimum cycle basis, decided
   * exactly or not at all, by the rule and work count of the native core. Written over the
   * bonds outside a breadth-first spanning forest, the cycles must be independent; each cycle
   * then has a dual bond set that it meets an odd number of times and every other cycle meets
   * an even number of times. They form a minimum basis exactly when no cycle shorter than a
   * member meets that member's dual set an odd number of times, which a breadth-first search
   * over atoms paired with a parity settles. Work is counted in words reduced and in atoms and
   * neighbours read; the answer is false when the cycles are not a basis or the work exceeds
   * {@code workLimit}.
   */
  static boolean certifiedMinimumBasis(int n, int[][] neighbors, int[][] cycles, int cycleRank, long workLimit) {
    if (cycleRank <= 0 || cycles.length != cycleRank) return false;
    Map<Long, Integer> edges = edgeIndex(n, neighbors);
    int edgeCount = edges.size();
    int[][] slotEdge = new int[n][];
    long bondEnds = 0;
    for (int v = 0; v < n; v++) {
      slotEdge[v] = new int[neighbors[v].length];
      for (int slot = 0; slot < neighbors[v].length; slot++) {
        int w = neighbors[v][slot];
        if (w == v || w < 0 || w >= n) return false;
        slotEdge[v][slot] = edges.get(MolGraph.bondKey(v, w));
      }
      bondEnds += neighbors[v].length;
    }
    if (bondEnds != 2L * edgeCount) return false;

    // Spanning forest in atom and neighbour order; the other bonds are the coordinates.
    boolean[] treeEdge = new boolean[edgeCount], reached = new boolean[n];
    int[] queue = new int[n];
    for (int root = 0; root < n; root++) {
      if (reached[root]) continue;
      reached[root] = true;
      int head = 0, tail = 0;
      queue[tail++] = root;
      while (head < tail) {
        int v = queue[head++];
        for (int slot = 0; slot < neighbors[v].length; slot++) {
          int w = neighbors[v][slot];
          if (reached[w]) continue;
          reached[w] = true;
          treeEdge[slotEdge[v][slot]] = true;
          queue[tail++] = w;
        }
      }
    }
    int[] coordinate = new int[edgeCount];
    int k = 0;
    for (int edge = 0; edge < edgeCount; edge++) coordinate[edge] = treeEdge[edge] ? -1 : k++;
    if (k != cycleRank) return false;
    int[] coordinateEdge = new int[k];
    for (int edge = 0; edge < edgeCount; edge++) if (coordinate[edge] >= 0) coordinateEdge[coordinate[edge]] = edge;
    int words = (k + 63) >>> 6;
    // The reduction keeps two k-by-k bit matrices, within a fixed share of the work bound.
    if ((long) k * 2 * words > (workLimit >>> 7)) return false;

    // transposed[r] has bit c when cycle c uses coordinate bond r; dual starts as the identity.
    long[][] transposed = new long[k][words], dual = new long[k][words];
    for (int c = 0; c < k; c++) {
      int[] cycle = cycles[c];
      if (cycle.length < 3) return false;
      for (int i = 0; i < cycle.length; i++) {
        int a = cycle[i], b = cycle[(i + 1) % cycle.length];
        if (a < 0 || a >= n || b < 0 || b >= n) return false;
        Integer edge = edges.get(MolGraph.bondKey(a, b));
        if (edge == null) return false;
        int r = coordinate[edge];
        if (r >= 0) transposed[r][c >>> 6] ^= 1L << (c & 63);
      }
    }
    for (int r = 0; r < k; r++) dual[r][r >>> 6] |= 1L << (r & 63);
    long work = 0;
    for (int c = 0; c < k; c++) {
      int word = c >>> 6;
      long bit = 1L << (c & 63);
      int pivot = c;
      while (pivot < k && (transposed[pivot][word] & bit) == 0) { pivot++; work++; }
      if (pivot == k) return false;
      long[] swap = transposed[pivot]; transposed[pivot] = transposed[c]; transposed[c] = swap;
      swap = dual[pivot]; dual[pivot] = dual[c]; dual[c] = swap;
      for (int r = 0; r < k; r++) {
        if (r == c || (transposed[r][word] & bit) == 0) continue;
        for (int w = 0; w < words; w++) {
          transposed[r][w] ^= transposed[c][w];
          dual[r][w] ^= dual[c][w];
        }
        work += 2L * words;
      }
      work += k;
      if (work > workLimit) return false;
    }

    // dual[i] meets cycle i once and every other cycle an even number of times.
    boolean[] inDual = new boolean[edgeCount], start = new boolean[n];
    int[] seen = new int[2 * n], depth = new int[2 * n], states = new int[2 * n];
    int stamp = 0;
    for (int i = 0; i < k; i++) {
      int length = cycles[i].length;
      if (length == 3) continue;
      Arrays.fill(inDual, false);
      Arrays.fill(start, false);
      for (int r = 0; r < k; r++)
        if ((dual[i][r >>> 6] & (1L << (r & 63))) != 0) inDual[coordinateEdge[r]] = true;
      // Every such closed walk crosses a dual bond, so it can start at one of their atoms.
      for (int v = 0; v < n; v++)
        for (int slot = 0; slot < slotEdge[v].length; slot++)
          if (inDual[slotEdge[v][slot]]) start[v] = true;
      for (int u = 0; u < n; u++) {
        if (!start[u]) continue;
        if (++stamp == 0) { Arrays.fill(seen, 0); stamp = 1; }
        int head = 0, tail = 0;
        states[tail++] = 2 * u;
        seen[2 * u] = stamp;
        depth[2 * u] = 0;
        while (head < tail) {
          int state = states[head++];
          int v = state >> 1, parity = state & 1;
          work += 1 + neighbors[v].length;
          if (work > workLimit) return false;
          int d = depth[state];
          if (d >= length - 1) continue;
          for (int slot = 0; slot < neighbors[v].length; slot++) {
            int next = 2 * neighbors[v][slot] + (parity ^ (inDual[slotEdge[v][slot]] ? 1 : 0));
            // A closed walk shorter than the member that meets its dual set oddly.
            if (next == 2 * u + 1) return false;
            if (seen[next] == stamp) continue;
            seen[next] = stamp;
            depth[next] = d + 1;
            states[tail++] = next;
          }
        }
      }
    }
    return true;
  }

  static int cycleRank(int n, int[][] neighbors) {
    int edgeCount = 0;
    for (int i = 0; i < n; i++) edgeCount += neighbors[i].length;
    edgeCount /= 2;

    boolean[] visited = new boolean[n];
    int components = 0;
    for (int i = 0; i < n; i++) {
      if (visited[i]) continue;
      components++;
      Deque<Integer> bfsQ = new ArrayDeque<>();
      bfsQ.add(i);
      visited[i] = true;
      while (!bfsQ.isEmpty()) {
        int u = bfsQ.poll();
        for (int v : neighbors[u]) if (!visited[v]) { visited[v] = true; bfsQ.add(v); }
      }
    }
    return edgeCount - n + components;
  }

  /** Fill the relevant set without replacing the already published basis object. */
  static Rings completeRelevant(int n, int[][] neighbors, Rings basis) {
    return completeRelevant(n, neighbors, basis, RingCandidateSearch.RingSearchOptions.DEFAULTS);
  }

  /**
   * As {@link #completeRelevant(int, int[][], Rings)}. After a completed basis the relevant cycles
   * still come from the basis the fallback chose, and stay marked incomplete.
   */
  static Rings completeRelevant(int n, int[][] neighbors, Rings basis, RingCandidateSearch.RingSearchOptions options) {
    Systems systems = systems(n, neighbors);
    int rank = cycleRank(n, neighbors);
    boolean[] cut = new boolean[1];
    List<int[]> relevant = relevantCycles(n, neighbors, systems, basis.baseSmallest, rank, cut, true);
    Rings rings = new Rings(basis.smallest, relevant.toArray(new int[0][]), basis.legacyRelevant,
                            basis.numericBasis, basis.basisTruncated,
                            cut[0] || basis.basisTruncated || basis.basisLifted, basis.truncated,
                            basis.baseSmallest, basis.basisLifted);
    return options != null && options.proveRelevantByFamilies ? proveRelevantComplete(n, neighbors, rings) : rings;
  }

  /**
   * The numeric SMARTS basis that {@link #selected} pairs with a basis: an equally weighted full
   * greedy basis, or the basis itself.
   */
  static int[][] numericBasis(int n, int[][] neighbors, Systems systems, int cycleRank, int[][] smallest) {
    if (smallest.length != cycleRank) return smallest;
    int[][] greedy = greedyBasis(n, neighbors, systems, cycleRank);
    long greedyWeight = 0, minimumWeight = 0;
    for (int[] cycle : greedy) greedyWeight += cycle.length;
    for (int[] cycle : smallest) minimumWeight += cycle.length;
    return greedy.length == cycleRank && greedyWeight == minimumWeight ? greedy : smallest;
  }

  /**
   * The rings with the relevant-cycle truncation flag cleared when the relevant-cycle families
   * prove the listed cycles complete: the basis is complete and proved without a completion, the
   * families are complete, and the listed cycles are distinct members of them, as many as all
   * their members. Otherwise the rings themselves; the list never changes.
   */
  static Rings proveRelevantComplete(int n, int[][] neighbors, Rings rings) {
    if (!rings.relevantTruncated || rings.basisTruncated || rings.basisLifted) return rings;
    RelevantCycleFamilies families = RelevantCycleFamilySearch.compute(n, neighbors);
    if (families.truncated() || !families.relevantCycleCount().equals(BigInteger.valueOf(rings.relevant.length)))
      return rings;
    Set<EdgeSet> seen = new HashSet<>();
    for (int[] cycle : rings.relevant) {
      if (cycle.length < 3 || !seen.add(new EdgeSet(cycle))) return rings;
      if (RelevantCycleFamilySearch.familyOf(n, neighbors, families, cycle) < 0) return rings;
    }
    return new Rings(rings.smallest, rings.relevant, rings.legacyRelevant, rings.numericBasis,
                     rings.basisTruncated, false, rings.truncated, rings.baseSmallest, rings.basisLifted);
  }

  /** The rings of a graph, found together by one candidate search. */
  static Rings perceive(int n, int[][] neighbors) {
    return perceive(n, neighbors, true);
  }

  static Rings perceive(int n, int[][] neighbors, boolean enumerateRelevant) {
    return perceive(n, neighbors, enumerateRelevant, RING_SEARCH_WORK_LIMIT);
  }

  /** As {@link #perceive(int, int[][], boolean)}, with the work bound a test may lower. */
  static Rings perceive(int n, int[][] neighbors, boolean enumerateRelevant, long workLimit) {
    return perceive(n, neighbors, enumerateRelevant, workLimit, RingCandidateSearch.RingSearchOptions.DEFAULTS, null);
  }

  /**
   * As {@link #perceive(int, int[][], boolean, long)}, with internal switches, all off by default,
   * and a record of the search in notes when notes is not null.
   */
  static Rings perceive(int n, int[][] neighbors, boolean enumerateRelevant, long workLimit,
                        RingCandidateSearch.RingSearchOptions options, RingCandidateSearch.RingSearchNotes notes) {
    if (n == 0) return new Rings(new int[0][], new int[0][]);

    int cycleRank = cycleRank(n, neighbors);
    if (cycleRank <= 0) return new Rings(new int[0][], new int[0][]);

    // With ring systems known, the candidates are read from one shortest-path tree per atom; the
    // pairwise search below remains for a graph whose atom lists itself or a neighbour twice.
    Systems systems = systems(n, neighbors);
    if (systems != null) {
      Rings rings = new RingCandidateSearch(n, neighbors, systems, cycleRank, workLimit, options, notes)
          .perceive(enumerateRelevant);
      return options != null && options.proveRelevantByFamilies ? proveRelevantComplete(n, neighbors, rings) : rings;
    }

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
          int[][][] bounded = boundedRings(n, neighbors, systems, cycleRank);
          if (!certifiedMinimumBasis(n, neighbors, bounded[0], cycleRank, workLimit))
            return new Rings(bounded[0], bounded[1], true);
          // The work bound stays recorded: fusion counts still come from the bounded candidates.
          return selected(n, neighbors, systems, cycleRank, bounded[0], bounded[1], enumerateRelevant, true);
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

    return selected(n, neighbors, systems, cycleRank, sssr.toArray(new int[0][]),
                    legacyRelevant.toArray(new int[0][]), enumerateRelevant, false);
  }

  /**
   * The rings of a selected basis: an equally weighted full greedy basis keeps the numeric
   * SMARTS tie-break, and the relevant cycles are enumerated from the basis when asked for.
   */
  static Rings selected(int n, int[][] neighbors, Systems systems, int cycleRank, int[][] smallest,
                        int[][] legacy, boolean enumerateRelevant, boolean searchLimited) {
    boolean basisTruncated = smallest.length != cycleRank;
    boolean truncated = basisTruncated || searchLimited;
    int[][] numericBasis = smallest;
    if (!basisTruncated) {
      int[][] greedy = greedyBasis(n, neighbors, systems, cycleRank);
      long greedyWeight = 0, minimumWeight = 0;
      for (int[] cycle : greedy) greedyWeight += cycle.length;
      for (int[] cycle : smallest) minimumWeight += cycle.length;
      if (greedy.length == cycleRank && greedyWeight == minimumWeight) numericBasis = greedy;
    }
    if (!enumerateRelevant)
      return new Rings(smallest, legacy, legacy, numericBasis, basisTruncated, basisTruncated, truncated);
    boolean[] relevantTruncated = new boolean[1];
    List<int[]> relevant = relevantCycles(n, neighbors, systems, smallest, cycleRank, relevantTruncated, true);
    return new Rings(smallest, relevant.toArray(new int[0][]), legacy, numericBasis,
                     basisTruncated, relevantTruncated[0] || basisTruncated, truncated);
  }

  /**
   * Breadth-first shortest paths avoiding one bond, with reused buffers and a running count of
   * the atoms and neighbours read. With ring systems known, the search stays in the ring
   * system of the target: its atoms are reached in the same order and with the same parents.
   */
  static final class PathSearch {
    final int[][] neighbors;
    final Systems systems;
    final int[] parent, visited, queue;
    int stamp;
    long work;

    PathSearch(int n, int[][] neighbors, Systems systems) {
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

  /**
   * Ring systems as connected ring atoms: atoms flagged as ring atoms and joined by a bond share
   * a system, numbered in atom order, and each system has a hash of its size, aromatic atoms,
   * elements and bond orders.
   */
  static final class AtomSystems {
    /** The system of each atom, or -1 for an atom outside every ring. */
    final int[] ofAtom;
    /** The hash of each system. */
    final long[] signature;

    AtomSystems(int[] ofAtom, long[] signature) {
      this.ofAtom = ofAtom;
      this.signature = signature;
    }
  }

  /** The ring systems of the ring atoms of {@code g}. */
  static AtomSystems atomSystems(MolGraph g) {
    int n = g.n;
    boolean[] ring = g.ring, aromatic = g.aromatic;
    int[] atomicNum = g.atomicNum;
    int[][] neighbors = g.neighbors;
    int[] uf = UnionFind.create(n);
    // Unite ring atoms connected by bonds where both endpoints are ring atoms
    for (int i = 0; i < n; i++) {
      if (!ring[i]) continue;
      for (int j : neighbors[i]) {
        if (j > i && ring[j]) UnionFind.join(uf, i, j);
      }
    }
    // Assign contiguous IDs
    int[] idMap = new int[n];
    Arrays.fill(idMap, -1);
    int count = 0;
    int[] sysId = new int[n];
    Arrays.fill(sysId, -1);
    for (int i = 0; i < n; i++) {
      if (!ring[i]) continue;
      int root = UnionFind.root(uf, i);
      if (idMap[root] < 0) idMap[root] = count++;
      sysId[i] = idMap[root];
    }
    // A 64-bit FNV-1a hash of each ring system
    long[] sig = new long[count];
    // Gather per-system data: size, aromCount, atomicNums, bondOrders
    int[] sysSize = new int[count];
    int[] sysArom = new int[count];
    @SuppressWarnings("unchecked")
    List<Integer>[] sysAtomicNums = new List[count];
    @SuppressWarnings("unchecked")
    List<Integer>[] sysBondOrds = new List[count];
    for (int s = 0; s < count; s++) {
      sysAtomicNums[s] = new ArrayList<>();
      sysBondOrds[s] = new ArrayList<>();
    }
    for (int i = 0; i < n; i++) {
      if (sysId[i] < 0) continue;
      int s = sysId[i];
      sysSize[s]++;
      if (aromatic[i]) sysArom[s]++;
      sysAtomicNums[s].add(atomicNum[i]);
      for (int j : neighbors[i]) {
        if (j > i && ring[j] && sysId[j] == s) {
          sysBondOrds[s].add(g.bondOrder(i, j));
        }
      }
    }
    for (int s = 0; s < count; s++) {
      long h = FingerprintEngine.FNV1A_SEED;
      h ^= sysSize[s]; h *= FingerprintEngine.FNV1A_PRIME;
      h ^= sysArom[s]; h *= FingerprintEngine.FNV1A_PRIME;
      Collections.sort(sysAtomicNums[s]);
      for (int v : sysAtomicNums[s]) { h ^= v; h *= FingerprintEngine.FNV1A_PRIME; }
      Collections.sort(sysBondOrds[s]);
      for (int v : sysBondOrds[s]) { h ^= v; h *= FingerprintEngine.FNV1A_PRIME; }
      sig[s] = h;
    }
    return new AtomSystems(sysId, sig);
  }

  /**
   * The relevant cycles of {@code g} grouped into ring families: two cycles of one size share a
   * family when one is the other with a single bond exchanged, or when every atom of both lies
   * in one atom orbit.
   */
  static int[][][] uniqueRingFamilies(MolGraph g) {
    g.ensureCanonical();
    int n = g.n;
    int[][] neighbors = g.neighbors;
    int[] orbit = g.orbit;
    int[][] rc = g.computeRelevantCycles();
    if (rc.length == 0) return new int[0][][];

    // Build edge index + O(1) inverse lookup
    Map<Long, Integer> edgeIndex = edgeIndex(n, neighbors);
    int numEdges = edgeIndex.size();
    int longWords = CommonBondGraphs.words(numEdges);

    // O(1) inverse lookup: edge index → (u, v) pair
    int[][] indexToEdge = new int[numEdges][2];
    for (Map.Entry<Long, Integer> e : edgeIndex.entrySet()) {
      long k = e.getKey();
      int idx = e.getValue();
      indexToEdge[idx] = new int[] { (int)(k >>> 32), (int)(k & 0xFFFFFFFFL) };
    }

    // Build edge-incidence vectors for each relevant cycle
    long[][] vectors = new long[rc.length][longWords];
    for (int ci = 0; ci < rc.length; ci++) {
      int[] c = rc[ci];
      for (int i = 0; i < c.length; i++) {
        int a = c[i], b = c[(i + 1) % c.length];
        long key = MolGraph.bondKey(Math.min(a, b), Math.max(a, b));
        Integer idx = edgeIndex.get(key);
        if (idx != null) vectors[ci][idx >>> 6] ^= 1L << (idx & 63);
      }
    }

    int nrc = rc.length;
    int[] family = UnionFind.create(nrc);

    // Generate orbit-based signatures for each cycle (lexicographically minimal)
    String[] sigs = new String[nrc];
    for (int i = 0; i < nrc; i++) sigs[i] = orbitSignature(rc[i], orbit);

    // Group cycles into URFs using two criteria (in order):
    //   1. Orbit-signature match: only when every atom in both cycles shares the
    //      same orbit value (vertex-transitive ring systems, e.g. cubane where
    //      all atoms have orbit=0). Without this guard, fused systems like
    //      naphthalene are incorrectly merged because both 6-rings happen to
    //      produce the same orbit sequence after rotation/reflection.
    //   2. Single-edge exchange: symmetric difference == 2 edges
    //      sharing a common vertex.
    for (int i = 0; i < nrc; i++) {
      for (int j = i + 1; j < nrc; j++) {
        if (rc[i].length != rc[j].length) continue;

        // 1. Orbit-signature match — guarded: only apply when all atoms in
        //    both cycles belong to the same orbit (vertex-transitive).
        if (sigs[i].equals(sigs[j])) {
          int orb0 = orbit[rc[i][0]];
          boolean allSame = true;
          for (int a : rc[i]) { if (orbit[a] != orb0) { allSame = false; break; } }
          if (allSame) for (int a : rc[j]) { if (orbit[a] != orb0) { allSame = false; break; } }
          if (allSame) {
            UnionFind.join(family, i, j);
            continue;
          }
        }

        // 2. Single-edge exchange (handles non-symmetric but interchangeable rings)
        long[] xor = new long[longWords];
        int popcount = 0;
        for (int w = 0; w < longWords; w++) {
          xor[w] = vectors[i][w] ^ vectors[j][w];
          popcount += Long.bitCount(xor[w]);
        }
        if (popcount == 2) {
          int[] diffEdges = new int[2];
          int di = 0;
          for (int w = 0; w < longWords && di < 2; w++) {
            long bits = xor[w];
            while (bits != 0 && di < 2) {
              int bit = Long.numberOfTrailingZeros(bits);
              diffEdges[di++] = (w << 6) + bit;
              bits &= bits - 1;
            }
          }
          int[] e1 = indexToEdge[diffEdges[0]]; // O(1) instead of O(E)
          int[] e2 = indexToEdge[diffEdges[1]];
          if (e1 != null && e2 != null && sharesVertex(e1, e2)) {
            UnionFind.join(family, i, j);
          }
        }
      }
    }

    // Group by family
    Map<Integer, List<Integer>> groups = new HashMap<>();
    for (int i = 0; i < nrc; i++) {
      groups.computeIfAbsent(UnionFind.root(family, i), k -> new ArrayList<>()).add(i);
    }

    int[][][] result = new int[groups.size()][][];
    int gi = 0;
    for (List<Integer> members : groups.values()) {
      result[gi] = new int[members.size()][];
      for (int m = 0; m < members.size(); m++) {
        result[gi][m] = rc[members.get(m)];
      }
      gi++;
    }
    return result;
  }

  private static boolean sharesVertex(int[] e1, int[] e2) {
    return e1[0] == e2[0] || e1[0] == e2[1] || e1[1] == e2[0] || e1[1] == e2[1];
  }

  /** The lexicographically least sequence of the atom orbits round a cycle, in either direction. */
  private static String orbitSignature(int[] c, int[] orbit) {
    int len = c.length;
    int[] bestSeq = null;
    for (int dir = 0; dir <= 1; dir++) {
      for (int start = 0; start < len; start++) {
        int[] seq = new int[len];
        for (int i = 0; i < len; i++) {
          int idx = dir == 0 ? (start + i) % len : (start - i + len) % len;
          seq[i] = orbit[c[idx]];
        }
        if (bestSeq == null || CanonicalLabels.compareIntArrays(seq, bestSeq) < 0) bestSeq = seq;
      }
    }
    return Arrays.toString(bestSeq);
  }

  /** The smallest set of smallest rings of {@code g}, sorted by ascending size. */
  static List<List<Integer>> smallestRingsBySize(MolGraph g) {
    int[][] raw = g.computeRings();
    List<List<Integer>> rings = new ArrayList<>(raw.length);
    for (int[] r : raw) {
      List<Integer> ring = new ArrayList<>(r.length);
      for (int a : r) ring.add(a);
      rings.add(ring);
    }
    rings.sort(Comparator.comparingInt(List::size));
    return rings;
  }

  /**
   * The rings of {@link #smallestRingsBySize} ordered for drawing: the fused system with most
   * atoms first, and within a system the smallest ring first and then its fused neighbours,
   * breadth-first and smallest first.
   */
  static List<List<Integer>> layoutOrder(MolGraph g) {
    List<List<Integer>> rings = smallestRingsBySize(g);
    if (rings.size() <= 1) return rings;

    int nr = rings.size();
    int[][] cycles = new int[nr][];
    for (int i = 0; i < nr; i++) cycles[i] = rings.get(i).stream().mapToInt(Integer::intValue).toArray();
    List<List<Integer>> adj = fusedNeighbours(cycles);
    int[] comp = fusedSystemOf(adj);
    int nComp = systemCount(comp);

    // Group by component, compute system size (unique atom count)
    List<List<Integer>> systems = new ArrayList<>(nComp);
    for (int i = 0; i < nComp; i++) systems.add(new ArrayList<>());
    for (int i = 0; i < nr; i++) systems.get(comp[i]).add(i);

    int[] sysSize = new int[nComp];
    for (int c = 0; c < nComp; c++) {
      Set<Integer> atoms = new HashSet<>();
      for (int ri : systems.get(c))
        atoms.addAll(rings.get(ri));
      sysSize[c] = atoms.size();
    }

    // Sort systems: largest first
    Integer[] sysOrder = new Integer[nComp];
    for (int i = 0; i < nComp; i++) sysOrder[i] = i;
    Arrays.sort(sysOrder, (a, b) -> Integer.compare(sysSize[b], sysSize[a]));

    List<List<Integer>> result = new ArrayList<>(nr);
    for (int ci : sysOrder) {
      List<Integer> sysRings = systems.get(ci);
      if (sysRings.size() == 1) {
        result.add(rings.get(sysRings.get(0)));
        continue;
      }

      sysRings.sort(Comparator.comparingInt(a -> rings.get(a).size()));

      boolean[] visited = new boolean[nr];
      Deque<Integer> bfs = new ArrayDeque<>();
      bfs.add(sysRings.get(0));
      visited[sysRings.get(0)] = true;
      while (!bfs.isEmpty()) {
        int u = bfs.poll();
        result.add(rings.get(u));
        List<Integer> nbrs = new ArrayList<>();
        for (int v : adj.get(u)) {
          if (!visited[v] && comp[v] == ci) nbrs.add(v);
        }
        nbrs.sort(Comparator.comparingInt(a -> rings.get(a).size()));
        for (int v : nbrs) {
          if (!visited[v]) { visited[v] = true; bfs.add(v); }
        }
      }
    }
    return result;
  }

  /** The atoms of each fused ring system of {@code g}: rings that share a bond share a system. */
  static List<Set<Integer>> fusedRingAtomSets(MolGraph g) {
    int[][] raw = g.computeRings();
    if (raw.length == 0) return Collections.emptyList();

    int nr = raw.length;
    int[] comp = fusedSystemOf(fusedNeighbours(raw));
    int nComp = systemCount(comp);

    List<Set<Integer>> systems = new ArrayList<>(nComp);
    for (int i = 0; i < nComp; i++) systems.add(new HashSet<>());
    for (int i = 0; i < nr; i++)
      for (int a : raw[i]) systems.get(comp[i]).add(a);
    return systems;
  }

  /** For each ring, the rings that share a bond with it, in ascending order. */
  private static List<List<Integer>> fusedNeighbours(int[][] rings) {
    int nr = rings.length;
    List<Set<Long>> edgeSets = new ArrayList<>(nr);
    for (int[] r : rings) {
      Set<Long> es = new HashSet<>();
      for (int j = 0; j < r.length; j++) es.add(MolGraph.bondKey(r[j], r[(j + 1) % r.length]));
      edgeSets.add(es);
    }

    List<List<Integer>> adj = new ArrayList<>(nr);
    for (int i = 0; i < nr; i++) adj.add(new ArrayList<>());
    for (int i = 0; i < nr; i++) {
      for (int j = i + 1; j < nr; j++) {
        for (Long e : edgeSets.get(i)) {
          if (edgeSets.get(j).contains(e)) {
            adj.get(i).add(j); adj.get(j).add(i); break;
          }
        }
      }
    }
    return adj;
  }

  /** The fused system of each ring, numbered in the order of the first ring of each system. */
  private static int[] fusedSystemOf(List<List<Integer>> adj) {
    int nr = adj.size();
    int[] comp = new int[nr];
    Arrays.fill(comp, -1);
    int nComp = 0;
    for (int i = 0; i < nr; i++) {
      if (comp[i] >= 0) continue;
      comp[i] = nComp;
      Deque<Integer> q = new ArrayDeque<>();
      q.add(i);
      while (!q.isEmpty()) {
        int u = q.poll();
        for (int v : adj.get(u)) {
          if (comp[v] < 0) { comp[v] = nComp; q.add(v); }
        }
      }
      nComp++;
    }
    return comp;
  }

  private static int systemCount(int[] systemOf) {
    int count = 0;
    for (int system : systemOf) count = Math.max(count, system + 1);
    return count;
  }
}
