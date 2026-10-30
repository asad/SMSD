/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
package com.bioinception.smsd.core;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

/**
 * The relevant cycles of a graph grouped into families, with the membership test and the member
 * enumeration of {@link MolGraph}, as in the native core. Every family is named by its root, the
 * highest-numbered atom of its members, and by the atoms opposite the root; its members are the
 * cycles made of two shortest paths from the root that stay among atoms numbered no higher than
 * the root, and its member count is exact.
 *
 * <p>Cost per root, for a ring system of a atoms and b bonds: O(a + b) for two breadth-first
 * searches, O(b ceil(a/64)) for the atom sets of the shortest paths, one union per bond of their
 * graph, with an exact addition of path counts along the same bonds, O(b ceil(a/64)) for families
 * of odd size and O(sum over atoms x of C(p_x, 2) ceil(a/64)) for families of even size, p_x being
 * the shortest-path parents of x. The relevance test costs O(prototypes * rank * ceil(|E|/64)).
 * Every atom of a ring system is a root, but only the highest atom of a system that is one simple
 * ring.
 * The work is counted against {@link #FAMILY_WORK_LIMIT}, tested after each root and after each
 * reduction; a ring system that exceeds it is left out with every later one, and the result is
 * marked truncated.
 */
final class RelevantCycleFamilySearch {
  /** Work allowance of one family enumeration. */
  static final long FAMILY_WORK_LIMIT = 1L << 29;
  /** The atom sets of one ring system's searches may take at most this many words. */
  static final long ATOM_SET_WORD_LIMIT = 1L << 23;

  private RelevantCycleFamilySearch() {}

  /** One family found for a root, before the relevance test. */
  private static final class Prototype {
    final int size, root, first, second, opposite;
    final int[] atoms;
    final BigInteger members;
    final long discovery;

    Prototype(int size, int root, int first, int second, int opposite, int[] atoms, BigInteger members,
              long discovery) {
      this.size = size;
      this.root = root;
      this.first = first;
      this.second = second;
      this.opposite = opposite;
      this.atoms = atoms;
      this.members = members;
      this.discovery = discovery;
    }
  }

  static RelevantCycleFamilies compute(int n, int[][] neighbors) {
    return compute(n, neighbors, FAMILY_WORK_LIMIT, null);
  }

  /** As {@link #compute(int, int[][])}, with the limit a test may lower; work[0] receives the count. */
  static RelevantCycleFamilies compute(int n, int[][] neighbors, long workLimit, long[] work) {
    Enumeration enumeration = new Enumeration(n, neighbors, workLimit);
    RelevantCycleFamilies families = enumeration.run();
    if (work != null) work[0] = enumeration.work;
    return families;
  }

  private static final class Enumeration {
    final int n;
    final int[][] neighbors;
    final long limit;
    long work, discovery;
    int[] ofAtom, local, distS, distR, parentCount, queue, order, pivotRow;
    int[][] parents;
    long[] reach = new long[0], small;
    BigInteger[] large;
    RingCandidateSearch.Bonds bonds;
    int edgeWords;

    Enumeration(int n, int[][] neighbors, long limit) {
      this.n = n;
      this.neighbors = neighbors;
      this.limit = limit;
    }

    RelevantCycleFamilies run() {
      List<RelevantCycleFamily> families = new ArrayList<>();
      if (n == 0 || RingPerception.cycleRank(n, neighbors) <= 0)
        return new RelevantCycleFamilies(families, BigInteger.ZERO, false);
      RingPerception.Systems systems = RingPerception.systems(n, neighbors);
      if (systems == null) return new RelevantCycleFamilies(families, BigInteger.ZERO, true);
      ofAtom = systems.ofAtom;
      bonds = new RingCandidateSearch.Bonds(n, neighbors);
      edgeWords = Math.max(1, (bonds.count + 63) >>> 6);
      int count = systems.atoms.length;
      int[][] members = new int[count][];
      for (int s = 0; s < count; s++) members[s] = new int[systems.atoms[s]];
      int[] fill = new int[count];
      local = new int[n];
      for (int atom = 0; atom < n; atom++) {
        int s = ofAtom[atom];
        if (s >= 0) { local[atom] = fill[s]; members[s][fill[s]++] = atom; }
      }
      distS = new int[n];
      distR = new int[n];
      Arrays.fill(distS, -1);
      Arrays.fill(distR, -1);
      parentCount = new int[n];
      parents = new int[n][];
      queue = new int[n];
      order = new int[n];
      small = new long[n];
      large = new BigInteger[n];
      pivotRow = new int[bonds.count];
      Arrays.fill(pivotRow, -1);
      List<Prototype> kept = new ArrayList<>();
      boolean truncated = false;
      for (int s = 0; s < count; s++) {
        List<Prototype> found = new ArrayList<>();
        if (!collect(s, members[s], systems.simpleRing[s], found)) { truncated = true; break; }
        List<Prototype> relevant = work > limit ? null : relevantOf(found);
        if (relevant == null) { truncated = true; break; }
        kept.addAll(relevant);
      }
      Collections.sort(kept, (p, q) -> p.size != q.size ? Integer.compare(p.size, q.size)
          : p.root != q.root ? Integer.compare(p.root, q.root) : Long.compare(p.discovery, q.discovery));
      BigInteger total = BigInteger.ZERO;
      for (Prototype p : kept) {
        families.add(new RelevantCycleFamily(p.size, p.root, p.first, p.second, p.opposite, p.atoms, p.members));
        total = total.add(p.members);
      }
      return new RelevantCycleFamilies(families, total, truncated);
    }

    /**
     * The prototypes of the roots of system s in root order, every atom or, for a system that is
     * one simple ring, its highest atom only; false when the work runs out.
     */
    private boolean collect(int s, int[] atomsOfSystem, boolean simpleRing, List<Prototype> found) {
      int a = atomsOfSystem.length, aw = (a + 63) >>> 6;
      if ((long) a * aw > ATOM_SET_WORD_LIMIT) return false;
      if (reach.length < a * aw) reach = new long[a * aw];
      for (int rootAt = simpleRing ? a - 1 : 0; rootAt < a; rootAt++) {
        int r = atomsOfSystem[rootAt];
        for (int u : atomsOfSystem) {
          distS[u] = -1;
          distR[u] = -1;
          parentCount[u] = 0;
        }
        int tail = 0;
        queue[tail++] = r;
        distS[r] = 0;
        for (int head = 0; head < tail; head++) {
          int u = queue[head];
          work += 1 + neighbors[u].length;
          for (int w : neighbors[u])
            if (ofAtom[w] == s && distS[w] < 0) { distS[w] = distS[u] + 1; queue[tail++] = w; }
        }
        // Over the atoms numbered up to r, keeping every shortest-path parent in dequeue order.
        int reached = 0;
        order[reached++] = r;
        distR[r] = 0;
        for (int head = 0; head < reached; head++) {
          int u = order[head];
          work += 1 + neighbors[u].length;
          for (int w : neighbors[u]) {
            if (ofAtom[w] != s || w > r) continue;
            if (distR[w] < 0) {
              distR[w] = distR[u] + 1;
              addParent(w, u);
              order[reached++] = w;
            } else if (distR[w] == distR[u] + 1) {
              addParent(w, u);
            }
          }
        }
        int top = (local[r] >>> 6) + 1;
        for (int i = 0; i < reached; i++) {
          int u = order[i], at = local[u] * aw;
          Arrays.fill(reach, at, at + top, 0L);
          reach[at + (local[u] >>> 6)] |= 1L << (local[u] & 63);
          large[u] = null;
          small[u] = u == r ? 1 : 0;
          for (int k = 0; k < parentCount[u]; k++) {
            int p = parents[u][k], pt = local[p] * aw;
            for (int word = 0; word < top; word++) reach[at + word] |= reach[pt + word];
            work += aw + Math.max(limbs(u), limbs(p));
            add(u, p);
          }
        }
        for (int i = 0; i < reached; i++) {
          int y = order[i];
          if (y == r || distR[y] != distS[y]) continue;
          for (int z : neighbors[y]) {
            if (z <= y || ofAtom[z] != s || z > r || distR[z] != distR[y] || distR[z] != distS[z]) continue;
            work += aw;
            if (disjoint(y, z, r, aw, top)) found.add(prototype(r, y, -1, z, 2 * distR[y] + 1));
          }
        }
        for (int i = 0; i < reached; i++) {
          int x = order[i];
          if (x == r || distR[x] < 2 || distR[x] != distS[x]) continue;
          int[] list = parents[x];
          for (int p = 0; p < parentCount[x]; p++)
            for (int q = p + 1; q < parentCount[x]; q++) {
              int y = Math.min(list[p], list[q]), z = Math.max(list[p], list[q]);
              work += aw;
              if (disjoint(y, z, r, aw, top)) found.add(prototype(r, y, x, z, 2 * distR[x]));
            }
        }
        if (work > limit) return false;
      }
      return true;
    }

    private void addParent(int atom, int parent) {
      int[] list = parents[atom];
      if (list == null) parents[atom] = list = new int[4];
      else if (parentCount[atom] == list.length) parents[atom] = list = Arrays.copyOf(list, 2 * list.length);
      list[parentCount[atom]++] = parent;
    }

    /** The 32-bit limbs of the path count of an atom. */
    private int limbs(int u) {
      int bits = large[u] != null ? large[u].bitLength() : 64 - Long.numberOfLeadingZeros(small[u]);
      return (bits + 31) >>> 5;
    }

    private void add(int u, int p) {
      if (large[u] == null && large[p] == null) {
        long sum = small[u] + small[p];
        if (sum >= 0) { small[u] = sum; return; }
      }
      large[u] = value(u).add(value(p));
    }

    private BigInteger value(int u) {
      return large[u] != null ? large[u] : BigInteger.valueOf(small[u]);
    }

    /** Whether the shortest paths to y and z meet only at r. */
    private boolean disjoint(int y, int z, int r, int aw, int top) {
      int ya = local[y] * aw, za = local[z] * aw, rootWord = local[r] >>> 6;
      long rootBit = 1L << (local[r] & 63);
      for (int i = 0; i < top; i++) {
        long both = reach[ya + i] & reach[za + i];
        if (i == rootWord) both &= ~rootBit;
        if (both != 0) return false;
      }
      return true;
    }

    /** The family of root r opposite y and z (and x for an even size), its prototype by first parents. */
    private Prototype prototype(int r, int y, int x, int z, int size) {
      int[] atoms = new int[size];
      int u = y;
      for (int i = distR[y]; i > 0; i--) { atoms[i] = u; u = parents[u][0]; }
      atoms[0] = r;
      int at = distR[y] + 1;
      if (x >= 0) atoms[at++] = x;
      for (u = z; u != r; u = parents[u][0]) atoms[at++] = u;
      work += (long) limbs(y) * limbs(z) + size;
      return new Prototype(size, r, y, z, x, atoms, value(y).multiply(value(z)), discovery++);
    }

    /**
     * The relevant prototypes of one system: those independent of all prototypes of strictly
     * smaller size, the prototypes of each size entering the elimination after that size is
     * tested. Null when the work runs out, which is tested after each reduction.
     */
    private List<Prototype> relevantOf(List<Prototype> found) {
      List<Prototype> sorted = new ArrayList<>(found);
      Collections.sort(sorted, (p, q) -> Integer.compare(p.size, q.size));
      List<Prototype> relevant = new ArrayList<>();
      List<long[]> rows = new ArrayList<>();
      int[] pivots = new int[8];
      boolean exceeded = false;
      for (int i = 0; i < sorted.size() && !exceeded; ) {
        int j = i;
        while (j < sorted.size() && sorted.get(j).size == sorted.get(i).size) j++;
        for (int k = i; k < j && !exceeded; k++) {
          if (reduce(mask(sorted.get(k)), rows) >= 0) relevant.add(sorted.get(k));
          exceeded = work > limit;
        }
        for (int k = i; k < j && !exceeded; k++) {
          long[] vector = mask(sorted.get(k));
          int pivot = reduce(vector, rows);
          if (pivot >= 0) {
            if (rows.size() == pivots.length) pivots = Arrays.copyOf(pivots, 2 * pivots.length);
            pivots[rows.size()] = pivot;
            pivotRow[pivot] = rows.size();
            rows.add(vector);
          }
          exceeded = work > limit;
        }
        i = j;
      }
      for (int row = 0; row < rows.size(); row++) pivotRow[pivots[row]] = -1;
      return exceeded ? null : relevant;
    }

    /** Reduces until a bond without a row is met; returns that bond, or -1 for zero. */
    private int reduce(long[] vector, List<long[]> rows) {
      for (int word = 0; word < edgeWords; word++)
        while (vector[word] != 0) {
          int bit = (word << 6) + Long.numberOfTrailingZeros(vector[word]);
          int row = pivotRow[bit];
          if (row < 0) return bit;
          long[] basis = rows.get(row);
          for (int k = word; k < edgeWords; k++) vector[k] ^= basis[k];
          work += edgeWords;
        }
      return -1;
    }

    private long[] mask(Prototype prototype) {
      long[] vector = new long[edgeWords];
      int[] atoms = prototype.atoms;
      for (int i = 0; i < atoms.length; i++) {
        int e = bonds.between(neighbors, atoms[i], atoms[(i + 1) % atoms.length]);
        vector[e >>> 6] ^= 1L << (e & 63);
      }
      return vector;
    }
  }

  /** Distances from a root over its system and over the atoms numbered up to it, with parents. */
  private static final class RootSearch {
    final int[] distS, distR, parentCount;
    final int[][] parents;

    RootSearch(int n, int[][] neighbors, int[] ofAtom, int r, boolean withSystem) {
      int s = ofAtom[r];
      distS = new int[n];
      distR = new int[n];
      parentCount = new int[n];
      parents = new int[n][];
      Arrays.fill(distS, -1);
      Arrays.fill(distR, -1);
      int[] queue = new int[n];
      if (withSystem) {
        int tail = 0;
        queue[tail++] = r;
        distS[r] = 0;
        for (int head = 0; head < tail; head++) {
          int u = queue[head];
          for (int w : neighbors[u])
            if (ofAtom[w] == s && distS[w] < 0) { distS[w] = distS[u] + 1; queue[tail++] = w; }
        }
      }
      int tail = 0;
      queue[tail++] = r;
      distR[r] = 0;
      for (int head = 0; head < tail; head++) {
        int u = queue[head];
        for (int w : neighbors[u]) {
          if (ofAtom[w] != s || w > r) continue;
          if (distR[w] < 0) {
            distR[w] = distR[u] + 1;
            add(w, u);
            queue[tail++] = w;
          } else if (distR[w] == distR[u] + 1) {
            add(w, u);
          }
        }
      }
    }

    private void add(int atom, int parent) {
      int[] list = parents[atom];
      if (list == null) parents[atom] = list = new int[4];
      else if (parentCount[atom] == list.length) parents[atom] = list = Arrays.copyOf(list, 2 * list.length);
      list[parentCount[atom]++] = parent;
    }

    boolean isParent(int atom, int parent) {
      for (int k = 0; k < parentCount[atom]; k++) if (parents[atom][k] == parent) return true;
      return false;
    }
  }

  /** The shortest paths from one atom back to the root, in a fixed order, one at a time. */
  private static final class Paths {
    final RootSearch search;
    final int[] atoms, choice;
    final int length;

    Paths(RootSearch search, int target) {
      this.search = search;
      length = search.distR[target];
      atoms = new int[length + 1];
      choice = new int[length + 1];
      atoms[0] = target;
      for (int i = 0; i < length; i++) atoms[i + 1] = search.parents[atoms[i]][0];
    }

    /** Moves to the next path; false when there is none. */
    boolean advance() {
      for (int i = length - 1; i >= 0; i--)
        if (choice[i] + 1 < search.parentCount[atoms[i]]) {
          choice[i]++;
          atoms[i + 1] = search.parents[atoms[i]][choice[i]];
          for (int j = i + 1; j < length; j++) {
            choice[j] = 0;
            atoms[j + 1] = search.parents[atoms[j]][0];
          }
          return true;
        }
      return false;
    }
  }

  private static boolean adjacent(int[][] neighbors, int u, int w) {
    for (int x : neighbors[u]) if (x == w) return true;
    return false;
  }

  /**
   * The index in families of the family that has ring as a member, or -1. The ring passes only
   * when, read both ways from its highest atom, the atom at every position up to the opposite
   * atoms is at that distance from the highest atom both in its ring system and among the atoms
   * numbered no higher, and a listed family is named by those atoms.
   */
  static int familyOf(int n, int[][] neighbors, RelevantCycleFamilies families, int[] ring) {
    if (families == null || ring == null || ring.length < 3) return -1;
    int l = ring.length, root = -1, start = -1;
    for (int i = 0; i < l; i++) {
      int u = ring[i];
      if (u < 0 || u >= n) return -1;
      if (u > root) { root = u; start = i; }
    }
    boolean[] seen = new boolean[n];
    for (int i = 0; i < l; i++) {
      int u = ring[i];
      if (seen[u] || !adjacent(neighbors, u, ring[(i + 1) % l])) return -1;
      seen[u] = true;
    }
    RingPerception.Systems systems = RingPerception.systems(n, neighbors);
    if (systems == null) return -1;
    int s = systems.ofAtom[root];
    if (s < 0) return -1;
    for (int u : ring) if (systems.ofAtom[u] != s) return -1;
    RootSearch search = new RootSearch(n, neighbors, systems.ofAtom, root, true);
    int k = l / 2;
    for (int i = 0; i <= k; i++) {
      int a = ring[(start + i) % l], b = ring[((start - i) % l + l) % l];
      if (search.distS[a] != i || search.distR[a] != i || search.distS[b] != i || search.distR[b] != i) return -1;
    }
    int y, z, opposite;
    if ((l & 1) == 1) {
      y = ring[(start + k) % l];
      z = ring[((start - k) % l + l) % l];
      opposite = -1;
    } else {
      opposite = ring[(start + k) % l];
      y = ring[(start + k - 1) % l];
      z = ring[((start - k + 1) % l + l) % l];
    }
    int first = Math.min(y, z), second = Math.max(y, z);
    List<RelevantCycleFamily> list = families.families();
    for (int index = 0; index < list.size(); index++) {
      RelevantCycleFamily family = list.get(index);
      if (family.size() == l && family.root() == root && family.first() == first && family.second() == second
          && family.opposite() == opposite) return index;
    }
    return -1;
  }

  /**
   * Up to limit members of family, always in the same order: paths to first() vary slowest, and
   * along each path the choice nearest the root varies fastest; the first member is the prototype.
   */
  static int[][] members(int n, int[][] neighbors, RelevantCycleFamily family, int limit) {
    int[][] none = new int[0][];
    if (family == null || limit <= 0) return none;
    int r = family.root(), y = family.first(), z = family.second(), x = family.opposite(), size = family.size();
    if (size < 3 || r < 0 || r >= n || y < 0 || y >= z || z >= n || x < -1 || x >= n) return none;
    RingPerception.Systems systems = RingPerception.systems(n, neighbors);
    if (systems == null) return none;
    int s = systems.ofAtom[r];
    if (s < 0 || systems.ofAtom[y] != s || systems.ofAtom[z] != s || (x >= 0 && systems.ofAtom[x] != s)) return none;
    RootSearch search = new RootSearch(n, neighbors, systems.ofAtom, r, true);
    boolean odd = (size & 1) == 1;
    int level = odd ? (size - 1) / 2 : size / 2 - 1;
    // The opposite atoms lie on shortest paths from the root, in the system and below the root.
    if (search.distR[y] != level || search.distR[z] != level || search.distS[y] != level || search.distS[z] != level)
      return none;
    if (odd ? x != -1 || !adjacent(neighbors, y, z)
            : x < 0 || search.distR[x] != level + 1 || search.distS[x] != level + 1
              || !search.isParent(x, y) || !search.isParent(x, z))
      return none;
    List<int[]> found = new ArrayList<>();
    Paths toFirst = new Paths(search, y);
    for (boolean more = true; more && found.size() < limit; more = toFirst.advance()) {
      Paths toSecond = new Paths(search, z);
      for (boolean next = true; next && found.size() < limit; next = toSecond.advance()) {
        int[] member = new int[size];
        int at = 0;
        for (int i = level; i >= 0; i--) member[at++] = toFirst.atoms[i];
        if (x >= 0) member[at++] = x;
        for (int i = 0; i < level; i++) member[at++] = toSecond.atoms[i];
        found.add(member);
      }
    }
    return found.toArray(new int[0][]);
  }
}
