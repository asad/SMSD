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
 * The candidate cycles of the exact ring basis, read from one shortest-path tree per atom of every
 * ring system that is not a single ring, as the native core reads them.
 *
 * <p>The pairwise search this stands for pairs every atom of a ring system with every bond of that
 * system and joins the shortest paths from the atom to both ends of the bond that avoid it. One
 * breadth-first tree per atom answers every such pair, with equal shortest paths broken as the
 * pairwise search breaks them, by discovery order and then by neighbour order. The candidate
 * list, its order and the atom order of every candidate are therefore the same, and so is every
 * result. The fallback decision is the pairwise search's own: it falls back exactly when the work
 * it counts before its last pair exceeds the limit. That work is found pair by pair, exactly or
 * within a bracket, and a bracketed pair is searched pairwise only when the brackets leave the
 * decision open.
 *
 * <p>Cost per root, for a system of a atoms and b bonds with cycle rank r, tree depth D,
 * w = ceil(b/64) and largest degree Delta: O(a + b) for the tree and its numbering, O(a + r) for
 * ordering the closing bonds, O((a + r) alpha(a)) for the bypasses of the tree bonds, O(a w) for
 * the tree-path bond sets when they are kept, and per bond O(1) for its charge and validity, O(w)
 * for its key and O(Delta D) for the climb of a valid tree bond; O(len + w) per stored candidate
 * and O(w) per stored key sharing the 64-bit hash. Selection costs O(w) per row reduced, with at
 * most r rows per candidate, and a pair searched pairwise costs that search. The worst case per
 * system is O(a b (w + Delta D)), against Theta(a b (a + b)) for the pairwise search.
 */
final class RingCandidateSearch {
  /** Pair classes recorded in {@link RingSearchNotes}. */
  static final int PAIR_NON_TREE = 0, PAIR_A = 1, PAIR_B1 = 2, PAIR_B1_PLUS = 3, PAIR_B2 = 4,
      PAIR_SIMPLE_RING = 5;
  /** Tree-path bond sets are kept for a root only when they take at most this many words. */
  static final long PATH_SET_WORD_LIMIT = 1L << 22;
  /** Stored key words allowed while a refused basis is completed. */
  static final long LIFT_KEY_WORD_LIMIT = 1L << 24;
  /** The upper end of a loosened bracket. */
  static final long UNBOUNDED = Long.MAX_VALUE;
  private static final int UNSEEN = -2;

  /** Internal switches of the ring search; every one is off in {@link #DEFAULTS}. */
  static final class RingSearchOptions {
    static final RingSearchOptions DEFAULTS = new RingSearchOptions(false, false, false, false);
    /** Complete a basis the fallback leaves unproved, keeping every other field of the result. */
    final boolean liftRefusedBasis;
    /** Clear the relevant-cycle truncation flag when the families prove the list complete. */
    final boolean proveRelevantByFamilies;
    /** Test only: every bracketed pair is bracketed by zero and no upper bound. */
    final boolean loosenBracketsForTesting;
    /** Test only: the work before the last pair is always found exactly. */
    final boolean replayAllForTesting;

    RingSearchOptions(boolean liftRefusedBasis, boolean proveRelevantByFamilies,
                      boolean loosenBracketsForTesting, boolean replayAllForTesting) {
      this.liftRefusedBasis = liftRefusedBasis;
      this.proveRelevantByFamilies = proveRelevantByFamilies;
      this.loosenBracketsForTesting = loosenBracketsForTesting;
      this.replayAllForTesting = replayAllForTesting;
    }
  }

  /** What one search did, for tests and diagnostics. */
  static final class RingSearchNotes {
    /** Atom-bond pairs of the pairwise search. */
    long pairs;
    /** Pairs whose charges were summed; the last pair never is. */
    long pairsSummed;
    /** Bounds on the work the pairwise search counts before its last pair. */
    long baseWorkLow, baseWorkHigh;
    /** The upper bound is not known. */
    boolean baseWorkHighUnbounded;
    /** The bounds are equal, so the work is known exactly. */
    boolean baseWorkExact;
    /** The search fell back to the bounded pipeline. */
    boolean fallback;
    /** Index of the pair after which the first pass stopped, or -1. */
    long stoppedAfterPair = -1;
    /** Work of this search, in atoms and neighbours read and in words read or written. */
    long searchWork;
    long nonTreePairs, classAPairs, classB1Pairs, classB1PlusPairs, classB2Pairs, simpleRingPairs;
    /** Roots searched again to settle the decision, and the pairs searched pairwise there. */
    long replayedRoots, replayedPairs;
    /** Distinct candidate cycles stored. */
    long candidates;
    /** The basis of the fallback was not proved minimum. */
    boolean baseBasisTruncated;
    /** That basis was completed from the full candidate list and proved. */
    boolean basisLifted;
    /** When not null, every pair of the first pass is recorded here. */
    List<PairRecord> trace;
  }

  /** One atom-bond pair, as recorded in {@link RingSearchNotes#trace}. */
  static final class PairRecord {
    int root, bond, kind;
    boolean valid;
    long bfsLow, bfsHigh, pathLength;
    int[] firstPath, secondPath;
  }

  /** The bonds of a graph in edge order: lower atom first, then neighbour order. */
  static final class Bonds {
    final int count;
    final int[] first, second;
    /** The bond of every neighbour slot. */
    final int[][] ofSlot;

    /** For graphs whose atoms list neither themselves nor a neighbour twice. */
    Bonds(int n, int[][] neighbors) {
      int total = 0;
      int[] start = new int[n + 1];
      for (int i = 0; i < n; i++)
        for (int j : neighbors[i])
          if (j > i) { total++; start[j + 1]++; }
      count = total;
      first = new int[total];
      second = new int[total];
      ofSlot = new int[n][];
      for (int j = 0; j < n; j++) start[j + 1] += start[j];
      int[] fill = Arrays.copyOf(start, n);
      int[] pending = new int[total];
      int bond = 0;
      for (int i = 0; i < n; i++) {
        int[] list = neighbors[i];
        ofSlot[i] = new int[list.length];
        Arrays.fill(ofSlot[i], -1);
        for (int slot = 0; slot < list.length; slot++) {
          int j = list[slot];
          if (j <= i) continue;
          first[bond] = i;
          second[bond] = j;
          ofSlot[i][slot] = bond;
          pending[fill[j]++] = bond;
          bond++;
        }
      }
      int[] mark = new int[n], slotOf = new int[n];
      Arrays.fill(mark, -1);
      for (int j = 0; j < n; j++) {
        int[] list = neighbors[j];
        for (int slot = 0; slot < list.length; slot++) {
          int i = list[slot];
          if (i < j) { mark[i] = j; slotOf[i] = slot; }
        }
        for (int at = start[j]; at < start[j + 1]; at++) {
          int b = pending[at], i = first[b];
          if (mark[i] == j) ofSlot[j][slotOf[i]] = b;
        }
      }
    }

    /** The bond joining two neighbours, or -1. */
    int between(int[][] neighbors, int u, int w) {
      int[] list = neighbors[u];
      for (int slot = 0; slot < list.length; slot++) if (list[slot] == w) return ofSlot[u][slot];
      return -1;
    }
  }

  /** The 64-bit key hash shared with the native core. */
  static long hash(long[] words, int from, int count) {
    long h = 0;
    for (int i = 0; i < count; i++) h = (Long.rotateLeft(h, 5) ^ words[from + i]) * 0x9E3779B97F4A7C15L;
    h ^= h >>> 33;
    h *= 0xFF51AFD7ED558CCDL;
    h ^= h >>> 33;
    return h;
  }

  private final int n;
  private final int[][] neighbors;
  private final int[] ofAtom;
  private final RingPerception.Systems systems;
  private final int cycleRank;
  private final long workLimit;
  private final RingSearchOptions options;
  private final RingSearchNotes notes;
  private final Bonds bonds;
  private final int totalEdges;
  /** The system-local number of every bond, -1 for a bond outside every system. */
  private final int[] localBond;
  private final int[][] systemEdges;
  private final int[] lowestAtom, words, systemRank;
  private final long pairCount;
  private final int lastRoot;
  private final long pathSetWordLimit;

  // The tree of the current root.
  private final int[] par, dep, rk, pre, subtreeSize, topChild, parentBond, order;
  private final long[] subtreeWork, prefixWork;
  private long[] levelWork = new long[2];
  private int orderSize, maxDepth;
  private int[] halfFrom, halfTo, halfBond, halfOrder;
  private int halfCount;
  private int[] lengthCount = new int[8];
  private final int[] setParent, setRank, setTop, coverFrom, coverTo, coverBond;
  private final int[] stackAtom, stackSlot;
  private long[] pathSets = new long[0];
  private boolean pathSetsKept;
  private final int[] climb, climbBond;
  private int climbLength;
  private long climbWork;

  // The candidates, in the order they are first found.
  private final List<int[]> candidateAtoms = new ArrayList<>();
  private int[] candidateSystem = new int[16], candidateKey = new int[16];
  private long[] candidateHash = new long[16];
  private final long[][] keyStore;
  private final int[] keyCount;
  private long storedKeyWords;
  private long[] tableHash = new long[64];
  private int[] tableEntry = new int[64];
  private final long[] key;

  private long searchWork;
  private RingPerception.PathSearch pathSearch;
  private int[] mark;
  private int markStamp;
  private int[] replayVisited, replayQueue;
  private int replayStamp;
  private long[] b2Low, b2High;
  private int[] b2Count;
  private final long[] classCounts = new long[6];

  // The pair just read.
  private int pairKind;
  private boolean pairValid;
  private long pairBfsLow, pairBfsHigh, pairPath;

  RingCandidateSearch(int n, int[][] neighbors, RingPerception.Systems systems, int cycleRank,
                      long workLimit, RingSearchOptions options, RingSearchNotes notes) {
    this(n, neighbors, systems, cycleRank, workLimit, options, notes, PATH_SET_WORD_LIMIT);
  }

  private RingCandidateSearch(int n, int[][] neighbors, RingPerception.Systems systems, int cycleRank,
                              long workLimit, RingSearchOptions options, RingSearchNotes notes,
                              long pathSetWordLimit) {
    this.pathSetWordLimit = pathSetWordLimit;
    this.n = n;
    this.neighbors = neighbors;
    this.systems = systems;
    this.ofAtom = systems.ofAtom;
    this.cycleRank = cycleRank;
    this.workLimit = workLimit;
    this.options = options != null ? options : RingSearchOptions.DEFAULTS;
    this.notes = notes;
    bonds = new Bonds(n, neighbors);
    totalEdges = bonds.count;
    int count = systems.atoms.length;
    int[] sizes = new int[count];
    localBond = new int[totalEdges];
    for (int e = 0; e < totalEdges; e++) {
      int s = ofAtom[bonds.first[e]];
      localBond[e] = s >= 0 && ofAtom[bonds.second[e]] == s ? sizes[s]++ : -1;
    }
    systemEdges = new int[count][];
    for (int s = 0; s < count; s++) systemEdges[s] = new int[sizes[s]];
    for (int e = 0; e < totalEdges; e++)
      if (localBond[e] >= 0) systemEdges[ofAtom[bonds.first[e]]][localBond[e]] = e;
    lowestAtom = new int[count];
    Arrays.fill(lowestAtom, -1);
    for (int atom = n - 1; atom >= 0; atom--) if (ofAtom[atom] >= 0) lowestAtom[ofAtom[atom]] = atom;
    words = new int[count];
    systemRank = new int[count];
    int widest = 1, mostBonds = 0;
    for (int s = 0; s < count; s++) {
      words[s] = (sizes[s] + 63) >>> 6;
      systemRank[s] = sizes[s] - systems.atoms[s] + 1;
      widest = Math.max(widest, words[s]);
      mostBonds = Math.max(mostBonds, sizes[s]);
    }
    long pairsFound = 0;
    int last = -1;
    for (int v = 0; v < n; v++) {
      int s = ofAtom[v];
      if (s < 0) continue;
      int pairs = !systems.simpleRing[s] ? sizes[s] : v == lowestAtom[s] && sizes[s] > 0 ? 1 : 0;
      if (pairs > 0) { pairsFound += pairs; last = v; }
    }
    pairCount = pairsFound;
    lastRoot = last;
    key = new long[widest];
    keyStore = new long[count][];
    keyCount = new int[count];
    par = new int[n];
    Arrays.fill(par, UNSEEN);
    dep = new int[n];
    rk = new int[n];
    pre = new int[n];
    subtreeSize = new int[n];
    topChild = new int[n];
    parentBond = new int[n];
    order = new int[n];
    subtreeWork = new long[n];
    prefixWork = new long[n];
    int halves = Math.max(2, 2 * mostBonds);
    halfFrom = new int[halves];
    halfTo = new int[halves];
    halfBond = new int[halves];
    halfOrder = new int[halves];
    setParent = new int[n];
    setRank = new int[n];
    setTop = new int[n];
    coverFrom = new int[n];
    coverTo = new int[n];
    coverBond = new int[n];
    stackAtom = new int[n];
    stackSlot = new int[n];
    climb = new int[n];
    climbBond = new int[n];
    Arrays.fill(tableEntry, -1);
  }

  private long weight(int atom) {
    return 1L + neighbors[atom].length;
  }

  private boolean ancestor(int u, int y) {
    return pre[u] <= pre[y] && pre[y] < pre[u] + subtreeSize[u];
  }

  private long levelWorkUpTo(long level) {
    if (level < 0) return 0;
    return levelWork[(int) Math.min(level, maxDepth)];
  }

  /**
   * The breadth-first tree of root v over system s: depths, ranks, the work of the atoms up to each
   * rank and each level, a depth-first numbering for ancestor tests, the closing half-bonds in
   * order of the cycle they close and then of the numbering walk, and for every atom but v the
   * first of those half-bonds that enters its subtree from outside.
   */
  private void buildTree(int v, int s, boolean withPathSets) {
    for (int i = 0; i < orderSize; i++) par[order[i]] = UNSEEN;
    order[0] = v;
    orderSize = 1;
    par[v] = -1;
    dep[v] = 0;
    parentBond[v] = -1;
    for (int head = 0; head < orderSize; head++) {
      int u = order[head];
      rk[u] = head;
      int[] list = neighbors[u];
      searchWork += 1 + list.length;
      for (int slot = 0; slot < list.length; slot++) {
        int w = list[slot];
        if (ofAtom[w] != s || par[w] != UNSEEN) continue;
        par[w] = u;
        dep[w] = dep[u] + 1;
        parentBond[w] = bonds.ofSlot[u][slot];
        order[orderSize++] = w;
      }
    }
    long total = 0;
    for (int i = 0; i < orderSize; i++) {
      total += weight(order[i]);
      prefixWork[i] = total;
    }
    maxDepth = dep[order[orderSize - 1]];
    if (levelWork.length < maxDepth + 1) levelWork = new long[Math.max(maxDepth + 1, 2 * levelWork.length)];
    Arrays.fill(levelWork, 0, maxDepth + 1, 0L);
    for (int i = 0; i < orderSize; i++) levelWork[dep[order[i]]] += weight(order[i]);
    for (int d = 1; d <= maxDepth; d++) levelWork[d] += levelWork[d - 1];

    // Depth-first numbering over the tree in neighbour order; every other slot inside the system
    // except the parent's is a closing half-bond, listed in walk order.
    int counter = 0, depth = 0;
    halfCount = 0;
    pre[v] = counter++;
    stackAtom[0] = v;
    stackSlot[0] = 0;
    searchWork += weight(v);
    while (depth >= 0) {
      int x = stackAtom[depth];
      int[] list = neighbors[x];
      if (stackSlot[depth] == list.length) { depth--; continue; }
      int slot = stackSlot[depth]++;
      int w = list[slot];
      if (ofAtom[w] != s || w == par[x]) continue;
      if (par[w] == x) {
        pre[w] = counter++;
        depth++;
        stackAtom[depth] = w;
        stackSlot[depth] = 0;
        searchWork += weight(w);
      } else {
        if (halfCount == halfFrom.length) growHalves();
        halfFrom[halfCount] = x;
        halfTo[halfCount] = w;
        halfBond[halfCount] = bonds.ofSlot[x][slot];
        halfCount++;
      }
    }
    for (int i = 0; i < orderSize; i++) {
      int u = order[i];
      subtreeSize[u] = 1;
      subtreeWork[u] = weight(u);
      topChild[u] = u == v ? -1 : par[u] == v ? u : topChild[par[u]];
    }
    for (int i = orderSize - 1; i > 0; i--) {
      int u = order[i];
      subtreeSize[par[u]] += subtreeSize[u];
      subtreeWork[par[u]] += subtreeWork[u];
    }

    // Stable counting sort of the half-bonds by the length of the cycle they close.
    int keys = 2 * maxDepth + 3;
    if (lengthCount.length < keys + 1) lengthCount = new int[keys + 1];
    Arrays.fill(lengthCount, 0, keys + 1, 0);
    for (int h = 0; h < halfCount; h++) lengthCount[dep[halfFrom[h]] + dep[halfTo[h]] + 2]++;
    for (int i = 1; i <= keys; i++) lengthCount[i] += lengthCount[i - 1];
    for (int h = 0; h < halfCount; h++) halfOrder[lengthCount[dep[halfFrom[h]] + dep[halfTo[h]] + 1]++] = h;
    searchWork += halfCount + 2L * maxDepth + 1;

    // Each atom takes the first half-bond, in that order, entering its subtree from outside. Sets
    // join covered atoms to their parents; the top of a set is its highest atom, still uncovered.
    for (int i = 0; i < orderSize; i++) {
      int u = order[i];
      setParent[u] = u;
      setRank[u] = 0;
      setTop[u] = u;
      coverFrom[u] = -1;
      coverTo[u] = -1;
      coverBond[u] = -1;
    }
    for (int j = 0; j < halfCount; j++) {
      int h = halfOrder[j], x = halfFrom[h], w = halfTo[h];
      int u = setTop[find(w)];
      while (!ancestor(u, x)) {
        coverFrom[u] = x;
        coverTo[u] = w;
        coverBond[u] = halfBond[h];
        int child = find(u), parent = find(par[u]);
        int top = setTop[parent];
        if (setRank[child] < setRank[parent]) { int swap = child; child = parent; parent = swap; }
        setParent[parent] = child;
        if (setRank[child] == setRank[parent]) setRank[child]++;
        setTop[child] = top;
        searchWork++;
        u = top;
      }
    }

    int width = words[s];
    pathSetsKept = withPathSets && (long) orderSize * width <= pathSetWordLimit;
    if (pathSetsKept) {
      int need = orderSize * width;
      if (pathSets.length < need)
        pathSets = new long[Math.max(need, (int) Math.min(PATH_SET_WORD_LIMIT, 2L * pathSets.length))];
      Arrays.fill(pathSets, 0, width, 0L);
      for (int i = 1; i < orderSize; i++) {
        int u = order[i];
        System.arraycopy(pathSets, rk[par[u]] * width, pathSets, i * width, width);
        int b = localBond[parentBond[u]];
        pathSets[i * width + (b >>> 6)] ^= 1L << (b & 63);
      }
      searchWork += (long) orderSize * width;
    }
  }

  private void growHalves() {
    int size = 2 * halfFrom.length;
    halfFrom = Arrays.copyOf(halfFrom, size);
    halfTo = Arrays.copyOf(halfTo, size);
    halfBond = Arrays.copyOf(halfBond, size);
    halfOrder = Arrays.copyOf(halfOrder, size);
  }

  private int find(int u) {
    while (setParent[u] != u) {
      setParent[u] = setParent[setParent[u]];
      u = setParent[u];
      searchWork++;
    }
    return u;
  }

  private void flip(int localIndex) {
    key[localIndex >>> 6] ^= 1L << (localIndex & 63);
  }

  private void flipTreePath(int u, int v) {
    for (int x = u; x != v; x = par[x]) flip(localBond[parentBond[x]]);
  }

  /** Reads pair (v, e): its class, validity, charge or bracket and, when valid, its candidate. */
  private void pair(int v, int s, int e) {
    int a = bonds.first[e], b = bonds.second[e];
    searchWork++;
    if (par[a] == b || par[b] == a) treeBond(v, s, e, a, b);
    else crossBond(v, s, e, a, b);
  }

  /** A bond outside the tree: the pairwise search finds the tree paths to both ends. */
  private void crossBond(int v, int s, int e, int a, int b) {
    pairKind = PAIR_NON_TREE;
    pairValid = topChild[a] != topChild[b];
    int later = rk[a] > rk[b] ? a : b;
    pairBfsLow = pairBfsHigh = prefixWork[rk[par[later]]];
    pairPath = dep[a] + dep[b] + 2;
    if (!pairValid) return;
    int width = words[s], size = dep[a] + dep[b] + 1;
    if (pathSetsKept) {
      int at = rk[a] * width, bt = rk[b] * width;
      for (int i = 0; i < width; i++) key[i] = pathSets[at + i] ^ pathSets[bt + i];
      searchWork += width;
    } else {
      Arrays.fill(key, 0, width, 0L);
      flipTreePath(a, v);
      flipTreePath(b, v);
      searchWork += size + width;
    }
    flip(localBond[e]);
    long h = hash(key, 0, width);
    if (lookup(s, h) >= 0) return;
    int[] atoms = new int[size];
    for (int u = a, i = dep[a]; i >= 0; u = par[u], i--) atoms[i] = u;
    int at = dep[a] + 1;
    for (int u = b; u != v; u = par[u]) atoms[at++] = u;
    store(s, h, atoms);
  }

  /**
   * A tree bond (p, c): the pairwise search reaches c through the first half-bond (x, w) that
   * enters the subtree of c, and then climbs to c one level a bond, always by the first listed
   * neighbour.
   */
  private void treeBond(int v, int s, int e, int a, int b) {
    int c = par[a] == b ? a : b, p = c == a ? b : a;
    int x = coverFrom[c], w = coverTo[c];
    int closing = dep[x] + dep[w] + 1;
    pairValid = p == v || topChild[x] != topChild[p];
    pairPath = closing + 1;
    if (w == c) {
      int step = dep[x] - dep[c];
      pairKind = step < 0 ? PAIR_A : step == 0 ? PAIR_B1 : PAIR_B1_PLUS;
      pairBfsLow = pairBfsHigh = prefixWork[rk[x]] - (rk[c] < rk[x] ? weight(c) : 0);
    } else {
      pairKind = PAIR_B2;
      long level = closing - dep[c] - 1;
      pairBfsHigh = levelWorkUpTo(level) - weight(c);
      pairBfsLow = Math.max(0L, levelWorkUpTo(level - 1) - subtreeWork[c]);
    }
    if (!pairValid) {
      loosen();
      return;
    }
    climbTo(w, c, s);
    searchWork += climbWork;
    if (pairKind == PAIR_B2) pairBfsLow += climbWork;
    loosen();
    int width = words[s];
    if (pathSetsKept) {
      int xt = rk[x] * width, pt = rk[p] * width;
      for (int i = 0; i < width; i++) key[i] = pathSets[xt + i] ^ pathSets[pt + i];
      searchWork += width;
    } else {
      Arrays.fill(key, 0, width, 0L);
      flipTreePath(x, v);
      flipTreePath(p, v);
      searchWork += closing + width;
    }
    flip(localBond[coverBond[c]]);
    for (int i = 0; i + 1 < climbLength; i++) flip(localBond[climbBond[i]]);
    flip(localBond[e]);
    long h = hash(key, 0, width);
    if (lookup(s, h) >= 0) return;
    int[] atoms = new int[closing];
    int at;
    if (c == a) {
      for (int u = x, i = dep[x]; i >= 0; u = par[u], i--) atoms[i] = u;
      at = dep[x] + 1;
      for (int i = 0; i < climbLength; i++) atoms[at++] = climb[i];
      for (int u = p; u != v; u = par[u]) atoms[at++] = u;
    } else {
      for (int u = p, i = dep[p]; i >= 0; u = par[u], i--) atoms[i] = u;
      at = dep[p] + 1;
      for (int i = climbLength - 1; i >= 0; i--) atoms[at++] = climb[i];
      for (int u = x; u != v; u = par[u]) atoms[at++] = u;
    }
    store(s, h, atoms);
  }

  private void loosen() {
    if (pairKind == PAIR_B2 && options.loosenBracketsForTesting) {
      pairBfsLow = 0;
      pairBfsHigh = UNBOUNDED;
    }
  }

  /** The climb from w up to c inside the subtree of c, with the work of every atom read but c. */
  private void climbTo(int w, int c, int s) {
    climbLength = 0;
    climbWork = 0;
    int u = w;
    climb[climbLength++] = u;
    while (u != c) {
      int[] list = neighbors[u];
      climbWork += weight(u);
      int next = -1, bond = -1;
      for (int slot = 0; slot < list.length; slot++) {
        int y = list[slot];
        if (ofAtom[y] == s && dep[y] == dep[u] - 1 && ancestor(c, y)) {
          next = y;
          bond = bonds.ofSlot[u][slot];
          break;
        }
      }
      climbBond[climbLength - 1] = bond;
      u = next;
      climb[climbLength++] = u;
    }
  }

  /**
   * The pair of a system that is one simple ring, its lowest atom and first bond, searched as the
   * pairwise search searches it. Returns the pair's charge.
   */
  private long simpleRingPair(int v, int s) {
    int e = systemEdges[s][0], a = bonds.first[e], b = bonds.second[e];
    if (pathSearch == null) pathSearch = new RingPerception.PathSearch(n, neighbors, systems);
    pathSearch.work = 0;
    List<Integer> path1 = new ArrayList<>(), path2 = new ArrayList<>();
    pathSearch.paths(v, a, b, a, b, path1, path2);
    long charge = pathSearch.work;
    pairKind = PAIR_SIMPLE_RING;
    pairBfsLow = pairBfsHigh = charge;
    pairPath = 0;
    pairValid = false;
    if (path1.isEmpty() || path2.isEmpty()) {
      searchWork += charge;
      return charge;
    }
    pairPath = path1.size() + path2.size();
    charge += pairPath;
    if (mark == null) mark = new int[n];
    if (++markStamp == 0) { Arrays.fill(mark, 0); markStamp = 1; }
    for (int i = 1; i < path1.size(); i++) mark[path1.get(i)] = markStamp;
    boolean disjoint = true;
    for (int i = 1; i < path2.size() && disjoint; i++) disjoint = mark[path2.get(i)] != markStamp;
    int size = path1.size() + path2.size() - 1;
    if (disjoint) {
      pairValid = true;
      charge += totalEdges;
      int width = words[s];
      Arrays.fill(key, 0, width, 0L);
      for (int i = 0; i + 1 < path1.size(); i++) flip(localBond[bonds.between(neighbors, path1.get(i), path1.get(i + 1))]);
      for (int i = 0; i + 1 < path2.size(); i++) flip(localBond[bonds.between(neighbors, path2.get(i), path2.get(i + 1))]);
      flip(localBond[e]);
      long h = hash(key, 0, width);
      if (lookup(s, h) < 0) {
        int[] atoms = new int[size];
        int at = 0;
        for (int atom : path1) atoms[at++] = atom;
        for (int i = path2.size() - 1; i >= 1; i--) atoms[at++] = path2.get(i);
        storeOnly(s, h, atoms);
      }
    }
    searchWork += charge + size;
    return charge;
  }

  /**
   * The stored candidate whose key is in {@link #key}, or -1. Keys sharing the hash are compared
   * in the order they were stored.
   */
  private int lookup(int s, long h) {
    int width = words[s], mask = tableEntry.length - 1;
    for (int at = home(h) & mask; ; at = (at + 1) & mask) {
      int entry = tableEntry[at];
      if (entry < 0) return -1;
      if (tableHash[at] != h || candidateSystem[entry] != s) continue;
      searchWork += width;
      long[] stored = keyStore[s];
      int from = candidateKey[entry] * width;
      boolean equal = true;
      for (int i = 0; i < width && equal; i++) equal = stored[from + i] == key[i];
      if (equal) return entry;
    }
  }

  private static int home(long h) {
    return (int) (h ^ (h >>> 32));
  }

  private void store(int s, long h, int[] atoms) {
    storeOnly(s, h, atoms);
    searchWork += atoms.length + words[s];
  }

  private void storeOnly(int s, long h, int[] atoms) {
    int width = words[s], c = candidateAtoms.size();
    candidateAtoms.add(atoms);
    if (c == candidateSystem.length) {
      candidateSystem = Arrays.copyOf(candidateSystem, 2 * c);
      candidateKey = Arrays.copyOf(candidateKey, 2 * c);
      candidateHash = Arrays.copyOf(candidateHash, 2 * c);
    }
    candidateSystem[c] = s;
    candidateHash[c] = h;
    long[] stored = keyStore[s];
    int need = (keyCount[s] + 1) * width;
    if (stored == null || stored.length < need)
      keyStore[s] = stored = stored == null ? new long[Math.max(need, 8 * width)]
                                            : Arrays.copyOf(stored, Math.max(need, 2 * stored.length));
    System.arraycopy(key, 0, stored, keyCount[s] * width, width);
    candidateKey[c] = keyCount[s]++;
    storedKeyWords += width;
    if (2 * (c + 1) > tableEntry.length) {
      int size = 2 * tableEntry.length;
      tableHash = new long[size];
      tableEntry = new int[size];
      Arrays.fill(tableEntry, -1);
      for (int other = 0; other < c; other++) place(other);
    }
    place(c);
  }

  private void place(int c) {
    int mask = tableEntry.length - 1;
    int at = home(candidateHash[c]) & mask;
    while (tableEntry[at] >= 0) at = (at + 1) & mask;
    tableEntry[at] = c;
    tableHash[at] = candidateHash[c];
  }

  /** The cost the pairwise search counts for (start, (eU, eV)), found by that search without its paths. */
  private long replayCost(int start, int eU, int eV, int system) {
    boolean firstReached = start == eU, secondReached = start == eV;
    if (firstReached && secondReached) return 0;
    if (replayVisited == null) {
      replayVisited = new int[n];
      replayQueue = new int[n];
    }
    if (++replayStamp == 0) { Arrays.fill(replayVisited, 0); replayStamp = 1; }
    replayVisited[start] = replayStamp;
    int head = 0, tail = 0;
    replayQueue[tail++] = start;
    long cost = 0;
    while (head < tail && !(firstReached && secondReached)) {
      int u = replayQueue[head++];
      cost += 1 + neighbors[u].length;
      for (int w : neighbors[u]) {
        if ((u == eU && w == eV) || (u == eV && w == eU)) continue;
        if (replayVisited[w] == replayStamp || ofAtom[w] != system) continue;
        replayVisited[w] = replayStamp;
        if (w == eU) firstReached = true;
        if (w == eV) secondReached = true;
        if (firstReached && secondReached) break;
        replayQueue[tail++] = w;
      }
    }
    return cost;
  }

  /** The rings, exactly as the pairwise search finds them. */
  RingPerception.Rings perceive(boolean enumerateRelevant) {
    boolean replayAll = options.replayAllForTesting;
    long low = 0, high = 0, unbounded = 0, summed = 0, index = 0, stoppedAfter = -1;
    long replayedRoots = 0, replayedPairs = 0;
    boolean decided = false, fallback = false;
    int stopRoot = -1, stopSlot = -1;
    pass:
    for (int v = 0; v < n; v++) {
      int s = ofAtom[v];
      if (s < 0) continue;
      if (systems.simpleRing[s]) {
        if (v != lowestAtom[s] || systemEdges[s].length == 0) continue;
        long k = index++;
        long charge = simpleRingPair(v, s);
        classCounts[PAIR_SIMPLE_RING]++;
        record(v, systemEdges[s][0], s);
        if (k + 1 < pairCount) {
          summed++;
          low += charge;
          high += charge;
          if (!replayAll && low > workLimit) {
            decided = fallback = true;
            stopRoot = v;
            stopSlot = 0;
            stoppedAfter = k;
            break;
          }
        }
        continue;
      }
      buildTree(v, s, true);
      int[] list = systemEdges[s];
      for (int i = 0; i < list.length; i++) {
        long k = index++;
        pair(v, s, list[i]);
        classCounts[pairKind]++;
        record(v, list[i], s);
        if (k + 1 >= pairCount) continue;
        long exact = pairPath + (pairValid ? totalEdges : 0);
        summed++;
        low += pairBfsLow + exact;
        high += exact;
        if (pairBfsHigh == UNBOUNDED) unbounded++;
        else high += pairBfsHigh;
        if (pairKind == PAIR_B2) {
          if (b2Count == null) {
            b2Low = new long[n];
            b2High = new long[n];
            b2Count = new int[n];
          }
          b2Low[v] += pairBfsLow;
          if (pairBfsHigh != UNBOUNDED) b2High[v] += pairBfsHigh;
          b2Count[v]++;
        }
        if (!replayAll && low > workLimit) {
          decided = fallback = true;
          stopRoot = v;
          stopSlot = i;
          stoppedAfter = k;
          break pass;
        }
      }
    }
    // Without pairs the pairwise search never tests its work, so it never falls back.
    if (pairCount == 0) decided = true;
    if (!decided && !replayAll) {
      if (unbounded == 0 && high <= workLimit) decided = true;
      else if (low > workLimit) decided = fallback = true;
    }
    if (!decided) {
      // Search the bracketed pairs pairwise, root by root, until the brackets settle the decision.
      for (int v = 0; v < n && !decided && b2Count != null; v++) {
        if (b2Count[v] == 0) continue;
        int s = ofAtom[v];
        // Only the replayed searches are counted, as in the native core, not the tree built again.
        long counted = searchWork;
        buildTree(v, s, false);
        searchWork = counted;
        int[] list = systemEdges[s];
        int end = v == lastRoot ? list.length - 1 : list.length;
        long exact = 0;
        for (int i = 0; i < end; i++) {
          int e = list[i], a = bonds.first[e], b = bonds.second[e];
          if (par[a] != b && par[b] != a) continue;
          int c = par[a] == b ? a : b;
          if (coverTo[c] == c) continue;
          long cost = replayCost(v, a, b, s);
          exact += cost;
          searchWork += cost;
          replayedPairs++;
        }
        replayedRoots++;
        low += exact - b2Low[v];
        high += exact - b2High[v];
        if (options.loosenBracketsForTesting) unbounded -= b2Count[v];
        if (!replayAll) {
          if (unbounded == 0 && high <= workLimit) decided = true;
          else if (low > workLimit) decided = fallback = true;
        }
      }
      if (!decided) {
        decided = true;
        fallback = low > workLimit;
      }
    }
    if (notes != null) {
      notes.pairs = pairCount;
      notes.pairsSummed = summed;
      notes.baseWorkLow = low;
      notes.baseWorkHigh = high;
      notes.baseWorkHighUnbounded = unbounded > 0 || stoppedAfter >= 0;
      notes.baseWorkExact = !notes.baseWorkHighUnbounded && low == high;
      notes.fallback = fallback;
      notes.stoppedAfterPair = stoppedAfter;
      notes.replayedRoots = replayedRoots;
      notes.replayedPairs = replayedPairs;
      notes.nonTreePairs = classCounts[PAIR_NON_TREE];
      notes.classAPairs = classCounts[PAIR_A];
      notes.classB1Pairs = classCounts[PAIR_B1];
      notes.classB1PlusPairs = classCounts[PAIR_B1_PLUS];
      notes.classB2Pairs = classCounts[PAIR_B2];
      notes.simpleRingPairs = classCounts[PAIR_SIMPLE_RING];
    }
    RingPerception.Rings result;
    if (fallback) {
      int[][][] bounded = RingPerception.boundedRings(n, neighbors, systems, cycleRank);
      if (!RingPerception.certifiedMinimumBasis(n, neighbors, bounded[0], cycleRank, workLimit)) {
        RingPerception.Rings refused = new RingPerception.Rings(bounded[0], bounded[1], true);
        if (notes != null) notes.baseBasisTruncated = true;
        result = options.liftRefusedBasis ? lift(refused, stopRoot, stopSlot) : refused;
      } else {
        result = RingPerception.selected(n, neighbors, systems, cycleRank, bounded[0], bounded[1],
                                         enumerateRelevant, true);
      }
    } else {
      int[][][] chosen = select();
      result = RingPerception.selected(n, neighbors, systems, cycleRank, chosen[0], chosen[1],
                                       enumerateRelevant, false);
    }
    if (notes != null) {
      notes.searchWork = searchWork;
      notes.candidates = candidateAtoms.size();
      notes.basisLifted = result.basisLifted;
    }
    return result;
  }

  /**
   * The refused basis completed from the full candidate list, when the search stays within the
   * work limit and the key-word guard and the certificate proves the basis; otherwise the refused
   * rings. Every field but the basis and its flag stays as refused.
   */
  private RingPerception.Rings lift(RingPerception.Rings refused, int stopRoot, int stopSlot) {
    if (stopRoot >= 0) {
      int s = ofAtom[stopRoot];
      if (!systems.simpleRing[s]) {
        buildTree(stopRoot, s, true);
        int[] list = systemEdges[s];
        for (int i = stopSlot + 1; i < list.length; i++) {
          if (overLiftLimits()) return refused;
          pair(stopRoot, s, list[i]);
        }
      }
      for (int v = stopRoot + 1; v < n; v++) {
        s = ofAtom[v];
        if (s < 0) continue;
        if (systems.simpleRing[s]) {
          if (v != lowestAtom[s] || systemEdges[s].length == 0) continue;
          if (overLiftLimits()) return refused;
          simpleRingPair(v, s);
          continue;
        }
        buildTree(v, s, true);
        for (int e : systemEdges[s]) {
          if (overLiftLimits()) return refused;
          pair(v, s, e);
        }
      }
    }
    if (overLiftLimits()) return refused;
    int[][][] chosen = select();
    int[][] basis = chosen[0];
    if (basis.length != cycleRank
        || !RingPerception.certifiedMinimumBasis(n, neighbors, basis, cycleRank, workLimit))
      return refused;
    int[][] numeric = RingPerception.numericBasis(n, neighbors, systems, cycleRank, basis);
    return new RingPerception.Rings(basis, refused.relevant, refused.legacyRelevant, numeric, false,
                                    refused.relevantTruncated, refused.truncated, refused.smallest, true);
  }

  private boolean overLiftLimits() {
    return searchWork > workLimit || storedKeyWords > LIFT_KEY_WORD_LIMIT;
  }

  /**
   * The two-phase selection of the pairwise search over the candidates in stable size order: a
   * candidate independent of strictly shorter rows is relevant, and a relevant one independent of
   * every row joins the basis. Each ring system keeps its own rows over its own bond numbers; once
   * a system has rows for its whole cycle rank, a candidate longer than all of them reduces to
   * zero and is passed over.
   */
  private int[][][] select() {
    int count = candidateAtoms.size(), largest = 0;
    for (int[] atoms : candidateAtoms) largest = Math.max(largest, atoms.length);
    int[] start = new int[largest + 2];
    for (int[] atoms : candidateAtoms) start[atoms.length + 1]++;
    for (int i = 1; i < start.length; i++) start[i] += start[i - 1];
    int[] sorted = new int[count];
    for (int c = 0; c < count; c++) sorted[start[candidateAtoms.get(c).length]++] = c;
    int systemsCount = systemEdges.length;
    long[][] rows = new long[systemsCount][];
    int[][] rowSize = new int[systemsCount][], rowOfPivot = new int[systemsCount][];
    int[] rowCount = new int[systemsCount], longestRow = new int[systemsCount];
    List<int[]> basis = new ArrayList<>(), relevant = new ArrayList<>();
    long[] vector = new long[key.length];
    for (int c : sorted) {
      int s = candidateSystem[c], width = words[s];
      int[] atoms = candidateAtoms.get(c);
      int size = atoms.length;
      if (rowCount[s] == systemRank[s] && size > longestRow[s]) continue;
      searchWork += width;
      if (rows[s] == null) {
        rows[s] = new long[systemRank[s] * width];
        rowSize[s] = new int[systemRank[s]];
        rowOfPivot[s] = new int[systemEdges[s].length];
        Arrays.fill(rowOfPivot[s], -1);
      }
      System.arraycopy(keyStore[s], candidateKey[c] * width, vector, 0, width);
      if (!reduceShorter(vector, width, rows[s], rowSize[s], rowOfPivot[s], size)) continue;
      relevant.add(atoms);
      int pivot = reduceAll(vector, width, rows[s], rowOfPivot[s]);
      if (pivot >= 0 && basis.size() < cycleRank) {
        int row = rowCount[s]++;
        System.arraycopy(vector, 0, rows[s], row * width, width);
        rowSize[s][row] = size;
        rowOfPivot[s][pivot] = row;
        longestRow[s] = size;
        basis.add(atoms);
      }
    }
    return new int[][][] {basis.toArray(new int[0][]), relevant.toArray(new int[0][])};
  }

  /** Reduces by the rows shorter than size; true when something is left. */
  private boolean reduceShorter(long[] vector, int width, long[] rows, int[] rowSize, int[] rowOfPivot, int size) {
    boolean left = false;
    for (int word = 0; word < width; word++) {
      long pending = vector[word];
      while (pending != 0) {
        int bit = Long.numberOfTrailingZeros(pending);
        int row = rowOfPivot[(word << 6) + bit];
        if (row >= 0 && rowSize[row] < size) {
          int at = row * width;
          for (int i = word; i < width; i++) vector[i] ^= rows[at + i];
          searchWork += width;
          pending = bit == 63 ? 0 : vector[word] & (-1L << (bit + 1));
        } else {
          pending &= pending - 1;
        }
      }
      if (vector[word] != 0) left = true;
    }
    return left;
  }

  /** Reduces by every row; returns the lowest bond left, or -1. */
  private int reduceAll(long[] vector, int width, long[] rows, int[] rowOfPivot) {
    int lowest = -1;
    for (int word = 0; word < width; word++) {
      long pending = vector[word];
      while (pending != 0) {
        int bit = Long.numberOfTrailingZeros(pending);
        int position = (word << 6) + bit;
        int row = rowOfPivot[position];
        if (row >= 0) {
          int at = row * width;
          for (int i = word; i < width; i++) vector[i] ^= rows[at + i];
          searchWork += width;
          pending = bit == 63 ? 0 : vector[word] & (-1L << (bit + 1));
        } else {
          if (lowest < 0) lowest = position;
          pending &= pending - 1;
        }
      }
    }
    return lowest;
  }

  /** Records the pair just read when the notes ask for a trace; counts no work. */
  private void record(int v, int e, int s) {
    if (notes == null || notes.trace == null) return;
    PairRecord record = new PairRecord();
    record.root = v;
    record.bond = e;
    record.kind = pairKind;
    record.valid = pairValid;
    record.bfsLow = pairBfsLow;
    record.bfsHigh = pairBfsHigh;
    record.pathLength = pairPath;
    if (pairKind != PAIR_SIMPLE_RING) {
      int a = bonds.first[e], b = bonds.second[e];
      if (par[a] == b || par[b] == a) {
        int c = par[a] == b ? a : b, p = c == a ? b : a;
        int x = coverFrom[c], w = coverTo[c];
        long work = climbWork;
        int length = climbLength;
        int[] saved = Arrays.copyOf(climb, length), savedBonds = Arrays.copyOf(climbBond, length);
        climbTo(w, c, s);
        int[] toC = new int[dep[x] + 1 + climbLength];
        for (int u = x, i = dep[x]; i >= 0; u = par[u], i--) toC[i] = u;
        System.arraycopy(climb, 0, toC, dep[x] + 1, climbLength);
        System.arraycopy(saved, 0, climb, 0, length);
        System.arraycopy(savedBonds, 0, climbBond, 0, length);
        climbLength = length;
        climbWork = work;
        int[] toP = treePath(p);
        record.firstPath = c == a ? toC : toP;
        record.secondPath = c == a ? toP : toC;
      } else {
        record.firstPath = treePath(a);
        record.secondPath = treePath(b);
      }
    }
    notes.trace.add(record);
  }

  private int[] treePath(int u) {
    int[] path = new int[dep[u] + 1];
    for (int x = u, i = dep[u]; i >= 0; x = par[x], i--) path[i] = x;
    return path;
  }

  /**
   * Test hook: the search with no tree-path bond sets kept, so that every key is read along its
   * paths, as for a system too large for them.
   */
  static RingPerception.Rings perceiveWalkingKeysForTesting(int n, int[][] neighbors, boolean enumerateRelevant,
                                                            long workLimit, RingSearchNotes notes) {
    RingPerception.Systems systems = RingPerception.systems(n, neighbors);
    int rank = RingPerception.cycleRank(n, neighbors);
    if (n == 0 || systems == null || rank <= 0) return RingPerception.perceive(n, neighbors, enumerateRelevant, workLimit);
    return new RingCandidateSearch(n, neighbors, systems, rank, workLimit, null, notes, 0).perceive(enumerateRelevant);
  }

  /**
   * Test hook: for every atom of the ring system of root but root itself, the closing half-bond
   * {from, to} the tree of root assigns to it; null elsewhere.
   */
  static int[][] coversForTesting(int n, int[][] neighbors, int root) {
    RingPerception.Systems systems = RingPerception.systems(n, neighbors);
    RingCandidateSearch search = new RingCandidateSearch(n, neighbors, systems, 1, Long.MAX_VALUE, null, null);
    search.buildTree(root, systems.ofAtom[root], false);
    int[][] covers = new int[n][];
    for (int i = 1; i < search.orderSize; i++) {
      int u = search.order[i];
      covers[u] = new int[] {search.coverFrom[u], search.coverTo[u]};
    }
    return covers;
  }
}
