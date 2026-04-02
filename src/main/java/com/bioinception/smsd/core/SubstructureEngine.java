/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms. */
package com.bioinception.smsd.core;

import java.util.*;
import org.openscience.cdk.interfaces.IAtomContainer;

/**
 * Substructure matching engine: VF2 and VF2PP matchers.
 * Package-private; accessed via SearchEngine facade.
 *
 * @since 6.5.3 — Domain space caching to avoid rebuilding O(Nq×Nt)
 *   compatibility matrices when the same graph pair is queried 6-18x
 *   per reaction mapping.
 */
final class SubstructureEngine {

  private SubstructureEngine() {}

  // ---- Domain space cache (v6.5.3 perf fix) ----
  // Keyed by the actual graph object identities, not identityHashCode alone.
  // Avoids O(Nq×Nt) atomsCompatFast rebuilds per matcher construction.
  private static final int DOMAIN_CACHE_MAX = 256;
  private static final java.util.concurrent.ConcurrentHashMap<DomainCacheKey, long[][]>
      domainCache = new java.util.concurrent.ConcurrentHashMap<>();

  static final class DomainCacheKey {
    final MolGraph gq, gt;
    final int hash;

    DomainCacheKey(MolGraph gq, MolGraph gt) {
      this.gq = gq;
      this.gt = gt;
      this.hash = 31 * System.identityHashCode(gq) + System.identityHashCode(gt);
    }

    @Override public boolean equals(Object o) {
      if (this == o) return true;
      if (!(o instanceof DomainCacheKey)) return false;
      DomainCacheKey other = (DomainCacheKey) o;
      return gq == other.gq && gt == other.gt;
    }

    @Override public int hashCode() { return hash; }
  }

  /** Clear domain cache (call when molecule objects are recycled). */
  static void clearDomainCache() { domainCache.clear(); }

  static boolean isSubstructure(MolGraph query, MolGraph target, ChemOptions C, long timeoutMs) {
    // Quick rejections before expensive VF2++
    if (query.n == 0) return true;
    if (query.n > target.n) return false;
    if (query == target) return true; // self-match
    // Fingerprint pre-screen: element frequency + degree check
    if (C.matchAtomType) {
      int[] qFreq = new int[120], tFreq = new int[120];
      for (int i = 0; i < query.n; i++) {
        int z = query.atomicNum[i];
        if (z >= 0 && z < 120) qFreq[z]++;
      }
      for (int i = 0; i < target.n; i++) {
        int z = target.atomicNum[i];
        if (z >= 0 && z < 120) tFreq[z]++;
      }
      for (int z = 0; z < 120; z++) {
        if (qFreq[z] > tFreq[z]) return false;
      }
    }
    return isSubstructureWithStats(query, target, C, timeoutMs).exists;
  }

  static SearchEngine.SubstructureResult isSubstructureWithStats(
      MolGraph query, MolGraph target, ChemOptions C, long timeoutMs) {
    SearchEngine.TimeBudget tb = new SearchEngine.TimeBudget(timeoutMs);
    Matcher m = makeMatcher(query, target, C, tb);
    long t0 = System.nanoTime();
    boolean ok = m.exists();
    long elapsed = (System.nanoTime() - t0) / 1_000_000L;
    return new SearchEngine.SubstructureResult(ok, Collections.emptyList(), m.buildStats(elapsed, ok ? 1 : 0));
  }

  static List<Map<Integer, Integer>> findAllSubstructures(
      MolGraph query, MolGraph target, ChemOptions C, int maxSolutions, long timeoutMs) {
    return findAllSubstructuresWithStats(query, target, C, maxSolutions, timeoutMs).mappings;
  }

  static SearchEngine.SubstructureResult findAllSubstructuresWithStats(
      MolGraph query, MolGraph target, ChemOptions C, int maxSolutions, long timeoutMs) {
    SearchEngine.TimeBudget tb = new SearchEngine.TimeBudget(timeoutMs);
    Matcher m = makeMatcher(query, target, C, tb);
    long t0 = System.nanoTime();
    List<Map<Integer, Integer>> out = new ArrayList<>();
    m.enumerate(maxSolutions, out);
    long elapsed = (System.nanoTime() - t0) / 1_000_000L;
    return new SearchEngine.SubstructureResult(!out.isEmpty(), out, m.buildStats(elapsed, out.size()));
  }

  static Matcher makeMatcher(MolGraph gq, MolGraph gt, ChemOptions C, SearchEngine.TimeBudget tb) {
    return C.matcherEngine == ChemOptions.MatcherEngine.VF2
        ? new VF2Matcher(gq, gt, C, tb) : new VF2PPMatcher(gq, gt, C, tb);
  }

  interface Matcher {
    boolean exists();
    void enumerate(int maxSolutions, List<Map<Integer, Integer>> out);
    SearchEngine.SubstructureStats buildStats(long elapsedMillis, int solutions);
  }

  @SuppressWarnings("unchecked")
  abstract static class AbstractVFMatcher implements Matcher {
    final MolGraph gq, gt;
    final ChemOptions C;
    final SearchEngine.TimeBudget tb;
    final int Nq, Nt;
    final int[] q2t, t2q;
    final int[][] qNLF1, tNLF1, qNLF2, tNLF2, qNLF3, tNLF3;
    final boolean useTwoHop, useThreeHop, useBitParallel, useRingOnly, useStereo;
    final boolean bitParallelSufficient;
    final long[][] domain;
    final long[] usedMask;
    final int tWords;
    final boolean singleWord;
    long nodesVisited = 0, backtracks = 0, candidatesTried = 0;
    long prunesAtom = 0, prunesBond = 0, prunesDegree = 0, prunesNLF = 0;
    boolean timedOut = false, found = false;
    boolean pivotEnabled = false; // VF3-style pivot: only enabled during exists()
    int[] qdeg, tdeg;
    final int[][] qNeighborsByDegDesc;
    // Pre-allocated buffers to avoid per-call allocations in selectCandidates/iterateBits
    final long[] availBuf;
    final int[] candBuf;
    // Per-depth candidate buffer stack: avoids Arrays.copyOf on every selectCandidates call
    final int[][] candBufByDepth;

    AbstractVFMatcher(MolGraph gq, MolGraph gt, ChemOptions C, SearchEngine.TimeBudget tb) {
      this.gq = gq; this.gt = gt;
      this.qdeg = gq.degree; this.tdeg = gt.degree;
      this.C = C; this.tb = tb;
      this.Nq = gq.n; this.Nt = gt.n;
      this.tWords = (Nt + 63) >>> 6;
      this.singleWord = (tWords <= 1);
      this.q2t = new int[Nq]; this.t2q = new int[Nt];
      Arrays.fill(q2t, -1); Arrays.fill(t2q, -1);
      this.useTwoHop = C.useTwoHopNLF && Nq > 12 && Nt > 12;
      this.useThreeHop = C.useThreeHopNLF && Nq > 20 && Nt > 20;
      this.useBitParallel = C.useBitParallelFeasibility;
      this.useRingOnly = C.ringMatchesRingOnly; this.useStereo = C.useBondStereo;
      this.bitParallelSufficient = useBitParallel
          && C.matchBondOrder == ChemOptions.BondOrderMode.ANY && !useStereo && !useRingOnly;
      // Use MolGraph's cached neighbor-by-degree sort (v6.5.3 perf fix)
      this.qNeighborsByDegDesc = gq.getNeighborsByDegDesc();
      // Use cached NLF tables from MolGraph (lazy-built, reused across calls on the same graph).
      this.qNLF1 = gq.getNLF1(); this.tNLF1 = gt.getNLF1();
      this.qNLF2 = useTwoHop ? gq.getNLF2() : null; this.tNLF2 = useTwoHop ? gt.getNLF2() : null;
      this.qNLF3 = useThreeHop ? gq.getNLF3() : null; this.tNLF3 = useThreeHop ? gt.getNLF3() : null;
      // Ensure expensive lazy fields are computed only when the options actually need them.
      if (C.tautomerAware) { gq.ensureTautomerClasses(); gt.ensureTautomerClasses(); }
      if (C.ringFusionMode != ChemOptions.RingFusionMode.IGNORE) { gq.ensureRingCounts(); gt.ensureRingCounts(); }
      // Domain space: try cache first (v6.5.3 perf fix).
      // Same graph pair queried 6-18x per reaction gets O(1) domain reuse.
      DomainCacheKey dKey = new DomainCacheKey(gq, gt);
      long[][] cached = domainCache.get(dKey);
      if (cached != null && cached.length == Nq
          && (Nq == 0 || cached[0].length == tWords)) {
        // Deep copy to avoid shared mutation across matchers
        this.domain = new long[Nq][tWords];
        for (int i = 0; i < Nq; i++)
          System.arraycopy(cached[i], 0, this.domain[i], 0, tWords);
      } else {
        this.domain = new long[Nq][tWords];
        if (Nq > 0 && Nt > 50) {
          Map<Integer, List<Integer>> targetByLabel = new HashMap<>();
          for (int j = 0; j < Nt; j++) {
            int key2 = C.matchAtomType ? gt.atomicNum[j] : 0;
            int ringBit = C.ringMatchesRingOnly ? (gt.ring[j] ? 1 : 0) : 0;
            key2 = (key2 << 2) | (gt.aromatic[j] ? 2 : 0) | ringBit;
            targetByLabel.computeIfAbsent(key2, k -> new ArrayList<>()).add(j);
          }
          for (int i = 0; i < Nq; i++) {
            for (Map.Entry<Integer, List<Integer>> entry : targetByLabel.entrySet()) {
              List<Integer> targets = entry.getValue();
              if (C.matchAtomType && gq.atomicNum[i] != gt.atomicNum[targets.get(0)]) { prunesAtom += targets.size(); continue; }
              for (int j : targets) {
                if (atomsCompatFast(gq, i, gt, j, C)) domain[i][j >>> 6] |= 1L << (j & 63);
                else prunesAtom++;
              }
            }
          }
        } else {
          for (int i = 0; i < Nq; i++)
            for (int j = 0; j < Nt; j++) {
              if (atomsCompatFast(gq, i, gt, j, C)) domain[i][j >>> 6] |= 1L << (j & 63);
              else prunesAtom++;
            }
        }
        // Cache for reuse by subsequent matchers on same graph pair
        if (domainCache.size() < DOMAIN_CACHE_MAX) {
          long[][] copy = new long[Nq][tWords];
          for (int i = 0; i < Nq; i++)
            System.arraycopy(this.domain[i], 0, copy[i], 0, tWords);
          domainCache.put(dKey, copy);
        }
      }
      this.usedMask = new long[tWords];
      this.availBuf = new long[tWords];
      this.candBuf = new int[Nt];
      this.candBufByDepth = new int[Nq][Nt];
    }

    static boolean atomsCompatFast(MolGraph gq, int qi, MolGraph gt, int tj, ChemOptions C) {
      // Tautomer-aware: C/N/O can interchange within tautomeric regions,
      // but all OTHER atom-level constraints still apply.
      boolean tautRelax = C.tautomerAware && gq.tautomerClass != null && gt.tautomerClass != null
          && gq.tautomerClass[qi] != -1 && gt.tautomerClass[tj] != -1;
      if (tautRelax) {
        int aq = gq.atomicNum[qi], at = gt.atomicNum[tj];
        if (!((aq == 6 || aq == 7 || aq == 8) && (at == 6 || at == 7 || at == 8)))
          tautRelax = false;
        // Degree guard: genuine tautomers shift at most 1 proton,
        // so degree difference should be <= 1 for same-element matches.
        if (tautRelax && aq == at && Math.abs(gq.degree[qi] - gt.degree[tj]) > 1)
          tautRelax = false;
      }
      if (!tautRelax && C.matchAtomType && gq.atomicNum[qi] != gt.atomicNum[tj]) return false;
      if (C.matchFormalCharge && gq.formalCharge[qi] != gt.formalCharge[tj]) return false;
      if (C.aromaticityMode == ChemOptions.AromaticityMode.STRICT && gq.aromatic[qi] != gt.aromatic[tj]) return false;
      if (C.ringMatchesRingOnly && gq.ring[qi] != gt.ring[tj]) return false;
      if (C.matchIsotope) {
        int qm = gq.massNumber[qi], tm = gt.massNumber[tj];
        if (qm != 0 && tm != 0 && qm != tm) return false;
      }
      if (C.useChirality && gq.tetraChirality != null && gt.tetraChirality != null) {
        int qs = gq.tetraChirality[qi], ts = gt.tetraChirality[tj];
        if (qs != 0 && ts != 0 && qs != ts) return false;
      }
      if (C.ringFusionMode != ChemOptions.RingFusionMode.IGNORE && gq.ring[qi] && gt.ring[tj]) {
        if (C.ringFusionMode == ChemOptions.RingFusionMode.STRICT) {
          if (gq.ringCount[qi] != gt.ringCount[tj]) return false;
        }
      }
      return true;
    }

    static int[] iterateBits(long[] words, int maxBits) {
      int count = 0; for (long w : words) count += Long.bitCount(w);
      int[] result = new int[count]; int idx = 0;
      for (int w = 0; w < words.length; w++) {
        long bits = words[w];
        while (bits != 0) { result[idx++] = (w << 6) | Long.numberOfTrailingZeros(bits); bits &= bits - 1; }
      }
      return result;
    }

    /** Write set bits into pre-allocated buffer, return count. Zero-allocation hot path. */
    int iterateBitsInto(long[] words, int[] buf) {
      int idx = 0;
      for (int w = 0; w < words.length; w++) {
        long bits = words[w];
        while (bits != 0 && idx < buf.length) { buf[idx++] = (w << 6) | Long.numberOfTrailingZeros(bits); bits &= bits - 1; }
      }
      return idx;
    }

    int[] domainCandidates(int qi) {
      for (int w = 0; w < tWords; w++) availBuf[w] = domain[qi][w] & ~usedMask[w];
      return iterateBits(availBuf, Nt);
    }

    boolean greedyProbe(int[] order) {
      int matched = 0;
      for (int qi : order) {
        int bestTj = -1, bestDist = Integer.MAX_VALUE; boolean found2 = false;
        for (int tj : domainCandidates(qi)) {
          int dist = Math.abs(tdeg[tj] - qdeg[qi]);
          if (dist < bestDist && feasible(qi, tj)) { bestTj = tj; bestDist = dist; found2 = true; if (dist == 0) break; }
        }
        if (!found2) {
          for (int k = matched - 1; k >= 0; k--) {
            int qk = order[k], tk = q2t[qk]; q2t[qk] = -1; t2q[tk] = -1;
            usedMask[tk >>> 6] &= ~(1L << (tk & 63));
          }
          return false;
        }
        q2t[qi] = bestTj; t2q[bestTj] = qi; usedMask[bestTj >>> 6] |= 1L << (bestTj & 63); matched++;
      }
      return true;
    }

    int[] fastisoOrder() {
      int[] order = new int[Nq]; boolean[] picked = new boolean[Nq];
      long[][] workDomain = new long[Nq][tWords];
      for (int i = 0; i < Nq; i++) System.arraycopy(domain[i], 0, workDomain[i], 0, tWords);
      long[] reachable = new long[tWords];
      for (int pos = 0; pos < Nq; pos++) {
        int bestQ = -1, bestSize = Integer.MAX_VALUE, bestDeg = -1;
        for (int i = 0; i < Nq; i++) {
          if (picked[i]) continue;
          int size = 0; for (int w = 0; w < tWords; w++) size += Long.bitCount(workDomain[i][w]);
          if (size < bestSize || (size == bestSize && gq.degree[i] > bestDeg)) { bestQ = i; bestSize = size; bestDeg = gq.degree[i]; }
        }
        order[pos] = bestQ; picked[bestQ] = true;
        for (int qk : gq.neighbors[bestQ]) {
          if (picked[qk]) continue;
          // Skip neighbors whose domain is already empty
          boolean hasAny = false;
          for (int w = 0; w < tWords && !hasAny; w++) if (workDomain[qk][w] != 0) hasAny = true;
          if (!hasAny) continue;
          Arrays.fill(reachable, 0, tWords, 0L);
          for (int w = 0; w < tWords; w++) {
            long bits = workDomain[bestQ][w];
            while (bits != 0) {
              int bit = Long.numberOfTrailingZeros(bits); int tj = (w << 6) | bit;
              if (tj < Nt) for (int rw = 0; rw < tWords; rw++) reachable[rw] |= gt.adjLong[tj][rw];
              bits &= bits - 1;
            }
          }
          for (int rw = 0; rw < tWords; rw++) workDomain[qk][rw] &= reachable[rw];
        }
      }
      return order;
    }

    public SearchEngine.SubstructureStats buildStats(long elapsedMillis, int solutions) {
      return new SearchEngine.SubstructureStats(nodesVisited, backtracks, candidatesTried,
          prunesAtom, prunesBond, prunesDegree, prunesNLF, elapsedMillis, timedOut, solutions);
    }

    void sortByDegreeProximity(int[] cands, int qi) {
      for (int i = 1; i < cands.length; i++) {
        int key = cands[i], keyDist = Math.abs(tdeg[key] - qdeg[qi]); int j = i - 1;
        while (j >= 0 && Math.abs(tdeg[cands[j]] - qdeg[qi]) > keyDist) { cands[j + 1] = cands[j]; j--; }
        cands[j + 1] = key;
      }
    }

    /** Sort first {@code len} entries of {@code buf} by degree proximity to query atom qi. */
    void sortByDegreeProximityBuf(int[] buf, int len, int qi) {
      for (int i = 1; i < len; i++) {
        int key = buf[i], keyDist = Math.abs(tdeg[key] - qdeg[qi]); int j = i - 1;
        while (j >= 0 && Math.abs(tdeg[buf[j]] - qdeg[qi]) > keyDist) { buf[j + 1] = buf[j]; j--; }
        buf[j + 1] = key;
      }
    }

    /** VF3-style pivot: prefer candidate with most mapped-neighbor overlap. */
    void applyPivotHeuristic(int[] buf, int len, int qi) {
      if (len <= 1) return;
      // Count mapped neighbors of qi -- pivot only valuable with >= 2
      int mappedNbrs = 0;
      for (int qk : gq.neighbors[qi]) if (q2t[qk] != -1) mappedNbrs++;
      if (mappedNbrs < 2) return;
      // Find candidate with most adjacency overlap to mapped neighbors
      int bestPivot = 0, bestOverlap = 0;
      for (int c = 0; c < len; c++) {
        int tj = buf[c]; int overlap = 0;
        for (int qk : gq.neighbors[qi]) {
          int tk = q2t[qk];
          if (tk != -1 && (gt.adjLong[tj][tk >>> 6] & (1L << (tk & 63))) != 0)
            overlap++;
        }
        if (overlap > bestOverlap) { bestOverlap = overlap; bestPivot = c; }
      }
      if (bestPivot != 0) { int tmp = buf[0]; buf[0] = buf[bestPivot]; buf[bestPivot] = tmp; }
    }

    boolean feasible(int qi, int tj) {
      if (gt.degree[tj] < gq.degree[qi]) { prunesDegree++; return false; }
      // Ring-only check (cheapest boolean test — run early for fast exit)
      if (useRingOnly && gq.ring[qi] != gt.ring[tj]) { prunesAtom++; return false; }
      if (!MolGraph.nlfOk(qNLF1[qi], tNLF1[tj])) { prunesNLF++; return false; }
      if (useTwoHop && !MolGraph.nlfOk(qNLF2[qi], tNLF2[tj])) { prunesNLF++; return false; }
      if (useThreeHop && !MolGraph.nlfOk(qNLF3[qi], tNLF3[tj])) { prunesNLF++; return false; }
      int[] nbSorted = qNeighborsByDegDesc[qi];
      if (useBitParallel) {
        long[] tAdj = gt.adjLong[tj];
        for (int qk : nbSorted) {
          int tk = q2t[qk];
          if (tk != -1 && (tAdj[tk >>> 6] & (1L << (tk & 63))) == 0) { prunesBond++; return false; }
        }
        if (bitParallelSufficient) return true;
      }
      for (int qk : nbSorted) {
        int tk = q2t[qk];
        if (tk != -1 && !MolGraph.ChemOps.bondsCompatible(gq, qi, qk, gt, tj, tk, C)) { prunesBond++; return false; }
      }
      return true;
    }

    int[] vf3LightOrder() {
      Map<Integer, Integer> labelFreq = new HashMap<>();
      for (int i = 0; i < Nq; i++) labelFreq.merge(gq.label[i], 1, Integer::sum);
      Integer[] order = new Integer[Nq];
      for (int i = 0; i < Nq; i++) order[i] = i;
      Arrays.sort(order, (a, b) -> {
        int freqA = labelFreq.getOrDefault(gq.label[a], 0), freqB = labelFreq.getOrDefault(gq.label[b], 0);
        if (freqA != freqB) return Integer.compare(freqA, freqB);
        if (gq.degree[a] != gq.degree[b]) return Integer.compare(gq.degree[b], gq.degree[a]);
        int domA = 0, domB = 0;
        for (int w = 0; w < tWords; w++) { domA += Long.bitCount(domain[a][w]); domB += Long.bitCount(domain[b][w]); }
        return Integer.compare(domA, domB);
      });
      int[] result = new int[Nq]; for (int i = 0; i < Nq; i++) result[i] = order[i]; return result;
    }

    /** Write candidates for query atom qi into candBufByDepth[pos] and return the count. */
    abstract int selectCandidates(int qi, int pos);
    void onMatch(int pos, int tj) {}
    void onUnmatch(int pos, int tj) {}

    void backtrack(int[] order, int pos) {
      if (tb.expired()) { timedOut = true; return; }
      nodesVisited++; if (found) return;
      if (pos == order.length) { found = true; return; }
      int qi = order[pos];
      int candCount = selectCandidates(qi, pos);
      for (int k = 0; k < candCount; k++) {
        int tj = candBufByDepth[pos][k];
        candidatesTried++; if (!feasible(qi, tj)) continue;
        q2t[qi] = tj; t2q[tj] = qi; usedMask[tj >>> 6] |= 1L << (tj & 63); onMatch(pos, tj);
        backtrack(order, pos + 1);
        if (found || timedOut) return;
        q2t[qi] = -1; t2q[tj] = -1; usedMask[tj >>> 6] &= ~(1L << (tj & 63)); onUnmatch(pos, tj); backtracks++;
      }
    }

    void enumerateRec(int[] order, int pos, List<Map<Integer, Integer>> out, int maxSolutions) {
      if (tb.expired()) { timedOut = true; return; }
      nodesVisited++; if (out.size() >= maxSolutions) return;
      if (pos == order.length) {
        Map<Integer, Integer> map = new LinkedHashMap<>(); for (int k : order) map.put(k, q2t[k]); out.add(map); return;
      }
      int qi = order[pos];
      int candCount = selectCandidates(qi, pos);
      for (int k = 0; k < candCount; k++) {
        int tj = candBufByDepth[pos][k];
        if (tb.expired()) { timedOut = true; return; }
        candidatesTried++; if (!feasible(qi, tj)) continue;
        q2t[qi] = tj; t2q[tj] = qi; usedMask[tj >>> 6] |= 1L << (tj & 63); onMatch(pos, tj);
        enumerateRec(order, pos + 1, out, maxSolutions);
        q2t[qi] = -1; t2q[tj] = -1; usedMask[tj >>> 6] &= ~(1L << (tj & 63)); onUnmatch(pos, tj); backtracks++;
        if (out.size() >= maxSolutions || timedOut) return;
      }
    }
  }

  static final class VF2Matcher extends AbstractVFMatcher {
    VF2Matcher(MolGraph gq, MolGraph gt, ChemOptions C, SearchEngine.TimeBudget tb) { super(gq, gt, C, tb); }
    public boolean exists() {
      if (Nq == 0) return true;
      for (int i = 0; i < Nq; i++) { boolean hasCandidate = false; for (int w = 0; w < tWords; w++) if (domain[i][w] != 0) { hasCandidate = true; break; } if (!hasCandidate) return false; }
      int[] order = fastisoOrder();
      if (Nq <= 50 && greedyProbe(order)) { found = true; return true; }
      backtrack(order, 0); return found;
    }
    public void enumerate(int maxSolutions, List<Map<Integer, Integer>> out) {
      if (Nq == 0) { out.add(new LinkedHashMap<>()); return; }
      enumerateRec(fastisoOrder(), 0, out, maxSolutions);
    }
    int selectCandidates(int qi, int pos) {
      int[] buf = candBufByDepth[pos];
      long[] termMask = new long[tWords]; boolean hasTerminal = false;
      for (int qk : gq.neighbors[qi]) { int tk = q2t[qk]; if (tk != -1) { for (int w = 0; w < tWords; w++) termMask[w] |= gt.adjLong[tk][w]; hasTerminal = true; } }
      if (hasTerminal) {
        for (int w = 0; w < tWords; w++) availBuf[w] = domain[qi][w] & ~usedMask[w] & termMask[w];
        int len = iterateBitsInto(availBuf, buf);
        if (len > 0) { sortByDegreeProximityBuf(buf, len, qi); return len; }
      }
      for (int w = 0; w < tWords; w++) availBuf[w] = domain[qi][w] & ~usedMask[w];
      int len = iterateBitsInto(availBuf, buf);
      sortByDegreeProximityBuf(buf, len, qi);
      return len;
    }
  }

  static final class VF2PPMatcher extends AbstractVFMatcher {
    long[] tFrontierMask; private final long[] candAvail; private final long[][] frontierStack;
    VF2PPMatcher(MolGraph gq, MolGraph gt, ChemOptions C, SearchEngine.TimeBudget tb) {
      super(gq, gt, C, tb);
      this.tFrontierMask = new long[tWords]; this.candAvail = new long[tWords]; this.frontierStack = new long[Nq][tWords];
    }
    public boolean exists() {
      if (Nq == 0) return true;
      pivotEnabled = true;
      int[] order = Nq > 30 ? vf3LightOrder() : fastisoOrder();
      if (Nq <= 50 && greedyProbe(order)) { found = true; pivotEnabled = false; return true; }
      Arrays.fill(tFrontierMask, 0); backtrack(order, 0); pivotEnabled = false; return found;
    }
    public void enumerate(int maxSolutions, List<Map<Integer, Integer>> out) {
      if (Nq == 0) { out.add(new LinkedHashMap<>()); return; }
      Arrays.fill(tFrontierMask, 0); enumerateRec(Nq > 30 ? vf3LightOrder() : fastisoOrder(), 0, out, maxSolutions);
    }
    int selectCandidates(int qi, int pos) {
      int[] buf = candBufByDepth[pos];
      boolean hasFrontier = false; for (int w = 0; w < tWords; w++) if (tFrontierMask[w] != 0) { hasFrontier = true; break; }
      // Only apply frontier if qi is connected to the mapped subgraph (prevents
      // false negatives on disconnected queries like pharmacophores [#6].[#8])
      if (hasFrontier) {
        boolean qiConn = false;
        for (int nb : gq.neighbors[qi]) if (q2t[nb] != -1) { qiConn = true; break; }
        if (!qiConn) hasFrontier = false;
      }
      if (hasFrontier) {
        for (int w = 0; w < tWords; w++) candAvail[w] = domain[qi][w] & ~usedMask[w] & tFrontierMask[w];
        boolean hasAny = false; for (int w = 0; w < tWords; w++) if (candAvail[w] != 0) { hasAny = true; break; }
        if (hasAny) {
          int len = iterateBitsInto(candAvail, buf);
          sortByDegreeProximityBuf(buf, len, qi);
          if (pivotEnabled) applyPivotHeuristic(buf, len, qi);
          return len;
        }
      }
      for (int w = 0; w < tWords; w++) candAvail[w] = domain[qi][w] & ~usedMask[w];
      int len = iterateBitsInto(candAvail, buf);
      sortByDegreeProximityBuf(buf, len, qi);
      if (pivotEnabled) applyPivotHeuristic(buf, len, qi);
      return len;
    }
    void onMatch(int pos, int tj) {
      System.arraycopy(tFrontierMask, 0, frontierStack[pos], 0, tWords);
      for (int u : gt.neighbors[tj]) if (t2q[u] == -1) tFrontierMask[u >>> 6] |= 1L << (u & 63);
      tFrontierMask[tj >>> 6] &= ~(1L << (tj & 63));
    }
    void onUnmatch(int pos, int tj) { System.arraycopy(frontierStack[pos], 0, tFrontierMask, 0, tWords); }
  }
}
