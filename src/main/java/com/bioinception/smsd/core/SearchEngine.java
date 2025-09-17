/*
 * SPDX-License-Identifier: Apache-2.0
 * © 2025 BioInception PVT LTD.
 */
package com.bioinception.smsd.core;

import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IBond;

import java.util.*;
import java.util.concurrent.*;
import java.util.stream.Collectors;

public final class SearchEngine {

    private SearchEngine() {}

    public static final class TimeBudget {
        private final long deadlineNanos;
        private final long checkEvery;
        private long counter;

        public TimeBudget(long timeoutMillis) {
            this.deadlineNanos = System.nanoTime() + Math.max(1, timeoutMillis) * 1_000_000L;
            this.checkEvery = 1024L;
            this.counter = 0L;
        }

        public boolean expired() {
            counter++;
            if ((counter & (checkEvery - 1)) == 0L) {
                return System.nanoTime() >= deadlineNanos;
            }
            return false;
        }
    }

    private static final class TimeoutExceeded extends RuntimeException {
        TimeoutExceeded() { super("Time budget exceeded"); }
    }

    public static boolean isSubstructure(final IAtomContainer query,
                                         final IAtomContainer target,
                                         final ChemOptions C,
                                         final long timeoutMillis) {
        final TimeBudget tb = new TimeBudget(timeoutMillis);
        try {
            return vf2Exists(query, target, C, tb);
        } catch (TimeoutExceeded te) {
            return false;
        }
    }

    public static List<Map<Integer, Integer>> findAllSubstructures(final IAtomContainer query,
                                                                   final IAtomContainer target,
                                                                   final ChemOptions C,
                                                                   final int maxMatches,
                                                                   final long timeoutMillis) {
        final TimeBudget tb = new TimeBudget(timeoutMillis);
        try {
            return vf2All(query, target, C, maxMatches, tb);
        } catch (TimeoutExceeded te) {
            return Collections.emptyList();
        }
    }

    private static boolean vf2Exists(IAtomContainer q, IAtomContainer t, ChemOptions C, TimeBudget tb) {
        VF2 st = new VF2(q, t, C, tb);
        return st.searchExists();
    }

    private static List<Map<Integer,Integer>> vf2All(IAtomContainer q, IAtomContainer t, ChemOptions C, int maxMatches, TimeBudget tb) {
        VF2 st = new VF2(q, t, C, tb);
        return st.searchAll(maxMatches);
    }

    private static final class VF2 {
        private final IAtomContainer q, t;
        private final ChemOptions C;
        private final TimeBudget tb;
        private final int Nq, Nt;
        private final int[] q2t, t2q;
        private final List<List<Integer>> QN, TN;
        private final int[] qdeg, tdeg;
        private final boolean[][] compat;
        private int depth;
        private final int[] qterm;
        private final int[] tterm;

        VF2(IAtomContainer q, IAtomContainer t, ChemOptions C, TimeBudget tb) {
            this.q = q; this.t = t; this.C = C; this.tb = tb;
            this.Nq = q.getAtomCount(); this.Nt = t.getAtomCount();
            this.q2t = new int[Nq]; this.t2q = new int[Nt];
            Arrays.fill(q2t, -1); Arrays.fill(t2q, -1);
            this.depth = 0;
            this.QN = new ArrayList<>(Nq); this.TN = new ArrayList<>(Nt);
            for (int i=0;i<Nq;i++) QN.add(neighbours(q, i));
            for (int j=0;j<Nt;j++) TN.add(neighbours(t, j));
            this.qdeg = new int[Nq]; this.tdeg = new int[Nt];
            for (int i=0;i<Nq;i++) qdeg[i] = QN.get(i).size();
            for (int j=0;j<Nt;j++) tdeg[j] = TN.get(j).size();
            this.compat = new boolean[Nq][Nt];
            precomputeCompatibility();
            this.qterm = new int[Nq]; this.tterm = new int[Nt];
        }

        private List<Integer> neighbours(IAtomContainer m, int i) {
            List<Integer> out = new ArrayList<>();
            IAtom a = m.getAtom(i);
            for (IAtom b : m.getConnectedAtomsList(a)) {
                out.add(m.getAtomNumber(b));
            }
            return out;
        }

        private void precomputeCompatibility() {
            for (int i=0;i<Nq;i++) {
                IAtom qa = q.getAtom(i);
                for (int j=0;j<Nt;j++) {
                    IAtom ta = t.getAtom(j);
                    if (ChemOps.atomsCompatible(qa, ta, C) && qdeg[i] <= tdeg[j] + C.degreeSlack) {
                        compat[i][j] = true;
                    }
                }
            }
        }

        boolean searchExists() {
            if (Nq == 0) return true;
            if (Nq > Nt) return false;
            try {
                return recExists();
            } catch (TimeoutExceeded te) {
                throw te;
            }
        }

        List<Map<Integer,Integer>> searchAll(int maxMatches) {
            List<Map<Integer,Integer>> out = new ArrayList<>();
            if (Nq == 0 || Nq > Nt) return out;
            try {
                recAll(maxMatches, out, new HashSet<>());
            } catch (TimeoutExceeded te) {
                // return what we have
            }
            return out;
        }

        private boolean recExists() {
            if (tb.expired()) throw new TimeoutExceeded();
            if (depth == Nq) return true;
            int qi = pickNextQuery(true);
            if (qi < 0) return false;
            for (int tj : candidateTargets(qi)) {
                if (!inducedCompatible(qi, tj)) continue;
                add(qi, tj);
                if (recExists()) return true;
                backtrack(qi, tj);
            }
            return false;
        }

        private void recAll(int maxMatches, List<Map<Integer,Integer>> out, Set<String> seen) {
            if (tb.expired()) throw new TimeoutExceeded();
            if (out.size() >= maxMatches) return;
            if (depth == Nq) {
                Map<Integer,Integer> m = new LinkedHashMap<>();
                for (int i=0;i<Nq;i++) m.put(i, q2t[i]);
                String key = m.values().toString();
                if (!seen.contains(key)) { seen.add(key); out.add(m); }
                return;
            }
            int qi = pickNextQuery(true);
            if (qi < 0) return;
            for (int tj : candidateTargets(qi)) {
                if (!inducedCompatible(qi, tj)) continue;
                add(qi, tj);
                recAll(maxMatches, out, seen);
                if (out.size() >= maxMatches) { backtrack(qi, tj); return; }
                backtrack(qi, tj);
            }
        }

        private boolean inducedCompatible(int qi, int tj) {
            for (int qn : QN.get(qi)) {
                int mapped = q2t[qn];
                if (mapped != -1) {
                    IBond qb = q.getBond(q.getAtom(qi), q.getAtom(qn));
                    IBond tb = t.getBond(t.getAtom(tj), t.getAtom(mapped));
                    if (tb == null || !ChemOps.bondsCompatible(qb, tb, C)) return false;
                }
            }
            return true;
        }

        private void add(int qi, int tj) {
            q2t[qi] = tj; t2q[tj] = qi; depth++;
            for (int n : QN.get(qi)) if (q2t[n]==-1 && qterm[n]==0) qterm[n] = depth;
            for (int n : TN.get(tj)) if (t2q[n]==-1 && tterm[n]==0) tterm[n] = depth;
        }

        private void backtrack(int qi, int tj) {
            for (int n : QN.get(qi)) if (qterm[n]==depth) qterm[n] = 0;
            for (int n : TN.get(tj)) if (tterm[n]==depth) tterm[n] = 0;
            q2t[qi] = -1; t2q[tj] = -1; depth--;
        }

        private Iterable<Integer> candidateTargets(int qi) {
            List<Integer> out = new ArrayList<>();
            for (int tj=0;tj<Nt;tj++) {
                if (!compat[qi][tj] || t2q[tj] != -1) continue;
                int qTerm=0,qNew=0,tTerm=0,tNew=0;
                for (int qn : QN.get(qi)) {
                    if (q2t[qn]==-1) { if (qterm[qn]>0) qTerm++; else qNew++; }
                }
                for (int tn : TN.get(tj)) {
                    if (t2q[tn]==-1) { if (tterm[tn]>0) tTerm++; else tNew++; }
                }
                if (qTerm <= tTerm && qNew <= tNew) out.add(tj);
            }
            out.sort((a,b) -> {
                IAtom A = t.getAtom(a), B = t.getAtom(b);
                int dA = tdeg[a], dB = tdeg[b];
                int inRingA = A.isInRing() ? 1 : 0, inRingB = B.isInRing() ? 1 : 0;
                int aromA = A.isAromatic() ? 1 : 0, aromB = B.isAromatic() ? 1 : 0;
                if (dA != dB) return Integer.compare(dB, dA);
                if (inRingA != inRingB) return Integer.compare(inRingB, inRingA);
                if (aromA != aromB) return Integer.compare(aromB, aromA);
                return Integer.compare(a, b);
            });
            return out;
        }

        private int pickNextQuery(boolean connectedOnly) {
            boolean hasTerminals = false;
            if (connectedOnly && depth > 0) {
                for (int i = 0; i < Nq; i++) {
                    if (q2t[i] == -1 && qterm[i] > 0) {
                        hasTerminals = true;
                        break;
                    }
                }
            }

            int best = -1;
            int bestCount = Integer.MAX_VALUE;
            int bestDeg = -1, bestRing = -1, bestArom = -1;
            for (int i = 0; i < Nq; i++) {
                if (q2t[i] != -1) continue;
                if (connectedOnly && depth > 0 && hasTerminals && qterm[i] == 0) continue;

                int cnt = 0;
                for (int tj = 0; tj < Nt; tj++)
                    if (compat[i][tj] && t2q[tj] == -1) cnt++;
                if (cnt == 0) return -1;
                int deg = qdeg[i];
                int ring = q.getAtom(i).isInRing() ? 1 : 0;
                int arom = q.getAtom(i).isAromatic() ? 1 : 0;
                if (cnt < bestCount ||
                        (cnt == bestCount && (deg > bestDeg ||
                                (deg == bestDeg && (ring > bestRing ||
                                        (ring == bestRing && (arom > bestArom ||
                                                (arom == bestArom && i < best)))))))) {
                    best = i;
                    bestCount = cnt;
                    bestDeg = deg;
                    bestRing = ring;
                    bestArom = arom;
                }
            }
            return best;
        }
    }

    public static final class McsOptions {
        public boolean connectedOnly = true;
        public boolean induced = false;
        public boolean maximizeBonds = false;
        public long timeoutMillis = 10_000L;
        public int cps = Runtime.getRuntime().availableProcessors();
    }

    public static Map<Integer, Integer> findMCS(final IAtomContainer m1,
                                                final IAtomContainer m2,
                                                final ChemOptions C,
                                                final McsOptions opt) {
        final TimeBudget tb = new TimeBudget(opt.timeoutMillis);
        try {
            return mcsImpl(m1, m2, C, opt, tb);
        } catch (TimeoutExceeded te) {
            return Collections.emptyMap();
        }
    }

    private static Map<Integer,Integer> mcsImpl(IAtomContainer m1, IAtomContainer m2, ChemOptions C, McsOptions opt, TimeBudget tb) {
        final List<Node> nodes = buildProduct(m1, m2, C);
        if (nodes.isEmpty()) return Collections.emptyMap();

        final int n = nodes.size();
        final BitSet[] adjEdge = new BitSet[n];
        final BitSet[] adjFull = new BitSet[n];
        for (int i = 0; i < n; i++) {
            adjEdge[i] = new BitSet(n);
            adjFull[i] = new BitSet(n);
        }

        for (int u=0; u<n; u++) {
            final Node a = nodes.get(u);
            for (int v=u+1; v<n; v++) {
                final Node b = nodes.get(v);
                if (a.q == b.q || a.t == b.t) continue;

                final IBond qb = m1.getBond(m1.getAtom(a.q), m1.getAtom(b.q));
                final IBond tbnd = m2.getBond(m2.getAtom(a.t), m2.getAtom(b.t));

                boolean bothBond = qb != null && tbnd != null && ChemOps.bondsCompatible(qb, tbnd, C);
                boolean bothNone = qb == null && tbnd == null;
                if (bothBond) {
                    adjEdge[u].set(v); adjEdge[v].set(u);
                    adjFull[u].set(v); adjFull[v].set(u);
                } else if (bothNone) {
                    adjFull[u].set(v); adjFull[v].set(u);
                }
            }
        }

        final BitSet[] ADJ = (opt.induced ? adjFull : adjEdge);
        bestMask = new BitSet(n);
        bestSize = 0;

        BitSet P = new BitSet(n);
        P.set(0, n);
        bbmc(new BitSet(n), P, new BitSet(n), ADJ, tb);

        Map<Integer,Integer> mapping = maskToMapping(nodes, bestMask);
        if (opt.connectedOnly) {
            mapping = largestConnectedInQuery(m1, mapping);
        }
        if (opt.maximizeBonds || C.matchBondOrder == ChemOptions.BondOrderMode.LOOSE || C.matchBondOrder == ChemOptions.BondOrderMode.ANY) {
            mapping = mcGregorExtend(m1, m2, mapping, C, tb, 2_000L);
            if (opt.connectedOnly) mapping = largestConnectedInQuery(m1, mapping);
        }
        return mapping;
    }

    private static final class Node {
        final int q, t;
        Node(int q, int t) { this.q=q; this.t=t; }
    }

    private static List<Node> buildProduct(IAtomContainer m1, IAtomContainer m2, ChemOptions C) {
        final List<Node> nodes = new ArrayList<>();
        for (int i=0;i<m1.getAtomCount();i++) {
            for (int j=0;j<m2.getAtomCount();j++) {
                if (ChemOps.atomsCompatible(m1.getAtom(i), m2.getAtom(j), C)) {
                    nodes.add(new Node(i, j));
                }
            }
        }
        return nodes;
    }

    private static Map<Integer,Integer> maskToMapping(List<Node> nodes, BitSet mask) {
        Map<Integer,Integer> m = new LinkedHashMap<>();
        for (int i = mask.nextSetBit(0); i >= 0; i = mask.nextSetBit(i + 1)) {
            Node n = nodes.get(i);
            m.put(n.q, n.t);
        }
        return m;
    }

    private static Map<Integer,Integer> largestConnectedInQuery(IAtomContainer q, Map<Integer,Integer> mapping) {
        if (mapping.isEmpty()) return mapping;
        List<Integer> keys = new ArrayList<>(mapping.keySet());
        Map<Integer,Integer> pos = new HashMap<>();
        for (int i=0;i<keys.size();i++) pos.put(keys.get(i), i);
        int n = keys.size();
        int[] parent = new int[n];
        for (int i=0;i<n;i++) parent[i]=i;
        java.util.function.IntUnaryOperator find = new java.util.function.IntUnaryOperator() {
            @Override public int applyAsInt(int x) {
                int r=x;
                while (parent[r]!=r) r=parent[r];
                while (parent[x]!=x) { int y=parent[x]; parent[x]=r; x=y; }
                return r;
            }
        };
        java.util.function.BiConsumer<Integer,Integer> union = (a,b) -> {
            int ra=find.applyAsInt(a), rb=find.applyAsInt(b);
            if (ra!=rb) parent[rb]=ra;
        };
        for (int i=0;i<keys.size();i++) {
            int qi = keys.get(i);
            for (IAtom nAtom : q.getConnectedAtomsList(q.getAtom(qi))) {
                int qj = q.getAtomNumber(nAtom);
                Integer pj = pos.get(qj);
                if (pj != null) union.accept(i, pj);
            }
        }
        int bestRoot=-1, bestCount=-1;
        Map<Integer,Integer> cnt = new HashMap<>();
        for (int i=0;i<n;i++) {
            int r = find.applyAsInt(i);
            int c = cnt.getOrDefault(r, 0)+1;
            cnt.put(r, c);
            if (c>bestCount) { bestCount=c; bestRoot=r; }
        }
        Set<Integer> keepIdx = new HashSet<>();
        for (int i=0;i<n;i++) if (find.applyAsInt(i)==bestRoot) keepIdx.add(i);
        Map<Integer,Integer> out = new LinkedHashMap<>();
        for (int i=0;i<n;i++) if (keepIdx.contains(i)) out.put(keys.get(i), mapping.get(keys.get(i)));
        return out;
    }

    private static BitSet bestMask;
    private static int bestSize = 0;

    private static void bbmc(BitSet R, BitSet P, BitSet X, BitSet[] ADJ, TimeBudget tb) {
        if (tb.expired()) throw new TimeoutExceeded();
        if (P.isEmpty() && X.isEmpty()) {
            int size = R.cardinality();
            if (size > bestSize) {
                bestSize = size;
                bestMask = (BitSet) R.clone();
            }
            return;
        }
        int bound = R.cardinality() + greedyColouringBound(P, ADJ);
        if (bound <= bestSize) return;

        BitSet UX = (BitSet) P.clone();
        UX.or(X);
        int u = -1;
        int bestDeg = -1;
        for (int i = UX.nextSetBit(0); i >= 0; i = UX.nextSetBit(i + 1)) {
            BitSet pIntersect = (BitSet) P.clone();
            pIntersect.and(ADJ[i]);
            int deg = pIntersect.cardinality();
            if (deg > bestDeg) {
                bestDeg = deg;
                u = i;
            }
        }
        BitSet candidates = (BitSet) P.clone();
        if (u != -1) {
            candidates.andNot(ADJ[u]);
        }

        for (int idx = candidates.nextSetBit(0); idx >= 0; idx = candidates.nextSetBit(idx + 1)) {
            BitSet nextR = (BitSet) R.clone();
            nextR.set(idx);
            BitSet nextP = (BitSet) P.clone();
            nextP.and(ADJ[idx]);
            BitSet nextX = (BitSet) X.clone();
            nextX.and(ADJ[idx]);
            bbmc(nextR, nextP, nextX, ADJ, tb);
            P.clear(idx);
            X.set(idx);
        }
    }

    private static int greedyColouringBound(BitSet P, BitSet[] ADJ) {
        int colours = 0;
        BitSet remaining = (BitSet) P.clone();
        while (!remaining.isEmpty()) {
            colours++;
            BitSet curr = new BitSet(ADJ.length);
            BitSet tmp = (BitSet) remaining.clone();
            for (int i = tmp.nextSetBit(0); i >= 0; i = tmp.nextSetBit(i + 1)) {
                BitSet currNeighbours = (BitSet) curr.clone();
                currNeighbours.and(ADJ[i]);
                if (currNeighbours.isEmpty()) {
                    curr.set(i);
                    remaining.clear(i);
                }
            }
        }
        return colours;
    }

    private static Map<Integer,Integer> mcGregorExtend(IAtomContainer m1, IAtomContainer m2,
                                                       Map<Integer,Integer> seed,
                                                       ChemOptions C, TimeBudget tb,
                                                       long localMillis) {
        final long localDeadline = System.nanoTime() + localMillis*1_000_000L;
        Map<Integer,Integer> best = new LinkedHashMap<>(seed);
        Set<Integer> qUn = new HashSet<>();
        Set<Integer> tUn = new HashSet<>();
        for (int i=0;i<m1.getAtomCount();i++) if (!seed.containsKey(i)) qUn.add(i);
        boolean[] tTaken = new boolean[m2.getAtomCount()];
        for (int j : seed.values()) tTaken[j] = true;
        for (int j=0;j<m2.getAtomCount();j++) if (!tTaken[j]) tUn.add(j);

        java.util.function.BiPredicate<Integer,Integer> feasible = (qi,tj) -> {
            if (!ChemOps.atomsCompatible(m1.getAtom(qi), m2.getAtom(tj), C)) return false;
            for (Map.Entry<Integer,Integer> e : seed.entrySet()) {
                int qk = e.getKey(), tl = e.getValue();
                IBond qb = m1.getBond(m1.getAtom(qi), m1.getAtom(qk));
                IBond tbnd = m2.getBond(m2.getAtom(tj), m2.getAtom(tl));
                if ((qb==null)!=(tbnd==null)) return false;
                if (qb!=null && !ChemOps.bondsCompatible(qb, tbnd, C)) return false;
            }
            return true;
        };

        Deque<Map<Integer,Integer>> dq = new ArrayDeque<>();
        dq.add(new LinkedHashMap<>(seed));
        while (!dq.isEmpty()) {
            if (System.nanoTime() >= localDeadline || tb.expired()) break;
            Map<Integer,Integer> cur = dq.pollFirst();
            if (cur.size() > best.size()) best = new LinkedHashMap<>(cur);
            
            Set<Integer> F = new HashSet<>();
            for (int qk : cur.keySet()) {
                for (IAtom nn : m1.getConnectedAtomsList(m1.getAtom(qk))) {
                    int qn = m1.getAtomNumber(nn);
                    if (!cur.containsKey(qn)) F.add(qn);
                }
            }
            if (F.isEmpty()) F.addAll(qUn);
            
            int bestQi = -1; List<Integer> candT = null;
            for (int qi : F) {
                List<Integer> ts = new ArrayList<>();
                for (int tj : tUn) if (feasible.test(qi, tj)) ts.add(tj);
                if (ts.isEmpty()) continue;
                if (bestQi == -1 || ts.size() < candT.size()) { bestQi = qi; candT = ts; }
            }
            if (bestQi == -1) continue;
            for (int tj : candT) {
                Map<Integer,Integer> nxt = new LinkedHashMap<>(cur);
                nxt.put(bestQi, tj);
                dq.addLast(nxt);
            }
        }
        return best;
    }

    private static final class ChemOps {

        static boolean atomsCompatible(IAtom q, IAtom t, ChemOptions C) {
            if (C.matchAtomType && q.getAtomicNumber() != 0 && q.getAtomicNumber() != t.getAtomicNumber())
                return false;
            if (C.matchFormalCharge && (q.getFormalCharge() != t.getFormalCharge()))
                return false;
            if (C.aromaticityMode == ChemOptions.AromaticityMode.STRICT && (q.isAromatic() != t.isAromatic()))
                return false;
            if (C.ringMatchesRingOnly && q.isInRing() && !t.isInRing())
                return false;
            if (C.useChirality) {
                Integer qp = q.getStereoParity();
                Integer tp = t.getStereoParity();
                if (qp != null && tp != null && qp != 0 && tp != 0 && !qp.equals(tp))
                    return false;
            }
            return true;
        }

        static boolean bondsCompatible(IBond qb, IBond tb, ChemOptions C) {
            if (C.useBondStereo) {
                IBond.Stereo qs = qb.getStereo();
                IBond.Stereo ts = tb.getStereo();
                if (qs != IBond.Stereo.NONE && ts != IBond.Stereo.NONE && qs != ts) {
                    return false;
                }
            }

            if (C.matchBondOrder == ChemOptions.BondOrderMode.ANY) {
                return true;
            }

            if (java.util.Objects.equals(qb.getOrder(), tb.getOrder())) {
                return true;
            }

            if (C.matchBondOrder == ChemOptions.BondOrderMode.LOOSE) {
                boolean qa = qb.isAromatic();
                boolean ta = tb.isAromatic();
                if ((qa && ta) ||
                    (qa && (tb.getOrder() == IBond.Order.SINGLE || tb.getOrder() == IBond.Order.DOUBLE)) ||
                    (ta && (qb.getOrder() == IBond.Order.SINGLE || qb.getOrder() == IBond.Order.DOUBLE))) {
                    return true;
                }
            }
            return false;
        }
    }
}