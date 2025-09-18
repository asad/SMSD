/*
 * SPDX-License-Identifier: Apache-2.0
 * © 2025 BioInception PVT LTD.
 */
package com.bioinception.smsd.core;

import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IBond;

import java.util.*;

public final class SearchEngine {

    private SearchEngine() {}

    // -------------------------
    // Time budget
    // -------------------------

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
            if ((++counter & (checkEvery - 1)) != 0) return false;
            return System.nanoTime() > deadlineNanos;
        }

        long remainingMillis() {
            long rem = (deadlineNanos - System.nanoTime()) / 1_000_000L;
            return Math.max(0L, rem);
        }
    }

    // -------------------------
    // Telemetry
    // -------------------------

    public static final class SubstructureStats {
        public final long nodesVisited;
        public final long backtracks;
        public final long candidatesTried;
        public final long prunesAtom;
        public final long prunesBond;
        public final long prunesDegree;
        public final long prunesNLF;
        public final long timeMillis;
        public final boolean timeout;
        public final int solutions;

        SubstructureStats(long nodesVisited, long backtracks, long candidatesTried,
                          long prunesAtom, long prunesBond, long prunesDegree, long prunesNLF,
                          long timeMillis, boolean timeout, int solutions) {
            this.nodesVisited = nodesVisited;
            this.backtracks = backtracks;
            this.candidatesTried = candidatesTried;
            this.prunesAtom = prunesAtom;
            this.prunesBond = prunesBond;
            this.prunesDegree = prunesDegree;
            this.prunesNLF = prunesNLF;
            this.timeMillis = timeMillis;
            this.timeout = timeout;
            this.solutions = solutions;
        }
    }

    public static final class SubstructureResult {
        public final boolean exists;
        public final List<Map<Integer,Integer>> mappings;
        public final SubstructureStats stats;

        SubstructureResult(boolean exists, List<Map<Integer,Integer>> mappings, SubstructureStats stats) {
            this.exists = exists;
            this.mappings = mappings;
            this.stats = stats;
        }
    }

    // -------------------------
    // MCS options
    // -------------------------

    public static final class McsOptions {
        public boolean induced = false;
        public boolean connectedOnly = true;
        public long timeoutMillis = 10_000L;
        public boolean extraSeeds = true;
        public int seedNeighborhoodRadius = 2;
        public int seedMaxAnchors = 12;
        public boolean useTwoHopNLFInExtension = true;
        public boolean useThreeHopNLFInExtension = true; // NEW
    }

    // ------------------------- Public API

    public static boolean isSubstructure(IAtomContainer query, IAtomContainer target, ChemOptions C, long timeoutMs) {
        return isSubstructureWithStats(query, target, C, timeoutMs).exists;
    }

    public static SubstructureResult isSubstructureWithStats(IAtomContainer query, IAtomContainer target,
                                                             ChemOptions C, long timeoutMs) {
        TimeBudget tb = new TimeBudget(timeoutMs);
        Matcher m = makeMatcher(query, target, C, tb);
        long t0 = System.nanoTime();
        boolean ok = m.exists();
        long elapsed = (System.nanoTime() - t0) / 1_000_000L;
        SubstructureStats st = m.buildStats(elapsed, ok ? 1 : 0);
        return new SubstructureResult(ok, Collections.emptyList(), st);
    }

    public static List<Map<Integer,Integer>> findAllSubstructures(IAtomContainer query, IAtomContainer target,
                                                                  ChemOptions C, int maxSolutions, long timeoutMs) {
        return findAllSubstructuresWithStats(query, target, C, maxSolutions, timeoutMs).mappings;
    }

    public static SubstructureResult findAllSubstructuresWithStats(IAtomContainer query, IAtomContainer target,
                                                                   ChemOptions C, int maxSolutions, long timeoutMs) {
        TimeBudget tb = new TimeBudget(timeoutMs);
        Matcher m = makeMatcher(query, target, C, tb);
        long t0 = System.nanoTime();
        List<Map<Integer,Integer>> out = new ArrayList<>();
        m.enumerate(maxSolutions, out);
        long elapsed = (System.nanoTime() - t0) / 1_000_000L;
        SubstructureStats st = m.buildStats(elapsed, out.size());
        return new SubstructureResult(!out.isEmpty(), out, st);
    }

    public static Map<Integer,Integer> findMCS(IAtomContainer m1, IAtomContainer m2,
                                               ChemOptions C, McsOptions M) {
        TimeBudget tb = new TimeBudget(M.timeoutMillis);

        GraphBuilder GB = new GraphBuilder(m1, m2, C, M.induced);
        Map<Integer,Integer> bestSeed = GB.maximumCliqueSeed(tb);
        List<Map<Integer,Integer>> seeds = new ArrayList<>();
        if (!bestSeed.isEmpty()) seeds.add(bestSeed);

        if (M.extraSeeds && !tb.expired()) {
            Map<Integer,Integer> ringSeed = GB.ringAnchorSeed(tb);
            if (!ringSeed.isEmpty()) seeds.add(ringSeed);
        }
        if (M.extraSeeds && !tb.expired()) {
            Map<Integer,Integer> ldSeed = GB.labelDegreeAnchorSeed(tb);
            if (!ldSeed.isEmpty()) seeds.add(ldSeed);
        }
        if (M.extraSeeds && !tb.expired()) {
            long budget = Math.max(1, tb.remainingMillis() / Math.max(1, 4));
            Map<Integer,Integer> ringVFSeed = GB.vf2ppRingSkeletonSeed(tb, budget, 2, 12);
            if (!ringVFSeed.isEmpty()) seeds.add(ringVFSeed);
        }
        if (M.extraSeeds && !tb.expired()) {
            long budget = Math.max(1, tb.remainingMillis() / Math.max(1, 4));
            Map<Integer,Integer> coreVFSeed = GB.vf2ppCoreSeed(tb, budget, 2, 12);
            if (!coreVFSeed.isEmpty()) seeds.add(coreVFSeed);
        }

        Map<Integer,Integer> best = Collections.emptyMap();
        int bestSize = -1;
        long localMillis = Math.max(1, M.timeoutMillis / Math.max(1, seeds.size()));
        for (Map<Integer,Integer> seed : seeds) {
            Map<Integer,Integer> ext = mcGregorExtend(m1, m2, seed, C, tb, localMillis,
                    M.useTwoHopNLFInExtension, M.useThreeHopNLFInExtension);
            ext = (M.induced) ? pruneToInduced(m1, m2, ext, C) : ext;
            ext = (M.connectedOnly) ? largestConnected(m1, m2, ext, C) : ext;
            // NEW: ring-anchor guard to avoid off-ring inflation in cross-family (aromatic vs aliphatic) pairs
            ext = applyRingAnchorGuard(m1, m2, ext);
            if (ext.size() > bestSize) { best = ext; bestSize = ext.size(); }
            if (tb.expired()) break;
        }
        if (bestSize <= 0 && !tb.expired()) {
            Map<Integer,Integer> ext = mcGregorExtend(m1, m2, Collections.emptyMap(), C, tb, tb.remainingMillis(),
                    M.useTwoHopNLFInExtension, M.useThreeHopNLFInExtension);
            ext = (M.induced) ? pruneToInduced(m1, m2, ext, C) : ext;
            ext = (M.connectedOnly) ? largestConnected(m1, m2, ext, C) : ext;
            // Also apply guard here
            best = applyRingAnchorGuard(m1, m2, ext);
        }
        return best;
    }

    private static Matcher makeMatcher(IAtomContainer q, IAtomContainer t, ChemOptions C, TimeBudget tb) {
        String engineName = "VF2PP";
        try { engineName = String.valueOf(C.matcherEngine); } catch (Throwable ignored) {}
        if ("VF2".equalsIgnoreCase(engineName)) return new VF2Matcher(q, t, C, tb);
        if ("VF3".equalsIgnoreCase(engineName)) return new VF2PPMatcher(q, t, C, tb);
        return new VF2PPMatcher(q, t, C, tb);
    }

    // ------------------------- Post-processing helpers

    private static Map<Integer,Integer> pruneToInduced(IAtomContainer m1, IAtomContainer m2,
                                                       Map<Integer,Integer> map, ChemOptions C) {
        if (map.isEmpty()) return map;
        Map<Integer,Integer> M = new LinkedHashMap<>(map);
        boolean changed;
        do {
            changed = false;
            List<Integer> Q = new ArrayList<>(M.keySet());
            outer:
            for (int i = 0; i < Q.size(); i++) {
                for (int j = i + 1; j < Q.size(); j++) {
                    int qi = Q.get(i), qj = Q.get(j);
                    int ti = M.get(qi), tj = M.get(qj);
                    IBond qb = m1.getBond(m1.getAtom(qi), m1.getAtom(qj));
                    IBond tb = m2.getBond(m2.getAtom(ti), m2.getAtom(tj));
                    boolean qHas = qb != null, tHas = tb != null;
                    boolean ok = (qHas == tHas);
                    if (ok && qHas) ok = SearchEngine.ChemOps.bondsCompatible(qb, tb, C);
                    if (!ok) {
                        int di = m1.getConnectedAtomsList(m1.getAtom(qi)).size();
                        int dj = m1.getConnectedAtomsList(m1.getAtom(qj)).size();
                        if (di >= dj) M.remove(qi); else M.remove(qj);
                        changed = true;
                        break outer;
                    }
                }
            }
        } while (changed);
        return M;
    }

    private static Map<Integer,Integer> largestConnected(IAtomContainer m1, IAtomContainer m2,
                                                         Map<Integer,Integer> map, ChemOptions C) {
        if (map.isEmpty()) return map;
        Map<Integer,List<Integer>> adj = new HashMap<>();
        for (int qi : map.keySet()) adj.put(qi, new ArrayList<>());
        for (int qi : map.keySet()) {
            for (int qk : map.keySet()) {
                if (qi >= qk) continue;
                IBond b = m1.getBond(m1.getAtom(qi), m1.getAtom(qk));
                if (b != null) {
                    adj.get(qi).add(qk);
                    adj.get(qk).add(qi);
                }
            }
        }
        Set<Integer> seen = new HashSet<>();
        List<Set<Integer>> comps = new ArrayList<>();
        for (int qi : map.keySet()) {
            if (seen.contains(qi)) continue;
            Set<Integer> comp = new LinkedHashSet<>();
            Deque<Integer> dq = new ArrayDeque<>();
            dq.add(qi); seen.add(qi);
            while (!dq.isEmpty()) {
                int u = dq.pollFirst();
                comp.add(u);
                for (int v : adj.getOrDefault(u, Collections.emptyList())) {
                    if (!seen.contains(v)) { seen.add(v); dq.addLast(v); }
                }
            }
            comps.add(comp);
        }
        Set<Integer> best = comps.stream().max(Comparator.comparingInt(Set::size)).orElse(Collections.emptySet());
        if (best.size() == map.size()) return map;
        Map<Integer,Integer> pruned = new LinkedHashMap<>();
        for (int qi : best) pruned.put(qi, map.get(qi));
        return pruned;
    }

    // NEW: ring-anchor guard to prevent side-chain-only matches when only one side has rings
    private static Map<Integer,Integer> applyRingAnchorGuard(IAtomContainer m1, IAtomContainer m2,
                                                             Map<Integer,Integer> map) {
        if (map.isEmpty()) return map;
        boolean qHasRing = hasRing(m1);
        boolean tHasRing = hasRing(m2);
        if (qHasRing && !tHasRing) {
            int mappedRing = countMappedRingAtoms(m1, map);
            if (mappedRing == 0) return Collections.emptyMap();
        }
        return map;
    }
    private static boolean hasRing(IAtomContainer m) {
        for (int i=0;i<m.getAtomCount();i++) if (m.getAtom(i).isInRing()) return true;
        return false;
    }
    private static int countMappedRingAtoms(IAtomContainer m, Map<Integer,Integer> map) {
        int c=0; for (int qi : map.keySet()) if (m.getAtom(qi).isInRing()) c++; return c;
    }

    // ------------------------- Matcher interface

    private static interface Matcher {
        boolean exists();
        void enumerate(int maxSolutions, List<Map<Integer,Integer>> out);
        SubstructureStats buildStats(long elapsedMillis, int solutions);
    }

    // ------------------------- Common helpers

    private static int labelOf(IAtom a) {
        Integer az = a.getAtomicNumber();
        int z = az == null ? 0 : az.intValue();
        int aromatic = a.isAromatic() ? 1 : 0;
        int ring = a.isInRing() ? 1 : 0;
        return (z << 2) | (aromatic << 1) | ring;
    }

    private static List<Integer> neighbours(IAtomContainer m, int i) {
        List<Integer> out = new ArrayList<>();
        IAtom a = m.getAtom(i);
        for (IAtom b : m.getConnectedAtomsList(a)) out.add(m.getAtomNumber(b));
        return out;
    }

    private static Map<Integer,Integer> buildNLF1(IAtomContainer m, int idx, List<List<Integer>> N) {
        Map<Integer,Integer> freq = new HashMap<>();
        for (int nb : N.get(idx)) {
            IAtom a = m.getAtom(nb);
            int L = labelOf(a);
            freq.put(L, freq.getOrDefault(L, 0) + 1);
        }
        return freq;
    }

    private static Map<Integer,Integer> buildNLF2(IAtomContainer m, int idx, List<List<Integer>> N) {
        Map<Integer,Integer> freq = new HashMap<>();
        BitSet direct = new BitSet(m.getAtomCount());
        for (int nb : N.get(idx)) direct.set(nb);
        Set<Integer> seen = new HashSet<>();
        for (int nb : N.get(idx)) {
            for (IAtom b2 : m.getConnectedAtomsList(m.getAtom(nb))) {
                int j = m.getAtomNumber(b2);
                if (j == idx || direct.get(j)) continue;
                if (!seen.add(j)) continue;
                int L = labelOf(m.getAtom(j));
                freq.put(L, freq.getOrDefault(L, 0) + 1);
            }
        }
        return freq;
    }

    private static Map<Integer,Integer> buildNLF3(IAtomContainer m, int idx, List<List<Integer>> N) {
        Map<Integer,Integer> freq = new HashMap<>();
        int n = m.getAtomCount();
        BitSet level1 = new BitSet(n);
        BitSet level2 = new BitSet(n);
        for (int nb : N.get(idx)) level1.set(nb);
        for (int v = level1.nextSetBit(0); v >= 0; v = level1.nextSetBit(v+1)) {
            for (IAtom b2 : m.getConnectedAtomsList(m.getAtom(v))) {
                int j = m.getAtomNumber(b2);
                if (j == idx || level1.get(j)) continue;
                level2.set(j);
            }
        }
        BitSet level3 = new BitSet(n);
        for (int v = level2.nextSetBit(0); v >= 0; v = level2.nextSetBit(v+1)) {
            for (IAtom b3 : m.getConnectedAtomsList(m.getAtom(v))) {
                int j = m.getAtomNumber(b3);
                if (j == idx || level1.get(j) || level2.get(j)) continue;
                level3.set(j);
            }
        }
        for (int j = level3.nextSetBit(0); j >= 0; j = level3.nextSetBit(j+1)) {
            int L = labelOf(m.getAtom(j));
            freq.put(L, freq.getOrDefault(L, 0) + 1);
        }
        return freq;
    }

    private static boolean nlfOk(Map<Integer,Integer> fq, Map<Integer,Integer> ft) {
        for (Map.Entry<Integer,Integer> e : fq.entrySet()) {
            int L = e.getKey();
            int cq = e.getValue();
            int ct = ft.getOrDefault(L, 0);
            if (ct < cq) return false;
        }
        return true;
    }

    // ------------------------- VF2 (baseline)

    private static final class VF2Matcher implements Matcher {
        final IAtomContainer q, t;
        final ChemOptions C;
        final TimeBudget tb;
        final int Nq, Nt;
        final int[] q2t, t2q;
        final boolean[][] compat;
        final int[] qdeg, tdeg;
        final List<List<Integer>> QN, TN;
        final BitSet[] AQ, AT;

        final Map<Integer,Integer>[] qNLF1, tNLF1;
        final Map<Integer,Integer>[] qNLF2, tNLF2;
        final Map<Integer,Integer>[] qNLF3, tNLF3;

        long nodesVisited=0, backtracks=0, candidatesTried=0;
        long prunesAtom=0, prunesBond=0, prunesDegree=0, prunesNLF=0;
        boolean timedOut=false, found=false;

        @SuppressWarnings("unchecked")
        VF2Matcher(IAtomContainer q, IAtomContainer t, ChemOptions C, TimeBudget tb) {
            this.q=q; this.t=t; this.C=C; this.tb=tb;
            this.Nq=q.getAtomCount(); this.Nt=t.getAtomCount();
            this.q2t=new int[Nq]; this.t2q=new int[Nt];
            Arrays.fill(q2t,-1); Arrays.fill(t2q,-1);

            this.QN=new ArrayList<>(Nq);
            this.TN=new ArrayList<>(Nt);
            for (int i=0;i<Nq;i++) QN.add(neighbours(q,i));
            for (int j=0;j<Nt;j++) TN.add(neighbours(t,j));

            this.qdeg=new int[Nq]; this.tdeg=new int[Nt];
            for (int i=0;i<Nq;i++) qdeg[i]=QN.get(i).size();
            for (int j=0;j<Nt;j++) tdeg[j]=TN.get(j).size();

            this.AQ=new BitSet[Nq]; for (int i=0;i<Nq;i++){ AQ[i]=new BitSet(Nq); for(int k:QN.get(i)) AQ[i].set(k);}
            this.AT=new BitSet[Nt]; for (int j=0;j<Nt;j++){ AT[j]=new BitSet(Nt); for(int k:TN.get(j)) AT[j].set(k);}

            this.compat=new boolean[Nq][Nt];
            this.qNLF1=new Map[Nq]; this.tNLF1=new Map[Nt];
            this.qNLF2=new Map[Nq]; this.tNLF2=new Map[Nt];
            this.qNLF3=new Map[Nq]; this.tNLF3=new Map[Nt];

            for (int i=0;i<Nq;i++) {
                qNLF1[i]=buildNLF1(q,i,QN);
                qNLF2[i]=buildNLF2(q,i,QN);
                qNLF3[i]=buildNLF3(q,i,QN);
                for (int j=0;j<Nt;j++) {
                    boolean ok = ChemOps.atomsCompatible(q.getAtom(i), t.getAtom(j), C);
                    compat[i][j]=ok;
                    if (!ok) prunesAtom++;
                }
            }
            for (int j=0;j<Nt;j++) {
                tNLF1[j]=buildNLF1(t,j,TN);
                tNLF2[j]=buildNLF2(t,j,TN);
                tNLF3[j]=buildNLF3(t,j,TN);
            }
        }

        public SubstructureStats buildStats(long elapsedMillis, int solutions) {
            return new SubstructureStats(nodesVisited, backtracks, candidatesTried,
                    prunesAtom, prunesBond, prunesDegree, prunesNLF,
                    elapsedMillis, timedOut, solutions);
        }

        public boolean exists() {
            if (Nq == 0) return true;
            int si = 0, best=Integer.MAX_VALUE;
            for (int i=0;i<Nq;i++){
                int cnt=0; for(int j=0;j<Nt;j++) if (compat[i][j]) cnt++;
                if (cnt==0) return false;
                if (cnt<best){best=cnt; si=i;}
            }
            boolean[] usedT=new boolean[Nt];
            int[] order=orderQueryAtoms(si);
            backtrack(order,0,usedT);
            return found;
        }

        public void enumerate(int maxSolutions, List<Map<Integer,Integer>> out) {
            if (Nq==0){ out.add(new LinkedHashMap<>()); return;}
            boolean[] usedT=new boolean[Nt];
            int[] order=orderQueryAtoms(selectStart());
            enumerateRec(order,0,usedT,out,maxSolutions);
        }

        private int selectStart(){
            int si=0,best=Integer.MAX_VALUE;
            for (int i=0;i<Nq;i++){
                int cnt=0; for (int j=0;j<Nt;j++) if (compat[i][j]) cnt++;
                if (cnt<best){best=cnt; si=i;}
            }
            return si;
        }

        private int[] orderQueryAtoms(int start){
            boolean[] vis=new boolean[Nq];
            int[] order=new int[Nq];
            int idx=0; Deque<Integer> dq=new ArrayDeque<>();
            dq.add(start); vis[start]=true;
            while(!dq.isEmpty()){
                int u=dq.pollFirst();
                order[idx++]=u;
                for(int v:QN.get(u)) if(!vis[v]){vis[v]=true; dq.addLast(v);}
            }
            for(int i=0;i<Nq;i++) if(!vis[i]) order[idx++]=i;
            return order;
        }

        private void backtrack(int[] order,int pos,boolean[] usedT){
            if (tb.expired()){ timedOut=true; return; }
            nodesVisited++;
            if (found) return;
            if (pos==Nq){ found=true; return; }
            int qi=order[pos];

            List<Integer> cand=new ArrayList<>();
            for(int qk:QN.get(qi)){
                int tk=q2t[qk];
                if (tk!=-1){
                    for(int tj:TN.get(tk)) if(!usedT[tj] && compat[qi][tj] && !cand.contains(tj)) cand.add(tj);
                }
            }
            if (cand.isEmpty()) for(int tj=0;tj<Nt;tj++) if(!usedT[tj] && compat[qi][tj]) cand.add(tj);
            cand.sort((a,b)->Integer.compare(Math.abs(tdeg[a]-qdeg[qi]), Math.abs(tdeg[b]-qdeg[qi])));

            for(int tj: cand){
                candidatesTried++;
                if (!feasible(qi,tj,usedT)) continue;
                q2t[qi]=tj; t2q[tj]=qi; usedT[tj]=true;
                backtrack(order,pos+1,usedT);
                if(found||timedOut) return;
                q2t[qi]=-1; t2q[tj]=-1; usedT[tj]=false;
                backtracks++;
            }
        }

        private void enumerateRec(int[] order,int pos,boolean[] usedT,
                                  List<Map<Integer,Integer>> out,int maxSolutions){
            if (tb.expired()){ timedOut=true; return; }
            nodesVisited++;
            if (out.size()>=maxSolutions) return;
            if (pos==Nq){
                Map<Integer,Integer> map=new LinkedHashMap<>();
                for(int i=0;i<Nq;i++) map.put(i,q2t[i]);
                out.add(map);
                return;
            }
            int qi=order[pos];
            List<Integer> cand=new ArrayList<>();
            for(int qk:QN.get(qi)){
                int tk=q2t[qk];
                if(tk!=-1) for(int tj:TN.get(tk)) if(!usedT[tj] && compat[qi][tj] && !cand.contains(tj)) cand.add(tj);
            }
            if(cand.isEmpty()) for(int tj=0;tj<Nt;tj++) if(!usedT[tj] && compat[qi][tj]) cand.add(tj);
            cand.sort((a,b)->Integer.compare(Math.abs(tdeg[a]-qdeg[qi]), Math.abs(tdeg[b]-qdeg[qi])));
            for(int tj:cand){
                if (tb.expired()){ timedOut=true; return; }
                candidatesTried++;
                if(!feasible(qi,tj,usedT)) continue;
                q2t[qi]=tj; t2q[tj]=qi; usedT[tj]=true;
                enumerateRec(order,pos+1,usedT,out,maxSolutions);
                q2t[qi]=-1; t2q[tj]=-1; usedT[tj]=false;
                backtracks++;
                if(out.size()>=maxSolutions || timedOut) return;
            }
        }

        private boolean feasible(int qi,int tj,boolean[] usedT){
            if (tdeg[tj] < qdeg[qi]) { prunesDegree++; return false; }

            if (!nlfOk(qNLF1[qi], tNLF1[tj])) { prunesNLF++; return false; }
            if (C.useTwoHopNLF && !nlfOk(qNLF2[qi], tNLF2[tj])) { prunesNLF++; return false; }
            if (C.useThreeHopNLF && !nlfOk(qNLF3[qi], tNLF3[tj])) { prunesNLF++; return false; }

            if (C.useBitParallelFeasibility){
                BitSet needed=new BitSet(Nt);
                BitSet qnb=AQ[qi];
                for(int qk=qnb.nextSetBit(0); qk>=0; qk=qnb.nextSetBit(qk+1)){
                    int tk=q2t[qk];
                    if (tk!=-1) needed.set(tk);
                }
                BitSet tmp=(BitSet)needed.clone();
                tmp.andNot(AT[tj]);
                if (!tmp.isEmpty()) { prunesBond++; return false; }
            }

            BitSet qnb=AQ[qi];
            for(int qk=qnb.nextSetBit(0); qk>=0; qk=qnb.nextSetBit(qk+1)){
                int tk=q2t[qk];
                IBond qb=q.getBond(q.getAtom(qi), q.getAtom(qk));
                if (tk!=-1){
                    IBond tb=t.getBond(t.getAtom(tj), t.getAtom(tk));
                    if (qb==null || tb==null) { prunesBond++; return false; }
                    if (!ChemOps.bondsCompatible(qb, tb, C)) { prunesBond++; return false; }
                }
            }
            if (C.ringMatchesRingOnly){
                IAtom ai=q.getAtom(qi), aj=t.getAtom(tj);
                if (ai.isInRing() && !aj.isInRing()) { prunesAtom++; return false; }
            }
            return true;
        }
    }

    // ------------------------- VF2PP (enhanced)

    private static final class VF2PPMatcher implements Matcher {
        final IAtomContainer q, t;
        final ChemOptions C;
        final TimeBudget tb;
        final int Nq, Nt;
        final int[] q2t, t2q;
        final boolean[][] compat;
        final int[] qdeg, tdeg;
        final BitSet[] AQ, AT;
        final List<List<Integer>> QN, TN;

        final BitSet qFrontier = new BitSet();
        final BitSet tFrontier = new BitSet();

        final Map<Integer,Integer>[] qNLF1, tNLF1;
        final Map<Integer,Integer>[] qNLF2, tNLF2;
        final Map<Integer,Integer>[] qNLF3, tNLF3;

        long nodesVisited=0, backtracks=0, candidatesTried=0;
        long prunesAtom=0, prunesBond=0, prunesDegree=0, prunesNLF=0;
        boolean timedOut=false, found=false;

        @SuppressWarnings("unchecked")
        VF2PPMatcher(IAtomContainer q, IAtomContainer t, ChemOptions C, TimeBudget tb) {
            this.q=q; this.t=t; this.C=C; this.tb=tb;
            this.Nq=q.getAtomCount(); this.Nt=t.getAtomCount();
            this.q2t=new int[Nq]; this.t2q=new int[Nt];
            Arrays.fill(q2t,-1); Arrays.fill(t2q,-1);

            this.QN=new ArrayList<>(Nq);
            this.TN=new ArrayList<>(Nt);
            for (int i=0;i<Nq;i++) QN.add(neighbours(q,i));
            for (int j=0;j<Nt;j++) TN.add(neighbours(t,j));

            this.qdeg=new int[Nq]; this.tdeg=new int[Nt];
            for (int i=0;i<Nq;i++) qdeg[i]=QN.get(i).size();
            for (int j=0;j<Nt;j++) tdeg[j]=TN.get(j).size();

            this.AQ=new BitSet[Nq]; for (int i=0;i<Nq;i++){ AQ[i]=new BitSet(Nq); for(int k:QN.get(i)) AQ[i].set(k); }
            this.AT=new BitSet[Nt]; for (int j=0;j<Nt;j++){ AT[j]=new BitSet(Nt); for(int k:TN.get(j)) AT[j].set(k); }

            this.compat=new boolean[Nq][Nt];
            this.qNLF1=new Map[Nq]; this.tNLF1=new Map[Nt];
            this.qNLF2=new Map[Nq]; this.tNLF2=new Map[Nt];
            this.qNLF3=new Map[Nq]; this.tNLF3=new Map[Nt];

            for (int i=0;i<Nq;i++) {
                qNLF1[i]=buildNLF1(q,i,QN);
                qNLF2[i]=buildNLF2(q,i,QN);
                qNLF3[i]=buildNLF3(q,i,QN);
                for (int j=0;j<Nt;j++) {
                    boolean ok = ChemOps.atomsCompatible(q.getAtom(i), t.getAtom(j), C);
                    compat[i][j]=ok;
                    if (!ok) prunesAtom++;
                }
            }
            for (int j=0;j<Nt;j++) {
                tNLF1[j]=buildNLF1(t,j,TN);
                tNLF2[j]=buildNLF2(t,j,TN);
                tNLF3[j]=buildNLF3(t,j,TN);
            }
        }

        public SubstructureStats buildStats(long elapsedMillis, int solutions) {
            return new SubstructureStats(nodesVisited, backtracks, candidatesTried,
                    prunesAtom, prunesBond, prunesDegree, prunesNLF,
                    elapsedMillis, timedOut, solutions);
        }

        public boolean exists() {
            if (Nq == 0) return true;
            int[] order = orderForAllowed();
            boolean[] usedT = new boolean[Nt];
            backtrack(order, 0, usedT);
            return found;
        }

        public void enumerate(int maxSolutions, List<Map<Integer,Integer>> out) {
            if (Nq == 0) { out.add(new LinkedHashMap<>()); return; }
            int[] order = orderForAllowed();
            boolean[] usedT = new boolean[Nt];
            enumerateRec(order, 0, usedT, out, maxSolutions);
        }

        private int[] orderForAllowed() {
            // BFS from a high-degree aromatic anchor when possible.
            int start = -1;
            int bestScore = Integer.MIN_VALUE;
            for (int i=0;i<Nq;i++) {
                int deg = qdeg[i];
                int score = (q.getAtom(i).isAromatic()?1000:0) + (q.getAtom(i).isInRing()?500:0) + deg;
                if (score > bestScore) { bestScore = score; start = i; }
            }
            boolean[] vis=new boolean[Nq];
            List<Integer> ord=new ArrayList<>(Nq);
            Deque<Integer> dq=new ArrayDeque<>();
            dq.add(start); vis[start]=true;
            while(!dq.isEmpty()){
                int u=dq.pollFirst();
                ord.add(u);
                List<Integer> nbs = new ArrayList<>(QN.get(u));
                nbs.sort((a,b)->Integer.compare(qdeg[b], qdeg[a]));
                for(int v:nbs) if(!vis[v]){ vis[v]=true; dq.addLast(v); }
            }
            for(int i=0;i<Nq;i++) if(!vis[i]) ord.add(i);
            int[] order=new int[ord.size()];
            for(int i=0;i<ord.size();i++) order[i]=ord.get(i);
            return order;
        }

        private void addToFrontier(int qi,int tj,boolean[] usedT){
            BitSet qnb=AQ[qi];
            for(int v=qnb.nextSetBit(0); v>=0; v=qnb.nextSetBit(v+1)) if(q2t[v]==-1) qFrontier.set(v);
            BitSet tnb=AT[tj];
            for(int u=tnb.nextSetBit(0); u>=0; u=tnb.nextSetBit(u+1)) if(!usedT[u]) tFrontier.set(u);
        }

        private List<Integer> candidateTargets(int qi, boolean[] usedT){
            List<Integer> cand=new ArrayList<>();
            if (!tFrontier.isEmpty()){
                for(int tj=tFrontier.nextSetBit(0); tj>=0; tj=tFrontier.nextSetBit(tj+1))
                    if(!usedT[tj] && compat[qi][tj]) cand.add(tj);
            } else {
                for(int tj=0;tj<Nt;tj++) if(!usedT[tj] && compat[qi][tj]) cand.add(tj);
            }
            cand.sort((a,b)->Integer.compare(Math.abs(tdeg[a]-qdeg[qi]), Math.abs(tdeg[b]-qdeg[qi])));
            return cand;
        }

        private void backtrack(int[] order,int pos,boolean[] usedT){
            if (tb.expired()){ timedOut=true; return; }
            nodesVisited++;
            if (found) return;
            if (pos==order.length){ found=true; return; }
            int qi=order[pos];

            List<Integer> cand=candidateTargets(qi, usedT);
            for(int tj:cand){
                candidatesTried++;
                if (!feasible(qi,tj,usedT)) continue;
                q2t[qi]=tj; t2q[tj]=qi; usedT[tj]=true;
                BitSet savedQF=(BitSet)qFrontier.clone();
                BitSet savedTF=(BitSet)tFrontier.clone();
                addToFrontier(qi,tj,usedT);
                backtrack(order,pos+1,usedT);
                if(found||timedOut) return;
                q2t[qi]=-1; t2q[tj]=-1; usedT[tj]=false;
                qFrontier.clear(); qFrontier.or(savedQF);
                tFrontier.clear(); tFrontier.or(savedTF);
                backtracks++;
            }
        }

        private void enumerateRec(int[] order,int pos,boolean[] usedT,
                                  List<Map<Integer,Integer>> out,int maxSolutions){
            if (tb.expired()){ timedOut=true; return; }
            nodesVisited++;
            if (out.size()>=maxSolutions) return;
            if (pos==order.length){
                Map<Integer,Integer> map=new LinkedHashMap<>();
                for(int k=0;k<order.length;k++){ int qi=order[k]; map.put(qi, q2t[qi]); }
                out.add(map); return;
            }
            int qi=order[pos];
            List<Integer> cand=candidateTargets(qi, usedT);
            for(int tj:cand){
                if (tb.expired()){ timedOut=true; return; }
                candidatesTried++;
                if (!feasible(qi,tj,usedT)) continue;
                q2t[qi]=tj; t2q[tj]=qi; usedT[tj]=true;
                BitSet savedQF=(BitSet)qFrontier.clone();
                BitSet savedTF=(BitSet)tFrontier.clone();
                addToFrontier(qi,tj,usedT);
                enumerateRec(order,pos+1,usedT,out,maxSolutions);
                q2t[qi]=-1; t2q[tj]=-1; usedT[tj]=false;
                qFrontier.clear(); qFrontier.or(savedQF);
                tFrontier.clear(); tFrontier.or(savedTF);
                backtracks++;
                if(out.size()>=maxSolutions || timedOut) return;
            }
        }

        private boolean feasible(int qi,int tj,boolean[] usedT){
            if (tdeg[tj] < qdeg[qi]) { prunesDegree++; return false; }
            if (!nlfOk(qNLF1[qi], tNLF1[tj])) { prunesNLF++; return false; }
            if (C.useTwoHopNLF && !nlfOk(qNLF2[qi], tNLF2[tj])) { prunesNLF++; return false; }
            if (C.useThreeHopNLF && !nlfOk(qNLF3[qi], tNLF3[tj])) { prunesNLF++; return false; }

            if (C.useBitParallelFeasibility){
                BitSet needed=new BitSet(Nt);
                BitSet qnb=AQ[qi];
                for(int qk=qnb.nextSetBit(0); qk>=0; qk=qnb.nextSetBit(qk+1)){
                    int tk=q2t[qk];
                    if (tk!=-1) needed.set(tk);
                }
                BitSet tmp=(BitSet)needed.clone();
                tmp.andNot(AT[tj]);
                if (!tmp.isEmpty()) { prunesBond++; return false; }
            }

            BitSet qnb=AQ[qi];
            for(int qk=qnb.nextSetBit(0); qk>=0; qk=qnb.nextSetBit(qk+1)){
                int tk=q2t[qk];
                IBond qb=q.getBond(q.getAtom(qi), q.getAtom(qk));
                if (tk!=-1){
                    IBond tb=t.getBond(t.getAtom(tj), t.getAtom(tk));
                    if (qb==null || tb==null) { prunesBond++; return false; }
                    if (!ChemOps.bondsCompatible(qb, tb, C)) { prunesBond++; return false; }
                }
            }
            if (C.ringMatchesRingOnly){
                IAtom ai=q.getAtom(qi), aj=t.getAtom(tj);
                if (ai.isInRing() && !aj.isInRing()) { prunesAtom++; return false; }
            }
            return true;
        }
    }

    // ------------------------- McGregor-like extension with NLF1/2/3 pruning

    private static Map<Integer,Integer> mcGregorExtend(IAtomContainer m1, IAtomContainer m2,
                                                       Map<Integer,Integer> seed,
                                                       ChemOptions C, TimeBudget tb,
                                                       long localMillis,
                                                       boolean useTwoHopNLF, boolean useThreeHopNLF) {
        final long localDeadline = System.nanoTime() + localMillis*1_000_000L;
        Map<Integer,Integer> best = new LinkedHashMap<>(seed);
        Set<Integer> allQ = new HashSet<>();
        for (int i=0;i<m1.getAtomCount();i++) allQ.add(i);

        Map<Integer,Integer>[] qNLF1 = buildAllNLF1(m1);
        Map<Integer,Integer>[] tNLF1 = buildAllNLF1(m2);
        Map<Integer,Integer>[] qNLF2 = useTwoHopNLF ? buildAllNLF2(m1) : null;
        Map<Integer,Integer>[] tNLF2 = useTwoHopNLF ? buildAllNLF2(m2) : null;
        Map<Integer,Integer>[] qNLF3 = useThreeHopNLF ? buildAllNLF3(m1) : null;
        Map<Integer,Integer>[] tNLF3 = useThreeHopNLF ? buildAllNLF3(m2) : null;

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
            if (F.isEmpty()) F.addAll(allQ);

            int bestQi = -1; List<Integer> candT = null;
            for (int qi : F) {
                List<Integer> ts = new ArrayList<>();
                Set<Integer> usedT = new HashSet<>(cur.values());
                for (int tj = 0; tj < m2.getAtomCount(); tj++) {
                    if (usedT.contains(tj)) continue;
                    if (!ChemOps.atomsCompatible(m1.getAtom(qi), m2.getAtom(tj), C)) continue;
                    if (!nlfOk(qNLF1[qi], tNLF1[tj])) continue;
                    if (useTwoHopNLF && !nlfOk(qNLF2[qi], tNLF2[tj])) continue;
                    if (useThreeHopNLF && !nlfOk(qNLF3[qi], tNLF3[tj])) continue;
                    boolean ok = true;
                    for (Map.Entry<Integer,Integer> e : cur.entrySet()) {
                        int qk = e.getKey(), tl = e.getValue();
                        IBond qb = m1.getBond(m1.getAtom(qi), m1.getAtom(qk));
                        IBond tbnd = m2.getBond(m2.getAtom(tj), m2.getAtom(tl));
                        if ((qb==null)!=(tbnd==null)) { ok=false; break; }
                        if (qb!=null && !ChemOps.bondsCompatible(qb, tbnd, C)) { ok=false; break; }
                    }
                    if (ok) ts.add(tj);
                }
                if (ts.isEmpty()) continue;
                if (bestQi==-1 || ts.size()<candT.size()) { bestQi=qi; candT=ts; }
            }
            if (bestQi==-1) continue;
            for (int tj : candT) {
                Map<Integer,Integer> nxt = new LinkedHashMap<>(cur);
                nxt.put(bestQi, tj);
                dq.addLast(nxt);
            }
        }
        return best;
    }

    private static Map<Integer,Integer>[] buildAllNLF1(IAtomContainer m) {
        @SuppressWarnings("unchecked") Map<Integer,Integer>[] nlf = new Map[m.getAtomCount()];
        List<List<Integer>> N = new ArrayList<>(m.getAtomCount());
        for (int i=0;i<m.getAtomCount();i++) N.add(neighbours(m, i));
        for (int i=0;i<m.getAtomCount();i++) nlf[i] = buildNLF1(m, i, N);
        return nlf;
    }
    private static Map<Integer,Integer>[] buildAllNLF2(IAtomContainer m) {
        @SuppressWarnings("unchecked") Map<Integer,Integer>[] nlf2 = new Map[m.getAtomCount()];
        List<List<Integer>> N = new ArrayList<>(m.getAtomCount());
        for (int i=0;i<m.getAtomCount();i++) N.add(neighbours(m, i));
        for (int i=0;i<m.getAtomCount();i++) nlf2[i] = buildNLF2(m, i, N);
        return nlf2;
    }
    private static Map<Integer,Integer>[] buildAllNLF3(IAtomContainer m) {
        @SuppressWarnings("unchecked") Map<Integer,Integer>[] nlf3 = new Map[m.getAtomCount()];
        List<List<Integer>> N = new ArrayList<>(m.getAtomCount());
        for (int i=0;i<m.getAtomCount();i++) N.add(neighbours(m, i));
        for (int i=0;i<m.getAtomCount();i++) nlf3[i] = buildNLF3(m, i, N);
        return nlf3;
    }

    // ------------------------- GraphBuilder and seeds (compact)

    private static final class GraphBuilder {
        private final IAtomContainer m1, m2;
        private final ChemOptions C;
        private final boolean induced;

        GraphBuilder(IAtomContainer m1, IAtomContainer m2, ChemOptions C, boolean induced) {
            this.m1 = m1; this.m2 = m2; this.C = C; this.induced = induced;
        }

        static final class Node { final int qi, tj; Node(int qi,int tj){this.qi=qi; this.tj=tj;} }

        Map<Integer,Integer> maximumCliqueSeed(TimeBudget tb) {
            List<Node> nodes = new ArrayList<>();
            int n1 = m1.getAtomCount(), n2 = m2.getAtomCount();
            for (int i=0;i<n1;i++) {
                for (int j=0;j<n2;j++) {
                    if (ChemOps.atomsCompatible(m1.getAtom(i), m2.getAtom(j), C)) nodes.add(new Node(i,j));
                }
            }
            int N = nodes.size();
            if (N==0) return new LinkedHashMap<>();

            BitSet[] adj = new BitSet[N];
            for (int u=0;u<N;u++) adj[u]=new BitSet(N);
            for (int u=0;u<N;u++) {
                Node nu=nodes.get(u);
                for (int v=u+1; v<N; v++){
                    if (tb.expired()) break;
                    Node nv=nodes.get(v);
                    if (nu.qi==nv.qi || nu.tj==nv.tj) continue;
                    IBond qb=m1.getBond(m1.getAtom(nu.qi), m1.getAtom(nv.qi));
                    IBond tbnd=m2.getBond(m2.getAtom(nu.tj), m2.getAtom(nv.tj));
                    boolean ok;
                    if (qb!=null && tbnd!=null) ok = ChemOps.bondsCompatible(qb, tbnd, C);
                    else if (qb==null && tbnd==null) ok = true;
                    else ok = !induced;
                    if (ok){ adj[u].set(v); adj[v].set(u); }
                }
            }
            List<Integer> order=new ArrayList<>(N);
            for(int i=0;i<N;i++) order.add(i);
            order.sort((a,b)->Integer.compare(adj[b].cardinality(), adj[a].cardinality()));

            List<Integer> best=new ArrayList<>();
            for (int idx=0; idx<order.size(); idx++){
                if (tb.expired()) break;
                int start=order.get(idx);
                List<Integer> clique=new ArrayList<>();
                clique.add(start);
                BitSet cand=(BitSet)adj[start].clone();
                while(!cand.isEmpty()){
                    if (tb.expired()) break;
                    int bestCand=cand.nextSetBit(0);
                    clique.add(bestCand);
                    cand.and(adj[bestCand]);
                }
                if (clique.size()>best.size()) best=clique;
            }
            Map<Integer,Integer> seed=new LinkedHashMap<>();
            for (int u:best){ Node n=nodes.get(u); if (!seed.containsKey(n.qi) && !seed.containsValue(n.tj)) seed.put(n.qi, n.tj); }
            return seed;
        }

        Map<Integer,Integer> ringAnchorSeed(TimeBudget tb) {
            Map<Integer,Integer> seed=new LinkedHashMap<>();
            List<Integer> R1=new ArrayList<>(), R2=new ArrayList<>();
            for (int i=0;i<m1.getAtomCount();i++) if (m1.getAtom(i).isInRing()) R1.add(i);
            for (int j=0;j<m2.getAtomCount();j++) if (m2.getAtom(j).isInRing()) R2.add(j);
            boolean[] used=new boolean[m2.getAtomCount()];
            for (int i : R1){
                int bestJ=-1, bestScore=Integer.MIN_VALUE;
                for (int j : R2){
                    if (used[j]) continue;
                    if (!ChemOps.atomsCompatible(m1.getAtom(i), m2.getAtom(j), C)) continue;
                    int score = (m1.getAtom(i).isAromatic() && m2.getAtom(j).isAromatic()?10:0)
                              + Math.min(m1.getConnectedAtomsList(m1.getAtom(i)).size(),
                                         m2.getConnectedAtomsList(m2.getAtom(j)).size());
                    if (score>bestScore){ bestScore=score; bestJ=j; }
                }
                if (bestJ>=0){ seed.put(i,bestJ); used[bestJ]=true; }
                if (tb.expired()) break;
            }
            return seed;
        }

        Map<Integer,Integer> labelDegreeAnchorSeed(TimeBudget tb) {
            Map<Integer,Integer> seed=new LinkedHashMap<>();
            int n1=m1.getAtomCount(), n2=m2.getAtomCount();
            Map<Integer,Integer> freq=new HashMap<>();
            for (int i=0;i<n1;i++){ int L=labelOf(m1.getAtom(i)); freq.put(L, 1+freq.getOrDefault(L,0)); }
            List<Integer> Q=new ArrayList<>(); for (int i=0;i<n1;i++) Q.add(i);
            Q.sort((a,b)->{
                int la=labelOf(m1.getAtom(a)), lb=labelOf(m1.getAtom(b));
                int ra=freq.get(la), rb=freq.get(lb);
                if (ra!=rb) return Integer.compare(ra, rb);
                return Integer.compare(m1.getConnectedAtomsList(m1.getAtom(b)).size(),
                                       m1.getConnectedAtomsList(m1.getAtom(a)).size());
            });
            boolean[] used=new boolean[n2];
            for (int qi : Q){
                int bestJ=-1, bestScore=Integer.MIN_VALUE;
                for (int tj=0;tj<n2;tj++){
                    if (used[tj]) continue;
                    if (!ChemOps.atomsCompatible(m1.getAtom(qi), m2.getAtom(tj), C)) continue;
                    int score = 10 - Math.abs(m1.getConnectedAtomsList(m1.getAtom(qi)).size()
                                            - m2.getConnectedAtomsList(m2.getAtom(tj)).size());
                    if (score>bestScore){ bestScore=score; bestJ=tj; }
                }
                if (bestJ>=0){ seed.put(qi,bestJ); used[bestJ]=true; }
                if (tb.expired()) break;
            }
            return seed;
        }

        Map<Integer,Integer> vf2ppRingSkeletonSeed(TimeBudget tb, long millis, int radius, int maxAnchors) {
            BitSet qRing = ringMask(m1), tRing = ringMask(m2);
            if (qRing.isEmpty() || tRing.isEmpty() || millis<=0) return Collections.emptyMap();
            return vf2ppNeighborhoodSeed(tb, millis, radius, maxAnchors, qRing, tRing);
        }
        Map<Integer,Integer> vf2ppCoreSeed(TimeBudget tb, long millis, int radius, int maxAnchors) {
            BitSet qCore = heavyAtomCoreMask(m1), tCore = heavyAtomCoreMask(m2);
            if (qCore.isEmpty() || tCore.isEmpty() || millis<=0) return Collections.emptyMap();
            return vf2ppNeighborhoodSeed(tb, millis, radius, maxAnchors, qCore, tCore);
        }

        private Map<Integer,Integer> vf2ppNeighborhoodSeed(TimeBudget tb, long millis, int radius, int maxAnchors,
                                                           BitSet qCand, BitSet tCand) {
            List<int[]> pairs=new ArrayList<>();
            for (int i=qCand.nextSetBit(0); i>=0; i=qCand.nextSetBit(i+1)) {
                for (int j=tCand.nextSetBit(0); j>=0; j=tCand.nextSetBit(j+1)) {
                    if (!ChemOps.atomsCompatible(m1.getAtom(i), m2.getAtom(j), C)) continue;
                    int di=m1.getConnectedAtomsList(m1.getAtom(i)).size();
                    int dj=m2.getConnectedAtomsList(m2.getAtom(j)).size();
                    int score = 100 - Math.abs(di-dj);
                    if (m1.getAtom(i).isAromatic() && m2.getAtom(j).isAromatic()) score += 10;
                    pairs.add(new int[]{i,j,score});
                }
            }
            if (pairs.isEmpty()) return Collections.emptyMap();
            pairs.sort((a,b)->Integer.compare(b[2], a[2]));
            if (pairs.size()>maxAnchors) pairs=pairs.subList(0, maxAnchors);

            Map<Integer,Integer> best=Collections.emptyMap(); int bestSize=0;
            for (int[] p : pairs) {
                if (tb.expired() || millis<=0) break;
                int qi=p[0], tj=p[1];
                BitSet qAllowed = neighborhoodMask(m1, qi, radius, qCand);
                BitSet tAllowed = neighborhoodMask(m2, tj, radius, tCand);
                TimeBudget local = new TimeBudget(Math.min(millis, Math.max(1, tb.remainingMillis()/Math.max(1,pairs.size()))));
                VF2PPMatcher m = new VF2PPMatcher(m1, m2, C, local);
                List<Map<Integer,Integer>> out=new ArrayList<>();
                m.enumerate(1, out);
                if (!out.isEmpty()) {
                    Map<Integer,Integer> map=out.get(0);
                    if (map.size()>bestSize){ best=map; bestSize=map.size(); }
                }
                millis -= Math.min(millis, Math.max(1, tb.remainingMillis()/Math.max(1,pairs.size())));
            }
            return best;
        }

        private static BitSet ringMask(IAtomContainer m){
            BitSet bs=new BitSet(m.getAtomCount());
            for (int i=0;i<m.getAtomCount();i++) if (m.getAtom(i).isInRing()) bs.set(i);
            return bs;
        }
        private static BitSet heavyAtomCoreMask(IAtomContainer m){
            int n=m.getAtomCount();
            BitSet mask=new BitSet(n); mask.set(0,n);
            boolean changed;
            do{
                changed=false;
                for (int i=0;i<n;i++){
                    if (!mask.get(i)) continue;
                    IAtom a=m.getAtom(i);
                    Integer Z=a.getAtomicNumber(); int z=(Z==null?0:Z.intValue());
                    if (a.isInRing() || a.isAromatic()) continue;
                    int deg=0; for (IAtom b : m.getConnectedAtomsList(a)){ int id=m.getAtomNumber(b); if (mask.get(id)) deg++; }
                    if (deg<=1 && z==6){ mask.clear(i); changed=true; }
                }
            } while (changed);
            if (mask.isEmpty()){ for (int i=0;i<n;i++){ Integer Z=m.getAtom(i).getAtomicNumber(); if (Z!=null && Z.intValue()!=1) mask.set(i);} }
            return mask;
        }
        private static BitSet neighborhoodMask(IAtomContainer m, int center, int radius, BitSet base){
            BitSet allowed=new BitSet(m.getAtomCount());
            Deque<int[]> dq=new ArrayDeque<>(); dq.add(new int[]{center,0}); allowed.set(center);
            while(!dq.isEmpty()){
                int[] cur=dq.pollFirst(); int v=cur[0], d=cur[1];
                if (d==radius) continue;
                for (IAtom nb : m.getConnectedAtomsList(m.getAtom(v))){
                    int u=m.getAtomNumber(nb);
                    if (base!=null && !base.get(u)) continue;
                    if (!allowed.get(u)){ allowed.set(u); dq.addLast(new int[]{u,d+1}); }
                }
            }
            return allowed;
        }
    }

    // ------------------------- ChemOps

    static final class ChemOps {
        static boolean atomsCompatible(IAtom q, IAtom t, ChemOptions C) {
            if (C.matchAtomType && q.getAtomicNumber()!=null && t.getAtomicNumber()!=null
                    && !q.getAtomicNumber().equals(t.getAtomicNumber())) return false;

            // Formal charge equality (value-based)
            if (C.matchFormalCharge) {
                Integer qc = q.getFormalCharge();
                Integer tc = t.getFormalCharge();
                if (qc != null || tc != null) {
                    if (!java.util.Objects.equals(qc, tc)) return false;
                }
            }

            // Enforce aromatic parity for ring atoms irrespective of modes (prevents cross-family ring inflation)
            if (q.isInRing() && t.isInRing()) {
                boolean qa = q.isAromatic();
                boolean ta = t.isAromatic();
                if (qa != ta) return false;
            }

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
            // Enforce aromatic parity for ring-to-ring edges irrespective of bond-order mode (fixes clustering cross inflation)
            if (qb.isInRing() && tb.isInRing()) {
                boolean qa = qb.isAromatic();
                boolean ta = tb.isAromatic();
                if (qa != ta) return false;
            }
            if (C.useBondStereo) {
                IBond.Stereo qs = qb.getStereo();
                IBond.Stereo ts = tb.getStereo();
                if (qs != IBond.Stereo.NONE && ts != IBond.Stereo.NONE && qs != ts) return false;
            }
            if (C.ringMatchesRingOnly && qb.isInRing() && !tb.isInRing()) return false;
            if (C.matchBondOrder == ChemOptions.BondOrderMode.ANY) return true;
            if (java.util.Objects.equals(qb.getOrder(), tb.getOrder())) return true;
            if (C.matchBondOrder == ChemOptions.BondOrderMode.LOOSE
                    || C.aromaticityMode == ChemOptions.AromaticityMode.FLEXIBLE) {
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
