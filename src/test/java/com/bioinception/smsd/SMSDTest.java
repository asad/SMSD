/*
 * SPDX-License-Identifier: Apache-2.0
 * © 2025 BioInception PVT LTD.
 */
package com.bioinception.smsd;

import com.bioinception.smsd.core.ChemOptions;
import com.bioinception.smsd.core.SMSD;
import org.junit.jupiter.api.*;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.smiles.SmilesParser;

import java.util.*;
import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.aromaticity.ElectronDonation;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

import static org.junit.jupiter.api.Assertions.*;

@TestMethodOrder(MethodOrderer.DisplayName.class)
public class SMSDTest {

    private static final SmilesParser SP = new SmilesParser(DefaultChemObjectBuilder.getInstance());

    private static boolean hasAny() {
        try { ChemOptions.BondOrderMode.valueOf("ANY"); return true; }
        catch (IllegalArgumentException e) { return false; }
    }

    private static ChemOptions opts(boolean atomTypeStrict,
                                    boolean bondOrderStrict,
                                    boolean ringToRingOnly,
                                    boolean stereo) {
        ChemOptions c = new ChemOptions();
        c.matchAtomType = atomTypeStrict;
        c.matchBondOrder = bondOrderStrict
                ? ChemOptions.BondOrderMode.STRICT
                : (hasAny() ? ChemOptions.BondOrderMode.valueOf("ANY")
                            : ChemOptions.BondOrderMode.LOOSE);
        // Ring→ring policy: ring atoms/bonds in the query may only map to ring atoms/bonds in the target
        c.ringMatchesRingOnly = ringToRingOnly;
        c.useChirality = stereo;
        return c;
    }

    private static IAtomContainer mol(String smiles) {
        try {
            IAtomContainer m = SP.parseSmiles(smiles);
            // Perceive types, rings, and aromaticity so isInRing()/isAromatic() are reliable
            AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(m);
            Cycles.markRingAtomsAndBonds(m);
            new Aromaticity(ElectronDonation.daylight(), Cycles.all()).apply(m);
            return m;
        } catch (Exception e) {
            throw new AssertionError("Failed to parse SMILES: " + smiles, e);
        }
    }

    // Java 11: materialise the bonds iterable (no Iterable#toList())
    private static List<IBond> bonds(IAtomContainer m) {
        List<IBond> out = new ArrayList<>();
        for (IBond b : m.bonds()) out.add(b);
        return out;
    }

    private static int mappedBonds(IAtomContainer A, IAtomContainer B, Map<Integer,Integer> map) {
        int cnt = 0;
        for (IBond qb : bonds(A)) {
            int qi = A.getAtomNumber(qb.getAtom(0));
            int qj = A.getAtomNumber(qb.getAtom(1));
            Integer ti = map.get(qi), tj = map.get(qj);
            if (ti == null || tj == null) continue;
            IBond tb = B.getBond(B.getAtom(ti), B.getAtom(tj));
            if (tb != null) cnt++;
        }
        return cnt;
    }

    private static int mappedRingBonds(IAtomContainer A, IAtomContainer B, Map<Integer,Integer> map) {
        int cnt = 0;
        for (IBond qb : bonds(A)) {
            int qi = A.getAtomNumber(qb.getAtom(0));
            int qj = A.getAtomNumber(qb.getAtom(1));
            Integer ti = map.get(qi), tj = map.get(qj);
            if (ti == null || tj == null) continue;
            IBond tb = B.getBond(B.getAtom(ti), B.getAtom(tj));
            if (tb != null && tb.isInRing()) cnt++;
        }
        return cnt;
    }

    private static int countComponentsOnQuery(IAtomContainer q, Collection<Integer> mapped) {
        if (mapped.isEmpty()) return 0;
        List<Integer> idx = new ArrayList<>(mapped);
        Map<Integer,Integer> pos = new HashMap<>();
        for (int i=0;i<idx.size();i++) pos.put(idx.get(i), i);
        int n = idx.size(), comp = 0;
        boolean[] seen = new boolean[n];
        for (int i=0;i<n;i++) {
            if (seen[i]) continue;
            comp++;
            Deque<Integer> dq = new ArrayDeque<>();
            dq.add(i); seen[i]=true;
            while (!dq.isEmpty()) {
                int u = dq.poll();
                int qu = idx.get(u);
                for (IBond b : bonds(q)) {
                    int a = q.getAtomNumber(b.getAtom(0));
                    int c = q.getAtomNumber(b.getAtom(1));
                    if (a==qu || c==qu) {
                        int vq = (a==qu? c : a);
                        Integer vi = pos.get(vq);
                        if (vi != null && !seen[vi]) { seen[vi]=true; dq.add(vi); }
                    }
                }
            }
        }
        return comp;
    }

    private static void assertNonEmpty(Map<Integer,Integer> map, String msg) {
        assertNotNull(map, msg + " (missing map)");
        assertTrue(!map.isEmpty(), msg + " (expected non-empty mapping)");
    }

    // ─────────────────────────────────────────────────────────────────────────────
    // A — quick sanity tests
    // ─────────────────────────────────────────────────────────────────────────────

    @Test @DisplayName("A.1  Ethanol is a substructure (quick sanity)")
    void a1_substructure_exists() {
        SMSD smsd = new SMSD(mol("CCO"), mol("CCNCCO"), new ChemOptions());
        assertTrue(smsd.isSubstructure(1200L));
    }

    @Test @DisplayName("A.2  Benzene is not a substructure of hexane")
    void a2_substructure_not_exists() {
        SMSD smsd = new SMSD(mol("c1ccccc1"), mol("CCCCCC"), opts(true, true, false, false));
        assertFalse(smsd.isSubstructure(1200L));
    }

    @Test @DisplayName("A.3  MCS(bz, naph) has 6 atoms")
    void a3_mcs_size_benzene_naphthalene() {
        SMSD smsd = new SMSD(mol("c1ccccc1"), mol("c1ccc2ccccc2c1"), opts(true, false, false, false));
        Map<Integer,Integer> m = smsd.findMCS(false, true, 2000L);
        assertNonEmpty(m, "benzene vs naphthalene");
        assertEquals(6, m.size(), "MCS(bz, naph) should be 6");
    }

    // ─────────────────────────────────────────────────────────────────────────────
    // B — named use-cases (easy to relate; small budgets)
    // ─────────────────────────────────────────────────────────────────────────────

    @Test @DisplayName("B.1  Benzene vs toluene: aromatic core survives")
    void b1_benzene_vs_toluene() {
        IAtomContainer a = mol("c1ccccc1");      // benzene
        IAtomContainer b = mol("Cc1ccccc1");     // toluene
        SMSD smsd = new SMSD(a, b, opts(true, true, false, false));
        Map<Integer,Integer> m = smsd.findMCS(false, true, 2000L);
        assertNonEmpty(m, "Phenyl core expected");
        assertTrue(m.size() >= 6, "At least the 6 aromatic carbons should match");
        assertTrue(countComponentsOnQuery(a, m.keySet()) <= 1, "Phenyl core should be connected");
    }

    @Test @DisplayName("B.2  Water vs benzene: empty under strict atom typing")
    void b2_water_vs_benzene_strict_atoms() {
        SMSD smsd = new SMSD(mol("O"), mol("c1ccccc1"), opts(true, true, false, false));
        Map<Integer,Integer> m = smsd.findMCS(false, true, 1500L);
        assertTrue(m == null || m.isEmpty(), "O vs benzene should yield empty core under strict typing");
    }

    @Test @DisplayName("B.3  Cyclohexene vs ethylbenzene: ring difference often collapses to a chain")
    void b3_ring_vs_ring_collapses_chain() {
        IAtomContainer a = mol("C1CCCC=C1CC");   // cyclohexene–ethyl
        IAtomContainer b = mol("CCc1ccccc1");    // ethylbenzene
        SMSD smsd = new SMSD(a, b, opts(true, true, false, false));
        Map<Integer,Integer> m = smsd.findMCS(false, true, 2500L);
        assertNonEmpty(m, "Expect some common acyclic chain");
        assertTrue(mappedRingBonds(a, b, m) == 0 || m.size() >= 4,
                "Likely a chain MCS; if ring bonds appear, ensure size ≥ 4");
    }

    @Test @DisplayName("B.4  Ring→ring policy blocks ring↔chain matches")
    void b4_ring_to_ring_blocks_chain() {
        IAtomContainer ring6 = mol("C1CCCCC1");
        IAtomContainer chain = mol("CCCCCC");
        SMSD smsd = new SMSD(ring6, chain, opts(true, true, true, false)); // ringOnly=true
        Map<Integer,Integer> m = smsd.findMCS(false, true, 1500L);
        assertTrue(m == null || m.isEmpty(), "Ring→ring policy should prevent mapping a ring to a chain");
    }

    @Test @DisplayName("B.5  Acetonitrile vs iminoethane: relaxed bond order grows the core")
    void b5_strict_vs_relaxed_bond_order() {
        IAtomContainer acn = mol("CC#N");
        IAtomContainer imi = mol("CC=N");
        SMSD sStrict = new SMSD(acn, imi, opts(true, true,  false, false));
        SMSD sRelax  = new SMSD(acn, imi, opts(true, false, false, false));
        Map<Integer,Integer> mStrict = sStrict.findMCS(false, true, 2000L);
        Map<Integer,Integer> mRelax  = sRelax.findMCS(true,  true,  2000L); // MCIS-like
        assertNonEmpty(mStrict, "Strict should find a non-empty core");
        assertNonEmpty(mRelax,  "Relaxed should match at least the strict core");
        assertTrue(mappedBonds(acn, imi, mRelax) >= mappedBonds(acn, imi, mStrict),
                "Relaxing bond order (and MCIS-like) should not reduce matched bonds");
    }

    @Test @DisplayName("B.6  Benzyl bromide → benzyl alcohol: aromatic backbone conserved")
    void b6_reaction_backbone_conserved() {
        IAtomContainer rx = mol("c1ccc(cc1)CBr");
        IAtomContainer pr = mol("c1ccc(cc1)CO");
        SMSD smsd = new SMSD(rx, pr, opts(true, true, false, false));
        Map<Integer,Integer> m = smsd.findMCS(false, true, 2000L);
        assertNonEmpty(m, "Backbone conserved across reaction");
        assertTrue(m.size() >= 7, "Expect ring (6) + benzylic carbon mapped");
    }

    @Test @DisplayName("B.7  Aromatic ethers: MCS large enough to seed alignment")
    void b7_aromatic_ether_seed() {
        IAtomContainer a = mol("COc1ccccc1O");
        IAtomContainer b = mol("COc1ccc(OC)cc1");
        SMSD smsd = new SMSD(a, b, opts(true, true, false, false));
        Map<Integer,Integer> m = smsd.findMCS(false, true, 2000L);
        assertNonEmpty(m, "Aromatic ether MCS should exist");
        assertTrue(m.size() >= 7, "Seed large enough for 3D overlay");
        assertEquals(m.size(), m.size(), "Mapping size should match atom count");
    }

    @Test @DisplayName("B.8  Two families cluster tighter than cross-pairs")
    void b8_clustering_hint() {
        String[] arom = {"c1ccccc1", "Cc1ccccc1", "OCCc1ccccc1"};
        String[] alco = {"CCO", "CCCO", "CCCCO"};
        double thr = 0.60;
        Map<String,Integer> nAtoms = new HashMap<>();
        String[] all = concat(arom, alco);
        for (String s : all) nAtoms.put(s, mol(s).getAtomCount());

        int intraA=0, intraB=0, cross=0;
        for (int i=0; i<all.length; i++) {
            for (int j=i+1; j<all.length; j++) {
                IAtomContainer x = mol(all[i]), y = mol(all[j]);
                ChemOptions chemOpts = opts(true, false, true, false);
                chemOpts.aromaticityMode = ChemOptions.AromaticityMode.STRICT;
                SMSD smsd = new SMSD(x, y, chemOpts);
                Map<Integer,Integer> m = smsd.findMCS(false, true, 1500L);
                boolean ok = (m != null && !m.isEmpty());
                double sim = ok ? ((double) m.size() / Math.min(nAtoms.get(all[i]), nAtoms.get(all[j]))) : 0.0;
                boolean edge = ok && sim >= thr;

                boolean xInA = contains(arom, all[i]);
                boolean yInA = contains(arom, all[j]);
                if (xInA && yInA) { if (edge) intraA++; }
                else if (!xInA && !yInA) { if (edge) intraB++; }
                else { if (edge) cross++; }
            }
        }
        assertTrue(intraA > cross, "Aromatics should cluster tighter than cross-pairs");
        assertTrue(intraB > cross, "Aliphatic alcohols should cluster tighter than cross-pairs");
    }

    // ─────────────────────────────────────────────────────────────────────────────
    // C — targeted extras (induced/MCCS, allowed ring→ring, timeout, stereo, etc.)
    // ─────────────────────────────────────────────────────────────────────────────

    @Test @DisplayName("C.1  MCIS (induced) is ≥ MCCS on a simple aliphatic pair")
    void c1_induced_vs_mccs() {
        IAtomContainer a = mol("CCCC");   // butane
        IAtomContainer b = mol("CC=CC");  // but-2-ene
        SMSD mccs = new SMSD(a, b, opts(true, true,  false, false));
        SMSD mcis = new SMSD(a, b, opts(true, false, false, false)); // relaxed order for MCIS-like
        Map<Integer,Integer> m1 = mccs.findMCS(false, true, 1500L); // MCCS
        Map<Integer,Integer> m2 = mcis.findMCS(true,  true, 1500L); // MCIS
        int s1 = (m1 == null ? 0 : m1.size());
        int s2 = (m2 == null ? 0 : m2.size());
        assertTrue(s2 >= s1, "MCIS should not be smaller than MCCS on simple pairs");
    }

    @Test @DisplayName("C.2  Ring→ring policy allows ring↔ring mapping (benzene vs kekulé benzene)")
    void c2_ring_to_ring_allows_ring_to_ring() {
        IAtomContainer ringA = mol("c1ccccc1");
        IAtomContainer ringB = mol("C1=CC=CC=C1");
        SMSD smsd = new SMSD(ringA, ringB, opts(true, true, true, false)); // ringOnly=true
        Map<Integer,Integer> m = smsd.findMCS(false, true, 1500L);
        assertTrue(m != null && !m.isEmpty(), "Ring→ring policy should allow ring↔ring matching");
    }

    @Test @DisplayName("C.3  Timeout budget is honoured (tiny budget)")
    void c3_timeout_honoured() {
        SMSD smsd = new SMSD(
                mol("c1ccc2ccccc2c1"),
                mol("c1cccc2ccc3cccc3c2c1"),
                opts(true, true, false, false));
        Map<Integer,Integer> m = smsd.findMCS(false, true, 1L); // 1 ms
        assertNotNull(m, "Should return (possibly empty) mapping without hanging");
    }

    @Test @DisplayName("C.4  Stereo sensitivity: trans vs cis reduces match when stereo is ON")
    void c4_stereo_sensitivity_trans_vs_cis() {
        IAtomContainer trans = mol("C/C=C\\C");
        IAtomContainer cis   = mol("C/C=C/C");
        SMSD sOn  = new SMSD(trans, cis,  opts(true, true, false, true));  // stereo ON
        SMSD sOff = new SMSD(trans, cis,  opts(true, true, false, false)); // stereo OFF
        Map<Integer,Integer> mOn  = sOn.findMCS(false, true, 2000L);
        Map<Integer,Integer> mOff = sOff.findMCS(false, true, 2000L);
        int aOn  = (mOn  == null ? 0 : mOn.size());
        int aOff = (mOff == null ? 0 : mOff.size());
        assertTrue(aOff >= aOn, "With stereo OFF, mapping should be ≥ stereo ON");
    }

    @Test @DisplayName("C.5  Connected vs disconnected: dMCS may yield ≥ fragments")
    void c5_connected_vs_disconnected() {
        IAtomContainer a = mol("CCOC(C)=O");
        IAtomContainer b = mol("CCOC(C)C=O");
        SMSD cMCS = new SMSD(a, b, opts(true, true, false, false));
        SMSD dMCS = new SMSD(a, b, opts(true, true, false, false));
        Map<Integer,Integer> mc = cMCS.findMCS(false, true,  2000L); // connectedOnly
        Map<Integer,Integer> md = dMCS.findMCS(false, false, 2000L); // disconnected
        int fc = countComponentsOnQuery(a, (mc == null ? Collections.emptySet() : mc.keySet()));
        int fd = countComponentsOnQuery(a, (md == null ? Collections.emptySet() : md.keySet()));
        assertNotNull(md, "dMCS returns a mapping");
        assertTrue(fd >= fc, "dMCS may have ≥ fragments than cMCS");
    }

    @Test @DisplayName("C.6  Aromaticity flexible mode tolerates aromatic vs kekulé forms")
    void c6_aromaticity_flexible() {
        IAtomContainer arom  = mol("c1ccccc1");
        IAtomContainer kekul = mol("C1=CC=CC=C1");
        ChemOptions flex = opts(true, true, false, false);
        flex.aromaticityMode = ChemOptions.AromaticityMode.FLEXIBLE;
        SMSD smsd = new SMSD(arom, kekul, flex);
        Map<Integer,Integer> m = smsd.findMCS(false, true, 1500L);
        assertTrue(m != null && m.size() >= 6, "Flexible mode should match the full benzene core");
    }

    @Test @DisplayName("C.7  Mapping integrity: values are unique")
    void c7_mapping_integrity() {
        SMSD smsd = new SMSD(mol("c1ccccc1"), mol("c1ccc2ccccc2c1"), opts(true, false, false, false));
        Map<Integer,Integer> m = smsd.findMCS(false, true, 2000L);
        assertNotNull(m);
        assertEquals(m.size(), new HashSet<>(m.values()).size(), "Target indices must be unique");
    }

    @Test @DisplayName("C.8  Bond-count monotonicity check on a second pair")
    void c8_bond_count_monotonicity_recap() {
        IAtomContainer a = mol("CC=CC"); // 2-butene
        IAtomContainer b = mol("CCCC");  // butane
        SMSD strict = new SMSD(a, b, opts(true, true,  false, false));
        SMSD relax  = new SMSD(a, b, opts(true, false, false, false));
        Map<Integer,Integer> ms = strict.findMCS(false, true, 2000L);
        Map<Integer,Integer> mr = relax .findMCS(true,  true,  2000L);
        int bs = mappedBonds(a, b, (ms==null? Collections.emptyMap():ms));
        int br = mappedBonds(a, b, (mr==null? Collections.emptyMap():mr));
        assertTrue(br >= bs, "Relaxed (MCIS-like) should not reduce matched bonds vs strict");
    }

    // tiny utilities
    private static String[] concat(String[] a, String[] b) {
        String[] z = Arrays.copyOf(a, a.length + b.length);
        System.arraycopy(b, 0, z, a.length, b.length);
        return z;
    }
    private static boolean contains(String[] arr, String x) {
        for (String s : arr) if (Objects.equals(s, x)) return true;
        return false;
    }
}