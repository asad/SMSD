/*
 * SPDX-License-Identifier: Apache-2.0
 * © 2025 BioInception PVT LTD.
 */
package com.bioinception.smsd.core;

/**
 * Matching options used by SearchEngine and SMSD facade.
 * Backward-compatible defaults; public fields for simplicity.
 */
public final class ChemOptions {

    public enum BondOrderMode { STRICT, LOOSE, ANY }
    public enum AromaticityMode { STRICT, FLEXIBLE }
    /** Pluggable matcher engine. VF3 is currently an alias of VF2PP inside SearchEngine. */
    public enum MatcherEngine { VF2, VF2PP, VF3 }

    // --- Matching knobs ---
    public boolean matchAtomType       = true;
    public boolean matchFormalCharge   = true;
    public boolean useChirality        = false;
    public boolean useBondStereo       = false;
    /** If true, ring atoms/bonds in query may only match ring atoms/bonds in target. */
    public boolean ringMatchesRingOnly = true;

    public BondOrderMode   matchBondOrder  = BondOrderMode.STRICT;
    public AromaticityMode aromaticityMode = AromaticityMode.FLEXIBLE;

    /** Engine selection (default = VF2PP). */
    public MatcherEngine matcherEngine = MatcherEngine.VF2PP;

    // --- Performance / pruning toggles ---
    /** Use 2-hop NLF (neighbor-of-neighbor label frequencies) pruning. */
    public boolean useTwoHopNLF = true;
    /** Use 3-hop NLF pruning (captures ring proximity at distance 3, e.g., benzylic contexts). */
    public boolean useThreeHopNLF = true;
    /** Use bit-parallel feasibility (BitSet adjacency subset tests). */
    public boolean useBitParallelFeasibility = true;

    public ChemOptions() {}

    /**
     * FIX: Added a static factory method to create pre-configured options from a profile name.
     * This resolves the compilation errors in SMSDCasesTest.
     *
     * @param name The name of the profile (e.g., "strict", "compat-substruct").
     * @return A configured ChemOptions instance.
     */
    public static ChemOptions profile(String name) {
        ChemOptions options = new ChemOptions();
        switch (name.toLowerCase()) {
            case "strict":
                options.matchBondOrder = BondOrderMode.STRICT;
                options.aromaticityMode = AromaticityMode.STRICT;
                break;
            case "compat-fmcs":
            case "compat-substruct":
            default:
                // A general-purpose flexible configuration
                options.matchBondOrder = BondOrderMode.LOOSE;
                options.aromaticityMode = AromaticityMode.FLEXIBLE;
                options.ringMatchesRingOnly = false;
                break;
        }
        return options;
    }

    public ChemOptions withMatcherEngine(MatcherEngine e) {
        this.matcherEngine = (e == null ? MatcherEngine.VF2PP : e);
        return this;
    }
    public ChemOptions withAromaticityMode(AromaticityMode m) {
        this.aromaticityMode = (m == null ? AromaticityMode.FLEXIBLE : m);
        return this;
    }
    public ChemOptions withBondOrderMode(BondOrderMode m) {
        this.matchBondOrder = (m == null ? BondOrderMode.STRICT : m);
        return this;
    }
    public ChemOptions withTwoHopNLF(boolean on) { this.useTwoHopNLF = on; return this; }
    public ChemOptions withThreeHopNLF(boolean on){ this.useThreeHopNLF = on; return this; }
    public ChemOptions withBitParallelFeasibility(boolean on) { this.useBitParallelFeasibility = on; return this; }
}