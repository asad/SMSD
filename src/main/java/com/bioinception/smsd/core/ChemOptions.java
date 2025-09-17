/*
 * SPDX-License-Identifier: Apache-2.0
 * © 2025 BioInception PVT LTD.
 */
package com.bioinception.smsd.core;

public class ChemOptions {

    public enum AromaticityMode { STRICT, FLEXIBLE }
    public enum BondOrderMode  { STRICT, LOOSE, ANY }
    public enum RingSizeMode   { EXACT, SUBSET, IGNORE }

    // ── Atom/Bond flags ───────────────────────────────────────────────────────
    public boolean matchAtomType       = true;
    public boolean matchIsotope        = false;
    public boolean matchFormalCharge   = true;

    // Stereochemistry
    /** Atom-centre chirality (R/S) enforcement. */
    public boolean useChirality        = false;
    /** Bond stereo (E/Z, cis/trans) enforcement. */
    public boolean useBondStereo       = false;

    // Flexibility toggles
    public AromaticityMode aromaticityMode = AromaticityMode.FLEXIBLE;
    public BondOrderMode  matchBondOrder   = BondOrderMode.LOOSE;
    public boolean        ringMatchesRingOnly = true;

    /** Allow mapping an atom of degree d_q to an atom of degree d_t when d_q ≤ d_t + degreeSlack. */
    public int            degreeSlack    = 1;

    // Optional ring-size handling (not currently used in SearchEngine, but kept for completeness)
    public RingSizeMode   ringSizeMode     = RingSizeMode.SUBSET;
    public int            ringSizeTolerance= 1;

    public ChemOptions() {}

    /**
     * Named presets to keep behaviour consistent across tests and CLI.
     *  - "compat-fmcs": FMCS-like (ignore stereo; flexible aromatics/bond order; ring→ring; degree slack 1)
     *  - "compat-substruct": substructure-friendly flexible defaults
     *  - "strict": enforce atom type, charge, aromaticity, exact bond order, and stereo; no degree slack
     */
    public static ChemOptions profile(String name) {
        ChemOptions c = new ChemOptions();
        String key = name == null ? "" : name.trim().toLowerCase();
        switch (key) {
            case "strict":
                c.matchAtomType = true;
                c.matchIsotope = false;
                c.matchFormalCharge = true;
                c.useChirality = true;          // atom stereo
                c.useBondStereo = true;         // bond stereo
                c.aromaticityMode = AromaticityMode.STRICT;
                c.matchBondOrder = BondOrderMode.STRICT;
                c.ringMatchesRingOnly = true;
                c.ringSizeMode = RingSizeMode.EXACT;
                c.ringSizeTolerance = 0;
                c.degreeSlack = 0;
                return c;

            case "compat-fmcs":
                c.useChirality = false;
                c.useBondStereo = false;
                c.aromaticityMode = AromaticityMode.FLEXIBLE;
                c.matchBondOrder = BondOrderMode.LOOSE;
                c.ringMatchesRingOnly = true;
                c.ringSizeMode = RingSizeMode.SUBSET;
                c.ringSizeTolerance = 1;
                c.degreeSlack = 1;
                return c;

            case "compat-substruct":
            default:
                c.useChirality = false;
                c.useBondStereo = false;
                c.aromaticityMode = AromaticityMode.FLEXIBLE;
                c.matchBondOrder = BondOrderMode.LOOSE;
                c.ringMatchesRingOnly = true;
                c.ringSizeMode = RingSizeMode.SUBSET;
                c.ringSizeTolerance = 1;
                c.degreeSlack = 1;
                return c;
        }
    }
}
