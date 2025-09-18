/*
 * SPDX-License-Identifier: Apache-2.0
 * © 2025 BioInception PVT LTD.
 */
/* * SPDX-License-Identifier: Apache-2.0
 * © 2025 BioInception PVT LTD.
 */
package com.bioinception.smsd.core;

import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.DefaultChemObjectBuilder;

/**
 * SMSD (Facade)
 * -------------
 * Standardises molecules on construction and delegates to SearchEngine.
 */
public final class SMSD {

    private final IAtomContainer query;
    private final IAtomContainer target;
    private final ChemOptions chem;
    // FIX: Add a field to hold the query if it's a SMARTS string.
    private final String querySmarts;

    /** Per-instance timeouts (ms). Can be overridden per call. */
    private long substructureTimeoutMs = 8_000L;
    private long mcsTimeoutMs = 10_000L;

    public SMSD(IAtomContainer query, IAtomContainer target, ChemOptions chem) {
        this.chem = (chem != null ? chem : new ChemOptions());
        this.querySmarts = null; // This constructor handles molecules
        IAtomContainer q = query, t = target;
        try {
            q = Standardiser.standardise(query, Standardiser.TautomerMode.NONE);
            t = Standardiser.standardise(target, Standardiser.TautomerMode.NONE);
        } catch (Exception ignored) {}
        this.query = q;
        this.target = t;
    }

    /**
     * FIX: Add a new constructor to explicitly handle SMARTS queries.
     * This avoids mis-parsing SMARTS as SMILES.
     *
     * @param querySmarts The SMARTS string query.
     * @param target The target molecule.
     * @param chem Chemistry options.
     */
    public SMSD(String querySmarts, IAtomContainer target, ChemOptions chem) {
        this.chem = (chem != null ? chem : new ChemOptions());
        this.querySmarts = querySmarts;
        this.query = null; // Query is SMARTS, not a molecule container
        IAtomContainer t = target;
        try {
            t = Standardiser.standardise(target, Standardiser.TautomerMode.NONE);
        } catch (Exception ignored) {}
        this.target = t;
    }

    public SMSD(String querySmiles, String targetSmiles, ChemOptions chem) throws Exception {
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer q = sp.parseSmiles(querySmiles);
        IAtomContainer t = sp.parseSmiles(targetSmiles);
        this.chem = (chem != null ? chem : new ChemOptions());
        this.querySmarts = null;
        try {
            q = Standardiser.standardise(q, Standardiser.TautomerMode.NONE);
            t = Standardiser.standardise(t, Standardiser.TautomerMode.NONE);
        } catch (Exception ignored) {}
        this.query = q;
        this.target = t;
    }

    public void setSubstructureTimeoutMs(long ms) { this.substructureTimeoutMs = ms; }
    public void setMcsTimeoutMs(long ms) { this.mcsTimeoutMs = ms; }

    // ---- Substructure

    public boolean isSubstructure() {
        return isSubstructure(substructureTimeoutMs);
    }

    public boolean isSubstructure(long timeoutMs) {
        // FIX: Delegate to SmartsUtil if the query is SMARTS, otherwise use SearchEngine.
        if (this.querySmarts != null) {
            try {
                // Note: The underlying SMARTS engine may not respect the timeout.
                return !SmartsUtil.matchAll(this.querySmarts, this.target).isEmpty();
            } catch (Exception e) {
                // Propagate exception or return false
                throw new RuntimeException("SMARTS matching failed", e);
            }
        }
        return SearchEngine.isSubstructure(query, target, chem, timeoutMs);
    }

    /** Enumerate up to maxSolutions substructure mappings (injective query->target). */
    public java.util.List<java.util.Map<Integer,Integer>> findAllSubstructures(int maxSolutions, long timeoutMs) {
        if (this.querySmarts != null) {
            throw new UnsupportedOperationException("findAllSubstructures is not supported for SMARTS queries in this facade.");
        }
        return SearchEngine.findAllSubstructures(query, target, chem, maxSolutions, timeoutMs);
    }

    /** Same as isSubstructure(...) but returns telemetry and no mappings. */
    public SearchEngine.SubstructureResult isSubstructureWithStats(long timeoutMs) {
         if (this.querySmarts != null) {
            throw new UnsupportedOperationException("isSubstructureWithStats is not supported for SMARTS queries in this facade.");
        }
        return SearchEngine.isSubstructureWithStats(query, target, chem, timeoutMs);
    }

    /** Same as findAllSubstructures(...) but returns telemetry and mappings. */
    public SearchEngine.SubstructureResult findAllSubstructuresWithStats(int maxSolutions, long timeoutMs) {
         if (this.querySmarts != null) {
            throw new UnsupportedOperationException("findAllSubstructuresWithStats is not supported for SMARTS queries in this facade.");
        }
        return SearchEngine.findAllSubstructuresWithStats(query, target, chem, maxSolutions, timeoutMs);
    }

    // ---- MCS (Maximum Common Substructure)

    public java.util.Map<Integer,Integer> findMCS() {
        SearchEngine.McsOptions opt = new SearchEngine.McsOptions();
        opt.timeoutMillis = mcsTimeoutMs;
        return SearchEngine.findMCS(query, target, chem, opt);
    }

    public java.util.Map<Integer,Integer> findMCS(long timeoutMs) {
        SearchEngine.McsOptions opt = new SearchEngine.McsOptions();
        opt.timeoutMillis = timeoutMs;
        return SearchEngine.findMCS(query, target, chem, opt);
    }

    public java.util.Map<Integer,Integer> findMCS(boolean induced, boolean connectedOnly, long timeoutMs) {
        SearchEngine.McsOptions opt = new SearchEngine.McsOptions();
        opt.induced = induced;
        opt.connectedOnly = connectedOnly;
        opt.timeoutMillis = timeoutMs;
        return SearchEngine.findMCS(query, target, chem, opt);
    }

    // Accessors for CLI/tests
    public org.openscience.cdk.interfaces.IAtomContainer getQuery() { return query; }
    public org.openscience.cdk.interfaces.IAtomContainer getTarget() { return target; }
    public ChemOptions getChem() { return chem; }
}