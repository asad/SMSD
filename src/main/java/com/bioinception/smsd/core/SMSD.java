/*
 * SPDX-License-Identifier: Apache-2.0
 * © 2025 BioInception PVT LTD.
 */
package com.bioinception.smsd.core;

import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.DefaultChemObjectBuilder;

import java.util.List;
import java.util.Map;

/**
 * SMSD (Facade)
 * -------------
 * Minimal, pure-Java facade used by the CLI and tests.
 * Wraps {@link SearchEngine} and exposes time-limited substructure and MCS.
 */
public final class SMSD {

    private final IAtomContainer query;
    private final IAtomContainer target;
    private final ChemOptions chem;

    /** Per-instance timeouts (ms). Can be overridden per call. */
    private long substructureTimeoutMs = 8_000L;
    private long mcsTimeoutMs = 10_000L;

    public SMSD(IAtomContainer query, IAtomContainer target, ChemOptions chem) {
        this.query = query;
        this.target = target;
        this.chem = (chem != null ? chem : new ChemOptions());
    }

    public SMSD(String querySmiles, String targetSmiles, ChemOptions chem) throws Exception {
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        this.query = sp.parseSmiles(querySmiles);
        this.target = sp.parseSmiles(targetSmiles);
        this.chem = (chem != null ? chem : new ChemOptions());
    }

    public void setSubstructureTimeoutMs(long ms) { this.substructureTimeoutMs = ms; }
    public void setMcsTimeoutMs(long ms) { this.mcsTimeoutMs = ms; }

    // ---- Substructure

    public boolean isSubstructure() {
        return SearchEngine.isSubstructure(query, target, chem, substructureTimeoutMs);
    }

    public boolean isSubstructure(long timeoutMs) {
        return SearchEngine.isSubstructure(query, target, chem, timeoutMs);
    }

    public List<Map<Integer,Integer>> findAllSubstructures(int maxMatches, long timeoutMs) {
        return SearchEngine.findAllSubstructures(query, target, chem, maxMatches, timeoutMs);
    }

    // ---- MCS

    /**
     * @param induced if true, compute MCIS (induced). If false, MCCS (edge preserving).
     * @param connectedOnly keep the largest connected component of the mapping
     * @return mapping qi->tj (possibly empty if timed out)
     */
    public Map<Integer,Integer> findMCS(boolean induced, boolean connectedOnly) {
        SearchEngine.McsOptions opt = new SearchEngine.McsOptions();
        opt.induced = induced;
        opt.connectedOnly = connectedOnly;
        opt.timeoutMillis = this.mcsTimeoutMs;
        return SearchEngine.findMCS(query, target, chem, opt);
    }

    public Map<Integer,Integer> findMCS(boolean induced, boolean connectedOnly, long timeoutMs) {
        SearchEngine.McsOptions opt = new SearchEngine.McsOptions();
        opt.induced = induced;
        opt.connectedOnly = connectedOnly;
        opt.timeoutMillis = timeoutMs;
        return SearchEngine.findMCS(query, target, chem, opt);
    }

    // Accessors for CLI/tests
    public IAtomContainer getQuery() { return query; }
    public IAtomContainer getTarget() { return target; }
    public ChemOptions getChem() { return chem; }
}
