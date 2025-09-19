/*
 * SPDX-License-Identifier: Apache-2.0
 * © 2025 BioInception PVT LTD.
 */

package com.bioinception.smsd.io;

import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.io.MDLV2000Writer;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;

import java.io.ByteArrayOutputStream;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.nio.charset.StandardCharsets;
import java.util.*;

/**
 * Output helpers: single molecule IO and mapping export.
 * NOTE: no CDK substructure/MCS utilities are used here.
 */
public final class OutputUtil {

    private OutputUtil() {}

    /** Output types accepted across helpers. */
    public enum OutType { JSON, SMI, SMARTS, MOL }

    /** Write a single molecule as SMILES (default) or MDL V2000 MOL. */
    public static void writeMolecule(IAtomContainer mol, OutType type, OutputStream os) throws Exception {
        switch (type) {
            case MOL:
                try (MDLV2000Writer w = new MDLV2000Writer(os)) { w.write(mol); }
                break;
            case SMI:
            default:
                SmilesGenerator sg = new SmilesGenerator(SmiFlavor.Unique);
                String s = sg.create(mol);
                try (OutputStreamWriter w = new OutputStreamWriter(os, StandardCharsets.UTF_8)) {
                    w.write(s); w.write(System.lineSeparator()); w.flush();
                }
        }
    }

    /**
     * Write substructure mappings in the requested format.
     * - JSON (default): an object with query/target sizes and a 'mappings' array of {index, pairs, query_atoms, target_atoms}.
     * - SMI: one line per mapping, SMILES of query and target with consistent atom-map numbers.
     * - SMARTS: same as SMI (SMILES is a subset of SMARTS) with atom-map numbers.
     * - MOL: one MOL block per mapping, extracted target subgraph induced by mapped atoms.
     *
     * No CDK substructure search or MCS is performed here; we format what was provided.
     */
    public static void writeMappings(IAtomContainer query,
                                     IAtomContainer target,
                                     java.util.List<java.util.Map<Integer,Integer>> mappings,
                                     OutType type,
                                     OutputStream os,
                                     boolean prettyJson) throws Exception {
        Objects.requireNonNull(target, "target");
        Objects.requireNonNull(mappings, "mappings");
        switch (type) {
            case JSON: {
                Map<String,Object> root = new LinkedHashMap<>();
                root.put("query_size", query != null ? query.getAtomCount() : null);
                root.put("target_size", target.getAtomCount());
                java.util.List<Map<String,Object>> out = new ArrayList<>();
                int idx = 0;
                for (java.util.Map<Integer,Integer> m : mappings) {
                    Map<String,Object> row = new LinkedHashMap<>();
                    row.put("index", idx++);
                    row.put("pairs", toPairs(m));
                    row.put("query_atoms", new ArrayList<>(m.keySet()));
                    row.put("target_atoms", new ArrayList<>(m.values()));
                    out.add(row);
                }
                root.put("mappings", out);
                String payload = (prettyJson
                        ? new com.fasterxml.jackson.databind.ObjectMapper()
                            .writerWithDefaultPrettyPrinter().writeValueAsString(root)
                        : new com.fasterxml.jackson.databind.ObjectMapper().writeValueAsString(root));
                try (OutputStreamWriter w = new OutputStreamWriter(os, StandardCharsets.UTF_8)) {
                    w.write(payload); w.flush();
                }
                break;
            }
            case SMI:
            case SMARTS: {
                // Assign atom-map numbers on both query and target so the two lines share numbering.
                final int flavor = SmiFlavor.Unique | SmiFlavor.UseAromaticSymbols | SmiFlavor.AtomAtomMap;
                SmilesGenerator sg = new SmilesGenerator(flavor);
                // Work on originals but restore properties after each mapping
                Map<IAtom, Object> qBackup = snapshotAtomMaps(query);
                Map<IAtom, Object> tBackup = snapshotAtomMaps(target);
                try (OutputStreamWriter w = new OutputStreamWriter(os, StandardCharsets.UTF_8)) {
                    int idx = 0;
                    for (java.util.Map<Integer,Integer> m : mappings) {
                        // assign 1..k labels by pair order
                        int n = 1;
                        for (Map.Entry<Integer,Integer> e : m.entrySet()) {
                            if (query != null) {
                                IAtom qa = query.getAtom(e.getKey());
                                qa.setProperty(CDKConstants.ATOM_ATOM_MAPPING, Integer.valueOf(n));
                            }
                            IAtom ta = target.getAtom(e.getValue());
                            ta.setProperty(CDKConstants.ATOM_ATOM_MAPPING, Integer.valueOf(n));
                            n++;
                        }
                        String qs = (query != null ? sg.create(query) : "");
                        String ts = sg.create(target);
                        if (type == OutType.SMARTS && qs.isEmpty()) {
                            // For SMARTS without a query container, just output target line
                            w.write(ts);
                        } else if (type == OutType.SMARTS) {
                            w.write(qs); w.write('	'); w.write(ts);
                        } else { // SMI
                            w.write(qs); w.write('	'); w.write(ts);
                        }
                        w.write(System.lineSeparator());
                        // clear maps for next round
                        restoreAtomMaps(query, qBackup);
                        restoreAtomMaps(target, tBackup);
                        idx++;
                    }
                    w.flush();
                }
                break;
            }
            case MOL: {
                // Extract induced subgraph of target by mapped atoms; write one MOL block per mapping
                int idx = 0;
                for (java.util.Map<Integer,Integer> m : mappings) {
                    IAtomContainer sub = inducedSubgraph(target, new java.util.HashSet<>(m.values()));
                    try (MDLV2000Writer w = new MDLV2000Writer(os)) {
                        w.write(sub);
                    }
                    // Separate MOL blocks
                    try (OutputStreamWriter w = new OutputStreamWriter(os, StandardCharsets.UTF_8)) {
                        w.write(System.lineSeparator());
                    }
                    idx++;
                }
                break;
            }
            default:
                throw new IllegalArgumentException("Unsupported mapping out type: " + type);
        }
    }

    // ---------- helpers ----------

    private static List<List<Integer>> toPairs(Map<Integer,Integer> m) {
        List<List<Integer>> out = new ArrayList<>();
        for (Map.Entry<Integer,Integer> e : m.entrySet()) {
            out.add(java.util.Arrays.asList(e.getKey(), e.getValue()));
        }
        return out;
    }

    /** Take a snapshot of CDKConstants.ATOM_ATOM_MAPPING for each atom (may be null). */
    private static Map<IAtom, Object> snapshotAtomMaps(IAtomContainer mol) {
        Map<IAtom, Object> snapshot = new IdentityHashMap<>();
        if (mol == null) return snapshot;
        for (IAtom a : mol.atoms()) snapshot.put(a, a.getProperty(CDKConstants.ATOM_ATOM_MAPPING));
        return snapshot;
    }

    /** Restore previously saved atom-map numbers. */
    private static void restoreAtomMaps(IAtomContainer mol, Map<IAtom, Object> snapshot) {
        if (mol == null || snapshot == null) return;
        for (Map.Entry<IAtom,Object> e : snapshot.entrySet()) e.getKey().setProperty(CDKConstants.ATOM_ATOM_MAPPING, e.getValue());
    }

    /** Build the induced subgraph on the given atom indices (target indices). */
    private static IAtomContainer inducedSubgraph(IAtomContainer target, java.util.Set<Integer> atomIdxs) {
        IAtomContainer sub = target.getBuilder().newInstance(IAtomContainer.class);
        // map old atom index -> new atom reference
        java.util.Map<Integer,IAtom> idx2new = new java.util.HashMap<>();
        for (Integer i : atomIdxs) {
            IAtom a = target.getAtom(i);
            IAtom b = sub.newAtom(a.getAtomicNumber() != null ? a.getAtomicNumber() : 0);
            // copy coordinates & charge if available
            b.setPoint2d(a.getPoint2d());
            b.setPoint3d(a.getPoint3d());
            b.setCharge(a.getCharge());
            b.setFormalCharge(a.getFormalCharge());
            idx2new.put(i, b);
        }
        // add bonds where both ends are in the subset
        for (IBond bo : target.bonds()) {
            int ai = target.getAtomNumber(bo.getAtom(0));
            int bi = target.getAtomNumber(bo.getAtom(1));
            if (atomIdxs.contains(ai) && atomIdxs.contains(bi)) {
                sub.newBond(idx2new.get(ai), idx2new.get(bi), bo.getOrder());
            }
        }
        return sub;
    }
}
