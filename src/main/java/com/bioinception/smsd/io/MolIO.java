/*
 * SPDX-License-Identifier: Apache-2.0
 * © 2025 BioInception PVT LTD.
 */

package com.bioinception.smsd.io;

import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.io.MDLV2000Reader;
import org.openscience.cdk.io.Mol2Reader;
import org.openscience.cdk.io.CMLReader;
import org.openscience.cdk.io.PDBReader;
import org.openscience.cdk.io.iterator.IteratingSDFReader;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;

import java.io.*;
import java.util.ArrayList;
import java.util.List;

public class MolIO {

    public static class Query {
        public final boolean isSmarts;
        public final String text;
        public final IAtomContainer container;
        public Query(boolean isSmarts, String text, IAtomContainer container) {
            this.isSmarts = isSmarts; this.text = text; this.container = container;
        }
    }

    public static Query loadQuery(String type, String value) throws Exception {
        switch (type.toUpperCase()) {
            case "SIG": // SMARTS
                return new Query(true, value, null);
            case "SMI":
                return new Query(false, value, new SmilesParser(SilentChemObjectBuilder.getInstance()).parseSmiles(value));
            case "MOL": {
                try (FileInputStream fis = new FileInputStream(value);
                     MDLV2000Reader r = new MDLV2000Reader(fis)) {
                    return new Query(false, value, r.read(DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainer.class)));
                }
            }
            case "ML2": {
                try (FileInputStream fis = new FileInputStream(value);
                     Mol2Reader r = new Mol2Reader(fis)) {
                    return new Query(false, value, r.read(DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainer.class)));
                }
            }
            case "PDB": {
                try (FileInputStream fis = new FileInputStream(value);
                     PDBReader r = new PDBReader(fis)) {
                    return new Query(false, value, r.read(DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainer.class)));
                }
            }
            case "CML": {
                try (FileInputStream fis = new FileInputStream(value);
                     CMLReader r = new CMLReader(fis)) {
                    return new Query(false, value, r.read(DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainer.class)));
                }
            }
            default:
                throw new IllegalArgumentException("Unsupported query type: " + type);
        }
    }

    public static IAtomContainer loadTarget(String type, String value) throws Exception {
        switch (type.toUpperCase()) {
            case "SMI":
                return new SmilesParser(SilentChemObjectBuilder.getInstance()).parseSmiles(value);
            case "MOL": {
                try (FileInputStream fis = new FileInputStream(value);
                     MDLV2000Reader r = new MDLV2000Reader(fis)) {
                    return r.read(DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainer.class));
                }
            }
            case "ML2": {
                try (FileInputStream fis = new FileInputStream(value);
                     Mol2Reader r = new Mol2Reader(fis)) {
                    return r.read(DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainer.class));
                }
            }
            case "PDB": {
                try (FileInputStream fis = new FileInputStream(value);
                     PDBReader r = new PDBReader(fis)) {
                    return r.read(DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainer.class));
                }
            }
            case "CML": {
                try (FileInputStream fis = new FileInputStream(value);
                     CMLReader r = new CMLReader(fis)) {
                    return r.read(DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainer.class));
                }
            }
            case "SDF":
                throw new IllegalArgumentException("Use loadTargetsSDF for multi-molecule SDF.");
            default:
                throw new IllegalArgumentException("Unsupported target type: " + type);
        }
    }

    public static List<IAtomContainer> loadTargetsSDF(String path) throws Exception {
        List<IAtomContainer> out = new ArrayList<IAtomContainer>();
        try (FileInputStream fis = new FileInputStream(path);
             IteratingSDFReader it = new IteratingSDFReader(fis, SilentChemObjectBuilder.getInstance())) {
            while (it.hasNext()) out.add(it.next());
        }
        return out;
    }
}