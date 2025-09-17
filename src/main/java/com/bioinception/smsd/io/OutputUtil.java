/*
 * SPDX-License-Identifier: Apache-2.0
 * © 2025 BioInception PVT LTD.
 */

package com.bioinception.smsd.io;

import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.io.MDLV2000Writer;

import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.nio.charset.StandardCharsets;

public class OutputUtil {
    public enum OutType { SMI, MOL }
    public static OutType parseOutType(String s) {
        if (s==null) return OutType.SMI;
        if ("MOL".equalsIgnoreCase(s)) return OutType.MOL;
        return OutType.SMI;
    }

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
}
