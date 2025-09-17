/*
 * SPDX-License-Identifier: Apache-2.0
 * © 2025 BioInception PVT LTD.
 */

package com.bioinception.smsd.core;

import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.aromaticity.ElectronDonation;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.tools.CDKHydrogenAdder;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/** Conservative standardisation; optional tautomer canonicalisation via InChI (if present). */
public class Standardiser {
    public enum TautomerMode { NONE, INCHI_CANONICAL }

    public static IAtomContainer standardise(IAtomContainer mol, TautomerMode mode) throws Exception {
        IAtomContainer m = mol.clone();
        AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(m);
        CDKHydrogenAdder.getInstance(SilentChemObjectBuilder.getInstance()).addImplicitHydrogens(m);
        Cycles.markRingAtomsAndBonds(m);
        new Aromaticity(ElectronDonation.daylight(), Cycles.or(Cycles.all(), Cycles.all(6))).apply(m);

        if (mode == TautomerMode.INCHI_CANONICAL) {
            try {
                List<IAtomContainer> tauts = generateTautomersInChI(m);
                if (!tauts.isEmpty()) {
                    SmilesGenerator sg = new SmilesGenerator(SmiFlavor.Unique);
                    List<String> s = new ArrayList<>(tauts.size());
                    for (IAtomContainer t : tauts) s.add(sg.create(t));
                    Collections.sort(s);
                    SmilesParser sp = new SmilesParser(SilentChemObjectBuilder.getInstance());
                    return sp.parseSmiles(s.get(0));
                }
            } catch (Throwable ignored) {}
        }
        return m;
    }

    @SuppressWarnings("unchecked")
    private static List<IAtomContainer> generateTautomersInChI(IAtomContainer m) throws Exception {
        Class<?> clazz = Class.forName("org.openscience.cdk.tautomers.InChITautomerGenerator");
        Object gen = clazz.getConstructor().newInstance();
        return (List<IAtomContainer>) clazz.getMethod("generateTautomers", IAtomContainer.class).invoke(gen, m);
    }
}
