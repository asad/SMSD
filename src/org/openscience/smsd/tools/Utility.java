/*
 *
 *
 * Copyright (C) 2009-2014  Syed Asad Rahman <asad@ebi.ac.uk>
 *
 * Contact: cdk-devel@lists.sourceforge.net
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 * All we ask is that proper credit is given for our work, which includes
 * - but is not limited to - adding the above copyright notice to the beginning
 * of your source code files, and to any copyright notice that you may distribute
 * with programs based on this work.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received query copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA.
 * 
 */
package org.openscience.smsd.tools;

import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.aromaticity.ElectronDonation;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.graph.CycleFinder;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

/**
 *
 * @author Syed Asad Rahman <asad@ebi.ac.uk>
 */
public class Utility {

    /**
     *
     * @param molecule
     * @throws CDKException
     */
    public static void aromatizeDayLight(IAtomContainer molecule) throws CDKException {
        ElectronDonation model = ElectronDonation.daylight();
        CycleFinder cycles = Cycles.or(Cycles.all(), Cycles.relevant());
        Aromaticity aromaticity = new Aromaticity(model, cycles);
        aromaticity.apply(molecule);
    }

    /**
     *
     * @param molecule
     * @throws CDKException
     */
    public static void aromatizeCDK(IAtomContainer molecule) throws CDKException {
        ElectronDonation model = ElectronDonation.cdk();
        CycleFinder cycles = Cycles.cdkAromaticSet();
        Aromaticity aromaticity = new Aromaticity(model, cycles);
        AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(molecule);
        aromaticity.apply(molecule);
    }
}
