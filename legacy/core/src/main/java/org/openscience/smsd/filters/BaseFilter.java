/*
 * SPDX-License-Identifier: Apache-2.0
 * © 2025 BioInception PVT LTD.
 */
package org.openscience.smsd.filters;

import java.util.logging.Level;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.isomorphism.matchers.IQueryAtomContainer;
import org.openscience.cdk.tools.ILoggingTool;
import org.openscience.cdk.tools.LoggingToolFactory;
import org.openscience.smsd.tools.ExtAtomContainerManipulator;

/**
 * @author Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>
 * @author maclean
 *
 */
public class BaseFilter {

    private final IAtomContainer mol1;
    private final IAtomContainer mol2;
    private final static ILoggingTool logger
            = LoggingToolFactory.createLoggingTool(BaseFilter.class);

    /**
     *
     * @param sourceMol
     * @param targetMol
     */
    public BaseFilter(IAtomContainer sourceMol, IAtomContainer targetMol) {
        try {
            ExtAtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(sourceMol);
        } catch (CDKException ex) {
            logger.error(Level.SEVERE, null, ex);
        }
        try {
            ExtAtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(targetMol);
        } catch (CDKException ex) {
            logger.error(Level.SEVERE, null, ex);
        }
        this.mol1 = sourceMol;
        this.mol2 = targetMol;

    }

    /**
     *
     * @param sourceMol
     * @param targetMol
     */
    public BaseFilter(IQueryAtomContainer sourceMol, IAtomContainer targetMol) {
        this.mol1 = sourceMol;
        this.mol2 = targetMol;

        try {
            ExtAtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(mol2);
        } catch (CDKException ex) {
            logger.error(Level.SEVERE, null, ex);
        }
    }

    /**
     * @return the mol1
     */
    public synchronized IAtomContainer getQuery() {
        return mol1;
    }

    /**
     * @return the mol2
     */
    public synchronized IAtomContainer getTarget() {
        return mol2;
    }
}
