/* Copyright (C) 2006-2011  Syed Asad Rahman <asad@ebi.ac.uk>
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
 * MERCHANTABILITY or FITNESS FOR sourceAtom PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA.
 */
package org.openscience.smsd.algorithm.single;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.logging.Level;
import org.openscience.cdk.annotations.TestClass;
import org.openscience.cdk.annotations.TestMethod;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.isomorphism.matchers.IQueryAtomContainer;
import org.openscience.cdk.tools.ILoggingTool;
import org.openscience.cdk.tools.LoggingToolFactory;
import org.openscience.smsd.AtomAtomMapping;
import org.openscience.smsd.helper.MoleculeInitializer;
import org.openscience.smsd.interfaces.IResults;

/**
 * This is a handler class for single atom mapping
 * ({@link org.openscience.cdk.smsd.algorithm.single.SingleMapping}).
 * @cdk.module smsd
 * @cdk.githash
 * @author Syed Asad Rahman <asad@ebi.ac.uk>
 */
@TestClass("org.openscience.cdk.smsd.algorithm.single.SingleMappingHandlerTest")
public class SingleMappingHandler extends MoleculeInitializer implements IResults {

    private final ILoggingTool Logger =
            LoggingToolFactory.createLoggingTool(SingleMappingHandler.class);
    private List<AtomAtomMapping> allAtomMCS = null;
    private List<Map<Integer, Integer>> allMCS = null;
    private final IAtomContainer source;
    private final IAtomContainer target;
    private final boolean shouldMatchRings;

    /**
     * 
     * @param source
     * @param target
     * @param bondTypeMatch 
     * @param shouldMatchRings  
     */
    @TestMethod("setMCSAlgorithm")
    public SingleMappingHandler(IAtomContainer source, IAtomContainer target, boolean bondTypeMatch, boolean shouldMatchRings) {
        allAtomMCS = new ArrayList<AtomAtomMapping>();
        allMCS = new ArrayList<Map<Integer, Integer>>();
        this.source = source;
        this.target = target;
        this.shouldMatchRings = shouldMatchRings;
        if (this.shouldMatchRings) {
            try {
                initializeMolecule(source);
                initializeMolecule(target);
            } catch (CDKException ex) {
            }
        }
        searchMCS();
    }

    /**
     * 
     * @param source
     * @param target  
     */
    @TestMethod("setMCSAlgorithm")
    public SingleMappingHandler(IQueryAtomContainer source, IAtomContainer target) {
        allAtomMCS = new ArrayList<AtomAtomMapping>();
        allMCS = new ArrayList<Map<Integer, Integer>>();
        this.source = source;
        this.target = target;
        this.shouldMatchRings = true;
        if (this.shouldMatchRings) {
            try {
                initializeMolecule(source);
                initializeMolecule(target);
            } catch (CDKException ex) {
            }
        }
        searchMCS();
    }

    //Function is called by the main program and serves as a starting point for the comparision procedure.
    /** {@inheritDoc}
     *
     */
    private synchronized void searchMCS() {
        SingleMapping singleMapping = new SingleMapping();
        List<Map<IAtom, IAtom>> mappings = null;
        try {
            if (!(source instanceof IQueryAtomContainer)) {
                mappings = singleMapping.getOverLaps(source, target);
            } else {
                mappings = singleMapping.getOverLaps((IQueryAtomContainer) source, target);
            }
        } catch (CDKException ex) {
            Logger.error(Level.SEVERE, null, ex);
        }
        setAllAtomMapping(mappings);
        setAllMapping(mappings);
        //setStereoScore();
    }

    /** {@inheritDoc}
     *
     * Set the mappings
     */
    private synchronized void setAllMapping(List<Map<IAtom, IAtom>> mappings) {
        try {
            int counter = 0;
            for (Map<IAtom, IAtom> solution : mappings) {
                Map<Integer, Integer> atomMappings = new TreeMap<Integer, Integer>();
                for (Map.Entry<IAtom, IAtom> map : solution.entrySet()) {
                    IAtom sourceAtom = map.getKey();
                    IAtom targetAtom = map.getValue();
                    atomMappings.put(source.getAtomNumber(sourceAtom), target.getAtomNumber(targetAtom));
                }
                allMCS.add(counter++, atomMappings);
            }
        } catch (Exception I) {
            I.getCause();
        }
    }

    private synchronized void setAllAtomMapping(List<Map<IAtom, IAtom>> mappings) {

        try {
            int counter = 0;
            for (Map<IAtom, IAtom> solution : mappings) {
                AtomAtomMapping atomMappings = new AtomAtomMapping(source, target);
                for (Map.Entry<IAtom, IAtom> map : solution.entrySet()) {

                    IAtom sourceAtom = map.getKey();
                    IAtom targetAtom = map.getValue();
                    atomMappings.put(sourceAtom, targetAtom);
                }
                allAtomMCS.add(counter++, atomMappings);
            }
        } catch (Exception I) {
            I.getCause();
        }
    }

    /** {@inheritDoc}
     */
    @Override
    @TestMethod("testGetAllMapping")
    public synchronized List<Map<Integer, Integer>> getAllMapping() {
        return allMCS;
    }

    /** {@inheritDoc}
     */
    @Override
    @TestMethod("testGetFirstMapping")
    public synchronized Map<Integer, Integer> getFirstMapping() {
        if (allMCS.iterator().hasNext()) {
            return allMCS.iterator().next();
        }
        return new TreeMap<Integer, Integer>();
    }

    /** {@inheritDoc}
     */
    @Override
    @TestMethod("testGetAllAtomMapping")
    public synchronized List<AtomAtomMapping> getAllAtomMapping() {
        return allAtomMCS;
    }

    /** {@inheritDoc}
     */
    @Override
    @TestMethod("testGetFirstAtomMapping")
    public synchronized AtomAtomMapping getFirstAtomMapping() {
        if (allAtomMCS.iterator().hasNext()) {
            return allAtomMCS.iterator().next();
        }
        return new AtomAtomMapping(source, target);
    }
}
