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
package org.openscience.smsd.algorithm.mcsplus;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import org.openscience.cdk.annotations.TestClass;
import org.openscience.cdk.annotations.TestMethod;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.isomorphism.matchers.IQueryAtomContainer;
import org.openscience.smsd.AtomAtomMapping;
import org.openscience.smsd.filters.PostFilter;
import org.openscience.smsd.helper.FinalMappings;
import org.openscience.smsd.interfaces.AbstractMCSAlgorithm;
import org.openscience.smsd.interfaces.IMCSBase;

/**
 * This class acts as a handler class for MCSPlus algorithm.
 * {@link org.openscience.cdk.smsd.algorithm.mcsplus.MCSPlus}
 * @cdk.module smsd
 * @cdk.githash
 * @author Syed Asad Rahman <asad@ebi.ac.uk>
 */
@TestClass("org.openscience.cdk.smsd.SMSDBondSensitiveTest")
public final class MCSPlusHandler extends AbstractMCSAlgorithm implements IMCSBase {

    private List<AtomAtomMapping> allAtomMCS = null;
    private List<Map<Integer, Integer>> allMCS = null;
    private IAtomContainer source = null;
    private IAtomContainer target = null;
    private boolean flagExchange = false;

    /**
     * Constructor for the MCS Plus algorithm class
     */
    public MCSPlusHandler() {
        allAtomMCS = Collections.synchronizedList(new ArrayList<AtomAtomMapping>());
        allMCS = Collections.synchronizedList(new ArrayList<Map<Integer, Integer>>());
    }

    /** {@inheritDoc}
     *
     * @param source
     * @param target
     */
    @Override
    @TestMethod("testSet_MolHandler_MolHandler")
    public synchronized void set(IAtomContainer source, IAtomContainer target) {
        this.source = source;
        this.target = target;
    }

    /** {@inheritDoc}
     *
     * @param source
     * @param target
     */
    @Override
    @TestMethod("testSet_IQueryAtomContainer_MolHandler")
    public synchronized void set(IQueryAtomContainer source, IAtomContainer target) {
        this.source = source;
        this.target = target;
    }

    /** {@inheritDoc}
     * Function is called by the main program and serves as a starting point for the comparison procedure.
     *
     * @param shouldMatchBonds 
     */
    @Override
    @TestMethod("testSearchMCS")
    public synchronized void searchMCS(boolean shouldMatchBonds, boolean shouldMatchRings) {
        List<List<Integer>> mappings = null;
        try {
            if (source.getAtomCount() > target.getAtomCount()) {
                mappings = Collections.synchronizedList(new MCSPlus().getOverlaps(source, target, shouldMatchBonds, false));
            } else {
                flagExchange = true;
                mappings = Collections.synchronizedList(new MCSPlus().getOverlaps(target, source, shouldMatchBonds, false));
            }
            PostFilter.filter(mappings);
            setAllMapping();
            setAllAtomMapping();
        } catch (CDKException e) {
            mappings = null;
        }
    }

    private synchronized void setAllMapping() {
        try {

            List<Map<Integer, Integer>> final_solution =
                    Collections.synchronizedList(FinalMappings.getInstance().getFinalMapping());
            int counter = 0;
            int bestSolSize = 0;
            for (Map<Integer, Integer> solution : final_solution) {
//                System.out.println("Number of MCS solution: " + solution);
                Map<Integer, Integer> validSolution = Collections.synchronizedSortedMap(new TreeMap<Integer, Integer>());
                if (!flagExchange) {
                    for (Map.Entry<Integer, Integer> map : solution.entrySet()) {
                        validSolution.put(map.getKey(), map.getValue());
                    }
                } else {
                    for (Map.Entry<Integer, Integer> map : solution.entrySet()) {
                        validSolution.put(map.getValue(), map.getKey());
                    }
//                    System.out.println("MCS solution: " + validSolution);
                }
                if (validSolution.size() > bestSolSize) {
                    bestSolSize = validSolution.size();
                    counter = 0;
                    allMCS.clear();
                }
                if (validSolution.size() == bestSolSize) {
                    allMCS.add(counter++, validSolution);
                }
            }

        } catch (Exception ex) {
            ex.printStackTrace();
        }
    }

    private synchronized void setAllAtomMapping() {
        try {

            int counter = 0;
            for (Map<Integer, Integer> solution : allMCS) {
                AtomAtomMapping atomMappings = new AtomAtomMapping(source, target);
                for (Map.Entry<Integer, Integer> map : solution.entrySet()) {

                    int IIndex = map.getKey();
                    int JIndex = map.getValue();

                    IAtom sourceAtom = null;
                    IAtom targetAtom = null;

                    sourceAtom = source.getAtom(IIndex);
                    targetAtom = target.getAtom(JIndex);
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
