/* Copyright (C) 2009-2011  Syed Asad Rahman <asad@ebi.ac.uk>
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
 * You should have received commonAtomList copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA.
 */
package org.openscience.smsd.algorithm.vflib;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.logging.Level;
import org.openscience.cdk.annotations.TestClass;
import org.openscience.cdk.annotations.TestMethod;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.isomorphism.matchers.IQueryAtomContainer;
import org.openscience.cdk.tools.ILoggingTool;
import org.openscience.cdk.tools.LoggingToolFactory;
import org.openscience.smsd.AtomAtomMapping;
import org.openscience.smsd.interfaces.IResults;

/**
 * This class should be used to find MCS between source
 * graph and target graph.
 *
 * First the algorithm runs VF lib {@link org.openscience.cdk.smsd.algorithm.vflib.map.VFMCSMapper}
 * and reports MCS between
 * run source and target graphs. Then these solutions are extended
 * using McGregor {@link org.openscience.cdk.smsd.algorithm.mcgregor.McGregor}
 * algorithm where ever required.
 *
 * @cdk.module smsd
 * @cdk.githash
 * @author Syed Asad Rahman <asad@ebi.ac.uk>
 */
@TestClass("org.openscience.cdk.smsd.algorithm.vflib.VFlibMCSHandlerTest")
public class VF2MCS extends BaseMCS implements IResults {

    private List<AtomAtomMapping> allAtomMCS = null;
    private List<Map<Integer, Integer>> allMCS = null;
    private final static ILoggingTool logger =
            LoggingToolFactory.createLoggingTool(VF2MCS.class);

    /**
     * Constructor for an extended VF Algorithm for the MCS search
     * @param source
     * @param target
     * @param shouldMatchBonds bond match
     * @param shouldMatchRings ring match 
     */
    public VF2MCS(IAtomContainer source, IAtomContainer target, boolean shouldMatchBonds, boolean shouldMatchRings) {
        super(source, target, shouldMatchBonds, shouldMatchRings);
        this.allAtomMCS = new ArrayList<AtomAtomMapping>();
        this.allMCS = new ArrayList<Map<Integer, Integer>>();
        searchMCS();
    }

    /**
     * Constructor for an extended VF Algorithm for the MCS search
     * @param source
     * @param target  
     */
    public VF2MCS(IQueryAtomContainer source, IQueryAtomContainer target) {
        super(source, target, true, true);
        this.allAtomMCS = new ArrayList<AtomAtomMapping>();
        this.allMCS = new ArrayList<Map<Integer, Integer>>();
        searchMCS();
    }

    /**
     *{@inheritDoc}
     *
     */
    private synchronized void searchMCS() {
        allMCS.clear();
        allAtomMCS.clear();

        try {
            addVFMatchesMappings();
//            System.out.println("\n\naddVFMatchesMappings() " + getLocalMCSSolution().size() + ", Solutions: " + getLocalMCSSolution());
            if (isExtensionRequired()) {
                addKochCliques();
//                System.out.println(" \nAfter addKochCliques() " + getLocalMCSSolution().size() + ", Solutions: " + getLocalMCSSolution());
                addUIT();
//                System.out.println(" \nAfter addUIT(); " + getLocalMCSSolution().size() + ", Solutions: " + getLocalMCSSolution());

                List<Map<Integer, Integer>> mcsSeeds = new ArrayList<Map<Integer, Integer>>();
                int solSize = 0;
                int counter = 0;
                if (!getLocalMCSSolution().isEmpty()) {
                    for (int i = 0; i < getLocalMCSSolution().size(); i++) {
                        Map<Integer, Integer> map = getLocalMCSSolution().get(i);
                        if (map.size() > solSize) {
                            solSize = map.size();
                            mcsSeeds.clear();
                            counter = 0;
                        }
                        if (!map.isEmpty()
                                && map.size() == solSize
                                && !hasClique(map, mcsSeeds)) {
                            mcsSeeds.add(counter, map);
                            counter++;
                        }
                    }
                }
                getLocalMCSSolution().clear();
                getLocalAtomMCSSolution().clear();
                /*
                 * Sort biggest clique to smallest
                 */
                Collections.sort(mcsSeeds, new Map2Comparator());
                extendCliquesWithMcGregor(mcsSeeds);
            }
        } catch (CDKException ex) {
            logger.error(Level.SEVERE, null, ex);
        } catch (IOException ex) {
            logger.error(Level.SEVERE, null, ex);
        }
        int solSize = 0;
        int counter = 0;
        if (!getLocalMCSSolution().isEmpty()) {
            for (int i = 0; i < getLocalMCSSolution().size(); i++) {
                Map<Integer, Integer> map = getLocalMCSSolution().get(i);
                AtomAtomMapping atomMCSMap = getLocalAtomMCSSolution().get(i);
                if (map.size() > solSize) {
                    solSize = map.size();
                    allMCS.clear();
                    allAtomMCS.clear();
                    counter = 0;
                }
                if (!map.isEmpty()
                        && map.size() == solSize) {
                    allMCS.add(counter, map);
                    allAtomMCS.add(counter, atomMCSMap);
                    counter++;
                }
            }
        }
    }

    /** {@inheritDoc}
     * @return 
     */
    @Override
    @TestMethod("testGetAllMapping")
    public synchronized List<Map<Integer, Integer>> getAllMapping() {
        return allMCS;
    }

    /** {@inheritDoc}
     * @return 
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
     * @return 
     */
    @Override
    @TestMethod("testGetAllAtomMapping")
    public synchronized List<AtomAtomMapping> getAllAtomMapping() {
        return allAtomMCS;
    }

    /** {@inheritDoc}
     * @return 
     */
    @Override
    @TestMethod("testGetFirstAtomMapping")
    public synchronized AtomAtomMapping getFirstAtomMapping() {
        if (allAtomMCS.iterator().hasNext()) {
            return allAtomMCS.iterator().next();
        }
        return new AtomAtomMapping(getReactantMol(), getProductMol());
    }
}
