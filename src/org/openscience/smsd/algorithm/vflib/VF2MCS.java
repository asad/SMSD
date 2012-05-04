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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received commonAtomList copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA.
 */
package org.openscience.smsd.algorithm.vflib;

import java.io.IOException;
import java.util.*;
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
 * This class should be used to find MCS between source graph and target graph.
 *
 * First the algorithm runs VF lib {@link org.openscience.cdk.smsd.algorithm.vflib.map.VFMCSMapper} and reports MCS
 * between run source and target graphs. Then these solutions are extended using McGregor
 * {@link org.openscience.cdk.smsd.algorithm.mcgregor.McGregor} algorithm where ever required.
 *
 * @cdk.module smsd@cdk.githash
 *
 * @author Syed Asad Rahman <asad@ebi.ac.uk>
 */
@TestClass("org.openscience.cdk.smsd.algorithm.vflib.VFlibMCSHandlerTest")
public final class VF2MCS extends BaseMCS implements IResults {

    private final List<AtomAtomMapping> allAtomMCS;
    private final static ILoggingTool logger =
            LoggingToolFactory.createLoggingTool(VF2MCS.class);

    /**
     * Constructor for an extended VF Algorithm for the MCS search
     *
     * @param source
     * @param target
     * @param shouldMatchBonds bond match
     * @param shouldMatchRings ring match
     */
    public VF2MCS(IAtomContainer source, IAtomContainer target, boolean shouldMatchBonds, boolean shouldMatchRings) {
        super(source, target, shouldMatchBonds, shouldMatchRings);
        this.allAtomMCS = new ArrayList<AtomAtomMapping>();

        try {
            addVFMatchesMappings();
            /*
             * An extension is triggered if its mcs solution is smaller than reactant and product. An enrichment is
             * triggered if its mcs solution is equal to reactant or product size.
             *
             *
             */
            if (isEnrichmentRequired() || isExtensionRequired()) {
                List<Map<Integer, Integer>> mcsSeeds = new ArrayList<Map<Integer, Integer>>();
                List<AtomAtomMapping> mcsKochCliques;
                mcsKochCliques = addKochCliques(isBondMatchFlag());
                for (AtomAtomMapping mapping : mcsKochCliques) {
                    Map<Integer, Integer> map = new TreeMap<Integer, Integer>();
                    map.putAll(mapping.getMappingsIndex());
                    mcsSeeds.add(map);
                }

                /*
                 * Copy UIT & Koch MCS solutions
                 */
                int solSize = 0;
                int counter = 0;
                List<Map<Integer, Integer>> cleanedMCSSeeds = new ArrayList<Map<Integer, Integer>>();
                if (!mcsSeeds.isEmpty()) {
                    for (Map<Integer, Integer> map : mcsSeeds) {
                        if (map.size() > solSize) {
                            solSize = map.size();
                            cleanedMCSSeeds.clear();
                            counter = 0;
                        }
                        if (!map.isEmpty()
                                && map.size() == solSize
                                && !hasClique(map, cleanedMCSSeeds)) {
                            cleanedMCSSeeds.add(counter, map);
                            counter++;
                        }
                    }
                }
                /*
                 * Copy VF MCS solutions
                 */
                for (Map<Integer, Integer> vfSolutions : allMCSCopy) {
                    Map<Integer, Integer> map = new TreeMap<Integer, Integer>();
                    map.putAll(vfSolutions);
                    if (!map.isEmpty()
                            && map.size() >= solSize
                            && !hasClique(map, cleanedMCSSeeds)) {
                        cleanedMCSSeeds.add(counter, map);
                        counter++;
                    }
                }

                allMCSCopy.clear();
                allLocalAtomAtomMapping.clear();
                /*
                 * Sort biggest clique to smallest
                 */
                Collections.sort(cleanedMCSSeeds, new Map2Comparator());
//                System.out.println("cleanedMCSSeeds " + cleanedMCSSeeds);
                extendCliquesWithMcGregor(cleanedMCSSeeds);

            }
        } catch (CDKException ex) {
            logger.error(Level.SEVERE, null, ex);
        } catch (IOException ex) {
            logger.error(Level.SEVERE, null, ex);
        }
        int solSize = 0;
        int counter = 0;
        if (!allLocalAtomAtomMapping.isEmpty()) {
            for (int i = 0; i < allLocalAtomAtomMapping.size(); i++) {
                AtomAtomMapping atomMCSMap = allLocalAtomAtomMapping.get(i);
                if (atomMCSMap.getCount() > solSize) {
                    solSize = atomMCSMap.getCount();
                    allAtomMCS.clear();
                    counter = 0;
                }
                if (!atomMCSMap.isEmpty()
                        && atomMCSMap.getCount() == solSize) {
                    allAtomMCS.add(counter, atomMCSMap);
                    counter++;
                }
            }
        }
    }

    /**
     * Constructor for an extended VF Algorithm for the MCS search
     *
     * @param source
     * @param target
     */
    public VF2MCS(IQueryAtomContainer source, IQueryAtomContainer target) {
        this(source, target, false, false);
    }

    /**
     * {@inheritDoc}
     *
     * @return
     */
    @Override
    @TestMethod("testGetAllAtomMapping")
    public synchronized List<AtomAtomMapping> getAllAtomMapping() {
        return Collections.unmodifiableList(allAtomMCS);
    }

    /**
     * {@inheritDoc}
     *
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
