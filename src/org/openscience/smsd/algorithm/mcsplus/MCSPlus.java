/* Copyright (C) 2009-2013  Syed Asad Rahman <asad@ebi.ac.uk>
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
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA.
 */
package org.openscience.smsd.algorithm.mcsplus;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Stack;
import java.util.TreeMap;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.openscience.cdk.annotations.TestClass;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.smsd.algorithm.mcgregor.McGregor;
import org.openscience.smsd.tools.IterationManager;

/**
 * This class handles MCS plus algorithm which is a combination of c-clique algorithm and McGregor algorithm.
 *
 * @cdk.module smsd
 * @cdk.githash
 *
 * @author Syed Asad Rahman <asad@ebi.ac.uk>
 */
@TestClass("org.openscience.cdk.smsd.SMSDBondSensitiveTest")
public class MCSPlus {

    private boolean shouldMatchRings;
    private boolean shouldMatchBonds;
    private boolean timeout = false;

    /**
     * Default constructor added
     */
    public MCSPlus() {
    }
    private IterationManager iterationManager = null;

    /**
     * @return the timeout
     */
    public synchronized boolean isTimeout() {
        return timeout;
    }

    /**
     * @return the iterationManager
     */
    private IterationManager getIterationManager() {
        return iterationManager;
    }

    /**
     * @param iterationManager the iterationManager to set
     */
    private void setIterationManager(IterationManager iterationManager) {
        this.iterationManager = iterationManager;
    }

    /**
     *
     * @param ac1
     * @param ac2
     * @param shouldMatchBonds
     * @param shouldMatchRings
     * @return
     * @throws CDKException
     */
    protected List<List<Integer>> getOverlaps(
            IAtomContainer ac1,
            IAtomContainer ac2,
            boolean shouldMatchBonds,
            boolean shouldMatchRings)
            throws CDKException {

        this.shouldMatchBonds = shouldMatchBonds;
        this.shouldMatchRings = shouldMatchRings;

        List<List<Integer>> extendMappings = null;

//        System.out.println("ac1 : " + ac1.getAtomCount());
//        System.out.println("ac2 : " + ac2.getAtomCount());
        setIterationManager(new IterationManager((ac1.getAtomCount() + ac2.getAtomCount())));
        try {
            GenerateCompatibilityGraph gcg = new GenerateCompatibilityGraph(ac1, ac2, isMatchBonds(), isMatchRings(), false);
            List<Integer> comp_graph_nodes = gcg.getCompGraphNodes();

            List<Integer> cEdges = gcg.getCEgdes();
            List<Integer> dEdges = gcg.getDEgdes();
//
//            System.out.println("**************************************************");
//            System.out.println("C_edges: " + cEdges.size());
//            System.out.println("D_edges: " + dEdges.size());
//            System.out.println("comp_graph_nodes: " + comp_graph_nodes);
            BKKCKCF init = new BKKCKCF(comp_graph_nodes, cEdges, dEdges);
            Stack<List<Integer>> maxCliqueSet = new Stack<>();
            maxCliqueSet.addAll(init.getMaxCliqueSet());

//            System.out.println("Max_Cliques_Set: " + maxCliqueSet);
//            System.out.println("Best Clique Size: " + init.getBestCliqueSize());
//            System.out.println("**************************************************");
            //clear all the compatibility graph content
            gcg.clear();
            List<Map<Integer, Integer>> mappings = new ArrayList<>();

            while (!maxCliqueSet.empty()) {
                Map<Integer, Integer> indexindexMapping;
                indexindexMapping = ExactMapping.extractMapping(comp_graph_nodes, maxCliqueSet.peek());
                if (indexindexMapping != null) {
                    mappings.add(indexindexMapping);
                }
                maxCliqueSet.pop();
            }
//            System.out.println("mappings: " + mappings.size());
            extendMappings = searchMcGregorMapping(ac1, ac2, mappings);
//            int size = !extendMappings.isEmpty() ? (extendMappings.size() / 2) : 0;
//            System.out.println("extendMappings: " + size);
        } catch (IOException ex) {
            Logger.getLogger(MCSPlus.class.getName()).log(Level.SEVERE, null, ex);
        }

        return extendMappings;
    }

    private List<List<Integer>> searchMcGregorMapping(
            IAtomContainer ac1,
            IAtomContainer ac2,
            List<Map<Integer, Integer>> allMCSCopy) throws IOException {

        List<List<Integer>> cliques = new ArrayList<>();

        boolean ROPFlag = true;
        for (Map<Integer, Integer> firstPassMappings : allMCSCopy) {
            Map<Integer, Integer> extendMapping = new TreeMap<>(firstPassMappings);
            McGregor mgit;
            if (ac1.getAtomCount() > ac2.getAtomCount()) {
                mgit = new McGregor(ac1, ac2, cliques, isMatchBonds(), isMatchRings());
            } else {
                extendMapping.clear();
                ROPFlag = false;
                for (Map.Entry<Integer, Integer> map : firstPassMappings.entrySet()) {
                    extendMapping.put(map.getValue(), map.getKey());
                }
                mgit = new McGregor(ac2, ac1, cliques, isMatchBonds(), isMatchRings());
            }
//            System.out.println("\nStart McGregor search");
            //Start McGregor search
            mgit.startMcGregorIteration(mgit.getMCSSize(), extendMapping);
            cliques = mgit.getMappings();
//            System.out.println("\nSol count after MG " + cliques.size());
            if (checkTimeout()) {
                break;
            }
        }
        List<List<Integer>> finalMappings = setMcGregorMappings(ROPFlag, cliques);
//        System.out.println("After set Sol count MG " + finalMappings.size());
        return finalMappings;
    }

    private List<List<Integer>> setMcGregorMappings(
            boolean RONP,
            List<List<Integer>> mappings) {
        int counter = 0;
        int mcsSize = 0;
        List<List<Integer>> finalMappings = new ArrayList<>();
        for (List<Integer> mapping : mappings) {
            List<Integer> indexindexMapping = new ArrayList<>();
            for (int index = 0; index < mapping.size(); index += 2) {
                Integer qIndex;
                Integer tIndex;

                if (RONP) {
                    qIndex = mapping.get(index);
                    tIndex = mapping.get(index + 1);
                } else {
                    qIndex = mapping.get(index + 1);
                    tIndex = mapping.get(index);
                }

                if (qIndex != null && tIndex != null) {
                    indexindexMapping.add(qIndex);
                    indexindexMapping.add(tIndex);
                }
            }
            if (!indexindexMapping.isEmpty() && indexindexMapping.size() > mcsSize) {
                mcsSize = indexindexMapping.size();
                finalMappings.clear();
                counter = 0;
            }
            if (!indexindexMapping.isEmpty() && !finalMappings.contains(indexindexMapping)
                    && (indexindexMapping.size()) == mcsSize) {
                finalMappings.add(counter, indexindexMapping);
                counter++;
            }
        }
        return finalMappings;
    }

    private synchronized boolean checkTimeout() {
        if (getIterationManager().isMaxIteration()) {
            this.timeout = true;
//            System.out.println("MCS+ iterations " + getIterationManager().getCounter());
            return true;
        }
        getIterationManager().increment();
        return false;
    }

    /**
     * @return the shouldMatchRings
     */
    public boolean isMatchRings() {
        return shouldMatchRings;
    }

    /**
     * @return the shouldMatchBonds
     */
    public boolean isMatchBonds() {
        return shouldMatchBonds;
    }
}
