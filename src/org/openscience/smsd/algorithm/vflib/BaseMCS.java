/* 
 * Copyright (C) 2009-2011  Syed Asad Rahman <asad@ebi.ac.uk>
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
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.Stack;
import java.util.TreeMap;
import java.util.logging.Level;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.isomorphism.matchers.IQueryAtomContainer;
import org.openscience.cdk.tools.ILoggingTool;
import org.openscience.cdk.tools.LoggingToolFactory;
import org.openscience.smsd.AtomAtomMapping;
import org.openscience.smsd.algorithm.mcgregor.McGregor;
import org.openscience.smsd.algorithm.mcsplus.BKKCKCF;
import org.openscience.smsd.algorithm.mcsplus.GenerateCompatibilityGraph;
import org.openscience.smsd.algorithm.rgraph.CDKRMapHandler;
import org.openscience.smsd.algorithm.vflib.interfaces.IMapper;
import org.openscience.smsd.algorithm.vflib.interfaces.INode;
import org.openscience.smsd.algorithm.vflib.interfaces.IQuery;
import org.openscience.smsd.algorithm.vflib.map.VFMCSMapper;
import org.openscience.smsd.algorithm.vflib.query.QueryCompiler;
import org.openscience.smsd.global.TimeOut;
import org.openscience.smsd.helper.FinalMappings;
import org.openscience.smsd.helper.MoleculeInitializer;
import org.openscience.smsd.tools.TimeManager;

/**
 * This class should be used to find MCS between source graph and target graph.
 *
 * First the algorithm runs VF lib {@link org.openscience.cdk.smsd.algorithm.vflib.map.VFMCSMapper} and reports MCS
 * between run source and target graphs. Then these solutions are extended using McGregor {@link org.openscience.cdk.smsd.algorithm.mcgregor.McGregor}
 * algorithm where ever required.
 *
 * @cdk.module smsd @cdk.githash
 *
 * @author Syed Asad Rahman <asad@ebi.ac.uk>
 */
public class BaseMCS extends MoleculeInitializer {

    private int countR = 0;
    private int countP = 0;
    private final IAtomContainer source;
    private final IAtomContainer target;
    private final TimeManager timeManager;
    private final boolean shouldMatchRings;
    private final boolean matchBonds;
    private final List<Map<INode, IAtom>> vfLibSolutions;
    final List<Map<Integer, Integer>> allMCSCopy;
    final List<AtomAtomMapping> allLocalAtomAtomMapping;
    private final static ILoggingTool Logger =
            LoggingToolFactory.createLoggingTool(BaseMCS.class);

    BaseMCS(IAtomContainer source, IAtomContainer target, boolean matchBonds, boolean shouldMatchRings) {
        this.allLocalAtomAtomMapping = new ArrayList<AtomAtomMapping>();
        this.timeManager = new TimeManager();
        this.allMCSCopy = new ArrayList<Map<Integer, Integer>>();
        this.shouldMatchRings = shouldMatchRings;
        this.matchBonds = matchBonds;
        this.vfLibSolutions = new ArrayList<Map<INode, IAtom>>();
        this.source = source;
        this.target = target;
        if (shouldMatchRings) {
            try {
                initializeMolecule(source);
                initializeMolecule(target);
            } catch (CDKException ex) {
                Logger.error(ex);
            }
        }

    }

    protected synchronized boolean hasClique(Map<Integer, Integer> cliqueMap, List<Map<Integer, Integer>> mapGlobal) {
        for (Map<Integer, Integer> storedMap : mapGlobal) {
            if (cliqueMap.size() < storedMap.size()) {
                return true;
            }
            if (cliqueMap.keySet().equals(storedMap.keySet())
                    && cliqueMap.values().equals(storedMap.values())) {
                return true;
            }
        }
        return false;
    }

    protected synchronized List<AtomAtomMapping> addKochCliques(boolean bondSensitiveFlag) {
//        System.out.println("addKochCliques ");
        List<AtomAtomMapping> allCliqueAtomMCS = new ArrayList<AtomAtomMapping>();
        try {
            GenerateCompatibilityGraph gcg = new GenerateCompatibilityGraph(source, target, bondSensitiveFlag, isMatchRings());
            List<Integer> comp_graph_nodes = gcg.getCompGraphNodes();

            List<Integer> cEdges = gcg.getCEgdes();
            List<Integer> dEdges = gcg.getDEgdes();

            BKKCKCF init = new BKKCKCF(comp_graph_nodes, cEdges, dEdges);
//                Koch init = new Koch(comp_graph_nodes, cEdges, dEdges);
            Stack<List<Integer>> maxCliqueSet = init.getMaxCliqueSet();
            //clear all the compatibility graph content
            gcg.clear();

            /*
             * Sort biggest clique to smallest
             */
            Collections.sort(maxCliqueSet, new Comparator<List<Integer>>() {

                @Override
                public int compare(List<Integer> a1, List<Integer> a2) {
                    return a2.size() - a1.size(); // assumes you want biggest to smallest
                }
            });

            while (!maxCliqueSet.empty()) {
                List<Integer> peek = maxCliqueSet.peek();
                AtomAtomMapping atomatomMapping = new AtomAtomMapping(source, target);

                for (Integer value : peek) {
                    int[] index = getIndex(value.intValue(), comp_graph_nodes);
                    Integer qIndex = index[0];
                    Integer tIndex = index[1];
                    if (qIndex != -1 && tIndex != -1) {
                        IAtom qAtom = source.getAtom(qIndex);
                        IAtom tAtom = target.getAtom(tIndex);
                        atomatomMapping.put(qAtom, tAtom);
                    } else {
                        try {
                            throw new CDKException("Atom index pointing to -1");
                        } catch (CDKException ex) {
                            Logger.error(Level.SEVERE, null, ex);
                        }
                    }
                }

                if (!atomatomMapping.isEmpty()) {
                    allCliqueAtomMCS.add(atomatomMapping);
                }
                maxCliqueSet.pop();
            }
        } catch (IOException ex) {
            Logger.error(Level.SEVERE, null, ex);
        }

        return allCliqueAtomMCS;
    }

    /**
     *
     * @param bondSensitiveFlag
     * @return
     */
    protected synchronized List<AtomAtomMapping> addUIT(boolean bondSensitiveFlag) {
        CDKRMapHandler rmap = new CDKRMapHandler();
        boolean rOnPFlag = false;
        List<Map<Integer, Integer>> sol = null;
        try {

            if (source.getAtomCount() > target.getAtomCount()) {
                rOnPFlag = true;
                rmap.calculateOverlapsAndReduce(source, target, bondSensitiveFlag, isMatchRings());
            } else {
                rOnPFlag = false;
                rmap.calculateOverlapsAndReduce(target, source, bondSensitiveFlag, isMatchRings());
            }
            sol = FinalMappings.getInstance().getFinalMapping();

        } catch (CDKException e) {
//            System.err.println("WARNING: graphContainer: most probably time out error ");
        }
        return setUITMappings(rOnPFlag, sol);
    }

    /*
     * Note: VF MCS will search for cliques which will match the types. Mcgregor will extend the cliques depending of
     * the bond type (sensitive and insensitive).
     */
    protected synchronized void addVFMatchesMappings() {
//        System.out.println("addVFMatchesMappings ");
        IQuery queryCompiler;
        IMapper mapper;

        if (!(source instanceof IQueryAtomContainer) && !(target instanceof IQueryAtomContainer)) {
            countR = getReactantMol().getAtomCount();
            countP = getProductMol().getAtomCount();
        }

        if (source instanceof IQueryAtomContainer) {
            queryCompiler = new QueryCompiler((IQueryAtomContainer) source).compile();
            mapper = new VFMCSMapper(queryCompiler);
            List<Map<INode, IAtom>> maps = mapper.getMaps(getProductMol());
            if (maps != null) {
                vfLibSolutions.addAll(maps);
            }
            setVFMappings(true, queryCompiler);

        } else if (countR <= countP) {
            queryCompiler = new QueryCompiler(this.source, isBondMatchFlag(), isMatchRings()).compile();
            mapper = new VFMCSMapper(queryCompiler);
            List<Map<INode, IAtom>> map = mapper.getMaps(this.target);
            if (map != null) {
                vfLibSolutions.addAll(map);
            }
            setVFMappings(true, queryCompiler);
        } else {
            queryCompiler = new QueryCompiler(this.target, isBondMatchFlag(), isMatchRings()).compile();
            mapper = new VFMCSMapper(queryCompiler);
            List<Map<INode, IAtom>> map = mapper.getMaps(this.source);
            if (map != null) {
                vfLibSolutions.addAll(map);
            }
            setVFMappings(false, queryCompiler);
        }
//        System.out.println("Sol count " + vfLibSolutions.size());
//        System.out.println("Sol size " + vfLibSolutions.iterator().next().size());
//        System.out.println("MCSSize " + vfMCSSize);
//        System.out.println("After Sol count " + allMCSCopy.size());

    }

    protected synchronized void extendCliquesWithMcGregor(List<Map<Integer, Integer>> refinedMCSSeeds) throws CDKException, IOException {
        List<List<Integer>> mappings = new ArrayList<List<Integer>>();
        boolean ROPFlag = true;
        for (Map<Integer, Integer> firstPassMappings : refinedMCSSeeds) {
            Map<Integer, Integer> extendMapping = new TreeMap<Integer, Integer>(firstPassMappings);
            McGregor mgit;
            if (source instanceof IQueryAtomContainer) {
                mgit = new McGregor((IQueryAtomContainer) source, target, mappings, isBondMatchFlag(), isMatchRings());
            } else {
                if (countR > countP) {
                    mgit = new McGregor(source, target, mappings, isBondMatchFlag(), isMatchRings());
                } else {
                    extendMapping.clear();
                    mgit = new McGregor(target, source, mappings, isBondMatchFlag(), isMatchRings());
                    ROPFlag = false;
                    for (Map.Entry<Integer, Integer> map : firstPassMappings.entrySet()) {
                        extendMapping.put(map.getValue(), map.getKey());
                    }
                }
            }
            //Start McGregor search
            mgit.startMcGregorIteration(mgit.getMCSSize(), extendMapping);
            mappings = mgit.getMappings();

            if (isTimeOut()) {
                System.err.println("\nMcGregorIteration timeout");
                break;
            }
        }
//        System.out.println("\nSol count after MG " + mappings.size());
        setMcGregorMappings(ROPFlag, mappings);
//        System.out.println("After set Sol count MG" + allMCS.size());
//        System.out.println("MCSSize " + vfMCSSize + "\n");
    }

    private synchronized void setVFMappings(boolean RONP, IQuery query) {
        /*
         * Sort biggest clique to smallest
         */
        Collections.sort(vfLibSolutions, new Map1Comparator());
        for (Map<INode, IAtom> solution : vfLibSolutions) {
            AtomAtomMapping atomatomMapping = new AtomAtomMapping(source, target);
            Map<Integer, Integer> indexindexMapping = new TreeMap<Integer, Integer>();

            for (INode node : solution.keySet()) {
                IAtom qAtom = null;
                IAtom tAtom = null;
                int qIndex = -1;
                int tIndex = -1;

                if (RONP) {
                    qAtom = query.getAtom(node);
                    tAtom = solution.get(node);
                    qIndex = source.getAtomNumber(qAtom);
                    tIndex = target.getAtomNumber(tAtom);
                } else {
                    tAtom = query.getAtom(node);
                    qAtom = solution.get(node);
                    qIndex = source.getAtomNumber(qAtom);
                    tIndex = target.getAtomNumber(tAtom);
                }

                if (qIndex != -1 && tIndex != -1) {
                    atomatomMapping.put(qAtom, tAtom);
                    indexindexMapping.put(qIndex, tIndex);
                } else {
                    try {
                        throw new CDKException("Atom index pointing to -1");
                    } catch (CDKException ex) {
                        Logger.error(Level.SEVERE, null, ex);
                    }
                }
            }

            if (!indexindexMapping.isEmpty()
                    && !hasClique(indexindexMapping, getLocalMCSSolution())) {
                getLocalAtomMCSSolution().add(atomatomMapping);
                getLocalMCSSolution().add(indexindexMapping);
            }
        }
    }

    private synchronized List<AtomAtomMapping> setUITMappings(boolean RONP, List<Map<Integer, Integer>> sol) {
        /*
         * Sort biggest clique to smallest
         */
        List<AtomAtomMapping> allCliqueAtomMCS = new ArrayList<AtomAtomMapping>();
        /*
         * Sort biggest clique to smallest
         */
        Collections.sort(sol, new Map2Comparator());
        for (Map<Integer, Integer> solution : sol) {
            AtomAtomMapping atomatomMapping = new AtomAtomMapping(source, target);

            for (Integer qAtomIndex : solution.keySet()) {
                IAtom qAtom = null;
                IAtom tAtom = null;
                Integer qIndex = -1;
                Integer tIndex = -1;

                if (RONP) {
                    qAtom = source.getAtom(qAtomIndex);
                    tAtom = target.getAtom(solution.get(qAtomIndex));
                    qIndex = source.getAtomNumber(qAtom);
                    tIndex = target.getAtomNumber(tAtom);
                } else {
                    tAtom = target.getAtom(qAtomIndex);
                    qAtom = source.getAtom(solution.get(qAtomIndex));
                    tIndex = target.getAtomNumber(tAtom);
                    qIndex = source.getAtomNumber(qAtom);
                }

                if (qIndex != -1 && tIndex != -1) {
                    atomatomMapping.put(qAtom, tAtom);
                } else {
                    try {
                        throw new CDKException("Atom index pointing to -1");
                    } catch (CDKException ex) {
                        Logger.error(Level.SEVERE, null, ex);
                    }
                }
            }

            if (!atomatomMapping.isEmpty()) {
                allCliqueAtomMCS.add(atomatomMapping);
            }
        }
        return allCliqueAtomMCS;
    }

    private synchronized void setMcGregorMappings(boolean RONP, List<List<Integer>> mappings) throws CDKException {
        int counter = 0;
        int solSize = 0;
        getLocalAtomMCSSolution().clear();
        getLocalMCSSolution().clear();
        for (List<Integer> mapping : mappings) {

            AtomAtomMapping atomatomMapping = new AtomAtomMapping(source, target);
            Map<Integer, Integer> indexindexMapping = new TreeMap<Integer, Integer>();
            for (int index = 0; index < mapping.size(); index += 2) {
                IAtom qAtom = null;
                IAtom tAtom = null;
                Integer qIndex = -1;
                Integer tIndex = -1;

                if (RONP) {
                    qAtom = source.getAtom(mapping.get(index));
                    tAtom = target.getAtom(mapping.get(index + 1));
                    qIndex = mapping.get(index);
                    tIndex = mapping.get(index + 1);
                } else {
                    qAtom = source.getAtom(mapping.get(index + 1));
                    tAtom = target.getAtom(mapping.get(index));
                    qIndex = mapping.get(index + 1);
                    tIndex = mapping.get(index);
                }

                if (qIndex != -1 && tIndex != -1) {
                    atomatomMapping.put(qAtom, tAtom);
                    indexindexMapping.put(qIndex, tIndex);
                } else {
                    throw new CDKException("Atom index pointing to NULL");
                }
            }
            if (indexindexMapping.size() > solSize) {
                solSize = indexindexMapping.size();
                getLocalAtomMCSSolution().clear();
                getLocalMCSSolution().clear();
                counter = 0;
            }
            if (!indexindexMapping.isEmpty()
                    && !hasClique(indexindexMapping, getLocalMCSSolution())
                    && indexindexMapping.size() == solSize) {
                getLocalAtomMCSSolution().add(counter, atomatomMapping);
                getLocalMCSSolution().add(counter, indexindexMapping);
                counter++;
            }
        }

    }

    private int[] getIndex(int cliqueIndex, List<Integer> comp_graph_nodes) {
        int[] v = new int[2];
        v[0] = -1;
        v[1] = -1;
        for (int i = 0; i < comp_graph_nodes.size(); i += 3) {
            if (cliqueIndex == comp_graph_nodes.get(i + 2)) {
                v[0] = comp_graph_nodes.get(i);
                v[1] = comp_graph_nodes.get(i + 1);
            }
        }
        return v;
    }

    protected synchronized IAtomContainer getReactantMol() {
        return source;
    }

    protected synchronized IAtomContainer getProductMol() {
        return target;
    }

    protected synchronized boolean isTimeOut() {
        if (getTimeout() > -1 && getTimeManager().getElapsedTimeInMinutes() > getTimeout()) {
            TimeOut.getInstance().setTimeOutFlag(true);
            return true;
        }
        return false;
    }

    /**
     * @return the shouldMatchRings
     */
    protected boolean isMatchRings() {
        return shouldMatchRings;
    }

    /**
     * @return the timeManager
     */
    protected synchronized TimeManager getTimeManager() {
        return timeManager;
    }

    /**
     * @return the timeout
     */
    protected synchronized double getTimeout() {
        return TimeOut.getInstance().getVFTimeout();
    }

    /**
     * @return the shouldMatchBonds
     */
    protected synchronized boolean isBondMatchFlag() {
        return matchBonds;
    }

    /**
     * @return the allMCSCopy
     */
    private synchronized List<Map<Integer, Integer>> getLocalMCSSolution() {
        return Collections.synchronizedList(allMCSCopy);
    }

    /**
     * @return the allLocalAtomAtomMapping
     */
    private synchronized List<AtomAtomMapping> getLocalAtomMCSSolution() {
        return Collections.synchronizedList(allLocalAtomAtomMapping);
    }

    protected synchronized boolean isExtensionRequired(List<Map<Integer, Integer>> mcsSeeds) {
        int maxSize = 0;
        for (Map<Integer, Integer> map : mcsSeeds) {
            if (map.size() > maxSize) {
                maxSize = map.size();
            }
        }
        return this.source.getAtomCount() > maxSize && this.target.getAtomCount() > maxSize ? true : false;
    }

    protected synchronized boolean isExtensionRequired() {
        int commonAtomCount = checkCommonAtomCount(getReactantMol(), getProductMol());
        int maxSize = 0;
        for (Map<Integer, Integer> map : allMCSCopy) {
            if (map.size() > maxSize) {
                maxSize = map.size();
            }
        }
        if (commonAtomCount > maxSize) {
            return true;
        }
        return false;
    }

    private synchronized int checkCommonAtomCount(IAtomContainer reactantMolecule, IAtomContainer productMolecule) {
        ArrayList<String> atoms = new ArrayList<String>();
        for (int i = 0; i < reactantMolecule.getAtomCount(); i++) {
            atoms.add(reactantMolecule.getAtom(i).getSymbol());
        }
        int common = 0;
        for (int i = 0; i < productMolecule.getAtomCount(); i++) {
            String symbol = productMolecule.getAtom(i).getSymbol();
            if (atoms.contains(symbol)) {
                atoms.remove(symbol);
                common++;
            }
        }
        return common;
    }

    protected synchronized boolean isEnrichmentRequired() {
        int maxSize = 0;
        for (Map<Integer, Integer> map : allMCSCopy) {
            if (map.size() > maxSize) {
                maxSize = map.size();
            }
        }
        return this.source.getAtomCount() == maxSize || this.target.getAtomCount() == maxSize ? true : false;
    }

    public class Map1Comparator implements Comparator<Map<INode, IAtom>> {

        /**
         *
         * @param object1
         * @param object2
         * @return
         */
        @Override
        public int compare(Map<INode, IAtom> object1, Map<INode, IAtom> object2) {

            Integer a1 = (Integer) ((Map) object1).size();
            Integer a2 = (Integer) ((Map) object2).size();
            return a2 - a1; // assumes you want biggest to smallest;

        }
    }

    public class Map2Comparator implements Comparator<Map<Integer, Integer>> {

        /**
         *
         * @param object1
         * @param object2
         * @return
         */
        @Override
        public int compare(Map<Integer, Integer> object1, Map<Integer, Integer> object2) {

            Integer a1 = (Integer) ((Map) object1).size();
            Integer a2 = (Integer) ((Map) object2).size();
            return a2 - a1; // assumes you want biggest to smallest;

        }
    }
}
