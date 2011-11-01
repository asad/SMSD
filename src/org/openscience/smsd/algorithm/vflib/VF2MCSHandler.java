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
import java.util.List;
import java.util.Map;
import java.util.Stack;
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
import org.openscience.smsd.algorithm.mcgregor.McGregor;
import org.openscience.smsd.algorithm.mcsplus.BKKCKCF;
import org.openscience.smsd.algorithm.mcsplus.GenerateCompatibilityGraph;
import org.openscience.smsd.algorithm.vflib.interfaces.IMapper;
import org.openscience.smsd.algorithm.vflib.interfaces.INode;
import org.openscience.smsd.algorithm.vflib.interfaces.IQuery;
import org.openscience.smsd.algorithm.vflib.map.VFMCSMapper;
import org.openscience.smsd.algorithm.vflib.query.QueryCompiler;
import org.openscience.smsd.global.TimeOut;
import org.openscience.smsd.interfaces.AbstractMCSAlgorithm;
import org.openscience.smsd.interfaces.IMCSBase;
import org.openscience.smsd.tools.TimeManager;

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
public class VF2MCSHandler extends AbstractMCSAlgorithm implements IMCSBase {

    private List<AtomAtomMapping> allAtomMCS = null;
    private List<AtomAtomMapping> allAtomMCSCopy = null;
    private List<Map<Integer, Integer>> allMCS = null;
    private List<Map<Integer, Integer>> allMCSCopy = null;
    private List<Map<INode, IAtom>> vfLibSolutions = null;
    private IAtomContainer source = null;
    private IAtomContainer target = null;
    private int vfMCSSize = -1;
    private boolean matchBonds;
    private int countR = 0;
    private int countP = 0;
    private final static ILoggingTool Logger =
            LoggingToolFactory.createLoggingTool(VF2MCSHandler.class);
    private TimeManager timeManager = null;
    private boolean shouldMatchRings;

    /**
     * @return the timeout
     */
    protected synchronized double getTimeout() {
        return TimeOut.getInstance().getVFTimeout();
    }

    /**
     * @return the timeManager
     */
    protected synchronized TimeManager getTimeManager() {
        return timeManager;
    }

    /**
     * @param aTimeManager the timeManager to set
     */
    public synchronized void setTimeManager(TimeManager aTimeManager) {
        TimeOut.getInstance().setTimeOutFlag(false);
        timeManager = aTimeManager;
    }

    /**
     * Constructor for an extended VF Algorithm for the MCS search
     */
    public VF2MCSHandler() {
        allAtomMCS = new ArrayList<AtomAtomMapping>();
        allAtomMCSCopy = new ArrayList<AtomAtomMapping>();
        allMCS = new ArrayList<Map<Integer, Integer>>();
        allMCSCopy = new ArrayList<Map<Integer, Integer>>();
    }

    /**
     *{@inheritDoc}
     *
     * @param matchBonds 
     */
    @Override
    @TestMethod("testSearchMCS")
    public synchronized void searchMCS(boolean matchBonds, boolean shouldMatchRings) {
        setTimeManager(new TimeManager());
        setBondMatchFlag(matchBonds);
        this.setMatchRings(shouldMatchRings);

        if (isMatchRings()) {
            try {
                initializeMolecule(source);
                initializeMolecule(target);
            } catch (CDKException ex) {
                Logger.error(Level.SEVERE, null, ex);
            }
        }

        searchVFMCSMappings();
        boolean flag = isExtensionFeasible();
        if (flag && !vfLibSolutions.isEmpty()) {
            try {
                GenerateCompatibilityGraph gcg = new GenerateCompatibilityGraph(source, target, true, isMatchRings());
                List<Integer> comp_graph_nodes = gcg.getCompGraphNodes();

                List<Integer> cEdges = gcg.getCEgdes();
                List<Integer> dEdges = gcg.getDEgdes();

                BKKCKCF init = new BKKCKCF(comp_graph_nodes, cEdges, dEdges);
//                Koch init = new Koch(comp_graph_nodes, cEdges, dEdges);
                Stack<List<Integer>> maxCliqueSet = init.getMaxCliqueSet();
                //clear all the compatibility graph content
                gcg.clear();
                int counter = 0;

                while (!maxCliqueSet.empty()) {
                    List<Integer> peek = maxCliqueSet.peek();
                    Map<Integer, Integer> mapping = new TreeMap<Integer, Integer>();
                    AtomAtomMapping atomatomMapping = new AtomAtomMapping(source, target);

                    for (Integer value : peek) {
                        int[] index = getIndex(value.intValue(), comp_graph_nodes);
                        Integer qIndex = index[0];
                        Integer tIndex = index[1];
                        if (qIndex != -1 && tIndex != -1) {
                            IAtom qAtom = source.getAtom(qIndex);
                            IAtom tAtom = target.getAtom(tIndex);
                            atomatomMapping.put(qAtom, tAtom);
                            mapping.put(qIndex, tIndex);
                        } else {
                            try {
                                throw new CDKException("Atom index pointing to -1");
                            } catch (CDKException ex) {
                                Logger.error(Level.SEVERE, null, ex);
                            }
                        }
                    }

                    if (!atomatomMapping.isEmpty() && !hasMap(mapping, allMCSCopy)) {
//                        System.out.println("\nvfMCSSize: " + vfMCSSize);
                        allAtomMCSCopy.add(counter, atomatomMapping);
                        allMCSCopy.add(counter, mapping);
                        counter++;
                    }
                    maxCliqueSet.pop();
                }

                searchMcGregorMapping();
            } catch (CDKException ex) {
                Logger.error(Level.SEVERE, null, ex);
            } catch (IOException ex) {
                Logger.error(Level.SEVERE, null, ex);
            }
        } else if (!allAtomMCSCopy.isEmpty()) {
            allAtomMCS.addAll(allAtomMCSCopy);
            allMCS.addAll(allMCSCopy);
        }
    }

    private synchronized boolean isExtensionFeasible() {
        int commonAtomCount = checkCommonAtomCount(getReactantMol(), getProductMol());
        if (commonAtomCount > vfMCSSize) {
            return true;
        }
        return false;
    }

    /** {@inheritDoc}
     *
     * Set the VFLib MCS software
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
    public void set(IQueryAtomContainer source, IAtomContainer target) {
        this.source = source;
        this.target = target;
    }

    private boolean hasMap(Map<Integer, Integer> maps, List<Map<Integer, Integer>> mapGlobal) {
        for (Map<Integer, Integer> test : mapGlobal) {
            if (test.size() > maps.size()) {
                return true;
            }
            if (test.equals(maps)) {
                return true;
            }
        }
        return false;
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

    /*
     * Note: VF MCS will search for cliques which will match the types.
     * Mcgregor will extend the cliques depending of the bond type 
     * (sensitive and insensitive).
     */
    private synchronized void searchVFMCSMappings() {
//        System.out.println("searchVFMCSMappings ");
        IQuery queryCompiler = null;
        IMapper mapper = null;

        if (!(source instanceof IQueryAtomContainer) && !(target instanceof IQueryAtomContainer)) {
            countR = getReactantMol().getAtomCount();
            countP = getProductMol().getAtomCount();
        }

        vfLibSolutions = new ArrayList<Map<INode, IAtom>>();
        if (source instanceof IQueryAtomContainer) {
            queryCompiler = new QueryCompiler((IQueryAtomContainer) source).compile();
            mapper = new VFMCSMapper(queryCompiler);
            List<Map<INode, IAtom>> maps = mapper.getMaps(getProductMol());
            if (maps != null) {
                vfLibSolutions.addAll(maps);
            }
            setVFMappings(true, queryCompiler);
        } else if (countR <= countP) {
//            queryCompiler = new QueryCompiler(this.source, isBondMatchFlag()).compile();
            queryCompiler = new QueryCompiler(this.source, true, isMatchRings()).compile();
            mapper = new VFMCSMapper(queryCompiler);
            List<Map<INode, IAtom>> maps = mapper.getMaps(getProductMol());
            if (maps != null) {
                vfLibSolutions.addAll(maps);
            }
            setVFMappings(true, queryCompiler);
        } else {
//            queryCompiler = new QueryCompiler(getProductMol(), isBondMatchFlag()).compile();
            queryCompiler = new QueryCompiler(getProductMol(), true, isMatchRings()).compile();
            mapper = new VFMCSMapper(queryCompiler);
            List<Map<INode, IAtom>> maps = mapper.getMaps(getReactantMol());
            if (maps != null) {
                vfLibSolutions.addAll(maps);
            }
            setVFMappings(false, queryCompiler);
        }
//        System.out.println("Sol count " + vfLibSolutions.size());
//        System.out.println("Sol size " + vfLibSolutions.iterator().next().size());
//        System.out.println("MCSSize " + vfMCSSize);
//        System.out.println("After Sol count " + allMCSCopy.size());

    }

    private synchronized void searchMcGregorMapping() throws CDKException, IOException {
        List<List<Integer>> mappings = new ArrayList<List<Integer>>();
        boolean ROPFlag = true;
        for (Map<Integer, Integer> firstPassMappings : allMCSCopy) {
            Map<Integer, Integer> extendMapping = new TreeMap<Integer, Integer>(firstPassMappings);
            McGregor mgit = null;
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
            mgit = null;

            if (isTimeOut()) {
                System.err.println("\nVFLibMCS hit by timeout in McGregor");
                break;
            }
        }
//        System.out.println("\nSol count after MG" + mappings.size());
        setMcGregorMappings(ROPFlag, mappings);
        vfMCSSize = vfMCSSize / 2;
//        System.out.println("After set Sol count MG" + allMCS.size());
//        System.out.println("MCSSize " + vfMCSSize + "\n");
    }

    private synchronized void setVFMappings(boolean RONP, IQuery query) {
        int counter = 0;
        for (Map<INode, IAtom> solution : vfLibSolutions) {
            AtomAtomMapping atomatomMapping = new AtomAtomMapping(source, target);
            Map<Integer, Integer> indexindexMapping = new TreeMap<Integer, Integer>();

            for (Map.Entry<INode, IAtom> mapping : solution.entrySet()) {
                IAtom qAtom = null;
                IAtom tAtom = null;
                Integer qIndex = -1;
                Integer tIndex = -1;

                if (RONP) {
                    qAtom = query.getAtom(mapping.getKey());
                    tAtom = mapping.getValue();
                    qIndex = getReactantMol().getAtomNumber(qAtom);
                    tIndex = getProductMol().getAtomNumber(tAtom);
                } else {
                    tAtom = query.getAtom(mapping.getKey());
                    qAtom = mapping.getValue();
                    qIndex = getReactantMol().getAtomNumber(qAtom);
                    tIndex = getProductMol().getAtomNumber(tAtom);
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
            if (indexindexMapping.size() > vfMCSSize) {
                vfMCSSize = indexindexMapping.size();
                allAtomMCSCopy.clear();
                allMCSCopy.clear();
                counter = 0;
            }
            if (!atomatomMapping.isEmpty() && !hasMap(indexindexMapping, allMCSCopy)
                    && indexindexMapping.size() == vfMCSSize) {
//                System.out.println("\nvfMCSSize: " + vfMCSSize);
                allAtomMCSCopy.add(counter, atomatomMapping);
                allMCSCopy.add(counter, indexindexMapping);
                counter++;
            }
        }
//        System.out.println("After set allMCSCopy " + allMCSCopy);
    }

    private synchronized void setMcGregorMappings(boolean RONP, List<List<Integer>> mappings) throws CDKException {
        int counter = 0;
        for (List<Integer> mapping : mappings) {

            AtomAtomMapping atomatomMapping = new AtomAtomMapping(source, target);
            Map<Integer, Integer> indexindexMapping = new TreeMap<Integer, Integer>();
            for (int index = 0; index < mapping.size(); index += 2) {
                IAtom qAtom = null;
                IAtom tAtom = null;
                Integer qIndex = -1;
                Integer tIndex = -1;

                if (RONP) {
                    qAtom = getReactantMol().getAtom(mapping.get(index));
                    tAtom = getProductMol().getAtom(mapping.get(index + 1));

                    qIndex = mapping.get(index);
                    tIndex = mapping.get(index + 1);
                } else {
                    qAtom = getReactantMol().getAtom(mapping.get(index + 1));
                    tAtom = getProductMol().getAtom(mapping.get(index));
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
            if (indexindexMapping.size() > vfMCSSize) {
                vfMCSSize = indexindexMapping.size();
                allAtomMCS.clear();
                allMCS.clear();
                counter = 0;
            }
            if (!atomatomMapping.isEmpty() && !hasMap(indexindexMapping, allMCS)
                    && (indexindexMapping.size()) == vfMCSSize) {
                allAtomMCS.add(counter, atomatomMapping);
                allMCS.add(counter, indexindexMapping);
                counter++;
            }
        }
    }

    /**
     * @return the shouldMatchBonds
     */
    public synchronized boolean isBondMatchFlag() {
        return matchBonds;
    }

    /**
     * @param shouldMatchBonds the shouldMatchBonds to set
     */
    public synchronized void setBondMatchFlag(boolean shouldMatchBonds) {
        this.matchBonds = shouldMatchBonds;
    }

    private synchronized IAtomContainer getReactantMol() {
        return source;
    }

    private synchronized IAtomContainer getProductMol() {
        return target;
    }

    public synchronized boolean isTimeOut() {
        if (getTimeout() > -1 && getTimeManager().getElapsedTimeInMinutes() > getTimeout()) {
            TimeOut.getInstance().setTimeOutFlag(true);
            return true;
        }
        return false;
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

    /**
     * @return the shouldMatchRings
     */
    public boolean isMatchRings() {
        return shouldMatchRings;
    }

    /**
     * @param shouldMatchRings the shouldMatchRings to set
     */
    public void setMatchRings(boolean shouldMatchRings) {
        this.shouldMatchRings = shouldMatchRings;
    }
}
