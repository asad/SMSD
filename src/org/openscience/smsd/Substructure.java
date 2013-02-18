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
package org.openscience.smsd;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.logging.Level;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.isomorphism.matchers.IQueryAtomContainer;
import org.openscience.cdk.tools.ILoggingTool;
import org.openscience.cdk.tools.LoggingToolFactory;
import org.openscience.smsd.algorithm.single.SingleMappingHandler;
import org.openscience.smsd.algorithm.vflib.VF2Sub;
import org.openscience.smsd.algorithm.vflib.substructure.VF2;

/**
 * This is an ultra fast method to report if query is a substructure for target
 * molecule. If this case is true then it returns only all mapping.
 *
 * This is much faster than {@link
 * org.openscience.cdk.smsd.algorithm.vflib.substructure} class as it only
 * reports first match and backtracks.
 *
 * This class should only be used to report if a query graph is a substructure
 * of the target graph.
 *
 *  *
 * <p>An example for <b>Substructure search</b>:</p> <font color="#003366">
 * <pre>
 * SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
 * IAtomContainer query = sp.parseSmiles("CC");
 * IAtomContainer target = sp.parseSmiles("C1CCC12CCCC2");
 * Substructure smsd = new Substructure(query, target, true, false, true);
 * Assert.assertTrue(smsd.isSubgraph());
 * Assert.assertEquals(18, smsd.getAllAtomMapping().size());
 *
 * IQueryAtomContainer queryContainer = QueryAtomContainerCreator.createSymbolAndBondOrderQueryContainer(query);
 * smsd = new Substructure(queryContainer, target, true);
 * Assert.assertTrue(smsd.isSubgraph());
 *
 *  </pre> </font>
 *
 *
 * @cdk.module smsd
 * @cdk.githash
 * @author Syed Asad Rahman <asad@ebi.ac.uk>
 */
public final class Substructure extends BaseMapping {

    private int vfMappingSize = -1;
    private final ILoggingTool Logger =
            LoggingToolFactory.createLoggingTool(Substructure.class);

    /**
     * Constructor for VF Substructure Algorithm
     *
     * @param query
     * @param target
     * @param shouldMatchBonds Match bond types (i.e. double to double etc)
     * @param matchRings Match ring atoms and ring size
     * @param findAllSubgraph report all subgraphs
     * @throws CDKException
     */
    public Substructure(
            IAtomContainer query,
            IAtomContainer target,
            boolean shouldMatchBonds,
            boolean matchRings,
            boolean findAllSubgraph) throws CDKException {
        super(shouldMatchBonds, matchRings, query, target);
        if (findAllSubgraph) {
            setSubgraph(findSubgraphs());
        } else {
            setSubgraph(findSubgraph());
        }
    }

    /**
     * Constructor for VF Substructure Algorithm
     *
     * @param query
     * @param target
     * @param findAllSubgraph report all subgraphs
     * @throws CDKException
     */
    public Substructure(
            IQueryAtomContainer query,
            IAtomContainer target,
            boolean findAllSubgraph) throws CDKException {
        super(true, true, query, target);
        if (findAllSubgraph) {
            setSubgraph(findSubgraphs());
        } else {
            setSubgraph(findSubgraph());
        }
    }

    private synchronized boolean hasMap(AtomAtomMapping map, List<AtomAtomMapping> mapGlobal) {
        for (AtomAtomMapping test : mapGlobal) {
            if (test.equals(map)) {
                return true;
            }
        }
        return false;
    }

    /**
     * Returns true if query is a subgraph of target molecule
     *
     * @return
     * @throws CDKException
     */
    private synchronized boolean findSubgraph() throws CDKException {
        boolean isSubgraph = false;

        if ((getTargetContainer() == null) || (getQueryContainer() == null)) {
            throw new CDKException("Query or Target molecule is not initialized (NULL)");
        }

        if (getQueryContainer().getAtomCount() == 1 || getTargetContainer().getAtomCount() == 1) {
            singleMapping(isMatchBonds());
        } else {
            if (getQueryContainer().getAtomCount() > getTargetContainer().getAtomCount()) {
                return false;
            }
            VF2 mapper = null;
            List<AtomAtomMapping> mappingsVF2 = new ArrayList<AtomAtomMapping>();
            if (getQueryContainer() instanceof IQueryAtomContainer) {
                mapper = new VF2((IQueryAtomContainer) getQueryContainer(), getTargetContainer());
            } else {
                mapper = new VF2(getQueryContainer(), getTargetContainer(), isMatchBonds(), isMatchRings());
            }
            isSubgraph = mapper.isSubgraph();
            List<AtomAtomMapping> atomMappings = mapper.getAllAtomMapping();
//            System.out.println("Mapping Size " + atomMappings.getCount());
            if (isSubgraph) {
                mappingsVF2.addAll(atomMappings);
            } else {
                return false;
            }
            setVFMappings(mappingsVF2);
        }
        return isSubgraph;
    }

    /**
     * Returns true if query is a subgraph of target molecule
     *
     * @return
     * @throws CDKException
     */
    private synchronized boolean findSubgraphs() throws CDKException {
        boolean isSubgraph = false;

        if ((getTargetContainer() == null) || (getQueryContainer() == null)) {
            throw new CDKException("Query or Target molecule is not initialized (NULL)");
        }

        if (getQueryContainer().getAtomCount() == 1 || getTargetContainer().getAtomCount() == 1) {
            singleMapping(isMatchBonds());
        } else {
            if (getQueryContainer().getAtomCount() > getTargetContainer().getAtomCount()) {
                return false;
            } else {
                List<AtomAtomMapping> mappingsVF2 = new ArrayList<AtomAtomMapping>();
                VF2Sub mapper = null;
                if (getQueryContainer() instanceof IQueryAtomContainer) {
                    mapper = new VF2Sub((IQueryAtomContainer) getQueryContainer(), getTargetContainer());
                } else {
                    mapper = new VF2Sub(getQueryContainer(), getTargetContainer(), isMatchBonds(), isMatchRings());
                }
                isSubgraph = mapper.isSubgraph();
                List<AtomAtomMapping> atomMappings = mapper.getAllAtomMapping();
//                    System.out.println("Mapping Size " + atomMapping.getCount());
                if (isSubgraph) {
                    mappingsVF2.addAll(atomMappings);
                } else {
                    return false;
                }
                setVFMappings(mappingsVF2);
            }
        }
        return isSubgraph;
    }

    private synchronized void setVFMappings(List<AtomAtomMapping> mappingsVF2) {
        int counter = 0;
        for (AtomAtomMapping solution : mappingsVF2) {
            AtomAtomMapping atomatomMapping = new AtomAtomMapping(getQueryContainer(), getTargetContainer());
            if (solution.getCount() > vfMappingSize) {
                this.vfMappingSize = solution.getCount();
                counter = 0;
            }
            for (Map.Entry<IAtom, IAtom> mapping : solution.getMappings().entrySet()) {
                IAtom qAtom = null;
                IAtom tAtom = null;

                qAtom = mapping.getKey();
                tAtom = mapping.getValue();

                if (qAtom != null && tAtom != null) {
                    atomatomMapping.put(qAtom, tAtom);
                } else {
                    try {
                        throw new CDKException("Atom index pointing to NULL");
                    } catch (CDKException ex) {
                        Logger.error(Level.SEVERE, null, ex);
                    }
                }
            }
            if (!atomatomMapping.isEmpty() && !hasMap(atomatomMapping, getMCSList())
                    && atomatomMapping.getCount() == vfMappingSize) {
                getMCSList().add(counter, atomatomMapping);
                counter++;
            }
        }
    }

    private synchronized void singleMapping(boolean shouldMatchBonds) {
        SingleMappingHandler mcs = null;
        if (!(getQueryContainer() instanceof IQueryAtomContainer) && !(getTargetContainer() instanceof IQueryAtomContainer)) {
            mcs = new SingleMappingHandler(getQueryContainer(), getTargetContainer(), shouldMatchBonds, isMatchRings());
        } else {
            mcs = new SingleMappingHandler((IQueryAtomContainer) getQueryContainer(), getTargetContainer());
        }
        getMCSList().addAll(mcs.getAllAtomMapping());
    }
}
