/*
 *
 *
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
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA.
 * 
 * 
 ** Copyright (C) 2009-2011 Kyle Lutz <kyle.r.lutz@gmail.com>
 **
 ** This file is part of chemkit. For more information see
 ** <http://www.chemkit.org>.
 **
 ** chemkit is free software: you can redistribute it and/or modify
 ** it under the terms of the GNU Lesser General Public License as published by
 ** the Free Software Foundation, either version 3 of the License, or
 ** (at your option) any later version.
 **
 ** chemkit is distributed in the hope that it will be useful,
 ** but WITHOUT ANY WARRANTY; without even the implied warranty of
 ** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 ** GNU Lesser General Public License for more details.
 **
 ** You should have received a copy of the GNU Lesser General Public License
 ** along with chemkit. If not, see <http://www.gnu.org/licenses/>.
 **
 ******************************************************************************/
package org.openscience.smsd.algorithm.vflib.substructure;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.logging.Level;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.isomorphism.matchers.IQueryAtomContainer;
import org.openscience.cdk.tools.ILoggingTool;
import org.openscience.cdk.tools.LoggingToolFactory;
import org.openscience.smsd.AtomAtomMapping;
import org.openscience.smsd.interfaces.AbstractSubGraph;
import org.openscience.smsd.interfaces.IMCSBase;

/**
 * This class finds mapping states between query and target
 * molecules.
 * 
 * @cdk.module smsd
 * @cdk.githash
 * @author Syed Asad Rahman <asad@ebi.ac.uk>
 */
public final class VF2 extends AbstractSubGraph implements IMCSBase {

    private List<AtomAtomMapping> allAtomMCS = null;
    private List<Map<Integer, Integer>> allMCS = null;
    private IAtomContainer source = null;
    private IAtomContainer target = null;
    private final boolean shouldMatchRings;
    private final ILoggingTool Logger =
            LoggingToolFactory.createLoggingTool(VF2.class);
    private final boolean shouldMatchBonds;

    /**
     * Constructor for an extended VF Algorithm for the MCS search
     * @param shouldMatchBonds
     * @param shouldMatchRings  
     */
    public VF2(boolean shouldMatchBonds, boolean shouldMatchRings) {
        this.shouldMatchRings = shouldMatchRings;
        this.shouldMatchBonds = shouldMatchBonds;
        allAtomMCS = new ArrayList<AtomAtomMapping>();
        allMCS = new ArrayList<Map<Integer, Integer>>();
    }

    /** 
     * The isomorphism method returns an isomorphism between two molecular
     *  graphs using the VF2Automorphism algorithm. This can be used for finding both
     *  graph-graph isomorphisms and graph-subgraph isomorphisms. In the latter
     *  case graph 'a' is the subgraph, implying a.size() < b.size(). In the case that
     *  no isomorphism is found an empty mapping is returned.
     * 
     * 
     * @param shouldMatchBonds 
     * @param shouldMatchRings 
     * @return
     */
    private synchronized void isomorphism() {

        if (shouldMatchRings) {
            try {
                initializeMolecule(source);
                initializeMolecule(target);
            } catch (CDKException ex) {
                Logger.error(Level.SEVERE, null, ex);
            }
        }

        if (!isDead(source, target) && testIsSubgraphHeuristics(source, target, shouldMatchBonds)) {
            State state = new State(source, target, shouldMatchBonds, shouldMatchRings);
            if (!state.isDead()) {
                state.matchFirst(state, allAtomMCS);
            }
        }
    }

    /** 
     * The isomorphism method returns an isomorphism between two molecular
     *  graphs using the VF2Automorphism algorithm. This can be used for finding both
     *  graph-graph isomorphisms and graph-subgraph isomorphisms. In the latter
     *  case graph 'a' is the subgraph, implying a.size() < b.size(). In the case that
     *  no isomorphism is found an empty mapping is returned.
     * 
     * 
     */
    private synchronized void isomorphisms() {

        if (shouldMatchRings) {
            try {
                initializeMolecule(source);
                initializeMolecule(target);
            } catch (CDKException ex) {
                Logger.error(Level.SEVERE, null, ex);
            }
        }

        if (!isDead(source, target) && testIsSubgraphHeuristics(source, target, shouldMatchBonds)) {
            State state = new State(source, target, shouldMatchBonds, shouldMatchRings);
            if (!state.isDead()) {
                state.matchAll(state, allAtomMCS);
            }
        }
    }

    // Returns true substructure is bigger than the target
    private synchronized boolean isDead(IAtomContainer a, IAtomContainer b) {
        return a.getAtomCount() > b.getAtomCount();
    }

    @Override
    public boolean isSubgraph() {
        isomorphism();
        return allAtomMCS.isEmpty() ? false : true;
    }

    public boolean isSubgraphs() {
        isomorphisms();
        return allAtomMCS.isEmpty() ? false : true;
    }

    @Override
    public void set(IAtomContainer source, IAtomContainer target) throws CDKException {
        this.source = source;
        this.target = target;
    }

    @Override
    public void set(IQueryAtomContainer source, IAtomContainer target) throws CDKException {
        this.source = source;
        this.target = target;
    }

    @Override
    public List<AtomAtomMapping> getAllAtomMapping() {
        return allAtomMCS;
    }

    @Override
    public List<Map<Integer, Integer>> getAllMapping() {
        return allMCS;
    }

    @Override
    public AtomAtomMapping getFirstAtomMapping() {
        if (allAtomMCS.iterator().hasNext()) {
            return allAtomMCS.iterator().next();
        }
        return new AtomAtomMapping(source, target);
    }

    @Override
    public Map<Integer, Integer> getFirstMapping() {
        if (allMCS.iterator().hasNext()) {
            return allMCS.iterator().next();
        }
        return new TreeMap<Integer, Integer>();
    }
}
