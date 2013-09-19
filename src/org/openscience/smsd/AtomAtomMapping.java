/*
 *
 *
 * Copyright (C) 2009-2013  Syed Asad Rahman <asad@ebi.ac.uk>
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
 * You should have received query copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA.
 * 
 * 
 */
package org.openscience.smsd;

import java.io.Serializable;
import java.util.*;
import org.openscience.cdk.graph.ConnectivityChecker;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;

/**
 * Holds atom-atom mappings information between source and target molecules
 *
 * @cdk.module smsd @cdk.githash
 *
 * @author Syed Asad Rahman <asad@ebi.ac.uk>
 */
public final class AtomAtomMapping implements Serializable {

    private static final long serialVersionUID = 1223637237262778L;
    private final IAtomContainer query;
    private final IAtomContainer target;
    private final Map<IAtom, IAtom> mapping;
    private final Map<Integer, Integer> mappingIndex;

    @Override
    public boolean equals(Object obj) {
        if (obj == null) {
            return false;
        }
        if (getClass() != obj.getClass()) {
            return false;
        }
        final AtomAtomMapping other = (AtomAtomMapping) obj;
        if (this.getQuery() != other.getQuery() && (this.getQuery() == null || !this.query.equals(other.query))) {
            return false;
        }
        if (this.getTarget() != other.getTarget() && (this.getTarget() == null || !this.target.equals(other.target))) {
            return false;
        }
        if (this.getMapping() != other.getMapping() && (this.getMapping() == null || !this.mapping.equals(other.mapping))) {
            return false;
        }
        if (this.getMappingByIndex() != other.getMappingByIndex() && (this.getMappingByIndex() == null || !this.mappingIndex.equals(other.mappingIndex))) {
            return false;
        }
        return true;
    }

    @Override
    public int hashCode() {
        int hash = 7;
        hash = 67 * hash + (this.getQuery() != null ? this.getQuery().hashCode() : 0);
        hash = 67 * hash + (this.getTarget() != null ? this.getTarget().hashCode() : 0);
        hash = 67 * hash + (this.getMapping() != null ? this.getMapping().hashCode() : 0);
        hash = 67 * hash + (this.getMappingByIndex() != null ? this.getMappingByIndex().hashCode() : 0);
        return hash;
    }

    /**
     *
     * @param query source molecule
     * @param target target molecule
     */
    public AtomAtomMapping(IAtomContainer query, IAtomContainer target) {
        this.query = query;
        this.target = target;
        this.mapping = Collections.synchronizedMap(new HashMap<IAtom, IAtom>());
        this.mappingIndex = Collections.synchronizedSortedMap(new TreeMap<Integer, Integer>());
    }

    /**
     *
     * @param atom1
     * @param atom2
     */
    public synchronized void put(IAtom atom1, IAtom atom2) {
        getMapping().put(atom1, atom2);
        getMappingByIndex().put(getQuery().getAtomNumber(atom1), getTarget().getAtomNumber(atom2));
    }

    /**
     * Returns String.
     *
     * @return string
     */
    @Override
    public synchronized String toString() {
        String s = "[";
        for (IAtom key : getMapping().keySet()) {
            int keyIndex = getQuery().getAtomNumber(key);
            int valueIndex = getTarget().getAtomNumber(getMapping().get(key));
            s += keyIndex + ":" + valueIndex + "|";
        }
        return s + "]";
    }

    /**
     * Returns true if 'query' is not isomorphic of 'target'.
     *
     * @return true if 'query' is not isomorphic of 'target'
     */
    public synchronized boolean isEmpty() {
        return getMapping().isEmpty();
    }

    /**
     *
     * Clear mappings
     */
    public synchronized void clear() {
        getMapping().clear();
        getMappingByIndex().clear();
    }

    /**
     *
     * Returns mapping size.
     *
     * @return mapping size
     */
    public synchronized int getCount() {
        return getMapping().isEmpty() ? 0 : getMapping().size();
    }

    /**
     * Returns atom-atom mappings
     *
     * @return atom-atom mappings
     */
    public synchronized Map<IAtom, IAtom> getMappingsByAtoms() {
        return Collections.unmodifiableMap(new HashMap<>(getMapping()));
    }

    /**
     * Returns atom-atom index mappings
     *
     * @return atom-atom index mappings
     */
    public synchronized Map<Integer, Integer> getMappingsByIndex() {
        return Collections.unmodifiableSortedMap(new TreeMap<>(getMappingByIndex()));
    }

    /**
     * Returns atom index of the given atom in the query molecule
     *
     * @param atom
     * @return
     */
    public synchronized int getQueryIndex(IAtom atom) {
        return getQuery().getAtomNumber(atom);
    }

    /**
     * Returns atom index of the given atom in the target molecule
     *
     * @param atom
     * @return
     */
    public synchronized int getTargetIndex(IAtom atom) {
        return getTarget().getAtomNumber(atom);
    }

    /**
     * Returns query molecule
     *
     * @return the query
     */
    public synchronized IAtomContainer getQuery() {
        return query;
    }

    /**
     * Returns target molecule
     *
     * @return the target
     */
    public synchronized IAtomContainer getTarget() {
        return target;
    }

    /**
     * Returns common mapped fragment in the query molecule.
     *
     * @return common mapped fragment in the query molecule
     * @throws CloneNotSupportedException
     */
    public synchronized IAtomContainer getCommonFragmentInQuery() throws CloneNotSupportedException {
        IAtomContainer ac = getQuery().clone();
        List<IAtom> uniqueAtoms = Collections.synchronizedList(new ArrayList<IAtom>());
        for (IAtom atom : getQuery().atoms()) {
            if (!mapping.containsKey(atom)) {
                uniqueAtoms.add(ac.getAtom(getQueryIndex(atom)));
            }
        }
        for (IAtom atom : uniqueAtoms) {
            ac.removeAtomAndConnectedElectronContainers(atom);
        }
        return ac;
    }

    /**
     * Returns common mapped fragment in the target molecule.
     *
     * @return common mapped fragment in the target molecule
     * @throws CloneNotSupportedException
     */
    public synchronized IAtomContainer getCommonFragmentInTarget() throws CloneNotSupportedException {
        IAtomContainer ac = getTarget().clone();
        List<IAtom> uniqueAtoms = Collections.synchronizedList(new ArrayList<IAtom>());
        for (IAtom atom : getTarget().atoms()) {
            if (!mapping.containsValue(atom)) {
                uniqueAtoms.add(ac.getAtom(getTargetIndex(atom)));
            }
        }
        for (IAtom atom : uniqueAtoms) {
            ac.removeAtomAndConnectedElectronContainers(atom);
        }
        return ac;
    }

    /**
     * Returns unique unmapped fragments in the query molecule.
     *
     * @return unique fragments in the query molecule
     * @throws CloneNotSupportedException
     */
    public synchronized IAtomContainerSet getUniqueFragmentsInQuery() throws CloneNotSupportedException {
        IAtomContainer ac = getQuery().clone();
        List<IAtom> commonAtoms = Collections.synchronizedList(new ArrayList<IAtom>());
        for (IAtom atom : getMapping().keySet()) {
            commonAtoms.add(ac.getAtom(getQueryIndex(atom)));
        }
        for (IAtom atom : commonAtoms) {
            ac.removeAtomAndConnectedElectronContainers(atom);
        }
        // now we probably have a set of disconnected components
        // so lets get a set of individual atom containers for
        // corresponding to each component
        return ConnectivityChecker.partitionIntoMolecules(ac);
    }

    /**
     * Returns unique unmapped fragments in the target molecule.
     *
     * @return unique fragments in the target molecule
     * @throws CloneNotSupportedException
     */
    public synchronized IAtomContainerSet getUniqueFragmentsInTarget() throws CloneNotSupportedException {
        IAtomContainer ac = getTarget().clone();
        List<IAtom> commonAtoms = Collections.synchronizedList(new ArrayList<IAtom>());
        for (IAtom atom : getMapping().values()) {
            commonAtoms.add(ac.getAtom(getTargetIndex(atom)));
        }
        for (IAtom atom : commonAtoms) {
            ac.removeAtomAndConnectedElectronContainers(atom);
        }
        // now we probably have a set of disconnected components
        // so lets get a set of individual atom containers for
        // corresponding to each component
        return ConnectivityChecker.partitionIntoMolecules(ac);
    }

    /**
     * @return the mapping
     */
    private Map<IAtom, IAtom> getMapping() {
        return mapping;
    }

    /**
     * @return the mappingIndex
     */
    private Map<Integer, Integer> getMappingByIndex() {
        return mappingIndex;
    }
}
