/*
 *
 *
 * Copyright (C) 2011  Syed Asad Rahman <asad@ebi.ac.uk>
 *                     Gilleain Torrance <gilleain.torrance@gmail.com>
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

import java.util.Collections;
import java.util.HashMap;
import java.util.Map;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;

/**
 * Holds atom-atom mappings between source and target atoms
 * 
 * @cdk.module smsd
 * @cdk.githash
 * @author Syed Asad Rahman <asad@ebi.ac.uk>
 */
public final class AtomAtomMapping {

    private final IAtomContainer query;
    private final IAtomContainer target;
    private final Map<IAtom, IAtom> mapping;

    /**
     * 
     * @param query source molecule
     * @param target target molecule
     */
    public AtomAtomMapping(IAtomContainer query, IAtomContainer target) {
        this.query = query;
        this.target = target;
        this.mapping = Collections.synchronizedMap(new HashMap<IAtom, IAtom>());
    }

    /**
     * 
     * @param atom1
     * @param atom2
     */
    public synchronized void put(IAtom atom1, IAtom atom2) {
        mapping.put(atom1, atom2);
    }

    @Override
    public synchronized String toString() {
        String s = "[";
        for (IAtom key : mapping.keySet()) {
            int keyIndex = getQuery().getAtomNumber(key);
            int valueIndex = getTarget().getAtomNumber(mapping.get(key));
            s += keyIndex + ":" + valueIndex + "|";
        }
        return s + "]";
    }

    /**
     * 
     * @return true if 'query' is not query subgraph of 'target'
     */
    public synchronized boolean isEmpty() {
        return mapping.isEmpty();
    }

    /**
     * 
     * clear mapping
     */
    public synchronized void clear() {
        mapping.clear();
    }

    /**
     * 
     * get mapping size 
     * @return
     */
    public synchronized int getCount() {
        return mapping.isEmpty() ? 0 : mapping.size();
    }

    /**
     * get atom-atom mappings
     * @return
     */
    public synchronized Map<IAtom, IAtom> getMappings() {
        return Collections.unmodifiableMap(new HashMap<IAtom, IAtom>(mapping));
    }

    /**
     * Returns atom index of the given atom
     * in the query molecule
     * @param atom
     * @return
     */
    public synchronized int getQueryIndex(IAtom atom) {
        return getQuery().getAtomNumber(atom);
    }

    /**
     * Returns atom index of the given atom
     * in the target molecule
     * @param atom
     * @return
     */
    public synchronized int getTargetIndex(IAtom atom) {
        return getTarget().getAtomNumber(atom);
    }

    /**
     * @return the query
     */
    public synchronized IAtomContainer getQuery() {
        return query;
    }

    /**
     * @return the target
     */
    public synchronized IAtomContainer getTarget() {
        return target;
    }
}