/*
 *
 *
 * Copyright (C) 2009-2011  Syed Asad Rahman <asad@ebi.ac.uk>
 *                          Gilleain Torrance <gilleain.torrance@gmail.com>
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
public class AtomMapping {

    private IAtomContainer a;
    private IAtomContainer b;
    private Map<IAtom, IAtom> mapping;

    /**
     * 
     * @param a source mol
     * @param b target mol
     */
    public AtomMapping(IAtomContainer a, IAtomContainer b) {
        this.a = a;
        this.b = b;
        this.mapping = new HashMap<IAtom, IAtom>();
    }

    /**
     * 
     * @param atom1
     * @param atom2
     */
    public void add(IAtom atom1, IAtom atom2) {
        mapping.put(atom1, atom2);
    }

    @Override
    public String toString() {
        String s = "[";
        for (IAtom key : mapping.keySet()) {
            int keyIndex = a.getAtomNumber(key);
            int valueIndex = b.getAtomNumber(mapping.get(key));
            s += keyIndex + ":" + valueIndex + "|";
        }
        return s + "]";
    }

    /**
     * 
     * @return true if 'a' is not a subgraph of 'b'
     */
    public boolean isEmpty() {
        return mapping.isEmpty();
    }

    /**
     * 
     * clear mapping
     */
    public void clear() {
        mapping.clear();
    }

    /**
     * 
     * get mapping size 
     * @return
     */
    public int getSize() {
        return mapping.isEmpty() ? 0 : mapping.size();
    }
    
    /**
     * get atom-atom mappings
     * @return
     */
    public Map<IAtom, IAtom> getAtomMapping() {
        return mapping;
    }
}