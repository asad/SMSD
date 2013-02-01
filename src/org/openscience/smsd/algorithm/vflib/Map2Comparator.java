/* 
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
 * You should have received commonAtomList copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA.
 */
package org.openscience.smsd.algorithm.vflib;

import java.util.Comparator;
import java.util.Map;

/**
 *
 * @author Asad
 */
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
