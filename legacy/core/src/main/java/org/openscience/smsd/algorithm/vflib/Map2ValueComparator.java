/*
 * SPDX-License-Identifier: Apache-2.0
 * © 2025 BioInception PVT LTD.
 */
package org.openscience.smsd.algorithm.vflib;

import java.util.Comparator;
import java.util.Map;
import org.openscience.cdk.interfaces.IAtom;

/*
 * 
 * 
 *
 * @author Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>
 */
public class Map2ValueComparator implements Comparator<Map<IAtom, IAtom>> {

    private final SortOrder sortOrder;

    public Map2ValueComparator(SortOrder sortOrder) {
        this.sortOrder = sortOrder;
    }

    /**
     *
     * @param object1
     * @param object2
     * @return
     */
    @Override
    public int compare(Map<IAtom, IAtom> object1, Map<IAtom, IAtom> object2) {
        int size1 = object1.size();
        int size2 = object2.size();
        int compare = Integer.signum(new Integer(size1).compareTo(size2));

        if (sortOrder == SortOrder.ASCENDING) {
            return compare;
        } else {
            return compare * (-1);
        }
        //return size2 - size1;  assumes you want biggest to smallest;
    }
}
