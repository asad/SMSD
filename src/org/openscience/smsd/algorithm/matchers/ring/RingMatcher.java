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
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA.
 */
package org.openscience.smsd.algorithm.matchers.ring;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.interfaces.IAtom;

/**
 * Checks if atom is in the ring, they are of equal size
 * then return true else returns false.
 * @cdk.module smsd
 * @cdk.githash
 * @author Syed Asad Rahman <asad@ebi.ac.uk>
 */
public class RingMatcher extends DefaultRingMatcher implements IRingMatcher {

    private static final long serialVersionUID = -7867686677884205622L;
    private final List<Integer> queryRingSizes;

    /**
     * Creates a new instance
     *
     * @param atom
     */
    @SuppressWarnings("unchecked")
    public RingMatcher(IAtom atom) {
        if (isRingAtom(atom) && hasRingSize(atom)) {
            this.queryRingSizes = (List<Integer>) atom.getProperty(CDKConstants.RING_SIZES);
        } else {
            this.queryRingSizes = new ArrayList<Integer>(0);
        }
    }

    /**
     * Returns true if query and target are atoms of same size.
     * @param targetAtom target atom of be checked 
     * @return true if query and target are atoms of same size
     */
    @Override
    public boolean matches(IAtom targetAtom) {
        if (queryRingSizes.isEmpty()) {
            return false;
        }
        if (isRingAtom(targetAtom) && hasRingSize(targetAtom)) {
            @SuppressWarnings("unchecked")
            List<Integer> targetRingSizes = (List<Integer>) targetAtom.getProperty(CDKConstants.RING_SIZES);
            if (targetRingSizes.isEmpty()) {
                return false;
            }
            return true;
        }
        return false;
    }

    private boolean isRingAtom(IAtom atom) {
        return atom.getFlag(CDKConstants.ISINRING) ? true : false;
    }

    private boolean hasRingSize(IAtom atom) {
        return atom.getProperty(CDKConstants.RING_SIZES) == null ? false : true;
    }
}
