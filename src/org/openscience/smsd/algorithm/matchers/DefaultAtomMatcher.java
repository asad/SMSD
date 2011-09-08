/* Copyright (C) 2011  Syed Asad Rahman <asad@ebi.ac.uk>
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
 */
package org.openscience.smsd.algorithm.matchers;

import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.annotations.TestClass;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.isomorphism.matchers.IQueryAtom;
import org.openscience.smsd.algorithm.matchers.ring.IRingMatcher;
import org.openscience.smsd.algorithm.matchers.ring.RingMatcher;


/**
 * Checks if atom is matching between query and target molecules.
 * @cdk.module smsd
 * @cdk.githash
 * @author Syed Asad Rahman <asad@ebi.ac.uk>
 */
@TestClass("org.openscience.cdk.smsd.algorithm.vflib.VFLibTest")
public final class DefaultAtomMatcher implements AtomMatcher {

    static final long serialVersionUID = -7861469841127327812L;
    private final String symbol;
    private final IAtom qAtom;
    private final IRingMatcher ringMatcher;
    private final boolean shouldMatchRings;


    /**
     * Constructor
     */
    public DefaultAtomMatcher() {
        this.qAtom = null;
        this.symbol = null;
        this.ringMatcher = null;
        this.shouldMatchRings = false;
    }

    /**
     * Constructor
     * @param atom query atom
     * @param shouldMatchRings ring matching flag
     */
    public DefaultAtomMatcher(IAtom atom, boolean shouldMatchRings) {
        super();
        this.qAtom = atom;
        this.symbol = atom.getSymbol();
        this.shouldMatchRings = shouldMatchRings;
        this.ringMatcher = new RingMatcher(atom);

//        System.out.println("Atom " + atom.getSymbol());
//        System.out.println("MAX allowed " + maximumNeighbors);
    }
    
    private boolean matchSymbol(IAtom atom) {
        if (symbol == null) {
            return false;
        }
        return symbol.equals(atom.getSymbol());
    }

    /** {@inheritDoc}
     *
     * @param targetAtom
     * @return
     */
    @Override
    public boolean matches(IAtom targetAtom) {
        if (targetAtom instanceof IQueryAtom) {
            return ((IQueryAtom) targetAtom).matches(qAtom);
        } else if (qAtom != null && qAtom instanceof IQueryAtom) {
            return ((IQueryAtom) qAtom).matches(targetAtom);
        } else {
            if (!matchSymbol(targetAtom)) {
                return false;
            }
            if (shouldMatchRings) {
                if (matchRingAtoms(targetAtom)) {
                    return ringMatcher.matches(targetAtom);
                } else {
                    return matchNonRingAtoms(targetAtom);
                }
            }
        }
        return true;
    }

    private boolean matchRingAtoms(IAtom tAtom) {
        return isRingAtom(qAtom) && isRingAtom(tAtom) ? true : false;
    }

    private boolean matchNonRingAtoms(IAtom tAtom) {
        return !isRingAtom(qAtom) && !isRingAtom(tAtom) ? true : false;
    }

    private boolean isRingAtom(IAtom atom) {
        return atom.getFlag(CDKConstants.ISINRING) ? true : false;
    }
}
