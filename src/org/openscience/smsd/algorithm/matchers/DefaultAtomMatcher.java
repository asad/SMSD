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
import org.openscience.smsd.interfaces.IMoleculeInitializer;

/**
 * Checks if atom is matching between query and target molecules. @cdk.module smsd @cdk.githash
 *
 * @author Syed Asad Rahman <asad@ebi.ac.uk>
 */
@TestClass("org.openscience.cdk.smsd.algorithm.vflib.VFLibTest")
public final class DefaultAtomMatcher implements AtomMatcher {

    static final long serialVersionUID = -7861469841127327812L;
    private final String symbol;
    private final IAtom qAtom;
    private final boolean shouldMatchRings;
    private int smallestRingSize;

    /**
     * Constructor
     *
     * @param qAtom
     * @param symbol
     * @param shouldMatchRings
     */
    public DefaultAtomMatcher(IAtom qAtom,
            String symbol,
            boolean shouldMatchRings) {
        this.qAtom = qAtom;
        this.symbol = symbol;
        this.shouldMatchRings = shouldMatchRings;
        if (shouldMatchRings) {
            this.smallestRingSize = isRingAtom(qAtom) ? getRingSize(qAtom).intValue() : 0;
        }
    }

    /**
     * Constructor
     *
     * @param atom query atom
     * @param shouldMatchRings ring matching flag
     */
    public DefaultAtomMatcher(IAtom atom, boolean shouldMatchRings) {
        this(atom, atom.getSymbol(), shouldMatchRings);
    }

    private boolean matchSymbol(IAtom atom) {
        if (symbol == null) {
            return false;
        }
        return symbol.equals(atom.getSymbol());
    }

    /**
     * {@inheritDoc}
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
            if (shouldMatchRings && (isRingAtom(qAtom) && isAliphaticAtom(targetAtom))) {
                return false;
            }
            if (shouldMatchRings && (isAliphaticAtom(qAtom) && isRingAtom(targetAtom))) {
                return false;
            }
            if (shouldMatchRings
                    && (isRingAtom(qAtom)
                    && isRingAtom(targetAtom)
                    && !isRingSizeEqual(targetAtom))) {
                return false;
            }
            if (shouldMatchRings
                    && !(isRingAtom(qAtom) && isRingAtom(targetAtom))
                    && qAtom.getHybridization() != null
                    && targetAtom.getHybridization() != null
                    && !qAtom.getHybridization().equals(targetAtom.getHybridization())) {
                return false;
            }
        }
        return true;
    }

    private boolean isAliphaticAtom(IAtom atom) {
        return atom.getFlag(CDKConstants.ISALIPHATIC) ? true : false;
    }

    private boolean isRingAtom(IAtom atom) {
        return atom.getFlag(CDKConstants.ISINRING) ? true : false;
    }

    private boolean isRingSizeEqual(IAtom atom) {

        if (isRingAtom(atom) && isRingAtom(qAtom)) {

            Integer ringSize = getRingSize(atom);
            if (ringSize == null) {
                return false;
            }

            if (ringSize.intValue() == this.smallestRingSize) {
                return true;
            }
        }
        return false;

    }

    private Integer getRingSize(IAtom atom) {
        return (Integer) atom.getProperty(IMoleculeInitializer.SMALLEST_RING_SIZE);
    }
}
