/*
 * SPDX-License-Identifier: Apache-2.0
 * © 2025 BioInception PVT LTD.
 */
package org.openscience.smsd.algorithm.vflib.vf2;

import java.util.List;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.isomorphism.matchers.IQueryAtom;

/**
 * Checks if atom is matching between query and target molecules.
 *
 *
 *
 *
 * @author Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>
 */
public final class DefaultAtomMatcher implements AtomMatcher {

    static final long serialVersionUID = -7861469841127327812L;
    private final boolean shouldMatchRings;
    private final boolean atomTypeMatcher;

    /**
     * Constructor
     *
     * @param shouldMatchRings match rings
     * @param atomTypeMatcher match atom types
     */
    public DefaultAtomMatcher(boolean shouldMatchRings, boolean atomTypeMatcher) {
        this.shouldMatchRings = shouldMatchRings;
        this.atomTypeMatcher = atomTypeMatcher;
    }

    private boolean matchSymbol(IAtom qAtom, IAtom atom) {
        if (qAtom == null) {
            return false;
        }
        return qAtom.getSymbol().equals(atom.getSymbol());
    }

    /**
     * {@inheritDoc}
     *
     * @param targetAtom
     * @return
     */
    @Override
    public boolean matches(IAtom queryAtom, IAtom targetAtom) {

        if (targetAtom instanceof IQueryAtom) {
            return ((IQueryAtom) targetAtom).matches(queryAtom);
        } else if (queryAtom != null && queryAtom instanceof IQueryAtom) {
            return ((IQueryAtom) queryAtom).matches(targetAtom);
        } else {
            if (!matchSymbol(queryAtom, targetAtom)) {
                return false;
            }

            if (isMatchRings()
                    && isRingAtom(queryAtom)
                    && !isRingAtom(targetAtom)) {
                return false;
            }

            if (isMatchRings()
                    && !isRingAtom(queryAtom)
                    && isRingAtom(targetAtom)) {
                return false;
            }

            if (isMatchRings()
                    && isRingAtom(queryAtom)
                    && isRingAtom(targetAtom)
                    && !isRingSizeMatch(queryAtom, targetAtom)) {
                return false;
            }

            if (isAtomTypeMatcher() && !matchAtomType(queryAtom, targetAtom)
                    && isAliphaticAtom(queryAtom)
                    && isAliphaticAtom(targetAtom)) {
                return false;
            }

        }
        return true;
    }

    private boolean isRingSizeMatch(IAtom qAtom, IAtom atom) {
        List<Integer> ringsizesQ = qAtom.getProperty(CDKConstants.RING_SIZES);
        List<Integer> ringsizesT = atom.getProperty(CDKConstants.RING_SIZES);
        if (ringsizesQ != null && ringsizesT != null) {
            if (ringsizesQ.stream().anyMatch((i) -> (ringsizesT.contains(i)))) {
                return true;
            }
        }
        return false;
    }

    private boolean isAliphaticAtom(IAtom atom) {
        return atom.getFlag(CDKConstants.ISALIPHATIC);
    }

    private boolean isRingAtom(IAtom atom) {
        return atom.getFlag(CDKConstants.ISINRING);
    }

    /**
     * @return the shouldMatchRings
     */
    public boolean isMatchRings() {
        return shouldMatchRings;
    }

    private boolean matchAtomType(IAtom qAtom, IAtom targetAtom) {
        String rAtom = qAtom.getAtomTypeName() == null
                ? qAtom.getSymbol() : qAtom.getAtomTypeName();
        String tAtom = targetAtom.getAtomTypeName() == null
                ? targetAtom.getSymbol() : targetAtom.getAtomTypeName();
        return rAtom.equals(tAtom);
    }

    /**
     * @return the atomTypeMatcher
     */
    public boolean isAtomTypeMatcher() {
        return atomTypeMatcher;
    }
}
