/*
 * SPDX-License-Identifier: Apache-2.0
 * © 2025 BioInception PVT LTD.
 */
package org.openscience.smsd.algorithm.matchers;

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
public final class DefaultAtomTypeMatcher implements AtomMatcher {

    static final long serialVersionUID = -7861469841127327812L;
    private final String symbol;
    private final IAtom qAtom;
    private final boolean shouldMatchRings;

    /**
     * Constructor
     *
     * @param qAtom
     * @param symbol
     * @param shouldMatchRings
     */
    public DefaultAtomTypeMatcher(IAtom qAtom,
            String symbol,
            boolean shouldMatchRings) {
        this.qAtom = qAtom;
        this.symbol = symbol;
        this.shouldMatchRings = shouldMatchRings;
    }

    /**
     * Constructor
     *
     * @param atom query atom
     * @param shouldMatchRings ring matching flag
     */
    public DefaultAtomTypeMatcher(IAtom atom, boolean shouldMatchRings) {
        this(atom, atom.getSymbol(), shouldMatchRings);
    }

    private boolean matchSymbol(IAtom atom) {
        if (getAtomSymbol() == null) {
            return false;
        }
        return getAtomSymbol().equals(atom.getSymbol());
    }

    /**
     * {@inheritDoc}
     *
     * @param targetAtom
     * @return true if condition meet else false
     */
    @Override
    public boolean matches(IAtom targetAtom) {

        if (targetAtom instanceof IQueryAtom) {
            return ((IQueryAtom) targetAtom).matches(getQueryAtom());
        } else if (getQueryAtom() != null && getQueryAtom() instanceof IQueryAtom) {
            return ((IQueryAtom) getQueryAtom()).matches(targetAtom);
        } else {
            if (!matchSymbol(targetAtom)) {
                return false;
            }

            if (isMatchRings()
                    && isRingAtom(getQueryAtom())
                    && !isRingAtom(targetAtom)) {
                return false;
            }

            if (isMatchRings()
                    && !isRingAtom(getQueryAtom())
                    && isRingAtom(targetAtom)) {
                return false;
            }

            if (isMatchRings()
                    && isRingAtom(getQueryAtom())
                    && isRingAtom(targetAtom)
                    && !isRingSizeMatch(targetAtom)) {
                return false;
            }

            if (!matchAtomType(targetAtom)
                    && isAliphaticAtom(getQueryAtom())
                    && isAliphaticAtom(targetAtom)) {
                return false;
            }

        }
        return true;
    }

    private boolean isRingSizeMatch(IAtom atom) {
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
     * @return the qAtom
     */
    @Override
    public IAtom getQueryAtom() {
        return qAtom;
    }

    /**
     * @return the symbol
     */
    public String getAtomSymbol() {
        return symbol;
    }

    /**
     * @return the shouldMatchRings
     */
    public boolean isMatchRings() {
        return shouldMatchRings;
    }

    private boolean matchAtomType(IAtom targetAtom) {
        String rAtom = qAtom.getAtomTypeName() == null
                ? qAtom.getSymbol() : qAtom.getAtomTypeName();
        String tAtom = targetAtom.getAtomTypeName() == null
                ? targetAtom.getSymbol() : targetAtom.getAtomTypeName();
        return rAtom.equals(tAtom);
    }
}
