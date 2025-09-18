/*
 * SPDX-License-Identifier: Apache-2.0
 * © 2025 BioInception PVT LTD.
 */
package org.openscience.smsd.algorithm.vflib.vf2;

import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.isomorphism.matchers.IQueryBond;

/**
 * Checks if a bond is matching between query and target molecules.
 *
 * 
 * 
 *
 * @author Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>
 */
public final class DefaultBondMatcher implements BondMatcher {

    static final long serialVersionUID = -7861469841127328812L;
    private final boolean shouldMatchBonds;
    private final boolean matchAtomTypes;
    private final boolean shouldMatchRings;

    /**
     * Constructor
     */
    public DefaultBondMatcher() {
        shouldMatchBonds = false;
        matchAtomTypes = false;
        shouldMatchRings = false;
    }

    /**
     * Constructor
     *
     * @param shouldMatchBonds bond match flag
     * @param shouldMatchRings bond in rings
     * @param matchAtomTypes atom types matched
     */
    public DefaultBondMatcher(
            boolean shouldMatchBonds,
            boolean shouldMatchRings,
            boolean matchAtomTypes) {
        this.shouldMatchBonds = shouldMatchBonds;
        this.shouldMatchRings = shouldMatchRings;
        this.matchAtomTypes = matchAtomTypes;
    }

    /**
     * {@inheritDoc}
     *
     * @param queryBond query GraphMolecule
     * @param targetBond target bond
     * @return true if bonds match
     */
    @Override
    public boolean matches(IBond queryBond, IBond targetBond) {

        if (queryBond != null && queryBond instanceof IQueryBond) {
            return ((IQueryBond) queryBond).matches(targetBond);
        } else if ((queryBond != null && targetBond != null)
                && isBondMatchFlag() && isBondTypeMatch(queryBond, targetBond)) {
            return true;
        } else if ((queryBond != null && targetBond != null)
                && !isBondMatchFlag() && isShouldMatchRings()) {
            if (queryBond.getFlag(CDKConstants.ISAROMATIC)
                    && targetBond.getFlag(CDKConstants.ISAROMATIC)) {
                return true;
            } else if (!queryBond.getFlag(CDKConstants.ISAROMATIC)
                    && !targetBond.getFlag(CDKConstants.ISAROMATIC)) {
                return true;
            }
        } else if ((queryBond != null && targetBond != null)
                && !isBondMatchFlag() && !isShouldMatchRings()) {
            return true;
        }
        return false;
    }

    /**
     * Return true if a bond is matched between query and target
     *
     * @param targetBond
     * @return
     */
    private boolean isBondTypeMatch(IBond queryBond, IBond targetBond) {

        if ((queryBond.getFlag(CDKConstants.ISAROMATIC) == targetBond.getFlag(CDKConstants.ISAROMATIC))
                && queryBond.getOrder().equals(targetBond.getOrder())) {
            return true;
        }

        if (queryBond.getFlag(CDKConstants.ISAROMATIC) && targetBond.getFlag(CDKConstants.ISAROMATIC)) {
            return true;
        }

        return !matchAtomTypes
                && queryBond.getFlag(CDKConstants.ISINRING)
                && targetBond.getFlag(CDKConstants.ISINRING)
                && (queryBond.getOrder() == IBond.Order.UNSET
                || targetBond.getOrder() == IBond.Order.UNSET);
    }

    /**
     * @return the shouldMatchBonds
     */
    public boolean isBondMatchFlag() {
        return shouldMatchBonds;
    }

    /**
     * @return the shouldMatchRings
     */
    public boolean isShouldMatchRings() {
        return shouldMatchRings;
    }
}
