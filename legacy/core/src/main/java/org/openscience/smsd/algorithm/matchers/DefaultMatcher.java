/*
 * SPDX-License-Identifier: Apache-2.0
 * © 2025 BioInception PVT LTD.
 */
package org.openscience.smsd.algorithm.matchers;

import org.openscience.cdk.interfaces.IBond;

/**
 * Checks if atom is matching between query and target molecules.
 *
 * 
 * 
 *
 * @author Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>
 */
public class DefaultMatcher {

    /**
     *
     * @param bondMatcher
     * @param bondA2
     * @param shouldMatchBonds
     * @return true if condition meet else false
     */
    private static boolean isBondMatch(BondMatcher queryBondMatcher, IBond targetBond) {
        return queryBondMatcher.matches(targetBond);
    }

    private static boolean isAtomMatch(IBond bondA1, IBond bondA2, boolean shouldMatchRings, boolean matchAtomTypes) {

        AtomMatcher atomMatcher1;
        AtomMatcher atomMatcher2;

        if (matchAtomTypes) {
            atomMatcher1 = new DefaultAtomTypeMatcher(bondA1.getAtom(0), shouldMatchRings);
            atomMatcher2 = new DefaultAtomTypeMatcher(bondA1.getAtom(1), shouldMatchRings);
        } else {
            atomMatcher1 = new DefaultAtomMatcher(bondA1.getAtom(0), shouldMatchRings);
            atomMatcher2 = new DefaultAtomMatcher(bondA1.getAtom(1), shouldMatchRings);
        }

        // ok, atoms match
        if (atomMatcher1.matches(bondA2.getAtom(0)) && atomMatcher2.matches(bondA2.getAtom(1))) {
//            System.out.println("Atom Matched");
            return true;
        }
        return atomMatcher1.matches(bondA2.getAtom(1)) && atomMatcher2.matches(bondA2.getAtom(0));
    }

    /**
     *
     * @param bondA1
     * @param bondA2
     * @param matchBond
     * @param shouldMatchRings
     * @param matchAtomTypes (atom type also matched and symbol matched)
     * @return true if condition meet else false
     */
    public static boolean matches(IBond bondA1, IBond bondA2,
            boolean matchBond, boolean shouldMatchRings, boolean matchAtomTypes) {

        if (!isAtomMatch(bondA1, bondA2, shouldMatchRings, matchAtomTypes)) {
            return false;
        }
        return isBondMatch(new DefaultBondMatcher(bondA1, matchBond, shouldMatchRings, matchAtomTypes), bondA2);
    }
}
