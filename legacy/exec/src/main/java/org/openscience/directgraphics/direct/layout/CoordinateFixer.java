/*
 * SPDX-License-Identifier: Apache-2.0
 * © 2025 BioInception PVT LTD.
 */
package org.openscience.directgraphics.direct.layout;

import javax.vecmath.Point2d;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IMapping;
import org.openscience.cdk.interfaces.IReaction;

public class CoordinateFixer {
    public static void fixCoordinates(IReaction reaction) {
        for (IMapping mapping : reaction.mappings()) {
            IAtom a0 = (IAtom)mapping.getChemObject(0);
            IAtom a1 = (IAtom)mapping.getChemObject(1);
            if (a0 == null || a1 == null) continue;
            a1.setPoint2d(new Point2d(a0.getPoint2d()));
        }
    }
}

