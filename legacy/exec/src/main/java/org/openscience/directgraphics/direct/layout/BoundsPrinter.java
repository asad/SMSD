/*
 * SPDX-License-Identifier: Apache-2.0
 * © 2025 BioInception PVT LTD.
 */
package org.openscience.directgraphics.direct.layout;

import java.awt.geom.Rectangle2D;

public class BoundsPrinter {

    public static String toString(Rectangle2D b) {
        return String.format("[(%2.0f, %2.0f), (%2.0f, %2.0f)] = (%2.0f x %2.0f) @ [%2.0f, %2.0f]", b.getMinX(), b.getMinY(), b.getMaxX(), b.getMaxY(), b.getWidth(), b.getHeight(), b.getCenterX(), b.getCenterY());
    }
}
