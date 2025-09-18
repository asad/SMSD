/*
 * SPDX-License-Identifier: Apache-2.0
 * © 2025 BioInception PVT LTD.
 */
package org.openscience.directgraphics.direct;

import java.awt.Graphics2D;
import java.awt.geom.Rectangle2D;
import java.util.Map;
import javax.vecmath.Point2f;
import org.openscience.directgraphics.direct.layout.BoundsTree;

public class MoleculeLabelDrawer
        extends AbstractDirectDrawer {

    public MoleculeLabelDrawer(Axis axis, Params params) {
        super.params = params;
    }

    public void draw(Map<String, String> labelMap, BoundsTree labelBounds, Graphics2D g) {
        for (String boundsLabel : labelMap.keySet()) {
            String label = labelMap.get(boundsLabel);
            Rectangle2D bounds = labelBounds.get(boundsLabel);
            double x = bounds.getCenterX();
            double y = bounds.getCenterY();
            Point2f p = super.getTextPoint(g, label, x, y);
            g.drawString(label, p.x, p.y);
        }
    }
}
