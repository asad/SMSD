/*
 * SPDX-License-Identifier: Apache-2.0
 * © 2025 BioInception PVT LTD.
 */
package org.openscience.directgraphics.direct.layout;

import java.awt.Dimension;
import java.awt.geom.Rectangle2D;
import java.util.HashMap;
import java.util.Map;
import javax.vecmath.Point2d;
import org.openscience.cdk.interfaces.IAtomContainer;

public abstract class AbstractCanvasGenerator
        implements CanvasGenerator {

    protected Map<IAtomContainer, Rectangle2D> canvasMap = new HashMap<>();

    @Override
    public Rectangle2D getCanvasForAtomContainer(IAtomContainer atomContainer) {
        return this.canvasMap.get(atomContainer);
    }

    public void createCanvas(IAtomContainer atomContainer, Point2d center, Dimension canvasDimensions) {
        double w = canvasDimensions.width;
        double h = canvasDimensions.height;
        double x = center.x - w / 2.0;
        double y = center.y - h / 2.0;
        Rectangle2D.Double canvas = new Rectangle2D.Double(x, y, w, h);
        this.canvasMap.put(atomContainer, canvas);
    }
}
