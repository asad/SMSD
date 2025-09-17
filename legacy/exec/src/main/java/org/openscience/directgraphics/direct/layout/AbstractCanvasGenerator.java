/*
 * Copyright (C) 2007-2018 Syed Asad Rahman <asad @ ebi.ac.uk>.
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301  USA
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
