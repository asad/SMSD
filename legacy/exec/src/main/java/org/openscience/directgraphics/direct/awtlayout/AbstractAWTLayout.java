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
package org.openscience.directgraphics.direct.awtlayout;

import java.awt.FontMetrics;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.geom.Rectangle2D;
import javax.vecmath.Point2f;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.directgraphics.direct.Params;
import org.openscience.directgraphics.direct.layout.BoundsTree;

public abstract class AbstractAWTLayout<T> {

    protected Graphics2D graphics;
    protected AbstractAWTLayout parent;
    protected T currentObject;
    protected Params params;
    protected BoundsTree boundsTree;

    public Params getParams() {
        return this.params;
    }

    public void setParams(Params params) {
        this.params = params;
    }

    public abstract BoundsTree layout(T var1, Graphics2D var2);

    public abstract BoundsTree layout(T var1, String var2, Graphics2D var3);

    public BoundsTree getBoundsTree() {
        return this.boundsTree;
    }

    public T getCurrentObject() {
        return this.currentObject;
    }

    public Graphics2D getGraphics() {
        return this.graphics;
    }

    public void setGraphics(Graphics2D graphics) {
        this.graphics = graphics;
    }

    public Point2f getTextPoint(Graphics g, String text, double cX, double cY) {
        FontMetrics metrics = g.getFontMetrics();
        Rectangle2D stringBounds = metrics.getStringBounds(text, g);
        double halfWidth = stringBounds.getWidth() / 2.0;
        double halfHeight = stringBounds.getHeight() / 2.0;
        double ascent = metrics.getAscent();
        float x = (float) (cX - halfWidth);
        float y = (float) (cY - halfHeight + ascent);
        return new Point2f(x, y);
    }

    public Rectangle2D getTextBounds(Graphics g, String text) {
        FontMetrics metrics = g.getFontMetrics();
        return metrics.getStringBounds(text, g);
    }

    public void translateTo(IAtomContainer ac, double x, double y, BoundsTree boundsTree) {
        Rectangle2D bounds = boundsTree.getRoot();
        double dx = x - bounds.getCenterX();
        double dy = y - bounds.getCenterY();
        for (IAtom atom : ac.atoms()) {
            atom.getPoint2d().x += dx;
            atom.getPoint2d().y += dy;
        }
        boundsTree.shift(dx, dy);
    }
}
