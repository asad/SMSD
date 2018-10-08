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
package org.openscience.directgraphics.direct;

import java.awt.Color;
import java.awt.FontMetrics;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.RenderingHints;
import java.awt.geom.Line2D;
import java.awt.geom.Rectangle2D;
import java.awt.image.BufferedImage;
import javax.vecmath.Point2d;
import javax.vecmath.Point2f;

public class AbstractDirectDrawer {

    protected Params params;

    public Params getParams() {
        return this.params;
    }

    public void setParams(Params params) {
        this.params = params;
    }

    public void drawLine(Point2d p1, Point2d p2, Graphics2D g) {
        g.draw(new Line2D.Double(p1.x, p1.y, p2.x, p2.y));
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

    public BufferedImage makeBlankImage(int w, int h) {
        return this.makeBlankImage(w, h, Color.WHITE);
    }

    public BufferedImage makeBlankImage(int w, int h, Color color) {
        BufferedImage image = new BufferedImage(w, h, 2);
        Graphics2D g = (Graphics2D) image.getGraphics();
        g.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
        g.setColor(color);
        g.fillRect(0, 0, w, h);
        return image;
    }
}
