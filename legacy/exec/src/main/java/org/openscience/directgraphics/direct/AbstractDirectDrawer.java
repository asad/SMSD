/*
 * SPDX-License-Identifier: Apache-2.0
 * © 2025 BioInception PVT LTD.
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
