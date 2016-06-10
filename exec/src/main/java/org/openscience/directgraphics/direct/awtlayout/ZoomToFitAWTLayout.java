/*
 * Copyright (C) 2007-2015 Syed Asad Rahman <asad @ ebi.ac.uk>.
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

import java.awt.Color;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.Graphics2D;
import java.awt.geom.AffineTransform;
import java.awt.geom.Rectangle2D;
import javax.vecmath.Vector2d;
import org.openscience.cdk.geometry.GeometryTools;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.directgraphics.direct.DirectMoleculeDrawer;
import org.openscience.directgraphics.direct.Params;
import org.openscience.directgraphics.direct.layout.AbstractDirectLayout;
import org.openscience.directgraphics.direct.layout.BoundsTree;

public class ZoomToFitAWTLayout
        extends AbstractDirectLayout<IAtomContainer> {

    private final DirectMoleculeDrawer drawer;

    public ZoomToFitAWTLayout(DirectMoleculeDrawer drawer) {
        this.drawer = drawer;
        super.setParams(drawer.getParams());
    }

    public void layout(IAtomContainer mol, Dimension cellCanvas, Graphics2D g) {
        double centerY;
        AffineTransform originalTransform = g.getTransform();
        double cW = cellCanvas.width;
        double cH = cellCanvas.height;
        if (this.shouldInvert) {
            super.invert(mol);
        }
        BoundsTree tree = this.getBoundsTree(mol, g);
        double tW = tree.getWidth();
        double tH = tree.getHeight();
        Rectangle2D stringBounds = null;
        String label = mol.getID();
        Font labelFont = new Font(this.params.labelPanelFont, 0, this.params.labelPanelFontSize);
        if (this.params.drawLabelPanel) {
            g.setFont(labelFont);
            FontMetrics metrics = g.getFontMetrics();
            stringBounds = metrics.getStringBounds(label, g);
            double labelHeight = stringBounds.getHeight();
            cH += labelHeight;
        }
        double zoom = this.calculateZoom(tW, tH, cW, cH);
        double centerX = cW / 2.0;
        Params paramsLocal = this.drawer.getParams();
        if (paramsLocal.drawMoleculeID) {
            centerY = cH / 2.0 - paramsLocal.labelYGap;
        } else if (paramsLocal.drawLabelPanel) {
            double labelHeight = stringBounds.getHeight();
            double scaledLabelHeight = labelHeight / 2.0 * (1.0 / zoom);
            centerY = (double) cellCanvas.height / 2.0 - scaledLabelHeight;
        } else {
            centerY = cH / 2.0;
        }
        g.translate(centerX, centerY);
        g.scale(zoom, zoom);
        this.drawer.drawMolecule(mol, g);
        g.setTransform(originalTransform);
        if (paramsLocal.drawLabelPanel) {
            double cX = cW / 2.0;
            double cY = cH / 2.0;
            g.setFont(labelFont);
            FontMetrics metrics = g.getFontMetrics();
            double halfWidth = stringBounds.getWidth() / 2.0;
            double halfHeight = stringBounds.getHeight() / 2.0;
            double halfScaledTreeWidth = tH * zoom / 2.0;
            double lY = cY + halfScaledTreeWidth - (double) paramsLocal.borderY;
            double ascent = metrics.getAscent();
            float x = (float) (cX - halfWidth);
            float y = (float) (lY - halfHeight + ascent);
            g.setColor(Color.BLACK);
            g.drawString(label, x, y);
        }
    }

    private BoundsTree getBoundsTree(IAtomContainer mol, Graphics2D g) {
        Rectangle2D bb = GeometryTools.getRectangle2D((IAtomContainer) mol);
        GeometryTools.translate2D((IAtomContainer) mol, (double) (-bb.getCenterX()), (double) (-bb.getCenterY()));
        GeometryTools.scaleMolecule((IAtomContainer) mol, (double) GeometryTools.getScaleFactor((IAtomContainer) mol, (double) this.params.bondLength));
        MoleculeLayout exactLayout = new MoleculeLayout(this.params);
        return exactLayout.layout(mol, g);
    }

    private double calculateZoom(double tw, double th, double cw, double ch) {
        Params paramsLocal = this.drawer.getParams();
        double borderX = paramsLocal.borderX;
        double borderY = paramsLocal.borderY;
        double rW = tw + borderX * 2.0;
        double rH = th + borderY * 2.0;
        return Math.min(cw / rW, ch / rH);
    }

    @Override
    public BoundsTree layout(IAtomContainer obj, Vector2d axis) {
        return null;
    }

    @Override
    public Vector2d getAxis() {
        return null;
    }

    @Override
    public double getAxisPosition() {
        return 0.0;
    }
}
