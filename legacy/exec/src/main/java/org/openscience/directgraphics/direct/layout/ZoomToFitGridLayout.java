/*
 * SPDX-License-Identifier: Apache-2.0
 * © 2025 BioInception PVT LTD.
 */
package org.openscience.directgraphics.direct.layout;

import java.awt.Dimension;
import java.awt.Graphics2D;
import java.awt.geom.AffineTransform;
import java.awt.geom.Rectangle2D;
import java.util.List;
import javax.vecmath.Point2d;
import org.openscience.cdk.geometry.GeometryTools;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.directgraphics.direct.DirectMoleculeDrawer;

public class ZoomToFitGridLayout {

    private DirectMoleculeDrawer drawer;
    private int rows;
    private int cols;

    public ZoomToFitGridLayout() {
        this.drawer = new DirectMoleculeDrawer();
    }

    public ZoomToFitGridLayout(int rows, int cols) {
        this(new DirectMoleculeDrawer(), rows, cols);
    }

    public ZoomToFitGridLayout(DirectMoleculeDrawer drawer, int rows, int cols) {
        this.drawer = drawer;
        this.rows = rows;
        this.cols = cols;
    }

    public void layout(List<IAtomContainer> mols, Dimension cellCanvas, Graphics2D g) {
        AffineTransform originalTransform = g.getTransform();
        double w = cellCanvas.width;
        double h = cellCanvas.height;
        double centerX = w / 2.0;
        double centerY = h / 2.0;
        int colCounter = 1;
        for (IAtomContainer mol : mols) {
            double zoom = this.calculateZoom(mol, cellCanvas);
            g.translate(centerX, centerY);
            g.scale(zoom, zoom);
            this.drawer.drawMolecule(mol, g);
            g.setTransform(originalTransform);
            if (colCounter < this.cols) {
                centerX += w;
                ++colCounter;
                continue;
            }
            centerY += h;
            centerX = w / 2.0;
            colCounter = 1;
        }
    }

    private double calculateZoom(IAtomContainer ac, Dimension canvas) {
        double scaleFactor = GeometryTools.getScaleFactor((IAtomContainer) ac, (double) this.drawer.getParams().bondLength);
        GeometryTools.translate2DCenterTo((IAtomContainer) ac, (Point2d) new Point2d(0.0, 0.0));
        GeometryTools.scaleMolecule((IAtomContainer) ac, (double) scaleFactor);
        Rectangle2D r2D = GeometryTools.getRectangle2D((IAtomContainer) ac);
        double canvasWidth = canvas.width;
        double canvasHeight = canvas.height;
        double borderX = this.drawer.getParams().borderX;
        double borderY = this.drawer.getParams().borderY;
        double objectWidth = r2D.getWidth() + borderX * 2.0;
        double objectHeight = r2D.getHeight() + borderY * 2.0;
        return Math.min(canvasWidth / objectWidth, canvasHeight / objectHeight);
    }
}
