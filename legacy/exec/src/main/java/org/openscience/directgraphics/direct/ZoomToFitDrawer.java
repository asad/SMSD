/*
 * SPDX-License-Identifier: Apache-2.0
 * © 2025 BioInception PVT LTD.
 */
package org.openscience.directgraphics.direct;

import java.awt.Dimension;
import java.awt.Graphics2D;
import java.awt.geom.AffineTransform;
import java.awt.geom.Rectangle2D;
import java.util.List;
import javax.vecmath.Point2d;
import org.openscience.cdk.geometry.GeometryTools;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.directgraphics.direct.layout.CanvasGenerator;
import org.openscience.directgraphics.direct.layout.GridCanvasGenerator;

public class ZoomToFitDrawer {
    private DirectMoleculeDrawer moleculeDrawer;
    private CanvasGenerator canvasGenerator;
    private Params params;

    public ZoomToFitDrawer() {
        this(new DirectMoleculeDrawer(), new GridCanvasGenerator());
    }

    public ZoomToFitDrawer(DirectMoleculeDrawer moleculeDrawer, CanvasGenerator canvasGenerator) {
        this.moleculeDrawer = moleculeDrawer;
        this.params = moleculeDrawer.getParams();
        this.canvasGenerator = canvasGenerator;
    }

    public void draw(List<IAtomContainer> mols, Dimension cellCanvas, Graphics2D g) {
        this.canvasGenerator.layout(mols, cellCanvas);
        AffineTransform originalTransform = g.getTransform();
        for (IAtomContainer mol : mols) {
            Rectangle2D canvas = this.canvasGenerator.getCanvasForAtomContainer(mol);
            g.translate(canvas.getCenterX(), canvas.getCenterY());
            double zoom = this.calculateZoom(mol, canvas);
            g.scale(zoom, zoom);
            this.moleculeDrawer.drawMolecule(mol, g);
            g.setTransform(originalTransform);
        }
    }

    private double calculateZoom(IAtomContainer ac, Rectangle2D canvas) {
        double scaleFactor = GeometryTools.getScaleFactor((IAtomContainer)ac, (double)this.params.bondLength);
        GeometryTools.translate2DCenterTo((IAtomContainer)ac, (Point2d)new Point2d(0.0, 0.0));
        GeometryTools.scaleMolecule((IAtomContainer)ac, (double)scaleFactor);
        Rectangle2D r2D = GeometryTools.getRectangle2D((IAtomContainer)ac);
        double canvasWidth = canvas.getWidth();
        double canvasHeight = canvas.getHeight();
        double objectWidth = r2D.getWidth() + (double)(this.params.borderX * 2);
        double objectHeight = r2D.getHeight() + (double)(this.params.borderY * 2);
        return Math.min(canvasWidth / objectWidth, canvasHeight / objectHeight);
    }
}

