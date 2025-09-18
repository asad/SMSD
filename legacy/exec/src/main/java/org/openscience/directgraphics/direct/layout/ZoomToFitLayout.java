/*
 * SPDX-License-Identifier: Apache-2.0
 * © 2025 BioInception PVT LTD.
 */
/*
 * Decompiled with CFR 0_114.
 * 
 * Could not load the following classes:
 *  javax.vecmath.Vector2d
 *  org.openscience.cdk.geometry.GeometryTools
 *  org.openscience.cdk.interfaces.IAtomContainer
 */
package org.openscience.directgraphics.direct.layout;

import java.awt.Dimension;
import java.awt.Graphics2D;
import java.awt.geom.AffineTransform;
import java.awt.geom.Rectangle2D;
import javax.vecmath.Vector2d;
import org.openscience.cdk.geometry.GeometryTools;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.directgraphics.direct.DirectMoleculeDrawer;
import org.openscience.directgraphics.direct.Params;

public class ZoomToFitLayout
extends AbstractDirectLayout<IAtomContainer> {
    private final DirectMoleculeDrawer drawer;

    public ZoomToFitLayout(DirectMoleculeDrawer drawer) {
        this.drawer = drawer;
    }

    public void layout(IAtomContainer mol, Dimension cellCanvas, Graphics2D g) {
        AffineTransform originalTransform = g.getTransform();
        double w = cellCanvas.width;
        double h = cellCanvas.height;
        double zoom = this.calculateZoom(mol, cellCanvas.width, cellCanvas.height);
        double centerX = w / 2.0;
        Params paramsLocal = this.drawer.getParams();
        double centerY = paramsLocal.drawMoleculeID ? h / 2.0 - paramsLocal.labelYGap : h / 2.0;
        g.translate(centerX, centerY);
        g.scale(zoom, zoom);
        this.drawer.drawMolecule(mol, g);
        g.setTransform(originalTransform);
    }

    private double calculateZoom(IAtomContainer ac, double w, double h) {
        double borderX = this.drawer.getParams().borderX;
        double borderY = this.drawer.getParams().borderY;
        double canvasWidth = w;
        double canvasHeight = h;
        double scaleFactor = GeometryTools.getScaleFactor((IAtomContainer)ac, (double)this.drawer.getParams().bondLength);
        Rectangle2D r2D = GeometryTools.getRectangle2D((IAtomContainer)ac);
        this.translateTo(ac, 0.0, 0.0, r2D);
        GeometryTools.scaleMolecule((IAtomContainer)ac, (double)scaleFactor);
        double objectWidth = r2D.getWidth() + borderX * 2.0;
        double objectHeight = r2D.getHeight() + borderY * 2.0;
        return Math.min(canvasWidth / objectWidth, canvasHeight / objectHeight);
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

