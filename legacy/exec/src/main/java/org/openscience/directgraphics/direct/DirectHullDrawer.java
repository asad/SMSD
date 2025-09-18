/*
 * SPDX-License-Identifier: Apache-2.0
 * © 2025 BioInception PVT LTD.
 */
package org.openscience.directgraphics.direct;

import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.Image;
import java.awt.RenderingHints;
import java.awt.geom.Ellipse2D;
import java.awt.image.BufferedImage;
import javax.vecmath.Point2d;
import javax.vecmath.Tuple2d;
import javax.vecmath.Vector2d;
import org.openscience.cdk.interfaces.IAtomContainer;

public class DirectHullDrawer
        extends AbstractDirectDrawer {

    private final DirectMoleculeDrawer moleculeDrawer = new DirectMoleculeDrawer();

    public DirectHullDrawer() {
        super.params = this.moleculeDrawer.getParams();
    }

    public Image drawHull(IAtomContainer atomContainer, int w, int h) {
        BufferedImage image = super.makeBlankImage(w, h);
        Graphics2D g = (Graphics2D) image.getGraphics();
        if (this.params.useAntialias) {
            g.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
        }
        for (int i = 0; i < atomContainer.getAtomCount(); ++i) {
            atomContainer.getAtom(i).setID(String.valueOf(i));
        }
        this.params.drawAtomID = true;
        this.drawHull(atomContainer, g);
        return image;
    }

    public void drawHull(IAtomContainer atomContainer, Graphics2D g) {
        DirectArrowDrawer arrowDrawer = new DirectArrowDrawer(this.getParams());
        ConvexHull hull = new ConvexHull(atomContainer);
        this.moleculeDrawer.drawMolecule(atomContainer, g);
        Point2d prev = null;
        Point2d first = null;
        for (Point2d hullPoint : hull) {
            if (prev == null) {
                first = prev = hullPoint;
                continue;
            }
            g.setColor(Color.RED);
            this.drawLine(prev, hullPoint, g);
            g.setColor(Color.BLACK);
            Point2d midPoint = new Point2d(prev);
            midPoint.interpolate((Tuple2d) hullPoint, 0.5);
            Vector2d direction = new Vector2d((Tuple2d) hullPoint);
            direction.sub((Tuple2d) prev);
            direction.normalize();
            arrowDrawer.drawArrow(g, midPoint, direction);
            prev = hullPoint;
        }
        g.setColor(Color.RED);
        this.drawLine(first, prev, g);
        g.setColor(Color.BLACK);
        Point2d midPoint = new Point2d(prev);
        midPoint.interpolate((Tuple2d) first, 0.5);
        Vector2d direction = new Vector2d((Tuple2d) first);
        direction.sub((Tuple2d) prev);
        direction.normalize();
        arrowDrawer.drawArrow(g, midPoint, direction);
        ConvexHull.Rectangle r = hull.getMinimumAreaBoundingRectangleBruteForce();
        Vector2d majorAxis = r.getMajorAxis();
        majorAxis.normalize();
        Point2d center = hull.getCenter();
        arrowDrawer.drawArrow(g, center, majorAxis);
        g.setColor(Color.BLACK);
        this.drawLine(r.cornerA, r.cornerB, g);
        this.drawLine(r.cornerB, r.cornerC, g);
        this.drawLine(r.cornerC, r.cornerD, g);
        this.drawLine(r.cornerD, r.cornerA, g);
        g.setColor(Color.BLUE);
        g.fill(new Ellipse2D.Double(r.cornerA.x - 3.0, r.cornerA.y - 3.0, 6.0, 6.0));
        g.setColor(Color.MAGENTA);
        g.fill(new Ellipse2D.Double(r.cornerB.x - 3.0, r.cornerB.y - 3.0, 6.0, 6.0));
        g.setColor(Color.YELLOW);
        g.fill(new Ellipse2D.Double(r.cornerC.x - 3.0, r.cornerC.y - 3.0, 6.0, 6.0));
        g.setColor(Color.CYAN);
        g.fill(new Ellipse2D.Double(r.cornerD.x - 3.0, r.cornerD.y - 3.0, 6.0, 6.0));
        g.setColor(Color.GREEN);
        g.fill(new Ellipse2D.Double(r.pointX.x - 2.0, r.pointX.y - 2.0, 4.0, 4.0));
        g.setColor(Color.PINK);
        g.fill(new Ellipse2D.Double(r.pointY.x - 2.0, r.pointY.y - 2.0, 4.0, 4.0));
        g.setColor(Color.ORANGE);
        g.fill(new Ellipse2D.Double(r.pointZ.x - 2.0, r.pointZ.y - 2.0, 4.0, 4.0));
    }
}
