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
package org.openscience.directgraphics.direct.layout;

import java.awt.Graphics2D;
import java.awt.geom.Rectangle2D;
import java.util.ArrayList;
import java.util.List;
import javax.vecmath.Point2d;
import javax.vecmath.Tuple2d;
import javax.vecmath.Vector2d;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.directgraphics.direct.DirectArrowDrawer;

public class ArrowWheel {

    private DirectArrowDrawer arrowDrawer;
    private List<Arrow> arrows;
    private IAtomContainer hub;
    private List<String> arrowLabels;

    public ArrowWheel(DirectArrowDrawer arrowDrawer, IAtomContainer hubMolecule, List<IAtomContainer> rimMolecules) {
        this(arrowDrawer, hubMolecule, rimMolecules, new ArrayList<String>());
    }

    public ArrowWheel(DirectArrowDrawer arrowDrawer, IAtomContainer hubMolecule, List<IAtomContainer> rimMolecules, List<String> arrowLabels) {
        this.arrowDrawer = arrowDrawer;
        this.arrows = new ArrayList<>();
        this.hub = hubMolecule;
        for (IAtomContainer molecule : rimMolecules) {
            this.arrows.add(new Arrow(this.hub, molecule));
        }
        this.arrowLabels = arrowLabels;
    }

    public void draw(CanvasGenerator canvasGenerator, Graphics2D g) {
        for (Arrow arrow : this.arrows) {
            Rectangle2D tailCanvas = canvasGenerator.getCanvasForAtomContainer(arrow.tail);
            Rectangle2D headCanvas = canvasGenerator.getCanvasForAtomContainer(arrow.head);
            Point2d tailCenter = new Point2d(tailCanvas.getCenterX(), tailCanvas.getCenterY());
            Point2d headCenter = new Point2d(headCanvas.getCenterX(), headCanvas.getCenterY());
            arrow.center = new Point2d(tailCenter);
            arrow.center.interpolate((Tuple2d) headCenter, 0.5);
            arrow.vector = new Vector2d((Tuple2d) headCenter);
            arrow.vector.sub((Tuple2d) tailCenter);
            arrow.vector.normalize();
        }
        int i = 0;
        for (Arrow arrow2 : this.arrows) {
            this.arrowDrawer.drawThinArrow(g, arrow2.center, arrow2.vector, this.arrowLabels.get(i));
            ++i;
        }
    }

    private class Arrow {

        public IAtomContainer tail;
        public IAtomContainer head;
        public Point2d center;
        public Vector2d vector;

        public Arrow(IAtomContainer tail, IAtomContainer head) {
            this.tail = tail;
            this.head = head;
        }
    }

}
