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

import java.awt.Font;
import java.awt.Graphics2D;
import java.awt.geom.Rectangle2D;
import javax.vecmath.Point2d;
import javax.vecmath.Tuple2d;
import javax.vecmath.Vector2d;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.directgraphics.direct.Params;
import org.openscience.directgraphics.direct.layout.BoundsTree;

public class LinearAtomContainerSetLayout
extends AbstractAWTLayout<IAtomContainerSet> {
    private Vector2d moleculeSetAxis;
    private MoleculeLayout moleculeLayout;

    public LinearAtomContainerSetLayout(Vector2d moleculeSetAxis) {
        this(moleculeSetAxis, new Params());
    }

    public LinearAtomContainerSetLayout(Vector2d moleculeSetAxis, Params params) {
        this.params = params;
        this.moleculeLayout = new MoleculeLayout(params);
        this.moleculeSetAxis = moleculeSetAxis;
    }

    @Override
    public BoundsTree layout(IAtomContainerSet atomContainerSet, Graphics2D graphics) {
        Font plusFont = new Font("ROMAN", 0, this.params.plusFontSize);
        graphics.setFont(plusFont);
        Rectangle2D plusBounds = super.getTextBounds(graphics, "+");
        double molGap = (double)(2 * this.params.plusGap) + plusBounds.getWidth();
        String atomContainerSetID = atomContainerSet.getID();
        this.boundsTree = new BoundsTree(atomContainerSetID);
        Point2d curr = new Point2d(0.0, 0.0);
        int moleculeCounter = 0;
        for (IAtomContainer molecule : atomContainerSet.atomContainers()) {
            String label = molecule.getID();
            label = label == null || label.equals("") ? "mol" + String.valueOf(moleculeCounter) : label + ":" + String.valueOf(moleculeCounter);
            BoundsTree molBounds = this.moleculeLayout.layout(molecule, label, graphics);
            double boundsWidth = molBounds.getWidth();
            double halfBoundsWidth = boundsWidth / 2.0;
            curr.scaleAdd(halfBoundsWidth, (Tuple2d)this.moleculeSetAxis, (Tuple2d)curr);
            this.translateTo(molecule, curr.x, curr.y, molBounds);
            curr.scaleAdd(halfBoundsWidth, (Tuple2d)this.moleculeSetAxis, (Tuple2d)curr);
            curr.scaleAdd(molGap, (Tuple2d)this.moleculeSetAxis, (Tuple2d)curr);
            this.boundsTree.add(atomContainerSetID, molBounds);
            ++moleculeCounter;
        }
        return this.boundsTree;
    }

    @Override
    public BoundsTree layout(IAtomContainerSet obj, String rootLabel, Graphics2D graphics) {
        return null;
    }
}

