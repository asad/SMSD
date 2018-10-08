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
package org.openscience.directgraphics.direct.layout;

import java.awt.geom.Rectangle2D;
import javax.vecmath.Vector2d;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.directgraphics.direct.Params;

public abstract class AbstractDirectLayout<T> {
    public static final String INVERTED = "Coordinates Inverted";
    protected Params params;
    protected BoundsTree boundsTree;
    public boolean shouldInvert;

    public abstract BoundsTree layout(T var1, Vector2d var2);

    public abstract Vector2d getAxis();

    public abstract double getAxisPosition();

    public AbstractDirectLayout() {
        this(true);
    }

    public AbstractDirectLayout(boolean shouldInvert) {
        this.shouldInvert = shouldInvert;
    }

    public Params getParams() {
        return this.params;
    }

    public void setParams(Params params) {
        this.params = params;
    }

    public void translateTo(IAtomContainer ac, double x, double y, Rectangle2D bounds) {
        double dx = x - bounds.getCenterX();
        double dy = y - bounds.getCenterY();
        for (IAtom atom : ac.atoms()) {
            atom.getPoint2d().x += dx;
            atom.getPoint2d().y += dy;
        }
        bounds.setRect(bounds.getMinX() + dx, bounds.getMinY() + dy, bounds.getWidth(), bounds.getHeight());
    }

    public void invert(IAtomContainer ac) {
        if (this.shouldInvert && ac.getProperty((Object)"Coordinates Inverted") == null || !((Boolean)ac.getProperty((Object)"Coordinates Inverted"))) {
            for (IAtom atom : ac.atoms()) {
                atom.getPoint2d().y *= -1.0;
            }
            ac.setProperty((Object)"Coordinates Inverted", (Object)Boolean.TRUE);
        }
        this.shouldInvert = false;
    }

    public void align(IAtomContainer atomContainer, Vector2d molAxis) {
        switch (this.params.moleculeAlignMethod) {
            case MAX_AXIS: {
                MoleculeAligner.alignToMaxWidth(atomContainer, molAxis);
            }
            case MIN_AREA: {
                MoleculeAligner.alignToMinAreaBox(atomContainer, molAxis);
            }
        }
        MoleculeAligner.alignToMaxWidth(atomContainer, molAxis);
    }

}

