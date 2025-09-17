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

import java.awt.Graphics2D;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.directgraphics.direct.LabelManager;
import org.openscience.directgraphics.direct.Params;
import org.openscience.directgraphics.direct.layout.BoundsTree;

public class MoleculeLayout
        extends AbstractAWTLayout<IAtomContainer> {

    private AtomLayout atomLayout;

    public MoleculeLayout(Params params) {
        this.atomLayout = new AtomLayout(this, params, new LabelManager());
    }

    public MoleculeLayout(AbstractAWTLayout parent, Params params) {
        this(params);
        this.parent = parent;
    }

    @Override
    public BoundsTree layout(IAtomContainer atomContainer, Graphics2D graphics) {
        return this.layout(atomContainer, atomContainer.getID(), graphics);
    }

    @Override
    public BoundsTree layout(IAtomContainer atomContainer, String rootLabel, Graphics2D graphics) {
        this.atomLayout.reset();
        this.setGraphics(graphics);
        this.currentObject = atomContainer;
        this.boundsTree = new BoundsTree(rootLabel);
        for (IAtom atom : atomContainer.atoms()) {
            this.boundsTree.add(rootLabel, this.atomLayout.layout(atom, graphics));
        }
        return this.boundsTree;
    }
}
