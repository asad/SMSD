/*
 * SPDX-License-Identifier: Apache-2.0
 * © 2025 BioInception PVT LTD.
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
