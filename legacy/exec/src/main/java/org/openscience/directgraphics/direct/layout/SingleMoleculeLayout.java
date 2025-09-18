/*
 * SPDX-License-Identifier: Apache-2.0
 * © 2025 BioInception PVT LTD.
 */
package org.openscience.directgraphics.direct.layout;

import java.awt.geom.Rectangle2D;
import javax.vecmath.Point2d;
import javax.vecmath.Tuple2d;
import javax.vecmath.Vector2d;
import org.openscience.cdk.AtomContainer;
import org.openscience.cdk.geometry.GeometryTools;
import org.openscience.cdk.graph.ConnectivityChecker;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.layout.StructureDiagramGenerator;
import org.openscience.directgraphics.direct.Params;

public class SingleMoleculeLayout
extends AbstractDirectLayout<IAtomContainer> {
    private StructureDiagramGenerator sdg;
    private boolean forceRelayout;

    public SingleMoleculeLayout(Params params) {
        this(params, false);
    }

    public SingleMoleculeLayout(Params params, boolean forceRelayout) {
        this.setParams(params);
        this.sdg = new StructureDiagramGenerator();
        this.forceRelayout = forceRelayout;
    }

    @Override
    public BoundsTree layout(IAtomContainer atomContainer, Vector2d axis) {
        Point2d center = new Point2d((Tuple2d)axis);
        if (this.forceRelayout || !GeometryTools.has2DCoordinates((IAtomContainer)atomContainer)) {
            this.sdg.setMolecule((IAtomContainer)new AtomContainer(atomContainer), false);
            try {
                if (ConnectivityChecker.isConnected((IAtomContainer)atomContainer)) {
                    this.sdg.generateCoordinates();
                } else {
                    System.err.println("Disconnected components needs to be layout separately");
                }
            }
            catch (Exception e) {
                e.printStackTrace();
            }
        }
        double scale = GeometryTools.getScaleFactor((IAtomContainer)atomContainer, (double)this.params.bondLength);
        Rectangle2D bounds = GeometryTools.getRectangle2D((IAtomContainer)atomContainer);
        GeometryTools.scaleMolecule((IAtomContainer)atomContainer, (double)scale);
        this.translateTo(atomContainer, center.x, center.y, bounds);
        String label = atomContainer.getID();
        return new BoundsTree(label, label, bounds);
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

