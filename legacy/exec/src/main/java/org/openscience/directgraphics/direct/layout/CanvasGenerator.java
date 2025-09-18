/*
 * SPDX-License-Identifier: Apache-2.0
 * © 2025 BioInception PVT LTD.
 */
package org.openscience.directgraphics.direct.layout;

import java.awt.Dimension;
import java.awt.geom.Rectangle2D;
import java.util.List;
import org.openscience.cdk.interfaces.IAtomContainer;

public interface CanvasGenerator {

    public void layout(List<IAtomContainer> var1, Dimension var2);

    public Rectangle2D getCanvasForAtomContainer(IAtomContainer var1);

    public Dimension getSize();
}
