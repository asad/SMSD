/*
 * SPDX-License-Identifier: Apache-2.0
 * © 2025 BioInception PVT LTD.
 */
package org.openscience.directgraphics.direct;

import java.awt.Color;
import java.awt.Graphics2D;
import java.util.List;
import java.util.Map;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;

public interface Highlighter {

    public void addHighlights(IAtomContainer var1, Color var2);

    public void addHighlights(List<IAtom> var1, List<IBond> var2);

    public void drawHighlights(IAtomContainer var1, Graphics2D var2);

    public void addToHighlights(Map<IAtom, Color> var1);

    public void clearHighlights();
}
