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
package org.openscience.directgraphics.direct;

import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.geom.Ellipse2D;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import javax.vecmath.Point2d;
import org.openscience.cdk.geometry.GeometryTools;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;

public class OutlineHighlighter
        extends AbstractHighlightDrawer
        implements Highlighter {

    private final Map<IAtomContainer, Color> colorMap = new HashMap<>();

    public OutlineHighlighter(Params params) {
        super(params);
    }

    @Override
    public void addHighlights(IAtomContainer highlightContainer, Color color) {
        this.colorMap.put(highlightContainer, color);
    }

    @Override
    public void addHighlights(List<IAtom> atoms, List<IBond> bonds) {
        IAtomContainer highlightContainer = null;
        if (atoms.size() > 0) {
            highlightContainer = (IAtomContainer) atoms.get(0).getBuilder().newInstance((Class) IAtomContainer.class, new Object[0]);
        } else if (bonds.size() > 0) {
            highlightContainer = (IAtomContainer) bonds.get(0).getBuilder().newInstance((Class) IAtomContainer.class, new Object[0]);
        } else {
            return;
        }
        for (IAtom atom : atoms) {
            highlightContainer.addAtom(atom);
        }
        for (IBond bond : bonds) {
            highlightContainer.addBond(bond);
        }
        this.addHighlights(highlightContainer, this.params.highlightColor);
    }

    @Override
    public void addToHighlights(Map<IAtom, Color> colorMap) {
    }

    @Override
    public void drawHighlights(IAtomContainer molecule, Graphics2D g) {
        ArrayList<IAtomContainer> highlightContainers;
        Point2d center = null;
        if (this.params.circularHighlightIsConcentric) {
            highlightContainers = new ArrayList<>(this.colorMap.keySet());
            Collections.sort(highlightContainers, (IAtomContainer ac0, IAtomContainer ac1) -> {
                if (ac0.getAtomCount() < ac1.getAtomCount()) {
                    return 1;
                }
                if (ac0.getAtomCount() > ac1.getAtomCount()) {
                    return -1;
                }
                return 0;
            });
            center = GeometryTools.get2DCenter((IAtomContainer) highlightContainers.get(highlightContainers.size() - 1));
        } else {
            highlightContainers = new ArrayList<>(this.colorMap.keySet());
        }
        for (int containerIndex = 0; containerIndex < highlightContainers.size(); ++containerIndex) {
            double x;
            double dim;
            double y;
            IAtomContainer highlightContainer = highlightContainers.get(containerIndex);
            Color savedColor = g.getColor();
            if (this.params.circularHighlightTransparentFilled) {
                g.setColor(this.getTranslucentColor(this.colorMap.get(highlightContainer)));
            } else {
                g.setColor(this.colorMap.get(highlightContainer));
            }
            if (!this.params.circularHighlightIsConcentric || center == null) {
                center = GeometryTools.get2DCenter((IAtomContainer) highlightContainer);
            }
            double maxDist = 0.0;
            for (IAtom highlightAtom : highlightContainer.atoms()) {
                double d;
                if (!molecule.contains(highlightAtom)) {
                    continue;
                }
                Point2d point = highlightAtom.getPoint2d();
                if (point != null && (d = center.distance(point)) > maxDist) {
                    maxDist = d;
                }
                if (!this.params.circularHighlightShowAtoms) {
                    continue;
                }
                double r = this.params.highlightRadius;
                g.fill(new Ellipse2D.Double(point.x - r, point.y - r, r * 2.0, r * 2.0));
            }
            if (highlightContainer.getAtomCount() == 1 && containerIndex == highlightContainers.size() - 1) {
                x = center.x - this.params.circularHighlightMinRadius;
                y = center.y - this.params.circularHighlightMinRadius;
                dim = 2.0 * this.params.circularHighlightMinRadius;
            } else {
                x = center.x - maxDist;
                y = center.y - maxDist;
                dim = 2.0 * maxDist;
            }
            if (this.params.circularHighlightTransparentFilled) {
                g.fill(new Ellipse2D.Double(x, y, dim, dim));
            } else {
                g.draw(new Ellipse2D.Double(x, y, dim, dim));
            }
            g.setColor(savedColor);
        }
    }

    @Override
    public void clearHighlights() {
        this.colorMap.clear();
    }

}
