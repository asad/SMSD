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
package org.openscience.directgraphics.direct;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.Stroke;
import java.awt.geom.Ellipse2D;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import javax.vecmath.Point2d;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;

public class SimpleHighlighter
extends AbstractHighlightDrawer
implements Highlighter {
    private Map<IAtom, Color> atomColorMap = new HashMap<>();
    private final Map<IBond, Color> bondColorMap = new HashMap<>();

    public SimpleHighlighter(Params params) {
        super(params);
    }

    @Override
    public void drawHighlights(IAtomContainer molecule, Graphics2D g) {
        Color color;
        for (IAtom atom : this.atomColorMap.keySet()) {
            if (!molecule.contains(atom)) continue;
            color = this.atomColorMap.get(atom);
            this.drawHighlight(atom, color, g);
        }
        for (IBond bond : this.bondColorMap.keySet()) {
            if (!molecule.contains(bond)) continue;
            color = this.bondColorMap.get(bond);
            this.drawHighlight(bond, color, g);
        }
    }

    @Override
    public void addHighlights(IAtomContainer highlightContainer, Color color) {
        this.registerColor(color);
        for (IAtom atom : highlightContainer.atoms()) {
            this.atomColorMap.put(atom, color);
        }
        for (IBond bond : highlightContainer.bonds()) {
            this.bondColorMap.put(bond, color);
        }
    }

    @Override
    public void addHighlights(List<IAtom> atoms, List<IBond> bonds) {
        for (IAtom atom : atoms) {
            this.atomColorMap.put(atom, this.params.highlightColor);
        }
        for (IBond bond : bonds) {
            this.bondColorMap.put(bond, this.params.highlightColor);
        }
    }

    @Override
    public void addToHighlights(Map<IAtom, Color> atomColorMap) {
        this.atomColorMap.putAll(atomColorMap);
    }

    public void setHighlights(Map<IAtom, Color> atomColorMap) {
        this.atomColorMap = atomColorMap;
    }

    public void drawHighlight(IAtom atom, Graphics2D g) {
        if (this.params.highlightsAbove) {
            this.drawHighlight(atom, this.translucentHighlightColor, g);
        } else {
            this.drawHighlight(atom, this.opaqueHighlightColor, g);
        }
    }

    public void drawHighlight(IAtom atom, Color color, Graphics2D g) {
        Color actualColor = this.params.highlightsAbove ? this.getTranslucentColor(color) : color;
        g.setColor(actualColor);
        double r = this.params.highlightRadius;
        double d = r * 2.0;
        Point2d p = atom.getPoint2d();
        g.fill(new Ellipse2D.Double(p.x - r, p.y - r, d, d));
    }

    public void drawHighlight(IBond bond, Color color, Graphics2D g) {
        Stroke stroke = g.getStroke();
        g.setStroke(new BasicStroke(this.params.highlightBondStroke));
        Point2d p0 = bond.getAtom(0).getPoint2d();
        Point2d p1 = bond.getAtom(1).getPoint2d();
        this.drawLine(p0, p1, g);
        g.setStroke(stroke);
    }

    public void drawHighlightContainer(IAtomContainer highlightContainer, Graphics2D g) {
        if (this.params.highlightsAbove) {
            this.drawHighlightContainer(highlightContainer, this.translucentHighlightColor, g);
        } else {
            this.drawHighlightContainer(highlightContainer, this.opaqueHighlightColor, g);
        }
    }

    public void drawHighlightContainer(IAtomContainer highlightContainer, Color color, Graphics2D g) {
        Color actualColor = this.params.highlightsAbove ? this.getTranslucentColor(color) : color;
        g.setColor(actualColor);
        double r = this.params.highlightRadius;
        double d = r * 2.0;
        for (IAtom atom : highlightContainer.atoms()) {
            Point2d p = atom.getPoint2d();
            g.fill(new Ellipse2D.Double(p.x - r, p.y - r, d, d));
        }
        Stroke stroke = g.getStroke();
        g.setStroke(new BasicStroke(this.params.highlightBondStroke));
        for (IBond bond : highlightContainer.bonds()) {
            Point2d p0 = bond.getAtom(0).getPoint2d();
            Point2d p1 = bond.getAtom(1).getPoint2d();
            this.drawLine(p0, p1, g);
        }
        g.setStroke(stroke);
    }

    @Override
    public void clearHighlights() {
        this.atomColorMap.clear();
        this.bondColorMap.clear();
    }
}

