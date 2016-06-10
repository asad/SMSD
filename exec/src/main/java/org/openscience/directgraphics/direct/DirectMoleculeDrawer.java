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
import java.awt.Font;
import java.awt.Graphics2D;
import java.awt.geom.Rectangle2D;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import javax.vecmath.Point2f;
import org.openscience.cdk.geometry.GeometryTools;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.directgraphics.stereo.IStereoAndConformation;

public class DirectMoleculeDrawer
        extends AbstractDirectDrawer {

    private Font moleculeIDFont;
    private List<Highlighter> highlightDrawers;
    private LabelManager labelManager;
    private DirectAtomDrawer atomDrawer;
    private DirectBondDrawer bondDrawer;
    private Map<IAtom, IStereoAndConformation> chiralMap;

    public DirectMoleculeDrawer(Params params) {
        super.setParams(params);
        params.bondLength = 20;
        this.highlightDrawers = new ArrayList<>();
        AbstractHighlightDrawer highlightDrawer = params.useCircularHighlight ? new OutlineHighlighter(params) : new SimpleHighlighter(params);
        this.highlightDrawers.add((Highlighter) ((Object) highlightDrawer));
        this.labelManager = new LabelManager();
        this.atomDrawer = new DirectAtomDrawer(params, this.labelManager);
        this.bondDrawer = new DirectBondDrawer(params, this.labelManager);
        this.chiralMap = new HashMap<>();
    }

    public DirectMoleculeDrawer() {
        this(new Params());
    }

    public void addToChiralMap(Map<IAtom, IStereoAndConformation> chirals) {
        this.chiralMap.putAll(chirals);
    }

    public Rectangle2D getDrawnBounds(List<IAtom> atoms) {
        return this.atomDrawer.getDrawnBounds(atoms);
    }

    public void clearHighlights() {
        for (Highlighter highlightDrawer : this.highlightDrawers) {
            highlightDrawer.clearHighlights();
        }
    }

    public Highlighter getFirstHighlighter() {
        Highlighter highlightDrawer;
        if (this.highlightDrawers.isEmpty()) {
            highlightDrawer = this.params.useCircularHighlight ? new OutlineHighlighter(this.params) : new SimpleHighlighter(this.params);
            this.highlightDrawers.add(highlightDrawer);
        } else {
            highlightDrawer = this.highlightDrawers.get(0);
        }
        return highlightDrawer;
    }

    public List<Highlighter> getHighlighters() {
        return this.highlightDrawers;
    }

    public void addHighlighter(Highlighter highlighter) {
        this.highlightDrawers.add(highlighter);
    }

    public void addHighlights(IAtomContainer highlightContainer, Color color) {
        Highlighter highlightDrawer = this.getFirstHighlighter();
        highlightDrawer.addHighlights(highlightContainer, color);
    }

    public void addHighlights(List<IAtom> atoms, Color color) {
        HashMap<IAtom, Color> atomColorMap = new HashMap<>();
        for (IAtom atom : atoms) {
            atomColorMap.put(atom, color);
        }
        Highlighter highlightDrawer = this.getFirstHighlighter();
        highlightDrawer.addToHighlights(atomColorMap);
    }

    public void addHighlights(IAtomContainer highlightContainer) {
        this.addHighlights(highlightContainer, this.params.highlightColor);
    }

    public void addHighlights(List<IAtom> atoms, List<IBond> bonds) {
        Highlighter highlightDrawer = this.getFirstHighlighter();
        highlightDrawer.addHighlights(atoms, bonds);
    }

    public void addHighlights(List<IAtom> atoms) {
        this.addHighlights(atoms, new ArrayList<IBond>());
    }

    public void addToHighlights(Map<IAtom, Color> colorMap) {
        Highlighter highlightDrawer = this.getFirstHighlighter();
        highlightDrawer.addToHighlights(colorMap);
    }

    public void drawMolecule(IAtomContainer molecule, Graphics2D g) {
        this.labelManager.reset();
        this.atomDrawer.setAtomSymbolFont(new Font("ROMAN", 0, this.params.atomSymbolFontSize));
        this.atomDrawer.setSubscriptFont(new Font("ROMAN", 0, this.params.subscriptTextSize));
        this.atomDrawer.setAtomIDFont(new Font("ROMAN", 0, this.params.atomIDFontSize));
        this.atomDrawer.setChiralSymbolFont(new Font("ROMAN", 0, this.params.chiralSymbolFontSize));
        this.moleculeIDFont = new Font("ROMAN", 0, this.params.moleculeLabelFontSize);
        Color savedColor = g.getColor();
        if (this.params.drawBounds) {
            Rectangle2D bounds = GeometryTools.getRectangle2D((IAtomContainer) molecule);
            g.draw(bounds);
        }
        if (this.params.drawHighlights && this.params.highlightsBelow) {
            this.drawHighlights(molecule, g);
        }
        this.atomDrawer.setChirals(this.chiralMap);
        this.bondDrawer.drawBonds(molecule, g);
        this.atomDrawer.drawAtoms(molecule, g);
        if (this.params.drawHighlights && this.params.highlightsAbove) {
            this.drawHighlights(molecule, g);
        }
        if (this.params.drawMoleculeID) {
            this.drawMoleculeID(molecule, g);
        }
        g.setColor(savedColor);
    }

    private void drawHighlights(IAtomContainer molecule, Graphics2D g) {
        for (Highlighter highlightDrawer : this.highlightDrawers) {
            highlightDrawer.drawHighlights(molecule, g);
        }
    }

    public Rectangle2D drawMoleculeID(IAtomContainer atomContainer, Graphics2D g) {
        String id = atomContainer.getID();
        if (id == null) {
            return null;
        }
        Rectangle2D moleculeBounds = GeometryTools.getRectangle2D((IAtomContainer) atomContainer);
        double labelCenterX = moleculeBounds.getCenterX();
        double labelCenterY = moleculeBounds.getMaxY() + this.params.labelYGap;
        Point2f textPoint = this.getTextPoint(g, id, labelCenterX, labelCenterY);
        g.setFont(this.moleculeIDFont);
        g.setColor(Color.BLACK);
        g.drawString(id, textPoint.x, textPoint.y);
        Rectangle2D textBounds = this.getTextBounds(g, id);
        return new Rectangle2D.Double(labelCenterX - textBounds.getWidth() / 2.0, labelCenterY - textBounds.getHeight() / 2.0, textBounds.getWidth(), textBounds.getHeight());
    }
}
