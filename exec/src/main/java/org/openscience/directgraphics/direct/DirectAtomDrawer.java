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
import java.awt.Font;
import java.awt.Graphics2D;
import java.awt.Stroke;
import java.awt.geom.Ellipse2D;
import java.awt.geom.Rectangle2D;
import java.util.BitSet;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import javax.vecmath.Point2d;
import javax.vecmath.Point2f;
import javax.vecmath.Tuple2d;
import javax.vecmath.Vector2d;
import org.openscience.cdk.PseudoAtom;
import org.openscience.cdk.geometry.GeometryTools;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.ILonePair;
import org.openscience.cdk.renderer.color.CDK2DAtomColors;
import org.openscience.cdk.renderer.color.IAtomColorer;
import org.openscience.directgraphics.stereo.IStereoAndConformation;

public class DirectAtomDrawer
        extends AbstractDirectDrawer {

    private Font atomSymbolFont;
    private Font subscriptFont;
    private Font atomIDFont;
    private Font chiralSymbolFont;
    private final IAtomColorer atomColorer;
    private final Map<IAtom, Rectangle2D> drawnAtomBounds;
    private final LabelManager labelManager;
    private Map<IAtom, IStereoAndConformation> chiralMap;

    /**
     *
     * @param params
     * @param labelManager
     */
    public DirectAtomDrawer(Params params, LabelManager labelManager) {
        super.params = params;
        this.labelManager = labelManager;
        this.atomColorer = new CDK2DAtomColors();
        this.drawnAtomBounds = new HashMap<>();
        this.chiralMap = new HashMap<>();
    }

    public void setChirals(Map<IAtom, IStereoAndConformation> chiralMap) {
        this.chiralMap = chiralMap;
    }

    public void setAtomSymbolFont(Font atomSymbolFont) {
        this.atomSymbolFont = atomSymbolFont;
    }

    public void setSubscriptFont(Font subscriptFont) {
        this.subscriptFont = subscriptFont;
    }

    public void setAtomIDFont(Font atomIDFont) {
        this.atomIDFont = atomIDFont;
    }

    public void setChiralSymbolFont(Font chiralSymbolFont) {
        this.chiralSymbolFont = chiralSymbolFont;
    }

    public Rectangle2D getDrawnBounds(List<IAtom> atoms) {
        Rectangle2D totalBounds = null;
        for (IAtom atom : atoms) {
            Rectangle2D bounds = this.drawnAtomBounds.get(atom);
            if (bounds == null) {
                continue;
            }
            if (totalBounds == null) {
                totalBounds = (Rectangle2D) bounds.clone();
            }
            totalBounds.add(bounds);
        }
        return totalBounds;
    }

    public void drawAtoms(IAtomContainer molecule, Graphics2D g) {
        Map<IAtom, Integer> lonePairMap = null;
        if (this.params.drawLonePairs) {
            lonePairMap = this.getLonePairCounts(molecule);
        }
        for (IAtom atom : molecule.atoms()) {
            int lonePairCount = 0;
            if (this.params.drawLonePairs) {
                Integer lonePairCountInteger = lonePairMap.get(atom);
                lonePairCount = lonePairCountInteger == null ? Integer.valueOf(0).intValue() : lonePairCountInteger.intValue();
            }
            this.drawnAtomBounds.put(atom, this.drawAtom(atom, molecule, lonePairCount, g));
        }
    }

    public Rectangle2D drawAtom(IAtom atom, IAtomContainer molecule, int lonePairCount, Graphics2D g) {
        Rectangle2D symbolBounds;
        Rectangle2D idBounds;
        if (this.shouldDraw(atom, molecule)) {
            Integer implicitHydrogenCount;
            symbolBounds = this.drawAtomSymbol(atom, g);
            if (this.isCharged(atom)) {
                Rectangle2D chargeBounds = this.drawCharge(atom, g);
                symbolBounds.add(chargeBounds);
            }
            if (this.params.drawImplicitHydrogens && (implicitHydrogenCount = atom.getImplicitHydrogenCount()) != null && implicitHydrogenCount > 0) {
                int align = GeometryTools.getBestAlignmentForLabel((IAtomContainer) molecule, (IAtom) atom);
                LabelManager.AnnotationPosition suggestedPosition = this.labelManager.alignmentToAnnotationPosition(align);
                if (atom.getSymbol().equals("O") && molecule.getConnectedAtomsCount(atom) == 0) {
                    suggestedPosition = LabelManager.AnnotationPosition.W;
                }
                if (this.labelManager.isUsed(atom, suggestedPosition)) {
                    suggestedPosition = this.labelManager.getNextSparePosition(atom);
                }
                this.labelManager.setUsedPosition(atom, suggestedPosition);
                Rectangle2D hBounds = this.drawImplicitHydrogens(atom, implicitHydrogenCount, suggestedPosition, g);
                if (hBounds != null) {
                    symbolBounds.add(hBounds);
                }
            }
        } else if (this.params.drawRS && this.chiralMap.containsKey(atom)) {
            symbolBounds = this.drawChiralSymbol(atom, this.chiralMap.get(atom), g);
        } else {
            Point2d p = atom.getPoint2d();
            symbolBounds = new Rectangle2D.Double(p.x, p.y, 0.0, 0.0);
        }
        if (this.params.drawAtomID && (idBounds = this.drawAtomID(atom, molecule, g)) != null) {
            symbolBounds.add(idBounds);
        }
        if (this.params.drawLonePairs) {
            Stroke stroke = g.getStroke();
            g.setStroke(new BasicStroke(0.05f));
            this.drawElectronPairs(atom, molecule, lonePairCount, g);
            g.setStroke(stroke);
        }
        return symbolBounds;
    }

    private Rectangle2D drawChiralSymbol(IAtom atom, IStereoAndConformation chirality, Graphics2D g) {
        Point2d p = atom.getPoint2d();
        if (chirality == IStereoAndConformation.NONE) {
            return new Rectangle2D.Double(p.x, p.y, 0.0, 0.0);
        }
        String text = chirality == IStereoAndConformation.R ? "(R)" : "(S)";
        g.setFont(this.chiralSymbolFont);
        Color color = Color.DARK_GRAY;
        return this.drawText(text, p, color, g);
    }

    public Rectangle2D drawImplicitHydrogens(IAtom atom, int implicitHydrogenCount, LabelManager.AnnotationPosition pos, Graphics2D g) {
        String text = atom.getSymbol();
        Point2d p = atom.getPoint2d();
        g.setFont(this.atomSymbolFont);
        Point2f pc = this.getTextPoint(g, text, p.x, p.y);
        Rectangle2D hBounds = this.getTextBounds(g, "H");
        double atomSymbolWidth = this.getTextBounds(g, text).getWidth();
        double hWidth = hBounds.getWidth();
        double hHeight = hBounds.getHeight();
        double subscriptWidth = 0.0;
        Rectangle2D.Double totalHBounds = null;
        g.setColor(Color.BLACK);
        if (pos == LabelManager.AnnotationPosition.E) {
            double cx = p.x + atomSymbolWidth / 2.0 + hWidth / 2.0;
            double cy = p.y;
            Point2f hP = this.getTextPoint(g, "H", cx, cy);
            String hString = "H";
            g.drawString(hString, hP.x, hP.y);
            totalHBounds = new Rectangle2D.Double(cx - hWidth / 2.0, cy - hHeight / 2.0, hWidth, hHeight);
            if (implicitHydrogenCount > 1) {
                g.setFont(this.subscriptFont);
                String hCount = String.valueOf(implicitHydrogenCount);
                Rectangle2D subscriptBounds = this.getTextBounds(g, hCount);
                subscriptWidth = subscriptBounds.getWidth();
                Point2f sP = this.getTextPoint(g, hCount, cx += hWidth / 2.0 + subscriptWidth / 2.0, cy += (double) this.params.subscriptHeight);
                double subscriptHeight = subscriptBounds.getHeight();
                Rectangle2D.Double finalHBounds = new Rectangle2D.Double(cx - subscriptWidth / 2.0, cy - subscriptHeight / 2.0, subscriptWidth, subscriptHeight);
                g.setColor(Color.WHITE);
                g.fill(finalHBounds);
                g.setColor(Color.BLACK);
                g.drawString(hCount, sP.x, sP.y);
                g.setFont(this.atomSymbolFont);
                totalHBounds.add(finalHBounds);
            }
        } else if (pos == LabelManager.AnnotationPosition.W) {
            float y;
            float x;
            if (implicitHydrogenCount > 1) {
                String hCount = String.valueOf(implicitHydrogenCount);
                g.setFont(this.subscriptFont);
                Rectangle2D subscriptBounds = this.getTextBounds(g, hCount);
                subscriptWidth = subscriptBounds.getWidth();
                x = (float) ((double) pc.x - subscriptWidth);
                y = pc.y + (float) this.params.subscriptHeight;
                g.drawString(hCount, x, y);
                g.setFont(this.atomSymbolFont);
                double subscriptHeight = subscriptBounds.getHeight();
                totalHBounds = new Rectangle2D.Double((double) x - subscriptWidth / 2.0, (double) y - subscriptHeight / 2.0, subscriptWidth, subscriptHeight);
            }
            x = (float) ((double) pc.x - atomSymbolWidth / 2.0 - subscriptWidth - hWidth / 2.0);
            y = pc.y;
            String hString = "H";
            Rectangle2D.Double hDrawnBounds = new Rectangle2D.Double(p.x - atomSymbolWidth / 2.0 - subscriptWidth - hWidth, p.y - hBounds.getHeight() / 2.0, hWidth, hHeight);
            g.setColor(Color.WHITE);
            g.fill(hDrawnBounds);
            g.setColor(Color.BLACK);
            g.drawString(hString, x, y);
            if (totalHBounds == null) {
                totalHBounds = hDrawnBounds;
            } else {
                totalHBounds.add(hDrawnBounds);
            }
        }
        return totalHBounds;
    }

    public Rectangle2D drawAtomSymbol(IAtom atom, Graphics2D g) {
        String text = atom.getSymbol();
        if (atom instanceof PseudoAtom) {
            text = ((PseudoAtom) atom).getLabel();
        }
        g.setFont(this.atomSymbolFont);
        Point2d p = atom.getPoint2d();
        return this.drawText(text, p, this.colorForAtom(atom), g);
    }

    private Rectangle2D drawText(String text, Point2d p, Color color, Graphics2D g) {
        Point2f pc = this.getTextPoint(g, text, p.x, p.y);
        Rectangle2D stringBounds = this.getTextBounds(g, text);
        double sW2 = stringBounds.getWidth() / 2.0;
        double sH2 = stringBounds.getHeight() / 2.0;
        double x = p.x - sW2;
        double y = p.y - sH2;
        g.setColor(Color.WHITE);
        Rectangle2D.Double bounds = new Rectangle2D.Double(x, y, sW2 * 2.0, sH2 * 2.0);
        g.fill(bounds);
        g.setColor(color);
        g.drawString(text, pc.x, pc.y);
        return bounds;
    }

    public Rectangle2D drawAtomID(IAtom atom, IAtomContainer container, Graphics2D g) {
        String atomID = atom.getID();
        if (atomID == null) {
            return null;
        }
        g.setFont(this.atomSymbolFont);
        Point2d p = atom.getPoint2d();
        Rectangle2D atomSymbolBounds = this.shouldDraw(atom, container) ? this.getTextBounds(g, atom.getSymbol()) : new Rectangle2D.Double(p.x, p.y, 1.0, 1.0);
        g.setFont(this.atomIDFont);
        Rectangle2D bounds = this.getTextBounds(g, atomID);
        Point2d pID = new Point2d(p);
        LabelManager.AnnotationPosition suggestedPosition = this.labelManager.alignmentToAnnotationPosition(GeometryTools.getBestAlignmentForLabelXY((IAtomContainer) container, (IAtom) atom));
        LabelManager.AnnotationPosition pos = this.labelManager.isUsed(atom, suggestedPosition) ? this.labelManager.getNextSparePosition(atom) : suggestedPosition;
        double aW2 = atomSymbolBounds.getWidth() / 2.0;
        double bW2 = bounds.getWidth() / 2.0;
        double aH2 = atomSymbolBounds.getHeight() / 2.0;
        double bH2 = bounds.getHeight() / 2.0;
        if (null != pos) {
            switch (pos) {
                case N:
                    pID.y -= aH2 + bH2;
                    break;
                case NE:
                    pID.x += aW2 + bW2;
                    pID.y -= aH2 + bH2;
                    break;
                case E:
                    pID.x += aW2 + bW2;
                    break;
                case SE:
                    pID.x += aW2 + bW2;
                    pID.y += aH2 + bH2;
                    break;
                case S:
                    pID.y += aH2 + bH2;
                    break;
                case SW:
                    pID.x -= aW2 + bW2;
                    pID.y += aH2 + bH2;
                    break;
                case W:
                    pID.x -= aW2 + bW2;
                    break;
                case NW:
                    pID.x -= aW2 + bW2;
                    pID.y -= aH2 + bH2;
                    break;
                default:
                    pID.x += aW2 + bW2;
                    break;
            }
        }
        if (pos != null) {
            this.labelManager.setUsedPosition(atom, pos);
        }
        Point2f tp = this.getTextPoint(g, atomID, pID.x, pID.y);
        g.setColor(Color.BLACK);
        g.drawString(atomID, tp.x, tp.y);
        g.setFont(this.atomSymbolFont);
        return new Rectangle2D.Double(pID.x - bounds.getWidth() / 2.0, pID.y - bounds.getHeight() / 2.0, bounds.getWidth(), bounds.getHeight());
    }

    public Rectangle2D drawElectronPairs(IAtom atom, IAtomContainer container, int lonePairCount, Graphics2D g) {
        if (lonePairCount == 0) {
            return null;
        }
        Point2d atomPoint = atom.getPoint2d();
        Rectangle2D atomSymbolBounds = this.getTextBounds(g, atom.getSymbol());
        BitSet positions = this.labelManager.getAtomAnnotationPositions(atom);
        double r = this.params.electronRadius;
        double d = r * 2.0;
        for (int i = 0; i < lonePairCount; ++i) {
            LabelManager.AnnotationPosition position = this.labelManager.getNextSparePosition(positions);
            Vector2d v = this.labelManager.getVectorFromPosition(position);
            Vector2d leftPerp = this.labelManager.getLeftPerpendicularFromPosition(position);
            Vector2d rightPerp = this.labelManager.getRightPerpendicularFromPosition(position);
            double dx = (atomSymbolBounds.getWidth() / 2.0 + d) * v.x;
            double dy = (atomSymbolBounds.getHeight() / 2.0 + d) * v.y;
            Point2d lp = new Point2d(atomPoint.x + dx, atomPoint.y + dy);
            Point2d llp = new Point2d(lp);
            llp.scaleAdd((double) (this.params.lonePairSeparation / 2), (Tuple2d) leftPerp, (Tuple2d) llp);
            Point2d rlp = new Point2d(lp);
            rlp.scaleAdd((double) (this.params.lonePairSeparation / 2), (Tuple2d) rightPerp, (Tuple2d) rlp);
            g.fill(new Ellipse2D.Double(llp.x - r, llp.y - r, d, d));
            g.fill(new Ellipse2D.Double(rlp.x - r, rlp.y - r, d, d));
            positions.set(position.ordinal());
        }
        return null;
    }

    private boolean shouldDraw(IAtom atom, IAtomContainer atomContainer) {
        String symbol = atom.getSymbol();
        if (symbol.equals("C")) {
            if (this.params.drawCarbons) {
                return true;
            }
            if (this.params.drawTerminalCarbons && this.isTerminal(atom, atomContainer)) {
                return true;
            }
            return this.getAttachedMultipleBondCount(atom, atomContainer) > 1;
        }
        if (symbol.equals("H")) {
            return this.params.drawExplicitHydrogens;
        }
        return true;
    }

    private int getAttachedMultipleBondCount(IAtom atom, IAtomContainer atomContainer) {
        int count = 0;
        for (IBond bond : atomContainer.getConnectedBondsList(atom)) {
            if (bond.getOrder() == IBond.Order.SINGLE) {
                continue;
            }
            ++count;
        }
        return count;
    }

    public boolean isCharged(IAtom atom) {
        Integer formalCharge = atom.getFormalCharge();
        return formalCharge != null && formalCharge != 0;
    }

    private boolean isTerminal(IAtom atom, IAtomContainer atomContainer) {
        int numberOfHeavyAtomsConnected = 0;
        for (IAtom connected : atomContainer.getConnectedAtomsList(atom)) {
            if (connected.getSymbol().equals("H")) {
                continue;
            }
            ++numberOfHeavyAtomsConnected;
        }
        return numberOfHeavyAtomsConnected < 2;
    }

    private Rectangle2D drawCharge(IAtom atom, Graphics2D g) {
        BitSet annotationPositions = this.labelManager.getAtomAnnotationPositions(atom);
        Integer formalCharge = atom.getFormalCharge();
        String chargeText = this.getChargeString(formalCharge);
        Rectangle2D atomBounds = this.getTextBounds(g, atom.getSymbol());
        Rectangle2D chargeBounds = this.getTextBounds(g, chargeText);
        g.setColor(Color.BLACK);
        Point2d atomPoint = atom.getPoint2d();
        Point2d chargePoint = new Point2d(atomPoint);
        double chargeDim = Math.min(chargeBounds.getWidth(), chargeBounds.getHeight());
        chargePoint.x += atomBounds.getWidth() / 2.0 + chargeDim / 2.0;
        chargePoint.y -= atomBounds.getHeight() / 2.0;
        annotationPositions.set(LabelManager.AnnotationPosition.NE.ordinal());
        Point2f sp = this.getTextPoint(g, chargeText, chargePoint.x, chargePoint.y);
        Rectangle2D.Double chargeBox = new Rectangle2D.Double(chargePoint.x - chargeBounds.getWidth() / 2.0, chargePoint.y - chargeBounds.getHeight() / 2.0, chargeBounds.getWidth(), chargeBounds.getHeight());
        g.setColor(Color.WHITE);
        g.fill(chargeBox);
        g.setColor(Color.BLACK);
        g.drawString(chargeText, sp.x, sp.y);
        return chargeBox;
    }

    private String getChargeString(Integer formalCharge) {
        if (formalCharge == 1) {
            return "+";
        }
        if (formalCharge == -1) {
            return "-";
        }
        if (formalCharge > 1) {
            return formalCharge + "+";
        }
        if (formalCharge < -1) {
            return formalCharge + "-";
        }
        return "";
    }

    private Map<IAtom, Integer> getLonePairCounts(IAtomContainer atomContainer) {
        HashMap<IAtom, Integer> lonePairMap = new HashMap<>();
        for (ILonePair lonePair : atomContainer.lonePairs()) {
            IAtom atom = lonePair.getAtom();
            int lonePairCount = lonePairMap.containsKey(atom) ? lonePairMap.get(atom) : 0;
            lonePairMap.put(atom, lonePairCount + 1);
        }
        return lonePairMap;
    }

    public Color colorForAtom(IAtom atom) {
        return this.atomColorer.getAtomColor(atom);
    }
}
