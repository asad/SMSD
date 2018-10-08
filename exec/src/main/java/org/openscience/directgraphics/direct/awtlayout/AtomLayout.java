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

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Font;
import java.awt.Graphics2D;
import java.awt.Stroke;
import java.awt.geom.Ellipse2D;
import java.awt.geom.Point2D;
import java.awt.geom.Rectangle2D;
import java.util.BitSet;
import java.util.HashMap;
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
import org.openscience.directgraphics.direct.LabelManager;
import org.openscience.directgraphics.direct.Params;
import org.openscience.directgraphics.direct.layout.BoundsTree;
import org.openscience.directgraphics.stereo.IStereoAndConformation;

public class AtomLayout
        extends AbstractAWTLayout<IAtom> {

    private Font atomSymbolFont;
    private Font subscriptFont;
    private Font atomIDFont;
    private Font chiralSymbolFont;
    private final IAtomColorer atomColorer;
    private final LabelManager labelManager;
    private Map<IAtom, IStereoAndConformation> chiralMap;

    public AtomLayout(AbstractAWTLayout parent, Params params, LabelManager labelManager) {
        this.parent = parent;
        this.setParams(params);
        this.labelManager = labelManager;
        this.atomColorer = new CDK2DAtomColors();
        this.chiralMap = new HashMap<>();
    }

    @Override
    public BoundsTree layout(IAtom atom, Graphics2D g) {
        Rectangle2D idBounds;
        this.currentObject = atom;
        String id = atom.getID();
        IAtomContainer molecule = null;
        if (this.parent != null) {
            molecule = (IAtomContainer) this.parent.getCurrentObject();
        }
        this.boundsTree = new BoundsTree(atom.getID());
        if (molecule == null || this.shouldDraw(atom, molecule)) {
            Integer implicitHydrogenCount;
            this.boundsTree.add(id + ":symbol", this.layoutAtomSymbol(atom, g));
            if (this.isCharged(atom)) {
                Rectangle2D chargeBounds = this.layoutCharge(atom, g);
                this.boundsTree.add(id + ":charge", chargeBounds);
            }
            if (this.params.drawImplicitHydrogens && (implicitHydrogenCount = atom.getImplicitHydrogenCount()) != null && implicitHydrogenCount > 0) {
                int align = 1;
                if (molecule != null) {
                    GeometryTools.getBestAlignmentForLabel((IAtomContainer) molecule, (IAtom) atom);
                }
                LabelManager.AnnotationPosition suggestedPosition = this.labelManager.alignmentToAnnotationPosition(align);
                if (atom.getSymbol().equals("O") && (molecule == null || molecule.getConnectedAtomsCount(atom) == 0)) {
                    suggestedPosition = LabelManager.AnnotationPosition.W;
                }
                if (this.labelManager.isUsed(atom, suggestedPosition)) {
                    suggestedPosition = this.labelManager.getNextSparePosition(atom);
                }
                this.labelManager.setUsedPosition(atom, suggestedPosition);
                Rectangle2D hBounds = this.layoutImplicitHydrogens(atom, implicitHydrogenCount, suggestedPosition, g);
                if (hBounds != null) {
                    this.boundsTree.add(id + ":hs", hBounds);
                }
            }
        } else if (this.params.drawRS && this.chiralMap.containsKey(atom)) {
            this.boundsTree.add(id + ":chiral", this.layoutChiralSymbol(atom, this.chiralMap.get(atom), g));
        } else {
            Point2d p = atom.getPoint2d();
            this.boundsTree.add(id + ":symbol", new Point2D.Double(p.x, p.y));
        }
        if (this.params.drawAtomID && molecule != null && (idBounds = this.layoutAtomID(atom, molecule, g)) != null) {
            this.boundsTree.add(id + ":id", idBounds);
        }
        if (this.params.drawLonePairs && molecule != null) {
            int lonePairCount = 0;
            for (ILonePair lonePair : molecule.lonePairs()) {
                if (!lonePair.contains(atom)) {
                    continue;
                }
                ++lonePairCount;
            }
            if (lonePairCount > 0) {
                Stroke stroke = g.getStroke();
                g.setStroke(new BasicStroke(0.05f));
                this.layoutElectronPairs(atom, molecule, lonePairCount, g);
                g.setStroke(stroke);
            }
        }
        return this.boundsTree;
    }

    public Rectangle2D layoutAtomSymbol(IAtom atom, Graphics2D g) {
        String text = atom.getSymbol();
        if (atom instanceof PseudoAtom) {
            text = ((PseudoAtom) atom).getLabel();
        }
        g.setFont(this.atomSymbolFont);
        Point2d p = atom.getPoint2d();
        return this.layoutText(text, p, g);
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

    private Rectangle2D layoutChiralSymbol(IAtom atom, IStereoAndConformation chirality, Graphics2D g) {
        Point2d p = atom.getPoint2d();
        if (chirality == IStereoAndConformation.NONE) {
            return new Rectangle2D.Double(p.x, p.y, 0.0, 0.0);
        }
        String text = chirality == IStereoAndConformation.R ? "(R)" : "(S)";
        g.setFont(this.chiralSymbolFont);
        return this.layoutText(text, p, g);
    }

    public Rectangle2D layoutImplicitHydrogens(IAtom atom, int implicitHydrogenCount, LabelManager.AnnotationPosition pos, Graphics2D g) {
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
        if (pos == LabelManager.AnnotationPosition.E) {
            double cx = p.x + atomSymbolWidth / 2.0 + hWidth / 2.0;
            double cy = p.y;
            totalHBounds = new Rectangle2D.Double(cx - hWidth / 2.0, cy - hHeight / 2.0, hWidth, hHeight);
            if (implicitHydrogenCount > 1) {
                g.setFont(this.subscriptFont);
                String hCount = String.valueOf(implicitHydrogenCount);
                Rectangle2D subscriptBounds = this.getTextBounds(g, hCount);
                subscriptWidth = subscriptBounds.getWidth();
                g.setFont(this.atomSymbolFont);
                double subscriptHeight = subscriptBounds.getHeight();
                totalHBounds.add(new Rectangle2D.Double((cx += hWidth / 2.0 + subscriptWidth / 2.0) - subscriptWidth / 2.0, (cy += (double) this.params.subscriptHeight) - subscriptHeight / 2.0, subscriptWidth, subscriptHeight));
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
                g.setFont(this.atomSymbolFont);
                double subscriptHeight = subscriptBounds.getHeight();
                totalHBounds = new Rectangle2D.Double((double) x - subscriptWidth / 2.0, (double) y - subscriptHeight / 2.0, subscriptWidth, subscriptHeight);
            }
            x = (float) ((double) pc.x - atomSymbolWidth / 2.0 - subscriptWidth - hWidth / 2.0);
            y = pc.y;
            Rectangle2D.Double hDrawnBounds = new Rectangle2D.Double((double) x - hWidth / 2.0, (double) y - hHeight / 2.0, hWidth, hHeight);
            if (totalHBounds == null) {
                totalHBounds = hDrawnBounds;
            } else {
                totalHBounds.add(hDrawnBounds);
            }
        }
        return totalHBounds;
    }

    private Rectangle2D layoutText(String text, Point2d p, Graphics2D g) {
        Rectangle2D stringBounds = this.getTextBounds(g, text);
        double sW2 = stringBounds.getWidth() / 2.0;
        double sH2 = stringBounds.getHeight() / 2.0;
        double x = p.x - sW2;
        double y = p.y - sH2;
        return new Rectangle2D.Double(x, y, sW2 * 2.0, sH2 * 2.0);
    }

    public Rectangle2D layoutAtomID(IAtom atom, IAtomContainer container, Graphics2D g) {
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
        if (null != pos) switch (pos) {
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
        if (pos != null) {
            this.labelManager.setUsedPosition(atom, pos);
        }
        g.setFont(this.atomSymbolFont);
        return new Rectangle2D.Double(pID.x - bounds.getWidth() / 2.0, pID.y - bounds.getHeight() / 2.0, bounds.getWidth(), bounds.getHeight());
    }

    public Rectangle2D layoutElectronPairs(IAtom atom, IAtomContainer container, int lonePairCount, Graphics2D g) {
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

    private Rectangle2D layoutCharge(IAtom atom, Graphics2D g) {
        BitSet annotationPositions = this.labelManager.getAtomAnnotationPositions(atom);
        Integer formalCharge = atom.getFormalCharge();
        String chargeText = this.getChargeString(formalCharge);
        Rectangle2D atomBounds = this.getTextBounds(g, atom.getSymbol());
        Rectangle2D chargeBounds = this.getTextBounds(g, chargeText);
        Point2d atomPoint = atom.getPoint2d();
        Point2d chargePoint = new Point2d(atomPoint);
        double chargeDim = Math.min(chargeBounds.getWidth(), chargeBounds.getHeight());
        chargePoint.x += atomBounds.getWidth() / 2.0 + chargeDim / 2.0;
        chargePoint.y -= atomBounds.getHeight() / 2.0;
        annotationPositions.set(LabelManager.AnnotationPosition.NE.ordinal());
        return new Rectangle2D.Double(chargePoint.x - chargeBounds.getWidth() / 2.0, chargePoint.y - chargeBounds.getHeight() / 2.0, chargeBounds.getWidth(), chargeBounds.getHeight());
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

    public Color colorForAtom(IAtom atom) {
        return this.atomColorer.getAtomColor(atom);
    }

    @Override
    public BoundsTree layout(IAtom obj, String rootLabel, Graphics2D graphics) {
        return null;
    }

    public void reset() {
        this.labelManager.reset();
    }
}
