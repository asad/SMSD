/*
 * SPDX-License-Identifier: Apache-2.0
 * © 2025 BioInception PVT LTD.
 */
package org.openscience.directgraphics.direct;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.Stroke;
import java.awt.geom.Ellipse2D;
import java.awt.geom.Path2D;
import java.awt.geom.Rectangle2D;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import javax.vecmath.Point2d;
import javax.vecmath.Tuple2d;
import javax.vecmath.Vector2d;
import org.openscience.cdk.geometry.GeometryTools;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IRing;
import org.openscience.cdk.interfaces.IRingSet;
import org.openscience.cdk.ringsearch.SSSRFinder;
import org.openscience.cdk.tools.manipulator.AtomContainerComparatorBy2DCenter;

public class DirectBondDrawer
        extends AbstractDirectDrawer {

    private final LabelManager labelManager;
    private final Stroke dashedWedgeStroke;
    private Stroke bondStroke;

    public DirectBondDrawer(Params params, LabelManager labelManager) {
        super.params = params;
        this.labelManager = labelManager;
        this.dashedWedgeStroke = new BasicStroke(params.dashedWedgeStroke);
    }

    private void setBondStroke() {
        int cap = 0;
        int join = 2;
        if (this.params.bondStrokeCap == Params.BondStrokeCap.ROUND) {
            cap = 1;
        } else if (this.params.bondStrokeCap == Params.BondStrokeCap.SQUARE) {
            cap = 2;
        }
        if (this.params.bondStrokeJoin == Params.BondStrokeJoin.BEVEL) {
            join = 2;
        } else if (this.params.bondStrokeJoin == Params.BondStrokeJoin.ROUND) {
            join = 1;
        }
        this.bondStroke = new BasicStroke(this.params.bondStrokeWidth, cap, join);
    }

    public void drawBonds(IAtomContainer molecule, Graphics2D g) {
        this.setBondStroke();
        g.setStroke(this.bondStroke);
        IRingSet ringSet = new SSSRFinder(molecule).findSSSR();
        ringSet.sortAtomContainers((Comparator) new AtomContainerComparatorBy2DCenter());
        this.addRingCentersToAtomAnnotationPositions(molecule, ringSet);
        Map<IBond, IAtomContainer> bondRingMap = this.fillBondRingMap(ringSet);
        g.setColor(Color.BLACK);
        for (IBond bond : molecule.bonds()) {
            if (this.shouldDraw(bond)) {
                this.drawBond(bond, bondRingMap, g);
            }
            this.labelManager.addBondToAtomAnnotationPositions(bond);
        }
        if (this.params.drawAromaticCircles) {
            for (IAtomContainer ring : ringSet.atomContainers()) {
                if (!this.ringIsAromatic(ring)) {
                    continue;
                }
                this.drawRingCircle(ring, g);
            }
        }
        ArrayList<IBond> drawnRingBonds = new ArrayList<>();
        for (IAtomContainer ring : ringSet.atomContainers()) {
            Point2d c = GeometryTools.get2DCenter((IAtomContainer) ring);
            for (IBond bond2 : ring.bonds()) {
                if (drawnRingBonds.contains(bond2)
                        || bond2.getFlag(32)
                        && this.params.drawAromaticCircles
                        || bond2.getOrder() == IBond.Order.SINGLE) {
                    continue;
                }
                Point2d p1 = bond2.getAtom(0).getPoint2d();
                Point2d p2 = bond2.getAtom(1).getPoint2d();
                this.drawOffsetBond(p1, p2, c, g);
                drawnRingBonds.add(bond2);
            }
        }
    }

    private void addRingCentersToAtomAnnotationPositions(IAtomContainer mol, IRingSet ringSet) {
        for (IAtomContainer ring : ringSet.atomContainers()) {
            for (IAtom atom : ring.atoms()) {
                List<IAtom> connectedAtoms = mol.getConnectedAtomsList(atom);
                ArrayList<IAtom> connectedAtomsInRing = new ArrayList<>();
                for (IAtom connectedAtom : connectedAtoms) {
                    if (!ring.contains(connectedAtom)) {
                        continue;
                    }
                    connectedAtomsInRing.add(connectedAtom);
                }
                this.labelManager.addRingCenterToAtomAnnotationPosition(atom, connectedAtomsInRing);
            }
        }
    }

    private Map<IBond, IAtomContainer> fillBondRingMap(IRingSet ringSet) {
        HashMap<IBond, IAtomContainer> bondRingMap = new HashMap<>();
        for (IAtomContainer ringAsAtomContainer : ringSet.atomContainers()) {
            IRing ring = (IRing) ringAsAtomContainer;
            for (IBond bond : ring.bonds()) {
                bondRingMap.put(bond, (IAtomContainer) ring);
            }
        }
        return bondRingMap;
    }

    public void drawBond(IBond bond, Map<IBond, IAtomContainer> bondRingMap, Graphics2D g) {
        Point2d p1 = bond.getAtom(0).getPoint2d();
        Point2d p2 = bond.getAtom(1).getPoint2d();
        IBond.Order order = bond.getOrder();
        IBond.Stereo stereo = bond.getStereo();
        if (stereo == IBond.Stereo.NONE && (order == IBond.Order.SINGLE || bond.getFlag(32))) {
            this.drawLine(p1, p2, g);
        } else if (order == IBond.Order.DOUBLE) {
            if (bondRingMap.containsKey(bond)) {
                this.drawLine(p1, p2, g);
            } else {
                this.drawDoubleBond(p1, p2, g);
            }
        } else if (order == IBond.Order.TRIPLE) {
            this.drawTripleBond(p1, p2, g);
        } else if (stereo != IBond.Stereo.NONE) {
            this.drawStereo(p1, p2, stereo, g);
        }
    }

    private void drawTripleBond(Point2d p1, Point2d p2, Graphics2D g) {
        Vector2d perpendicular = this.makePerpendicular(p1, p2);
        perpendicular.scale(this.params.tripleBondGap);
        Vector2d negativePerp = new Vector2d(perpendicular);
        negativePerp.negate();
        this.drawLine(this.displace(p1, perpendicular), this.displace(p2, perpendicular), g);
        this.drawLine(p1, p2, g);
        this.drawLine(this.displace(p1, negativePerp), this.displace(p2, negativePerp), g);
    }

    private void drawStereo(Point2d p1, Point2d p2, IBond.Stereo stereo, Graphics2D g) {
        if (null != stereo) {
            switch (stereo) {
                case UP_OR_DOWN:
                    this.drawWigglyLine(p1, p2, g);
                    break;
                case DOWN:
                    this.drawWedge(p1, p2, false, g);
                    break;
                case DOWN_INVERTED:
                    this.drawWedge(p2, p1, false, g);
                    break;
                case UP:
                    this.drawWedge(p1, p2, true, g);
                    break;
                case UP_INVERTED:
                    this.drawWedge(p2, p1, true, g);
                    break;
                default:
                    break;
            }
        }
    }

    private void drawWedge(Point2d p1, Point2d p2, boolean isFilled, Graphics2D g) {
        Vector2d halfWidthVector = new Vector2d(p2.y - p1.y, p1.x - p2.x);
        halfWidthVector.normalize();
        halfWidthVector.scale((double) (this.params.filledWedgeWidth / 2));
        Vector2d negHalfWidthVector = new Vector2d(halfWidthVector);
        negHalfWidthVector.negate();
        Point2d p2a = this.displace(p2, halfWidthVector);
        Point2d p2b = this.displace(p2, negHalfWidthVector);
        if (isFilled) {
            this.drawFilledWedge(p1, p2a, p2b, g);
        } else {
            this.drawDashedWedge2(p1, p2a, p2b, g);
        }
    }

    public void drawDashedWedge(Point2d a, Point2d b, Point2d c, Graphics2D g) {
        Stroke savedStroke = g.getStroke();
        g.setStroke(this.dashedWedgeStroke);
        double distance = b.distance(a);
        double gapFactor = this.params.dashedGapFactor;
        double gap = distance * gapFactor;
        double numberOfDashes = distance / gap;
        double d = 0.0;
        int i = 0;
        while ((double) i < numberOfDashes) {
            Point2d p1 = new Point2d();
            p1.interpolate((Tuple2d) a, (Tuple2d) b, d);
            Point2d p2 = new Point2d();
            p2.interpolate((Tuple2d) a, (Tuple2d) c, d);
            this.drawLine(p1, p2, g);
            if (distance * (d + gapFactor) >= distance) {
                break;
            }
            d += gapFactor;
            ++i;
        }
        g.setStroke(savedStroke);
    }

    public void drawDashedWedge2(Point2d a, Point2d b, Point2d c, Graphics2D g) {
        Stroke savedStroke = g.getStroke();
        g.setStroke(this.dashedWedgeStroke);
        double distance = b.distance(a);
        double gapFactor = this.params.dashedGapFactor;
        double gap = distance * gapFactor;
        double numberOfDashes = distance / gap;
        double currentDistance = 0.0;
        Point2d d = new Point2d(b);
        d.interpolate((Tuple2d) c, 0.5);
        Vector2d perp = this.makePerpendicular(a, d);
        Vector2d nPerp = new Vector2d(perp);
        nPerp.negate();
        double maxWidth = this.params.dashedWedgeWidth / 4.0;
        double currentWidth = maxWidth * this.params.dashedWidthFactor;
        int i = 0;
        while ((double) i < numberOfDashes) {
            Point2d rungCenter = new Point2d(a);
            rungCenter.interpolate((Tuple2d) d, currentDistance);
            Point2d p1 = new Point2d(rungCenter);
            p1.scaleAdd(currentWidth, (Tuple2d) perp, (Tuple2d) p1);
            Point2d p2 = new Point2d(rungCenter);
            p2.scaleAdd(currentWidth, (Tuple2d) nPerp, (Tuple2d) p2);
            this.drawLine(p1, p2, g);
            if (distance * (currentDistance + gapFactor) >= distance) {
                break;
            }
            currentDistance += gapFactor;
            currentWidth += maxWidth * this.params.dashedWidthFactor;
            ++i;
        }
        g.setStroke(savedStroke);
    }

    private void drawFilledWedge(Point2d a, Point2d b, Point2d c, Graphics2D g) {
        Path2D.Double path = new Path2D.Double();
        path.moveTo(a.x, a.y);
        path.lineTo(b.x, b.y);
        path.lineTo(c.x, c.y);
        path.closePath();
        g.fill(path);
    }

    public void drawWigglyLine(Point2d p1, Point2d p2, Graphics2D g) {
        double gapProportion = 0.1;
        double wiggleWidth = this.params.wiggleLineWidth;
        Vector2d line = new Vector2d((Tuple2d) p2);
        line.sub((Tuple2d) p1);
        double length = line.length();
        double gap = length * gapProportion;
        int numberOfSegments = 10;
        line.normalize();
        Vector2d perpendicular = this.makePerpendicular(line);
        Vector2d negPerp = new Vector2d(perpendicular);
        negPerp.negate();
        Point2d centerLinePoint = new Point2d(p1);
        Path2D.Double path = new Path2D.Double();
        path.moveTo(p1.x, p1.y);
        centerLinePoint.scaleAdd(gap / 2.0, (Tuple2d) line, (Tuple2d) centerLinePoint);
        Point2d tipPoint = new Point2d(centerLinePoint);
        tipPoint.scaleAdd(wiggleWidth / 2.0, (Tuple2d) perpendicular, (Tuple2d) tipPoint);
        for (int i = 0; i < numberOfSegments - 1; ++i) {
            centerLinePoint.scaleAdd(gap / 2.0, (Tuple2d) line, (Tuple2d) centerLinePoint);
            path.quadTo(tipPoint.x, tipPoint.y, centerLinePoint.x, centerLinePoint.y);
            centerLinePoint.scaleAdd(gap / 2.0, (Tuple2d) line, (Tuple2d) centerLinePoint);
            tipPoint = new Point2d(centerLinePoint);
            if (i % 2 == 0) {
                tipPoint.scaleAdd(wiggleWidth / 2.0, (Tuple2d) negPerp, (Tuple2d) tipPoint);
                continue;
            }
            tipPoint.scaleAdd(wiggleWidth / 2.0, (Tuple2d) perpendicular, (Tuple2d) tipPoint);
        }
        g.draw(path);
    }

    private Vector2d makePerpendicular(Point2d p1, Point2d p2) {
        Vector2d line = new Vector2d((Tuple2d) p1);
        line.sub((Tuple2d) p2);
        line.normalize();
        return this.makePerpendicular(line);
    }

    private void drawDoubleBond(Point2d p1, Point2d p2, Graphics2D g) {
        Vector2d perpendicular = this.makePerpendicular(p1, p2);
        perpendicular.scale(this.params.doubleBondGap);
        Vector2d negativePerp = new Vector2d(perpendicular);
        negativePerp.negate();
        this.drawLine(this.displace(p1, perpendicular), this.displace(p2, perpendicular), g);
        this.drawLine(this.displace(p1, negativePerp), this.displace(p2, negativePerp), g);
    }

    private void drawOffsetBond(Point2d p1, Point2d p2, Point2d c, Graphics2D g) {
        double distanceProportion = this.params.offsetBondDistanceProportion;
        Point2d w = new Point2d();
        w.interpolate((Tuple2d) c, (Tuple2d) p1, distanceProportion);
        Point2d u = new Point2d();
        u.interpolate((Tuple2d) c, (Tuple2d) p2, distanceProportion);
        this.drawLine(w, u, g);
    }

    private Vector2d makePerpendicular(Vector2d line) {
        Vector2d perp = new Vector2d(-line.y, line.x);
        perp.normalize();
        return perp;
    }

    private Point2d displace(Point2d point, Vector2d vector) {
        Point2d displacedPoint = new Point2d(point);
        displacedPoint.add((Tuple2d) vector);
        return displacedPoint;
    }

    private void drawRingCircle(IAtomContainer ring, Graphics2D g) {
        Point2d center = GeometryTools.get2DCenter((IAtomContainer) ring);
        Rectangle2D bounds = GeometryTools.getRectangle2D((IAtomContainer) ring);
        double diameter = Math.min(bounds.getWidth(), bounds.getHeight());
        double radius = (diameter *= this.params.ringProportion) / 2.0;
        g.draw(new Ellipse2D.Double(center.x - radius, center.y - radius, diameter, diameter));
    }

    private boolean ringIsAromatic(IAtomContainer ring) {
        for (IAtom atom : ring.atoms()) {
            if (atom.getFlag(32)) {
                continue;
            }
            return false;
        }
        for (IBond b : ring.bonds()) {
            if (b.getFlag(32)) {
                continue;
            }
            return false;
        }
        return true;
    }

    private boolean shouldDraw(IBond bond) {
        boolean neitherAreH;
        boolean symbol0IsH = bond.getAtom(0).getSymbol().equals("H");
        boolean symbol1IsH = bond.getAtom(1).getSymbol().equals("H");
        boolean bothAreH = symbol0IsH && symbol1IsH;
        boolean atLeastOneIsH = symbol0IsH || symbol1IsH;
        boolean bl = neitherAreH = !symbol0IsH && !symbol1IsH;
        if (bothAreH || neitherAreH) {
            return true;
        }
        if (atLeastOneIsH) {
            return this.params.drawExplicitHydrogens;
        }
        return true;
    }
}
