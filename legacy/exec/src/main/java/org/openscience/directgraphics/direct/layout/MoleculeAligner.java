/*
 * SPDX-License-Identifier: Apache-2.0
 * © 2025 BioInception PVT LTD.
 */
package org.openscience.directgraphics.direct.layout;

import javax.vecmath.Point2d;
import javax.vecmath.Tuple2d;
import javax.vecmath.Vector2d;
import org.openscience.cdk.geometry.GeometryTools;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.directgraphics.direct.ConvexHull;

public class MoleculeAligner {
    public static final Vector2d X_AXIS = new Vector2d(1.0, 0.0);
    public static final Vector2d Y_AXIS = new Vector2d(0.0, 1.0);

    public static void alignToMinAreaBox(IAtomContainer atomContainer, Vector2d axis) {
        ConvexHull hull = new ConvexHull(atomContainer);
        MoleculeAligner.alignToAxis(atomContainer, hull.getMajorAxis(), axis, hull.getCenter());
    }

    public static Vector2d getMaxWidthVector(IAtomContainer atomContainer) {
        int nAtoms = atomContainer.getAtomCount();
        Vector2d widthVector = null;
        IAtom maxI = null;
        IAtom maxJ = null;
        double maxDistance = 0.0;
        for (int indexI = nAtoms - 1; indexI >= 0; --indexI) {
            IAtom atomI = atomContainer.getAtom(indexI);
            Point2d pointI = atomI.getPoint2d();
            if (pointI == null) continue;
            for (int indexJ = indexI - 1; indexJ >= 0; --indexJ) {
                double distance;
                IAtom atomJ = atomContainer.getAtom(indexJ);
                Point2d pointJ = atomJ.getPoint2d();
                if (pointJ == null || (distance = pointI.distance(pointJ)) <= maxDistance) continue;
                maxDistance = distance;
                maxI = atomI;
                maxJ = atomJ;
            }
        }
        if (maxI == null || maxJ == null) {
            return new Vector2d(0.0, 0.0);
        }
        widthVector = new Vector2d((Tuple2d)maxI.getPoint2d());
        widthVector.sub((Tuple2d)maxJ.getPoint2d());
        return widthVector;
    }

    public static void alignToMaxWidth(IAtomContainer atomContainer, Vector2d axis) {
        Vector2d widthVector = MoleculeAligner.getMaxWidthVector(atomContainer);
        Point2d center = GeometryTools.get2DCenter((IAtomContainer)atomContainer);
        MoleculeAligner.alignToAxis(atomContainer, widthVector, axis, center);
    }

    private static double getPolarAngle(Vector2d vector) {
        double x = vector.x;
        double y = vector.y;
        if (x > 0.0) {
            return Math.atan(y / x);
        }
        if (x < 0.0) {
            if (y >= 0.0) {
                return Math.atan(y / x) + 3.141592653589793;
            }
            return Math.atan(y / x) - 3.141592653589793;
        }
        if (y > 0.0) {
            return 1.5707963267948966;
        }
        if (y < 0.0) {
            return -1.5707963267948966;
        }
        return 0.0;
    }

    public static double getMinAngle(Vector2d axisFrom, Vector2d axisTo) {
        double polarAngleForwardFrom = Math.atan2(axisFrom.y, axisFrom.x);
        double polarAngleBackwardFrom = Math.atan2(- axisFrom.y, - axisFrom.x);
        double polarAngleTo = Math.atan2(axisTo.y, axisTo.x);
        double forwardDiff = polarAngleForwardFrom - polarAngleTo;
        double backwardDiff = polarAngleBackwardFrom - polarAngleTo;
        double minAngleDiff = Math.abs(forwardDiff) < Math.abs(backwardDiff) ? forwardDiff : backwardDiff;
        return - minAngleDiff;
    }

    public static void alignToAxis(IAtomContainer atomContainer, Vector2d axisFrom, Vector2d axisTo, Point2d center) {
        double angle = MoleculeAligner.getMinAngle(axisFrom, axisTo);
        double cosA = Math.cos(angle);
        double sinA = Math.sin(angle);
        double minCosA = 1.0 - cosA;
        for (IAtom atom : atomContainer.atoms()) {
            Point2d p = atom.getPoint2d();
            double x = cosA * p.x - sinA * p.y + center.x * minCosA + center.y * sinA;
            double y = sinA * p.x + cosA * p.y + center.y * minCosA - center.x * sinA;
            p.x = x;
            p.y = y;
            atom.setPoint2d(p);
        }
    }
}

