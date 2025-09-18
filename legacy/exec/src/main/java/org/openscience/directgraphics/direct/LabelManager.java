/*
 * SPDX-License-Identifier: Apache-2.0
 * © 2025 BioInception PVT LTD.
 */
package org.openscience.directgraphics.direct;

import java.util.ArrayList;
import java.util.BitSet;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import javax.vecmath.Point2d;
import javax.vecmath.Tuple2d;
import javax.vecmath.Vector2d;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IBond;

public class LabelManager {

    private static final Vector2d POS_X = new Vector2d(1.0, 0.0);
    private static final Vector2d POS_Y = new Vector2d(0.0, 1.0);
    private static final Vector2d vN = new Vector2d(0.0, -1.0);
    private static final Vector2d vNE = new Vector2d(1.0, -1.0);
    private static final Vector2d vE = new Vector2d(1.0, 0.0);
    private static final Vector2d vSE = new Vector2d(1.0, 1.0);
    private static final Vector2d vS = new Vector2d(0.0, 1.0);
    private static final Vector2d vSW = new Vector2d(-1.0, 1.0);
    private static final Vector2d vW = new Vector2d(-1.0, 0.0);
    private static final Vector2d vNW = new Vector2d(-1.0, -1.0);
    private final Map<IAtom, BitSet> atomAnnotationPositions = new HashMap<>();

    public String getAnnotationPositionsAsString(IAtom atom) {
        StringBuilder sb = new StringBuilder("|");
        BitSet positions = this.getAtomAnnotationPositions(atom);
        AnnotationPosition[] values = AnnotationPosition.values();
        for (int i = 0; i < values.length; ++i) {
            if (!positions.get(i)) {
                continue;
            }
            sb.append((Object) values[i]);
            sb.append("|");
        }
        return sb.toString();
    }

    public AnnotationPosition getNextSparePosition(IAtom atom) {
        return this.getNextSparePosition(this.getAtomAnnotationPositions(atom));
    }

    public AnnotationPosition getNextSparePosition(BitSet positions) {
        for (int i = 0; i < AnnotationPosition.values().length; ++i) {
            if (positions.get(i)) {
                continue;
            }
            return AnnotationPosition.values()[i];
        }
        return null;
    }

    public Vector2d getVectorFromPosition(AnnotationPosition position) {
        switch (position) {
            case N: {
                return vN;
            }
            case NE: {
                return vNE;
            }
            case E: {
                return vE;
            }
            case SE: {
                return vSE;
            }
            case S: {
                return vS;
            }
            case SW: {
                return vSW;
            }
            case W: {
                return vW;
            }
            case NW: {
                return vNW;
            }
        }
        return vN;
    }

    public Vector2d getLeftPerpendicularFromPosition(AnnotationPosition position) {
        switch (position) {
            case N: {
                return vW;
            }
            case NE: {
                return vNW;
            }
            case E: {
                return vN;
            }
            case SE: {
                return vNE;
            }
            case S: {
                return vE;
            }
            case SW: {
                return vSE;
            }
            case W: {
                return vS;
            }
            case NW: {
                return vSW;
            }
        }
        return vN;
    }

    public Vector2d getRightPerpendicularFromPosition(AnnotationPosition position) {
        switch (position) {
            case N: {
                return vE;
            }
            case NE: {
                return vSE;
            }
            case E: {
                return vS;
            }
            case SE: {
                return vSW;
            }
            case S: {
                return vW;
            }
            case SW: {
                return vNW;
            }
            case W: {
                return vN;
            }
            case NW: {
                return vNE;
            }
        }
        return vS;
    }

    public BitSet getAtomAnnotationPositions(IAtom atom) {
        if (this.atomAnnotationPositions.containsKey(atom)) {
            return this.atomAnnotationPositions.get(atom);
        }
        BitSet positions = new BitSet();
        this.atomAnnotationPositions.put(atom, positions);
        return positions;
    }

    public void setUsedPosition(IAtom atom, AnnotationPosition position) {
        BitSet pos = this.getAtomAnnotationPositions(atom);
        pos.set(position.ordinal());
    }

    public AnnotationPosition alignmentToAnnotationPosition(int align) {
        switch (align) {
            case 1: {
                return AnnotationPosition.E;
            }
            case -1: {
                return AnnotationPosition.W;
            }
            case -2: {
                return AnnotationPosition.N;
            }
            case 2: {
                return AnnotationPosition.S;
            }
        }
        return AnnotationPosition.E;
    }

    public void addBondToAtomAnnotationPositions(IBond bond) {
        IAtom atom0 = bond.getAtom(0);
        IAtom atom1 = bond.getAtom(1);
        BitSet positions = this.getAtomAnnotationPositions(atom0);
        AnnotationPosition bondPosition = this.calculateBondPosition(atom0, atom1);
        positions.set(bondPosition.ordinal());
        positions = this.getAtomAnnotationPositions(atom1);
        bondPosition = this.calculateBondPosition(atom1, atom0);
        positions.set(bondPosition.ordinal());
    }

    public AnnotationPosition calculateBondPosition(IAtom atomFrom, IAtom atomTo) {
        AnnotationPosition pos = this.calculateRelativePosition(atomFrom.getPoint2d(), atomTo.getPoint2d());
        return pos;
    }

    public AnnotationPosition calculateRelativePosition(Point2d fromPoint, Point2d toPoint) {
        Vector2d bondVector = new Vector2d((Tuple2d) toPoint);
        bondVector.sub((Tuple2d) fromPoint);
        bondVector.normalize();
        double xAng = Math.toDegrees(bondVector.angle(POS_X));
        double yAng = Math.toDegrees(bondVector.angle(POS_Y));
        if (xAng < 22.5 && yAng > 67.5 && yAng < 115.5) {
            return AnnotationPosition.E;
        }
        if (xAng > 22.5 && xAng < 67.5 && yAng > 115.5 && yAng < 155.5) {
            return AnnotationPosition.NE;
        }
        if (xAng > 67.5 && xAng < 115.5 && yAng > 155.5) {
            return AnnotationPosition.N;
        }
        if (xAng > 115.5 && xAng < 155.5 && yAng > 115.5 && yAng < 155.5) {
            return AnnotationPosition.NW;
        }
        if (xAng > 155.5 && yAng > 67.5 && yAng < 115.5) {
            return AnnotationPosition.W;
        }
        if (xAng > 115.5 && xAng < 155.5 && yAng > 22.5 && yAng < 67.5) {
            return AnnotationPosition.SW;
        }
        if (xAng > 67.5 && xAng < 115.5 && yAng < 22.5) {
            return AnnotationPosition.S;
        }
        if (xAng > 22.5 && xAng < 67.5 && yAng > 22.5 && yAng < 67.5) {
            return AnnotationPosition.SE;
        }
        return AnnotationPosition.E;
    }

    private void blockRingSegment(IAtom atom, List<AnnotationPosition> ringPositions) {
        AnnotationPosition b;
        BitSet positions = this.getAtomAnnotationPositions(atom);
        if (ringPositions.size() != 2) {
            return;
        }
        AnnotationPosition a = ringPositions.get(0);
        if (this.positionsEqual(a, b = ringPositions.get(1), AnnotationPosition.N, AnnotationPosition.SW)) {
            positions.set(AnnotationPosition.NW.ordinal());
            positions.set(AnnotationPosition.W.ordinal());
        } else if (this.positionsEqual(a, b, AnnotationPosition.N, AnnotationPosition.SE)) {
            positions.set(AnnotationPosition.NE.ordinal());
            positions.set(AnnotationPosition.E.ordinal());
        } else if (this.positionsEqual(a, b, AnnotationPosition.NW, AnnotationPosition.S)) {
            positions.set(AnnotationPosition.W.ordinal());
            positions.set(AnnotationPosition.SW.ordinal());
        } else if (this.positionsEqual(a, b, AnnotationPosition.NE, AnnotationPosition.S)) {
            positions.set(AnnotationPosition.E.ordinal());
            positions.set(AnnotationPosition.SE.ordinal());
        } else if (this.positionsEqual(a, b, AnnotationPosition.W, AnnotationPosition.SE)) {
            positions.set(AnnotationPosition.SW.ordinal());
            positions.set(AnnotationPosition.S.ordinal());
        } else if (this.positionsEqual(a, b, AnnotationPosition.E, AnnotationPosition.SW)) {
            positions.set(AnnotationPosition.SE.ordinal());
            positions.set(AnnotationPosition.S.ordinal());
        } else if (this.positionsEqual(a, b, AnnotationPosition.NW, AnnotationPosition.E)) {
            positions.set(AnnotationPosition.N.ordinal());
            positions.set(AnnotationPosition.NE.ordinal());
        } else if (this.positionsEqual(a, b, AnnotationPosition.NE, AnnotationPosition.W)) {
            positions.set(AnnotationPosition.NW.ordinal());
            positions.set(AnnotationPosition.N.ordinal());
        } else if (this.positionsEqual(a, b, AnnotationPosition.NW, AnnotationPosition.NE)) {
            positions.set(AnnotationPosition.N.ordinal());
        } else if (this.positionsEqual(a, b, AnnotationPosition.SW, AnnotationPosition.SE)) {
            positions.set(AnnotationPosition.S.ordinal());
        } else if (this.positionsEqual(a, b, AnnotationPosition.NW, AnnotationPosition.SW)) {
            positions.set(AnnotationPosition.W.ordinal());
        } else if (this.positionsEqual(a, b, AnnotationPosition.NE, AnnotationPosition.SE)) {
            positions.set(AnnotationPosition.E.ordinal());
        }
    }

    private boolean positionsEqual(AnnotationPosition a, AnnotationPosition b, AnnotationPosition c, AnnotationPosition d) {
        return a == c && b == d || a == d && b == c;
    }

    public void addRingCenterToAtomAnnotationPosition(IAtom atom, List<IAtom> connectedAtomsInRing) {
        Point2d p1 = atom.getPoint2d();
        ArrayList<AnnotationPosition> ringPositions = new ArrayList<>();
        for (IAtom connectedAtom : connectedAtomsInRing) {
            Point2d p2 = connectedAtom.getPoint2d();
            ringPositions.add(this.calculateRelativePosition(p1, p2));
        }
        this.blockRingSegment(atom, ringPositions);
    }

    public boolean isUsed(IAtom atom, AnnotationPosition suggestedPosition) {
        int index = suggestedPosition.ordinal();
        return this.getAtomAnnotationPositions(atom).get(index);
    }

    public void reset() {
        this.atomAnnotationPositions.clear();
    }

    public static enum AnnotationPosition {
        N,
        W,
        S,
        E,
        NW,
        NE,
        SW,
        SE;

        private AnnotationPosition() {
        }
    }

}
