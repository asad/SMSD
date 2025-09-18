/*
 * SPDX-License-Identifier: Apache-2.0
 * © 2025 BioInception PVT LTD.
 */
package org.openscience.directgraphics.direct;

import java.awt.geom.Rectangle2D;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import javax.vecmath.Point2d;
import javax.vecmath.Tuple2d;
import javax.vecmath.Vector2d;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;

public class ConvexHull
        implements Iterable<Point2d> {

    private Point2d[] hull;
    private String[] hullIDs;
    private final Vector2d X_AXIS = new Vector2d(1.0, 0.0);

    public ConvexHull(IAtomContainer atomContainer) {
        Point2d[] points = new Point2d[atomContainer.getAtomCount()];
        int i = 0;
        for (IAtom atom : atomContainer.atoms()) {
            points[i] = atom.getPoint2d();
            ++i;
        }
        if (i < atomContainer.getAtomCount()) {
            Point2d[] nonNullPoints = new Point2d[i];
            int k = 0;
            for (int j = 0; j < points.length; ++j) {
                if (points[j] == null) {
                    continue;
                }
                nonNullPoints[k] = points[j];
                ++k;
            }
            points = nonNullPoints;
        }
        this.makeFromPoints(points);
        this.hullIDs = new String[this.hull.length];
        block2:
        for (IAtom atom2 : atomContainer.atoms()) {
            if (atom2.getPoint2d() == null || atom2.getID() == null) {
                continue;
            }
            Point2d point = atom2.getPoint2d();
            String id = atom2.getID();
            int hullIndex = 0;
            for (Point2d hullPoint : this.hull) {
                if (hullPoint == point) {
                    this.hullIDs[hullIndex] = id;
                    continue block2;
                }
                ++hullIndex;
            }
        }
    }

    public ConvexHull(Point2d[] points) {
        this.makeFromPoints(points);
    }

    public Vector2d getMajorAxis() {
        Rectangle minimumAreaBoundingRectangle = this.getMinimumAreaBoundingRectangleBruteForce();
        return minimumAreaBoundingRectangle.getMajorAxis();
    }

    public Point2d getCenter() {
        Point2d center = new Point2d();
        for (Point2d hullPoint : this.hull) {
            center.x += hullPoint.x;
            center.y += hullPoint.y;
        }
        center.x /= (double) this.hull.length;
        center.y /= (double) this.hull.length;
        return center;
    }

    public Rectangle getMinimumAreaBoundingRectangleBruteForce() {
        Rectangle minRect = null;
        double minArea = -1.0;
        int winnerIndex = -1;
        for (int index = 0; index < this.hull.length - 1; ++index) {
            Vector2d edge = this.edgeVector(this.hull[index], this.hull[index + 1]);
            Rectangle rect = this.getRectangleBrute(edge, index, index + 1);
            double area = rect.area();
            if (minRect != null && area >= minArea) {
                continue;
            }
            minRect = rect;
            minArea = area;
            winnerIndex = index;
        }
        Vector2d edge = this.edgeVector(this.hull[this.hull.length - 1], this.hull[0]);
        Rectangle rect = this.getRectangleBrute(edge, this.hull.length - 1, 0);
        double area = rect.area();
        if (minRect == null || area < minArea) {
            minRect = rect;
            minArea = area;
            winnerIndex = this.hull.length;
        }
        return minRect;
    }

    public Rectangle getMinimumAreaBoundingRectangle() {
        double minAngle;
        assert (this.hull != null);
        Point2d minY = null;
        Point2d maxY = null;
        int indexA = -1;
        int indexB = -1;
        for (int index = 0; index < this.hull.length; ++index) {
            Point2d point = this.hull[index];
            if (minY == null || point.y < minY.y) {
                minY = point;
                indexA = index;
            }
            if (maxY != null && point.y <= maxY.y) {
                continue;
            }
            maxY = point;
            indexB = index;
        }
        Vector2d caliperA = new Vector2d(1.0, 0.0);
        Vector2d caliperB = new Vector2d(-1.0, 0.0);
        double minArea = Double.MAX_VALUE;
        Rectangle minRect = null;
        for (double rotatedAngle = 0.0; rotatedAngle < 3.141592653589793; rotatedAngle += minAngle) {
            if (indexA == this.hull.length - 1) {
                indexA = 0;
            }
            if (indexB == this.hull.length - 1) {
                indexB = 0;
            }
            Vector2d edgeA = this.edgeVector(this.hull[indexA], this.hull[indexA + 1]);
            Vector2d edgeB = this.edgeVector(this.hull[indexB], this.hull[indexB + 1]);
            double angleA = edgeA.angle(caliperA);
            double angleB = edgeB.angle(caliperB);
            minAngle = Math.min(angleA, angleB);
            caliperA = this.rotate(caliperA, minAngle);
            caliperB = this.rotate(caliperB, minAngle);
            Rectangle rectangle = angleA < angleB ? this.getRectangle(edgeA, indexA, ++indexA + 1) : this.getRectangle(edgeB, indexB, ++indexB + 1);
            double area = rectangle.area();
            if (area >= minArea) {
                continue;
            }
            minArea = area;
            minRect = rectangle;
        }
        return minRect;
    }

    private Rectangle getRectangleBrute(Vector2d vector, int tailPointIndex, int headPointIndex) {
        Point2d headPoint = this.hull[headPointIndex];
        Point2d tailPoint = this.hull[tailPointIndex];
        int index = headPointIndex;
        int min = 0;
        Point2d vMax = null;
        Point2d thirdPoint = null;
        double thirdPointDist = 0.0;
        double minAngle = 6.283185307179586;
        double maxAngle = 0.0;
        Vector2d vN = new Vector2d(vector);
        vN.normalize();
        int max = 0;
        for (int visited = 0; visited < this.hull.length; ++visited) {
            if (index == this.hull.length) {
                index = 0;
            }
            if (vMax == null) {
                vMax = this.hull[index];
            } else {
                double angle = this.prj(tailPoint, headPoint, this.hull[index]);
                if (angle < minAngle) {
                    min = index;
                    minAngle = angle;
                }
                if (angle > maxAngle) {
                    vMax = this.hull[index];
                    max = index;
                    maxAngle = angle;
                }
            }
            if (thirdPoint == null) {
                thirdPoint = this.hull[index];
            } else {
                double d = this.pointLineDistance(tailPoint, headPoint, this.hull[index]);
                if (d > thirdPointDist) {
                    thirdPointDist = d;
                    thirdPoint = this.hull[index];
                }
            }
            ++index;
        }
        Point2d vMin = this.hull[min];
        Point2d tailProj = this.project(tailPoint, headPoint, vMax, true);
        Point2d headProj = this.project(tailPoint, headPoint, vMin, true);
        Rectangle r = new Rectangle(thirdPoint, tailProj, headProj, thirdPointDist);
        r.pointY = vMin;
        r.pointZ = vMax;
        return r;
    }

    private Rectangle getRectangle(Vector2d vector, int tailPointIndex, int headPointIndex) {
        Point2d headPoint = this.hull[headPointIndex];
        Point2d tailPoint = this.hull[tailPointIndex];
        int tailExPtIndex = tailPointIndex;
        Point2d tailExPt = this.hull[tailExPtIndex];
        boolean increasing = true;
        double proj = this.pointLineDistance(tailPoint, headPoint, tailExPt);
        while (increasing) {
            int nextIndex = tailExPtIndex > 0 ? tailExPtIndex - 1 : this.hull.length - 1;
            Point2d nextPoint = this.hull[nextIndex];
            double nextProj = this.pointLineDistance(tailPoint, headPoint, nextPoint);
            if (nextProj > proj) {
                proj = nextProj;
                tailExPtIndex = nextIndex;
                tailExPt = nextPoint;
                continue;
            }
            increasing = false;
        }
        Vector2d negV = new Vector2d(vector);
        negV.negate();
        Point2d projTail = this.project(tailPoint, headPoint, tailExPt);
        int headExPtIndex = headPointIndex;
        Point2d headExPt = this.hull[headExPtIndex];
        increasing = true;
        proj = this.pointLineDistance(tailPoint, headPoint, headExPt);
        while (increasing) {
            int nextIndex = headExPtIndex < this.hull.length - 1 ? headExPtIndex + 1 : 0;
            Point2d nextPoint = this.hull[nextIndex];
            double nextProj = this.pointLineDistance(tailPoint, headPoint, nextPoint);
            if (nextProj > proj) {
                proj = nextProj;
                headExPtIndex = nextIndex;
                headExPt = nextPoint;
                continue;
            }
            increasing = false;
        }
        Point2d projHead = this.project(tailPoint, headPoint, headExPt);
        int remainExPtIndex = headExPtIndex;
        Point2d remainExPoint = this.hull[remainExPtIndex];
        increasing = true;
        double dist = this.pointLineDistance(tailPoint, headPoint, remainExPoint);
        while (increasing) {
            int nextIndex = remainExPtIndex < this.hull.length - 1 ? remainExPtIndex + 1 : 0;
            Point2d nextPoint = this.hull[nextIndex];
            double nextDistance = this.pointLineDistance(tailPoint, headPoint, nextPoint);
            if (nextDistance > dist) {
                dist = nextDistance;
                remainExPtIndex = nextIndex;
                remainExPoint = nextPoint;
                continue;
            }
            increasing = false;
        }
        return new Rectangle(remainExPoint, projTail, projHead, this.pointLineDistance(tailPoint, headPoint, remainExPoint));
    }

    public /* varargs */ String toString(Point2d... points) {
        String str = "[";
        for (Point2d point : points) {
            str = str + String.format("(%2.0f, %2.0f)", point.x, point.y);
        }
        return str + "]";
    }

    private Point2d project(Point2d p1, Point2d p2, Point2d p3) {
        return this.project(p1, p2, p3, false);
    }

    private Point2d project(Point2d p1, Point2d p2, Point2d p3, boolean outSeg) {
        double dx = p2.x - p1.x;
        double dy = p2.y - p1.y;
        if (dx == 0.0 && dy == 0.0) {
            return new Point2d(p1);
        }
        double t = ((p3.x - p1.x) * dx + (p3.y - p1.y) * dy) / (dx * dx + dy * dy);
        Point2d p = outSeg && t > 0.0 && t < 1.0 ? (t > 0.5 ? p2 : p1) : new Point2d(p1.x + t * dx, p1.y + t * dy);
        return p;
    }

    private double prj(Point2d p1, Point2d p2, Point2d p3) {
        double dx = p2.x - p1.x;
        double dy = p2.y - p1.y;
        return ((p3.x - p1.x) * dx + (p3.y - p1.y) * dy) / (dx * dx + dy + dy);
    }

    private double pointLineDistance(Point2d p1, Point2d p2, Point2d p3) {
        Point2d p = this.project(p1, p2, p3);
        return p3.distance(p);
    }

    private Vector2d rotate(Vector2d vector, double angle) {
        Vector2d rotatedVector = new Vector2d();
        double cosTh = Math.cos(angle);
        double sinTh = Math.sin(angle);
        rotatedVector.x = cosTh * vector.x - sinTh * vector.y;
        rotatedVector.y = sinTh * vector.x + cosTh * vector.y;
        return rotatedVector;
    }

    private Vector2d edgeVector(Point2d fromPoint, Point2d toPoint) {
        Vector2d edge = new Vector2d((Tuple2d) fromPoint);
        edge.sub((Tuple2d) toPoint);
        return edge;
    }

    public Rectangle2D getAxisAlignedMinimumBoundingRectangle() {
        double minX = Double.MAX_VALUE;
        double minY = Double.MAX_VALUE;
        double maxX = Double.MIN_VALUE;
        double maxY = Double.MIN_VALUE;
        for (Point2d point : this.hull) {
            if (point.x < minX) {
                minX = point.x;
            }
            if (point.y < minY) {
                minY = point.y;
            }
            if (point.x > maxX) {
                maxX = point.x;
            }
            if (point.y <= maxY) {
                continue;
            }
            maxY = point.y;
        }
        return new Rectangle2D.Double(minX, minY, maxX, maxY);
    }

    private void makeFromPoints(Point2d[] points) {
        if (points.length < 4) {
            this.hull = points;
            return;
        }
        int indexOfLowPoint = -1;
        Point2d lowPoint = null;
        for (int index = 0; index < points.length; ++index) {
            Point2d current = points[index];
            if (indexOfLowPoint != -1 && current.y <= lowPoint.y) {
                continue;
            }
            lowPoint = current;
            indexOfLowPoint = index;
        }
        this.swap(points, 0, indexOfLowPoint);
        points = this.sortByPolarAngle(points);
        int m = 2;
        for (int i = 3; i < points.length; ++i) {
            while (this.ccw(points[m - 1], points[m], points[i]) <= 0.0) {
                if (m == 2) {
                    this.swap(points, m, i);
                    ++i;
                    continue;
                }
                --m;
            }
            this.swap(points, ++m, i);
        }
        this.hull = new Point2d[m];
        System.arraycopy(points, 0, this.hull, 0, m);
    }

    private Point2d[] sortByPolarAngle(Point2d[] points) {
        Point2d ref = points[0];
        final HashMap<Point2d, Double> angles = new HashMap<>();
        angles.put(ref, 0.0);
        for (int pointIndex = 1; pointIndex < points.length; ++pointIndex) {
            Point2d point = points[pointIndex];
            double angle = this.getAngle(ref, point);
            angles.put(point, angle);
        }
        Arrays.sort(points, new Comparator<Point2d>() {

            @Override
            public int compare(Point2d p0, Point2d p1) {
                return (angles.get(p0)).compareTo(angles.get(p1));
            }
        });
        Point2d[] sortedPoints = new Point2d[points.length + 1];
        sortedPoints[0] = points[points.length - 1];
        System.arraycopy(points, 0, sortedPoints, 1, points.length);
        return sortedPoints;
    }

    private double getAngle(Point2d ref, Point2d point) {
        Vector2d rp = new Vector2d((Tuple2d) ref);
        rp.sub((Tuple2d) point);
        rp.normalize();
        return this.X_AXIS.angle(rp);
    }

    private void swap(Point2d[] points, int i, int j) {
        Point2d tmp = points[i];
        points[i] = points[j];
        points[j] = tmp;
    }

    private double ccw(Point2d p1, Point2d p2, Point2d p3) {
        return (p2.x - p1.x) * (p3.y - p1.y) - (p2.y - p1.y) * (p3.x - p1.x);
    }

    @Override
    public Iterator<Point2d> iterator() {
        return Arrays.asList(this.hull).iterator();
    }

    public class Rectangle {

        public Point2d pointX;
        public Point2d pointY;
        public Point2d pointZ;
        public Point2d cornerA;
        public Point2d cornerB;
        public Point2d cornerC;
        public Point2d cornerD;

        public Rectangle(Point2d pointOnAB, Point2d cornerC, Point2d cornerD, double distToCD) {
            this.pointX = new Point2d(pointOnAB);
            this.cornerC = new Point2d(cornerC);
            this.cornerD = new Point2d(cornerD);
            Vector2d cdVec = new Vector2d((Tuple2d) cornerD);
            cdVec.sub((Tuple2d) cornerC);
            Vector2d cdVecNormalized = new Vector2d(cdVec);
            if (cdVec.x != 0.0 && cdVec.y != 0.0) {
                cdVecNormalized.normalize();
            }
            Vector2d perp = new Vector2d(cdVecNormalized.y, -cdVecNormalized.x);
            this.cornerA = new Point2d(cornerD);
            this.cornerA.scaleAdd(distToCD, (Tuple2d) perp, (Tuple2d) this.cornerA);
            this.cornerB = new Point2d(cornerC);
            this.cornerB.scaleAdd(distToCD, (Tuple2d) perp, (Tuple2d) this.cornerB);
        }

        public double area() {
            return new Point2d(this.cornerA).distance(new Point2d(this.cornerC)) * new Point2d(this.cornerC).distance(new Point2d(this.cornerD));
        }

        @Override
        public String toString() {
            return String.format("[(%2.0f, %2.0f), (%2.0f, %2.0f), (%2.0f, %2.0f), (%2.0f, %2.0f)]", this.cornerA.x, this.cornerA.y, this.cornerB.x, this.cornerB.y, this.cornerC.x, this.cornerC.y, this.cornerD.x, this.cornerD.y);
        }

        public double getWidth() {
            Vector2d cd = new Vector2d((Tuple2d) this.cornerC);
            cd.sub((Tuple2d) this.cornerD);
            return cd.length();
        }

        public Vector2d getMajorAxis() {
            Vector2d cd = new Vector2d((Tuple2d) this.cornerC);
            cd.sub((Tuple2d) this.cornerD);
            double cdLen = cd.length();
            Vector2d ad = new Vector2d((Tuple2d) this.cornerA);
            ad.sub((Tuple2d) this.cornerD);
            double adLen = ad.length();
            if (adLen > cdLen) {
                return ad;
            }
            return cd;
        }

        public double getHeight() {
            Vector2d ac = new Vector2d((Tuple2d) this.cornerA);
            ac.sub((Tuple2d) this.cornerC);
            return ac.length();
        }
    }

}
