/*
 * SPDX-License-Identifier: Apache-2.0
 * © 2025 BioInception PVT LTD.
 */
package org.openscience.directgraphics.direct;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Font;
import java.awt.Graphics2D;
import java.awt.Stroke;
import java.awt.font.FontRenderContext;
import java.awt.font.GlyphVector;
import java.awt.geom.AffineTransform;
import java.awt.geom.Line2D;
import java.awt.geom.Path2D;
import java.awt.geom.Rectangle2D;
import javax.vecmath.Point2d;
import javax.vecmath.Tuple2d;
import javax.vecmath.Vector2d;

public class DirectArrowDrawer
extends AbstractDirectDrawer {
    private static final Vector2d X_AXIS = new Vector2d(1.0, 0.0);
    private static final Vector2d Y_AXIS = new Vector2d(0.0, 1.0);
    private final Params paramsLocal;

    public DirectArrowDrawer(Params params) {
        this.paramsLocal = params;
    }

    public void drawArrow(Graphics2D g, Point2d c, Vector2d v) {
        Stroke savedStroke = g.getStroke();
        g.setStroke(new BasicStroke());
        if (this.paramsLocal.drawFatArrow) {
            if (this.paramsLocal.arrowType == Params.ArrowType.BIDIRECTIONAL) {
                this.drawDoubleHeadedFatArrow(g, c, v, null);
            } else {
                this.drawFatArrow(g, c, v, null);
            }
        } else {
            this.drawThinArrow(g, c, v, null);
        }
        g.setStroke(savedStroke);
    }

    public void drawFatArrow(Graphics2D g, Point2d c, Vector2d v, String text) {
        int arrowLength = this.paramsLocal.arrowLength;
        int arrowHeadLength = this.paramsLocal.arrowHeadLength;
        int arrowHeadIndent = this.paramsLocal.arrowHeadIndent;
        int arrowBodyWidth = this.paramsLocal.arrowBodyWidth;
        double arrowHeadAngleRad = Math.toRadians(this.paramsLocal.arrowHeadAngle);
        double arrowHeadAngleRadPrime = Math.toRadians(360.0 - this.paramsLocal.arrowHeadAngle);
        double cosA = Math.cos(arrowHeadAngleRad);
        double sinA = Math.sin(arrowHeadAngleRad);
        double cosAPrime = Math.cos(arrowHeadAngleRadPrime);
        double sinAPrime = Math.sin(arrowHeadAngleRadPrime);
        int halfLength = arrowLength / 2;
        int halfBodyWidth = arrowBodyWidth / 2;
        g.setColor(Color.BLACK);
        Vector2d nV = new Vector2d(v);
        nV.negate();
        Vector2d p = new Vector2d(v.y, - v.x);
        Vector2d nP = new Vector2d(- v.y, v.x);
        Point2d tail = new Point2d(c.x, c.y);
        tail.scaleAdd((double)halfLength, (Tuple2d)nV, (Tuple2d)c);
        Point2d upperTail = new Point2d(tail);
        upperTail.scaleAdd((double)halfBodyWidth, (Tuple2d)nP, (Tuple2d)upperTail);
        Point2d lowerTail = new Point2d(tail);
        lowerTail.scaleAdd((double)halfBodyWidth, (Tuple2d)p, (Tuple2d)lowerTail);
        Point2d head = new Point2d(c);
        head.scaleAdd((double)halfLength, (Tuple2d)v, (Tuple2d)c);
        Vector2d ccwVec = new Vector2d(cosA * nV.x + sinA * nV.y, cosA * nV.y - sinA * nV.x);
        Vector2d cwVec = new Vector2d(cosAPrime * nV.x + sinAPrime * nV.y, cosAPrime * nV.y - sinAPrime * nV.x);
        Point2d headCCW = new Point2d(head.x, head.y);
        headCCW.scaleAdd((double)arrowHeadLength, (Tuple2d)ccwVec, (Tuple2d)head);
        Point2d indentCCW = new Point2d(headCCW);
        indentCCW.scaleAdd((double)arrowHeadIndent, (Tuple2d)p, (Tuple2d)indentCCW);
        Point2d headCW = new Point2d(head.x, head.y);
        headCW.scaleAdd((double)arrowHeadLength, (Tuple2d)cwVec, (Tuple2d)head);
        Point2d indentCW = new Point2d(headCW);
        indentCW.scaleAdd((double)arrowHeadIndent, (Tuple2d)nP, (Tuple2d)indentCW);
        Path2D.Double polygon = new Path2D.Double();
        polygon.moveTo(head.x, head.y);
        polygon.lineTo(headCCW.x, headCCW.y);
        polygon.lineTo(indentCCW.x, indentCCW.y);
        polygon.lineTo(upperTail.x, upperTail.y);
        polygon.lineTo(lowerTail.x, lowerTail.y);
        polygon.lineTo(indentCW.x, indentCW.y);
        polygon.lineTo(headCW.x, headCW.y);
        polygon.closePath();
        if (this.paramsLocal.drawArrowFilled) {
            g.fill(polygon);
        } else {
            g.draw(polygon);
        }
        if (text != null) {
            this.drawText(g, text, c, v, nV);
        }
    }

    public void drawDoubleHeadedFatArrow(Graphics2D g, Point2d c, Vector2d v, String text) {
        int arrowLength = this.paramsLocal.arrowLength;
        int arrowHeadLength = this.paramsLocal.arrowHeadLength;
        int arrowHeadIndent = this.paramsLocal.arrowHeadIndent;
        double arrowHeadAngleRad = Math.toRadians(this.paramsLocal.arrowHeadAngle);
        double arrowHeadAngleRadPrime = Math.toRadians(360.0 - this.paramsLocal.arrowHeadAngle);
        double cosA = Math.cos(arrowHeadAngleRad);
        double sinA = Math.sin(arrowHeadAngleRad);
        double cosAPrime = Math.cos(arrowHeadAngleRadPrime);
        double sinAPrime = Math.sin(arrowHeadAngleRadPrime);
        int halfLength = arrowLength / 2;
        g.setColor(Color.BLACK);
        Vector2d nV = new Vector2d(v);
        nV.negate();
        Vector2d p = new Vector2d(v.y, - v.x);
        Vector2d nP = new Vector2d(- v.y, v.x);
        Point2d tail = new Point2d(c.x, c.y);
        tail.scaleAdd((double)halfLength, (Tuple2d)nV, (Tuple2d)c);
        Point2d head = new Point2d(c);
        head.scaleAdd((double)halfLength, (Tuple2d)v, (Tuple2d)c);
        Vector2d ccwVec = new Vector2d(cosA * nV.x + sinA * nV.y, cosA * nV.y - sinA * nV.x);
        Vector2d nCCWVec = new Vector2d(ccwVec);
        nCCWVec.negate();
        Vector2d cwVec = new Vector2d(cosAPrime * nV.x + sinAPrime * nV.y, cosAPrime * nV.y - sinAPrime * nV.x);
        Vector2d nCWVec = new Vector2d(cwVec);
        nCWVec.negate();
        Point2d headCCW = new Point2d(head.x, head.y);
        headCCW.scaleAdd((double)arrowHeadLength, (Tuple2d)ccwVec, (Tuple2d)head);
        Point2d headIndentCCW = new Point2d(headCCW);
        headIndentCCW.scaleAdd((double)arrowHeadIndent, (Tuple2d)p, (Tuple2d)headIndentCCW);
        Point2d headCW = new Point2d(head.x, head.y);
        headCW.scaleAdd((double)arrowHeadLength, (Tuple2d)cwVec, (Tuple2d)head);
        Point2d headIndentCW = new Point2d(headCW);
        headIndentCW.scaleAdd((double)arrowHeadIndent, (Tuple2d)nP, (Tuple2d)headIndentCW);
        Point2d tailCCW = new Point2d(tail);
        tailCCW.scaleAdd((double)arrowHeadLength, (Tuple2d)nCWVec, (Tuple2d)tailCCW);
        Point2d tailCW = new Point2d(tail);
        tailCW.scaleAdd((double)arrowHeadLength, (Tuple2d)nCCWVec, (Tuple2d)tailCW);
        Point2d upperTail = new Point2d(tailCCW);
        upperTail.scaleAdd((double)arrowHeadIndent, (Tuple2d)p, (Tuple2d)upperTail);
        Point2d lowerTail = new Point2d(tailCW);
        lowerTail.scaleAdd((double)arrowHeadIndent, (Tuple2d)nP, (Tuple2d)lowerTail);
        Path2D.Double polygon = new Path2D.Double();
        polygon.moveTo(head.x, head.y);
        polygon.lineTo(headCCW.x, headCCW.y);
        polygon.lineTo(headIndentCCW.x, headIndentCCW.y);
        polygon.lineTo(upperTail.x, upperTail.y);
        polygon.lineTo(tailCCW.x, tailCCW.y);
        polygon.lineTo(tail.x, tail.y);
        polygon.lineTo(tailCW.x, tailCW.y);
        polygon.lineTo(lowerTail.x, lowerTail.y);
        polygon.lineTo(headIndentCW.x, headIndentCW.y);
        polygon.lineTo(headCW.x, headCW.y);
        polygon.closePath();
        if (this.paramsLocal.drawArrowFilled) {
            g.fill(polygon);
        } else {
            g.draw(polygon);
        }
        if (text != null) {
            this.drawText(g, text, c, v, nV);
        }
    }

    public void drawThinArrow(Graphics2D g, Point2d c, Vector2d v, String text) {
        int arrowLength = this.paramsLocal.arrowLength;
        int arrowHeadLength = this.paramsLocal.arrowHeadLength;
        double arrowHeadAngleRad = Math.toRadians(this.paramsLocal.arrowHeadAngle);
        double arrowHeadAngleRadPrime = Math.toRadians(360.0 - this.paramsLocal.arrowHeadAngle);
        double cosA = Math.cos(arrowHeadAngleRad);
        double sinA = Math.sin(arrowHeadAngleRad);
        double cosAPrime = Math.cos(arrowHeadAngleRadPrime);
        double sinAPrime = Math.sin(arrowHeadAngleRadPrime);
        int halfLength = arrowLength / 2;
        g.setColor(Color.BLACK);
        Vector2d nV = new Vector2d(v);
        nV.negate();
        Point2d tail = new Point2d(c.x, c.y);
        tail.scaleAdd((double)halfLength, (Tuple2d)nV, (Tuple2d)c);
        Point2d head = new Point2d(c);
        head.scaleAdd((double)halfLength, (Tuple2d)v, (Tuple2d)c);
        Vector2d ccwVec = new Vector2d(cosA * nV.x + sinA * nV.y, cosA * nV.y - sinA * nV.x);
        Vector2d cwVec = new Vector2d(cosAPrime * nV.x + sinAPrime * nV.y, cosAPrime * nV.y - sinAPrime * nV.x);
        Point2d headCCW = new Point2d(head.x, head.y);
        headCCW.scaleAdd((double)arrowHeadLength, (Tuple2d)ccwVec, (Tuple2d)head);
        Point2d headCW = new Point2d(head.x, head.y);
        headCW.scaleAdd((double)arrowHeadLength, (Tuple2d)cwVec, (Tuple2d)head);
        g.draw(new Line2D.Double(tail.x, tail.y, head.x, head.y));
        g.draw(new Line2D.Double(head.x, head.y, headCCW.x, headCCW.y));
        g.draw(new Line2D.Double(head.x, head.y, headCW.x, headCW.y));
        if (text != null) {
            this.drawText(g, text, c, v, nV);
        }
    }

    private void drawText(Graphics2D g, String text, Point2d c, Vector2d v, Vector2d nV) {
        AffineTransform originalTransform = g.getTransform();
        double angle = this.getAngle(v);
        Rectangle2D textBounds = this.getTextBounds(g, text);
        double distance = textBounds.getWidth() / 2.0;
        Point2d start = new Point2d(c.x, c.y);
        if (angle < Math.toRadians(90.0) || angle > Math.toRadians(270.0)) {
            start.scaleAdd(distance, (Tuple2d)nV, (Tuple2d)c);
        } else {
            start.scaleAdd(distance, (Tuple2d)v, (Tuple2d)c);
            double angDeg = (180.0 + Math.toDegrees(angle)) % 360.0;
            angle = Math.toRadians(angDeg);
        }
        g.translate(start.x, start.y);
        g.rotate(angle);
        Font font = g.getFont();
        FontRenderContext frc = g.getFontRenderContext();
        GlyphVector gv = font.createGlyphVector(frc, text);
        int length = gv.getNumGlyphs();
        for (int i = 0; i < length; ++i) {
            g.fill(gv.getGlyphOutline(i));
        }
        g.rotate(6.283185307179586 - angle);
        g.setTransform(originalTransform);
    }

    private double getAngle(Vector2d v) {
        double xAngle = X_AXIS.angle(v);
        double yAngle = Y_AXIS.angle(v);
        if (xAngle < Math.toRadians(90.0)) {
            if (yAngle < Math.toRadians(90.0)) {
                return xAngle;
            }
            return Math.toRadians(360.0) - xAngle;
        }
        return Math.toRadians(90.0) + yAngle;
    }
}

