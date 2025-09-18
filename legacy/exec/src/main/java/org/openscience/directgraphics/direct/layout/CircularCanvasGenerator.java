/*
 * SPDX-License-Identifier: Apache-2.0
 * © 2025 BioInception PVT LTD.
 */
package org.openscience.directgraphics.direct.layout;

import java.awt.Dimension;
import java.util.List;
import javax.vecmath.Point2d;
import javax.vecmath.Tuple2d;
import javax.vecmath.Vector2d;
import org.openscience.cdk.interfaces.IAtomContainer;

public class CircularCanvasGenerator
extends AbstractCanvasGenerator
implements CanvasGenerator {
    private Vector2d vectorToStart;
    private boolean putFirstInCenter;
    private Dimension size;

    public CircularCanvasGenerator() {
        this(new Vector2d(-1.0, 0.0));
    }

    public CircularCanvasGenerator(Vector2d vectorToStart) {
        this.vectorToStart = vectorToStart;
    }

    public CircularCanvasGenerator(boolean putFirstInCenter) {
        this(new Vector2d(-1.0, 0.0), putFirstInCenter);
    }

    public CircularCanvasGenerator(Vector2d vectorToStart, boolean putFirstInCenter) {
        this.vectorToStart = vectorToStart;
        this.putFirstInCenter = putFirstInCenter;
    }

    @Override
    public void layout(List<IAtomContainer> atomContainers, Dimension cellCanvas) {
        int index;
        int n = this.putFirstInCenter ? atomContainers.size() - 1 : atomContainers.size();
        if (n < 1) {
            return;
        }
        double maxDim = Math.max(cellCanvas.width, cellCanvas.height);
        double alpha = Math.toRadians(360 / n);
        double cosA = Math.cos(alpha);
        double sinA = Math.sin(alpha);
        double circleRadius = maxDim / 2.0 / Math.sin(alpha / 2.0);
        double totalDim = 2.0 * circleRadius + maxDim;
        this.size = new Dimension((int)totalDim, (int)totalDim);
        Point2d center = new Point2d(totalDim / 2.0, totalDim / 2.0);
        Vector2d v = new Vector2d(this.vectorToStart);
        v.normalize();
        if (this.putFirstInCenter) {
            this.createCanvas(atomContainers.get(0), center, cellCanvas);
            index = 1;
        } else {
            index = 0;
        }
        while (index < atomContainers.size()) {
            IAtomContainer atomContainer = atomContainers.get(index);
            Point2d canvasCenter = new Point2d(center);
            canvasCenter.scaleAdd(circleRadius, (Tuple2d)v, (Tuple2d)canvasCenter);
            this.createCanvas(atomContainer, canvasCenter, cellCanvas);
            Vector2d w = new Vector2d();
            w.x = cosA * v.x + sinA * v.y;
            w.y = (- sinA) * v.x + cosA * v.y;
            v = w;
            ++index;
        }
    }

    @Override
    public Dimension getSize() {
        return this.size;
    }
}

