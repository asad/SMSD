/*
 * SPDX-License-Identifier: Apache-2.0
 * © 2025 BioInception PVT LTD.
 */
package org.openscience.directgraphics.direct.layout;

import java.awt.Dimension;
import java.util.List;
import javax.vecmath.Point2d;
import org.openscience.cdk.interfaces.IAtomContainer;

public class GridCanvasGenerator
extends AbstractCanvasGenerator
implements CanvasGenerator {
    private int rows;
    private int cols;
    private Dimension size;

    public GridCanvasGenerator() {
        this(1, 1);
    }

    public GridCanvasGenerator(int rows, int cols) {
        this.rows = rows;
        this.cols = cols;
    }

    @Override
    public void layout(List<IAtomContainer> atomContainers, Dimension cellCanvas) {
        double w = cellCanvas.width;
        double h = cellCanvas.height;
        double centerX = w / 2.0;
        double centerY = h / 2.0;
        int colCounter = 0;
        int rowCounter = 0;
        for (IAtomContainer atomContainer : atomContainers) {
            this.createCanvas(atomContainer, new Point2d(centerX, centerY), cellCanvas);
            if (++colCounter < this.cols) {
                centerX += w;
            } else {
                centerY += h;
                centerX = w / 2.0;
                colCounter = 0;
                ++rowCounter;
            }
            if (rowCounter <= this.rows) continue;
            System.err.println("WARNING : Row limit exceeded");
        }
        this.size = new Dimension(this.cols * cellCanvas.width, this.rows * cellCanvas.height);
    }

    @Override
    public Dimension getSize() {
        return this.size;
    }
}

