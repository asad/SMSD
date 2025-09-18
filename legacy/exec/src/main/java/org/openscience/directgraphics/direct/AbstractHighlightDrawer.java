/*
 * SPDX-License-Identifier: Apache-2.0
 * © 2025 BioInception PVT LTD.
 */
package org.openscience.directgraphics.direct;

import java.awt.Color;
import java.util.HashMap;
import java.util.Map;

public class AbstractHighlightDrawer
extends AbstractDirectDrawer {
    protected Color opaqueHighlightColor;
    protected Color translucentHighlightColor;
    private final Map<Color, Color> opaqueToTranslucentColorMap;

    public AbstractHighlightDrawer(Params params) {
        super.setParams(params);
        this.opaqueToTranslucentColorMap = new HashMap<>();
        this.opaqueHighlightColor = params.highlightColor;
        this.translucentHighlightColor = this.getTranslucentColor(this.opaqueHighlightColor);
    }

    public void registerColor(Color color) {
        if (this.opaqueToTranslucentColorMap.containsKey(color)) {
            return;
        }
        this.opaqueToTranslucentColorMap.put(color, this.makeTranslucentColor(color));
    }

    protected Color getTranslucentColor(Color color) {
        if (this.opaqueToTranslucentColorMap.containsKey(color)) {
            return this.opaqueToTranslucentColorMap.get(color);
        }
        Color translucentColor = this.makeTranslucentColor(color);
        this.opaqueToTranslucentColorMap.put(color, translucentColor);
        return translucentColor;
    }

    private Color makeTranslucentColor(Color color) {
        float[] c = color.getColorComponents(null);
        return new Color(c[0], c[1], c[2], this.params.highlightAlpha);
    }
}

