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

