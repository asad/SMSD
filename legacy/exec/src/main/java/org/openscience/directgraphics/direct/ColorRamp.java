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
import java.util.ArrayList;
import java.util.List;

public class ColorRamp {

    public static List<Color> getColors(int number) {
        ArrayList<Color> colors = new ArrayList<>();
        for (int i = 0; i < number; ++i) {
            colors.add(ColorRamp.colorRamp(i, 0, number));
        }
        return colors;
    }

    public static Color colorRamp(int v, int vmin, int vmax) {
        double r = 1.0;
        double g = 1.0;
        double b = 1.0;
        if (v < vmin) {
            v = vmin;
        }
        if (v > vmax) {
            v = vmax;
        }
        int dv = vmax - vmin;
        try {
            if ((double) v < (double) vmin + 0.25 * (double) dv) {
                r = 0.0;
                g = 4.0 * (double) (v - vmin) / (double) dv;
            } else if ((double) v < (double) vmin + 0.5 * (double) dv) {
                r = 0.0;
                b = 1.0 + 4.0 * ((double) vmin + 0.25 * (double) dv - (double) v) / (double) dv;
            } else if ((double) v < (double) vmin + 0.75 * (double) dv) {
                r = 4.0 * ((double) (v - vmin) - 0.5 * (double) dv) / (double) dv;
                b = 0.0;
            } else {
                g = 1.0 + 4.0 * ((double) vmin + 0.75 * (double) dv - (double) v) / (double) dv;
                b = 0.0;
            }
            float[] hsb = Color.RGBtoHSB((int) (r * 255.0), (int) (g * 255.0), (int) (b * 255.0), null);
            return Color.getHSBColor(hsb[0], hsb[1], hsb[2]);
        } catch (ArithmeticException zde) {
            float[] hsb = Color.RGBtoHSB(0, 0, 0, null);
            return Color.getHSBColor(hsb[0], hsb[1], hsb[2]);
        }
    }
}
