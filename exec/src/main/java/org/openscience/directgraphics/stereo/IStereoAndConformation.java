/*
 * Copyright (C) 2007-2015 Syed Asad Rahman <asad @ ebi.ac.uk>.
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
package org.openscience.directgraphics.stereo;

public enum IStereoAndConformation implements Comparable<IStereoAndConformation> {
    NONE(0, "CHIRALITY NONE"),
    R(1, "CHIRALITY Rectus"),
    S(2, "CHIRALITY Sinister"),
    EITHER(3, "CHIRALITY R or S"),
    M(4, "CHIRALITY M Configuration"),
    P(5, "CHIRALITY P Configuration"),
    Z(6, "TOGETHER atom Configuration"),
    E(7, "OPPOSITE atom Configuration");

    private final int type;
    private final String description;

    private IStereoAndConformation(int aStatus, String desc) {
        this.type = aStatus;
        this.description = desc;
    }

    public int type() {
        return this.type;
    }

    public String description() {
        return this.description;
    }
}
