/*
 * SPDX-License-Identifier: Apache-2.0
 * © 2025 BioInception PVT LTD.
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
