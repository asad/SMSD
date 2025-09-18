/*
 * SPDX-License-Identifier: Apache-2.0
 * © 2025 BioInception PVT LTD.
 */
package org.openscience.smsd.interfaces;


/**
 * 
 * This class represents various algorithm type supported by SMSD.
 * Presently SMSD supports 5 different kinds of algorithms:
 * 
 * <OL>
 * <lI>0: default,
 * <lI>1: MCSPlus,
 * <lI>2: VFLibMCS,
 * <lI>3: CDKMCS,
 * <lI>4: SubStructure
 * <lI>5: TurboSubStructure
 * </OL>
 *
 * 
 * 
 * @author Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>
 */
public enum Algorithm implements Comparable<Algorithm> {

    /**
     * Default SMSD algorithm.
     */
    DEFAULT(0, "Default SMSD algorithm"),
    /**
     * MCS Plus algorithm.
     */
    MCSPlus(2, "MCS Plus algorithm"),
    /**
     * VF-Koch-McGregor Lib based MCS algorithm.
     */
    VFLibMCS(3, "VF-Koch-McGregor Lib based MCS algorithm"),
    /**
     * CDK UIT MCS.
     */
    CDKMCS(4, "CDK UIT MCS");
    private final int type;
    private final String description;

    Algorithm(int aStatus, String desc) {
        this.type = aStatus;
        this.description = desc;
    }

    /**
     * Returns type of algorithm.
     * @return type of algorithm
     */
    public int type() {
        return this.type;
    }

    /**
     * Returns short description of the algorithm.
     * @return description of the algorithm
     */
    public String description() {
        return this.description;
    }
}
