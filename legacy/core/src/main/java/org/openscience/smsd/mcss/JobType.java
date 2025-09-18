/*
 * SPDX-License-Identifier: Apache-2.0
 * © 2025 BioInception PVT LTD.
 */
package org.openscience.smsd.mcss;

/**
 * 
 * 
 * @author Syed Asad Rahman <asad @ ebi.ac.uk>
 *
 */
public enum JobType implements Comparable<JobType> {

    /**
     * Default MULTIPLE algorithm.
     */
    MULTIPLE(0, "Multiple Fragments"),
    /**
     * SINGLE search algorithm.
     */
    SINGLE(1, "Single Fragment");
    private final int type;
    private final String description;

    JobType(int aStatus, String desc) {
        this.type = aStatus;
        this.description = desc;
    }

    /**
     * Returns type of algorithm.
     *
     * @return type of algorithm
     */
    public int type() {
        return this.type;
    }

    /**
     * Returns short description of the algorithm.
     *
     * @return description of the algorithm
     */
    public String description() {
        return this.description;
    }
}
