/*
 * SPDX-License-Identifier: Apache-2.0
 * © 2025 BioInception PVT LTD.
 */
package org.openscience.smsd.tools;

import java.io.Serializable;

/**
 * Class that handles execution time of the MCS search.
 *
 *
 * 
 * 
 * @author Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>
 */
public class IterationManager implements Serializable {

    private static final long serialVersionUID = 396239639826981L;
    private int max;
    private int counter;
    private int coverage;
    private final int limit;

    /**
     * Constructor for storing execution time
     */
    public IterationManager() {
        this(Integer.MAX_VALUE);
    }

    /**
     * Constructor for storing execution time
     *
     * @param maxIteration
     */
    public IterationManager(int maxIteration) {
        this.counter = 0;
        this.coverage = 250;
        this.max = maxIteration;
        this.limit = this.max * this.coverage;
    }

    /**
     * Returns Number of iterations
     *
     * @return Number of iterations
     */
    public synchronized int getCounter() {
        return counter;
    }

    /**
     * increment the counter
     *
     *
     */
    public synchronized void increment() {
        counter++;
    }

    /**
     * decrement the counter
     *
     *
     */
    public synchronized void decrement() {
        counter--;
    }

    public synchronized boolean isMaxIteration() {
        return getCounter() > limit;
    }

    /**
     * @return the coverage
     */
    public synchronized int getCoverage() {
        return coverage;
    }

    /**
     * @param coverage the coverage to set
     */
    public synchronized void setCoverage(int coverage) {
        this.coverage = coverage;
    }

    /**
     * Returns max allowed iterations (upper limit)
     *
     * @return
     */
    public int getIterationLimit() {
        return limit;
    }
}
