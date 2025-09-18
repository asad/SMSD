/*
 * SPDX-License-Identifier: Apache-2.0
 * © 2025 BioInception PVT LTD.
 */
package org.openscience.smsd.algorithm.vflib.substructure;

/**
 * Holds source and target atoms
 * 
 * 
 * 
 * @author Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>
 */
/**
 * @param <T>
 * @param <S>
 */
public class Pair<T, S> {

    private T source;
    private S target;

    public Pair(T a, S b) {
        this.source = a;
        this.target = b;
    }

    @Override
    public synchronized String toString() {
        return "(" + getSourceAtom() + ", " + getTargetAtom() + ")";
    }

    /**
     * @return the source
     */
    public synchronized T getSourceAtom() {
        return source;
    }

    /**
     * @param source the source to set
     */
    public synchronized void setSourceAtom(T first) {
        this.source = first;
    }

    /**
     * @return the target
     */
    public synchronized S getTargetAtom() {
        return target;
    }

    /**
     * @param target the target to set
     */
    public synchronized void setTargetAtom(S second) {
        this.target = second;
    }
}
