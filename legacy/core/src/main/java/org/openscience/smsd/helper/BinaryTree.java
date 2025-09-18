/*
 * SPDX-License-Identifier: Apache-2.0
 * © 2025 BioInception PVT LTD.
 */
package org.openscience.smsd.helper;

/**
 * Class to construct a Binary tree for McGregor search.
 *
 *
 *
 * @author Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>
 */
public class BinaryTree {

    /**
     * Creates a new instance of BinaryTree.
     *
     * @param value node value
     */
    public BinaryTree(int value) {
        this.value = value;
    }
    /**
     * not equal is initialized as null
     */
    private BinaryTree equal = null;
    private BinaryTree notEqual = null;
    private int value = -1;

    /**
     * Return value of the node
     *
     * @return get the value of the current node
     */
    public synchronized int getValue() {
        return this.value;
    }

    /**
     * Returns equal node
     *
     * @return the equal
     */
    public synchronized BinaryTree getEqual() {
        return equal;
    }

    /**
     * Set equal node
     *
     * @param equal the equal to set
     */
    public synchronized void setEqual(BinaryTree equal) {
        this.equal = equal;
    }

    /**
     * Returns not equal node
     *
     * @return the notEqual
     */
    public synchronized BinaryTree getNotEqual() {
        return notEqual;
    }

    /**
     * Set not equal node
     *
     * @param notEqual the notEqual to set
     */
    public synchronized void setNotEqual(BinaryTree notEqual) {
        this.notEqual = notEqual;
    }
}
