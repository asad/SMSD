/**
 *
 * Copyright (C) 2009-2018 Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>
 *
 * Contact: cdk-devel@lists.sourceforge.net
 *
 * This program is free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation; either version 2.1 of the License, or (at your option)
 * any later version. All we ask is that proper credit is given for our work,
 * which includes - but is not limited to - adding the above copyright notice to
 * the beginning of your source code files, and to any copyright notice that you
 * may distribute with programs based on this work.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
 * details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA.
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
