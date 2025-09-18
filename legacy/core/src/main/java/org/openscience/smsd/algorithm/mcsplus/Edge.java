/*
 * SPDX-License-Identifier: Apache-2.0
 * © 2025 BioInception PVT LTD.
 */
package org.openscience.smsd.algorithm.mcsplus;

import java.io.Serializable;
import java.util.Comparator;

/**
 * This class generates compatibility graph between query and target molecule. It also marks edges in the compatibility
 * graph as c-edges or d-edges.
 *
 * 
 * 
 *
 * @author Syed Asad Rahman <asad at ebi.ac.uk>
 */
public class Edge implements Comparable<Edge>, Comparator<Edge>, Serializable {

    private static final long serialVersionUID = 52343464641L;

    @Override
    public int hashCode() {
        int hash = 5;
        hash = 23 * hash + this.i;
        hash = 23 * hash + this.j;
        return hash;
    }

    @Override
    public boolean equals(Object obj) {
        if (obj == null) {
            return false;
        }
        if (getClass() != obj.getClass()) {
            return false;
        }
        final Edge other = (Edge) obj;
        if (this.i != other.i && this.j == other.j) {
            return false;
        }
        if (this.i == other.i && this.j != other.j) {
            return false;
        }
        return this.i == other.i || this.j == other.j;
    }

    private final int i;
    private final int j;

    /**
     *
     * @param i
     * @param j
     */
    public Edge(int i, int j) {
        this.i = i;
        this.j = j;
    }

    @Override
    public int compareTo(Edge o) {
        String a = this.getSource() + "_" + this.getSink();
        String b = o.getSource() + "_" + o.getSink();

        return a.compareTo(b);
    }

    @Override
    public int compare(Edge o1, Edge o2) {
        String a = o1.getSource() + "_" + o1.getSink();
        String b = o2.getSource() + "_" + o2.getSink();

        return a.compareTo(b);
    }

    /**
     * @return the i
     */
    public int getSource() {
        return i;
    }

    /**
     * @return the j
     */
    public int getSink() {
        return j;
    }
}
