/*
 * SPDX-License-Identifier: Apache-2.0
 * © 2025 BioInception PVT LTD.
 */
package org.openscience.smsd.interfaces;

import java.util.Iterator;
import java.util.List;
import java.util.Map;

/**
 * Interface for mappings.
 * 
 * 
 * @author Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>
 */

public interface IFinalMapping {

   
    /**
     * Adds mapping to the mapping list
     * @param mapping List of all MCS mapping between a given
     * reactant and product 
     */
    public void add(Map<Integer, Integer> mapping);

    /**
     * Sets mapping list
     * @param list List of all MCS mapping between a given
     * reactant and product 
     */
    public void set(List<Map<Integer, Integer>> list);

    /**
     * Returns a mapping Iterator
     * @return Iterator of mappings
     */
    public Iterator<Map<Integer, Integer>> getIterator();

    /**
     * clear the mapping
     */
    public void clear();

    /**
     * Returns the stored mappings
     * @return get of MCS mapping List
     */
    public List<Map<Integer, Integer>> getFinalMapping();

    /**
     * Returns number of stored mappings
     * @return size of the mapping
     */
    public int getSize();
}
