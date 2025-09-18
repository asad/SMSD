/*
 * SPDX-License-Identifier: Apache-2.0
 * © 2025 BioInception PVT LTD.
 */
package org.openscience.smsd.interfaces;

import java.util.List;
import org.openscience.smsd.AtomAtomMapping;

/**
 * Interface that holds basic core interface for all MCS algorithm.  
 *
 * @author Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>
 */
public interface IResults {

    /**
     * Returns all plausible mappings between query and target molecules. Each map in the list has atom-atom equivalence
     * of the mappings between query and target molecule i.e. map.getKey() for the query and map.getValue() for the
     * target molecule
     *
     * @return All possible MCS atom Mappings
     */
    public abstract List<AtomAtomMapping> getAllAtomMapping();

    /**
     * Returns one of the best matches with atoms mapped.
     *
     * @return Best Atom Mapping
     */
    public abstract AtomAtomMapping getFirstAtomMapping();
}
