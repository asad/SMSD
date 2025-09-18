/*
 * SPDX-License-Identifier: Apache-2.0
 * © 2025 BioInception PVT LTD.
 */
package org.openscience.smsd.filters;

import java.util.List;
import java.util.Map;

import org.openscience.cdk.exception.CDKException;
import org.openscience.smsd.AtomAtomMapping;

/**
 * A filter on SMSD results.
 *
 * @param <T>
 *
 * @author Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>
 * @author maclean
 *
 */
public interface IChemicalFilter<T> {

    /**
     * Calculates a score for each MCS, and sorts the results on that score,
     * returning the best.
     *
     * @param allAtomMCS
     * @param selectionMap
     * @return
     * @throws CDKException
     */
    public T sortResults(
            Map<Integer, AtomAtomMapping> allAtomMCS,
            Map<Integer, T> selectionMap) throws CDKException;

    public List<T> getScores();

    public void clearScores();

    public void addScore(int counter, T value);

    public void fillMap(Map<Integer, T> map);
}
