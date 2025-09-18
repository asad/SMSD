/*
 * SPDX-License-Identifier: Apache-2.0
 * © 2025 BioInception PVT LTD.
 */
package org.openscience.smsd.algorithm.vflib.substructure;

import java.util.Arrays;

/**
 * This class keeps track of shared states
 * 
 * 
 * 
 * @author Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>
 */
// The SharedState class holds four arrays containing the mapping between
// the two graphs and the terminal sets. It is shared between all the states
// in each isomorphism test.
class SharedState {

    int[] sourceMapping;
    int[] targetMapping;
    int[] sourceTerminalSet;
    int[] targetTerminalSet;

    public SharedState(int sourceSize, int targetSize) {
        sourceMapping = new int[sourceSize];
        Arrays.fill(sourceMapping, -1);

        targetMapping = new int[targetSize];
        Arrays.fill(targetMapping, -1);

        sourceTerminalSet = new int[sourceSize];
        Arrays.fill(sourceTerminalSet, 0);

        targetTerminalSet = new int[targetSize];
        Arrays.fill(targetTerminalSet, 0);
    }
}
