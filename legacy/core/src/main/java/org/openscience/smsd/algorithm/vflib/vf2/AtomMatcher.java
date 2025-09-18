/*
 * SPDX-License-Identifier: Apache-2.0
 * © 2025 BioInception PVT LTD.
 */
package org.openscience.smsd.algorithm.vflib.vf2;

import org.openscience.cdk.interfaces.IAtom;

/**
 * Interface for the AtomMatcher (atoms) in graph.
 *
 * 
 * 
 * @author Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>
 */
public interface AtomMatcher {

    boolean matches(IAtom queryAtom, IAtom targetAtom);
}
