/*
 * SPDX-License-Identifier: Apache-2.0
 * © 2025 BioInception PVT LTD.
 */
package org.openscience.smsd.algorithm.matchers;

import org.openscience.cdk.interfaces.IAtom;

/**
 * Interface for the AtomMatcher (atoms) in graph.
 *
 * 
 * 
 * @author Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>
 */
public interface AtomMatcher {

    boolean matches(IAtom atom);

    IAtom getQueryAtom();
}
