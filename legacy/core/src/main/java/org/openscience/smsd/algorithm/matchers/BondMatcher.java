/*
 * SPDX-License-Identifier: Apache-2.0
 * © 2025 BioInception PVT LTD.
 */
package org.openscience.smsd.algorithm.matchers;

import org.openscience.cdk.interfaces.IBond;

/**
 * Interface for the BondMatcher (bonds) in graph.
 *
 * 
 * 
 * @author Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>
 */
public interface BondMatcher {

    boolean matches(IBond bond);
}
