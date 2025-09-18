/*
 * SPDX-License-Identifier: Apache-2.0
 * © 2025 BioInception PVT LTD.
 */
package org.openscience.smsd.algorithm.rgraph;

import org.junit.Assert;
import org.junit.Test;

/**
 *  test-smsd
 * @author     Syed Asad Rahman
 *  java1.8+
 */
public class CDKRGraphTest {

    @Test
    public void testRGraph() {
        CDKRGraph graph = new CDKRGraph();
        Assert.assertNotNull(graph);
    }
}
