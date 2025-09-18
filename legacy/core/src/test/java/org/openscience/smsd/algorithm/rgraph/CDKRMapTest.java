/*
 * SPDX-License-Identifier: Apache-2.0
 * © 2025 BioInception PVT LTD.
 */
package org.openscience.smsd.algorithm.rgraph;

import org.junit.Assert;
import org.junit.Test;

/**
 * test-smsd
 *
 * @author Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>
 * java1.8+
 */
public class CDKRMapTest {

    @Test
    public void testRMap_int_int() {
        CDKRMap node = new CDKRMap(1, 2);
        Assert.assertNotNull(node);
        Assert.assertEquals(1, node.getId1());
        Assert.assertEquals(2, node.getId2());
    }
}
