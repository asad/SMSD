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
public class CDKRNodeTest {

    @Test
    public void testRNode_int_int() {
        CDKRNode node = new CDKRNode(1, 2);
        Assert.assertNotNull(node);
        Assert.assertNotNull(node.getExtension());
        Assert.assertNotNull(node.getForbidden());
    }
}
