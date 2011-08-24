/* Copyright (C) 2009-2011  Syed Asad Rahman <asad@ebi.ac.uk>
 *
 * Contact: cdk-devel@lists.sourceforge.net
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 * All we ask is that proper credit is given for our work, which includes
 * - but is not limited to - adding the above copyright notice to the beginning
 * of your source code files, and to any copyright notice that you may distribute
 * with programs based on this work.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received commonAtomList copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA.
 */
package org.openscience.smsd.algorithm.vflib;

import java.util.List;
import java.util.Map;
import junit.framework.Assert;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.isomorphism.matchers.IQueryAtomContainer;
import org.openscience.smsd.AtomAtomMapping;
import org.openscience.smsd.tools.TimeManager;

/**
 * This class should be used to find MCS between source
 * graph and target graph.
 *
 * First the algorithm runs VF lib {@link org.openscience.cdk.smsd.algorithm.vflib.map.VFMCSMapper}
 * and reports MCS between
 * run source and target graphs. Then these solutions are extended
 * using McGregor {@link org.openscience.cdk.smsd.algorithm.mcgregor.McGregor}
 * algorithm where ever required.
 *
 * @cdk.module smsd-test
 * @cdk.githash
 * @author Syed Asad Rahman <asad@ebi.ac.uk>
 */
public class VF2SubTest {

    public VF2SubTest() {
    }

    @BeforeClass
    public static void setUpClass() throws Exception {
    }

    @AfterClass
    public static void tearDownClass() throws Exception {
    }

    @Before
    public void setUp() {
    }

    @After
    public void tearDown() {
    }

    /**
     * Test of setTimeManager method, of class VF2Sub.
     */
    @Test
    public void testSetTimeManager() {
        System.out.println("setTimeManager");
        TimeManager aTimeManager = null;
        VF2Sub instance = new VF2Sub();
        instance.setTimeManager(aTimeManager);
        // TODO review the generated test code and remove the default call to fail.
        Assert.fail("The test case is a prototype.");
    }

    /**
     * Test of isSubgraph method, of class VF2Sub.
     */
    @Test
    public void testIsSubgraph() {
        System.out.println("isSubgraph");
        boolean shouldMatchBonds = false;
        boolean shouldMatchRings = false;
        VF2Sub instance = new VF2Sub();
        boolean expResult = false;
        boolean result = instance.isSubgraph(shouldMatchBonds, shouldMatchRings);
        Assert.assertEquals(expResult, result);
        // TODO review the generated test code and remove the default call to fail.
        Assert.fail("The test case is a prototype.");
    }

    /**
     * Test of set method, of class VF2Sub.
     */
    @Test
    public void testSet_IAtomContainer_IAtomContainer() {
        System.out.println("set");
        IAtomContainer source = null;
        IAtomContainer target = null;
        VF2Sub instance = new VF2Sub();
        instance.set(source, target);
        // TODO review the generated test code and remove the default call to fail.
        Assert.fail("The test case is a prototype.");
    }

    /**
     * Test of set method, of class VF2Sub.
     */
    @Test
    public void testSet_IQueryAtomContainer_IAtomContainer() {
        System.out.println("set");
        IQueryAtomContainer source = null;
        IAtomContainer target = null;
        VF2Sub instance = new VF2Sub();
        instance.set(source, target);
        // TODO review the generated test code and remove the default call to fail.
        Assert.fail("The test case is a prototype.");
    }

    /**
     * Test of getAllMapping method, of class VF2Sub.
     */
    @Test
    public void testGetAllMapping() {
        System.out.println("getAllMapping");
        VF2Sub instance = new VF2Sub();
        List expResult = null;
        List result = instance.getAllMapping();
        Assert.assertEquals(expResult, result);
        // TODO review the generated test code and remove the default call to fail.
        Assert.fail("The test case is a prototype.");
    }

    /**
     * Test of getFirstMapping method, of class VF2Sub.
     */
    @Test
    public void testGetFirstMapping() {
        System.out.println("getFirstMapping");
        VF2Sub instance = new VF2Sub();
        Map expResult = null;
        Map result = instance.getFirstMapping();
        Assert.assertEquals(expResult, result);
        // TODO review the generated test code and remove the default call to fail.
        Assert.fail("The test case is a prototype.");
    }

    /**
     * Test of getAllAtomMapping method, of class VF2Sub.
     */
    @Test
    public void testGetAllAtomMapping() {
        System.out.println("getAllAtomMapping");
        VF2Sub instance = new VF2Sub();
        List expResult = null;
        List result = instance.getAllAtomMapping();
        Assert.assertEquals(expResult, result);
        // TODO review the generated test code and remove the default call to fail.
        Assert.fail("The test case is a prototype.");
    }

    /**
     * Test of getFirstAtomMapping method, of class VF2Sub.
     */
    @Test
    public void testGetFirstAtomMapping() {
        System.out.println("getFirstAtomMapping");
        VF2Sub instance = new VF2Sub();
        AtomAtomMapping expResult = null;
        AtomAtomMapping result = instance.getFirstAtomMapping();
        Assert.assertEquals(expResult, result);
        // TODO review the generated test code and remove the default call to fail.
        Assert.fail("The test case is a prototype.");
    }

    /**
     * Test of isBondMatchFlag method, of class VF2Sub.
     */
    @Test
    public void testIsBondMatchFlag() {
        System.out.println("isBondMatchFlag");
        VF2Sub instance = new VF2Sub();
        boolean expResult = false;
        boolean result = instance.isBondMatchFlag();
        Assert.assertEquals(expResult, result);
        // TODO review the generated test code and remove the default call to fail.
        Assert.fail("The test case is a prototype.");
    }

    /**
     * Test of setBondMatchFlag method, of class VF2Sub.
     */
    @Test
    public void testSetBondMatchFlag() {
        System.out.println("setBondMatchFlag");
        boolean shouldMatchBonds = false;
        VF2Sub instance = new VF2Sub();
        instance.setBondMatchFlag(shouldMatchBonds);
        // TODO review the generated test code and remove the default call to fail.
        Assert.fail("The test case is a prototype.");
    }

    /**
     * Test of isTimeOut method, of class VF2Sub.
     */
    @Test
    public void testIsTimeOut() {
        System.out.println("isTimeOut");
        VF2Sub instance = new VF2Sub();
        boolean expResult = false;
        boolean result = instance.isTimeOut();
        Assert.assertEquals(expResult, result);
        // TODO review the generated test code and remove the default call to fail.
        Assert.fail("The test case is a prototype.");
    }

    /**
     * Test of isMatchRings method, of class VF2Sub.
     */
    @Test
    public void testIsMatchRings() {
        System.out.println("isMatchRings");
        VF2Sub instance = new VF2Sub();
        boolean expResult = false;
        boolean result = instance.isMatchRings();
        Assert.assertEquals(expResult, result);
        // TODO review the generated test code and remove the default call to fail.
        Assert.fail("The test case is a prototype.");
    }

    /**
     * Test of setMatchRings method, of class VF2Sub.
     */
    @Test
    public void testSetMatchRings() {
        System.out.println("setMatchRings");
        boolean shouldMatchRings = false;
        VF2Sub instance = new VF2Sub();
        instance.setMatchRings(shouldMatchRings);
        // TODO review the generated test code and remove the default call to fail.
        Assert.fail("The test case is a prototype.");
    }
}
