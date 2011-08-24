/* Copyright (C) 2010  Egon Willighagen <egonw@users.sf.net>
 *
 * Contact: cdk-devel@lists.sourceforge.net
 *
 * This program is free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation; either version 2.1 of the License, or (at your option)
 * any later version. All we ask is that proper credit is given for our work,
 * which includes - but is not limited to - adding the above copyright notice to
 * the beginning of your source code files, and to any copyright notice that you
 * may distribute with programs based on this work.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
 * details.
 *
 * You should have received rAtomCount copy of the GNU Lesser General Public
 * License along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA.
 */
package org.openscience.smsd.algorithm.vflib.substructure;

import java.util.List;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import static org.junit.Assert.*;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.isomorphism.matchers.IQueryAtomContainer;
import org.openscience.smsd.AtomAtomMapping;

/**
 * Interface class for reporting only substructure searches.
 * @cdk.module smsd
 * @cdk.githash
 * @author Syed Asad Rahman <asad@ebi.ac.uk>
 */
public class VF2Test {

    public VF2Test() {
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
     * Test of isomorphism method, of class VF2.
     */
    @Test
    public void testIsomorphism_4args() {
        System.out.println("isomorphism");
        IAtomContainer a = null;
        IAtomContainer b = null;
        boolean shouldMatchBonds = false;
        boolean shouldMatchRings = false;
        VF2 instance = new VF2();
        AtomAtomMapping expResult = null;
        AtomAtomMapping result = instance.isomorphism(a, b, shouldMatchBonds, shouldMatchRings);
        assertEquals(expResult, result);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of isomorphism method, of class VF2.
     */
    @Test
    public void testIsomorphism_IQueryAtomContainer_IAtomContainer() {
        System.out.println("isomorphism");
        IQueryAtomContainer a = null;
        IAtomContainer b = null;
        VF2 instance = new VF2();
        AtomAtomMapping expResult = null;
        AtomAtomMapping result = instance.isomorphism(a, b);
        assertEquals(expResult, result);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of isomorphisms method, of class VF2.
     */
    @Test
    public void testIsomorphisms() {
        System.out.println("isomorphisms");
        IAtomContainer a = null;
        IAtomContainer b = null;
        boolean shouldMatchBonds = false;
        boolean shouldMatchRings = false;
        VF2 instance = new VF2();
        List expResult = null;
        List result = instance.isomorphisms(a, b, shouldMatchBonds, shouldMatchRings);
        assertEquals(expResult, result);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }
}
