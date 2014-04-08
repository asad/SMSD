
/* Copyright (C) 2009-2014 Syed Asad Rahman <asad@ebi.ac.uk>
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
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software
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
import org.openscience.smsd.AtomAtomMapping;

/**
 * @author Syed Asad Rahman <asad@ebi.ac.uk>
 *
 * @cdk.module test-smsd @cdk.require java1.6+
 */
public class VF2Test {
    
    public VF2Test() {
    }
    
    @BeforeClass
    public static void setUpClass() {
    }
    
    @AfterClass
    public static void tearDownClass() {
    }
    
    @Before
    public void setUp() {
    }
    
    @After
    public void tearDown() {
    }

    /**
     * Test of getAllAtomMapping method, of class VF2.
     */
    @Test
    public void testGetAllAtomMapping() {
        //////System.out.println("getAllAtomMapping");
        VF2 instance = null;
        List<AtomAtomMapping> expResult = null;
        List<AtomAtomMapping> result = instance.getAllAtomMapping();
        assertEquals(expResult, result);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of getFirstAtomMapping method, of class VF2.
     */
    @Test
    public void testGetFirstAtomMapping() {
        //////System.out.println("getFirstAtomMapping");
        VF2 instance = null;
        AtomAtomMapping expResult = null;
        AtomAtomMapping result = instance.getFirstAtomMapping();
        assertEquals(expResult, result);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of isSubgraph method, of class VF2.
     */
    @Test
    public void testIsSubgraph() {
        //////System.out.println("isSubgraph");
        VF2 instance = null;
        boolean expResult = false;
        boolean result = instance.isSubgraph();
        assertEquals(expResult, result);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }
    
}
