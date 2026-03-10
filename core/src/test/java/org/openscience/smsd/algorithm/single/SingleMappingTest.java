
/* Copyright (C) 2009-2014 Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>
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
package org.openscience.smsd.algorithm.single;

import org.junit.jupiter.api.AfterEach;
import org.junit.jupiter.api.AfterAll;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;
import org.openscience.cdk.Atom;
import org.openscience.cdk.AtomContainer;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;

/**
 * @author Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>
 *
 *  test-smsd
 *  java1.8+
 */
public class SingleMappingTest {

    public SingleMappingTest() {
    }

    @BeforeAll
    public static void setUpClass() throws Exception {
    }

    @AfterAll
    public static void tearDownClass() throws Exception {
    }

    @BeforeEach 

    public void setUp() {
    }

    @AfterEach 

    public void tearDown() {
    }

    /**
     * Test of getOverLaps method, of class SingleMapping.
     */
    @Test
    public void testGetOverLaps() throws CDKException {
        ////////System.out.println("getOverLaps");
        IAtom atomSource = new Atom("R");
        IAtom atomTarget = new Atom("R");
        IAtomContainer source = new AtomContainer();
        source.addAtom(atomSource);
        IAtomContainer target = new AtomContainer();
        target.addAtom(atomTarget);
        SingleMapping instance = new SingleMapping();
        Assertions.assertNotNull(instance.getOverLaps(source, target));
    }
}




