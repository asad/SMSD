
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
package org.openscience.smsd;

import org.junit.jupiter.api.AfterEach;
import org.junit.jupiter.api.AfterAll;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.BeforeAll;
import org.junit.platform.suite.api.SelectClasses;
import org.junit.platform.suite.api.Suite;


/**
 * @author Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>
 *
 *  test-smsd
 *  java1.8+
 */
@Suite
@SelectClasses({org.openscience.smsd.IsomorphismTest.class,
    org.openscience.smsd.helper.HelperSuite.class,
    //        org.openscience.smsd.interfaces.InterfacesSuite.class,
    org.openscience.smsd.filters.FiltersSuite.class,
    org.openscience.smsd.SubstructureTest.class,
    org.openscience.smsd.algorithm.AlgorithmSuite.class,
    org.openscience.smsd.tools.ToolsSuite.class})
public class SmsdSuite {

    @BeforeAll
    public static void setUpClass() throws Exception {
    }

    @AfterAll
    public static void tearDownClass() throws Exception {
    }

    @BeforeEach 

    public void setUp() throws Exception {
    }

    @AfterEach 

    public void tearDown() throws Exception {
    }
}





