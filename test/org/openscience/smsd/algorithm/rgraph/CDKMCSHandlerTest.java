
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
package org.openscience.smsd.algorithm.rgraph;

import java.util.logging.Level;
import java.util.logging.Logger;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.smiles.SmilesParser;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 * @author Syed Asad Rahman <asad@ebi.ac.uk>
 *
 * @cdk.module test-smsd @cdk.require java1.6+
 */
public class CDKMCSHandlerTest {

    public CDKMCSHandlerTest() {
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
     * Test of searchMCS method, of class CDKMCSHandler.
     */
    @Test
    public void testSearchMCS() {
        try {
            ////////System.out.println("searchMCS");
            SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
            IAtomContainer target = null;
            target = sp.parseSmiles("C\\C=C/Nc1cccc(c1)N(O)\\C=C\\C\\C=C\\C=C/C");
            IAtomContainer query = null;
            query = sp.parseSmiles("Nc1ccccc1");
            CDKMCSHandler smsd1 = new CDKMCSHandler(query, target, true, false,false);
            assertNotNull(smsd1.getFirstAtomMapping());
        } catch (InvalidSmilesException ex) {
            Logger.getLogger(CDKMCSHandlerTest.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

    /**
     * Test of set method, of class CDKMCSHandler.
     *
     * @throws Exception
     */
    @Test
    public void testSet_IMolecule_IMolecule() throws Exception {
        ////////System.out.println("set");
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer target = sp.parseSmiles("C\\C=C/Nc1cccc(c1)N(O)\\C=C\\C\\C=C\\C=C/C");
        IAtomContainer query = sp.parseSmiles("Nc1ccccc1");

        CDKMCSHandler smsd1 = new CDKMCSHandler(query, target, true, false,false);
        assertNotNull(smsd1.getFirstAtomMapping());
    }

    /**
     * Test of set method, of class CDKMCSHandler.
     *
     * @throws InvalidSmilesException
     */
    @Test
    public void testSet_AtomContainer_AtomContainer() throws InvalidSmilesException {
        ////////System.out.println("set");
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());

        IAtomContainer target = sp.parseSmiles("C\\C=C/Nc1cccc(c1)N(O)\\C=C\\C\\C=C\\C=C/C");
        IAtomContainer query = sp.parseSmiles("Nc1ccccc1");
        CDKMCSHandler instance = new CDKMCSHandler(query, target, true, false,false);
        assertNotNull(instance.getFirstAtomMapping());
    }

    /**
     * Test of getAllAtomMapping method, of class CDKMCSHandler.
     *
     * @throws InvalidSmilesException
     */
    @Test
    public void testGetAllAtomMapping() throws InvalidSmilesException {
        ////////System.out.println("getAllAtomMapping");
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer target = sp.parseSmiles("C\\C=C/Nc1cccc(c1)N(O)\\C=C\\C\\C=C\\C=C/C");
        IAtomContainer query = sp.parseSmiles("Nc1ccccc1");

        CDKMCSHandler smsd1 = new CDKMCSHandler(query, target, true, false,false);
        assertNotNull(smsd1.getFirstAtomMapping());
        assertEquals(4, smsd1.getAllAtomMapping().size());
    }

    /**
     * Test of getAllAtomMapping method, of class CDKMCSHandler.
     *
     * @throws InvalidSmilesException
     */
    @Test
    public void testgetAllAtomMapping() throws InvalidSmilesException {
        ////////System.out.println("getAllAtomMapping");
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer target = sp.parseSmiles("C\\C=C/Nc1cccc(c1)N(O)\\C=C\\C\\C=C\\C=C/C");
        IAtomContainer query = sp.parseSmiles("Nc1ccccc1");

        CDKMCSHandler smsd1 = new CDKMCSHandler(query, target, true, false,false);
        assertNotNull(smsd1.getFirstAtomMapping());
        assertEquals(4, smsd1.getAllAtomMapping().size());
    }

    /**
     * Test of getFirstAtomMapping method, of class CDKMCSHandler.
     *
     * @throws InvalidSmilesException
     */
    @Test
    public void testGetFirstAtomMapping() throws InvalidSmilesException {
        ////////System.out.println("getFirstAtomMapping");
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer target = sp.parseSmiles("C\\C=C/Nc1cccc(c1)N(O)\\C=C\\C\\C=C\\C=C/C");
        IAtomContainer query = sp.parseSmiles("Nc1ccccc1");

        CDKMCSHandler smsd1 = new CDKMCSHandler(query, target, true, false,false);
        assertNotNull(smsd1.getFirstAtomMapping());
        assertEquals(7, smsd1.getFirstAtomMapping().getCount());
    }

    /**
     * Test of getFirstAtomMapping method, of class CDKMCSHandler.
     *
     * @throws InvalidSmilesException
     */
    @Test
    public void testgetFirstAtomMapping() throws InvalidSmilesException {
        ////////System.out.println("getFirstAtomMapping");
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer target = sp.parseSmiles("C\\C=C/Nc1cccc(c1)N(O)\\C=C\\C\\C=C\\C=C/C");
        IAtomContainer query = sp.parseSmiles("Nc1ccccc1");

        CDKMCSHandler smsd1 = new CDKMCSHandler(query, target, true, false,false);
        assertNotNull(smsd1.getFirstAtomMapping());
        assertEquals(7, smsd1.getFirstAtomMapping().getCount());
    }
}
