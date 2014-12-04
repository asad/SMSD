
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
package org.openscience.smsd.algorithm.mcsplus;

import java.util.logging.Level;
import java.util.logging.Logger;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.normalize.SMSDNormalizer;
import org.openscience.cdk.smiles.SmilesParser;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import static org.junit.Assert.*;
import org.openscience.smsd.tools.ExtAtomContainerManipulator;

/**
 * @author Syed Asad Rahman <asad@ebi.ac.uk>
 *
 * @cdk.module test-smsd
 * @cdk.require java1.6+
 */
public class MCSPlusHandlerTest {

    public MCSPlusHandlerTest() {
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
     * Test of searchMCS method, of class MCSPlusHandler.
     */
    @Test
    public void testSearchMCS() {
        try {
            ////////System.out.println("1");
            SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
            IAtomContainer target;
            target = sp.parseSmiles("C\\C=C/Nc1cccc(c1)N(O)\\C=C\\C\\C=C\\C=C/C");
            IAtomContainer queryac = null;
            queryac = sp.parseSmiles("Nc1ccccc1");
            MCSPlusHandler smsd1 = new MCSPlusHandler(queryac, target, true, false, false);
            assertNotNull(smsd1.getFirstAtomMapping());
        } catch (InvalidSmilesException ex) {
            Logger.getLogger(MCSPlusHandlerTest.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

    /**
     * Test of set method, of class MCSPlusHandler.
     *
     * @throws InvalidSmilesException
     */
    @Test
    public void testSet_IAtomContainer_IAtomContainer() throws InvalidSmilesException {
        ////////System.out.println("2");
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer target = sp.parseSmiles("C\\C=C/Nc1cccc(c1)N(O)\\C=C\\C\\C=C\\C=C/C");
        IAtomContainer queryac = sp.parseSmiles("Nc1ccccc1");

        MCSPlusHandler smsd1 = new MCSPlusHandler(queryac, target, true, false, false);
        assertNotNull(smsd1.getFirstAtomMapping());

    }

    /**
     * Test of set method, of class MCSPlusHandler.
     *
     * @throws Exception
     */
    @Test
    public void testSet_IMolecule_IMolecule() throws Exception {
        ////////System.out.println("3");
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer target = sp.parseSmiles("C\\C=C/Nc1cccc(c1)N(O)\\C=C\\C\\C=C\\C=C/C");
        IAtomContainer queryac = sp.parseSmiles("Nc1ccccc1");

        MCSPlusHandler smsd1 = new MCSPlusHandler(queryac, target, true, false, false);
        assertNotNull(smsd1.getFirstAtomMapping());
    }

    /**
     * Test of set method, of class MCSPlusHandler.
     *
     * @throws InvalidSmilesException
     */
    @Test
    public void testSet_AtomContainer_AtomContainer() throws InvalidSmilesException {
        ////////System.out.println("4");
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());

        IAtomContainer target = sp.parseSmiles("C\\C=C/Nc1cccc(c1)N(O)\\C=C\\C\\C=C\\C=C/C");
        IAtomContainer queryac = sp.parseSmiles("Nc1ccccc1");
        MCSPlusHandler instance = new MCSPlusHandler(queryac, target, true, false, false);
        assertNotNull(instance.getFirstAtomMapping());
    }

    /**
     * Test of getAllAtomMapping method, of class MCSPlusHandler.
     *
     * @throws InvalidSmilesException
     */
    @Test
    public void testGetAllAtomMapping() throws InvalidSmilesException {
        ////////System.out.println("5");
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer query = sp.parseSmiles("Nc1cccc(NO)c1");
        IAtomContainer target = sp.parseSmiles("Nc1ccccc1");

        MCSPlusHandler smsd1 = new MCSPlusHandler(query, target, true, false, false);
        assertNotNull(smsd1.getFirstAtomMapping());
        assertEquals(4, smsd1.getAllAtomMapping().size());
    }

    /**
     * Test of getAllAtomMapping method, of class MCSPlusHandler.
     *
     * @throws InvalidSmilesException
     */
    @Test
    public void testgetAllAtomMapping() throws InvalidSmilesException, CDKException {
        ////////System.out.println("6");
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer query = sp.parseSmiles("Nc1cccc(NO)c1");
        IAtomContainer target = sp.parseSmiles("Nc1ccccc1");
        MCSPlusHandler comparison = new MCSPlusHandler(query, target, true, true, false);
        assertNotNull(comparison.getFirstAtomMapping());
        assertEquals(4, comparison.getAllAtomMapping().size());
    }

    /**
     * Test of getFirstAtomMapping method, of class MCSPlusHandler.
     *
     * @throws InvalidSmilesException
     */
    @Test
    public void testGetFirstAtomMapping() throws InvalidSmilesException {
        ////////System.out.println("7");
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer query = sp.parseSmiles("Nc1cccc(NO)c1");
        IAtomContainer target = sp.parseSmiles("Nc1ccccc1");

        MCSPlusHandler smsd1 = new MCSPlusHandler(query, target, true, false, false);
        assertNotNull(smsd1.getFirstAtomMapping());

        assertEquals(7, smsd1.getFirstAtomMapping().getCount());
    }

    /**
     * Test of getFirstAtomMapping method, of class MCSPlusHandler.
     *
     * @throws InvalidSmilesException
     * @throws CDKException
     */
    @Test
    public void testgetFirstAtomMapping() throws InvalidSmilesException, CDKException {
        ////////System.out.println("8");
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer query = sp.parseSmiles("Nc1cccc(NO)c1");
        IAtomContainer target = sp.parseSmiles("Nc1ccccc1");

        SMSDNormalizer.percieveAtomTypesAndConfigureAtoms(target);
        ExtAtomContainerManipulator.aromatizeDayLight(target);

        SMSDNormalizer.percieveAtomTypesAndConfigureAtoms(query);
        ExtAtomContainerManipulator.aromatizeDayLight(query);

        MCSPlusHandler smsd1 = new MCSPlusHandler(query, target, true, false, false);
        assertNotNull(smsd1.getFirstAtomMapping());
        assertEquals(7, smsd1.getFirstAtomMapping().getCount());
    }
}
