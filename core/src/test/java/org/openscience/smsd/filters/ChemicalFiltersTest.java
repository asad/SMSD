
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
package org.openscience.smsd.filters;

import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.smsd.Isomorphism;
import org.openscience.smsd.interfaces.Algorithm;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 * @author Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>
 *
 *  test-smsd
 *  java1.8+
 */
public class ChemicalFiltersTest {

    public ChemicalFiltersTest() {
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
     * Test of sortResultsByStereoAndBondMatch method, of class ChemicalFilters.
     *
     * @throws Exception
     */
    @Test
    public void testSortResultsByStereoAndBondMatch() throws Exception {
        //////System.out.println("sortResultsByStereoAndBondMatch");
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer target = sp.parseSmiles("C\\C=C/OCC=C");
        IAtomContainer queryac = sp.parseSmiles("CCCOCC(C)=C");

        Isomorphism smsd = new Isomorphism(queryac, target, Algorithm.DEFAULT, false, false, false);
        smsd.setChemFilters(false, false, false);
        assertEquals(4, smsd.getAllAtomMapping().size());

        Isomorphism smsd1 = new Isomorphism(queryac, target, Algorithm.DEFAULT, false, false, false);
        smsd1.setChemFilters(true, false, false);
        assertEquals(1, smsd1.getAllAtomMapping().size());
    }

    /**
     * Test of sortResultsByFragments method, of class ChemicalFilters.
     *
     * @throws InvalidSmilesException
     * @throws CDKException
     */
    @Test
    public void testSortResultsByFragments() throws InvalidSmilesException, CDKException {
        //////System.out.println("sortResultsByFragments");
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer target = sp.parseSmiles("C\\C=C/Nc1cccc(c1)N(O)\\C=C\\C\\C=C\\C=C/C");
        IAtomContainer queryac = sp.parseSmiles("Nc1ccccc1");

        Isomorphism smsd = new Isomorphism(queryac, target, Algorithm.CDKMCS, false, false, false);
        smsd.setChemFilters(false, false, false);
        assertEquals(4, smsd.getAllAtomMapping().size());

        Isomorphism smsd1 = new Isomorphism(queryac, target, Algorithm.CDKMCS, false, false, false);
        smsd1.setChemFilters(false, true, false);
        assertEquals(2, smsd1.getAllAtomMapping().size());
    }

    /**
     * Test of sortResultsByEnergies method, of class ChemicalFilters.
     *
     * @throws Exception
     */
    @Test
    public void testSortResultsByEnergies() throws Exception {
        //////System.out.println("sortResultsByEnergies");
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer target = sp.parseSmiles("C\\C=C/Nc1cccc(c1)N(O)\\C=C\\C\\C=C\\C=C/C");
        IAtomContainer queryac = sp.parseSmiles("Nc1ccccc1");

        Isomorphism smsd = new Isomorphism(queryac, target, Algorithm.DEFAULT, true, false, false);
        smsd.setChemFilters(false, false, false);
        assertEquals(4, smsd.getAllAtomMapping().size());

        Isomorphism smsd1 = new Isomorphism(queryac, target, Algorithm.DEFAULT, true, false, false);
        smsd1.setChemFilters(false, false, true);
        assertEquals(2, smsd1.getAllAtomMapping().size());
    }

    /**
     * Test of getSortedEnergy method, of class ChemicalFilters.
     *
     * @throws InvalidSmilesException
     * @throws CDKException
     */
    @Test
    public void testGetSortedEnergy() throws InvalidSmilesException, CDKException {
        //////System.out.println("getSortedEnergy");
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer target = sp.parseSmiles("C\\C=C/Nc1cccc(c1)N(O)\\C=C\\C\\C=C\\C=C/C");
        IAtomContainer queryac = sp.parseSmiles("Nc1ccccc1");

        Isomorphism smsd1 = new Isomorphism(queryac, target, Algorithm.DEFAULT, true, false, false);
        smsd1.setChemFilters(false, false, true);
        Double score = 610.0;
        assertEquals(score, smsd1.getEnergyScore(0));
    }

    /**
     * Test of getSortedFragment method, of class ChemicalFilters.
     *
     * @throws InvalidSmilesException
     * @throws CDKException
     */
    @Test
    public void testGetSortedFragment() throws InvalidSmilesException, CDKException {
        //////System.out.println("getSortedFragment");
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer target = sp.parseSmiles("C\\C=C/Nc1cccc(c1)N(O)\\C=C\\C\\C=C\\C=C/C");
        IAtomContainer queryac = sp.parseSmiles("Nc1ccccc1");

        Isomorphism smsd1 = new Isomorphism(queryac, target, Algorithm.DEFAULT, true, false, false);
        smsd1.setChemFilters(false, true, false);
        Integer score = 2;
        assertEquals(score, smsd1.getFragmentSize(0));
    }

    /**
     * Test of getStereoMatches method, of class ChemicalFilters.
     *
     * @throws InvalidSmilesException
     * @throws CDKException
     */
    @Test
    public void testGetStereoMatches() throws InvalidSmilesException, CDKException {
        //////System.out.println("getStereoMatches");
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer target = sp.parseSmiles("C\\C=C/OCC=C");
        IAtomContainer queryac = sp.parseSmiles("CCCOCC(C)=C");

        Isomorphism smsd1 = new Isomorphism(queryac, target, Algorithm.DEFAULT, false, false, false);
        smsd1.setChemFilters(true, false, false);
        Integer score = 77;//1048 to 77 as per new score
        assertEquals(score, smsd1.getStereoScore(0));
    }
}
