
/* 
 * Copyright (C) 2009-2014 Syed Asad Rahman <asad@ebi.ac.uk>
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

import java.io.IOException;
import java.util.logging.Level;
import java.util.logging.Logger;
import junit.framework.Assert;
import org.junit.*;
import org.openscience.cdk.AtomContainer;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.isomorphism.matchers.IQueryAtomContainer;
import org.openscience.cdk.isomorphism.matchers.QueryAtomContainer;
import org.openscience.cdk.isomorphism.matchers.QueryAtomContainerCreator;
import org.openscience.cdk.normalize.SMSDNormalizer;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.smsd.algorithm.mcsplus.MCSPlusHandlerTest;
import org.openscience.smsd.interfaces.Algorithm;
import org.openscience.smsd.tools.ExtAtomContainerManipulator;

/**
 * Unit testing for the {@link Isomorphism} class.
 *
 * @author Syed Asad Rahman <asad@ebi.ac.uk>
 *
 * @cdk.module test-smsd
 * @cdk.require java1.6+
 */
public class IsomorphismTest {

    /**
     * Tests if the CDKMCS can be instantiated without throwing exceptions.
     */
    @Test
    public void IsomorphismTest() {
        Assert.assertNotNull(
                new Isomorphism(new AtomContainer(), new AtomContainer(),
                        Algorithm.DEFAULT,
                        true, false, true));
        Assert.assertNotNull(
                new Isomorphism(new AtomContainer(), new AtomContainer(),
                        Algorithm.DEFAULT,
                        false, false, false));
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
     * Test of init method, of class Isomorphism.
     *
     * @throws CDKException
     */
    @Test
    public void testInit_3args_1() throws CDKException {
        ////System.out.println("1");
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer target = sp.parseSmiles("C\\C=C/OCC=C");
        IAtomContainer query = sp.parseSmiles("CCCOCC(C)=C");

        Isomorphism smsd1 = new Isomorphism(query, target, Algorithm.DEFAULT, false, false, false);
        smsd1.setChemFilters(true, false, false);
        Assert.assertNotNull(smsd1.getQuery());
        Assert.assertNotNull(smsd1.getTarget());
    }

    /**
     * Test of init method, of class Isomorphism.
     *
     * @throws CDKException
     */
    @Test
    public void testInit_3args_2() throws CDKException {
        ////System.out.println("2");
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer target = sp.parseSmiles("C\\C=C/OCC=C");
        IAtomContainer query = sp.parseSmiles("CCCOCC(C)=C");

        Isomorphism smsd1 = new Isomorphism(query, target, Algorithm.DEFAULT, false, false, false);
        smsd1.setChemFilters(true, false, false);
        Assert.assertNotNull(smsd1.getQuery());
        Assert.assertNotNull(smsd1.getTarget());
    }

    /**
     * Test of searchMCS method, of class Isomorphism.
     *
     * @throws CDKException
     */
    @Test
    public void testSearchMCS() throws CDKException {
        try {
            ////System.out.println("3");
            SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
            IAtomContainer target;
            target = sp.parseSmiles("C\\C=C/Nc1cccc(c1)N(O)\\C=C\\C\\C=C\\C=C/C");
            //("NC1=CC=C(N)C=C1");
            IAtomContainer query = sp.parseSmiles("Nc1ccccc1");
            Isomorphism smsd1 = new Isomorphism(query, target, Algorithm.CDKMCS, true, false, false);
            smsd1.setChemFilters(false, true, true);
            Assert.assertEquals(7, smsd1.getFirstAtomMapping().getCount());
            Assert.assertEquals(2, smsd1.getAllAtomMapping().size());
            Assert.assertNotNull(smsd1.getFirstAtomMapping());
        } catch (InvalidSmilesException ex) {
            Logger.getLogger(MCSPlusHandlerTest.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

    /**
     * Test of set method, of class Isomorphism.
     *
     * @throws CDKException
     */
    @Test
    public void testSet_IAtomContainer_IAtomContainer() throws CDKException {
        ////System.out.println("4");
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer target = sp.parseSmiles("C\\C=C/Nc1cccc(c1)N(O)\\C=C\\C\\C=C\\C=C/C");
        IAtomContainer query = sp.parseSmiles("Nc1ccccc1");

        Isomorphism smsd1 = new Isomorphism(query, target, Algorithm.DEFAULT, true, false, false);
        smsd1.setChemFilters(true, true, true);
        Assert.assertNotNull(smsd1.getFirstAtomMapping());

    }

    /**
     * Test of set method, of class Isomorphism.
     *
     * @throws Exception
     */
    @Test
    public void testSet_IMolecule_IMolecule() throws Exception {
        ////System.out.println("5");
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer target = sp.parseSmiles("C\\C=C/Nc1cccc(c1)N(O)\\C=C\\C\\C=C\\C=C/C");
        IAtomContainer query = sp.parseSmiles("Nc1ccccc1");

        Isomorphism smsd1 = new Isomorphism(query, target, Algorithm.DEFAULT, true, false, false);
        smsd1.setChemFilters(true, true, true);
        Assert.assertNotNull(smsd1.getFirstAtomMapping());
    }

    /**
     * Test of set method, of class Isomorphism.
     *
     * @throws CDKException
     * @throws IOException
     */
    @Test
    public void testSet_String_String() throws CDKException, IOException {
        ////System.out.println("6");
        String molfile = "C1CCC2CCCCC2C1";//decalin
        String queryfile = "C1CCC2CCCCC2C1";//decalin
        IAtomContainer query;
        IAtomContainer target;
        QueryAtomContainer query1 = null;
        QueryAtomContainer query2 = null;

        query = new SmilesParser(DefaultChemObjectBuilder.getInstance()).parseSmiles(molfile);
        target = new SmilesParser(DefaultChemObjectBuilder.getInstance()).parseSmiles(queryfile);

        Isomorphism smsd1 = new Isomorphism(query, target, Algorithm.DEFAULT, true, false, false);
        smsd1.setChemFilters(true, true, true);
        double score = 1.0;
        Assert.assertEquals(score, smsd1.getTanimotoSimilarity(), 0.0001);
    }

    /**
     * Test of set method, of class Isomorphism.
     *
     * @throws CDKException
     */
    @Test
    public void testSet_AtomContainer_AtomContainer() throws CDKException {
        ////System.out.println("7");
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());

        IAtomContainer target = sp.parseSmiles("C\\C=C/Nc1cccc(c1)N(O)\\C=C\\C\\C=C\\C=C/C");
        IAtomContainer query = sp.parseSmiles("Nc1ccccc1");
        Isomorphism smsd1 = new Isomorphism(query, target, Algorithm.DEFAULT, true, false, false);
        smsd1.setChemFilters(true, true, true);
        Assert.assertNotNull(smsd1.getFirstAtomMapping());
    }

    /**
     * Test of getAllAtomMapping method, of class Isomorphism.
     *
     * @throws CDKException
     */
    @Test
    public void testGetAllAtomMapping() throws CDKException {
        ////System.out.println("8");
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer target = sp.parseSmiles("C\\C=C/Nc1cccc(c1)N(O)\\C=C\\C\\C=C\\C=C/C");
        IAtomContainer query = sp.parseSmiles("Nc1ccccc1");

        SMSDNormalizer.percieveAtomTypesAndConfigureAtoms(query);
        SMSDNormalizer.percieveAtomTypesAndConfigureAtoms(target);

        Isomorphism smsd1 = new Isomorphism(query, target, Algorithm.DEFAULT, false, false, false);
        smsd1.setChemFilters(false, true, true);
        Assert.assertNotNull(smsd1.getFirstAtomMapping());
        Assert.assertEquals(2, smsd1.getAllAtomMapping().size());
    }

    /**
     * Test of getAllAtomMapping method, of class Isomorphism.
     *
     * @throws CDKException
     */
    @Test
    public void testgetAllAtomMapping() throws CDKException {
        ////System.out.println("9");
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer target = sp.parseSmiles("C\\C=C/Nc1cccc(c1)N(O)\\C=C\\C\\C=C\\C=C/C");
        IAtomContainer query = sp.parseSmiles("Nc1ccccc1");

        Isomorphism smsd1 = new Isomorphism(query, target, Algorithm.DEFAULT, true, false, false);
        smsd1.setChemFilters(false, true, true);
        Assert.assertNotNull(smsd1.getFirstAtomMapping());
        Assert.assertEquals(2, smsd1.getAllAtomMapping().size());
    }

    /**
     * Test of getFirstAtomMapping method, of class Isomorphism.
     *
     * @throws CDKException
     */
    @Test
    public void testGetFirstAtomMapping() throws CDKException {
        ////System.out.println("10");
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer target = sp.parseSmiles("C\\C=C/Nc1cccc(c1)N(O)\\C=C\\C\\C=C\\C=C/C");
        IAtomContainer query = sp.parseSmiles("Nc1ccccc1");

        Isomorphism smsd1 = new Isomorphism(query, target, Algorithm.DEFAULT, true, false, false);
        smsd1.setChemFilters(true, true, true);
        Assert.assertNotNull(smsd1.getFirstAtomMapping());
        Assert.assertEquals(7, smsd1.getFirstAtomMapping().getCount());
    }

    /**
     * Test of getFirstAtomMapping method, of class Isomorphism.
     *
     * @throws CDKException
     */
    @Test
    public void testGetFirstMapping() throws CDKException {
        ////System.out.println("11");
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer target = sp.parseSmiles("C\\C=C/Nc1cccc(c1)N(O)\\C=C\\C\\C=C\\C=C/C");
        IAtomContainer query = sp.parseSmiles("Nc1ccccc1");

        Isomorphism smsd1 = new Isomorphism(query, target, Algorithm.DEFAULT, true, false, false);
        smsd1.setChemFilters(true, true, true);
        Assert.assertNotNull(smsd1.getFirstAtomMapping());
        Assert.assertEquals(7, smsd1.getFirstAtomMapping().getCount());
    }

    /**
     * Test of setChemFilters method, of class Isomorphism.
     *
     * @throws CDKException
     */
    @Test
    public void testSetChemFilters() throws CDKException {
        ////System.out.println("12");
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer target = sp.parseSmiles("C\\C=C/OCC=C");
        IAtomContainer query = sp.parseSmiles("CCCOCC(C)=C");

        Isomorphism smsd1 = new Isomorphism(query, target, Algorithm.DEFAULT, true, true, false);
        smsd1.setChemFilters(true, true, true);
        Assert.assertEquals(1, smsd1.getAllAtomMapping().size());
        Assert.assertEquals(5, smsd1.getAllAtomMapping().listIterator().next().getCount());
    }

    /**
     * Test of getFragmentSize method, of class Isomorphism.
     *
     * @throws CDKException
     */
    @Test
    public void testGetFragmentSize() throws CDKException {
        ////System.out.println("13");
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer target = sp.parseSmiles("C\\C=C/Nc1cccc(c1)N(O)\\C=C\\C\\C=C\\C=C/C");
        IAtomContainer query = sp.parseSmiles("Nc1ccccc1");

        Isomorphism smsd1 = new Isomorphism(query, target, Algorithm.DEFAULT, true, false, true);
        smsd1.setChemFilters(false, true, false);
        Integer score = 2;
        Assert.assertEquals(score, smsd1.getFragmentSize(0));
    }

    /**
     * Test of getStereoScore method, of class Isomorphism.
     *
     * @throws CDKException
     */
    @Test
    public void testGetStereoScore() throws CDKException {
        ////System.out.println("14");
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer target = sp.parseSmiles("C\\C=C/OCC=C");
        IAtomContainer query = sp.parseSmiles("CCCOCC(C)=C");

        Isomorphism smsd1 = new Isomorphism(query, target, Algorithm.DEFAULT, false, false, false);
        smsd1.setChemFilters(true, false, false);
        Integer score = 77; //1024 to 77 due to new scoring fuction
        Assert.assertEquals(score, smsd1.getStereoScore(0));
    }

    /**
     * Test of getEnergyScore method, of class Isomorphism.
     *
     * @throws CDKException
     */
    @Test
    public void testGetEnergyScore() throws CDKException {
        ////System.out.println("15");
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer target = sp.parseSmiles("C\\C=C/Nc1cccc(c1)N(O)\\C=C\\C\\C=C\\C=C/C");
        IAtomContainer query = sp.parseSmiles("Nc1ccccc1");

        Isomorphism smsd1 = new Isomorphism(query, target, Algorithm.DEFAULT, true, false, false);
        smsd1.setChemFilters(false, false, true);
        Double score = 610.0;
        Assert.assertEquals(score, smsd1.getEnergyScore(0));
    }

    /**
     * Test of getQuery method, of class Isomorphism.
     *
     * @throws CDKException
     */
    @Test
    public void testGetReactantMolecule() throws CDKException {
        ////System.out.println("16");
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer target = sp.parseSmiles("C\\C=C/Nc1cccc(c1)N(O)\\C=C\\C\\C=C\\C=C/C");
        IAtomContainer query = sp.parseSmiles("Nc1ccccc1");

        Isomorphism smsd1 = new Isomorphism(query, target, Algorithm.DEFAULT, true, false, true);
        smsd1.setChemFilters(true, true, true);
        Assert.assertEquals(7, smsd1.getQuery().getAtomCount());
    }

    /**
     * Test of getTarget method, of class Isomorphism.
     *
     * @throws CDKException
     */
    @Test
    public void testGetProductMolecule() throws CDKException {
        ////System.out.println("17");
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer target = sp.parseSmiles("C\\C=C/Nc1cccc(c1)N(O)\\C=C\\C\\C=C\\C=C/C");
        IAtomContainer query = sp.parseSmiles("Nc1ccccc1");

        Isomorphism smsd1 = new Isomorphism(query, target, Algorithm.DEFAULT, true, false, true);
        smsd1.setChemFilters(true, true, true);
        Assert.assertEquals(20, smsd1.getTarget().getAtomCount());
    }

    /**
     * Test of getTanimotoSimilarity method, of class Isomorphism.
     *
     * @throws Exception
     */
    @Test
    public void testGetTanimotoSimilarity() throws Exception {
        ////System.out.println("18");
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer target = sp.parseSmiles("C\\C=C/Nc1cccc(c1)N(O)\\C=C\\C\\C=C\\C=C/C");
        IAtomContainer query = sp.parseSmiles("Nc1ccccc1");

        Isomorphism smsd1 = new Isomorphism(query, target, Algorithm.DEFAULT, true, false, true);
        smsd1.setChemFilters(true, true, true);

        double score = 0.35;
        Assert.assertEquals(score, smsd1.getTanimotoSimilarity(), 0);
    }

    /**
     * Test of isStereoMisMatch method, of class Isomorphism.
     *
     * @throws CDKException
     */
    @Test
    public void testIsStereoMisMatch() throws CDKException {
        ////System.out.println("19");
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer target = sp.parseSmiles("C\\C=C/Nc1cccc(c1)N(O)\\C=C\\C\\C=C\\C=C/C");
        IAtomContainer query = sp.parseSmiles("Nc1ccccc1");

        Isomorphism smsd1 = new Isomorphism(query, target, Algorithm.DEFAULT, false, false, true);
        smsd1.setChemFilters(true, true, true);
        Assert.assertFalse(smsd1.isStereoMisMatch());
    }

    @Test
    public void testQueryAtomContainerDefaultVF2() throws CDKException {
        ////System.out.println("20");
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer query = sp.parseSmiles("CC");
        IAtomContainer target = sp.parseSmiles("C1CCC12CCCC2");

        QueryAtomContainer queryContainer = QueryAtomContainerCreator.createSymbolAndBondOrderQueryContainer(query);
        Isomorphism smsd = new Isomorphism(queryContainer, target, Algorithm.VFLibMCS);
        boolean foundMatches = smsd.isSubgraph();
        Assert.assertTrue(foundMatches);
    }

    @Test
    public void testQueryAtomContainerMCSPLUS() throws CDKException {
        ////System.out.println("21");
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer query = sp.parseSmiles("CC");
        IAtomContainer target = sp.parseSmiles("C1CCC12CCCC2");
        Isomorphism smsd = new Isomorphism(query, target, Algorithm.MCSPlus, true, false, true);

        Assert.assertEquals(18, smsd.getAllAtomMapping().size());
        Assert.assertTrue(smsd.isSubgraph());

        IQueryAtomContainer queryContainer = QueryAtomContainerCreator.createSymbolAndBondOrderQueryContainer(query);
        smsd = new Isomorphism(queryContainer, target, Algorithm.MCSPlus);
        Assert.assertTrue(smsd.isSubgraph());
    }

    @Test
    public void testMatchAtomCount() throws CDKException {
        ////System.out.println("22");
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer query = sp.parseSmiles("CC");
        IAtomContainer target = sp.parseSmiles("C1CCC12CCCC2");
        Isomorphism smsd = new Isomorphism(query, target, Algorithm.CDKMCS, true, false, true);

        Assert.assertEquals(18, smsd.getAllAtomMapping().size());
        Assert.assertTrue(smsd.isSubgraph());

        IQueryAtomContainer queryContainer = QueryAtomContainerCreator.createSymbolAndBondOrderQueryContainer(query);
        smsd = new Isomorphism(queryContainer, target, Algorithm.DEFAULT);
        Assert.assertTrue(smsd.isSubgraph());
    }

    @Test
    public void testMatchCountCDKMCS() throws CDKException {
        ////System.out.println("23");
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer query = sp.parseSmiles("CC");
        IAtomContainer target = sp.parseSmiles("C1CCC12CCCC2");
        Isomorphism smsd = new Isomorphism(query, target, Algorithm.CDKMCS, true, false, true);
        Assert.assertEquals(18, smsd.getAllAtomMapping().size());
        Assert.assertTrue(smsd.isSubgraph());
    }

    /**
     * Test Tanimoto NAD+ & NADH for Bond Sensitive.
     *
     * @throws Exception
     */
    @Test
    public void testNADPlusNADHBondSensitive() throws Exception {
        ////System.out.println("24");
        SmilesParser smilesParser = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer molecule1 = smilesParser.parseSmiles("NC(=O)c1ccc[n+](c1)[C@@H]1O[C@H](COP(O)(=O)OP(O)(=O)OC[C@H]2O[C@H]([C@H](O)[C@@H]2O)n2cnc3c(N)ncnc23)[C@@H](O)[C@H]1O");
        IAtomContainer molecule2 = smilesParser.parseSmiles("NC(=O)C1=CN(C=CC1)[C@@H]1O[C@H](COP(O)(=O)OP(O)(=O)OC[C@H]2O[C@H]([C@H](O)[C@@H]2O)n2cnc3c(N)ncnc23)[C@@H](O)[C@H]1O");

        double score = 0.6923;
        Isomorphism comparison = new Isomorphism(molecule1, molecule2, Algorithm.DEFAULT, true, false, true);
        comparison.setChemFilters(true, true, true);
//        ////System.out.println("SMILES :" + comparison.getFirstAtomMapping().getCommonFragmentAsSMILES());
        Assert.assertEquals(score, comparison.getTanimotoSimilarity(), 0.001);
    }

    /**
     * Test ring match using MCS VF2Plus
     *
     * @throws Exception
     */
    @Test
    public void testCDKMCS() throws Exception {
        ////System.out.println("25");
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        // Benzene
        IAtomContainer query = sp.parseSmiles("C1=CC=CC=C1");
        // Napthalene
        IAtomContainer target = sp.parseSmiles("C1=CC2=C(C=C1)C=CC=C2");
        Isomorphism comparison = new Isomorphism(query, target, Algorithm.CDKMCS, true, true, true);
        // set chemical filter true
        comparison.setChemFilters(true, true, true);

        Assert.assertEquals(0.6, comparison.getTanimotoSimilarity());
        Assert.assertEquals(12, comparison.getAllAtomMapping().size());
    }

    /**
     * Test ring match using MCS VF2Plus
     *
     * @throws Exception
     */
    @Test
    public void testVF2MCS() throws Exception {
        ////System.out.println("26");
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        // Benzene
        IAtomContainer query = sp.parseSmiles("C1=CC=CC=C1");
        // Napthalene
        IAtomContainer target = sp.parseSmiles("C1=CC2=C(C=C1)C=CC=C2");
        //{ 0: Default Isomorphism Algorithm, 1: MCSPlus Algorithm, 2: VFLibMCS Algorithm, 3: CDKMCS Algorithm}
        //Algorithm is VF2MCS
        //Bond Sensitive is set True
        //Ring Match is set True

        Isomorphism comparison = new Isomorphism(query, target, Algorithm.VFLibMCS, true, true, true);
        // set chemical filter true
        comparison.setChemFilters(true, true, true);

        //Get similarity score
//        //////System.out.println("Tanimoto coefficient:  " + comparison.getTanimotoSimilarity());
        Assert.assertEquals(0.6, comparison.getTanimotoSimilarity());
        Assert.assertEquals(12, comparison.getAllAtomMapping().size());
        // Print the mapping between molecules
//        //////System.out.println(" Mappings: ");
//        for (AtomAtomMapping atomatomMapping : comparison.getAllAtomMapping()) {
//            for (Map.Entry<IAtom, IAtom> mapping : atomatomMapping.getMappingsByAtoms().entrySet()) {
//                IAtom sourceAtom = mapping.getKey();
//                IAtom targetAtom = mapping.getValue();
//                //////System.out.println(sourceAtom.getSymbol() + " " + targetAtom.getSymbol());
//                //////System.out.println(atomatomMapping.getQueryIndex(sourceAtom) + " " + atomatomMapping.getTargetIndex(targetAtom));
//            }
//            //////System.out.println("");
//        }
    }

    /**
     * Test ring match using MCS MCSPlus
     *
     * @throws Exception
     */
    @Test
    public void testMCSPlus() throws Exception {
        ////System.out.println("27");
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        // Benzene
        IAtomContainer query = sp.parseSmiles("C1=CC=CC=C1");
        // Napthalene
        IAtomContainer target = sp.parseSmiles("C1=CC2=C(C=C1)C=CC=C2");
        //{ 0: Default Isomorphism Algorithm, 1: MCSPlus Algorithm, 2: VFLibMCS Algorithm, 3: CDKMCS Algorithm}
        //Algorithm is MCSPlus
        //Bond Sensitive is set True
        //Ring Match is set True

        Isomorphism comparison = new Isomorphism(query, target, Algorithm.MCSPlus, false, false, true);
        // set chemical filter true
        comparison.setChemFilters(true, true, true);
        Assert.assertEquals(0.6, comparison.getTanimotoSimilarity());
        Assert.assertEquals(12, comparison.getAllAtomMapping().size());
    }

    /**
     * Test Tanimoto NAD+ & NADH for Bond InSensitive. Time taking cases
     *
     * @throws Exception
     */
    @Test
    public void testNADPlusNADHBondInSensitive() throws Exception {
        ////System.out.println("28");
        SmilesParser smilesParser = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer query = smilesParser.parseSmiles("NC(=O)c1ccc[n+](c1)[C@@H]1O[C@H](COP(O)(=O)OP(O)(=O)OC[C@H]2O[C@H]([C@H](O)[C@@H]2O)n2cnc3c(N)ncnc23)[C@@H](O)[C@H]1O");
        IAtomContainer target = smilesParser.parseSmiles("NC(=O)C1=CN(C=CC1)[C@@H]1O[C@H](COP(O)(=O)OP(O)(=O)OC[C@H]2O[C@H]([C@H](O)[C@@H]2O)n2cnc3c(N)ncnc23)[C@@H](O)[C@H]1O");

        IAtomContainer ac1 = ExtAtomContainerManipulator.removeHydrogens(query);
        IAtomContainer ac2 = ExtAtomContainerManipulator.removeHydrogens(target);
        double score = 1.0;
        Isomorphism smsd1 = new Isomorphism(ac1, ac2, Algorithm.DEFAULT, false, false, false);
        smsd1.setChemFilters(true, true, true);
        Assert.assertEquals(score, smsd1.getTanimotoSimilarity(), 0.001);
    }

//    @Test
//    public void testComplex1() throws Exception {
//        ////System.out.println("29");
//        SmilesParser smilesParser = new SmilesParser(DefaultChemObjectBuilder.getInstance());
//        IAtomContainer query = smilesParser.parseSmiles("C[C@@H](NC(=O)[C@@H](C)NC(=O)[C@H](CCCCN)NC(=O)CC[C@@H](NC(=O)[C@H](C)NC(=O)[C@@H](C)O[C@@H]1[C@@H](NC(C)=O)[C@@H](OP(O)(=O)OP(O)(=O)OC\\C=C(\\C)CC\\C=C(\\C)CC\\C=C(\\C)CC\\C=C(\\C)CC\\C=C(\\C)CC\\C=C(\\C)CC\\C=C(\\C)CC\\C=C(\\C)CC\\C=C(/C)CC\\C=C(/C)CCC=C(C)C)O[C@H](CO)[C@H]1O[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1NC(C)=O)C(O)=O)C(O)=O");
//        IAtomContainer target = smilesParser.parseSmiles("C[C@@H](NC(=O)[C@@H](C)NC(=O)[C@H](CCCCN)NC(=O)CC[C@@H](NC(=O)[C@H](C)NC(=O)[C@@H](C)O[C@@H]1[C@@H](NC(C)=O)[C@@H](OP(O)(=O)OP(O)(=O)OC\\C=C(\\C)CC\\C=C(\\C)CC\\C=C(\\C)CC\\C=C(\\C)CC\\C=C(\\C)CC\\C=C(\\C)CC\\C=C(\\C)CC\\C=C(\\C)CC\\C=C(/C)CC\\C=C(/C)CCC=C(C)C)O[C@H](CO)[C@H]1O[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1NC(C)=O)C(N)=O)C(O)=O");
//        IAtomContainer ac1 = ExtAtomContainerManipulator.removeHydrogens(query);
//        IAtomContainer ac2 = ExtAtomContainerManipulator.removeHydrogens(target);
//        double score = 0.9847;
//        Isomorphism smsd1 = new Isomorphism(ac1, ac2, Algorithm.VFLibMCS, false, false, false);
//        smsd1.setChemFilters(true, true, true);
////        ////System.out.println(" overlap " + smsd1.getFirstAtomMapping().getCommonFragmentAsSMILES());
//        Assert.assertEquals(score, smsd1.getTanimotoSimilarity(), 0.001);
//    }
//
//    @Test
//    public void testComplex2() throws Exception {
//        ////System.out.println("30");
//        SmilesParser smilesParser = new SmilesParser(DefaultChemObjectBuilder.getInstance());
//        IAtomContainer query = smilesParser.parseSmiles("NC1=NC=NC2=C1N=CN2[C@@H]1O[C@H](COP(O)(=O)OP(O)(O)=O)[C@@H](O)[C@H]1O");
//        IAtomContainer target = smilesParser.parseSmiles("C[C@@H](NC(=O)[C@@H](C)NC(=O)[C@H](CCCCN)NC(=O)CC[C@@H](NC(=O)[C@H](C)NC(=O)[C@@H](C)O[C@@H]1[C@@H](NC(C)=O)[C@@H](OP(O)(=O)OP(O)(=O)OC\\C=C(\\C)CC\\C=C(\\C)CC\\C=C(\\C)CC\\C=C(\\C)CC\\C=C(\\C)CC\\C=C(\\C)CC\\C=C(\\C)CC\\C=C(\\C)CC\\C=C(/C)CC\\C=C(/C)CCC=C(C)C)O[C@H](CO)[C@H]1O[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1NC(C)=O)C(N)=O)C(O)=O");
//        IAtomContainer ac1 = ExtAtomContainerManipulator.removeHydrogens(query);
//        IAtomContainer ac2 = ExtAtomContainerManipulator.removeHydrogens(target);
//        double score = 0.1214;
//        Isomorphism smsd1 = new Isomorphism(ac1, ac2, Algorithm.VFLibMCS, false, false, false);
//        smsd1.setChemFilters(true, true, true);
////        ////System.out.println(" overlap " + smsd1.getFirstAtomMapping().getCommonFragmentAsSMILES());
//        Assert.assertEquals(score, smsd1.getTanimotoSimilarity(), 0.001);
//    }
//
//    @Test
//    public void testComplex3() throws Exception {
//        ////System.out.println("31");
//        SmilesParser smilesParser = new SmilesParser(DefaultChemObjectBuilder.getInstance());
//        IAtomContainer query = smilesParser.parseSmiles("COc1cccc(C\\C=C(/C)CC\\C=C(/C)CC\\C=C(/C)CC\\C=C(/C)CC\\C=C(/C)CC\\C=C(/C)CC\\C=C(/C)CCC=C(C)C)c1O");
//        IAtomContainer target = smilesParser.parseSmiles("COC1=CC(=O)C=C(C\\C=C(/C)CC\\C=C(/C)CC\\C=C(/C)CC\\C=C(/C)CC\\C=C(/C)CC\\C=C(/C)CC\\C=C(/C)CCC=C(C)C)C1=O");
//        IAtomContainer ac1 = ExtAtomContainerManipulator.removeHydrogens(query);
//        IAtomContainer ac2 = ExtAtomContainerManipulator.removeHydrogens(target);
//        double score = 0.98;
//        Isomorphism smsd1 = new Isomorphism(ac1, ac2, Algorithm.VFLibMCS, false, false, false);
//        smsd1.setChemFilters(true, true, true);
////        ////System.out.println(" overlap " + smsd1.getFirstAtomMapping().getCommonFragmentAsSMILES());
//        Assert.assertEquals(score, smsd1.getTanimotoSimilarity(), 0.001);
//        Assert.assertEquals(2, smsd1.getAllAtomMapping().size());
//    }
    /**
     * Time taking cases Test ring match using MCS VF2Plus
     *
     * @throws Exception
     */
    @Test
    public void testOpenCloseRing() throws Exception {
    ////System.out.println("31");
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        //(S)-2,3-epoxysqualene (CHEBI:15441)
        IAtomContainer query = sp.parseSmiles("CC(C)=CCC\\C(C)=C\\CC\\C(C)=C\\CC\\C=C(/C)CC\\C=C(/C)CC[C@@H]1OC1(C)C");
        // lupeol (CHEBI:6570)
        IAtomContainer target = sp.parseSmiles("[H][C@]12[C@@H](CC[C@]1(C)CC[C@]1(C)[C@]2([H])CC[C@]2([H])[C@@]3(C)CC[C@H](O)C(C)(C)[C@]3([H])CC[C@@]12C)C(C)=C");

        IAtomContainer ac1 = ExtAtomContainerManipulator.removeHydrogens(query);
        IAtomContainer ac2 = ExtAtomContainerManipulator.removeHydrogens(target);

        Isomorphism comparison = new Isomorphism(ac1, ac2, Algorithm.VFLibMCS, false, false, false);
        // set chemical filter true
        comparison.setChemFilters(true, true, true);
        SmilesGenerator aromatic = SmilesGenerator.unique().aromatic();
        ////System.out.println("SMILES :" + aromatic.create(comparison.getFirstAtomMapping().getMapCommonFragmentOnQuery()));
        ////System.out.println("SMILES :" + aromatic.create(comparison.getFirstAtomMapping().getMapCommonFragmentOnTarget()));
        ////System.out.println("SMILES :" + comparison.getFirstAtomMapping().getCommonFragmentAsSMILES());

        Assert.assertEquals(0.8235, comparison.getTanimotoSimilarity(), .001);
        Assert.assertEquals(2, comparison.getAllAtomMapping().size());
    }
    /**
     * Test Common Fragment Expected SMLIES is Given the two molecules CCCNCC &
     * CNCCS the MCS returned is N([CH2])CC which is correct. But it could also
     * have been [CH2]CN[CH2]
     *
     * @throws Exception
     */
    @Test
    public void testCommonFragment() throws Exception {
        ////System.out.println("32");
        SmilesParser smilesParser = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer query = smilesParser.parseSmiles("CCCNCC");//("CNCCS");
        IAtomContainer target = smilesParser.parseSmiles("CNCCS");//("CCCNCC");

        IAtomContainer ac1 = ExtAtomContainerManipulator.removeHydrogens(query);
        IAtomContainer ac2 = ExtAtomContainerManipulator.removeHydrogens(target);
        double score = 0.5714;
        Isomorphism overlap = new Isomorphism(ac1, ac2, Algorithm.VFLibMCS, true, true, true);
        overlap.setChemFilters(true, true, true);
        Assert.assertEquals(score, overlap.getTanimotoSimilarity(), 0.001);
        SmilesGenerator aromatic = SmilesGenerator.unique().aromatic();
//        ////System.out.println("SMILES Q :" + aromatic.create(overlap.getFirstAtomMapping().getMapCommonFragmentOnQuery()));
//        ////System.out.println("SMILES T :" + aromatic.create(overlap.getFirstAtomMapping().getMapCommonFragmentOnTarget()));
//        ////System.out.println("SMILES Common:" + overlap.getFirstAtomMapping().getCommonFragmentAsSMILES());
//
    }

    /**
     * Test Common Fragment Expected SMLIES is Given the two molecules CCCNCC &
     * CNCCS the MCS returned is N([CH2])CC which is correct. But it could also
     * have been [CH2]CN[CH2]
     *
     * @throws Exception
     */
    @Test
    public void testAliphaticToRing() throws Exception {
        ////System.out.println("32");
        SmilesParser smilesParser = new SmilesParser(DefaultChemObjectBuilder.getInstance());
//        IAtomContainer query = smilesParser.parseSmiles("CC(C)=CCC\\C(C)=C\\CC\\C(C)=C\\CC\\C=C(/C)CC\\C=C(/C)CC[C@@H]1OC1(C)C");//("CNCCS");
//        IAtomContainer target = smilesParser.parseSmiles("[H][C@]12[C@@H](CC[C@]1(C)CC[C@]1(C)[C@]2([H])CC[C@]2([H])[C@@]3(C)CC[C@H](O)C(C)(C)[C@]3([H])CC[C@@]12C)C(C)=C");//("CCCNCC");

        IAtomContainer query = smilesParser.parseSmiles("O1C(CCC(=CCCC(=CCCC=C(C)CCC=C(C)CCC=C(C)C)C)C)C1(C)C");//("CNCCS");
        IAtomContainer target = smilesParser.parseSmiles("OC1CCC2(C)C(CCC3(C)C2CCC4C5C(C(=C)C)CCC5(C)CCC43C)C1(C)C");

        IAtomContainer ac1 = ExtAtomContainerManipulator.removeHydrogens(query);
        IAtomContainer ac2 = ExtAtomContainerManipulator.removeHydrogens(target);
        double score = 0.8235;
        Isomorphism overlap = new Isomorphism(ac1, ac2, Algorithm.VFLibMCS, false, false, false);
        overlap.setChemFilters(true, true, true);
        SmilesGenerator aromatic = SmilesGenerator.unique().aromatic();
        System.out.println("SMILES Q :" + aromatic.create(overlap.getFirstAtomMapping().getMapCommonFragmentOnQuery()));
        System.out.println("SMILES T :" + aromatic.create(overlap.getFirstAtomMapping().getMapCommonFragmentOnTarget()));
        System.out.println(" SMILES Common: " + overlap.getFirstAtomMapping().getCount() + ", " + overlap.getFirstAtomMapping().getCommonFragmentAsSMILES());
        Assert.assertEquals(score, overlap.getTanimotoSimilarity(), 0.001);

    }

    @Test
    public void testjg1() throws Exception {
        ////System.out.println("32");
        SmilesParser smilesParser = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer query = smilesParser.parseSmiles("COc1ccccc1C=CC=c2c(=O)n3c(s2)nc(n3)c4cccs4");//("CNCCS");
        IAtomContainer target = smilesParser.parseSmiles("COc1ccccc1C=CC=c2c(=O)n3c(s2)nc(n3)c4cccs4");

        IAtomContainer ac1 = ExtAtomContainerManipulator.removeHydrogens(query);
        IAtomContainer ac2 = ExtAtomContainerManipulator.removeHydrogens(target);
        double score = 1.0;
        Isomorphism overlap = new Isomorphism(ac1, ac2, Algorithm.VFLibMCS, false, false, false);
        overlap.setChemFilters(true, true, true);
        SmilesGenerator aromatic = SmilesGenerator.unique().aromatic();
        System.out.println("SMILES Q :" + aromatic.create(overlap.getFirstAtomMapping().getMapCommonFragmentOnQuery()));
        System.out.println("SMILES T :" + aromatic.create(overlap.getFirstAtomMapping().getMapCommonFragmentOnTarget()));
        System.out.println(" SMILES Common: " + overlap.getFirstAtomMapping().getCount() + ", " + overlap.getFirstAtomMapping().getCommonFragmentAsSMILES());
        Assert.assertEquals(score, overlap.getTanimotoSimilarity(), 0.001);

    }

    @Test
    public void testjg2() throws Exception {
        ////System.out.println("32");
        /*
         CDK RMAP bug, equals is || not &&
         */
        SmilesParser smilesParser = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer query = smilesParser.parseSmiles("O=C(c1cn(c2c1cccc2)CCN3CCOCC3)C4C(C4(C)C)(C)C");//("CNCCS");
        setID(query);
        IAtomContainer target = smilesParser.parseSmiles("C(ON=O)C(C)C");
        setID(target);

        IAtomContainer ac1 = ExtAtomContainerManipulator.removeHydrogens(query);
        IAtomContainer ac2 = ExtAtomContainerManipulator.removeHydrogens(target);
        double score = 0.1379;
        Isomorphism overlap = new Isomorphism(ac1, ac2, Algorithm.CDKMCS, true, false, false);
        overlap.setChemFilters(true, true, true);
        SmilesGenerator aromatic = SmilesGenerator.unique().aromatic();
        System.out.println("SMILES Q :" + aromatic.create(overlap.getFirstAtomMapping().getMapCommonFragmentOnQuery()));
        System.out.println("SMILES T :" + aromatic.create(overlap.getFirstAtomMapping().getMapCommonFragmentOnTarget()));
        System.out.println(" SMILES Common: " + overlap.getFirstAtomMapping().getCount() + ", " + overlap.getFirstAtomMapping().getCommonFragmentAsSMILES());
        Assert.assertEquals(score, overlap.getTanimotoSimilarity(), 0.001);

    }

    @Test
    public void testjg3() throws Exception {
        ////System.out.println("32");
        /*
         Now check ring comparision, CDK RMAP bug, equals is || not &&
         */
        SmilesParser smilesParser = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer query = smilesParser.parseSmiles("C(NC)(C)CC1=CCC=CC1");//("CNCCS");
        setID(query);
        IAtomContainer target = smilesParser.parseSmiles("c12OCOc1cc(cc2)C(=O)C(N(C)C)CC");
        setID(target);

        IAtomContainer ac1 = ExtAtomContainerManipulator.removeHydrogens(query);
        IAtomContainer ac2 = ExtAtomContainerManipulator.removeHydrogens(target);
        double score = 0.6471;
        Isomorphism overlap = new Isomorphism(ac1, ac2, Algorithm.DEFAULT, false, false, false);
        overlap.setChemFilters(true, true, true);
        SmilesGenerator aromatic = SmilesGenerator.unique().aromatic();
        System.out.println("SMILES Q :" + aromatic.create(overlap.getFirstAtomMapping().getMapCommonFragmentOnQuery()));
        System.out.println("SMILES T :" + aromatic.create(overlap.getFirstAtomMapping().getMapCommonFragmentOnTarget()));
        System.out.println(" SMILES Common: " + overlap.getFirstAtomMapping().getCount() + ", " + overlap.getFirstAtomMapping().getCommonFragmentAsSMILES());
        Assert.assertEquals(score, overlap.getTanimotoSimilarity(), 0.001);

    }

    /*
     g1: NC1CCCCC1
     g2: NC1CCC(N)CC1
     g3: CNC1CCC(N)CC1
     g4: CNC1CCC(CC1)NC

     The correct number of all MCS mappings should be

     |MCS(g1, g1)| = 2
     |MCS(g1, g2)| = 4
     |MCS(g1, g3)| = 4
     |MCS(g1, g4)| = 4
     */
    @Test
    public void testSolutionCount() throws Exception {
        ////System.out.println("33");
        String g1 = "NC1CCCCC1";
        String g2 = "NC1CCC(N)CC1";
        String g3 = "CNC1CCC(N)CC1";
        String g4 = "CNC1CCC(CC1)NC";

//        String g1 = "OC1CCCCC1";
//        String g2 = "OC1CCC(O)CC1";
//        String g3 = "COC1CCC(O)CC1";
//        String g4 = "COC1CCC(CC1)OC";
        SmilesParser smilesParser = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer ac1 = smilesParser.parseSmiles(g1);
        setID(ac1);
        IAtomContainer ac2 = smilesParser.parseSmiles(g2);
        setID(ac2);
        IAtomContainer ac3 = smilesParser.parseSmiles(g3);
        setID(ac3);
        IAtomContainer ac4 = smilesParser.parseSmiles(g4);
        setID(ac4);

        Isomorphism overlap = new Isomorphism(ac1, ac2, Algorithm.DEFAULT, true, false, false);
        overlap.setChemFilters(false, false, true);
        Assert.assertEquals(4, overlap.getAllAtomMapping().size());
        overlap = new Isomorphism(ac1, ac3, Algorithm.DEFAULT, true, false, false);
        overlap.setChemFilters(false, false, true);
        Assert.assertEquals(4, overlap.getAllAtomMapping().size());
        overlap = new Isomorphism(ac1, ac4, Algorithm.DEFAULT, true, false, false);
        overlap.setChemFilters(false, false, true);
        Assert.assertEquals(4, overlap.getAllAtomMapping().size());
        overlap = new Isomorphism(ac1, ac1, Algorithm.DEFAULT, true, false, false);
        overlap.setChemFilters(false, false, true);
        Assert.assertEquals(2, overlap.getAllAtomMapping().size());

//        SmilesGenerator aromatic = SmilesGenerator.unique().aromatic();
//        ////System.out.println("SMILES Q :" + aromatic.create(overlap.getFirstAtomMapping().getMapCommonFragmentOnQuery()));
//        ////System.out.println("SMILES T :" + aromatic.create(overlap.getFirstAtomMapping().getMapCommonFragmentOnTarget()));
//        ////System.out.println("SMILES Common:" + overlap.getFirstAtomMapping().getCommonFragmentAsSMILES());
//
    }

    private void setID(IAtomContainer ac) {
        for (IAtom a : ac.atoms()) {
            a.setID(ac.getAtomNumber(a) + "");
        }
    }
}
