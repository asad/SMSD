
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
import java.io.InputStream;
import java.util.logging.Level;
import java.util.logging.Logger;
import junit.framework.Assert;
import org.junit.*;
import org.openscience.cdk.AtomContainer;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.smsd.tools.Utility;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.io.IChemObjectReader.Mode;
import org.openscience.cdk.io.MDLV2000Reader;
import org.openscience.cdk.isomorphism.matchers.IQueryAtomContainer;
import org.openscience.cdk.isomorphism.matchers.QueryAtomContainerCreator;
import org.openscience.cdk.normalize.SMSDNormalizer;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import org.openscience.smsd.algorithm.mcsplus.MCSPlusHandlerTest;
import org.openscience.smsd.interfaces.Algorithm;

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
//        //System.out.println("init");
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
//        //System.out.println("init");
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
//            //System.out.println("searchMCS");
            SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
            IAtomContainer target;
            target = sp.parseSmiles("C\\C=C/Nc1cccc(c1)N(O)\\C=C\\C\\C=C\\C=C/C");
            IAtomContainer query = sp.parseSmiles("Nc1ccccc1");
            Isomorphism smsd1 = new Isomorphism(query, target, Algorithm.DEFAULT, true, false, false);
            smsd1.setChemFilters(true, true, true);
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
//        //System.out.println("set");
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
//        //System.out.println("set");
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
//        //System.out.println("set");
        String molfile = "data/mdl/decalin.mol";
        String queryfile = "data/mdl/decalin.mol";
        AtomContainer query = new AtomContainer();
        AtomContainer target = new AtomContainer();

        InputStream ins = this.getClass().getClassLoader().getResourceAsStream(molfile);
        MDLV2000Reader reader = new MDLV2000Reader(ins, Mode.STRICT);
        reader.read(query);
        ins = this.getClass().getClassLoader().getResourceAsStream(queryfile);
        reader = new MDLV2000Reader(ins, Mode.STRICT);
        reader.read(target);

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
//        //System.out.println("set");
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
//        //System.out.println("getAllAtomMapping");
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer target = sp.parseSmiles("C\\C=C/Nc1cccc(c1)N(O)\\C=C\\C\\C=C\\C=C/C");
        IAtomContainer query = sp.parseSmiles("Nc1ccccc1");

        SMSDNormalizer.percieveAtomTypesAndConfigureAtoms(query);
        SMSDNormalizer.percieveAtomTypesAndConfigureAtoms(target);

//	Calling the main algorithm to perform MCS cearch
        Utility.aromatizeDayLight(query);
        Utility.aromatizeDayLight(target);

        Isomorphism smsd1 = new Isomorphism(query, target, Algorithm.DEFAULT, true, false, false);
        smsd1.setChemFilters(true, true, true);
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
//        //System.out.println("getAllAtomMapping");
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer target = sp.parseSmiles("C\\C=C/Nc1cccc(c1)N(O)\\C=C\\C\\C=C\\C=C/C");
        IAtomContainer query = sp.parseSmiles("Nc1ccccc1");

        Isomorphism smsd1 = new Isomorphism(query, target, Algorithm.DEFAULT, true, false, false);
        smsd1.setChemFilters(true, true, true);
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
//        //System.out.println("getFirstAtomMapping");
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
//        //System.out.println("getFirstAtomMapping");
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
//        //System.out.println("setChemFilters");
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
//        //System.out.println("getFragmentSize");
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
//        //System.out.println("getStereoScore");
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
//        //System.out.println("getEnergyScore");
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
//        //System.out.println("getQuery");
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
//        //System.out.println("getTarget");
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
//        //System.out.println("getTanimotoSimilarity");
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
//        //System.out.println("isStereoMisMatch");
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer target = sp.parseSmiles("C\\C=C/Nc1cccc(c1)N(O)\\C=C\\C\\C=C\\C=C/C");
        IAtomContainer query = sp.parseSmiles("Nc1ccccc1");

        Isomorphism smsd1 = new Isomorphism(query, target, Algorithm.DEFAULT, false, false, true);
        smsd1.setChemFilters(true, true, true);
        Assert.assertFalse(smsd1.isStereoMisMatch());
    }

    @Test
    public void testQueryAtomContainerDefaultVF2() throws CDKException {
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer query = sp.parseSmiles("CC");
        IAtomContainer target = sp.parseSmiles("C1CCC12CCCC2");
        Isomorphism smsd = new Isomorphism(query, target, Algorithm.VFLibMCS, true, false, true);

        boolean foundMatches = smsd.isSubgraph();
        Assert.assertEquals(18, smsd.getAllAtomMapping().size());
        Assert.assertTrue(foundMatches);

        IQueryAtomContainer queryContainer = QueryAtomContainerCreator.createSymbolAndBondOrderQueryContainer(query);
        smsd = new Isomorphism(queryContainer, target, Algorithm.VFLibMCS);
        foundMatches = smsd.isSubgraph();
        Assert.assertTrue(foundMatches);
    }

    @Test
    public void testQueryAtomContainerMCSPLUS() throws CDKException {
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
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer query = sp.parseSmiles("CC");
        IAtomContainer target = sp.parseSmiles("C1CCC12CCCC2");
        Isomorphism smsd = new Isomorphism(query, target, Algorithm.DEFAULT, true, false, true);

        Assert.assertEquals(18, smsd.getAllAtomMapping().size());
        Assert.assertTrue(smsd.isSubgraph());

        IQueryAtomContainer queryContainer = QueryAtomContainerCreator.createSymbolAndBondOrderQueryContainer(query);
        smsd = new Isomorphism(queryContainer, target, Algorithm.DEFAULT);
        Assert.assertTrue(smsd.isSubgraph());
    }

    @Test
    public void testMatchCountCDKMCS() throws CDKException {
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
//        //System.out.println("getTanimoto for NAD+ and NADH");
        SmilesParser smilesParser = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer molecule1 = smilesParser.parseSmiles("NC(=O)c1ccc[n+](c1)[C@@H]1O[C@H](COP(O)(=O)OP(O)(=O)OC[C@H]2O[C@H]([C@H](O)[C@@H]2O)n2cnc3c(N)ncnc23)[C@@H](O)[C@H]1O");
        IAtomContainer molecule2 = smilesParser.parseSmiles("NC(=O)C1=CN(C=CC1)[C@@H]1O[C@H](COP(O)(=O)OP(O)(=O)OC[C@H]2O[C@H]([C@H](O)[C@@H]2O)n2cnc3c(N)ncnc23)[C@@H](O)[C@H]1O");

        double score = 0.7;
        Isomorphism smsd1 = new Isomorphism(molecule1, molecule2, Algorithm.MCSPlus, true, false, true);
        //smsd1.setChemFilters(true, true, true);
        Assert.assertEquals(score, smsd1.getTanimotoSimilarity(), 0.0999);
    }

    /**
     * Test ring match using MCS VF2Plus
     *
     * @throws Exception
     */
    @Test
    public void testCDKMCS() throws Exception {
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        // Benzene
        IAtomContainer query = sp.parseSmiles("C1=CC=CC=C1");
        // Napthalene
        IAtomContainer target = sp.parseSmiles("C1=CC2=C(C=C1)C=CC=C2");
        //{ 0: Default Isomorphism Algorithm, 1: MCSPlus Algorithm, 2: VFLibMCS Algorithm, 3: CDKMCS Algorithm}
        //Algorithm is VF2MCS
        //Bond Sensitive is set True
        //Ring Match is set True

        Isomorphism comparison = new Isomorphism(query, target, Algorithm.CDKMCS, true, true, true);
        // set chemical filter true
        comparison.setChemFilters(true, true, true);

        //Get similarity score
//        //System.out.println("Tanimoto coefficient:  " + comparison.getTanimotoSimilarity());
        Assert.assertEquals(0.6, comparison.getTanimotoSimilarity());
        Assert.assertEquals(12, comparison.getAllAtomMapping().size());
        // Print the mapping between molecules
//        //System.out.println(" Mappings: ");
//        for (AtomAtomMapping atomatomMapping : comparison.getAllAtomMapping()) {
//            for (Map.Entry<IAtom, IAtom> mapping : atomatomMapping.getMappingsByAtoms().entrySet()) {
//                IAtom sourceAtom = mapping.getKey();
//                IAtom targetAtom = mapping.getValue();
//                //System.out.println(sourceAtom.getSymbol() + " " + targetAtom.getSymbol());
//                //System.out.println(atomatomMapping.getQueryIndex(sourceAtom) + " " + atomatomMapping.getTargetIndex(targetAtom));
//            }
//            //System.out.println("");
//        }
    }

    /**
     * Test ring match using MCS VF2Plus
     *
     * @throws Exception
     */
    @Test
    public void testVF2MCS() throws Exception {
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
//        //System.out.println("Tanimoto coefficient:  " + comparison.getTanimotoSimilarity());
        Assert.assertEquals(0.6, comparison.getTanimotoSimilarity());
        Assert.assertEquals(12, comparison.getAllAtomMapping().size());
        // Print the mapping between molecules
//        //System.out.println(" Mappings: ");
//        for (AtomAtomMapping atomatomMapping : comparison.getAllAtomMapping()) {
//            for (Map.Entry<IAtom, IAtom> mapping : atomatomMapping.getMappingsByAtoms().entrySet()) {
//                IAtom sourceAtom = mapping.getKey();
//                IAtom targetAtom = mapping.getValue();
//                //System.out.println(sourceAtom.getSymbol() + " " + targetAtom.getSymbol());
//                //System.out.println(atomatomMapping.getQueryIndex(sourceAtom) + " " + atomatomMapping.getTargetIndex(targetAtom));
//            }
//            //System.out.println("");
//        }
    }

    /**
     * Test ring match using MCS MCSPlus
     *
     * @throws Exception
     */
    @Test
    public void testMCSPlus() throws Exception {
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

//        SmilesGenerator sg = new SmilesGenerator();
//        String createSMILES = sg.create(comparison.getFirstAtomMapping().getCommonFragmentInTarget());
//        System.out.println("createSMILES " + createSMILES);
        //Get similarity score
        //System.out.println("Tanimoto coefficient:  " + comparison.getTanimotoSimilarity());
        Assert.assertEquals(0.6, comparison.getTanimotoSimilarity());
        Assert.assertEquals(12, comparison.getAllAtomMapping().size());
    }

    /**
     * Time taking cases Test ring match using MCS VF2Plus
     *
     * @throws Exception
     */
    @Test
    public void testOpenCloseRing() throws Exception {
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        //(S)-2,3-epoxysqualene (CHEBI:15441)
        IAtomContainer query = sp.parseSmiles("CC(C)=CCC\\C(C)=C\\CC\\C(C)=C\\CC\\C=C(/C)CC\\C=C(/C)CC[C@@H]1OC1(C)C");
        // lupeol (CHEBI:6570)
        IAtomContainer target = sp.parseSmiles("[H][C@]12[C@@H](CC[C@]1(C)CC[C@]1(C)[C@]2([H])CC[C@]2([H])[C@@]3(C)CC[C@H](O)C(C)(C)[C@]3([H])CC[C@@]12C)C(C)=C");

        IAtomContainer ac1 = AtomContainerManipulator.removeHydrogens(query);
        IAtomContainer ac2 = AtomContainerManipulator.removeHydrogens(target);

        Isomorphism comparison = new Isomorphism(ac1, ac2, Algorithm.DEFAULT, false, false, false);
        // set chemical filter true
        comparison.setChemFilters(true, true, true);
//        SmilesGenerator aromatic = SmilesGenerator.unique().aromatic();
//        System.out.println("SMILES :" + aromatic.create(comparison.getFirstAtomMapping().getCommonFragmentInQuery()));
//        System.out.println("SMILES :" + aromatic.create(comparison.getFirstAtomMapping().getCommonFragmentInTarget()));

        Assert.assertEquals(0.8235, comparison.getTanimotoSimilarity(), .001);
        Assert.assertEquals(2, comparison.getAllAtomMapping().size());
    }

    /**
     * Test Tanimoto NAD+ & NADH for Bond InSensitive. Time taking cases
     *
     * @throws Exception
     */
    @Test
    public void testNADPlusNADHBondInSensitive() throws Exception {
//        //System.out.println("getTanimoto for NAD+ and NADH");
        SmilesParser smilesParser = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer query = smilesParser.parseSmiles("NC(=O)c1ccc[n+](c1)[C@@H]1O[C@H](COP(O)(=O)OP(O)(=O)OC[C@H]2O[C@H]([C@H](O)[C@@H]2O)n2cnc3c(N)ncnc23)[C@@H](O)[C@H]1O");
        IAtomContainer target = smilesParser.parseSmiles("NC(=O)C1=CN(C=CC1)[C@@H]1O[C@H](COP(O)(=O)OP(O)(=O)OC[C@H]2O[C@H]([C@H](O)[C@@H]2O)n2cnc3c(N)ncnc23)[C@@H](O)[C@H]1O");

        IAtomContainer ac1 = AtomContainerManipulator.removeHydrogens(query);
        IAtomContainer ac2 = AtomContainerManipulator.removeHydrogens(target);
        double score = 1.0;
        Isomorphism smsd1 = new Isomorphism(ac1, ac2, Algorithm.DEFAULT, false, false, false);
        smsd1.setChemFilters(true, true, true);
        Assert.assertEquals(score, smsd1.getTanimotoSimilarity(), 0.001);
    }
}
