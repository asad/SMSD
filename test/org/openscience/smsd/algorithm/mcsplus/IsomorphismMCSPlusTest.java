
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
package org.openscience.smsd.algorithm.mcsplus;

import org.openscience.smsd.*;
import java.util.logging.Level;
import java.util.logging.Logger;
import junit.framework.Assert;
import org.junit.*;
import org.openscience.cdk.AtomContainer;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.isomorphism.matchers.IQueryAtomContainer;
import org.openscience.cdk.isomorphism.matchers.QueryAtomContainerCreator;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import org.openscience.smsd.interfaces.Algorithm;

/**
 * Unit testing for the {@link Isomorphism} class.
 *
 * @author Syed Asad Rahman <asad@ebi.ac.uk>
 *
 * @cdk.module test-smsd
 * @cdk.require java1.6+
 */
public class IsomorphismMCSPlusTest extends ImageUtility {

    /**
     * Tests if the CDKMCS can be instantiated without throwing exceptions.
     */
    @Test
    public void IsomorphismTest() {
        Assert.assertNotNull(
                new Isomorphism(new AtomContainer(), new AtomContainer(),
                        Algorithm.DEFAULT,
                        true, false, false));
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
            Isomorphism smsd1 = new Isomorphism(query, target, Algorithm.MCSPlus, true, false, false);
            smsd1.setChemFilters(true, true, true);
            Assert.assertEquals(7, smsd1.getFirstAtomMapping().getCount());
            Assert.assertEquals(2, smsd1.getAllAtomMapping().size());
            Assert.assertNotNull(smsd1.getFirstAtomMapping());
            double score = 0.35;
            Assert.assertEquals(score, smsd1.getTanimotoSimilarity(), 0);
        } catch (InvalidSmilesException ex) {
            Logger.getLogger(MCSPlusHandlerTest.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

    @Test
    public void testQueryAtomContainerMCSPLUS() throws CDKException {
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer query = sp.parseSmiles("CC");
        IAtomContainer target = sp.parseSmiles("C1CCC12CCCC2");
        Isomorphism smsd = new Isomorphism(query, target, Algorithm.MCSPlus, true, false, false);

        Assert.assertEquals(18, smsd.getAllAtomMapping().size());
        Assert.assertTrue(smsd.isSubgraph());

        IQueryAtomContainer queryContainer = QueryAtomContainerCreator.createSymbolAndBondOrderQueryContainer(query);
        smsd = new Isomorphism(queryContainer, target, Algorithm.MCSPlus);
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

        double score = 0.7959;
//        IAtomContainer molecule1 = AtomContainerManipulator.removeHydrogens(molecule1);
//        IAtomContainer molecule2 = AtomContainerManipulator.removeHydrogens(molecule2);
        Isomorphism smsd1 = new Isomorphism(molecule1, molecule2, Algorithm.MCSPlus, true, false, false);
        smsd1.setChemFilters(true, true, true);
        Assert.assertEquals(score, smsd1.getTanimotoSimilarity(), 0.001);
    }

    /**
     * Test Tanimoto NAD+ & NADH for Bond InSensitive.
     *
     * @throws Exception
     */
    @Test
    public void testNADPlusNADHBondInSensitive() throws Exception {
//        //System.out.println("getTanimoto for NAD+ and NADH");
        SmilesParser smilesParser = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer molecule1 = smilesParser.parseSmiles("NC(=O)c1ccc[n+](c1)[C@@H]1O[C@H](COP(O)(=O)OP(O)(=O)OC[C@H]2O[C@H]([C@H](O)[C@@H]2O)n2cnc3c(N)ncnc23)[C@@H](O)[C@H]1O");
        IAtomContainer molecule2 = smilesParser.parseSmiles("NC(=O)C1=CN(C=CC1)[C@@H]1O[C@H](COP(O)(=O)OP(O)(=O)OC[C@H]2O[C@H]([C@H](O)[C@@H]2O)n2cnc3c(N)ncnc23)[C@@H](O)[C@H]1O");

        double score = 1.0;

        IAtomContainer ac1 = AtomContainerManipulator.removeHydrogens(molecule1);
        IAtomContainer ac2 = AtomContainerManipulator.removeHydrogens(molecule2);
        Isomorphism smsd1 = new Isomorphism(ac1, ac2, Algorithm.MCSPlus, false, false, false);
        smsd1.setChemFilters(true, true, true);
        Assert.assertEquals(score, smsd1.getTanimotoSimilarity(), 0.001);
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

//        for (AtomAtomMapping m : comparison.getAllAtomMapping()) {
//            SmilesGenerator sg = new SmilesGenerator();
//            String createSMILESQ = sg.createSMILES(m.getMapCommonFragmentOnQuery());
//            String createSMILEST = sg.createSMILES(m.getMapCommonFragmentOnTarget());
//            System.out.println("createSMILES " + createSMILESQ + "." + createSMILEST);
//            System.out.println("Map " + m.getMappingsByIndex());
//        }
//        RenderedImage generateImage = generateImage(query, target, comparison);
//        boolean write = ImageIO.write(generateImage, "png", new File("MCSPLUS_C1=CC=CC=C1.C1=CC2=C(C=C1)C=CC=C2.png"));
        Assert.assertEquals(0.6, comparison.getTanimotoSimilarity());
        Assert.assertEquals(12, comparison.getAllAtomMapping().size());
    }

    /**
     * Test ring match using MCS VF2Plus
     *
     * @throws Exception
     */
    @Test
    public void testOpenRing() throws Exception {
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        // Benzene
        IAtomContainer query = sp.parseSmiles("O=P(O)(O)OCC1OC(O)(COP(=O)(O)O)C(O)C1(O)");
        // Napthalene
        IAtomContainer target = sp.parseSmiles("O=C(CO)COP(=O)(O)O");
        IAtomContainer ac1 = AtomContainerManipulator.removeHydrogens(query);
        IAtomContainer ac2 = AtomContainerManipulator.removeHydrogens(target);
        Isomorphism comparison = new Isomorphism(ac1, ac2, Algorithm.MCSPlus, false, false, false);
        System.out.println("SMILES  " + comparison.getFirstAtomMapping().getCommonFragmentAsSMILES());
        // set chemical filter true
        comparison.setChemFilters(true, true, true);
//        RenderedImage generateImage = generateImage(ac1, ac2, comparison);
//        boolean write = ImageIO.write(generateImage, "png", new File("MCSPLUS_Ring_Open.png"));
        Assert.assertEquals(0.5, comparison.getTanimotoSimilarity(), .09);
        Assert.assertEquals(4, comparison.getAllAtomMapping().size());
    }

    /**
     * Test ring match using MCS VF2Plus
     *
     * @throws Exception
     */
    @Test
    public void testC06006_C14463() throws Exception {
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        // C06006
        IAtomContainer query = sp.parseSmiles("O=C(O)C(O)(C(=O)C)CC");

        // C14463
        IAtomContainer target = sp.parseSmiles("O=C(O)C(=O)C(O)(C)CC");

        IAtomContainer ac1 = AtomContainerManipulator.removeHydrogens(query);
        IAtomContainer ac2 = AtomContainerManipulator.removeHydrogens(target);
        Isomorphism comparison = new Isomorphism(ac1, ac2, Algorithm.MCSPlus, false, false, false);
        // set chemical filter true
        comparison.setChemFilters(false, true, true);

        System.out.println("SMILES  " + comparison.getFirstAtomMapping().getCommonFragmentAsSMILES());
//        RenderedImage generateImage = generateImage(ac1, ac2, comparison);
//        boolean write = ImageIO.write(generateImage, "png", new File("MCSPLUS_C06006_C14463.png"));
//        Assert.assertEquals(1.0, comparison.getTanimotoSimilarity(), .09);
        Assert.assertEquals(2, comparison.getAllAtomMapping().size());
    }

    /**
     * Test ring match using MCS VF2Plus Exhaustive
     *
     * @throws Exception
     */
    @Test
    public void testNADPNNADPH() throws Exception {
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        // NAD
        IAtomContainer query = sp.parseSmiles("NC(=O)c1ccc[n+](c1)[C@@H]1O[C@H](COP(O)(=O)OP(O)(=O)OC[C@H]2O[C@H]([C@H](O)[C@@H]2O)n2cnc3c(N)ncnc23)[C@@H](O)[C@H]1O");

        // NADH
        IAtomContainer target = sp.parseSmiles("NC(=O)C1=CN(C=CC1)[C@@H]1O[C@H](COP(O)(=O)OP(O)(=O)OC[C@H]2O[C@H]([C@H](O)[C@@H]2O)n2cnc3c(N)ncnc23)[C@@H](O)[C@H]1O");

        IAtomContainer ac1 = AtomContainerManipulator.removeHydrogens(query);
        IAtomContainer ac2 = AtomContainerManipulator.removeHydrogens(target);
        Isomorphism comparison = new Isomorphism(ac1, ac2, Algorithm.MCSPlus, false, false, false);
        // set chemical filter true
        comparison.setChemFilters(true, true, true);
//        RenderedImage generateImage = generateImage(ac1, ac2, comparison);
//        boolean write = ImageIO.write(generateImage, "png", new File("MCSPLUS_NADH_NADP.png"));
        Assert.assertEquals(1.0, comparison.getTanimotoSimilarity(), .09);
        Assert.assertEquals(1, comparison.getAllAtomMapping().size());
    }

}
//c-edges 2288
//d-edges 395230
