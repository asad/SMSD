/* 
 * Copyright (C) 2009-2014 Syed Asad Rahman <asad@ebi.ac.uk>
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
package org.openscience.smsd.algorithm.vflib;

import java.util.List;
import java.util.Map;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Assert;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.isomorphism.matchers.IQueryAtomContainer;
import org.openscience.cdk.isomorphism.matchers.QueryAtomContainerCreator;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.smsd.algorithm.vflib.vf2.Pattern;
import org.openscience.smsd.algorithm.vflib.vf2.VF;

/**
 * Unit testing for the {@link VF} class.
 *
 * @author Syed Asad Rahman
 * @author egonw
 * @cdk.module test-smsd
 */
public class VFTest {

    public VFTest() {
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

    @Test
    public void testQueryAtomContainerDefault() throws CDKException {
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer query = sp.parseSmiles("CC");
        IAtomContainer target = sp.parseSmiles("C1CCC12CCCC2");
        Pattern findIdentical = VF.findSubstructure(query, true, true, true);

        boolean foundMatches = findIdentical.matches(target);
        Assert.assertTrue(foundMatches);

        IQueryAtomContainer queryContainer = QueryAtomContainerCreator.createSymbolAndBondOrderQueryContainer(query);
        findIdentical = VF.findSubstructure(queryContainer);
        foundMatches = findIdentical.matches(target);
        Assert.assertTrue(foundMatches);
    }

    @Test
    /*
     * Imp test case
     */
    public void testMatchAtomCount() throws CDKException {
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer query = sp.parseSmiles("CC");
        IAtomContainer target = sp.parseSmiles("C1CCC12CCCC2");
        int i = 1;
        for (IAtom a : query.atoms()) {
            a.setID(String.valueOf(i));
            //System.out.print(i + ",");
            i++;
        }
        //////System.out.println();
        i = 1;
        for (IAtom a : target.atoms()) {
            a.setID(String.valueOf(i));
            //System.out.print(i + ",");
            i++;
        }

        Pattern findIdentical = VF.findSubstructure(query, true, true, true);
        List<Map<IAtom, IAtom>> matchAll = findIdentical.matchAll(target);
        Assert.assertTrue(!matchAll.isEmpty());
        Assert.assertEquals(18, matchAll.size());

        IQueryAtomContainer queryContainer = QueryAtomContainerCreator.createSymbolAndBondOrderQueryContainer(query);
        findIdentical = VF.findSubstructure(queryContainer);
        matchAll = findIdentical.matchAll(target);
        Assert.assertEquals(18, matchAll.size());
    }

    @Test
    public void testIsSubgraph() throws CDKException {
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer query = sp.parseSmiles("[#8]P([#8])(=O)[#8]-[#6@H]-1-[#6@H](-[#8]P([#8])([#8])=O)-[#6@@H](-[#8]P([#8])([#8])=O)-[#6@H](-[#8]P([#8])([#8])=O)-[#6@H](-[#8]P([#8])([#8])=O)-[#6@@H]-1-[#8]P([#8])([#8])=O");
        IAtomContainer target = sp.parseSmiles("[#8]P([#8])(=O)[#8]-[#6@@H]-1-[#6@H](-[#8]P([#8])([#8])=O)-[#6@@H](-[#8]P([#8])([#8])=O)-[#6@H](-[#8]P([#8])(=O)[#8]P([#8])([#8])=O)-[#6@@H](-[#8]P([#8])([#8])=O)-[#6@@H]-1-[#8]P([#8])([#8])=O");
        int i = 1;
        for (IAtom a : query.atoms()) {
            a.setID(String.valueOf(i));
            //System.out.print(i + ",");
            i++;
        }
        //////System.out.println();
        i = 1;
        for (IAtom a : target.atoms()) {
            a.setID(String.valueOf(i));
            //System.out.print(i + ",");
            i++;
        }

        //////System.out.println();
//        Pattern pattern = Ullmann.findSubstructure(query);
        Pattern pattern = VF.findSubstructure(query, true, true, true);
        int hits = 0;
        List<Map<IAtom, IAtom>> matchAll = pattern.matchAll(target);
        for (Map<IAtom, IAtom> mapping : matchAll) {
            hits++;
        }
        System.out.println("HITS " + hits);

        Assert.assertEquals(768, matchAll.size());

        IQueryAtomContainer queryContainer = QueryAtomContainerCreator.createSymbolAndBondOrderQueryContainer(query);
        pattern = VF.findSubstructure(queryContainer);
        Assert.assertTrue(pattern.matches(target));
    }

    @Test
    /*
     * Imp test case
     */
    public void testSeeds() throws CDKException {
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer query = sp.parseSmiles("CCCCCn1c2c(cccc2)c(c1)C(=O)c3ccc(c4c3cccc4)Cl");
        IAtomContainer target = sp.parseSmiles("CCCCCn1c2c(cccc2)c(c1)C(=O)c3cccc4c3cccc4Cl");
        int i = 1;
        for (IAtom a : query.atoms()) {
            a.setID(String.valueOf(i));
            //System.out.print(i + ",");
            i++;
        }
        //////System.out.println();
        i = 1;
        for (IAtom a : target.atoms()) {
            a.setID(String.valueOf(i));
            //System.out.print(i + ",");
            i++;
        }

        Pattern find = VF.findSeeds(query, true, true, true);
        List<Map<IAtom, IAtom>> matchAll = find.matchAll(target);
        Assert.assertTrue(!matchAll.isEmpty());
        Assert.assertEquals(186, matchAll.size());

        IQueryAtomContainer queryContainer = QueryAtomContainerCreator.createSymbolAndBondOrderQueryContainer(query);
        find = VF.findSeeds(queryContainer);
        matchAll = find.matchAll(target);
        Assert.assertEquals(246, matchAll.size());
    }

}
