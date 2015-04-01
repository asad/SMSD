/**
 * Copyright (C) 2009-2013 Syed Asad Rahman <asad@ebi.ac.uk>
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
 * FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
 * details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA.
 */
package tools;

import java.util.ArrayList;
import java.util.Calendar;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.junit.Assert;
import org.junit.Test;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.smsd.mcss.MCSSJF;

/**
 *
 * @author Syed Asad Rahman <asad@ebi.ac.uk>
 *
 */
public class TestMCSSJF {

    @Test
    public void case1() {
        HashMap<String, String> map = new HashMap<>();
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        //Data

        map.put("30331701", "Cc1ccc(C2N(C(c3c(Cl)cccc3)=O)c(c4N=C(C2=C5O)CC(c6c(Cl)cccc6)C5)cccc4)cc1");
        map.put("4170207", "CCC(N1C(c2ccc(F)cc2)C(C3=Nc(c14)cccc4)=C(O)CC(c5ccc(C(C)(C)C)cc5)C3)=O");
        map.put("6368123", "CC(N1C(c2ccc(F)cc2)C(C3=Nc(c14)cccc4)=C(O)CC(c5ccc(Cl)cc5)C3)=O");
        map.put("MCS", "O=C(N4c5ccccc5(N=C2C(=C(O)CC(c1ccccc1)C2)C4(c3ccccc3)))C");
        map.put("MMV011436", "CC(c1ccc(C2CC(C3=C(O)C2)=Nc(c4N(C(c5c(Cl)cccc5)=O)C3c6ccc(F)cc6)cccc4)cc1)(C)C");
        map.put("4196500", "CC(c1ccc(C2CC(C3=C(O)C2)=Nc(c4N(C(c5c(Cl)cccc5)=O)C3c6ccc(F)cc6)cccc4)cc1)(C)C");

       List<IAtomContainer> jobs = new ArrayList<>();
        for (String s : map.values()) {
            try {
                jobs.add(sp.parseSmiles(s));
            } catch (InvalidSmilesException ex) {
                Logger.getLogger(TestMCSSJF.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
        long startTime = Calendar.getInstance().getTimeInMillis();
        MCSSJF mcss = new MCSSJF(jobs, true, true, true);
        for (IAtomContainer ac : mcss.getSolutions()) {
            ////System.out.println("Result MCS " + getMCSSSmiles(ac));
            Assert.assertEquals(31, ac.getAtomCount());
        }
        long endCalcTime = Calendar.getInstance().getTimeInMillis();

        ////System.out.println("Total time: " + (endCalcTime - startTime) + "ms");
    }

    @Test
    public void case2() {
        HashMap<String, String> map = new HashMap<>();
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        // Data

        map.put("1", "CCCCC=NNN");
        map.put("2", "CCCC-NNN");
        map.put("3", "CCCSNNN");
        map.put("4", "CCCC=NNN");

        List<IAtomContainer> jobs = new ArrayList<>();
        for (String s : map.values()) {
            try {
                jobs.add(sp.parseSmiles(s));
            } catch (InvalidSmilesException ex) {
                Logger.getLogger(TestMCSSJF.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
        long startTime = Calendar.getInstance().getTimeInMillis();
        MCSSJF mcss = new MCSSJF(jobs, true, true, true);
        for (IAtomContainer ac : mcss.getSolutions()) {
            ////System.out.println("Result MCS " + getMCSSSmiles(ac));
            Assert.assertEquals(3, ac.getAtomCount());
        }
        long endCalcTime = Calendar.getInstance().getTimeInMillis();

        ////System.out.println("Total time: " + (endCalcTime - startTime) + "ms");
    }

    @Test
    public void case3() {
        HashMap<String, String> map = new HashMap<>();
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        // Data

        map.put("1", "CCC=NNN");
        map.put("2", "CCC-NNN");
        map.put("3", "CCCSNNN");

        List<IAtomContainer> jobs = new ArrayList<>();
        for (String s : map.values()) {
            try {
                jobs.add(sp.parseSmiles(s));
            } catch (InvalidSmilesException ex) {
                Logger.getLogger(TestMCSSJF.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
        long startTime = Calendar.getInstance().getTimeInMillis();
        MCSSJF mcss = new MCSSJF(jobs, true, true, false);
        for (IAtomContainer ac : mcss.getSolutions()) {
            //System.out.println("Result MCS " + getMCSSSmiles(ac));
            Assert.assertEquals(3, ac.getAtomCount());
        }
        long endCalcTime = Calendar.getInstance().getTimeInMillis();
        ////System.out.println("Total time: " + (endCalcTime - startTime) + "ms");
    }

    @Test
    public void case4() {
        HashMap<String, String> map = new HashMap<>();
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        // Data

        map.put("1", "CCC=NNN");
        map.put("2", "CCC-NNN");
        map.put("3", "CCCSNNN");
        map.put("4", "SSSCCC");
        map.put("5", "CSSSCSNNN");

        List<IAtomContainer> jobs = new ArrayList<>();
        for (String s : map.values()) {
            try {
                jobs.add(sp.parseSmiles(s));
            } catch (InvalidSmilesException ex) {
                Logger.getLogger(TestMCSSJF.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
        long startTime = Calendar.getInstance().getTimeInMillis();
        MCSSJF mcss = new MCSSJF(jobs, true, true, true);
        Collection<IAtomContainer> solutions = mcss.getSolutions();
        System.out.println("Number of Solutions: " + solutions.size());
        for (IAtomContainer ac : solutions) {
            try {
                System.out.println("Result MCS " + new SmilesGenerator().create(ac));
            } catch (CDKException ex) {
                Logger.getLogger(TestMCSSJF.class.getName()).log(Level.SEVERE, null, ex);
            }
            Assert.assertEquals(1, ac.getAtomCount());
        }
        long endCalcTime = Calendar.getInstance().getTimeInMillis();

        ////System.out.println("Total time: " + (endCalcTime - startTime) + "ms");
    }

    @Test
    public void case5() {
        HashMap<String, String> map = new HashMap<>();
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        // Data

        map.put("1", "O=C(OC1OC(CO)C(O)C(O)C1(O))Cc3c[nH]c2ccccc23");
        map.put("2", "O=C(OC1C(O)C(O)C(O)C(O)C1(O))Cc3c[nH]c2ccccc23");

        List<IAtomContainer> jobs = new ArrayList<>();
        for (String s : map.values()) {
            try {
                IAtomContainer parseSmiles = sp.parseSmiles(s);
                System.out.println(s + ", size:" + parseSmiles.getAtomCount());
                jobs.add(parseSmiles);
            } catch (InvalidSmilesException ex) {
                Logger.getLogger(TestMCSSJF.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
        long startTime = Calendar.getInstance().getTimeInMillis();
        MCSSJF mcss = new MCSSJF(jobs, true, true, false);
        for (IAtomContainer ac : mcss.getSolutions()) {
            try {
                System.out.println("Result MCS " + new SmilesGenerator().create(ac));
            } catch (CDKException ex) {
                Logger.getLogger(TestMCSSJF.class.getName()).log(Level.SEVERE, null, ex);
            }
            Assert.assertEquals(21, ac.getAtomCount());
        }
        long endCalcTime = Calendar.getInstance().getTimeInMillis();

        ////System.out.println("Total time: " + (endCalcTime - startTime) + "ms");
    }
}
