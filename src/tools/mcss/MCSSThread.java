/**
 * Copyright (C) 2009-2013 Syed Asad Rahman <asad@ebi.ac.uk>
 *
 * Contact: cdk-devel@lists.sourceforge.net
 *
 * This program is free software; you can redistribute it and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either version 2.1 of the License, or (at your option)
 * any later version. All we ask is that proper credit is given for our work, which includes - but is not limited to -
 * adding the above copyright notice to the beginning of your source code files, and to any copyright notice that you
 * may distribute with programs based on this work.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
 * warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
 * details.
 *
 * You should have received a copy of the GNU Lesser General Public License along with this program; if not, write to
 * the Free Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA.
 */
package tools.mcss;

import java.util.*;
import java.util.concurrent.Callable;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import org.openscience.smsd.AtomAtomMapping;
import org.openscience.smsd.BaseMapping;
import org.openscience.smsd.Isomorphism;
import org.openscience.smsd.Substructure;
import org.openscience.smsd.interfaces.Algorithm;

/**
 *
 * @author Syed Asad Rahman <asad@ebi.ac.uk>
 *
 */
final public class MCSSThread implements Callable<List<IAtomContainer>> {

    private final List<IAtomContainer> mcssList;
    private final JobType jobType;
    private int taskNumber = 0;

    /**
     *
     * @param mcssList
     * @param jobType MCS/Substructure
     * @param taskNumber
     */
    public MCSSThread(List<IAtomContainer> mcssList, JobType jobType, int taskNumber) {
        this.mcssList = mcssList;
        this.jobType = jobType;
        this.taskNumber = taskNumber;
    }

    @Override
    public synchronized List<IAtomContainer> call() {

//        System.out.println("Calling MCSSTask " + taskNumber + " with " + mcssList.size() + " items");
        List<IAtomContainer> resultsList = new ArrayList<IAtomContainer>();
        long startTime = Calendar.getInstance().getTimeInMillis();
        IAtomContainer querySeed = AtomContainerManipulator.removeHydrogens(mcssList.get(0));
        long calcTime = startTime;


        try {
            for (int index = 1; index < mcssList.size(); index++) {
                IAtomContainer target = AtomContainerManipulator.removeHydrogens(mcssList.get(index));
                Collection<Fragment> fragmentsFomMCS;
                if (this.jobType == JobType.MCS) {
                    Isomorphism comparison = new Isomorphism(querySeed, target, Algorithm.DEFAULT, true, true);
                    comparison.setChemFilters(true, true, true);
                    fragmentsFomMCS = getMCSS(comparison);
                } else {
                    Substructure comparison = new Substructure(querySeed, target, true, true, false);
                    comparison.setChemFilters(true, true, true);
                    fragmentsFomMCS = getMCSS(comparison);
                }
//                System.out.println("MCSS for task " + taskNumber + " has " + querySeed.getAtomCount() + " atoms, and " + querySeed.getBondCount() + " bonds");
//                System.out.println("Target for task " + taskNumber + " has " + target.getAtomCount() + " atoms, and " + target.getBondCount() + " bonds");


                long endCalcTime = Calendar.getInstance().getTimeInMillis();
//                System.out.println("Task " + taskNumber + " index " + index + " took " + (endCalcTime - calcTime) + "ms");
                calcTime = endCalcTime;

                if (querySeed == null || querySeed.isEmpty()) {
                    break;
                }
//                System.out.println("comparison for task " + taskNumber + " has " + fragmentsFomMCS.size()
//                        + " unique matches of size " + comparison.getFirstAtomMapping().getCount());
                querySeed = fragmentsFomMCS.iterator().next().getContainer();
            }

        } catch (Exception e) {
            Logger.getLogger(MCSSThread.class.getName()).log(Level.SEVERE, null, e);
        }
        if (resultsList != null) {
            resultsList.add(querySeed);
        }

        long endTime = Calendar.getInstance().getTimeInMillis();
//        System.out.print("Done: task " + taskNumber + " took " + (endTime - startTime) + "ms");
//        System.out.println(" and mcss has " + querySeed.getAtomCount() + " atoms, and " + querySeed.getBondCount() + " bonds");
        return resultsList;
    }

    private synchronized Collection<Fragment> getMCSS(BaseMapping comparison) {
        Set<Fragment> matchList = new HashSet<Fragment>();
        // System.out.println("Found "+comparison.getAllAtomMapping().size()+" mappings");
        for (AtomAtomMapping mapping : comparison.getAllAtomMapping()) {
            //IAtomContainer match = getMatchedSubgraph(mol1, mapping.getMappings());
            IAtomContainer match;
            try {
                match = mapping.getCommonFragmentInQuery();
                try {
                    matchList.add(new Fragment(match));
                } catch (CDKException ex) {
                    Logger.getLogger(MCSSThread.class.getName()).log(Level.SEVERE, null, ex);
                }
            } catch (CloneNotSupportedException ex) {
                Logger.getLogger(MCSSThread.class.getName()).log(Level.SEVERE, null, ex);
            }
            // System.out.println("match has "+match.getAtomCount()+" atoms, and "+match.getBondCount()+" bonds");
        }

        return matchList;
    }

    /**
     * Return SMILES
     *
     * @param ac
     * @return
     */
    public synchronized String getMCSSSmiles(IAtomContainer ac) {
        SmilesGenerator g = new SmilesGenerator();
        g.setUseAromaticityFlag(true);
        return g.createSMILES(ac);
    }
}
