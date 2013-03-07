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
package org.openscience.smsd.mcss;

import java.util.*;
import java.util.concurrent.Callable;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.tools.ILoggingTool;
import org.openscience.cdk.tools.LoggingToolFactory;
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

    private final static ILoggingTool logger =
            LoggingToolFactory.createLoggingTool(MCSSThread.class);
    private final List<IAtomContainer> mcssList;
    private final JobType jobType;
    private final int taskNumber;
    private final boolean matchBonds;
    private final boolean matchRings;

    /**
     *
     * @param mcssList
     * @param jobType MCS/Substructure
     * @param taskNumber
     */
    public MCSSThread(List<IAtomContainer> mcssList, JobType jobType, int taskNumber) {
        this(mcssList, jobType, taskNumber, true, true);
    }

    MCSSThread(List<IAtomContainer> mcssList, JobType jobType, int taskNumber, boolean matchBonds, boolean matchRings) {
        this.mcssList = mcssList;
        this.jobType = jobType;
        this.taskNumber = taskNumber;
        this.matchBonds = matchBonds;
        this.matchRings = matchRings;
    }

    @Override
    public synchronized List<IAtomContainer> call() {

//        System.out.println("Calling MCSSTask " + taskNumber + " with " + mcssList.size() + " items");
        List<IAtomContainer> resultsList = new ArrayList<IAtomContainer>();
//        long startTime = Calendar.getInstance().getTimeInMillis();
        IAtomContainer querySeed = AtomContainerManipulator.removeHydrogens(mcssList.get(0));
//        long calcTime = startTime;


        try {
            for (int index = 1; index < mcssList.size(); index++) {
                IAtomContainer target = AtomContainerManipulator.removeHydrogens(mcssList.get(index));
                Collection<Fragment> fragmentsFomMCS;
                BaseMapping comparison;
                if (this.jobType.equals(JobType.MCS)) {
                    comparison = new Isomorphism(querySeed, target, Algorithm.DEFAULT, matchBonds, matchRings);
                    comparison.setChemFilters(true, true, true);
                    fragmentsFomMCS = getMCSS(comparison);
                    querySeed = null;
                } else {
                    comparison = new Substructure(querySeed, target, matchBonds, matchRings, false);
                    comparison.setChemFilters(true, true, true);
                    fragmentsFomMCS = getMCSS(comparison);
                    querySeed = null;
                }
//                System.out.println("comparison for task " + taskNumber + " has " + fragmentsFomMCS.size()
//                        + " unique matches of size " + comparison.getFirstAtomMapping().getCount());
//                System.out.println("MCSS for task " + taskNumber + " has " + querySeed.getAtomCount() + " atoms, and " + querySeed.getBondCount() + " bonds");
//                System.out.println("Target for task " + taskNumber + " has " + target.getAtomCount() + " atoms, and " + target.getBondCount() + " bonds");
//                long endCalcTime = Calendar.getInstance().getTimeInMillis();
//                System.out.println("Task " + taskNumber + " index " + index + " took " + (endCalcTime - calcTime) + "ms");
//                calcTime = endCalcTime;

                if (fragmentsFomMCS == null || fragmentsFomMCS.isEmpty()) {
                    break;
                }
                querySeed = fragmentsFomMCS.iterator().next().getContainer();
            }

        } catch (Exception e) {
            logger.error("ERROR IN MCS Thread: ", e);
        }
        if (querySeed != null) {
            resultsList.add(querySeed);
        }

//        long endTime = Calendar.getInstance().getTimeInMillis();
//        System.out.println("Done: task " + taskNumber + " took " + (endTime - startTime) + "ms");
//        System.out.println(" and mcss has " + querySeed.getAtomCount() + " atoms, and " + querySeed.getBondCount() + " bonds");
        return resultsList;
    }

    private synchronized Collection<Fragment> getMCSS(BaseMapping comparison) {
        Set<Fragment> matchList = new HashSet<Fragment>();
        for (AtomAtomMapping mapping : comparison.getAllAtomMapping()) {
            IAtomContainer match;
            try {
                match = mapping.getCommonFragmentInQuery();
                try {
                    matchList.add(new Fragment(match));
                } catch (CDKException ex) {
                    logger.error("ERROR IN MCS Thread: ", ex);
                }
            } catch (CloneNotSupportedException ex) {
                logger.error("ERROR IN MCS Thread: ", ex);
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

    /**
     * @return the taskNumber
     */
    public int getTaskNumber() {
        return taskNumber;
    }
}
