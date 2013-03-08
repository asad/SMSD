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
import java.util.concurrent.ConcurrentLinkedQueue;
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
        /*
         * Store final solution here
         */
        List<IAtomContainer> mcss = new ArrayList<IAtomContainer>();

//        System.out.println("Calling MCSSTask " + taskNumber + " with " + mcssList.size() + " items");
//        long startTime = Calendar.getInstance().getTimeInMillis();
        IAtomContainer querySeed = AtomContainerManipulator.removeHydrogens(mcssList.get(0));
//        long calcTime = startTime;

        ConcurrentLinkedQueue<IAtomContainer> seeds = new ConcurrentLinkedQueue<IAtomContainer>();
        try {
            /*
             * Local Seeds
             */
            Set<Fragment> localSeeds = new TreeSet<Fragment>();
            int minSeedSize = querySeed.getAtomCount();

            for (int index = 1; index < mcssList.size(); index++) {
                IAtomContainer target = mcssList.get(index);
                Collection<Fragment> fragmentsFromMCS;
                BaseMapping comparison;
                if (this.jobType.equals(JobType.MCS)) {
                    comparison = new Isomorphism(querySeed, target, Algorithm.DEFAULT, matchBonds, matchRings);
                    comparison.setChemFilters(true, true, true);
                    fragmentsFromMCS = getMCSS(comparison);
                } else {
                    comparison = new Substructure(querySeed, target, matchBonds, matchRings, false);
                    comparison.setChemFilters(true, true, true);
                    fragmentsFromMCS = getMCSS(comparison);
                }
//                System.out.println("comparison for task " + taskNumber + " has " + fragmentsFromMCS.size()
//                        + " unique matches of size " + comparison.getFirstAtomMapping().getCount());
//                System.out.println("MCSS for task " + taskNumber + " has " + querySeed.getAtomCount() + " atoms, and " + querySeed.getBondCount() + " bonds");
//                System.out.println("Target for task " + taskNumber + " has " + target.getAtomCount() + " atoms, and " + target.getBondCount() + " bonds");
//                long endCalcTime = Calendar.getInstance().getTimeInMillis();
//                System.out.println("Task " + taskNumber + " index " + index + " took " + (endCalcTime - calcTime) + "ms");
//                calcTime = endCalcTime;

                if (fragmentsFromMCS == null || fragmentsFromMCS.isEmpty()) {
                    localSeeds.clear();
                    break;
                }
                Iterator<Fragment> iterator = fragmentsFromMCS.iterator();
                /*
                 * Store rest of the unique hits
                 */
                while (iterator.hasNext()) {
                    Fragment fragment = iterator.next();
                    if (minSeedSize > fragment.getContainer().getAtomCount()) {
                        localSeeds.clear();
                        minSeedSize = fragment.getContainer().getAtomCount();
                    }
                    if (minSeedSize == fragment.getContainer().getAtomCount()) {
                        localSeeds.add(fragment);
                    }
                }
            }
            /*
             * Add all the Maximum Unique Substructures
             */
            if (!localSeeds.isEmpty()) {
                for (Fragment f : localSeeds) {
                    seeds.add(f.getContainer());
                }
                localSeeds.clear();
            }

//        System.out.println("No of Potential MCS " + seeds.size());

            /*
             * Choose only cleaned MCS Substructures
             */
            minSeedSize = Integer.MAX_VALUE;

            while (!seeds.isEmpty()) {
                IAtomContainer fragmentMCS = seeds.poll();
                localSeeds = new TreeSet<Fragment>();
//            System.out.println("Potential MCS " + getMCSSSmiles(fragmentMCS));
                Collection<Fragment> fragmentsFromMCS;
                for (int index = 0; index < mcssList.size(); index++) {
                    IAtomContainer target = mcssList.get(index);
                    Isomorphism comparison = new Isomorphism(fragmentMCS, target, Algorithm.DEFAULT, matchBonds, matchRings);
                    comparison.setChemFilters(true, true, true);
                    fragmentsFromMCS = getMCSS(comparison);

                    /*
                     * Only true MCSS is added
                     */
                    if (fragmentsFromMCS == null || fragmentsFromMCS.isEmpty()) {
                        localSeeds.clear();
                        break;
                    }
                    Iterator<Fragment> iterator = fragmentsFromMCS.iterator();
                    /*
                     * Store rest of the unique hits
                     */
                    while (iterator.hasNext()) {
                        Fragment fragment = iterator.next();
                        if (minSeedSize > fragment.getContainer().getAtomCount()) {
                            localSeeds.clear();
                            minSeedSize = fragment.getContainer().getAtomCount();
                        }
                        if (minSeedSize == fragment.getContainer().getAtomCount()) {
                            localSeeds.add(fragment);
                        }
                    }
                    /*
                     * Top solution
                     */
                    fragmentMCS = localSeeds.iterator().next().getContainer();
                }

                /*
                 * Add all the Maximum Unique Substructures
                 */
                if (!localSeeds.isEmpty()) {
                    for (Fragment f : localSeeds) {
                        mcss.add(f.getContainer());
                    }
                    localSeeds.clear();
                }

            }
        } catch (Exception e) {
            logger.error("ERROR IN MCS Thread: ", e);
        }
//        long endTime = Calendar.getInstance().getTimeInMillis();
//        System.out.println("Done: task " + taskNumber + " took " + (endTime - startTime) + "ms");
//        System.out.println(" and mcss has " + querySeed.getAtomCount() + " atoms, and " + querySeed.getBondCount() + " bonds");
        return mcss;
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
