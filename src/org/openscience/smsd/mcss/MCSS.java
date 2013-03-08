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
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.LinkedBlockingQueue;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.tools.ILoggingTool;
import org.openscience.cdk.tools.LoggingToolFactory;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import org.openscience.smsd.tools.AtomContainerComparator;

/**
 *
 * @author Syed Asad Rahman <asad@ebi.ac.uk>
 *
 */
final public class MCSS {

    private final static ILoggingTool logger =
            LoggingToolFactory.createLoggingTool(MCSS.class);
    private final Collection<IAtomContainer> calculateMCSS;
    private final boolean matchBonds;
    private final boolean matchRings;

    /**
     *
     * @param jobList
     * @param jobType
     * @param numberOfThreads
     */
    public MCSS(List<IAtomContainer> jobList, JobType jobType, int numberOfThreads) {
        this(jobList, jobType, numberOfThreads, true, true);
    }

    /**
     *
     * @param jobList
     * @param jobType
     * @param numberOfThreads
     * @param matchBonds
     * @param matchRings
     */
    public MCSS(List<IAtomContainer> jobList, JobType jobType, int numberOfThreads, boolean matchBonds, boolean matchRings) {
        int threadsAvailable = Runtime.getRuntime().availableProcessors() - 1;

        logger.debug("Demand threads: " + numberOfThreads);
        logger.debug(", Available threads: " + threadsAvailable);
        if (numberOfThreads > 0 && threadsAvailable >= numberOfThreads) {
            threadsAvailable = numberOfThreads;
        } else if (threadsAvailable <= 0) {
            threadsAvailable = 1;
        }
        logger.debug(", Assigned threads: " + threadsAvailable + "\n");
        /*
         * Remove hydrogen from the molecules
         **/
        List<IAtomContainer> selectedJobs = new ArrayList<IAtomContainer>(jobList.size());
        for (IAtomContainer ac : jobList) {
            selectedJobs.add(AtomContainerManipulator.removeHydrogens(ac));
        }
        /*
         * Sort the molecules in the ascending order of atom size and atom type
         */
        Comparator<IAtomContainer> comparator = new AtomContainerComparator();
        Collections.sort(selectedJobs, comparator);
        this.matchBonds = matchBonds;
        this.matchRings = matchRings;
        /*
         * Call the MCS
         */
        calculateMCSS = calculateMCSS(selectedJobs, jobType, threadsAvailable);
        selectedJobs.clear();
    }

    private synchronized Collection<IAtomContainer> calculateMCSS(List<IAtomContainer> mcssList, JobType jobType, int nThreads) {
        List<IAtomContainer> newMCSSList;
        if (nThreads == 1) {
            newMCSSList = new LinkedList<IAtomContainer>(submitSingleThreadedJob(mcssList, jobType, nThreads));
        } else {
            /*
             * Calling recursive MCS
             */
            newMCSSList = new LinkedList<IAtomContainer>(submitMultiThreadedJob(mcssList, jobType, nThreads));
            while (newMCSSList.size() > 1) {
                if (newMCSSList.size() > 2) {
                    newMCSSList = new LinkedList<IAtomContainer>(submitMultiThreadedJob(newMCSSList, jobType, nThreads));
                } else {
                    newMCSSList = new LinkedList<IAtomContainer>(submitMultiThreadedJob(newMCSSList, jobType, 1));
                }
            }
        }
        if (mcssList.get(mcssList.size() - 1) == newMCSSList.get(0)) {
            return new LinkedBlockingQueue<IAtomContainer>();
        }
        return newMCSSList;
    }

    /**
     * @return the calculateMCSS
     */
    public synchronized Collection<IAtomContainer> getCalculateMCSS() {
        return Collections.unmodifiableCollection(calculateMCSS);
    }

    private synchronized LinkedBlockingQueue<IAtomContainer> submitSingleThreadedJob(List<IAtomContainer> mcssList, JobType jobType, int nThreads) {
        LinkedBlockingQueue<IAtomContainer> solutions = new LinkedBlockingQueue<IAtomContainer>();
        MCSSThread task = new MCSSThread(mcssList, jobType, 1);
        LinkedBlockingQueue<IAtomContainer> results = task.call();
        if (results != null) {
            solutions.addAll(results);
        }
        return solutions;
    }

    private synchronized LinkedBlockingQueue<IAtomContainer> submitMultiThreadedJob(List<IAtomContainer> mcssList, JobType jobType, int nThreads) {
        int taskNumber = 1;
        LinkedBlockingQueue<IAtomContainer> solutions = new LinkedBlockingQueue<IAtomContainer>();
        LinkedBlockingQueue<Callable<LinkedBlockingQueue<IAtomContainer>>> callablesQueue = new LinkedBlockingQueue<Callable<LinkedBlockingQueue<IAtomContainer>>>();
        ExecutorService threadPool = Executors.newFixedThreadPool(nThreads);
        int step = (int) Math.ceil((double) mcssList.size() / (double) nThreads);
        if (step < 2) {
            step = 2; // Can't have a step size of less than 2
        }
        for (int i = 0; i < mcssList.size(); i += step) {
            int endPoint = i + step;
            if (endPoint > mcssList.size()) {
                endPoint = mcssList.size();
            }
            List<IAtomContainer> subList = new ArrayList<IAtomContainer>(mcssList.subList(i, endPoint));
            if (subList.size() > 1) {
                MCSSThread mcssJobThread = new MCSSThread(subList, jobType, taskNumber, matchBonds, matchRings);
                callablesQueue.add(mcssJobThread);
                taskNumber++;
            } else {
                solutions.add(subList.get(0));
            }
        }
        try {
            /*
             * Wait all the threads to finish
             */
            List<Future<LinkedBlockingQueue<IAtomContainer>>> futureList = threadPool.invokeAll(callablesQueue);
            /*
             * Collect the results
             */
            for (Iterator<Future<LinkedBlockingQueue<IAtomContainer>>> it = futureList.iterator(); it.hasNext();) {
                Future<LinkedBlockingQueue<IAtomContainer>> callable = it.next();
                LinkedBlockingQueue<IAtomContainer> mapping = callable.get();
                if (callable.isDone() && mapping != null) {
                    solutions.addAll(mapping);
                } else {
                    logger.warn("WARNING: InComplete job in AtomMappingTool: ");
                }
            }
            threadPool.shutdown();
            // Wait until all threads are finish
            while (!threadPool.isTerminated()) {
            }
            System.gc();
        } catch (Exception e) {
            logger.debug("ERROR: in AtomMappingTool: " + e.getMessage());
            logger.error(e);
        } finally {
            threadPool.shutdown();
        }

        return solutions;
    }

    public synchronized String getTitle() {
        return "Calculating Maximum Commmon Substrutures (MCSS) using SMSD";
    }
}
