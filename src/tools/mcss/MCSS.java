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

import cmd.AtomContainerComparator;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.openscience.cdk.interfaces.IAtomContainer;

/**
 *
 * @author Syed Asad Rahman <asad@ebi.ac.uk>
 *
 */
final public class MCSS {

    private final List<IAtomContainer> calculateMCSS;

    public MCSS(List<IAtomContainer> jobList, JobType jobType, int numberOfThreads) {
        int threadsAvailable = Runtime.getRuntime().availableProcessors() - 1;
        if (numberOfThreads > 0) {
            threadsAvailable = numberOfThreads;
        }
        Comparator<IAtomContainer> comparator = new AtomContainerComparator();
        Collections.sort(jobList, comparator);
        calculateMCSS = calculateMCSS(jobList, jobType, threadsAvailable);
    }

    private synchronized List<IAtomContainer> calculateMCSS(List<IAtomContainer> mcssList, JobType jobType, int nThreads) {
        List<IAtomContainer> newMCSSList = Collections.synchronizedList(new ArrayList<IAtomContainer>(nThreads));
        if (nThreads == 1) {
            MCSSThread task = new MCSSThread(mcssList, jobType, 1);
            List<IAtomContainer> results = task.call();
            if (results != null) {
                newMCSSList.addAll(results);
            }
            return newMCSSList;
        } else {
            /*
             * Calling recursive MCS
             */
            newMCSSList = submitMultiThreadedJob(mcssList, jobType, nThreads);
            while (newMCSSList.size() > 1) {
                if (newMCSSList.size() > 2) {
                    newMCSSList = submitMultiThreadedJob(newMCSSList, jobType, nThreads);
                } else {
                    newMCSSList = submitMultiThreadedJob(newMCSSList, jobType, 1);
                }
            }
        }
        return newMCSSList;
    }

    /**
     * @return the calculateMCSS
     */
    public synchronized List<IAtomContainer> getCalculateMCSS() {
        return Collections.unmodifiableList(calculateMCSS);
    }

    private synchronized List<IAtomContainer> submitMultiThreadedJob(List<IAtomContainer> mcssList, JobType jobType, int nThreads) {
        int taskNumber = 1;
        List<IAtomContainer> newMCSSList = Collections.synchronizedList(new ArrayList<IAtomContainer>(nThreads));
        List<Future<List<IAtomContainer>>> futureList = new ArrayList<Future<List<IAtomContainer>>>();
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
                MCSSThread mcssJobThread = new MCSSThread(subList, jobType, taskNumber++);
                Future<List<IAtomContainer>> callMCSSThread = threadPool.submit(mcssJobThread);
                futureList.add(callMCSSThread);
            } else {
                newMCSSList.add(subList.get(0));
            }
        }

        for (Future<List<IAtomContainer>> results : futureList) {

            if (results == null) {
                continue;
            }
            try {
                List<IAtomContainer> f = results.get();
                if (f != null) {
                    newMCSSList.addAll(results.get());
                }
            } catch (Exception e) {
                Logger.getLogger(MCSSThread.class.getName()).log(Level.WARNING, "Execution exception: {0}", e);
            }

        }

        threadPool.shutdown();
        return newMCSSList;
    }

    public synchronized String getTitle() {
        return "Calculating Maximum Commmon Substrutures (MCSS)";
    }
}
