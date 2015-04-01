/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.openscience.smsd.mcss;

import java.util.ArrayList;
import java.util.Calendar;
import java.util.Collection;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;
import java.util.concurrent.RecursiveTask;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.tools.ILoggingTool;
import org.openscience.cdk.tools.LoggingToolFactory;
import org.openscience.smsd.AtomAtomMapping;
import org.openscience.smsd.Isomorphism;
import org.openscience.smsd.interfaces.Algorithm;

/**
 *
 * @author Asad
 */
public class ForkAndJoinMCSS extends RecursiveTask<List<Set<Fragment>>> {

    public static final long THRESHOLD = 2;
    private final static ILoggingTool logger
            = LoggingToolFactory.createLoggingTool(ForkAndJoinMCSS.class);
    private static final long serialVersionUID = 1322393982398293L;
    private final List<IAtomContainer> calculateMCSS;
    private IAtomContainer root;
    private final boolean matchBonds;
    private final boolean matchRings;
    private final boolean matchAtomType;
    private final int start;
    private final int end;

    /**
     *
     * @param root
     * @param calculateMCSS
     * @param matchBonds
     * @param matchRings
     * @param matchAtomType
     */
    public ForkAndJoinMCSS(IAtomContainer root, List<IAtomContainer> calculateMCSS, boolean matchBonds, boolean matchRings, boolean matchAtomType) {
        this(root, calculateMCSS, matchBonds, matchRings, matchAtomType, 0, calculateMCSS.size());
    }

    private ForkAndJoinMCSS(IAtomContainer root, List<IAtomContainer> calculateMCSS, boolean matchBonds, boolean matchRings, boolean matchAtomType, int start, int end) {
        this.root = root;
        this.calculateMCSS = calculateMCSS;
        this.matchBonds = matchBonds;
        this.matchRings = matchRings;
        this.matchAtomType = matchAtomType;
        this.start = start;
        this.end = end;
    }

    @Override
    protected List<Set<Fragment>> compute() {
        int length = end - start;
        if (length <= THRESHOLD) {
            Set<Fragment> mcs = ComputeSimilarity(root, calculateMCSS, matchBonds, matchRings, matchAtomType, start, end);
            ArrayList<Set<Fragment>> arrayList = new ArrayList<>();
            arrayList.add(mcs);
            return arrayList;
        }
        ForkAndJoinMCSS leftTask = new ForkAndJoinMCSS(root, calculateMCSS, matchBonds, matchRings, matchAtomType, start, start + length / 2);
        leftTask.fork();
        ForkAndJoinMCSS rightTask = new ForkAndJoinMCSS(root, calculateMCSS, matchBonds, matchRings, matchAtomType, start + length / 2, end);

        List<Set<Fragment>> rightResult = rightTask.compute();
        List<Set<Fragment>> leftResult = leftTask.join();

        List<Set<Fragment>> joinArrayList = new ArrayList<>(rightResult);
        joinArrayList.addAll(leftResult);
        return joinArrayList;
    }

    /*
     * MULTIPLE Fragments of MCS are returned if present
     */
    private synchronized Set<Fragment>
            ComputeSimilarity(IAtomContainer root,
                    List<IAtomContainer> jobs,
                    boolean matchBonds,
                    boolean matchRings,
                    boolean matchAtomType,
                    int start,
                    int end) {
        /*
         * Store final solution here
         */

        logger.debug("Calling MCSSTask " + start + " with " + end + " items");
        long startTime = Calendar.getInstance().getTimeInMillis();
        IAtomContainer query = root;
        long calcTime = startTime;

        /*
         * Local Seeds
         */
        Set<Fragment> localSeeds = new TreeSet<>();
        int minSeedSize = query.getAtomCount();

        for (int j = start; j < end; j++) {
            IAtomContainer target = jobs.get(j);
            Collection<Fragment> fragmentsFromMCS;
            Isomorphism comparison = null;

            try {
                comparison = new Isomorphism(query, target, Algorithm.DEFAULT, matchBonds, matchRings, matchAtomType);
                //comparison.setChemFilters(true, true, true);
            } catch (Exception ex) {
                Logger.getLogger(ForkAndJoinMCSS.class.getName()).log(Level.SEVERE, null, ex);
            }
            fragmentsFromMCS = getMCSS(comparison);
            long endCalcTime = Calendar.getInstance().getTimeInMillis();
            calcTime = endCalcTime;
            if (fragmentsFromMCS.isEmpty()) {
                localSeeds.clear();
                break;
            }
            Iterator<Fragment> iterator = fragmentsFromMCS.iterator();
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
        Set<Fragment> seeds = new HashSet<>();
        if (!localSeeds.isEmpty()) {
            seeds.addAll(localSeeds);
            localSeeds.clear();
        }

        logger.debug("No of Potential MULTIPLE " + seeds.size());
//        System.out.println("No of Potential MULTIPLE " + seeds.size());
        long endTime = Calendar.getInstance().getTimeInMillis();
        logger.debug("Done: task " + start + " took " + (endTime - startTime) + "ms");
        logger.debug(" and mcss has " + query.getAtomCount() + " atoms, and " + query.getBondCount() + " bonds");
        return seeds;
    }

    private synchronized Set<Fragment> getMCSS(Isomorphism comparison) {
        Set<Fragment> matchList = new HashSet<>();
        for (AtomAtomMapping mapping : comparison.getAllAtomMapping()) {
            IAtomContainer match;
            try {
                match = mapping.getCommonFragment();
                try {
                    Fragment fragment = new Fragment(match);
                    if (!matchList.contains(fragment)) {
                        matchList.add(fragment);
                    }
                } catch (CDKException ex) {
                    logger.error("ERROR IN MCS Thread: ", ex);
                }
            } catch (CloneNotSupportedException ex) {
                logger.error("ERROR IN MCS Thread: ", ex);
            }
        }
        return matchList;
    }

}
