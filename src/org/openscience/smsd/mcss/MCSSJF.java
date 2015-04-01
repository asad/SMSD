/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.openscience.smsd.mcss;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.concurrent.ForkJoinPool;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.tools.ILoggingTool;
import org.openscience.cdk.tools.LoggingToolFactory;
import org.openscience.smsd.Isomorphism;
import org.openscience.smsd.interfaces.Algorithm;

/**
 *
 * @author Asad
 */
public class MCSSJF {

    private final static ILoggingTool logger
            = LoggingToolFactory.createLoggingTool(MCSSJF.class);
    private final boolean matchBonds;
    private final boolean matchRings;
    private final boolean matchAtomType;
    private final List<IAtomContainer> solutions;

    private synchronized List<IAtomContainer> submitFJJob(List<IAtomContainer> mcssList) {
        final ForkJoinPool forkJoinPool = new ForkJoinPool();
        ForkAndJoinMCSS mcssJobThread = new ForkAndJoinMCSS(mcssList.get(0), new ArrayList<>(mcssList), matchBonds, matchRings, matchAtomType);
        List<Set<Fragment>> solutionsJF = forkJoinPool.invoke(mcssJobThread);
        forkJoinPool.shutdown();

        System.out.println("First Round " + solutionsJF.size());

        /*
         Test all the solutions from the first round
         */
        for (Set<Fragment> l : solutionsJF) {
            if (l.isEmpty()) {
                /*
                 This means atleast one of the solutions had no MCS
                 Hence we need not look further
                 */
                return new ArrayList<>();
            }
        }

        /*
         Find the root list as there maybe more than one unique solution
         */
        Set<Fragment> rootSet = new HashSet<>();
        if (!solutionsJF.isEmpty()) {
            rootSet.addAll(solutionsJF.get(0));
        }

        List<IAtomContainer> uniqueSolutions = new ArrayList<>();


        /*
         Now compute MCS for all the unique solutions
         */
        for (Fragment root : rootSet) {
            try {
                System.out.println("sub-root: " + root.getFragmentSMILES());
            } catch (CDKException ex) {
                Logger.getLogger(MCSSJF.class.getName()).log(Level.SEVERE, null, ex);
            }
            int solutionSize = Integer.MAX_VALUE;
            IAtomContainer mcs = null;
            for (IAtomContainer target : mcssList) {
                Isomorphism s = new Isomorphism(root.getContainer(), target, Algorithm.VFLibMCS, matchBonds, matchRings, matchAtomType);
                if (s.getMappingCount() > 0) {
                    int count = s.getFirstAtomMapping().getCount();
                    if (count < solutionSize) {
                        solutionSize = count;
                        try {
                            mcs = s.getFirstAtomMapping().getCommonFragment();
                        } catch (CloneNotSupportedException ex) {
                            Logger.getLogger(MCSSJF.class.getName()).log(Level.SEVERE, null, ex);
                        }
                    }
                } else {
                    solutionSize = 0;
                }

                try {
                    System.out.println("Q: " + root.getFragmentSMILES()
                            + ", mcs " + s.getFirstAtomMapping().getCommonFragmentAsSMILES()
                            + ", target " + new SmilesGenerator().create(target));
                } catch (CDKException | CloneNotSupportedException ex) {
                    Logger.getLogger(MCSSJF.class.getName()).log(Level.SEVERE, null, ex);
                }
            }
            if (solutionSize > 0 && mcs != null) {
                uniqueSolutions.add(mcs);
            }
        }

        System.out.println("Second Round " + uniqueSolutions.size());
        return uniqueSolutions;
    }

    /**
     *
     * @param calculateMCSS
     * @param matchBonds
     * @param matchRings
     * @param matchAtomType
     */
    public MCSSJF(List<IAtomContainer> calculateMCSS, boolean matchBonds, boolean matchRings, boolean matchAtomType) {
        this.matchBonds = matchBonds;
        this.matchRings = matchRings;
        this.matchAtomType = matchAtomType;
        this.solutions = submitFJJob(calculateMCSS);
    }

    public synchronized String getTitle() {
        return "Calculating Maximum Commmon Substrutures (MCSS) using SMSD";
    }

    /**
     * @return the solutions
     */
    public Collection<IAtomContainer> getSolutions() {
        return solutions;
    }
}
