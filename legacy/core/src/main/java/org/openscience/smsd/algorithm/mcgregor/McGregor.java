
/* Copyright (C) 2005-2006 Markus Leber
 *               2006-2014 Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>
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
package org.openscience.smsd.algorithm.mcgregor;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Stack;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.isomorphism.matchers.IQueryAtomContainer;
import org.openscience.smsd.helper.BinaryTree;
import org.openscience.smsd.tools.IterationManager;

/**
 * Class which reports MCS solutions based on the McGregor algorithm published
 * in 1982.
 *
 * <p>
 * The SMSD algorithm is described in this paper. <font color="#FF0000">please
 * refer Rahman <i>et.al. 2009</i></font> {
 *
 * @cdk.cite SMSD2009}. </p>
 *
 *
 *
 *
 * @author Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>
 */
public final class McGregor {

    private final boolean shouldMatchRings;
    private final boolean bondMatch;
    private final boolean matchAtomType;
    private IterationManager iterationManager = null;
    private boolean timeout = false;

    /**
     * @return the timeout
     */
    public synchronized boolean isTimeout() {
        return timeout;
    }

    private synchronized boolean checkTimeout() {
        if (getIterationManager().isMaxIteration()) {
            this.timeout = true;
//            System.err.println("McGregor MCS has hit the iteration limits " + getIterationManager().getCounter());
            return true;
        }
        getIterationManager().increment();
        return false;
    }

    /**
     * @return the iterationManager
     */
    public IterationManager getIterationManager() {
        return iterationManager;
    }

    /**
     * @param iterationManager the iterationManager to set
     */
    public void setIterationManager(IterationManager iterationManager) {
        this.iterationManager = iterationManager;
    }
    /*
     *
     * McGregor starts
     */
    private final IAtomContainer target;
    private BinaryTree last = null;
    private BinaryTree first = null;
    private final Stack<List<Integer>> bestARCS;
    private final List<Integer> modifiedARCS;
    private int bestarcsleft;
    private int globalMCSSize;
    private final List<List<Integer>> mappings;
    /*
     * This should be more or equal to all the atom types
     */
    private final String[] SIGNS = {
        "$1", "$2", "$3", "$4", "$5", "$6", "$7", "$8", "$9", "$10", "$11", "$12",
        "$13", "$15", "$16", "$17", "$18", "$19", "$20", "$21", "$22", "$23", "$24",
        "$25", "$26", "$27", "$28", "$29", "$30", "$31", "$32", "$33", "$34", "$35",
        "$36", "$37", "$38", "$39", "$40", "$41", "$42", "$43", "$44", "$45", "$46",
        "$47", "$48", "$49", "$50", "$51", "$52", "$53", "$54", "$55"
    };
    private boolean newMatrix = false;

    /**
     * Constructor for the McGregor algorithm.
     *
     * @param source
     * @param target
     * @param mappings
     * @param shouldMatchBonds
     * @param shouldMatchRings
     * @param matchAtomType
     */
    public McGregor(IAtomContainer source,
            IAtomContainer target,
            List<List<Integer>> mappings,
            boolean shouldMatchBonds,
            boolean shouldMatchRings,
            boolean matchAtomType) {
        this.shouldMatchRings = shouldMatchRings;
        this.bondMatch = shouldMatchBonds;
        this.matchAtomType = matchAtomType;
        this.target = target;
        this.mappings = Collections.synchronizedList(mappings);
        this.bestarcsleft = 0;
        setIterationManager(new IterationManager((source.getAtomCount() + this.target.getAtomCount()) * 1000));

        if (!mappings.isEmpty()) {
            this.globalMCSSize = mappings.get(0).size();
        } else {
            this.globalMCSSize = 0;
        }
//        System.out.println("globalMCSSize " + globalMCSSize);
        this.modifiedARCS = Collections.synchronizedList(new ArrayList<>());
        this.bestARCS = new Stack<>();
        this.newMatrix = false;
    }

    /**
     * Constructor for the McGregor algorithm.
     *
     * @param source
     * @param target
     * @param mappings
     */
    public McGregor(IQueryAtomContainer source, IAtomContainer target, List<List<Integer>> mappings) {
        this.shouldMatchRings = true;
        this.bondMatch = true;
        this.matchAtomType = true;
        this.target = target;
        this.mappings = Collections.synchronizedList(mappings);
        this.bestarcsleft = 0;
        setIterationManager(new IterationManager((source.getAtomCount() + this.target.getAtomCount()) * 1000));

        if (!mappings.isEmpty()) {
            this.globalMCSSize = mappings.get(0).size();
        } else {
            this.globalMCSSize = 0;
        }
        this.modifiedARCS = Collections.synchronizedList(new ArrayList<Integer>());
        this.bestARCS = new Stack<>();
        this.newMatrix = false;
    }

    /**
     * Start McGregor search and extend the mappings if possible.
     *
     * @param source
     * @param largestMappingSize
     * @param present_Mapping
     * @throws IOException
     */
    public synchronized void startMcGregorIteration(IAtomContainer source, int largestMappingSize, Map<Integer, Integer> present_Mapping) throws IOException {

        this.globalMCSSize = (largestMappingSize / 2);
//        System.out.println("globalMCSSize " + globalMCSSize);
        List<String> c_tab1_copy = McGregorChecks.generateCTabCopy(source);
        List<String> c_tab2_copy = McGregorChecks.generateCTabCopy(target);

        //find mapped atoms of both molecules and store these in mappedAtoms
        List<Integer> mapped_atoms = new ArrayList<>();
//        System.out.println("\nMapped Atoms");
        present_Mapping.entrySet().stream().map((map) -> {
            //            System.out.println("i:" + map.getKey() + " j:" + map.getValue());
            mapped_atoms.add(map.getKey());
            return map;
        }).forEachOrdered((map) -> {
            mapped_atoms.add(map.getValue());
        });
        int mapping_size = present_Mapping.size();

        List<Integer> i_bond_neighborsA = new ArrayList<>();
        List<String> c_bond_neighborsA = new ArrayList<>();

        List<Integer> i_bond_setA = new ArrayList<>();
        List<String> c_bond_setA = new ArrayList<>();

        List<Integer> i_bond_neighborsB = new ArrayList<>();
        List<Integer> i_bond_setB = new ArrayList<>();
        List<String> c_bond_neighborsB = new ArrayList<>();
        List<String> c_bond_setB = new ArrayList<>();

        //find unmapped atoms of molecule A
        List<Integer> unmapped_atoms_molA = McGregorChecks.markUnMappedAtoms(true, source, present_Mapping);
        int counter = 0;
        int gSetBondNumA = 0;
        int gSetBondNumB = 0;
        int gNeighborBondnumA = 0; //number of remaining molecule A bonds after the clique search, which are neighbors of the MCS_1
        int gNeighborBondNumB = 0; //number of remaining molecule B bonds after the clique search, which are neighbors of the MCS_1

        QueryProcessor queryProcess = new QueryProcessor(
                c_tab1_copy,
                c_tab2_copy,
                SIGNS,
                gNeighborBondnumA,
                gSetBondNumA,
                i_bond_neighborsA,
                c_bond_neighborsA,
                mapping_size,
                i_bond_setA,
                c_bond_setA);

        if (!(source instanceof IQueryAtomContainer)) {
            queryProcess.process(
                    source,
                    target,
                    unmapped_atoms_molA,
                    mapped_atoms,
                    counter);
        } else {
            queryProcess.process(
                    (IQueryAtomContainer) source,
                    target,
                    unmapped_atoms_molA,
                    mapped_atoms,
                    counter);
        }

        gSetBondNumA = queryProcess.getBondNumA();
        gNeighborBondnumA = queryProcess.getNeighborBondNumA();
        i_bond_neighborsA = queryProcess.getIBondNeighboursA();
        c_bond_neighborsA = queryProcess.getCBondNeighborsA();

        //find unmapped atoms of molecule B
        List<Integer> unmapped_atoms_molB = McGregorChecks.markUnMappedAtoms(false, target, present_Mapping);

//        System.out.println("unmapped_atoms_molB: " + unmapped_atoms_molB.size());
        //Extract bonds which are related with unmapped atoms of molecule B.
        //In case that unmapped atoms are connected with already mapped atoms, the mapped atoms are labelled with
        //new special signs -> the result are two vectors: cBondNeighborsA and int_bonds_molB, which contain those
        //bonds of molecule B, which are relevant for the McGregorBondTypeInSensitive algorithm.
        //The special signs must be transfered to the corresponding atoms of molecule A
        TargetProcessor targetProcess = new TargetProcessor(
                c_tab2_copy,
                SIGNS,
                gNeighborBondNumB,
                gSetBondNumB,
                i_bond_neighborsB,
                c_bond_neighborsB,
                gNeighborBondnumA,
                i_bond_neighborsA,
                c_bond_neighborsA);

        targetProcess.process(
                target,
                unmapped_atoms_molB,
                mapping_size,
                i_bond_setB,
                c_bond_setB,
                mapped_atoms,
                counter);

        gSetBondNumB = targetProcess.getBondNumB();
        gNeighborBondNumB = targetProcess.getNeighborBondNumB();
        i_bond_neighborsB = targetProcess.getIBondNeighboursB();
        c_bond_neighborsB = targetProcess.getCBondNeighborsB();

        boolean dummy = false;

        McgregorHelper mcGregorHelper = new McgregorHelper(dummy,
                present_Mapping.size(),
                mapped_atoms,
                gNeighborBondnumA,
                gNeighborBondNumB,
                i_bond_neighborsA,
                i_bond_neighborsB,
                c_bond_neighborsA,
                c_bond_neighborsB,
                gSetBondNumA,
                gSetBondNumB,
                i_bond_setA,
                i_bond_setB,
                c_bond_setA,
                c_bond_setB);
        iterator(source, mcGregorHelper);
    }

    private synchronized int iterator(IAtomContainer source, McgregorHelper mcGregorHelper) throws IOException {

        boolean mappingCheckFlag = mcGregorHelper.isMappingCheckFlag();
        int mappedAtomCount = mcGregorHelper.getMappedAtomCount();
        List<Integer> mappedAtoms = new ArrayList<>(mcGregorHelper.getMappedAtomsOrg());
        int neighborBondNumA = mcGregorHelper.getNeighborBondNumA();
        int neighborBondNumB = mcGregorHelper.getNeighborBondNumB();

//        //check possible mappings:
        boolean furtherMappingFlag;
        if (source instanceof IQueryAtomContainer) {
            furtherMappingFlag = McGregorChecks.isFurtherMappingPossible((IQueryAtomContainer) source, target, mcGregorHelper, isBondMatch(), isMatchRings(), isMatchAtomType());
        } else {
            furtherMappingFlag = McGregorChecks.isFurtherMappingPossible(source, target, mcGregorHelper, isBondMatch(), isMatchRings(), isMatchAtomType());
        }

        if (neighborBondNumA == 0 || neighborBondNumB == 0 || mappingCheckFlag || !furtherMappingFlag) {
            setFinalMappings(mappedAtoms, mappedAtomCount);
            return 0;
        }

        modifiedARCS.clear();
        int size = neighborBondNumA * neighborBondNumB;
        for (int i = 0; i < size; i++) {
            modifiedARCS.add(i, 0);
        }
        setModifedArcs(source, mcGregorHelper);
        first = new BinaryTree(-1);
        last = first;
        last.setEqual(null);
        last.setNotEqual(null);
        bestarcsleft = 0;

        startsearch(mcGregorHelper);
        Stack<List<Integer>> bestARCSClone = new Stack<>();

        bestARCSClone.addAll(bestARCS);
        while (!bestARCS.empty()) {
            bestARCS.pop();
        }
        searchAndExtendMappings(source, bestARCSClone, mcGregorHelper);

        //System.out.println("In the iterator Termination");
        //System.out.println("============+++++++++==============");
        //System.out.println("Mapped Atoms before iterator Over: " + mappedAtoms);
        return 0;
    }

    private synchronized void searchAndExtendMappings(
            IAtomContainer source,
            Stack<List<Integer>> bestARCSClone,
            McgregorHelper mcGregorHelper) throws IOException {
        int mappedAtomCount = mcGregorHelper.getMappedAtomCount();
        int setNumA = mcGregorHelper.getSetNumA();
        int setNumB = mcGregorHelper.getsetNumB();
        List<Integer> i_bond_setA = mcGregorHelper.getIBondSetA();
        List<Integer> i_bond_setB = mcGregorHelper.getIBondSetB();
        List<String> c_bond_setA = mcGregorHelper.getCBondSetA();
        List<String> c_bond_setB = mcGregorHelper.getCBondSetB();

        while (!bestARCSClone.empty()) {

            List<Integer> MARCS_vector = new ArrayList<>(bestARCSClone.peek());
            List<Integer> new_Mapping = findMcGregorMapping(source, MARCS_vector, mcGregorHelper);

            int newMapingSize = new_Mapping.size() / 2;
            boolean no_further_MAPPINGS = false;
            if (mappedAtomCount == newMapingSize) {
                no_further_MAPPINGS = true;
            }

            List<Integer> new_i_neighborsA = new ArrayList<>(); //instead of iBondNeighborAtomsA
            List<Integer> new_i_neighborsB = new ArrayList<>(); //instead of iBondNeighborAtomsB
            List<String> new_c_neighborsA = new ArrayList<>(); //instead of cBondNeighborsA
            List<String> new_c_neighborsB = new ArrayList<>(); //instead of cBondNeighborsB
            List<Integer> new_i_bond_setA = new ArrayList<>(); //instead of i_bond_setA
            List<Integer> new_i_bond_setB = new ArrayList<>(); //instead of i_bond_setB
            List<String> new_c_bond_setA = new ArrayList<>(); //instead of c_bond_setA
            List<String> new_c_bond_setB = new ArrayList<>(); //instead of c_bond_setB
            //new values for setNumA + setNumB
            //new arrays for i_bond_setA + i_bond_setB + c_bond_setB + c_bond_setB

            List<String> c_setA_copy = McGregorChecks.generateCSetCopy(setNumA, c_bond_setA);
            List<String> c_setB_copy = McGregorChecks.generateCSetCopy(setNumB, c_bond_setB);

            //find unmapped atoms of molecule A
            List<Integer> unmapped_atoms_molA = new ArrayList<>();
            int unmapped_numA = 0;
            boolean atomA_is_unmapped = true;

            for (int a = 0; a < source.getAtomCount(); a++) {
                for (int b = 0; b < newMapingSize; b++) {
                    if (a == new_Mapping.get(b * 2 + 0)) {
                        atomA_is_unmapped = false;
                    }
                }
                if (atomA_is_unmapped) {
                    unmapped_atoms_molA.add(unmapped_numA, a);
                    unmapped_numA += 1;
                }
                atomA_is_unmapped = true;
            }

            //The special signs must be transfered to the corresponding atoms of molecule B
            int counter = 0;
            //number of remaining molecule A bonds after the clique search, which aren't neighbors
            int newSetBondNumA = 0; //instead of setNumA
            int newNeighborNumA = 0; //instead of localNeighborBondnumA

            QueryProcessor queryProcess
                    = new QueryProcessor(
                            c_setA_copy,
                            c_setB_copy,
                            SIGNS,
                            newNeighborNumA,
                            newSetBondNumA,
                            new_i_neighborsA,
                            new_c_neighborsA,
                            newMapingSize,
                            new_i_bond_setA,
                            new_c_bond_setA);

            queryProcess.process(
                    setNumA,
                    setNumB,
                    i_bond_setA,
                    i_bond_setB,
                    unmapped_atoms_molA,
                    new_Mapping,
                    counter);

            newSetBondNumA = queryProcess.getBondNumA();
            newNeighborNumA = queryProcess.getNeighborBondNumA();
            new_i_neighborsA = queryProcess.getIBondNeighboursA();
            new_c_neighborsA = queryProcess.getCBondNeighborsA();

            //find unmapped atoms of molecule B
            List<Integer> unmapped_atoms_molB = new ArrayList<>();
            int unmapped_numB = 0;
            boolean atomB_is_unmapped = true;

            for (int a = 0; a < target.getAtomCount(); a++) {
                for (int b = 0; b < newMapingSize; b++) {
                    if (a == new_Mapping.get(b * 2 + 1)) {
                        atomB_is_unmapped = false;
                    }
                }
                if (atomB_is_unmapped) {
                    unmapped_atoms_molB.add(unmapped_numB, a);
                    unmapped_numB += 1;
                }
                atomB_is_unmapped = true;
            }

            //number of remaining molecule B bonds after the clique search, which aren't neighbors
            int newSetBondNumB = 0; //instead of setNumB
            int newNeighborNumB = 0; //instead of localNeighborBondNumB

            TargetProcessor targetProcess = new TargetProcessor(
                    c_setB_copy,
                    SIGNS,
                    newNeighborNumB,
                    newSetBondNumB,
                    new_i_neighborsB,
                    new_c_neighborsB,
                    newNeighborNumA,
                    new_i_neighborsA,
                    new_c_neighborsA);

            targetProcess.process(
                    setNumB,
                    unmapped_atoms_molB,
                    newMapingSize,
                    i_bond_setB,
                    c_bond_setB,
                    new_Mapping,
                    counter,
                    new_i_bond_setB,
                    new_c_bond_setB);

            newSetBondNumB = targetProcess.getBondNumB();
            newNeighborNumB = targetProcess.getNeighborBondNumB();
            new_i_neighborsB = targetProcess.getIBondNeighboursB();
            new_c_neighborsB = targetProcess.getCBondNeighborsB();

//             System.out.println("Mapped Atoms before Iterator2: " + mappedAtoms);
            McgregorHelper newMH = new McgregorHelper(no_further_MAPPINGS,
                    newMapingSize,
                    new_Mapping,
                    newNeighborNumA,
                    newNeighborNumB,
                    new_i_neighborsA,
                    new_i_neighborsB,
                    new_c_neighborsA,
                    new_c_neighborsB,
                    newSetBondNumA,
                    newSetBondNumB,
                    new_i_bond_setA,
                    new_i_bond_setB,
                    new_c_bond_setA,
                    new_c_bond_setB);

            iterator(source, newMH);
            bestARCSClone.pop();
//            System.out.println("End of the iterator!!!!");
        }
    }

    private synchronized List<Integer> findMcGregorMapping(IAtomContainer source, List<Integer> MARCS, McgregorHelper mcGregorHelper) {

        int neighborBondNumA = mcGregorHelper.getNeighborBondNumA();
        int neighborBondNumB = mcGregorHelper.getNeighborBondNumB();
        List<Integer> currentMapping = new ArrayList<>(mcGregorHelper.getMappedAtomsOrg());
        List<Integer> additional_mapping = new ArrayList<>();

        for (int x = 0; x < neighborBondNumA; x++) {
            for (int y = 0; y < neighborBondNumB; y++) {
                if (MARCS.get(x * neighborBondNumB + y) == 1) {
                    extendMapping(source, x, y, mcGregorHelper, additional_mapping, currentMapping);
                }
            }
        }

        int additionalMappingSize = additional_mapping.size();
        //add McGregorBondTypeInSensitive mapping to the Clique mapping
        for (int a = 0; a < additionalMappingSize; a += 2) {
            currentMapping.add(additional_mapping.get(a + 0));
            currentMapping.add(additional_mapping.get(a + 1));
        }

//        remove recurring mappings from currentMapping
        List<Integer> unique_MAPPING = McGregorChecks.removeRecurringMappings(currentMapping);
        return unique_MAPPING;
    }

    private synchronized void setModifedArcs(IAtomContainer source, McgregorHelper mcGregorHelper) {
        int neighborBondNumA = mcGregorHelper.getNeighborBondNumA();
        int neighborBondNumB = mcGregorHelper.getNeighborBondNumB();
        List<Integer> iBondNeighborAtomsA = mcGregorHelper.getiBondNeighborAtomsA();
        List<Integer> iBondNeighborAtomsB = mcGregorHelper.getiBondNeighborAtomsB();
        List<String> cBondNeighborsA = mcGregorHelper.getcBondNeighborsA();
        List<String> cBondNeighborsB = mcGregorHelper.getcBondNeighborsB();
        for (int row = 0; row < neighborBondNumA; row++) {
            for (int column = 0; column < neighborBondNumB; column++) {

                String G1A = cBondNeighborsA.get(row * 4 + 0);
                String G2A = cBondNeighborsA.get(row * 4 + 1);
                String G1B = cBondNeighborsB.get(column * 4 + 0);
                String G2B = cBondNeighborsB.get(column * 4 + 1);

                if (!(source instanceof IQueryAtomContainer) && matchGAtoms(G1A, G2A, G1B, G2B)) {
                    int Index_I = iBondNeighborAtomsA.get(row * 3 + 0);
                    int Index_IPlus1 = iBondNeighborAtomsA.get(row * 3 + 1);

                    IAtom R1_A = source.getAtom(Index_I);
                    IAtom R2_A = source.getAtom(Index_IPlus1);
                    IBond reactantBond = source.getBond(R1_A, R2_A);

                    int Index_J = iBondNeighborAtomsB.get(column * 3 + 0);
                    int Index_JPlus1 = iBondNeighborAtomsB.get(column * 3 + 1);

                    IAtom P1_B = target.getAtom(Index_J);
                    IAtom P2_B = target.getAtom(Index_JPlus1);
                    IBond productBond = target.getBond(P1_B, P2_B);
                    if (McGregorChecks.isMatchFeasible(reactantBond, productBond, isBondMatch(), isMatchRings(), isMatchAtomType())) {
                        modifiedARCS.set(row * neighborBondNumB + column, 1);
                    }
                } else if (source instanceof IQueryAtomContainer) {
                    int Index_I = iBondNeighborAtomsA.get(row * 3 + 0);
                    int Index_IPlus1 = iBondNeighborAtomsA.get(row * 3 + 1);

                    IAtom R1_A = source.getAtom(Index_I);
                    IAtom R2_A = source.getAtom(Index_IPlus1);
                    IBond reactantBond = source.getBond(R1_A, R2_A);

                    int Index_J = iBondNeighborAtomsB.get(column * 3 + 0);
                    int Index_JPlus1 = iBondNeighborAtomsB.get(column * 3 + 1);

                    IAtom P1_B = target.getAtom(Index_J);
                    IAtom P2_B = target.getAtom(Index_JPlus1);
                    IBond productBond = target.getBond(P1_B, P2_B);
                    if (McGregorChecks.isMatchFeasible(reactantBond, productBond, isBondMatch(), isMatchRings(), isMatchAtomType())) {
                        modifiedARCS.set(row * neighborBondNumB + column, 1);
                    }
                }
            }
        }
    }

    private synchronized void partsearch(int xstart, int ystart, List<Integer> TEMPMARCS_ORG, McgregorHelper mcGregorHelper) {

        if (checkTimeout()) {
            return;
        }

        int neighborBondNumA = mcGregorHelper.getNeighborBondNumA();
        int neighborBondNumB = mcGregorHelper.getNeighborBondNumB();

        int xIndex = xstart;
        int yIndex = ystart;

        List<Integer> TEMPMARCS = new ArrayList<>(TEMPMARCS_ORG);

        if (TEMPMARCS.get(xstart * neighborBondNumB + ystart) == 1) {

            McGregorChecks.removeRedundantArcs(xstart, ystart, TEMPMARCS, mcGregorHelper);
            int arcsleft = McGregorChecks.countArcsLeft(TEMPMARCS, neighborBondNumA, neighborBondNumB);

            //test best arcs left and skip rest if needed
            if (arcsleft >= bestarcsleft) {
                setArcs(xIndex, yIndex, arcsleft, TEMPMARCS, mcGregorHelper);
            }
        } else {
            do {
                yIndex++;
                if (yIndex == neighborBondNumB) {
                    yIndex = 0;
                    xIndex++;
                }

            } while ((xIndex < neighborBondNumA) && (TEMPMARCS.get(xIndex * neighborBondNumB + yIndex) != 1)); //Correction by ASAD set value minus 1

            if (xIndex < neighborBondNumA) {

                partsearch(xIndex, yIndex, TEMPMARCS, mcGregorHelper);
                TEMPMARCS.set(xIndex * neighborBondNumB + yIndex, 0);
                partsearch(xIndex, yIndex, TEMPMARCS, mcGregorHelper);
            } else {
                int arcsleft = McGregorChecks.countArcsLeft(TEMPMARCS, neighborBondNumA, neighborBondNumB);
                if (arcsleft >= bestarcsleft) {
                    popBestArcs(arcsleft);
                    if (checkMARCS(TEMPMARCS, neighborBondNumA, neighborBondNumB)) {
                        bestARCS.push(TEMPMARCS);
                    }
                }
            }
        }
    }

//The function is called in function partsearch. The function is given indexZ temporary matrix.
//The function checks whether the temporary matrix is already found by calling the function
//"verifyNodes". If the matrix already exists the function returns false which means that
//the matrix will not be stored. Otherwise the function returns true which means that the
//matrix will be stored in function partsearch.
    private synchronized boolean checkMARCS(List<Integer> MARCS_T, int neighborBondNumA, int neighborBondNumB) {

        int size = neighborBondNumA * neighborBondNumA;
        List<Integer> posnum_list = new ArrayList<>(size);

        for (int i = 0; i < posnum_list.size(); i++) {
            posnum_list.add(i, 0);
        }

        int yCounter = 0;
        int count_entries = 0;
        for (int x = 0; x < (neighborBondNumA * neighborBondNumB); x++) {
            if (MARCS_T.get(x) == 1) {
                posnum_list.add(yCounter, x);
                yCounter += 1;
                count_entries++;
            }
        }
        boolean flag = false;

        verifyNodes(posnum_list, first, 0, count_entries);
        if (isNewMatrix()) {
            flag = true;
        }

        return flag;

    }

    private synchronized boolean verifyNodes(List<Integer> matrix, BinaryTree currentStructure, int index, int fieldLength) {
        if (index < fieldLength) {
            if (matrix.get(index) == currentStructure.getValue() && currentStructure.getEqual() != null) {
                setNewMatrix(false);
                verifyNodes(matrix, currentStructure.getEqual(), index + 1, fieldLength);
            }
            if (matrix.get(index) != currentStructure.getValue()) {
                if (currentStructure.getNotEqual() != null) {
                    verifyNodes(matrix, currentStructure.getNotEqual(), index, fieldLength);
                }

                if (currentStructure.getNotEqual() == null) {
                    currentStructure.setNotEqual(new BinaryTree(matrix.get(index)));
                    currentStructure.getNotEqual().setNotEqual(null);
                    int yIndex = 0;

                    BinaryTree last_one = currentStructure.getNotEqual();

                    while ((yIndex + index + 1) < fieldLength) {
                        last_one.setEqual(new BinaryTree(matrix.get(yIndex + index + 1)));
                        last_one = last_one.getEqual();
                        last_one.setNotEqual(null);
                        yIndex++;

                    }
                    last_one.setEqual(null);
                    setNewMatrix(true);
                }

            }
        }
        return true;
    }

    private synchronized void startsearch(McgregorHelper mcGregorHelper) {

        int neighborBondNumA = mcGregorHelper.getNeighborBondNumA();
        int neighborBondNumB = mcGregorHelper.getNeighborBondNumB();

        int size = neighborBondNumA * neighborBondNumB;
        List<Integer> FIXARCS = new ArrayList<>(size);//  Initialize FIXARCS with 0
        for (int i = 0; i < size; i++) {
            FIXARCS.add(i, 0);
        }

        int xIndex = 0;
        int yIndex = 0;

        while ((xIndex < neighborBondNumA) && (modifiedARCS.get(xIndex * neighborBondNumB + yIndex) != 1)) {
            yIndex++;
            if (yIndex == neighborBondNumB) {
                yIndex = 0;
                xIndex++;
            }
        }

        if (xIndex == neighborBondNumA) {
            yIndex = neighborBondNumB - 1;
            xIndex -= 1;
        }

        if (modifiedARCS.get(xIndex * neighborBondNumB + yIndex) == 0) {
            partsearch(xIndex, yIndex, modifiedARCS, mcGregorHelper);
        }

        if (modifiedARCS.get(xIndex * neighborBondNumB + yIndex) != 0) {
            partsearch(xIndex, yIndex, modifiedARCS, mcGregorHelper);
            modifiedARCS.set(xIndex * neighborBondNumB + yIndex, 0);
            partsearch(xIndex, yIndex, modifiedARCS, mcGregorHelper);
        }

    }

    /**
     * Returns computed mappings.
     *
     * @return mappings
     */
    public synchronized List<List<Integer>> getMappings() {
        return Collections.synchronizedList(mappings);
    }

    /**
     * Returns MCS size.
     *
     * @return MCS size
     */
    public synchronized int getMCSSize() {
        return this.globalMCSSize;
    }

    private synchronized void setFinalMappings(List<Integer> mapped_atoms, int mappedAtomCount) {
        if (mappedAtomCount >= globalMCSSize) {
//            System.out.println("Hello-1");
            if (mappedAtomCount > globalMCSSize) {
//                System.out.println("Hello-2");
                this.globalMCSSize = mappedAtomCount;
//                System.out.println("best_MAPPING_size: " + globalMCSSize);
                mappings.clear();
            }
            mappings.add(mapped_atoms);
//            System.out.println("mappings " + mappings);
        }
    }

    private synchronized void setArcs(int xIndex, int yIndex, int arcsleft, List<Integer> TEMPMARCS, McgregorHelper mcGregorHelper) {
        int neighborBondNumA = mcGregorHelper.getNeighborBondNumA();
        int neighborBondNumB = mcGregorHelper.getNeighborBondNumB();
        do {
            yIndex += 1;
            if (yIndex == neighborBondNumB) {
                yIndex = 0;
                xIndex += 1;

            }
        } //Correction by ASAD set value minus 1
        while ((xIndex < neighborBondNumA) && (TEMPMARCS.get(xIndex * neighborBondNumB + yIndex) != 1));
        if (xIndex < neighborBondNumA) {

            partsearch(xIndex, yIndex, TEMPMARCS, mcGregorHelper);
            TEMPMARCS.set(xIndex * neighborBondNumB + yIndex, 0);
            partsearch(xIndex, yIndex, TEMPMARCS, mcGregorHelper);

        } else {
            popBestArcs(arcsleft);
            if (checkMARCS(TEMPMARCS, neighborBondNumA, neighborBondNumB)) {
                bestARCS.push(TEMPMARCS);
            }
        }
    }

    private synchronized void popBestArcs(int arcsleft) {
        if (arcsleft > bestarcsleft) {
            McGregorChecks.removeTreeStructure(first);
            first = new BinaryTree(-1);
            last = new BinaryTree(-1);
            last.setEqual(null);
            last.setNotEqual(null);
            while (!bestARCS.empty()) {
                bestARCS.pop();
            }
        }
        bestarcsleft = arcsleft;
    }

    private synchronized void extendMapping(IAtomContainer source, int xIndex, int yIndex, McgregorHelper mcGregorHelper, List<Integer> additional_mapping, List<Integer> currentMapping) {

        int Atom1_moleculeA = mcGregorHelper.getiBondNeighborAtomsA().get(xIndex * 3 + 0);
        int Atom2_moleculeA = mcGregorHelper.getiBondNeighborAtomsA().get(xIndex * 3 + 1);
        int Atom1_moleculeB = mcGregorHelper.getiBondNeighborAtomsB().get(yIndex * 3 + 0);
        int Atom2_moleculeB = mcGregorHelper.getiBondNeighborAtomsB().get(yIndex * 3 + 1);

        IAtom R1_A = source.getAtom(Atom1_moleculeA);
        IAtom R2_A = source.getAtom(Atom2_moleculeA);
        IBond reactantBond = source.getBond(R1_A, R2_A);

        IAtom P1_B = target.getAtom(Atom1_moleculeB);
        IAtom P2_B = target.getAtom(Atom2_moleculeB);
        IBond productBond = target.getBond(P1_B, P2_B);

//      Bond Order Check Introduced by Asad
        if (McGregorChecks.isMatchFeasible(reactantBond, productBond, isBondMatch(), isMatchRings(), isMatchAtomType())) {

            for (int indexZ = 0; indexZ < mcGregorHelper.getMappedAtomCount(); indexZ++) {

                int Mapped_Atom_1 = currentMapping.get(indexZ * 2 + 0);
                int Mapped_Atom_2 = currentMapping.get(indexZ * 2 + 1);

                if ((Mapped_Atom_1 == Atom1_moleculeA) && (Mapped_Atom_2 == Atom1_moleculeB)) {
                    additional_mapping.add(Atom2_moleculeA);
                    additional_mapping.add(Atom2_moleculeB);
                } else if ((Mapped_Atom_1 == Atom1_moleculeA) && (Mapped_Atom_2 == Atom2_moleculeB)) {
                    additional_mapping.add(Atom2_moleculeA);
                    additional_mapping.add(Atom1_moleculeB);
                } else if ((Mapped_Atom_1 == Atom2_moleculeA) && (Mapped_Atom_2 == Atom1_moleculeB)) {
                    additional_mapping.add(Atom1_moleculeA);
                    additional_mapping.add(Atom2_moleculeB);
                } else if ((Mapped_Atom_1 == Atom2_moleculeA) && (Mapped_Atom_2 == Atom2_moleculeB)) {
                    additional_mapping.add(Atom1_moleculeA);
                    additional_mapping.add(Atom1_moleculeB);
                }
            }//for loop
        }
    }

    private synchronized boolean matchGAtoms(String G1A, String G2A, String G1B, String G2B) {
        return (G1A.compareToIgnoreCase(G1B) == 0
                && G2A.compareToIgnoreCase(G2B) == 0)
                || (G1A.compareToIgnoreCase(G2B) == 0
                && G2A.compareToIgnoreCase(G1B) == 0);
    }

    /**
     * Checks if its a new Matrix.
     *
     * @return the newMatrix
     */
    public synchronized boolean isNewMatrix() {
        return newMatrix;
    }

    /**
     * set a new Matrix.
     *
     * @param newMatrix the newMatrix to set
     */
    public synchronized void setNewMatrix(boolean newMatrix) {
        this.newMatrix = newMatrix;
    }

    /**
     * Should bonds match
     *
     * @return the bondMatch
     */
    private synchronized boolean isBondMatch() {
        return bondMatch;
    }

    /**
     * @return the shouldMatchRings
     */
    public boolean isMatchRings() {
        return shouldMatchRings;
    }

    /**
     * @return the matchAtomType
     */
    public boolean isMatchAtomType() {
        return matchAtomType;
    }
}
