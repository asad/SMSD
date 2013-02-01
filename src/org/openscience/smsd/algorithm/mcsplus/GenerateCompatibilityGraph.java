/*
 *
 * Copyright (C) 2009-2013  Syed Asad Rahman <asad@ebi.ebi.ac.uk>
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
 * You should have received iIndex copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA.
 */
package org.openscience.smsd.algorithm.mcsplus;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import org.openscience.cdk.annotations.TestClass;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.isomorphism.matchers.IQueryAtom;
import org.openscience.cdk.isomorphism.matchers.IQueryBond;
import org.openscience.smsd.algorithm.matchers.DefaultMatcher;
import org.openscience.smsd.helper.LabelContainer;

/**
 * This class generates compatibility graph between query and target molecule. It also marks edges in the compatibility
 * graph as c-edges or d-edges. @cdk.module smsd @cdk.githash
 *
 * @author Syed Asad Rahman <asad@ebi.ac.uk>
 */
@TestClass("org.openscience.cdk.smsd.SMSDBondSensitiveTest")
public final class GenerateCompatibilityGraph {

    private List<Integer> compGraphNodes = null;
    private List<Integer> compGraphNodesCZero = null;
    private List<Integer> cEdges = null;
    private List<Integer> dEdges = null;
    private int cEdgesSize = 0;
    private int dEdgesSize = 0;
    private final IAtomContainer source;
    private final IAtomContainer target;
    private final boolean shouldMatchBonds;
    private final boolean shouldMatchRings;

    /**
     * Generates a compatibility graph between two molecules
     *
     * @param source
     * @param target
     * @param shouldMatchBonds
     * @param shouldMatchRings
     * @throws java.io.IOException
     */
    public GenerateCompatibilityGraph(
            IAtomContainer source,
            IAtomContainer target,
            boolean shouldMatchBonds,
            boolean shouldMatchRings) throws IOException {
        this.shouldMatchRings = shouldMatchRings;
        this.shouldMatchBonds = shouldMatchBonds;
        this.source = source;
        this.target = target;
        compGraphNodes = new ArrayList<Integer>();
        compGraphNodesCZero = new ArrayList<Integer>();
        cEdges = Collections.synchronizedList(new ArrayList<Integer>());
        dEdges = Collections.synchronizedList(new ArrayList<Integer>());
        compatibilityGraphNodes();
        compatibilityGraph();


        if (getCEdgesSize() == 0) {
            clearCompGraphNodes();

            clearCEgdes();
            clearDEgdes();

            resetCEdgesSize();
            resetDEdgesSize();

            compatibilityGraphNodesIfCEdgeIsZero();
            compatibilityGraphCEdgeZero();
            clearCompGraphNodesCZero();
        }
    }

    private List<List<Integer>> labelAtoms(IAtomContainer atomCont) {
        List<List<Integer>> label_list = new ArrayList<List<Integer>>();
        LabelContainer labelContainer = LabelContainer.getInstance();

        for (int i = 0; i < atomCont.getAtomCount(); i++) {
            List<Integer> label = new ArrayList<Integer>(7);
            for (int a = 0; a < 7; a++) {
                label.add(a, 0);
            }
            IAtom refAtom = atomCont.getAtom(i);
            String referenceAtom = refAtom.getSymbol();

            label.set(0, labelContainer.getLabelID(referenceAtom));
            List<IAtom> connAtoms = atomCont.getConnectedAtomsList(refAtom);

            int counter = 1;

            for (IAtom negAtom : connAtoms) {
                String neighbouringAtom = negAtom.getSymbol();
                label.set(counter, labelContainer.getLabelID(neighbouringAtom));
                counter += 1;
            }
            bubbleSort(label);
            label_list.add(label);
        }
        return label_list;
    }

    private void bubbleSort(List<Integer> label) {

        boolean flag = true; // set flag to 1 to begin initial pass
        int temp; // holding variable
        for (int i = 0; i < 7 && flag; i++) {
            flag = false;
            for (int j = 0; j < 6; j++) {
                if (label.get(i) > label.get(j + 1)) {
                    // descending order simply changes to >
                    temp = label.get(i); // swap elements

                    label.set(i, label.get(j + 1));
                    label.set(j + 1, temp);
                    flag = true; // indicates that iIndex swap occurred.
                }
            }
        }
    }

    private List<IAtom> reduceAtomSet(IAtomContainer atomCont) {
        List<IAtom> atoms = new ArrayList<IAtom>();
        for (IAtom atom : atomCont.atoms()) {
            atoms.add(atom);
        }
        return atoms;
    }

    /*
     * labelA [0, 2, 2, 2, 0, 0, 0] labelB [0, 2, 2, 2, 0, 0, 0]
     *
     * return true
     *
     * labelA [0, 2, 2, 0, 0, 0, 0] labelB [0, 2, 2, 2, 0, 0, 0]
     *
     * return true
     *
     * labelA [0, 2, 2, 3, 0, 0, 0] labelB [0, 2, 2, 2, 0, 0, 0]
     *
     * return false
     *
     * labelA [0, 2, 2, 2, 0, 0, 0] labelB [0, 2, 2, 0, 0, 0, 0]
     *
     * return false
     *
     */
    private boolean isSubset(List<Integer> labelA, List<Integer> labelB) {
        boolean flag = true;
        for (int i = 0; i < labelA.size(); i++) {
            if (labelA.get(i) != labelB.get(i)) {
                if (labelA.get(i) > 0) {
                    flag = false;
                }
            }
        }
        return flag;
    }

    /**
     * Generate Compatibility Graph Nodes
     *
     * @return
     * @throws IOException
     */
    protected int compatibilityGraphNodes() throws IOException {

        compGraphNodes.clear();
        List<IAtom> sourceAtoms = null;
        List<IAtom> targetAtoms = null;

        sourceAtoms = reduceAtomSet(source);
        targetAtoms = reduceAtomSet(target);

        List<List<Integer>> label_list_molA = labelAtoms(source);
        List<List<Integer>> label_list_molB = labelAtoms(target);

        int sourceNodes = 0;
        int nodeCount = 1;

        for (List<Integer> labelA : label_list_molA) {
            int targetNodes = 0;
            for (List<Integer> labelB : label_list_molB) {
//                if (labelA.equals(labelB)) {
                if (isSubset(labelA, labelB)) {

                    compGraphNodes.add(source.getAtomNumber(sourceAtoms.get(sourceNodes)));
                    compGraphNodes.add(target.getAtomNumber(targetAtoms.get(targetNodes)));
                    compGraphNodes.add(nodeCount);
                    nodeCount += 1;
                }
                targetNodes++;
            }
            sourceNodes++;
        }
        return 0;
    }

    /**
     * Generate Compatibility Graph Nodes Bond Insensitive
     *
     * @return
     * @throws IOException
     */
    protected int compatibilityGraph() throws IOException {

        int comp_graph_nodes_List_size = compGraphNodes.size();

//        System.out.println("Vector_Size: " + compGraphNodes);

//        System.out.println("compGraphNodes size " + comp_graph_nodes_List_size);

        cEdges = new ArrayList<Integer>(); //Initialize the cEdges List
        dEdges = new ArrayList<Integer>(); //Initialize the dEdges List

        for (int a = 0; a < comp_graph_nodes_List_size; a += 3) {
            for (int b = a + 3; b < comp_graph_nodes_List_size; b += 3) {
                if ((a != b)
                        && (compGraphNodes.get(a) != compGraphNodes.get(b))
                        && (compGraphNodes.get(a + 1) != compGraphNodes.get(b + 1))) {

                    boolean molecule1_pair_connected = false;
                    boolean molecule2_pair_connected = false;

                    IBond reactantBond = null;
                    IBond productBond = null;

                    //exists a bond in molecule 2, so that molecule 1 pair is connected?
                    for (int x = 0; x < source.getBondCount(); x++) {

                        reactantBond = source.getBond(x);

                        int q1 = source.getAtomNumber(reactantBond.getAtom(0));
                        int q2 = source.getAtomNumber(reactantBond.getAtom(1));

                        if (((compGraphNodes.get(a) == q1)
                                && (compGraphNodes.get(b) == q2))
                                || ((compGraphNodes.get(a) == q2)
                                && (compGraphNodes.get(b) == q1))) {
                            molecule1_pair_connected = true;
                            break;
                        }
                    }
                    //exists a bond in molecule 2, so that molecule 2 pair is connected?
                    for (int y = 0; y < target.getBondCount(); y++) {

                        productBond = target.getBond(y);

                        int t1 = target.getAtomNumber(productBond.getAtom(0));
                        int t2 = target.getAtomNumber(productBond.getAtom(1));

                        if (((compGraphNodes.get(a + 1) == t1)
                                && (compGraphNodes.get(b + 1) == t2))
                                || ((compGraphNodes.get(a + 1) == t2)
                                && (compGraphNodes.get(b + 1) == t1))) {
                            molecule2_pair_connected = true;
                            break;
                        }
                    }

                    if (reactantBond != null && productBond != null) {
                        addEdges(reactantBond, productBond, a, b, molecule1_pair_connected, molecule2_pair_connected);
                    }
                }
            }
        }
        cEdgesSize = cEdges.size();
        dEdgesSize = dEdges.size();
        return 0;
    }

    private void addEdges(IBond reactantBond, IBond productBond, int iIndex, int jIndex,
            boolean molecule1_pair_connected, boolean molecule2_pair_connected) {

        if (molecule1_pair_connected && molecule2_pair_connected) {
            if (isMatchFeasible(reactantBond, productBond, isMatchBond(), isMatchRings())) {
                cEdges.add((iIndex / 3) + 1);
                cEdges.add((jIndex / 3) + 1);
            } else {
                dEdges.add((iIndex / 3) + 1);
                dEdges.add((jIndex / 3) + 1);
            }
        } else if ((!molecule1_pair_connected && !molecule2_pair_connected)) {
            dEdges.add((iIndex / 3) + 1);
            dEdges.add((jIndex / 3) + 1);
        }
    }

    /**
     * compGraphNodesCZero is used to build up of the edges of the compatibility graph
     *
     * @return
     * @throws IOException
     */
    protected Integer compatibilityGraphNodesIfCEdgeIsZero() throws IOException {

        int count_nodes = 1;
        List<String> map = new ArrayList<String>();
        compGraphNodesCZero = new ArrayList<Integer>(); //Initialize the compGraphNodesCZero List
        LabelContainer labelContainer = LabelContainer.getInstance();
        compGraphNodes.clear();

        for (int i = 0; i < source.getAtomCount(); i++) {
            for (int j = 0; j < target.getAtomCount(); j++) {
                IAtom atom1 = source.getAtom(i);
                IAtom atom2 = target.getAtom(j);

                //You can also check object equal or charge, hydrogen count etc

                if ((atom1 instanceof IQueryAtom)
                        && ((IQueryAtom) atom1).matches(atom2)
                        && !map.contains(i + "_" + j)) {
                    compGraphNodesCZero.add(i);
                    compGraphNodesCZero.add(j);
                    compGraphNodesCZero.add(labelContainer.getLabelID(atom2.getSymbol())); //i.e C is label 1
                    compGraphNodesCZero.add(count_nodes);
                    compGraphNodes.add(i);
                    compGraphNodes.add(j);
                    compGraphNodes.add(count_nodes);
                    count_nodes += 1;
                    map.add(i + "_" + j);
                } else if (atom1.getSymbol().equalsIgnoreCase(atom2.getSymbol())
                        && !map.contains(i + "_" + j)) {
                    compGraphNodesCZero.add(i);
                    compGraphNodesCZero.add(j);
                    compGraphNodesCZero.add(labelContainer.getLabelID(atom1.getSymbol())); //i.e C is label 1
                    compGraphNodesCZero.add(count_nodes);
                    compGraphNodes.add(i);
                    compGraphNodes.add(j);
                    compGraphNodes.add(count_nodes);
                    count_nodes += 1;
                    map.add(i + "_" + j);
                }
            }
        }
        map.clear();
        return count_nodes;
    }

    /**
     * compatibilityGraphCEdgeZero is used to build up of the edges of the compatibility graph BIS
     *
     * @return
     * @throws IOException
     */
    protected int compatibilityGraphCEdgeZero() throws IOException {

        int compGraphNodesCZeroListSize = compGraphNodesCZero.size();
        cEdges = new ArrayList<Integer>(); //Initialize the cEdges List
        dEdges = new ArrayList<Integer>(); //Initialize the dEdges List

        for (int a = 0; a < compGraphNodesCZeroListSize; a += 4) {
            int index_a = compGraphNodesCZero.get(a);
            int index_aPlus1 = compGraphNodesCZero.get(a + 1);
            for (int b = a + 4; b < compGraphNodesCZeroListSize; b += 4) {
                int index_b = compGraphNodesCZero.get(b);
                int index_bPlus1 = compGraphNodesCZero.get(b + 1);

                // if element atomCont !=jIndex and atoms on the adjacent sides of the bonds are not equal
                if ((a != b) && (index_a != index_b)
                        && (index_aPlus1 != index_bPlus1)) {

                    IBond reactantBond = null;
                    IBond productBond = null;

                    reactantBond = source.getBond(source.getAtom(index_a), source.getAtom(index_b));
                    productBond = target.getBond(target.getAtom(index_aPlus1), target.getAtom(index_bPlus1));

                    if (reactantBond != null && productBond != null) {
                        addZeroEdges(reactantBond, productBond, a, b);
                    }

                }
            }
        }

        //Size of C and D edges of the compatibility graph
        cEdgesSize = cEdges.size();
        dEdgesSize = dEdges.size();
        return 0;
    }

    private void addZeroEdges(IBond reactantBond, IBond productBond, int indexI, int indexJ) {
        if (isMatchFeasible(reactantBond, productBond, isMatchBond(), isMatchRings())) {
            cEdges.add((indexI / 4) + 1);
            cEdges.add((indexJ / 4) + 1);
        }
        if (reactantBond == null && productBond == null) {
            dEdges.add((indexI / 4) + 1);
            dEdges.add((indexJ / 4) + 1);
        }
    }

    /**
     *
     * @param bondA1
     * @param bondA2
     * @param shouldMatchBonds
     * @param shouldMatchRings
     * @return
     */
    protected static boolean isMatchFeasible(
            IBond bondA1,
            IBond bondA2,
            boolean shouldMatchBonds,
            boolean shouldMatchRings) {

        if (bondA1 instanceof IQueryBond) {
            if (((IQueryBond) bondA1).matches(bondA2)) {
                IQueryAtom atom1 = (IQueryAtom) (bondA1.getAtom(0));
                IQueryAtom atom2 = (IQueryAtom) (bondA1.getAtom(1));
                // ok, bonds match
                if (atom1.matches(bondA2.getAtom(0)) && atom2.matches(bondA2.getAtom(1))
                        || atom1.matches(bondA2.getAtom(1)) && atom2.matches(bondA2.getAtom(0))) {
                    // ok, atoms match in either order
                    return true;
                }
                return false;
            }
            return false;
        } else {
            return DefaultMatcher.matches(bondA1, bondA2, shouldMatchBonds, shouldMatchRings) ? true : false;
        }
    }

    public List<Integer> getCEgdes() {
        return Collections.unmodifiableList(cEdges);
    }

    public List<Integer> getDEgdes() {
        return Collections.unmodifiableList(dEdges);
    }

    public List<Integer> getCompGraphNodes() {
        return Collections.unmodifiableList(compGraphNodes);
    }

    public int getCEdgesSize() {
        return cEdgesSize;
    }

    public int getDEdgesSize() {
        return dEdgesSize;
    }

    protected List<Integer> getCompGraphNodesCZero() {
        return Collections.unmodifiableList(compGraphNodesCZero);
    }

    protected void clearCEgdes() {
        cEdges.clear();
    }

    protected void clearDEgdes() {
        dEdges.clear();
    }

    protected void clearCompGraphNodes() {
        compGraphNodes.clear();
    }

    protected void clearCompGraphNodesCZero() {
        compGraphNodesCZero.clear();
    }

    protected void resetCEdgesSize() {
        cEdgesSize = 0;
    }

    protected void resetDEdgesSize() {
        dEdgesSize = 0;
    }

    public synchronized void clear() {
        cEdges = null;
        dEdges = null;
        compGraphNodes = null;
        compGraphNodesCZero = null;
    }

    /**
     * @return the shouldMatchBonds
     */
    public boolean isMatchBond() {
        return shouldMatchBonds;
    }

    /**
     * @return the shouldMatchRings
     */
    public boolean isMatchRings() {
        return shouldMatchRings;
    }
}
