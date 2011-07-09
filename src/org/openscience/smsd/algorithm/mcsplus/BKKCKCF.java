/**
 *
 * Copyright (C) 2006-2011  Syed Asad Rahman <asad@ebi.ac.uk>
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
 * You should have received index copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA.
 */
package org.openscience.smsd.algorithm.mcsplus;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Stack;
import org.openscience.cdk.annotations.TestClass;

/**
 * This class implements Bron-Kerbosch clique detection algorithm as it is
 * described in [F. Cazals, C. Karande: An Algorithm for reporting maximal c-cliques;
 * processedVertex.Comp. Sc. (2005); vol 349; pp.
 * 484-490]
 *
 *
 * BronKerboschCazalsKarandeKochCliqueFinder.java
 *
 * @cdk.githash
 * @cdk.module smsd
 * @author Syed Asad Rahman <asad@ebi.ac.uk>
 */
@TestClass("org.openscience.cdk.smsd.BKKCKCFTest")
public class BKKCKCF {

    private List<List<Integer>> maxCliquesSet = null;
    /***********************************************************************/
    private List<Integer> cEdges = null;
    private List<Integer> dEdges = null;
    private int bestCliqueSize = 0;
    private List<Integer> compGraphNodes = null;
    private double dEdgeIterationSize = 0;
    private double cEdgeIterationSize = 0;
    /*
     *V: stored all the vertices for the Graph G
     * V[G]
     *nodes of vector comp_graph_nodes are stored in V
     */
    private Stack<Integer> V;

    /**
     * Creates index new instance of Bron Kerbosch Cazals Karande Koch Clique Finder
     * This class implements Bron-Kerbosch clique detection algorithm as it is
     * described in [F. Cazals, C. Karande: An Algorithm for reporting maximal c-cliques;
     * processedVertex.Comp. Sc. (2005); vol 349; pp.
     * 484-490]
     * @param compGraphNodes
     * @param edgesC C-Edges set of allowed edges
     * @param edgesD D-Edges set of prohibited edges
     */
    public BKKCKCF(List<Integer> compGraphNodes, List<Integer> edgesC, List<Integer> edgesD) {
        this.compGraphNodes = compGraphNodes;
        this.cEdges = edgesC;
        this.dEdges = edgesD;
        bestCliqueSize = 0;
        //Orignal assignment as per paper
        dEdgeIterationSize = dEdges.size() / 2;

        //Orignal assignment as per paper
        cEdgeIterationSize = cEdges.size() / 2;

        //Initialization maxCliquesSet

        maxCliquesSet = Collections.synchronizedList(new ArrayList<List<Integer>>());

//        System.out.println("C-edges " + edgesC.size());
//        System.out.println("d-edges " + edgesD.size());

        init();

    }

    /*
     * Call the wrapper for ENUMERATE_CLIQUES
     *
     */
    private synchronized void init() {


        /********************************************************************/
        /*
         * vertex: stored all the vertices for the Graph G
         * vertex[G]
         * nodes of vector compGraphNodes are stored in vertex
         */
        V = new Stack<Integer>(); //Initialization of ArrayList vertex

        int vertexCount = compGraphNodes.size() / 3;

        //System.out.println("ArrayList vertex is initialized");
        for (int a = 0; a < vertexCount; a++) {
            V.add(compGraphNodes.get(a * 3 + 2));
            //System.out.print("vertex[" + index + "]: " + compGraphNodes.get(index * 3 + 2) + " ");
        }
        //System.out.println();

        V.add(0);
//        System.out.println("Vector V :" + V);

        /*
         *processedVertex: is index set of vertices which have already been used
         */
        List<Integer> processedVertex = new ArrayList<Integer>();
        /*
         * Let processedVertex be the set of Nodes already been used in the initialization
         *
         */
        initIterator(processedVertex);
        processedVertex.clear();
        //System.out.println("maxCliquesSet: " + maxCliquesSet);

    }

    private synchronized int enumerateCliques(List<Integer> C, Stack<Integer> P,
            List<Integer> D, List<Integer> S, List<Integer> excludedCVertex) {
        List<Integer> potentialVertex = new ArrayList<Integer>();//Defined as P' in the paper

        for (Integer I : P) {
            potentialVertex.add(I);
        }

//        System.out.print("P " + P);
//        System.out.print("S " + S);

        if ((P.size() == 1) && (S.isEmpty())) {

//            System.out.print("Clique found!!!!!!!!!!!!!!!!");
//            System.out.print("\n");
//            int c_num = C.size();
//            System.out.print("C Vector Size: ");
//            System.out.print(c_num);
//            System.out.print("\n");
//            for (int cl = 0; cl < c_num; cl++) {
//                System.out.print(C.get(cl));
//                System.out.print(" ");
//            }
//            System.out.print("\n");

            //store best solutions in stack maxCliquesSet
            int clique_size = C.size();

            if (clique_size >= bestCliqueSize) {
                if (clique_size > bestCliqueSize) {

                    maxCliquesSet.clear();
                    bestCliqueSize = clique_size;

                }
                if (clique_size == bestCliqueSize) {
//                    System.out.println("C-Clique " + C);
                    Collections.sort(C);
                    if (!maxCliquesSet.contains(C)) {
                        maxCliquesSet.add(C);
                    }
                }
            }
            return 0;
        }
        findCliques(
                potentialVertex,
                C,
                P,
                D,
                S,
                excludedCVertex);
        return 0;
    }

    private synchronized void findCliques(List<Integer> potentialVertex, List<Integer> C,
            Stack<Integer> P, List<Integer> D, List<Integer> S,
            List<Integer> excludedCVertex) {
        int index = 0;
        List<Integer> neighbourVertex = new ArrayList<Integer>(); ////Initialization ArrayList N

        while (potentialVertex.get(index) != 0) {
            int potentialVertexIndex = potentialVertex.get(index);

            P.removeElement(potentialVertexIndex);

            List<Integer> R_copy = new ArrayList<Integer>(C);
            Stack<Integer> P_copy = new Stack<Integer>();
            List<Integer> Q_copy = new ArrayList<Integer>(D);
            List<Integer> X_copy = new ArrayList<Integer>(S);
            List<Integer> Y_copy = new ArrayList<Integer>(excludedCVertex);

            neighbourVertex.clear();


            for (Integer obj : P) {
                P_copy.add(obj);
            }

            P_copy.pop();
            //find the neighbors of the central node from P
            //System.out.println("potentialVertex.elementAt(index): " + potentialVertex.elementAt(index));

            neighbourVertex = findNeighbors(potentialVertexIndex);
            groupNeighbors(index,
                    P_copy,
                    Q_copy,
                    X_copy,
                    Y_copy,
                    neighbourVertex,
                    D,
                    potentialVertex,
                    S,
                    excludedCVertex);
            Stack<Integer> P_copy_N_intersec = new Stack<Integer>();
            List<Integer> Q_copy_N_intersec = new ArrayList<Integer>();
            List<Integer> X_copy_N_intersec = new ArrayList<Integer>();
            List<Integer> Y_copy_N_intersec = new ArrayList<Integer>();

            copyVertex(neighbourVertex,
                    P_copy_N_intersec,
                    P_copy,
                    Q_copy_N_intersec,
                    Q_copy,
                    X_copy_N_intersec,
                    X_copy,
                    Y_copy_N_intersec,
                    Y_copy);

            P_copy_N_intersec.push(0);
            R_copy.add(potentialVertexIndex);
            enumerateCliques(R_copy, P_copy_N_intersec, Q_copy_N_intersec, X_copy_N_intersec, Y_copy_N_intersec);
            S.add(potentialVertexIndex);
            index++;
        }
    }

    private synchronized void groupNeighbors(int index,
            Stack<Integer> P_copy,
            List<Integer> Q_copy,
            List<Integer> X_copy,
            List<Integer> Y_copy,
            List<Integer> neighbourVertex,
            List<Integer> potentialDVertex,
            List<Integer> potentialVertex,
            List<Integer> excludedVertex,
            List<Integer> excludedCVertex) {

        int N_size = neighbourVertex.size();

//        System.out.println("Neighbors: " + neighbourVertex);

        for (int b = 0; b < N_size; b += 2) {
            // N[index] is node v
            //Grouping of the neighbors:


            Integer Nelement_at_b = neighbourVertex.get(b);

            if (neighbourVertex.get(b + 1) == 1) {
                //u and v are adjacent via index C-edge

                if (potentialDVertex.contains(Nelement_at_b)) {

                    P_copy.push(Nelement_at_b);
                    //delete N[index] bzw. D[c] from set Q_copy, remove C-edges
                    Q_copy.remove(Nelement_at_b);

                }
                if (excludedCVertex.contains(Nelement_at_b)) {
                    if (excludedVertex.contains(Nelement_at_b)) {
                        X_copy.add(Nelement_at_b);
                    }
                    Y_copy.remove(Nelement_at_b);
                }
            }

            //find respective neighbor position in potentialVertex, which is needed for the deletion from potentialVertex

            if (potentialVertex.indexOf(Nelement_at_b) <= index && potentialVertex.indexOf(Nelement_at_b) > -1) {
                --index;
            }
            potentialVertex.remove(Nelement_at_b);
        }
    }

    private synchronized List<Integer> findNeighbors(int central_node) {

        List<Integer> neighborVertex = new ArrayList<Integer>();

        for (int a = 0; a < cEdgeIterationSize; a++) {
            if (cEdges.get(a * 2 + 0) == central_node) {
                //          System.out.println( cEdges.get(index*2+0) + " " + cEdges.get(index*2+1));
                neighborVertex.add(cEdges.get(a * 2 + 1));
                neighborVertex.add(1); // 1 means: is connected via C-edge
            } else if (cEdges.get(a * 2 + 1) == central_node) {
                //           System.out.println(cEdges.get(index*2+0) + " " + cEdges.get(index*2+1));
                neighborVertex.add(cEdges.get(a * 2 + 0));
                neighborVertex.add(1); // 1 means: is connected via C-edge
            }

        }
        for (int a = 0; a < dEdgeIterationSize; a++) {
            if (dEdges.get(a * 2 + 0) == central_node) {
                //       System.out.println( dEdges.get(index*2+0) + " " + dEdges.get(index*2+1));
                neighborVertex.add(dEdges.get(a * 2 + 1));
                neighborVertex.add(2); // 2 means: is connected via D-edge
            } else if (dEdges.get(a * 2 + 1) == central_node) {
                //        System.out.println(dEdges.get(index*2+0) + " " + dEdges.get(index*2+1));
                neighborVertex.add(dEdges.get(a * 2 + 0));
                neighborVertex.add(2); // 2 means: is connected via D-edge
            }
        }
        return neighborVertex;
    }

    protected synchronized int getBestCliqueSize() {
        return bestCliqueSize;
    }

    /**
     * Maximum clique set
     * @return
     */
    public synchronized Stack<List<Integer>> getMaxCliqueSet() {
        Stack<List<Integer>> solution = new Stack<List<Integer>>();
        for (List<Integer> list : maxCliquesSet) {
            Collections.sort(list);
            if (!solution.contains(list)) {
                solution.add(list);
            }
        }
        return solution;
    }

    private synchronized void copyVertex(List<Integer> neighbourVertex, Stack<Integer> P_copy_N_intersec, Stack<Integer> P_copy,
            List<Integer> Q_copy_N_intersec, List<Integer> Q_copy, List<Integer> X_copy_N_intersec,
            List<Integer> X_copy, List<Integer> Y_copy_N_intersec, List<Integer> Y_copy) {
        int nElement = -1;
        int N_size = neighbourVertex.size();

        for (int sec = 0; sec < N_size; sec += 2) {

            nElement = neighbourVertex.get(sec);

            if (P_copy.contains(nElement)) {
                P_copy_N_intersec.push(nElement);
            }
            if (Q_copy.contains(nElement)) {
                Q_copy_N_intersec.add(nElement);
            }
            if (X_copy.contains(nElement)) {
                X_copy_N_intersec.add(nElement);
            }
            if (Y_copy.contains(nElement)) {
                Y_copy_N_intersec.add(nElement);
            }
        }
    }

    private synchronized void initIterator(List<Integer> T) {
        /*
         * C: set of vertices belonging to the current clique
         */
        List<Integer> C = new ArrayList<Integer>();
        /*
         *P: is index set of vertices which <index>can</index> be added
         *to C, because they are
         * neighbours of V u via <i>c-edges</i>
         */
        Stack<Integer> P = new Stack<Integer>();
        /*
         *D: is index set of vertices which <index>cannot</index> be added to
         *C, because they are
         * neighbours of V u via <i>d-edges</i>
         */

        List<Integer> D = new ArrayList<Integer>();
        /*
         *S: set of vertices which are not allowed to be added
         * to C
         */
        List<Integer> S = new ArrayList<Integer>();


        /*
         *excludedCVertex: set of vertices which are not allowed to be added
         * to C
         */

        List<Integer> excludedCVertex = new ArrayList<Integer>();

        /*
         * N[u]: set of neighbours of V u in Graph G
         *
         */

        List<Integer> N = new ArrayList<Integer>();

        int index = 0;
        while (V.get(index) != 0) {
            int central_node = V.get(index);
            P.clear();
            D.clear();
            S.clear();
            C.clear();

            //find the neighbors of the central node from V
            N = findNeighbors(central_node);

            for (int c = 0; c < N.size(); c = c + 2) {
                /*
                 * u and v are adjacent via index C-edge
                 */
                Integer neighbourVertexOfC = N.get(c);

                //find respective neighbor position in P, which is needed for the deletion from V
                //delete neighbor from set V

                if (N.get(c + 1) == 1) {
                    // u and v are adjacent via index D-edge
//                    System.out.println("u and v are adjacent via index C-edge: " + neighbourVertexOfC);
                    if (T.contains(neighbourVertexOfC)) {
                        S.add(neighbourVertexOfC);
                    } else {
                        P.push(neighbourVertexOfC);
                    }
                } else if (N.get(c + 1) == 2) {
                    // u and v are adjacent via index D-edge
//                    System.out.println("u and v are adjacent via index D-edge: " + neighbourVertexOfC);

                    if (T.contains(neighbourVertexOfC)) {
                        excludedCVertex.add(neighbourVertexOfC);
                    } else {
                        D.add(neighbourVertexOfC);
                    }
                }

                if (V.indexOf(neighbourVertexOfC) <= index && V.indexOf(neighbourVertexOfC) > -1) {
                    --index;
                }
                V.remove(neighbourVertexOfC);
                //System.out.println("Elements Removed from V:" + neighbourVertexOfC);
            }

            P.add(0);
            C.add(central_node);

//            System.out.println(" Calling Enumerate_Cliques, Vector Size in CliquesGenerator: ");
//
//            System.out.println("C: " + C.size() + ", P: " + P.size() + ", D: " + D.size() + ", S: " + S.size());

            enumerateCliques(C, P, D, S, excludedCVertex);
            T.add(V.get(index));
            index++;
        }
    }
}
