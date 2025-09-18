/*
 * SPDX-License-Identifier: Apache-2.0
 * © 2025 BioInception PVT LTD.
 */
package org.openscience.smsd.algorithm.mcsplus;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.TreeMap;

/**
 * This class handles MCS between two identical molecules. Hence they generate
 * am MCS where all atoms are mapped.
 *
 *
 *
 * @author Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>
 */
public class ExactMapping {

    /**
     *
     * Extract atom mapping from the cliques and stores it in a List
     *
     * @param comp_graph_nodes
     * @param clique_List_org
     * @return
     */
    private static List<Integer> extractCliqueMapping(List<Integer> comp_graph_nodes, List<Integer> clique_List_org) {

        List<Integer> clique_mapping = Collections.synchronizedList(new ArrayList<>());
        List<Integer> clique_List = new ArrayList<>(clique_List_org);
        int clique_siz = clique_List.size();
        int vec_size = comp_graph_nodes.size();
//        System.out.println("VEC  SIZE " + vec_size);
        for (int a = 0; a < clique_siz; a++) {
            for (int b = 0; b < vec_size; b += 3) {
                if (Objects.equals(clique_List.get(a), comp_graph_nodes.get(b + 2))) {
                    clique_mapping.add(comp_graph_nodes.get(b));
                    clique_mapping.add(comp_graph_nodes.get(b + 1));
                }
            }
        }

        return clique_mapping;
    }

    //extract atom mapping from the clique List and print it on the screen
    /**
     *
     * @param _mappings
     * @param comp_graph_nodes
     * @param clique_List_org
     * @return
     */
    public static List<List<Integer>> extractMapping(List<List<Integer>> _mappings, List<Integer> comp_graph_nodes,
            List<Integer> clique_List_org) {
        try {
            List<Integer> clique_List = extractCliqueMapping(comp_graph_nodes, clique_List_org);
            _mappings.add(clique_List);
        } catch (Exception e) {
            System.err.println("Error in FinalMapping List: " + e.getCause());
            e.printStackTrace();
            System.exit(1);
        }
        return _mappings;
    }

    //extract atom mapping from the clique List and print it on the screen
    /**
     *
     * @param comp_graph_nodes
     * @param clique_List_org
     * @return
     */
    public static Map<Integer, Integer> extractMapping(List<Integer> comp_graph_nodes, List<Integer> clique_List_org) {
        Map<Integer, Integer> clique_mapping = Collections.synchronizedSortedMap(new TreeMap<Integer, Integer>());

        try {
            List<Integer> clique_List = new ArrayList<>(clique_List_org);

//        System.out.println("VEC  SIZE " + vec_size);
            for (int a = 0; a < clique_List.size(); a++) {
                for (int b = 0; b < comp_graph_nodes.size(); b += 3) {
                    if (Objects.equals(clique_List.get(a), comp_graph_nodes.get(b + 2))) {
                        clique_mapping.put(comp_graph_nodes.get(b), comp_graph_nodes.get(b + 1));
                    }
                }
            }

//            System.out.println("clique_mapping  SIZE " + clique_mapping.size());
        } catch (Exception e) {
            System.err.println("Error in FinalMapping List: " + e.getCause());
            e.printStackTrace();
            System.exit(1);
        }
        return clique_mapping;
    }
}
