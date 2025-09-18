/*
 * SPDX-License-Identifier: Apache-2.0
 * © 2025 BioInception PVT LTD.
 */
package org.openscience.smsd.filters;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

/**
 * Class that cleans redundant mappings from the solution set.
 * <OL>
 *
 * <lI>1: Stereo match, bond type, ring etc,
 * <lI>2: Fragment size,
 * <lI>3: Bond breaking energy
 *
 * </OL>
 *
 *
 *
 * @author Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>
 */
public class PostFilter {

    /**
     *
     * Creates a new instance of Post Filter and removes redundant mapping(s).
     *
     * @param mappings
     * @return Filtered non-redundant mappings
     */
    public synchronized static List<Map<Integer, Integer>> filter(List<List<Integer>> mappings) {
        List<Map<Integer, Integer>> final_MAPPINGS = new ArrayList<>();
        if (mappings != null && !mappings.isEmpty()) {
            List<Map<Integer, Integer>> removeRedundantMapping = removeRedundantMapping(mappings);
            final_MAPPINGS.addAll(removeRedundantMapping);
        }
        return final_MAPPINGS;
    }

    private synchronized static boolean hasMap(Map<Integer, Integer> newMap, List<Map<Integer, Integer>> nonRedundantMapping) {
        for (Map<Integer, Integer> storedMap : nonRedundantMapping) {
            if (storedMap.equals(newMap)) {
                return true;
            }
        }
        return false;
    }

    /**
     *
     * @param mapping_org
     * @return
     */
    private synchronized static List<Map<Integer, Integer>> removeRedundantMapping(List<List<Integer>> mapping_org) {
        List<Map<Integer, Integer>> nonRedundantMapping = Collections.synchronizedList(new ArrayList<>());
        mapping_org.stream().map((M) -> getMappingMapFromList(M)).filter((newMap) -> (!hasMap(newMap, nonRedundantMapping))).forEachOrdered((newMap) -> {
            nonRedundantMapping.add(newMap);
        });
        return nonRedundantMapping;
    }

    private synchronized static Map<Integer, Integer> getMappingMapFromList(List<Integer> list) {
        Map<Integer, Integer> newMap = Collections.synchronizedSortedMap(new TreeMap<>());
        for (int index = 0; index < list.size(); index += 2) {
            newMap.put(list.get(index), list.get(index + 1));
        }
        return newMap;
    }
}
