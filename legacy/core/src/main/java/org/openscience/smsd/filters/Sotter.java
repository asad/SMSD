/*
 * SPDX-License-Identifier: Apache-2.0
 * © 2025 BioInception PVT LTD.
 */
package org.openscience.smsd.filters;

import java.util.Collections;
import java.util.Comparator;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

/**
 * @author Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>
 * @author maclean
 *
 */
public class Sotter {

    public synchronized static Map<Integer, Double> sortMapByValueInAscendingOrder(Map<Integer, Double> map) {
        List<Map.Entry<Integer, Double>> list = new LinkedList<>(map.entrySet());
        // Sort the list using an annonymous inner class implementing Comparator for the compare method
        Collections.sort(list, (Map.Entry<Integer, Double> entry, Map.Entry<Integer, Double> entry1) -> (entry.getValue().equals(entry1.getValue()) ? 0 : (entry.getValue() > entry1.getValue() ? 1 : -1)) // Return 0 for eAtom match, -1 for less than and +1 for more then (Aceending Order Sort)
        );
        // logger.info(list);
        Map<Integer, Double> result = new LinkedHashMap<>();
        list.forEach((entry) -> {
            result.put(entry.getKey(), entry.getValue());
        });
        return result;
    }

    public synchronized static Map<Integer, Double> sortMapByValueInDescendingOrder(Map<Integer, Double> map) {
        List<Map.Entry<Integer, Double>> list = new LinkedList<>(map.entrySet());
        // Sort the list using an annonymous inner class implementing Comparator for the compare method
        Collections.sort(list, (Map.Entry<Integer, Double> entry, Map.Entry<Integer, Double> entry1) -> (entry.getValue().equals(entry1.getValue()) ? 0
                : (entry.getValue() < entry1.getValue() ? 1 : -1)) // Return 0 for eAtom match, -1 for less than and +1 for more then (Decending Order Sort)
        );
        // logger.info(list);
        Map<Integer, Double> result = new LinkedHashMap<>();
        list.forEach((entry) -> {
            result.put(entry.getKey(), entry.getValue());
        });
        return result;
    }
}
