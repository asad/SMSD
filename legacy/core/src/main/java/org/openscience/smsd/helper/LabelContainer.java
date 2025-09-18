/*
 * SPDX-License-Identifier: Apache-2.0
 * © 2025 BioInception PVT LTD.
 */
package org.openscience.smsd.helper;

import java.util.ArrayList;
import java.util.List;

/**
 * Class that handles atoms and assigns an integer label to them.
 *
 *
 *
 * @author Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>
 */
public class LabelContainer {

    private final List<String> labelMap;
    private int labelCounter = 0;
    private static LabelContainer instance = null;

    protected LabelContainer() {

        // System.err.println("List Initialized");
        labelMap = new ArrayList<>();
        labelMap.add(labelCounter++, "X");
        labelMap.add(labelCounter++, "R");
    }

    /**
     * Create ids from atom labels
     *
     * @return instance of this object
     */
    synchronized public static LabelContainer getInstance() {
        if (instance == null) {
            instance = new LabelContainer();
        }
        return instance;
    }

    /**
     * Add label if its not present
     *
     * @param label
     */
    synchronized public void addLabel(String label) {
        if (!labelMap.contains(label)) {
            labelMap.add(labelCounter++, label);
        }
    }

    /**
     * Returns label ID
     *
     * @param label
     * @return labelID
     */
    synchronized public Integer getLabelID(String label) {
        addLabel(label);
        return labelMap.indexOf(label);
    }

    /**
     * Returns Label of a given ID
     *
     * @param labelID
     * @return label
     */
    synchronized public String getLabel(Integer labelID) {
        return labelMap.get(labelID);
    }

    /**
     * Returns label count
     *
     * @return size of the labels
     */
    synchronized public int getSize() {
        return labelMap.size();
    }
}
