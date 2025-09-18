/*
 * SPDX-License-Identifier: Apache-2.0
 * © 2025 BioInception PVT LTD.
 */
package org.openscience.smsd.algorithm.single;

import java.util.*;
import java.util.logging.Level;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.isomorphism.matchers.IQueryAtomContainer;
import org.openscience.cdk.tools.ILoggingTool;
import org.openscience.cdk.tools.LoggingToolFactory;
import org.openscience.smsd.AtomAtomMapping;
import org.openscience.smsd.interfaces.IResults;

/**
 * This is a handler class for single atom mapping
 * ({@link org.openscience.smsd.algorithm.single.SingleMapping}).
 *
 *  
 *
 * @author Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>
 */
public class SingleMappingHandler implements IResults {

    private final ILoggingTool Logger
            = LoggingToolFactory.createLoggingTool(SingleMappingHandler.class);
    private List<AtomAtomMapping> allAtomMCS = null;
    private final IAtomContainer source;
    private final IAtomContainer target;
    private final boolean shouldMatchRings;

    /**
     *
     * @param source
     * @param target
     * @param shouldMatchRings
     */
    public SingleMappingHandler(
            IAtomContainer source, 
            IAtomContainer target, 
            boolean shouldMatchRings) {
        allAtomMCS = new ArrayList<>();
        this.source = source;
        this.target = target;
        this.shouldMatchRings = shouldMatchRings;
        searchMCS();
    }

    /**
     *
     * @param source
     * @param target
     */
    public SingleMappingHandler(
            IQueryAtomContainer source, 
            IAtomContainer target) {
        allAtomMCS = new ArrayList<>();
        this.source = source;
        this.target = target;
        this.shouldMatchRings = true;
        searchMCS();
    }

    /**
     * Function is called by the main program and serves as a starting point for
     * the comparison procedure. {@inheritDoc}
     *
     */
    private synchronized void searchMCS() {
        SingleMapping singleMapping = new SingleMapping();
        List<Map<IAtom, IAtom>> mappings = null;
        try {
            if (target instanceof IQueryAtomContainer) {
                throw new CDKException("Target can't be IQueryAtomContainer");
            } else if (!(source instanceof IQueryAtomContainer)) {
                mappings = singleMapping.getOverLaps(source, target);
            } else {
                mappings = singleMapping.getOverLaps((IQueryAtomContainer) source, target);
            }
        } catch (CDKException ex) {
            Logger.error(Level.SEVERE, null, ex);
        }
        setAllAtomMapping(mappings);
        //setStereoScore();
    }

    private synchronized void setAllAtomMapping(List<Map<IAtom, IAtom>> mappings) {

        try {
            int counter = 0;
            for (Map<IAtom, IAtom> solution : mappings) {
                AtomAtomMapping atomMappings = new AtomAtomMapping(source, target);
                solution.entrySet().forEach((map) -> {
                    IAtom sourceAtom = map.getKey();
                    IAtom targetAtom = map.getValue();
                    atomMappings.put(sourceAtom, targetAtom);
                });
                allAtomMCS.add(counter, atomMappings);
                counter++;
            }
        } catch (Exception I) {
            I.getCause();
        }
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public synchronized List<AtomAtomMapping> getAllAtomMapping() {
        return Collections.unmodifiableList(allAtomMCS);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public synchronized AtomAtomMapping getFirstAtomMapping() {
        if (allAtomMCS.iterator().hasNext()) {
            return allAtomMCS.iterator().next();
        }
        return new AtomAtomMapping(source, target);
    }
}
