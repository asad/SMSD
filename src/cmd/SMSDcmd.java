/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package cmd;

import java.io.IOException;
import java.util.*;

import org.apache.commons.cli.MissingOptionException;
import org.apache.commons.cli.ParseException;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.AtomContainer;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.aromaticity.CDKHueckelAromaticityDetector;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.graph.ConnectivityChecker;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.io.iterator.IIteratingChemObjectReader;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import org.openscience.smsd.AtomAtomMapping;
import org.openscience.smsd.BaseMapping;
import org.openscience.smsd.Isomorphism;
import org.openscience.smsd.Substructure;
import org.openscience.smsd.interfaces.Algorithm;

/**
 *
 * @author sar
 */
public class SMSDcmd {

    /**
     * @param args the command line arguments
     */
    @SuppressWarnings("unchecked")
    public static void main(String[] args) {
        ArgumentHandler argumentHandler = new ArgumentHandler();
        try {
            argumentHandler.parseCommandLineOptions(args);
            InputHandler inputHandler = new InputHandler(argumentHandler);
            if (argumentHandler.isHelp()) {
                argumentHandler.printHelp();
                inputHandler.printDataTypeHelp();
            } else if (argumentHandler.isImageOptionHelp()) {
                OutputHandler outputHandler = new OutputHandler(argumentHandler);
                outputHandler.printImageOptionsHelp();
            } else {
                run(argumentHandler, inputHandler);
            }
        } catch (ParseException pe) {
            System.err.println("Problem with the arguments : " + pe.getMessage());
        }
    }

    public static void run(ArgumentHandler argumentHandler) {
        run(argumentHandler, new InputHandler(argumentHandler));
    }

    public static void run(ArgumentHandler argumentHandler, InputHandler inputHandler) {
        OutputHandler outputHandler = new OutputHandler(argumentHandler);
        try {
            InputHandler.MatchType matchType = inputHandler.validateInput();
            switch (matchType) {
                case SINGLE_QUERY_SINGLE_TARGET:
                    runSingleQuerySingleTarget(inputHandler, outputHandler, argumentHandler);
                    break;
                case SINGLE_QUERY_MULTIPLE_TARGET:
                    runSingleQueryMultipleTarget(inputHandler, outputHandler, argumentHandler);
                    break;
                case NMCS:
                    runNMCS(inputHandler, outputHandler, argumentHandler);
                    break;
                case UNKNOWN:
                default:
                    throw new IOException("Unknown types " + argumentHandler.getQueryType() + " " + argumentHandler.getTargetType());
            }

        } catch (IOException ioe) {
            System.err.println("IO Problem : " + ioe.getMessage());
//            ioe.printStackTrace();
        } catch (CDKException e) {
            System.err.println("CDK Problem : " + e.getMessage());
//            e.printStackTrace();
        } catch (CloneNotSupportedException e) {
            System.err.println(e.toString());
        } catch (MissingOptionException e) {
            System.err.println("Missing argument : " + e.getMessage());
        }

    }

    public static void runNMCS(
            InputHandler inputHandler,
            OutputHandler outputHandler,
            ArgumentHandler argumentHandler) throws IOException, CDKException, CloneNotSupportedException {
        IIteratingChemObjectReader reader = inputHandler.getAllTargets();
        String targetType = argumentHandler.getTargetType();
        if (reader == null) {
            throw new IOException("Unknown input type " + targetType);
        }
        List<IAtomContainer> atomContainerSet = new ArrayList<IAtomContainer>();
        while (reader.hasNext()) {
            IAtomContainer target = (IAtomContainer) reader.next();
            atomContainerSet.add(target);
        }

        Comparator<IAtomContainer> comparator = new AtomContainerComparator();
        Collections.sort(atomContainerSet, comparator);

        IAtomContainer mcsAtomContainer = null;
        boolean matchBonds = argumentHandler.isMatchBondType();
        boolean matchRings = argumentHandler.isMatchRingType();
        boolean removeHydrogens = argumentHandler.isApplyHRemoval();
        int filter = argumentHandler.getChemFilter();
        List<IAtomContainer> targets = new ArrayList<IAtomContainer>();

        for (IAtomContainer target : atomContainerSet) {
            boolean flag = ConnectivityChecker.isConnected(target);
            if (!flag) {
                System.err.println("WARNING : Skipping target AtomContainer "
                        + target.getProperty(CDKConstants.TITLE) + " as it is not connected.");
                continue;
            } else {
                if (target.getProperty(CDKConstants.TITLE) != null) {
                    target.setID((String) target.getProperty(CDKConstants.TITLE));
                    argumentHandler.setTargetMolOutName(target.getID() == null ? "Target" : target.getID());
                }
            }
            if (removeHydrogens) {
                target = new AtomContainer(AtomContainerManipulator.removeHydrogens(target));
            }

            if (mcsAtomContainer != null) {
                flag = ConnectivityChecker.isConnected(mcsAtomContainer);
                if (!flag) {
                    System.err.println("WARNING : Skipping file "
                            + mcsAtomContainer.getProperty(CDKConstants.TITLE) + " not connected ");
                    return;
                } else if (mcsAtomContainer.getProperty(CDKConstants.TITLE) != null) {
                    String mcsFilenName = mcsAtomContainer.getProperty(CDKConstants.TITLE).equals("untitled")
                            ? "mcs" : (String) mcsAtomContainer.getProperty(CDKConstants.TITLE);
                    mcsAtomContainer.setID(mcsFilenName);
                    argumentHandler.setQueryMolOutName(mcsAtomContainer.getID());
                }
                if (removeHydrogens) {
                    mcsAtomContainer = new AtomContainer(AtomContainerManipulator.removeHydrogens(mcsAtomContainer));
                }
            }

            inputHandler.configure(target, targetType);

            if (mcsAtomContainer == null) {
                mcsAtomContainer = target;
                targets.add(target);
            } else {
                BaseMapping smsd = run(mcsAtomContainer, target, filter, matchBonds, matchRings);
                target = target.getBuilder().newInstance(IAtomContainer.class, smsd.getFirstAtomMapping().getTarget());
                targets.add(target);
                Map<Integer, Integer> mapping = getIndexMapping(smsd.getFirstAtomMapping());
                IAtomContainer subgraph = getSubgraph(target, mapping);
                mcsAtomContainer = new AtomContainer(subgraph);
            }
        }
        inputHandler.configure(mcsAtomContainer, targetType);

        if (argumentHandler.shouldOutputSubgraph()) {
            String outpath = argumentHandler.getOutputFilepath();
            String outtype = argumentHandler.getOutputFiletype();
            outputHandler.writeMol(outtype, mcsAtomContainer, outpath);
        }
        if (mcsAtomContainer != null && argumentHandler.isImage()) {
            // now that we have the N-MCS, remap
            List<Map<Integer, Integer>> mappings = new ArrayList<Map<Integer, Integer>>();
            List<IAtomContainer> secondRoundTargets = new ArrayList<IAtomContainer>();
            IChemObjectBuilder builder = DefaultChemObjectBuilder.getInstance();
            for (IAtomContainer target : targets) {
                BaseMapping smsd = run(mcsAtomContainer, target, filter, matchBonds, matchRings);
                mappings.add(getIndexMapping(smsd.getFirstAtomMapping()));
                secondRoundTargets.add(
                        builder.newInstance(IAtomContainer.class, smsd.getFirstAtomMapping().getTarget()));
            }

            String name = inputHandler.getTargetName();
            outputHandler.writeCircleImage(mcsAtomContainer, secondRoundTargets, name, mappings);
        }
    }

    public static void runSingleQueryMultipleTarget(
            InputHandler inputHandler,
            OutputHandler outputHandler,
            ArgumentHandler argumentHandler) throws IOException, CDKException, CloneNotSupportedException {
        IAtomContainer query = inputHandler.getQuery();
        boolean removeHydrogens = argumentHandler.isApplyHRemoval();

        /*
         * check connectivity
         */
        boolean flag = ConnectivityChecker.isConnected(query);
        if (!flag) {
            System.err.println("WARNING : Skipping file " + inputHandler.getQueryName() + " not connected ");
            return;
        }
        if (removeHydrogens) {
            query = new AtomContainer(AtomContainerManipulator.removeHydrogens(query));
        }

        outputHandler.writeQueryMol(query);

        String out = ".out";
        outputHandler.startAppending(out);

        long startTime = System.currentTimeMillis();
        BaseMapping smsd = null;
        boolean matchBonds = argumentHandler.isMatchBondType();
        boolean matchRings = argumentHandler.isMatchRingType();

        int targetNumber = 0;
        IIteratingChemObjectReader reader = inputHandler.getAllTargets();
        String targetType = argumentHandler.getTargetType();
        if (reader == null) {
            throw new IOException("Unknown input type " + targetType);
        }
        while (reader.hasNext()) {
            IAtomContainer target = (IAtomContainer) reader.next();
            flag = ConnectivityChecker.isConnected(target);
            if (!flag) {
                System.err.println("WARNING : Skipping target AtomContainer "
                        + target.getProperty(CDKConstants.TITLE) + " as it is not connected.");
                continue;
            }

            /*
             * remove target hydrogens
             */
            if (removeHydrogens) {
                target = new AtomContainer(AtomContainerManipulator.removeHydrogens(target));
            }
            if (target.getProperty(CDKConstants.TITLE) != null) {
                target.setID((String) target.getProperty(CDKConstants.TITLE));
                argumentHandler.setTargetMolOutName(target.getID() == null ? "Target" : target.getID());
            }

            inputHandler.configure(target, targetType);

            if (argumentHandler.isSubstructureMode()) {
                smsd = runSubstructure(query, target, argumentHandler.getChemFilter(), matchBonds, matchRings);
            } else {
                smsd = run(query, target, argumentHandler.getChemFilter(), matchBonds, matchRings);
            }


            long endTime = System.currentTimeMillis();
            long executionTime = endTime - startTime;
            outputHandler.writeTargetMol(smsd.getTargetContainer());

            String queryPath = argumentHandler.getQueryFilepath();
            String targetPath = argumentHandler.getTargetFilepath();

            query = query.getBuilder().newInstance(IAtomContainer.class, smsd.getFirstAtomMapping().getQuery());
            target = target.getBuilder().newInstance(IAtomContainer.class, smsd.getFirstAtomMapping().getTarget());


            Map<IAtom, IAtom> mcs = smsd.getFirstAtomMapping().getMappings();
            int nAtomsMatched = (mcs == null) ? 0 : mcs.size();
            double tanimotoSimilarity = smsd.getTanimotoSimilarity();

            //print out all mappings
            if (mcs != null && argumentHandler.isAllMapping()) {
                outputHandler.printHeader(queryPath, targetPath, nAtomsMatched);
                for (AtomAtomMapping aam : smsd.getAllAtomMapping()) {
                    Map<Integer, Integer> mapping = getIndexMapping(aam);
                    int counter = 0;
                    if (argumentHandler.isImage()) {
                        double stereoScore = smsd.getStereoScore(counter);
                        String label = outputHandler.makeLabel(tanimotoSimilarity, stereoScore);
                        outputHandler.addImage(query, target, label, mapping);
                    }
                    outputHandler.printMapping((counter + 1), aam.getMappings());
                    counter += 1;
                }
            } //print out top one
            else if (mcs != null && !argumentHandler.isAllMapping()) {
                Map<Integer, Integer> mcsNumber = getIndexMapping(smsd.getFirstAtomMapping());
                double stereoScore = smsd.getStereoScore(0);
                outputHandler.printHeader(queryPath, targetPath, nAtomsMatched);
                String qrefName = inputHandler.getQRefName();
                String trefName = inputHandler.getTRefName();
                outputHandler.printTopMapping(
                        nAtomsMatched, mcs, mcsNumber, qrefName, trefName);
                if (argumentHandler.isImage()) {
                    String label = outputHandler.makeLabel(tanimotoSimilarity, stereoScore);
                    outputHandler.makeImage(query, target, label, mcsNumber);
                }
            }

            double tanimotoGraph = smsd.getTanimotoSimilarity();
//            double tanimotoAtom = smsd.getTanimotoAtomSimilarity();
//            double tanimotoBond = smsd.getTanimotoBondSimilarity();
            double euclidianGraph = smsd.getEuclideanDistance();
//            outputHandler.writeResults(query, target, tanimotoGraph, tanimotoAtom, tanimotoBond, euclidianGraph, nAtomsMatched, executionTime);
            outputHandler.writeResults(query, target, tanimotoGraph, euclidianGraph, nAtomsMatched, executionTime);

            if (mcs != null && argumentHandler.isImage()) {
                String qName = inputHandler.getQueryName();
                String tName = inputHandler.getTargetName() + "_" + targetNumber;
                outputHandler.writeImage(qName, tName);
            }
            targetNumber++;
        }
        outputHandler.closeFiles();
    }

    public static void runSingleQuerySingleTarget(
            InputHandler inputHandler,
            OutputHandler outputHandler,
            ArgumentHandler argumentHandler)
            throws IOException, CDKException, CloneNotSupportedException {
        IAtomContainer query = inputHandler.getQuery();
        IAtomContainer target = inputHandler.getTarget();

        boolean removeHydrogens = argumentHandler.isApplyHRemoval();

        /*
         * check connectivity
         */
        boolean flag = ConnectivityChecker.isConnected(query);
        if (!flag) {
            System.err.println("WARNING : Skipping file " + inputHandler.getQueryName() + " not connectted ");
            return;
        }
        flag = ConnectivityChecker.isConnected(target);
        if (!flag) {
            System.err.println("WARNING : Skipping target AtomContainer "
                    + inputHandler.getTargetName() + " as it is not connected.");
            return;
        }
        /*
         * remove hydrogens
         */
        if (removeHydrogens) {
            query = new AtomContainer(AtomContainerManipulator.removeHydrogens(query));
            target = new AtomContainer(AtomContainerManipulator.removeHydrogens(target));
        }

        if (target.getProperty(CDKConstants.TITLE) != null) {
            String fileName = target.getProperty(CDKConstants.TITLE) == null ? "Target"
                    : (String) target.getProperty(CDKConstants.TITLE);
            target.setID(fileName);
            argumentHandler.setTargetMolOutName(target.getID());
        }
        if (query.getProperty(CDKConstants.TITLE) != null) {
            String fileName = query.getProperty(CDKConstants.TITLE) == null ? "Query"
                    : (String) query.getProperty(CDKConstants.TITLE);
            query.setID(fileName);
            argumentHandler.setQueryMolOutName(query.getID());
        }

        String out = ".out";
        if (!argumentHandler.isAppendMode()) {
            outputHandler.startAppending(out);
        } else {
            outputHandler.startNew(out);
        }

        // XXX - this is also done in the InputHandler!
        CDKHueckelAromaticityDetector.detectAromaticity(query);
        CDKHueckelAromaticityDetector.detectAromaticity(target);

        // XXX - this is also done in the InputHandler!
        AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(query);
        AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(target);

        if (argumentHandler.isApplyHAdding()) {
            AtomContainerManipulator.convertImplicitToExplicitHydrogens(query);
            AtomContainerManipulator.convertImplicitToExplicitHydrogens(target);
        }

        long startTime = System.currentTimeMillis();
        BaseMapping smsd;
        boolean matchBonds = argumentHandler.isMatchBondType();
        boolean matchRings = argumentHandler.isMatchRingType();


        if (argumentHandler.isSubstructureMode()) {
            smsd = runSubstructure(query, target, argumentHandler.getChemFilter(), matchBonds, matchRings);
        } else {
            smsd = run(query, target, argumentHandler.getChemFilter(), matchBonds, matchRings);
        }

        query = query.getBuilder().newInstance(IAtomContainer.class, smsd.getFirstAtomMapping().getQuery());
        target = target.getBuilder().newInstance(IAtomContainer.class, smsd.getFirstAtomMapping().getTarget());

        long endTime = System.currentTimeMillis();
        long executionTime = endTime - startTime;

        // write out the input AtomContainers to files
        outputHandler.writeQueryMol(smsd.getFirstAtomMapping().getQuery());
        outputHandler.writeTargetMol(smsd.getFirstAtomMapping().getTarget());

        String queryPath = argumentHandler.getQueryFilepath();
        String targetPath = argumentHandler.getTargetFilepath();

        Map<IAtom, IAtom> mcs = smsd.getFirstAtomMapping().getMappings();
        int nAtomsMatched = (mcs == null) ? 0 : mcs.size();
        double tanimotoSimilarity = smsd.getTanimotoSimilarity();

        //print out all mappings
        if (mcs != null && argumentHandler.isAllMapping()) {
            outputHandler.printHeader(queryPath, targetPath, nAtomsMatched);
            for (Iterator<AtomAtomMapping> it = smsd.getAllAtomMapping().iterator(); it.hasNext();) {
                AtomAtomMapping aam = it.next();
                Map<Integer, Integer> mapping = getIndexMapping(aam);
                int counter = 0;
                if (argumentHandler.isImage()) {
                    double stereoScore = smsd.getStereoScore(counter);
                    String label = outputHandler.makeLabel(tanimotoSimilarity, stereoScore);
                    outputHandler.addImage(query, target, label, mapping);
                }
                outputHandler.printMapping((counter + 1), aam.getMappings());
                counter += 1;
            }
        } //print out top one
        else if (mcs != null && !argumentHandler.isAllMapping()) {
            Map<Integer, Integer> mcsNumber = getIndexMapping(smsd.getFirstAtomMapping());
            double stereoScore = smsd.getStereoScore(0);
            outputHandler.printHeader(queryPath, targetPath, nAtomsMatched);
            String qrefName = inputHandler.getQRefName();
            String trefName = inputHandler.getTRefName();
            outputHandler.printTopMapping(
                    nAtomsMatched, mcs, mcsNumber, qrefName, trefName);
            if (argumentHandler.isImage()) {
                String label = outputHandler.makeLabel(tanimotoSimilarity, stereoScore);
                outputHandler.makeImage(query, target, label, mcsNumber);
            }
        }

        double tanimotoGraph = smsd.getTanimotoSimilarity();
//        double tanimotoAtom = smsd.getTanimotoAtomSimilarity();
//        double tanimotoBond = smsd.getTanimotoBondSimilarity();
        double euclidianGraph = smsd.getEuclideanDistance();
//        outputHandler.writeResults(query, target, tanimotoGraph, tanimotoAtom, tanimotoBond, euclidianGraph, nAtomsMatched, executionTime);
        outputHandler.writeResults(query, target, tanimotoGraph, euclidianGraph, nAtomsMatched, executionTime);

        if (mcs != null && argumentHandler.isImage()) {
            String qName = inputHandler.getQueryName();
            String tName = inputHandler.getTargetName();
            outputHandler.writeImage(qName, tName);
        }

        if (argumentHandler.shouldOutputSubgraph()) {
            Map<Integer, Integer> mapping = getIndexMapping(smsd.getFirstAtomMapping());
            IAtomContainer subgraph = getSubgraph(target, mapping);
            String outpath = argumentHandler.getOutputFilepath();
            String outtype = argumentHandler.getOutputFiletype();
            outputHandler.writeMol(outtype, subgraph, outpath);
        }
        outputHandler.closeFiles();
    }

    private static IAtomContainer getSubgraph(
            IAtomContainer container, Map<Integer, Integer> mapping) throws CloneNotSupportedException {
        Collection<Integer> values = mapping.values();
        List<IAtom> subgraphAtoms = new ArrayList<IAtom>();
        IAtomContainer subgraph = (IAtomContainer) container.clone();
        for (Integer index : values) {
            subgraphAtoms.add(subgraph.getAtom(index));
        }
        List<IAtom> atoms = new ArrayList<IAtom>();
        for (IAtom atom : subgraph.atoms()) {
            atoms.add(atom);
        }
        for (IAtom atom : atoms) {
            if (!subgraphAtoms.contains(atom)) {
                subgraph.removeAtomAndConnectedElectronContainers(atom);
            }
        }
        return subgraph;
    }

    private static BaseMapping run(
            IAtomContainer query,
            IAtomContainer target,
            int filter,
            boolean matchBonds,
            boolean matchRings) throws CDKException {
        // XXX - if clean and configure is 'true', is that not duplicate configuring?
        BaseMapping smsd = new Isomorphism(query, target, Algorithm.VFLibMCS, matchBonds, matchRings);

        if (filter == 0) {
            smsd.setChemFilters(false, false, false);
        }
        if (filter == 1) {
            smsd.setChemFilters(true, false, false);
        }
        if (filter == 2) {
            smsd.setChemFilters(true, true, false);
        }
        if (filter == 3) {
            smsd.setChemFilters(true, true, true);
        }
        return smsd;
    }

    private static BaseMapping runSubstructure(
            IAtomContainer query,
            IAtomContainer target,
            int filter,
            boolean matchBonds,
            boolean matchRings) throws CDKException {
        // XXX - if clean and configure is 'true', is that not duplicate configuring?
        BaseMapping smsd = new Substructure(query, target, matchBonds, matchRings, false);

        if (smsd.isSubgraph()) {
            if (filter == 0) {
                smsd.setChemFilters(false, false, false);
            }
            if (filter == 1) {
                smsd.setChemFilters(true, false, false);
            }
            if (filter == 2) {
                smsd.setChemFilters(true, true, false);
            }
            if (filter == 3) {
                smsd.setChemFilters(true, true, true);
            }
        }

        return smsd;
    }

    private static Map<Integer, Integer> getIndexMapping(AtomAtomMapping aam) {
        Map<IAtom, IAtom> mappings = aam.getMappings();
        Map<Integer, Integer> mapping = new TreeMap<Integer, Integer>();
        for (IAtom keys : mappings.keySet()) {
            mapping.put(aam.getQueryIndex(keys), aam.getTargetIndex(mappings.get(keys)));
        }
        return mapping;
    }
}
