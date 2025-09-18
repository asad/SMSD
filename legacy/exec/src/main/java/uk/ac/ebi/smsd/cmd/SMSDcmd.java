/*
 * SPDX-License-Identifier: Apache-2.0
 * © 2025 BioInception PVT LTD.
 */
package uk.ac.ebi.smsd.cmd;

import java.io.IOException;
import java.util.*;

import org.apache.commons.cli.MissingOptionException;
import org.apache.commons.cli.ParseException;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.AtomContainer;
import org.openscience.cdk.AtomContainerSet;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.graph.ConnectivityChecker;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.tools.ILoggingTool;
import org.openscience.cdk.tools.LoggingToolFactory;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import org.openscience.smsd.AtomAtomMapping;
import org.openscience.smsd.BaseMapping;
import org.openscience.smsd.Isomorphism;
import org.openscience.smsd.Substructure;
import org.openscience.smsd.interfaces.Algorithm;
import org.openscience.smsd.mcss.JobType;
import org.openscience.smsd.mcss.MCSS;
import org.openscience.smsd.tools.AtomContainerComparator;
import org.openscience.smsd.tools.ExtAtomContainerManipulator;

/**
 *
 * @author Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>
 */
public class SMSDcmd {

    private final static ILoggingTool logger
            = LoggingToolFactory.createLoggingTool(InputHandler.class);

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

    /**
     *
     * @param argumentHandler
     * @param inputHandler
     */
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
            logger.error("IO Problem : " + ioe.getMessage());
//            ioe.printStackTrace();
        } catch (CDKException e) {
            logger.error("CDK Problem : " + e.getMessage());
//            e.printStackTrace();
        } catch (CloneNotSupportedException e) {
            logger.error(e.toString());
        } catch (MissingOptionException e) {
            logger.error("Missing argument : " + e.getMessage());
        }

    }

    /**
     *
     * @param inputHandler
     * @param outputHandler
     * @param argumentHandler
     * @throws IOException
     * @throws CDKException
     * @throws CloneNotSupportedException
     */
    public static void runNMCS(
            InputHandler inputHandler,
            OutputHandler outputHandler,
            ArgumentHandler argumentHandler) throws IOException, CDKException, CloneNotSupportedException {
        List<IAtomContainer> atomContainerSet = inputHandler.getAllTargets();
        String targetType = argumentHandler.getTargetType();
        if (atomContainerSet == null) {
            throw new IOException("Unknown input type " + targetType);
        }

        Comparator<IAtomContainer> comparator = new AtomContainerComparator();
        Collections.sort(atomContainerSet, comparator);

        boolean matchBonds = argumentHandler.isMatchBondType();
        boolean matchRings = argumentHandler.isMatchRingType();
        boolean matchAtomTypes = argumentHandler.isMatchAtomType();
        int filter = argumentHandler.getChemFilter();

        /*
         * Configure the targets
         */
        for (IAtomContainer target : atomContainerSet) {
            inputHandler.configure(target, targetType);
        }


        /*
         * Run N MULTIPLE on targets
         */
        MCSS mcss = new MCSS(atomContainerSet, JobType.MULTIPLE, 0,
                matchBonds, matchRings, matchAtomTypes);
        Collection<IAtomContainer> calculatedMCSS = mcss.getCalculateMCSS();
        IAtomContainerSet solutions = new AtomContainerSet();
        for (IAtomContainer mcsAtomContainer : calculatedMCSS) {
            if (mcsAtomContainer != null && mcsAtomContainer.getAtomCount() > 0) {
                boolean flag = ConnectivityChecker.isConnected(mcsAtomContainer);
                if (!flag) {
                    System.err.println("WARNING : Skipping file "
                            + mcsAtomContainer.getProperty(CDKConstants.TITLE) + " not connected ");
                    return;
                } else if (mcsAtomContainer.getProperty(CDKConstants.TITLE) != null) {
                    String mcsFilenName = mcsAtomContainer.getProperty(CDKConstants.TITLE).equals("untitled")
                            ? "mcs" : (String) mcsAtomContainer.getProperty(CDKConstants.TITLE);
                    mcsAtomContainer.setID(mcsFilenName);
                    argumentHandler.setQueryMolOutName(mcsAtomContainer.getID());
                } else if (mcsAtomContainer.getProperty(CDKConstants.TITLE) == null) {
                    String mcsFilenName = "Fragment";
                    mcsAtomContainer.setID(mcsFilenName);
                    argumentHandler.setQueryMolOutName(mcsAtomContainer.getID());
                }

                inputHandler.configure(mcsAtomContainer, targetType);
                solutions.addAtomContainer(mcsAtomContainer);
            }
        }

        if (argumentHandler.shouldOutputSubgraph()) {
            String outpath = argumentHandler.getOutputFilepath();
            String outtype = argumentHandler.getOutputFiletype();
            outputHandler.writeMol(outtype, solutions, outpath);
        }

        /*
         * For image generation RE-RUN the MULTIPLE with the common fragment
         */
        if (argumentHandler.isImage() && !solutions.isEmpty()) {
            int index = 1;
            for (IAtomContainer ac : solutions.atomContainers()) {
                if (ac != null && ac.getAtomCount() > 0) {
                    IAtomContainer mcsAtomContainer = ac.clone();
                    // now that we have the N-MULTIPLE, remap
                    List<Map<Integer, Integer>> mappings = new ArrayList<>();
                    List<IAtomContainer> secondRoundTargets = new ArrayList<>();
                    IChemObjectBuilder builder = DefaultChemObjectBuilder.getInstance();
                    for (IAtomContainer target : atomContainerSet) {
                        BaseMapping smsd = run(mcsAtomContainer, target, filter, matchBonds, matchRings, matchBonds);
                        mappings.add(getIndexMapping(smsd.getFirstAtomMapping()));
                        secondRoundTargets.add(
                                builder.newInstance(IAtomContainer.class, smsd.getFirstAtomMapping().getTarget()));
                    }

                    String name = inputHandler.getTargetName() + "_" + String.valueOf(index);
                    outputHandler.writeCircleImage(mcsAtomContainer, secondRoundTargets, name, mappings);
                }
                index++;
            }
        }
    }

    /**
     *
     * @param inputHandler
     * @param outputHandler
     * @param argumentHandler
     * @throws IOException
     * @throws CDKException
     * @throws CloneNotSupportedException
     */
    public static void runSingleQueryMultipleTarget(
            InputHandler inputHandler,
            OutputHandler outputHandler,
            ArgumentHandler argumentHandler) throws IOException, CDKException, CloneNotSupportedException {
        IAtomContainer query = inputHandler.getQuery();
        String name = (String) query.getProperty(CDKConstants.TITLE);
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
            query.setProperty(CDKConstants.TITLE, name);
            query.setID(name);
        }

        outputHandler.writeQueryMol(query);

        String out = ".out";
        outputHandler.startAppending(out);

        long startTime = System.currentTimeMillis();
        BaseMapping smsd;
        boolean matchBonds = argumentHandler.isMatchBondType();
        boolean matchRings = argumentHandler.isMatchRingType();
        boolean matchAtomTypes = argumentHandler.isMatchAtomType();

        int targetNumber = 0;
        List<IAtomContainer> allTargets = inputHandler.getAllTargets();
        String targetType = argumentHandler.getTargetType();
        if (allTargets == null) {
            throw new IOException("Unknown input type " + targetType);
        }
        for (IAtomContainer target : allTargets) {
            flag = ConnectivityChecker.isConnected(target);
            if (!flag) {
                logger.error("WARNING : Skipping target AtomContainer "
                        + target.getProperty(CDKConstants.TITLE) + " as it is not connected.");
                continue;
            }

            inputHandler.configure(target, targetType);

            if (argumentHandler.isSubstructureMode()) {
                smsd = runSubstructure(query, target, argumentHandler.getChemFilter(), matchBonds, matchRings, matchAtomTypes);
            } else {
                smsd = run(query, target, argumentHandler.getChemFilter(), matchBonds, matchRings, matchAtomTypes);
            }

            long endTime = System.currentTimeMillis();
            long executionTime = endTime - startTime;
            outputHandler.writeTargetMol(smsd.getTarget());

            String queryPath = argumentHandler.getQueryFilepath();
            String targetPath = argumentHandler.getTargetFilepath();

            IAtomContainer queryLocal = query.getBuilder().newInstance(IAtomContainer.class, smsd.getFirstAtomMapping().getQuery());
            IAtomContainer targetLocal = target.getBuilder().newInstance(IAtomContainer.class, smsd.getFirstAtomMapping().getTarget());
            Map<IAtom, IAtom> mcs = smsd.getFirstAtomMapping().getMappingsByAtoms();
            int nAtomsMatched = (mcs == null) ? 0 : mcs.size();
            double tanimotoSimilarity = smsd.getTanimotoSimilarity();
            //print out all mappings
            if (mcs != null && !mcs.isEmpty() && argumentHandler.isAllMapping()) {
                outputHandler.printHeader(queryPath, targetPath, nAtomsMatched);
                int counter = 0;
                for (AtomAtomMapping aam : smsd.getAllAtomMapping()) {
                    Map<Integer, Integer> mapping = aam.getMappingsByIndex();
                    if (argumentHandler.isImage() && !mapping.isEmpty()) {
                        double stereoScore = smsd.getStereoScore(counter);
                        String label = outputHandler.makeLabel(tanimotoSimilarity, stereoScore);
                        outputHandler.addImage(queryLocal, targetLocal, label, mapping);
                    }
                    outputHandler.printMapping((counter + 1), mapping);
                    counter += 1;
                }
            } //print out top one
            else if (mcs != null && !mcs.isEmpty() && !argumentHandler.isAllMapping()) {
                Map<Integer, Integer> mcsNumber = smsd.getFirstAtomMapping().getMappingsByIndex();
                double stereoScore = smsd.getStereoScore(0);
                outputHandler.printHeader(queryPath, targetPath, nAtomsMatched);
                String qrefName = inputHandler.getQRefName();
                String trefName = inputHandler.getTRefName();
                outputHandler.printTopMapping(
                        nAtomsMatched, mcs, mcsNumber, qrefName, trefName);
                if (argumentHandler.isImage() && !mcs.isEmpty()) {
                    String label = outputHandler.makeLabel(tanimotoSimilarity, stereoScore);
                    outputHandler.makeImage(queryLocal, targetLocal, label, mcsNumber);
                }
            }
            double tanimotoGraph = smsd.getTanimotoSimilarity();
//            double tanimotoAtom = smsd.getTanimotoAtomSimilarity();
//            double tanimotoBond = smsd.getTanimotoBondSimilarity();
            double euclidianGraph = smsd.getEuclideanDistance();
//            outputHandler.writeResults(queryLocal, targetLocal, tanimotoGraph, tanimotoAtom, tanimotoBond, euclidianGraph, nAtomsMatched, executionTime);
//            outputHandler.writeResults(queryLocal, targetLocal, tanimotoGraph, euclidianGraph, nAtomsMatched, executionTime);
            outputHandler.writeResults(queryLocal, targetLocal, tanimotoGraph, euclidianGraph, nAtomsMatched, executionTime, smsd.getFirstAtomMapping());

            if (mcs != null && !mcs.isEmpty() && argumentHandler.isImage()) {
                String qName = inputHandler.getQueryName();
                String tName = inputHandler.getTargetName() + "_" + targetNumber;
                outputHandler.writeImage(qName, tName);
            }
            targetNumber++;
        }
        outputHandler.closeFiles();
    }

    /**
     *
     * @param inputHandler
     * @param outputHandler
     * @param argumentHandler
     * @throws IOException
     * @throws CDKException
     * @throws CloneNotSupportedException
     */
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
        if (argumentHandler.isImage()) {
            boolean flag = ConnectivityChecker.isConnected(query);
            if (!flag) {
                logger.error("WARNING : Skipping file " + inputHandler.getQueryName() + " not connectted ");
                return;
            }
            flag = ConnectivityChecker.isConnected(target);
            if (!flag) {
                logger.error("WARNING : Skipping target AtomContainer "
                        + inputHandler.getTargetName() + " as it is not connected.");
                return;
            }
        }

        String fileNameQ = "Query";
        String fileNameT = "Target";

        if (target.getProperty(CDKConstants.TITLE) != null) {
            fileNameQ = target.getProperty(CDKConstants.TITLE) == null ? fileNameT
                    : (String) target.getProperty(CDKConstants.TITLE);
            target.setID(fileNameQ);
            argumentHandler.setTargetMolOutName(target.getID());
        }
        if (query.getProperty(CDKConstants.TITLE) != null) {
            fileNameT = query.getProperty(CDKConstants.TITLE) == null ? fileNameQ
                    : (String) query.getProperty(CDKConstants.TITLE);
            query.setID(fileNameT);
            argumentHandler.setQueryMolOutName(query.getID());
        }

        /*
         * remove hydrogens
         */
        if (removeHydrogens) {
            query = new AtomContainer(AtomContainerManipulator.removeHydrogens(query));
            query.setID(fileNameQ);
            target = new AtomContainer(AtomContainerManipulator.removeHydrogens(target));
            target.setID(fileNameT);
        }

        String out = ".out";
        if (!argumentHandler.isAppendMode()) {
            outputHandler.startAppending(out);
        } else {
            outputHandler.startNew(out);
        }

        ExtAtomContainerManipulator.aromatizeDayLight(query);
        ExtAtomContainerManipulator.aromatizeDayLight(target);

        if (argumentHandler.isApplyHAdding()) {
            AtomContainerManipulator.convertImplicitToExplicitHydrogens(query);
            AtomContainerManipulator.convertImplicitToExplicitHydrogens(target);
        }

        long startTime = System.currentTimeMillis();
        BaseMapping smsd;
        boolean matchBonds = argumentHandler.isMatchBondType();
        boolean matchRings = argumentHandler.isMatchRingType();
        boolean matchAtomTypes = argumentHandler.isMatchAtomType();

        if (argumentHandler.isSubstructureMode()) {
            smsd = runSubstructure(query, target, argumentHandler.getChemFilter(), matchBonds, matchRings, matchAtomTypes);
        } else {
            smsd = run(query, target, argumentHandler.getChemFilter(), matchBonds, matchRings, matchAtomTypes);
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
        Map<IAtom, IAtom> mcs = smsd.getFirstAtomMapping().getMappingsByAtoms();
        int nAtomsMatched = (mcs == null) ? 0 : mcs.size();
        double tanimotoSimilarity = smsd.getTanimotoSimilarity();

        //print out all mappings
        if (mcs != null && !mcs.isEmpty() && argumentHandler.isAllMapping()) {
            outputHandler.printHeader(queryPath, targetPath, nAtomsMatched);
            int counter = 0;
            for (AtomAtomMapping aam : smsd.getAllAtomMapping()) {
                Map<Integer, Integer> mapping = aam.getMappingsByIndex();
                if (argumentHandler.isImage() && !mapping.isEmpty()) {
                    double stereoScore = smsd.getStereoScore(counter);
                    String label = outputHandler.makeLabel(tanimotoSimilarity, stereoScore);
                    outputHandler.addImage(query, target, label, mapping);
                }
                outputHandler.printMapping((counter + 1), mapping);
                counter += 1;
            }
        } //print out top one
        else if (mcs != null && !mcs.isEmpty() && !argumentHandler.isAllMapping()) {
            Map<Integer, Integer> mcsNumber = smsd.getFirstAtomMapping().getMappingsByIndex();
            double stereoScore = smsd.getStereoScore(0);
            outputHandler.printHeader(queryPath, targetPath, nAtomsMatched);
            String qrefName = inputHandler.getQRefName();
            String trefName = inputHandler.getTRefName();
            outputHandler.printTopMapping(nAtomsMatched, mcs, mcsNumber, qrefName, trefName);
            if (argumentHandler.isImage() && !mcsNumber.isEmpty()) {
                String label = outputHandler.makeLabel(tanimotoSimilarity, stereoScore);
                outputHandler.makeImage(query, target, label, mcsNumber);
            }
        }

        double tanimotoGraph = smsd.getTanimotoSimilarity();
//        double tanimotoAtom = smsd.getTanimotoAtomSimilarity();
//        double tanimotoBond = smsd.getTanimotoBondSimilarity();
        double euclidianGraph = smsd.getEuclideanDistance();
//        outputHandler.writeResults(queryLocal, targetLocal, tanimotoGraph, tanimotoAtom, tanimotoBond, euclidianGraph, nAtomsMatched, executionTime);
//        outputHandler.writeResults(query, target, tanimotoGraph, euclidianGraph, nAtomsMatched, executionTime);

        outputHandler.writeResults(query, target, tanimotoGraph, euclidianGraph, nAtomsMatched, executionTime, smsd.getFirstAtomMapping());
        if (mcs != null && !mcs.isEmpty() && argumentHandler.isImage()) {
            String qName = inputHandler.getQueryName();
            String tName = inputHandler.getTargetName();
            outputHandler.writeImage(qName, tName);
        }

        if (argumentHandler.shouldOutputSubgraph()) {
            IAtomContainer subgraph = smsd.getFirstAtomMapping().getCommonFragment();
            String outpath = argumentHandler.getOutputFilepath();
            String outtype = argumentHandler.getOutputFiletype();
            outputHandler.writeMol(outtype, subgraph, outpath);
        }

        outputHandler.closeFiles();
    }

    private static BaseMapping run(
            IAtomContainer query,
            IAtomContainer target,
            int filter,
            boolean matchBonds,
            boolean matchRings,
            boolean matchAtomType) throws CDKException {
        // XXX - if clean and configure is 'true', is that not duplicate configuring?
        BaseMapping smsd = new Isomorphism(query, target, Algorithm.DEFAULT, matchBonds, matchRings, matchAtomType);

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
            boolean matchRings,
            boolean matchAtomTypes) throws CDKException {
        // XXX - if clean and configure is 'true', is that not duplicate configuring?
        BaseMapping smsd = new Substructure(query, target, matchBonds, matchRings, matchAtomTypes, true);

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
        return aam.isEmpty() ? new TreeMap<>() : aam.getMappingsByIndex();
    }
}
