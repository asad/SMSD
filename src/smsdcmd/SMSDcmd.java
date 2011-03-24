/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package smsdcmd;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;

import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.openscience.cdk.Molecule;
import org.openscience.cdk.aromaticity.CDKHueckelAromaticityDetector;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.interfaces.IMolecule;
import org.openscience.cdk.io.iterator.IIteratingChemObjectReader;
import org.openscience.cdk.nonotify.NoNotificationChemObjectBuilder;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import org.openscience.smsd.IAtomAtomMapping;
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
            } else {
                run(argumentHandler, inputHandler);
            }
        } catch (ParseException pe) {
            System.err.println("Problem with the arguments : " + pe.getMessage());
        }
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
        IMolecule mcsMolecule = null;
        boolean matchBonds = argumentHandler.isMatchBondType();
        boolean removeHydrogens = argumentHandler.isApplyHRemoval();
        int filter = argumentHandler.getChemFilter();
//        Isomorphism smsd = new Isomorphism(Algorithm.DEFAULT, matchBonds);
        List<IAtomContainer> targets = new ArrayList<IAtomContainer>();
        while (reader.hasNext()) {
            IMolecule target = (IMolecule) reader.next();
            inputHandler.configure(target, targetType);
            if (mcsMolecule == null) {
                mcsMolecule = target;
                targets.add(target);
            } else {
                Isomorphism smsd = new Isomorphism(Algorithm.DEFAULT, matchBonds);
                run(smsd, mcsMolecule, target, filter);
                target = target.getBuilder().newInstance(IMolecule.class, smsd.getProductMolecule());
                targets.add(target);
                Map<Integer, Integer> mapping = smsd.getFirstMapping();
                IAtomContainer subgraph = getSubgraph(target, mapping);
                mcsMolecule = new Molecule(subgraph);
            }
        }
        inputHandler.configure(mcsMolecule, targetType);

        if (argumentHandler.shouldOutputSubgraph()) {
            String outpath = argumentHandler.getOutputFilepath();
            String outtype = argumentHandler.getOutputFiletype();
            outputHandler.writeMol(outtype, mcsMolecule, outpath);
        }
        if (mcsMolecule != null && argumentHandler.isImage()) {
            // now that we have the N-MCS, remap
            List<Map<Integer, Integer>> mappings = new ArrayList<Map<Integer, Integer>>();
            List<IAtomContainer> secondRoundTargets = new ArrayList<IAtomContainer>();
            IChemObjectBuilder builder = NoNotificationChemObjectBuilder.getInstance();
            for (IAtomContainer target : targets) {
                Isomorphism smsd = new Isomorphism(Algorithm.DEFAULT, matchBonds);
                run(smsd, mcsMolecule, (IMolecule) target, filter);
                mappings.add(smsd.getFirstMapping());
                secondRoundTargets.add(
                        builder.newInstance(IAtomContainer.class, smsd.getProductMolecule()));
            }

            String name = inputHandler.getTargetName();
            outputHandler.writeCircleImage(mcsMolecule, secondRoundTargets, name, mappings);
        }
    }

    public static void runSingleQueryMultipleTarget(
            InputHandler inputHandler,
            OutputHandler outputHandler,
            ArgumentHandler argumentHandler) throws IOException, CDKException, CloneNotSupportedException {
        IMolecule query = inputHandler.getQuery();
        outputHandler.writeQueryMol(query);

        String out = ".out";
        outputHandler.startAppending(out);

        long startTime = System.currentTimeMillis();
        IAtomAtomMapping smsd = null;
        Substructure substructure = null;
        Isomorphism isomorphism = null;
        boolean matchBonds = argumentHandler.isMatchBondType();
        if (argumentHandler.isSubstructureMode()) {
            substructure = new Substructure();
        } else {
            isomorphism = new Isomorphism(Algorithm.DEFAULT, matchBonds);
        }

        boolean removeHydrogens = argumentHandler.isApplyHRemoval();
        int targetNumber = 0;
        IIteratingChemObjectReader reader = inputHandler.getAllTargets();
        String targetType = argumentHandler.getTargetType();
        if (reader == null) {
            throw new IOException("Unknown input type " + targetType);
        }
        while (reader.hasNext()) {
            IMolecule target = (IMolecule) reader.next();
            inputHandler.configure(target, targetType);
            if (argumentHandler.isSubstructureMode()) {
                runSubstructure(substructure, query, target, argumentHandler.getChemFilter(), matchBonds);
                smsd = substructure;
            } else {
                run(isomorphism, query, target, argumentHandler.getChemFilter());
                smsd = isomorphism;
            }
            long endTime = System.currentTimeMillis();
            long executionTime = endTime - startTime;
            outputHandler.writeTargetMol(smsd.getProductMolecule());

            String queryPath = argumentHandler.getQueryFilepath();
            String targetPath = argumentHandler.getTargetFilepath();

            query = query.getBuilder().newInstance(IMolecule.class, smsd.getReactantMolecule());
            target = target.getBuilder().newInstance(IMolecule.class, smsd.getProductMolecule());


            Map<IAtom, IAtom> mcs = smsd.getFirstAtomMapping();
            int nAtomsMatched = (mcs == null) ? 0 : mcs.size();
            double tanimotoSimilarity = smsd.getTanimotoSimilarity();

            //print out all mappings
            if (mcs != null && argumentHandler.isAllMapping()) {
                outputHandler.printHeader(queryPath, targetPath, nAtomsMatched);
                List<Map<Integer, Integer>> allMappings = smsd.getAllMapping();
                int counter = 1;
                for (Map<Integer, Integer> mapping : allMappings) {
                    if (argumentHandler.isImage()) {
                        double stereoScore = smsd.getStereoScore(counter);
                        String label = outputHandler.makeLabel(tanimotoSimilarity, stereoScore);
                        outputHandler.addImage(query, target, label, mapping);
                    }
                    outputHandler.printMapping(counter, mcs);
                }
            } //print out top one
            else if (mcs != null && !argumentHandler.isAllMapping()) {
                Map<Integer, Integer> mcsNumber = smsd.getFirstMapping();
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
        IMolecule query = inputHandler.getQuery();
        IMolecule target = inputHandler.getTarget();

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
        IAtomAtomMapping smsd = null;
        Substructure substructure = null;
        Isomorphism isomorphism = null;
        boolean matchBonds = argumentHandler.isMatchBondType();
        if (argumentHandler.isSubstructureMode()) {
            substructure = new Substructure();
        } else {
            isomorphism = new Isomorphism(Algorithm.DEFAULT, matchBonds);
        }
        boolean removeHydrogens = argumentHandler.isApplyHRemoval();
        if (argumentHandler.isSubstructureMode()) {
            runSubstructure(substructure, query, target, argumentHandler.getChemFilter(), matchBonds);
            smsd = substructure;
        } else {
            run(isomorphism, query, target, argumentHandler.getChemFilter());
            smsd = isomorphism;
        }

        query = query.getBuilder().newInstance(IMolecule.class, smsd.getReactantMolecule());
        target = target.getBuilder().newInstance(IMolecule.class, smsd.getProductMolecule());

        long endTime = System.currentTimeMillis();
        long executionTime = endTime - startTime;

        // write out the input molecules to files
        outputHandler.writeQueryMol(smsd.getReactantMolecule());
        outputHandler.writeTargetMol(smsd.getProductMolecule());

        String queryPath = argumentHandler.getQueryFilepath();
        String targetPath = argumentHandler.getTargetFilepath();

        Map<IAtom, IAtom> mcs = smsd.getFirstAtomMapping();
        int nAtomsMatched = (mcs == null) ? 0 : mcs.size();
        double tanimotoSimilarity = smsd.getTanimotoSimilarity();

        //print out all mappings
        if (mcs != null && argumentHandler.isAllMapping()) {
            outputHandler.printHeader(queryPath, targetPath, nAtomsMatched);
            List<Map<Integer, Integer>> allMappings = smsd.getAllMapping();
            int counter = 1;
            for (Map<Integer, Integer> mapping : allMappings) {
                if (argumentHandler.isImage()) {
                    double stereoScore = smsd.getStereoScore(counter);
                    String label = outputHandler.makeLabel(tanimotoSimilarity, stereoScore);
                    outputHandler.addImage(query, target, label, mapping);
                }
                outputHandler.printMapping(counter, mcs);
            }
        } //print out top one
        else if (mcs != null && !argumentHandler.isAllMapping()) {
            Map<Integer, Integer> mcsNumber = smsd.getFirstMapping();
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
            Map<Integer, Integer> mapping = smsd.getFirstMapping();
            IAtomContainer subgraph = getSubgraph(target, mapping);
            String outpath = argumentHandler.getOutputFilepath();
            String outtype = argumentHandler.getOutputFiletype();
            outputHandler.writeMol(outtype, subgraph, outpath);
        }
        outputHandler.closeFiles();
    }

    private static IAtomContainer getSubgraph(
            IMolecule container, Map<Integer, Integer> mapping) throws CloneNotSupportedException {
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

    private static void run(Isomorphism smsd,
            IMolecule query,
            IMolecule target,
            int filter) throws CDKException {
        // XXX - if clean and configure is 'true', is that not duplicate configuring?
        smsd.init(query, target);
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

    private static void runSubstructure(Substructure smsd,
            IMolecule query,
            IMolecule target,
            int filter,
            boolean matchBonds) throws CDKException {
        // XXX - if clean and configure is 'true', is that not duplicate configuring?
        smsd.init(query, target);
        boolean findSubgraph = smsd.findSubgraph(matchBonds);
        if (findSubgraph) {
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
    }
}
