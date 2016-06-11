/*
 *
 *
 * Copyright (C) 2009-2015  Syed Asad Rahman <asad@ebi.ac.uk>
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
 * You should have received query copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA.
 * 
 */
package uk.ac.ebi.smsd.cmd;

import uk.ac.ebi.smsd.cmd.pdb.LigandHelper;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import org.apache.commons.cli.MissingOptionException;
import org.openscience.cdk.AtomContainer;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.ChemFile;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.aromaticity.Kekulization;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.graph.ConnectivityChecker;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IChemFile;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.io.CMLReader;
import org.openscience.cdk.io.IChemObjectReader;
import org.openscience.cdk.io.ISimpleChemObjectReader;
import org.openscience.cdk.io.MDLV2000Reader;
import org.openscience.cdk.io.Mol2Reader;
import org.openscience.cdk.io.PDBReader;
import org.openscience.cdk.io.ReaderFactory;
import org.openscience.cdk.io.SMILESReader;
import org.openscience.cdk.io.iterator.IteratingSDFReader;
import org.openscience.cdk.layout.StructureDiagramGenerator;
import org.openscience.cdk.signature.MoleculeSignature;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.tools.CDKHydrogenAdder;
import org.openscience.cdk.tools.ILoggingTool;
import org.openscience.cdk.tools.LoggingToolFactory;
import org.openscience.cdk.tools.manipulator.ChemFileManipulator;
import org.openscience.smsd.tools.ExtAtomContainerManipulator;

/**
 *
 * @author Syed Asad Rahman <asad@ebi.ac.uk>
 */
public class InputHandler {

    private final static ILoggingTool logger
            = LoggingToolFactory.createLoggingTool(InputHandler.class);
    private final ArgumentHandler argumentHandler;
    private final StructureDiagramGenerator sdg;
    private final Map<String, String> singularDataTypes;
    private final Map<String, String> multipleDataTypes;
    private final Map<String, String> stringDataTypes;
    private boolean isSingleFileQuery;
    private boolean isSingleFileTarget;
    private boolean isMultipleTarget;
    private boolean isStringQuery;
    private boolean isStringTarget;

    public enum MatchType {

        SINGLE_QUERY_SINGLE_TARGET,
        SINGLE_QUERY_MULTIPLE_TARGET,
        NMCS,
        UNKNOWN
    };
    private MatchType matchType;

    public InputHandler(ArgumentHandler argumentHandler) {
        this.argumentHandler = argumentHandler;
        sdg = new StructureDiagramGenerator();

        singularDataTypes = new HashMap<>();
        singularDataTypes.put("CML", "Chemical Markup Language");
        singularDataTypes.put("MOL", "MDL V2000 format");
        singularDataTypes.put("ML2", "MOL2 Tripos format");
        singularDataTypes.put("PDB", "Protein Databank Format");

        multipleDataTypes = new HashMap<>();
        multipleDataTypes.put("SDF", "SD file format");
        multipleDataTypes.put("SMIF", "SMILES file format");

        stringDataTypes = new HashMap<>();
        stringDataTypes.put("SMI", "SMILES string format");
        stringDataTypes.put("SIG", "Signature string format");
    }

    public void printDataTypeHelp() {
        System.out.println("Allowed types for single-molecules (query or target):");
        for (String singularType : singularDataTypes.keySet()) {
            String description = singularDataTypes.get(singularType);
            System.out.println(String.format("%s\t%s", singularType, description));
        }
        for (String stringType : stringDataTypes.keySet()) {
            String description = stringDataTypes.get(stringType);
            System.out.println(String.format("%s\t%s", stringType, description));
        }
        System.out.println();
        System.out.println("Allowed types for multiple-molecules (targets only):");
        for (String multipleType : multipleDataTypes.keySet()) {
            String description = multipleDataTypes.get(multipleType);
            System.out.println(String.format("%s\t%s", multipleType, description));
        }
    }

    public Map<String, String> getStringDataTypes() {
        return Collections.unmodifiableMap(stringDataTypes);
    }

    public Map<String, String> getSingularDataTypes() {
        return Collections.unmodifiableMap(singularDataTypes);
    }

    public Map<String, String> getMultipleDataTypes() {
        return Collections.unmodifiableMap(multipleDataTypes);
    }

    public String getQRefName() {
        String suffix = argumentHandler.getSuffix();
        String fileName = argumentHandler.getQueryMolOutName() == null ? "Query" : argumentHandler.getQueryMolOutName();
        if (!fileName.equals("Query")) {
            fileName = argumentHandler.getQueryMolOutName().equals("untitled") ? "Query" : argumentHandler.getQueryMolOutName();
        }
        return fileName + suffix + ".mol";
    }

    public String getTRefName() {
        String suffix = argumentHandler.getSuffix();
        String fileName = argumentHandler.getTargetMolOutName() == null ? "Target" : argumentHandler.getTargetMolOutName();
        if (!fileName.equals("Target")) {
            fileName = argumentHandler.getTargetMolOutName().equals("untitled") ? "Target" : argumentHandler.getTargetMolOutName();
        }
        return fileName + suffix + ".mol";
    }

    public MatchType validateInput() throws MissingOptionException {
        validateQueryType();
        validateTargetType();
        if ((isSingleFileQuery && isSingleFileTarget)
                || (isStringQuery && isSingleFileTarget)
                || (isSingleFileQuery && isStringTarget)
                || (isStringQuery && isStringTarget)) {
            matchType = MatchType.SINGLE_QUERY_SINGLE_TARGET;
        } else if ((isSingleFileQuery || isStringQuery) && isMultipleTarget) {
            matchType = MatchType.SINGLE_QUERY_MULTIPLE_TARGET;
        } else if (!isSingleFileQuery && isMultipleTarget) {
            if (argumentHandler.isNMCS()) {
                matchType = MatchType.NMCS;
            } else {
                throw new MissingOptionException("Set N-MCS to true");
            }
        } else {
            matchType = MatchType.UNKNOWN;
        }
        return matchType;
    }

    private void validateQueryType() {
        String queryType = argumentHandler.getQueryType();
        if (queryType != null) {
            queryType = queryType.toUpperCase();
            if (singularDataTypes.containsKey(queryType)) {
                isSingleFileQuery = true;
            } else if (multipleDataTypes.containsKey(queryType)) {
                // TODO : throw error!
            } else if (stringDataTypes.containsKey(queryType)) {
                isStringQuery = true;
                isSingleFileQuery = false;
            }
        } else {
            isSingleFileQuery = false;
            isStringQuery = false;
        }
    }

    private void validateTargetType() {
        String targetType = argumentHandler.getTargetType().toUpperCase();
        if (singularDataTypes.containsKey(targetType)) {
            isSingleFileTarget = true;
        } else if (multipleDataTypes.containsKey(targetType)) {
            isMultipleTarget = true;
        } else if (stringDataTypes.containsKey(targetType)) {
            isSingleFileTarget = false;
            isMultipleTarget = false;
            isStringTarget = true;
        } else {
            // TODO : throw error! - must have either a target
        }
    }

    public boolean isNMCSInput() {
        return false;   // TODO
    }

    public boolean isSingleFileQueryInput() {
        return isSingleFileQuery;
    }

    public boolean isSingleFileTargetInput() {
        return isSingleFileTarget;
    }

    public boolean isMultiTargetInput() {
        return isMultipleTarget;
    }

    private ISimpleChemObjectReader getReader(String type, String filename) throws IOException {
        File input = new File(filename);
        if (input.isDirectory()) {
            throw new IOException(
                    "Input path " + filename + " is a directory, not a file");
        }
        switch (type) {
            case "MOL":
                return new MDLV2000Reader(
                        new FileReader(input), IChemObjectReader.Mode.RELAXED);
            case "CML":
                return new CMLReader(new FileInputStream(input));
            case "ML2":
                return new Mol2Reader(new FileReader(input));
            case "PDB":
                return new PDBReader(new FileReader(input));
        }
        return null;
    }

    /**
     *
     * @param molecule
     * @param type
     * @throws CDKException
     */
    public void configure(IAtomContainer molecule, String type) throws CDKException {
        IAtomContainer mol = molecule;
        String id = "";
        switch (type) {
            case "PDB":
                LigandHelper.addMissingBondOrders(mol);
                break;
            case "SDF":
                id = (String) mol.getProperty(CDKConstants.TITLE);
                break;
        }
        ExtAtomContainerManipulator.aromatizeMolecule(mol);
        mol = new AtomContainer(mol);
        mol.setID(id);

        if (argumentHandler.isImage()) {
            sdg.setMolecule(mol, false);
            sdg.generateCoordinates();
        }
        setAtomID(mol);
    }

    private IAtomContainer getMolFromString(String stringData, String type) throws CDKException {
        switch (type) {
            case "SMI":
                return getMolFromSmiles(stringData);
            case "SIG":
                return getMolFromSignature(stringData);
            default:
                return null;
        }
    }

    private IAtomContainer getMolFromSmiles(String smiles) throws CDKException {
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer atomContainer = sp.parseSmiles(smiles);
        IAtomContainer mol = new AtomContainer(atomContainer);
        configure(mol, "SMI");
        return mol;
    }

    private IAtomContainer getMolFromSignature(String signatureString) throws CDKException {
        IAtomContainer atomContainer = MoleculeSignature.fromSignatureString(
                signatureString, DefaultChemObjectBuilder.getInstance());
        IAtomContainer mol = new AtomContainer(atomContainer);
        configure(mol, "SIG");
        return mol;
    }

    public String getQueryName() {
        String filename = argumentHandler.getQueryFilepath();
        File input = new File(filename);
        return input.getName().split("\\.")[0];
    }

    public String getTargetName() {
        String filename = argumentHandler.getTargetFilepath();
        File input = new File(filename);
        return input.getName().split("\\.")[0];
    }

    /**
     * Return Query molecule
     *
     * @return
     * @throws IOException
     * @throws CDKException
     */
    public IAtomContainer getQuery() throws IOException, CDKException {
        String filenameOrData = argumentHandler.getQueryFilepath();
        String type = argumentHandler.getQueryType();
        if (isSingleFileQuery) {
            ISimpleChemObjectReader reader = getReader(type, filenameOrData);
            IChemFile chemFile = reader.read(new ChemFile());
            List<IAtomContainer> allAtomContainers = ChemFileManipulator.getAllAtomContainers(chemFile);
            IAtomContainer molecule = null;
            for (IAtomContainer frag : allAtomContainers) {
                if (molecule == null || frag.getAtomCount() > molecule.getAtomCount()) {
                    molecule = frag;
                }
            }
            configure(molecule, type);
            return molecule;
        } else {
            return getMolFromString(filenameOrData, type);
        }
    }

    /**
     * Returns Target molecule
     *
     * @return
     * @throws IOException
     * @throws CDKException
     */
    public IAtomContainer getTarget() throws IOException, CDKException {
        String filenameOrData = argumentHandler.getTargetFilepath();
        String type = argumentHandler.getTargetType();
        if (isSingleFileTarget) {
            ISimpleChemObjectReader reader = getReader(type, filenameOrData);
            IChemFile chemFile = reader.read(new ChemFile());
            List<IAtomContainer> allAtomContainers = ChemFileManipulator.getAllAtomContainers(chemFile);
            IAtomContainer molecule = null;
            for (IAtomContainer frag : allAtomContainers) {
                if (molecule == null || frag.getAtomCount() > molecule.getAtomCount()) {
                    molecule = frag;
                }
            }
            configure(molecule, type);
            return molecule;
        } else {
            return getMolFromString(filenameOrData, type);
        }
    }

    /**
     * Returns an SDF files iterator
     *
     * @return
     * @throws FileNotFoundException
     * @throws IOException
     * @throws CDKException
     */
    public List<IAtomContainer> getAllTargets() throws FileNotFoundException, IOException, CDKException {
        String type = argumentHandler.getTargetType();
        IChemObjectBuilder builder = DefaultChemObjectBuilder.getInstance();
        ISimpleChemObjectReader reader;
        boolean deducebonds = false;

        ReaderFactory readerFactory = new ReaderFactory();
        IChemFile emptyChemFile;
        IChemFile chemFile;

        String infileName = argumentHandler.getTargetFilepath();
        File inputFile = new File(infileName);

        if (!inputFile.isFile()) {
            throw new FileNotFoundException("ERROR: Input File Not Found " + infileName);
        }

        List<IAtomContainer> allAtomContainers = new ArrayList<>();
        switch (type) {
            case "SDF":
                IteratingSDFReader iteratingSDFReader
                        = new IteratingSDFReader(
                                new FileReader(inputFile), DefaultChemObjectBuilder.getInstance());
                while (iteratingSDFReader.hasNext()) {
                    IAtomContainer mol = iteratingSDFReader.next();
                    String id = (String) mol.getProperty(CDKConstants.TITLE);
                    mol.setID(id);
                    allAtomContainers.add(mol);
                }
                iteratingSDFReader.close();
                break;
            case "SMIF":
                reader = new SMILESReader(new FileReader(inputFile));
                deducebonds = true;
                emptyChemFile = builder.newInstance(IChemFile.class);
                chemFile = reader.read(emptyChemFile);
                allAtomContainers = ChemFileManipulator.getAllAtomContainers(chemFile);
                break;
            default:
                reader = readerFactory.createReader(new FileReader(inputFile));
                emptyChemFile = builder.newInstance(IChemFile.class);
                chemFile = reader.read(emptyChemFile);
                allAtomContainers = ChemFileManipulator.getAllAtomContainers(chemFile);
                break;
        }
        if (!allAtomContainers.isEmpty()) {
            // Get Molecules
            List<IAtomContainer> atomContainerList = new ArrayList<>(allAtomContainers.size());

            CDKHydrogenAdder adder = CDKHydrogenAdder.getInstance(DefaultChemObjectBuilder.getInstance());
            for (int atomContainerNr = 0; atomContainerNr < allAtomContainers.size();) {
                IAtomContainer temp = allAtomContainers.get(atomContainerNr);
                IAtomContainer atomcontainerHFree = ExtAtomContainerManipulator.removeHydrogens(temp);
                ExtAtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(atomcontainerHFree);

                if (deducebonds) {
//                    DeduceBondSystemTool dbst = new DeduceBondSystemTool();
//                    atomcontainerHFree = dbst.fixAromaticBondOrders(atomcontainerHFree);
                    Kekulization.kekulize(atomcontainerHFree);
                }

                adder.addImplicitHydrogens(atomcontainerHFree);
                String index = String.valueOf((atomContainerNr + 1));
                boolean flag = ConnectivityChecker.isConnected(atomcontainerHFree);
                String title = atomcontainerHFree.getProperty(CDKConstants.TITLE) != null
                        ? (String) atomcontainerHFree.getProperty(CDKConstants.TITLE) : index;
                atomcontainerHFree.setProperty(CDKConstants.TITLE, index);
                if (!flag) {
                    System.err.println("WARNING : Skipping target AtomContainer "
                            + title + " as it is not connected.");
                    continue;
                } else {
                    if (title != null) {
                        atomcontainerHFree.setID(title);
                    }
                    argumentHandler.setTargetMolOutName(atomcontainerHFree.getID());
                }

                atomContainerList.add(atomContainerNr, atomcontainerHFree);
                atomContainerNr++;
            }
            allAtomContainers.clear();
            return atomContainerList;
        } else {
            return null;
        }
    }

    private static void setAtomID(IAtomContainer mol) {
        int index = 1;
        for (IAtom atom : mol.atoms()) {
            atom.setID(String.valueOf(index));
            index++;
        }
    }
}
