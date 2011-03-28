package cmd;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

import org.apache.commons.cli.MissingOptionException;
import org.openscience.cdk.ChemFile;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.Molecule;
import org.openscience.cdk.aromaticity.CDKHueckelAromaticityDetector;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IChemFile;
import org.openscience.cdk.interfaces.IMolecule;
import org.openscience.cdk.io.CMLReader;
import org.openscience.cdk.io.IChemObjectReader;
import org.openscience.cdk.io.ISimpleChemObjectReader;
import org.openscience.cdk.io.MDLReader;
import org.openscience.cdk.io.Mol2Reader;
import org.openscience.cdk.io.PDBReader;
import org.openscience.cdk.io.iterator.IIteratingChemObjectReader;
import org.openscience.cdk.io.iterator.IteratingMDLReader;
import org.openscience.cdk.layout.StructureDiagramGenerator;
import org.openscience.cdk.nonotify.NoNotificationChemObjectBuilder;
import org.openscience.cdk.signature.MoleculeSignature;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import org.openscience.cdk.tools.manipulator.ChemFileManipulator;

import cmd.pdb.LigandHelper;

public class InputHandler {
    
    private ArgumentHandler argumentHandler;
    private StructureDiagramGenerator sdg;
    private Map<String, String> singularDataTypes;
    private Map<String, String> multipleDataTypes;
    private Map<String, String> stringDataTypes;
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
        
        singularDataTypes = new HashMap<String, String>();
        singularDataTypes.put("CML", "Chemical Markup Language");
        singularDataTypes.put("MOL", "MDL V2000 format");
        singularDataTypes.put("ML2", "MOL2 Tripos format");
        singularDataTypes.put("PDB", "Protein Databank Format");
        
        multipleDataTypes = new HashMap<String, String>();
        multipleDataTypes.put("SDF", "SD file format");
        
        stringDataTypes = new HashMap<String, String>();
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
        return stringDataTypes;
    }
    
    public Map<String, String> getSingularDataTypes() {
        return singularDataTypes;
    }
    
    public Map<String, String> getMultipleDataTypes() {
        return multipleDataTypes;
    }
    
    public String getQRefName() {
        String suffix = argumentHandler.getSuffix(); 
        return argumentHandler.getQueryMolOutName() + suffix + ".mol";
    }
    
    public String getTRefName() {
        String suffix = argumentHandler.getSuffix(); 
        return argumentHandler.getTargetMolOutName() + suffix + ".mol";
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
//        System.err.println("MatchType " + matchType + 
//                " isStringQuery " + isStringQuery +
//                " isSingleFileQuery " + isSingleFileQuery +
//                " isSingleFileTarget " + isSingleFileTarget +
//                " isMultipleTarget " + isMultipleTarget +
//                " queryType " + argumentHandler.getQueryType() +
//                " targetType " + argumentHandler.getTargetType()
//                );
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
        if (type.equals("MOL")) {
            return new MDLReader(
                    new FileReader(input), IChemObjectReader.Mode.RELAXED);
        } else if (type.equals("CML")) {
            return new CMLReader(new FileInputStream(input));
        } else if (type.equals("ML2")) {
            return new Mol2Reader(new FileReader(input));
        } else if (type.equals("PDB")) {
            return new PDBReader(new FileReader(input));
        }
        return null;
    }
    
    public void configure(IMolecule mol, String type) throws CDKException {
        if (type.equals("PDB")) {
            LigandHelper.addMissingBondOrders(mol);
        }
        AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(mol);
        CDKHueckelAromaticityDetector.detectAromaticity(mol);
        sdg.setMolecule(new Molecule(mol), false);
        sdg.generateCoordinates();
        setID(mol);
    }
    
    private IMolecule getMolFromString(String stringData, String type) throws CDKException {
        if (type.equals("SMI")) {
            return getMolFromSmiles(stringData);
        } else if (type.equals("SIG")) {
            return getMolFromSignature(stringData);
        } else {
            return null;
        }
    }
    
    private IMolecule getMolFromSmiles(String smiles) throws CDKException {
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer atomContainer = sp.parseSmiles(smiles);
        IMolecule mol = new Molecule(atomContainer);
        configure(mol, "SMI");
        return mol;
    }
    
    private IMolecule getMolFromSignature(String signatureString) throws CDKException {
        IAtomContainer atomContainer = MoleculeSignature.fromSignatureString(
                signatureString, DefaultChemObjectBuilder.getInstance());
        IMolecule mol = new Molecule(atomContainer);
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
    
    public IMolecule getQuery() throws IOException, CDKException {
        String filenameOrData = argumentHandler.getQueryFilepath();
        String type = argumentHandler.getQueryType(); 
        if (isSingleFileQuery) {
            ISimpleChemObjectReader reader = getReader(type, filenameOrData);
            IChemFile chemFile = reader.read(new ChemFile());
            IMolecule molecule = 
                (IMolecule)ChemFileManipulator.getAllAtomContainers(chemFile).get(0);
            configure(molecule, type);
            return molecule;
        } else {
            return getMolFromString(filenameOrData, type);
        }
    }
    
    public IMolecule getTarget() throws IOException, CDKException {
        String filenameOrData = argumentHandler.getTargetFilepath();
        String type = argumentHandler.getTargetType();
        if (isSingleFileTarget) {
            ISimpleChemObjectReader reader = getReader(type, filenameOrData);
            IChemFile chemFile = reader.read(new ChemFile());
            IMolecule molecule = 
                (IMolecule)ChemFileManipulator.getAllAtomContainers(chemFile).get(0);
            configure(molecule, type);
            return molecule;
        } else {
            return getMolFromString(filenameOrData, type);
        }
    }
    
    public IIteratingChemObjectReader getAllTargets() throws FileNotFoundException {
        String type = argumentHandler.getTargetType();
        if (type.equals("SDF")) {
            FileReader in = new FileReader(argumentHandler.getTargetFilepath());
            return new IteratingMDLReader(
                        in, NoNotificationChemObjectBuilder.getInstance());
        }
        return null;
    }
    
    private static void setID(IMolecule mol) {
        int index = 1;

        for (IAtom atom : mol.atoms()) {
            atom.setID(String.valueOf(index++));
        }
    }

}
