/*
 * SPDX-License-Identifier: Apache-2.0
 * © 2025 BioInception PVT LTD.
 */
package uk.ac.ebi.smsd.cmd;

import java.io.Writer;
import java.util.Properties;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

/**
 *
 * @author Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>
 */
public class ArgumentHandler {

    final static String NEW_LINE = System.getProperty("line.separator");

    /**
     * Get any arguments that were left over.
     *
     * @return the remaining arguments
     */
    public String[] getRemainingArgs() {
        return remainingArgs;
    }

    /**
     * @return the matchFile
     */
    public String getMatchFile() {
        return matchFile;
    }

    /**
     * @return the fingerFile
     */
    public String getFingerFile() {
        return fingerFile;
    }

    /**
     * @return the descriptorFile
     */
    public String getDescriptorFile() {
        return descriptorFile;
    }
    /*
     Boolean Options
     */
    private boolean applyHAdding = false;
    private boolean applyHRemoval = false;
    private boolean applyTest = false;
    private boolean applySuffix = false;
    private boolean appendMode = false;
    private boolean matchBondType = false;
    private boolean matchRingType = false;

    private boolean allMapping = false;
    private boolean image = false;
    private boolean substructureMode = false;
    private boolean isNMCS = false;
    private boolean outputSubgraph = false;
    private boolean matchAtomType = false;

    /*
    
     */
    private final String matchFile = "mcs";
    private final String fingerFile = "finger";
    private final String descriptorFile = "molDescriptors";
    private String queryOutfileName;
    private String targetOutfileName;
    private String suffix = "";
    private int chemFilter = 3;
    private String outputFilepath;
    private Writer outputWriter;
    private String outputFiletype;
    private String queryFilepath;
    private String targetFilepath;
    private String queryType;
    private String targetType;
    private int imageWidth = -1;
    private int imageHeight = -1;
    private boolean helpRequested = false;
    private Options options;
    private String[] remainingArgs;
    private Properties imageProperties;
    private boolean isImageOptionHelp = false;

    protected void printError(String errmessg) {
        System.err.println(errmessg);
        System.exit(-1);
    }

    /**
     * Parses the options in the command line arguments and returns an array of
     * strings corresponding to the filenames given as arguments only
     *
     * @param args
     * @throws org.apache.commons.cli.ParseException
     */
    @SuppressWarnings("static-access")
    public void parseCommandLineOptions(String[] args) throws ParseException {

        options = new Options();

        options.addOption("h", "help", false, "Help page for command usage");

        options.addOption("s", false, "SubStructure detection");

        options.addOption("a", false, "Add Hydrogen");

        options.addOption("x", false, "Match Atom Type");

        options.addOption("r", false, "Remove Hydrogen");

        options.addOption("z", false, "Ring matching");

        options.addOption("b", false, "Match Bond types (Single, Double etc)");

        options.addOption(
                Option.builder("q")
                .argName("query_file")
                .hasArg()
                .desc("Query filename")
                .build());

        options.addOption(
                Option.builder("t")
                .argName("target_file")
                .hasArg()
                .desc("Target filename")
                .build());
        options.addOption(
                Option.builder("S")
                .argName("suffix")
                .hasArg()
                .desc("Add suffix to the files")
                .build());

        options.addOption("g", false, "create png of the mapping");

        options.addOption(
                Option.builder("d")
                .argName("WIDTHxHEIGHT")
                .hasArg()
                .desc("Dimension of the image in pixels")
                .build());

        options.addOption("m", false, "Report all Mappings");

        String filterDescr = "Default: 0, Stereo: 1, "
                + "Stereo+Fragment: 2, Stereo+Fragment+Energy: 3";

        options.addOption(
                Option.builder("f")
                .argName("filter_number")
                .hasArg()
                .desc(filterDescr)
                .build());

        options.addOption("A", false,
                "Appends output to existing files, else creates new files");

        options.addOption(
                Option.builder("N")
                .desc("Do N-way MCS on the target SD file")
                .build());

        options.addOption(
                Option.builder("Q")
                .argName("query_type")
                .hasArg()
                .desc("Query type (MOL, SMI, etc.)")
                .build());

        options.addOption(
                Option.builder("T")
                .argName("target_type")
                .hasArg()
                .desc("Target type (MOL, SMI, SMIF, etc.)")
                .build());

        options.addOption(
                Option.builder("o")
                .argName("filename")
                .hasArg()
                .desc("Output the substructure to a file")
                .build());

        options.addOption(
                Option.builder("O")
                .argName("type")
                .hasArg()
                .desc("Output type (SMI, MOL)")
                .build());

        options.addOption(
                Option.builder("I")
                .argName("option=value")
                .hasArg()
                .desc("Image options")
                .build());

        DefaultParser parser = new DefaultParser();
        CommandLine line = parser.parse(options, args, true);

        if (line.hasOption('Q')) {
            queryType = line.getOptionValue("Q");
        } //else {
//            queryType = "MOL";
//        } //XXX default type?

        if (line.hasOption('T')) {
            targetType = line.getOptionValue("T");
        } else {
            targetType = "MOL";
        }

        if (line.hasOption('a')) {
            this.setApplyHAdding(true);
        }

        if (line.hasOption('r')) {
            this.setApplyHRemoval(true);
        }

        if (line.hasOption('m')) {
            this.setAllMapping(true);
        }

        if (line.hasOption('s')) {
            this.setSubstructureMode(true);
        }

        if (line.hasOption('g')) {
            this.setImage(true);
        }

        if (line.hasOption('b')) {
            this.setMatchBondType(true);
        }

        if (line.hasOption('z')) {
            this.setMatchRingType(true);
        }

        if (line.hasOption('x')) {
            this.setMatchAtomType(true);
        }

        remainingArgs = line.getArgs();

        if (line.hasOption('h') || line.getOptions().length == 0) {
//            System.out.println("Hello");
            helpRequested = true;
        }

        if (line.hasOption('S')) {
            String[] suffix_reader = line.getOptionValues('S');
            if (suffix_reader.length < 1) {
                System.out.println("Suffix required!");
                helpRequested = true;
            }
            setSuffix(suffix_reader[0]);
            setApplySuffix(true);
        }

        if (line.hasOption('f')) {
            String[] filters = line.getOptionValues('f');
            if (filters.length < 1) {
                System.out.println("Chemical filter required (Ranges: 0 to 3)!");
                helpRequested = true;
            }
            setChemFilter((int) new Integer(filters[0]));
        }

        if (line.hasOption('q')) {
            queryFilepath = line.getOptionValue('q');
        }

        if (line.hasOption('t')) {
            targetFilepath = line.getOptionValue('t');
        }

        if (line.hasOption("A")) {
            this.setAppendMode(true);
        }

        if (line.hasOption("N")) {
            setNMCS(true);
        }

        if (line.hasOption("o")) {
            outputSubgraph = true;
            outputFilepath = line.getOptionValue("o");
        }

        if (line.hasOption("O")) {
            outputFiletype = line.getOptionValue("O");
        } else {
            outputFiletype = "MOL";
        }

        if (line.hasOption("d")) {
            String dimensionString = line.getOptionValue("d");
            if (dimensionString.contains("x")) {
                String[] parts = dimensionString.split("x");
                try {
                    setImageWidth(Integer.parseInt(parts[0]));
                    setImageHeight(Integer.parseInt(parts[1]));
                    System.out.println("set image dim to "
                            + getImageWidth() + "x" + getImageHeight());
                } catch (NumberFormatException nfe) {
                    throw new ParseException("Malformed dimension string " + dimensionString);
                }
            } else {
                throw new ParseException("Malformed dimension string " + dimensionString);
            }
        }

        if (line.hasOption("I")) {
            imageProperties = line.getOptionProperties("I");
            if (imageProperties.isEmpty()) {
                // used just "-I" by itself
                isImageOptionHelp = true;
            }
        }
    }

    public Properties getImageProperties() {
        return imageProperties;
    }

    public void setImageProperties(Properties imageProperties) {
        this.imageProperties = imageProperties;
    }

    public Options getOptions() {
        return options;
    }

    public void printHelp() {
        printHelp(options);
    }

    public void printHelp(Options options) {
        HelpFormatter formatter = new HelpFormatter();

        System.out.println(NEW_LINE + "++++++++++++++++++++++++++++++++++++++++++++++"
                + NEW_LINE + "Small Molecule Similarity Detector (SMSD)"
                + NEW_LINE + "+++++++++++++++++++++++++++++++++++++++++++++"
                + NEW_LINE + "Contact: Syed Asad Rahman,"
                + NEW_LINE + "\t EMBL-EBI, Hinxton "
                + NEW_LINE + "\t Cambridge CB10 1SD"
                + NEW_LINE + "\t United Kingdom "
                + NEW_LINE + "e-mail:  asad.rahman@bioinceptionlabs.com"
                + NEW_LINE + "++++++++++++++++++++++++++++++++++++++++++++++" + NEW_LINE
                + NEW_LINE + "SMSD can calculate similarity between small molecules"
                + NEW_LINE + "by using an inhouse algorithm developed at EMBL-EBI. It"
                + NEW_LINE + "reports similarity scores and highlights the matched "
                + NEW_LINE + "substructure." + NEW_LINE
                + NEW_LINE + "It uses CDK library for handling chemistry. "
                + NEW_LINE + "++++++++++++++++++++++++++++++++++++++++++++++"
                + NEW_LINE + "Publication:"
                + NEW_LINE + "   S.A. Rahman et.al, (2009), Journal of Cheminformatics"
                + NEW_LINE + "   (http://www.jcheminf.com/content/1/1/12)"
                + NEW_LINE + "++++++++++++++++++++++++++++++++++++++++++++++" + NEW_LINE);

        formatter.printHelp(NEW_LINE, options);
        System.out.println(getExamples());
        System.out.println(NEW_LINE + "++++++++++++++++++++++++++++++++++++++++++++++" + NEW_LINE);
    }

    private StringBuilder getExamples() {

        StringBuilder sb = new StringBuilder();
        sb.append(NEW_LINE);
        sb.append("Usage examples: ");
        sb.append(NEW_LINE);
        sb.append("To start the GUI: java -jar SMSD.jar ").append(NEW_LINE).append(NEW_LINE);
        sb.append("Java Command: ").append("java -Xms500M -Xmx512M -cp smsd.jar: uk.ac.ebi.smsd.cmd.SMSDcmd").append(NEW_LINE).append(NEW_LINE);
        sb.append("Help Command: ").append("java -Xms500M -Xmx512M -cp smsd.jar: uk.ac.ebi.smsd.cmd.SMSDcmd -h").append(NEW_LINE).append(NEW_LINE);
        sb.append("Shell Command Line Options (recommended flags):")
                .append("\n -r -z -b will remove hydrogens, match rings, match bond types respectively.").append(NEW_LINE).append(NEW_LINE);
        sb.append("a) Find MCS between two molecules and write the output as a MOL file:").append(NEW_LINE)
                .append("\tsh SMSD.sh -Q SMI -q \"CCC\" -T PDB -t ATP.pdb -O MOL -o subgraph.mol -r -z -b").append(NEW_LINE);
        sb.append("b) Find MCS between two molecules and print the output:").append(NEW_LINE)
                .append("\tsh SMSD.sh -Q SMI -q \"CCC\" -T SMI -t \"CCN\" -O SMI -o - -r -z -b").append(NEW_LINE);
        sb.append("c) Find MCS between two molecules and generate an image highlighting the matched parts:").append(NEW_LINE)
                .append("\tsh SMSD.sh -Q MOL -q ADP.mol -T MOL -t ATP.mol -g -r -z -b").append(NEW_LINE);
        sb.append("d) Find MCS between N-molecules and print the common substructure between them:").append(NEW_LINE)
                .append("\tsh SMSD.sh -T SDF -t arom.sdf -N -O SMI -o - -r -z -b").append(NEW_LINE);
        sb.append("e) Find MCS between N-molecules and highlighting the common substructure between them:").append(NEW_LINE);
        sb.append("\tWARNING: This option might require large virtual machine memory allocation").append(NEW_LINE)
                .append("\tsh SMSD.sh -T SDF -t arom.sdf -N -O SMI -o - -g -r -z -b").append(NEW_LINE).append(NEW_LINE);
        sb.append("Note: You could use various file formats").append(NEW_LINE);
        return sb;
    }

    public boolean isHelp() {
        return helpRequested;
    }

    public void setOutputWriter(Writer writer) {
        outputWriter = writer;
    }

    public Writer getOutputWriter() {
        return outputWriter;
    }

    public boolean shouldOutputSubgraph() {
        return outputSubgraph;
    }

    public String getOutputFilepath() {
        return outputFilepath;
    }

    public String getOutputFiletype() {
        return outputFiletype;
    }

    public String getQueryFilepath() {
        return queryFilepath;
    }

    public void setQueryFilepath(String filepath) {
        queryFilepath = filepath;
    }

    public String getTargetFilepath() {
        return targetFilepath;
    }

    public void setTargetFilepath(String filepath) {
        targetFilepath = filepath;
    }

    /**
     * Use N-way MCS.
     *
     * @return True if condition meet else False
     */
    public boolean isNMCS() {
        return isNMCS;
    }

    /**
     * Set the use of N-MCS for finding the MCS of N molecules from an SDF file.
     *
     * @param value true if N-MCS should be used
     */
    public void setNMCS(boolean value) {
        isNMCS = value;
    }

    /**
     * @return the applyHAdding
     */
    public boolean isApplyHAdding() {
        return applyHAdding;
    }

    /**
     * @param applyHAdding the applyHAdding to set
     */
    public void setApplyHAdding(boolean applyHAdding) {
        this.applyHAdding = applyHAdding;
    }

    /**
     * @return the applyHRemoval
     */
    public boolean isApplyHRemoval() {
        return applyHRemoval;
    }

    /**
     * @param applyHRemoval the applyHRemoval to set
     */
    public void setApplyHRemoval(boolean applyHRemoval) {
        this.applyHRemoval = applyHRemoval;
    }

    /**
     * @return the applyTest
     */
    public boolean isApplyTest() {
        return applyTest;
    }

    /**
     * @param applyTest the applyTest to set
     */
    public void setApplyTest(boolean applyTest) {
        this.applyTest = applyTest;
    }

    /**
     * @return the applySuffix
     */
    public boolean isApplySuffix() {
        return applySuffix;
    }

    /**
     * @param applySuffix the applySuffix to set
     */
    public void setApplySuffix(boolean applySuffix) {
        this.applySuffix = applySuffix;
    }

    /**
     * @return the appendMode
     */
    public boolean isAppendMode() {
        return appendMode;
    }

    /**
     * @param appendMode the appendMode to set
     */
    public void setAppendMode(boolean appendMode) {
        this.appendMode = appendMode;
    }

    /**
     * @return the matchBondType
     */
    public boolean isMatchBondType() {
        return matchBondType;
    }

    /**
     * @return the matchRingType
     */
    public boolean isMatchRingType() {
        return matchRingType;
    }

    /**
     * @param matchBondType the matchBondType to set
     */
    public void setMatchBondType(boolean matchBondType) {
        this.matchBondType = matchBondType;
    }

    /**
     * @param matchRingType the match ring size and ring atom
     */
    public void setMatchRingType(boolean matchRingType) {
        this.matchRingType = matchRingType;
    }

    /**
     * @return the suffix
     */
    public String getSuffix() {
        return suffix;
    }

    /**
     * @param suffix the suffix to set
     */
    public void setSuffix(String suffix) {
        this.suffix = suffix;
    }

    /**
     * @return the chemFilter
     */
    public int getChemFilter() {
        return chemFilter;
    }

    /**
     * @return the Target
     */
    public String getTargetMolOutName() {
        return targetOutfileName;
    }

    /**
     * @return the Query
     */
    public String getQueryMolOutName() {
        return queryOutfileName;
    }

    /**
     * @param name
     */
    public void setTargetMolOutName(String name) {
        this.targetOutfileName = name;
    }

    /**
     * @param name
     */
    public void setQueryMolOutName(String name) {
        this.queryOutfileName = name;
    }

    /**
     * @param chemFilter the chemFilter to set
     */
    public void setChemFilter(int chemFilter) {
        this.chemFilter = chemFilter;
    }

    /**
     * @return the allMapping
     */
    public boolean isAllMapping() {
        return allMapping;
    }

    /**
     * @param allMapping the allMapping to set
     */
    public void setAllMapping(boolean allMapping) {
        this.allMapping = allMapping;
    }

    /**
     * @return the image
     */
    public boolean isImage() {
        return image;
    }

    /**
     * @param image the image to set
     */
    public void setImage(boolean image) {
        this.image = image;
    }

    /**
     * @return the substructureMode
     */
    public boolean isSubstructureMode() {
        return substructureMode;
    }

    /**
     * @param substructureMode the substructureMode to set
     */
    public void setSubstructureMode(boolean substructureMode) {
        this.substructureMode = substructureMode;
    }

    public String getQueryType() {
        return queryType;
    }

    public String getTargetType() {
        return targetType;
    }

    public void setQueryType(String queryType) {
        this.queryType = queryType;
    }

    public void setTargetType(String targetType) {
        this.targetType = targetType;
    }

    public void setOutputSubgraph(boolean value) {
        outputSubgraph = value;
    }

    public void setOutputFilepath(String filepath) {
        outputFilepath = filepath;
    }

    public void setOutputFiletype(String filetype) {
        outputFiletype = filetype;
    }

    public int getImageWidth() {
        return imageWidth;
    }

    public void setImageWidth(int imageWidth) {
        this.imageWidth = imageWidth;
    }

    public int getImageHeight() {
        return imageHeight;
    }

    public void setImageHeight(int imageHeight) {
        this.imageHeight = imageHeight;
    }

    public boolean isImageOptionHelp() {
        return isImageOptionHelp;
    }

    private void setMatchAtomType(boolean x) {
        this.matchAtomType = x;
    }

    /**
     * @return the matchAtomType
     */
    public boolean isMatchAtomType() {
        return matchAtomType;
    }
}
