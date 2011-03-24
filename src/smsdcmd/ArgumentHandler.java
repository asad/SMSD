/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package smsdcmd;

import java.io.Writer;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.cli.PosixParser;

/**
 *
 * @author sar
 */
public class ArgumentHandler {

   
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
     * @return the graphFile
     */
    public String getGraphFile() {
        return graphFile;
    }

    /**
     * @return the descriptorFile
     */
    public String getDescriptorFile() {
        return descriptorFile;
    }
    private boolean applyHAdding = false;
    private boolean applyHRemoval = false;
    private boolean applyTest = false;
    private boolean applySuffix = false;
    private boolean appendMode = false;
    private boolean matchBondType = false;
    private String matchFile = "mcs";
    private String fingerFile = "finger";
    private String graphFile = "graph";
    private String descriptorFile = "molDescriptors";
    private String queryOutfileName = "Query";
    private String targetOutfileName = "Target";
    private String suffix = "";
    private int chemFilter = 3;
    private boolean allMapping = false;
    private boolean image = false;
    private boolean substructureMode = false;
    private boolean isNMCS = false;
    private boolean outputSubgraph;
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

    protected void printError(String errmessg) {
        System.err.println(errmessg);
        System.exit(-1);
    }

    /**
     * Parses the options in the command line arguments and returns
     * an array of strings corresponding to the filenames given as arguments only
     * @param args
     * @return
     * @throws org.apache.commons.cli.ParseException 
     */
    @SuppressWarnings("static-access")
    protected void parseCommandLineOptions(String[] args) throws ParseException {

        options = new Options();

        options.addOption("h", "help", false, "Help page for command usage");

        options.addOption("s", false, "SubStructure detection");
        
        options.addOption("a", false, "Add Hydrogen");

        options.addOption("r", false, "Remove Hydrogen");

        options.addOption("b", false, "Match Bond types (Single, Double etc)");

        options.addOption(
                OptionBuilder.hasArg().withDescription("Query filename").withArgName("filepath").create("q")
        );

        options.addOption(
                OptionBuilder.hasArg().withDescription("Target filename").withArgName("filepath").create("t")
        );

        options.addOption(
                OptionBuilder.hasArg().withDescription("Add suffix to the files").withArgName("suffix").create("S")
        );

        options.addOption("g", false, "create png of the mapping");
        
        options.addOption(
                OptionBuilder.hasArg().withDescription("Dimension of the image in pixels").withArgName("WIDTHxHEIGHT").create("d")
        );

        options.addOption("m", false, "Report all Mappings");

        String filterDescr = "Default: 0, Stereo: 1, "
                                 + "Stereo+Fragment: 2, Stereo+Fragment+Energy: 3";
        options.addOption(
                OptionBuilder.hasArg().withDescription(filterDescr).withArgName("number").create("f")
        );

        options.addOption("A", false, 
             "Appends output to existing files, else creates new files");
        
        options.addOption(
                OptionBuilder.withDescription("Do N-way MCS on the target SD file").create("N")
        );
        
        options.addOption(
                OptionBuilder.hasArg().withDescription("Query type (MOL, SMI, etc)").withArgName("type").create("Q")
        );
        
        options.addOption(
                OptionBuilder.hasArg().withDescription("Target type (MOL, SMI, etc)").withArgName("type").create("T")
        );
        
        options.addOption(
                OptionBuilder.hasArg().withDescription("Output the substructure to a file").withArgName("filename").create("o")
        );
        
        options.addOption(
                OptionBuilder.hasArg().withDescription("Output type").withArgName("type").create("O")
        );

        PosixParser parser = new PosixParser();
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
        
        remainingArgs = line.getArgs();

        if (line.hasOption('h') || line.getOptions().length == 0) {
            System.out.println("Hello");
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
       
    }
    
    public Options getOptions() {
        return options;
    }
    
    public void printHelp() {
        printHelp(options);
    }
    
    public void printHelp(Options options) {
        HelpFormatter formatter = new HelpFormatter();

        System.out.println("\n++++++++++++++++++++++++++++++++++++++++++++++"
                + "\nSMSD (Small Molecule Similarity Detector)"
                + "\n++++++++++++++++++++++++++++++++++++++++++++++"
                + "\nContact: Syed Asad Rahman,"
                + "\n\t EMBL-EBI, Hinxton "
                + "\n\t Cambridge CB10 1SD"
                + "\n\t United Kingdom "
                + "\ne-mail: asad@ebi.ac.uk"
                + "\n++++++++++++++++++++++++++++++++++++++++++++++\n"
                + "\nSMSD software can calculate the similarity between"
                + "\ntwo small molecules by using an inhouse algorithm"
                + "\ndeveloped at EMBL-EBI. "
                + "It also uses CDK based"
                + "\nFingerprints and Maximum Common Subgraph (MCS) "
                + "\nmatching for finding similarity. "
                + "It create four"
                + "\nfiles which contains similarity scores between"
                + "\nmatched atoms.\n"
                + "\n++++++++++++++++++++++++++++++++++++++++++++++"
                + "\nNOTE: The graph matching is performed by removing"
                + "\nthe Hydrogens"
                + "\n++++++++++++++++++++++++++++++++++++++++++++++\n");

        formatter.printHelp("\n", options);
        System.out.println("\n++++++++++++++++++++++++++++++++++++++++++++++\n");
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
     * @return
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
     * @param matchBondType the matchBondType to set
     */
    public void setMatchBondType(boolean matchBondType) {
        this.matchBondType = matchBondType;
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
        return this.targetOutfileName;
    }

    /**
     * @return the Query
     */
    public String getQueryMolOutName() {
        return this.queryOutfileName;
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
    
}
