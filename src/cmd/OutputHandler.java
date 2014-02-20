/*
 *
 *
 * Copyright (C) 2009-2014  Syed Asad Rahman <asad@ebi.ac.uk>
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
package cmd;

import java.awt.image.RenderedImage;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.Writer;
import java.lang.reflect.Field;
import java.text.NumberFormat;
import java.util.List;
import java.util.Map;

import javax.imageio.ImageIO;

import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.io.MDLV2000Writer;
import org.openscience.cdk.io.SDFWriter;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.tools.ILoggingTool;
import org.openscience.cdk.tools.LoggingToolFactory;
import org.openscience.reactionblast.graphics.direct.Params;

/**
 * Writes the results of SMSD to text files and images.
 *
 * @author maclean
 *
 */
public class OutputHandler {

    private final static ILoggingTool logger
            = LoggingToolFactory.createLoggingTool(OutputHandler.class);
    private final ArgumentHandler argumentHandler;
    private BufferedWriter outGFile = null;
    private BufferedWriter outMFile = null;
    private BufferedWriter outDescriptorFile = null;
    private final ImageGenerator imageGenerator;
    private final NumberFormat nf;
    final static String NEW_LINE = System.getProperty("line.separator");

    OutputHandler(ArgumentHandler argumentHandler) {
        this.argumentHandler = argumentHandler;
        imageGenerator = new ImageGenerator(argumentHandler);

        ////set the format right for the Tanimoto score (only two digits printed)
        nf = NumberFormat.getInstance();
        nf.setMaximumFractionDigits(2);
        nf.setMinimumFractionDigits(2);
    }

    void printImageOptionsHelp() {
        Params params = imageGenerator.getParams();
        System.out.println("Image options, and default values");
        for (Field field : params.getClass().getFields()) {
            try {
                System.out.println(field.getName() + "=" + field.get(params));
            } catch (IllegalArgumentException e) {
                logger.error("ERROR: ", e);
            } catch (IllegalAccessException e) {
                logger.error("ERROR: ", e);
            }
        }
    }

    /**
     *
     * @param mol
     * @throws IllegalArgumentException
     * @throws IOException
     * @throws CDKException
     */
    void writeQueryMol(IAtomContainer mol) throws IllegalArgumentException, IOException, CDKException {
        String suffix = argumentHandler.getSuffix();
        String fileName = argumentHandler.getQueryMolOutName() == null ? "Query" : argumentHandler.getQueryMolOutName();
        if (!fileName.equals("Query")) {
            fileName = argumentHandler.getQueryMolOutName().equals("untitled") ? "Query" : argumentHandler.getQueryMolOutName();
        }
        String qRefName = fileName + suffix + ".mol";
        writeMolToMolfile(mol, qRefName);
    }

    /**
     *
     * @param mol
     * @throws IllegalArgumentException
     * @throws IOException
     * @throws CDKException
     */
    public void writeTargetMol(IAtomContainer mol) throws IllegalArgumentException, IOException, CDKException {
        String suffix = argumentHandler.getSuffix();
        String fileName = argumentHandler.getQueryMolOutName() == null ? "Target" : argumentHandler.getQueryMolOutName();
        if (!fileName.equals("Target")) {
            fileName = argumentHandler.getQueryMolOutName().equals("untitled") ? "Target" : argumentHandler.getQueryMolOutName();
        }
        String tRefName = fileName + suffix + ".mol";
        writeMolToMolfile(mol, tRefName);
    }

    /**
     *
     * @param out
     * @throws IOException
     */
    void startAppending(String out) throws IOException {
        String suffix = argumentHandler.getSuffix();
        outGFile = new BufferedWriter(new FileWriter(argumentHandler.getGraphFile() + suffix + out));
        outMFile = new BufferedWriter(new FileWriter(argumentHandler.getMatchFile() + suffix + out));
        outDescriptorFile = new BufferedWriter(new FileWriter(argumentHandler.getDescriptorFile() + suffix + out));
    }

    /**
     *
     * @param out
     * @throws IOException
     */
    void startNew(String out) throws IOException {
        String suffix = argumentHandler.getSuffix();
        outGFile = new BufferedWriter(new FileWriter(argumentHandler.getGraphFile() + suffix + out, true));
        outMFile = new BufferedWriter(new FileWriter(argumentHandler.getMatchFile() + suffix + out, true));

        File f = new File(argumentHandler.getDescriptorFile() + out);
        if (!f.exists()) {
            outDescriptorFile = new BufferedWriter(new FileWriter(argumentHandler.getDescriptorFile() + suffix + out, true));

            outDescriptorFile.write("Query" + "\t");
            outDescriptorFile.write("Target" + "\t");
            outDescriptorFile.write("Tanimoto (Sim.)" + "\t");
            outDescriptorFile.write("Tanimoto (Bond Sim.)" + "\t");
            outDescriptorFile.write("Tanimoto (Atom Sim.)" + "\t");
            outDescriptorFile.write("Euclidian (Dist.)" + "\t");
            outDescriptorFile.write("Cosine (Sim.)" + "\t");
            outDescriptorFile.write("Soergel (Dist.)" + "\t");
            outDescriptorFile.write("Query (Atom Count)" + "\t");
            outDescriptorFile.write("Target (Atom Count)" + "\t");
            outDescriptorFile.write("Query (Bond Count)" + "\t");
            outDescriptorFile.write("Target (Bond Count)" + "\t");
            outDescriptorFile.write("Match (Size)" + "\t");
            outDescriptorFile.write("Query (Wt.)" + "\t");
            outDescriptorFile.write("Target (Wt.)");
            outDescriptorFile.newLine();
        } else {
            outDescriptorFile = new BufferedWriter(new FileWriter(argumentHandler.getDescriptorFile() + suffix + out, true));
        }
    }

    /**
     *
     * @param queryMolInput
     * @param targetMolInput
     * @param tanimotoGraph
     * @throws IOException
     */
    void writeGraphScores(String queryMolInput, String targetMolInput, double tanimotoGraph) throws IOException {
        String graphScoresL = queryMolInput + "\t" + targetMolInput + "\t" + nf.format(tanimotoGraph);
        outGFile.write(graphScoresL);
        outGFile.newLine();
    }

//    public void writeResults(IAtomContainer mol1, IAtomContainer mol2,
//            double tanimotoGraph, double tanimotoAtom, double tanimotoBond, double euclidianGraph, int nAtomsMatched,
//            long executionTime) throws IOException, CloneNotSupportedException {
    /**
     *
     * @param mol1
     * @param mol2
     * @param tanimotoGraph
     * @param euclidianGraph
     * @param nAtomsMatched
     * @param executionTime
     * @throws IOException
     * @throws CloneNotSupportedException
     */
    void writeResults(IAtomContainer mol1, IAtomContainer mol2,
            double tanimotoGraph, double euclidianGraph, int nAtomsMatched,
            long executionTime) throws IOException, CloneNotSupportedException {
        String queryMolInput = argumentHandler.getQueryFilepath();
        String targetMolInput = argumentHandler.getTargetFilepath();

        double cosineGraph;
        double SoergelGraph;

        int mol1Size = mol1.getAtomCount();
        int mol2Size = mol2.getAtomCount();
        int mol1BondSize = mol1.getBondCount();
        int mol2BondSize = mol2.getBondCount();

        if (nAtomsMatched != 0) {
            cosineGraph = nAtomsMatched / Math.sqrt((double) mol1Size * (double) mol2Size);
            SoergelGraph = ((double) mol1Size + (double) mol2Size - 2 * nAtomsMatched)
                    / ((double) mol1Size + (double) mol2Size - nAtomsMatched);
        } else {
            tanimotoGraph = 0.0;
            euclidianGraph = 0.0;
            cosineGraph = 0.0;
            SoergelGraph = 0.0;
        }

        //to hold the graph matching descriptor
        String graphDescrptorL = queryMolInput + "\t" + targetMolInput;

        if (!argumentHandler.isAppendMode()) {
            outDescriptorFile.write(graphDescrptorL + " ");
            outDescriptorFile.write("Tanimoto (Sim.)= " + nf.format(tanimotoGraph) + " ");
            outDescriptorFile.write("Euclidian (Dist.)= " + nf.format(euclidianGraph) + " ");
            outDescriptorFile.write("Cosine (Sim.)= " + nf.format(cosineGraph) + " ");
            outDescriptorFile.write("Soergel (Dist.)= " + nf.format(SoergelGraph) + " ");
            outDescriptorFile.write("Query (Atom Count)= " + mol1Size + " ");
            outDescriptorFile.write("Target (Atom Count)= " + mol2Size + " ");
            outDescriptorFile.write("Query (Bond Count)= " + mol1BondSize + " ");
            outDescriptorFile.write("Target (Bond Count)= " + mol2BondSize + " ");
            outDescriptorFile.write("Match (Size)= " + nAtomsMatched + " ");
            outDescriptorFile.write("Time (ms):" + executionTime + " ");
            outDescriptorFile.newLine();
        } else {
            outDescriptorFile.write(graphDescrptorL + "\t");
            outDescriptorFile.write(nf.format(tanimotoGraph) + "\t");
            outDescriptorFile.write(nf.format(euclidianGraph) + "\t");
            outDescriptorFile.write(nf.format(cosineGraph) + "\t");
            outDescriptorFile.write(nf.format(SoergelGraph) + "\t");
            outDescriptorFile.write(mol1Size + "\t");
            outDescriptorFile.write(mol2Size + "\t");
            outDescriptorFile.write(mol1BondSize + "\t");
            outDescriptorFile.write(mol2BondSize + "\t");
            outDescriptorFile.write(nAtomsMatched + "\t");
            outDescriptorFile.write("Time (ms):" + executionTime + " ");
            outDescriptorFile.newLine();
            outGFile.flush();
            outMFile.flush();
            outDescriptorFile.flush();
        }
    }

    /**
     *
     * @param queryMolInput
     * @param targetMolInput
     * @param nAtomsMatched
     * @throws IOException
     */
    void printHeader(
            String queryMolInput, String targetMolInput, int nAtomsMatched) throws IOException {
        outMFile.write("AtomContainer 1=\t" + queryMolInput);
        outMFile.newLine();
        outMFile.write("AtomContainer 2=\t" + targetMolInput);
        outMFile.newLine();
        outMFile.write("Max atoms matched=\t" + nAtomsMatched);
        outMFile.newLine();
    }

    /**
     *
     * @param solutionIndex
     * @param mcs
     * @throws IOException
     */
    void printMapping(int solutionIndex, Map<Integer, Integer> mcs) throws IOException {
        outMFile.newLine();
        outMFile.write("Solution=\t" + solutionIndex);
        outMFile.newLine();
        for (Integer queryIndex : mcs.keySet()) {
            Integer targetIndex = mcs.get(queryIndex);
            outMFile.write((queryIndex + 1) + "\t" + (targetIndex + 1));
            outMFile.newLine();
        }

        outMFile.newLine();
        outMFile.write("//");
        outMFile.newLine();
    }

    /**
     *
     * @param tanimotoSimilarity
     * @param stereoScore
     * @return
     */
    String makeLabel(double tanimotoSimilarity, double stereoScore) {
        String tanimoto = nf.format(tanimotoSimilarity);
        String stereo = nf.format(stereoScore);
        return "Scores [" + "Tanimoto: " + tanimoto + ", Stereo: " + stereo + "]";
    }

    /**
     *
     * @param nAtomsMatched
     * @param mcs
     * @param mcsNumber
     * @param qrefName
     * @param trefName
     * @throws IOException
     */
    void printTopMapping(
            int nAtomsMatched, Map<IAtom, IAtom> mcs, Map<Integer, Integer> mcsNumber,
            String qrefName, String trefName) throws IOException {

        for (Map.Entry<IAtom, IAtom> map : mcs.entrySet()) {
            String queryIndex = map.getKey().getID();
            String targetIndex = map.getValue().getID();
            outMFile.write(queryIndex + "\t" + targetIndex);
            outMFile.newLine();
        }
        outMFile.newLine();
        outMFile.write("------------------------------------");
        outMFile.newLine();
        outMFile.write("Query =" + qrefName);
        outMFile.newLine();
        outMFile.write("Target = " + trefName);
        outMFile.newLine();
        outMFile.write("Max atoms matched=\t" + nAtomsMatched);
        outMFile.newLine();
        for (Map.Entry<Integer, Integer> map : mcsNumber.entrySet()) {
            int queryIndex = map.getKey();
            int targetIndex = map.getValue();
            outMFile.write((queryIndex + 1) + "\t" + (targetIndex + 1));
            outMFile.newLine();
        }
    }

    void closeFiles() throws IOException {
        outGFile.close();
        outMFile.close();
        outDescriptorFile.close();
    }

    /**
     *
     * @param mol1
     * @param mol2
     * @param label
     * @param mcsNumber
     * @throws IOException
     * @throws CloneNotSupportedException
     */
    void makeImage(IAtomContainer mol1, IAtomContainer mol2, String label, Map<Integer, Integer> mcsNumber) throws IOException, CloneNotSupportedException {
        imageGenerator.addImages(mol1, mol2, label, mcsNumber);
    }

    /**
     *
     * @param mol1
     * @param mol2
     * @param label
     * @param mcsNumber
     * @throws CloneNotSupportedException
     * @throws IOException
     */
    void addImage(IAtomContainer mol1, IAtomContainer mol2, String label, Map<Integer, Integer> mcsNumber) throws CloneNotSupportedException, IOException {
        imageGenerator.addImages(mol1, mol2, label, mcsNumber);
    }

    /**
     *
     * @param qName
     * @param tName
     */
    void writeImage(String qName, String tName) {
        String suffix = argumentHandler.getSuffix();
        String outImageFileName = qName + "_" + tName + suffix;
        int w = argumentHandler.getImageWidth();
        int h = argumentHandler.getImageHeight();
        if (w != - 1 && h != -1) {
            imageGenerator.createImage(outImageFileName, qName, tName, w, h);
        } else {
            imageGenerator.createImage(outImageFileName, qName, tName);
        }

    }

    /**
     *
     * @param hub
     * @param rim
     * @param name
     * @param mappings
     * @throws IOException
     */
    void writeCircleImage(IAtomContainer hub, List<IAtomContainer> rim, String name, List<Map<Integer, Integer>> mappings) throws IOException {
        int w = argumentHandler.getImageWidth();
        int h = argumentHandler.getImageHeight();
        RenderedImage image;
        if (w != -1 && h != -1) {
            // these are the width and heights of the individual canvases
            image = imageGenerator.createHubWheelImage(hub, rim, mappings, w, h);
        } else {
            image = imageGenerator.createHubWheelImage(hub, rim, mappings);
        }
        ImageIO.write(image, "PNG", new File(name + ".png"));
    }

    /**
     *
     * @param mol
     * @param filepath
     * @throws IOException
     * @throws IllegalArgumentException
     * @throws CDKException
     */
    void writeMolToMolfile(IAtomContainer mol, String filepath) throws IOException, IllegalArgumentException, CDKException {
        Writer out;
        if (filepath.equals("--")) {
            out = new PrintWriter(System.out);
        } else {
            out = new FileWriter(filepath);
        }
        writeMolToMolfile(mol, out);
    }

    /**
     *
     * @param mol
     * @param out
     * @throws IOException
     * @throws IllegalArgumentException
     * @throws CDKException
     */
    void writeMolToMolfile(IAtomContainer mol, Writer out) throws IOException, IllegalArgumentException, CDKException {
        try (MDLV2000Writer writer = new MDLV2000Writer(out)) {
            writer.write(DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainer.class, mol));
            writer.close();
        }
    }

    /**
     *
     * @param mols
     * @param out
     * @throws IOException
     * @throws IllegalArgumentException
     * @throws CDKException
     */
    void writeMolsToMolfile(IAtomContainerSet mols, Writer out) throws IOException, IllegalArgumentException, CDKException {
        try (SDFWriter writer = new SDFWriter(out)) {
            writer.write(mols);
            writer.close();
        }
    }

    /**
     *
     * @param mol
     * @param out
     * @throws IOException
     */
    void writeMolToSmiles(IAtomContainer mol, Writer out) throws IOException, CDKException {
        SmilesGenerator smilesGenerator = new SmilesGenerator().aromatic();
        String smiles = smilesGenerator.create(mol);
        out.write(smiles);
        out.write(NEW_LINE);
        out.close();
    }

    /**
     *
     * @param mols
     * @param out
     * @throws IOException
     */
    void writeMolsToSmiles(IAtomContainerSet mols, Writer out) throws IOException, CDKException {
        SmilesGenerator smilesGenerator = new SmilesGenerator().aromatic();
        for (IAtomContainer mol : mols.atomContainers()) {
            String smiles = smilesGenerator.create(mol);
            out.write(smiles);
            out.write(NEW_LINE);
        }
        out.close();
    }

    void writeMol(String outputType, IAtomContainer mol, String filepath) throws IOException, IllegalArgumentException, CDKException {
        Writer out;
        if (filepath.equals("--")) {
            Writer outWriter = argumentHandler.getOutputWriter();
            if (outWriter == null) {
                out = new PrintWriter(System.out);
            } else {
                out = outWriter;
            }
        } else {
            out = new FileWriter(filepath);
        }
        switch (outputType) {
            case "MOL":
                writeMolToMolfile(mol, out);
                break;
            case "SMI":
                writeMolToSmiles(mol, out);
                break;
        }
    }

    void writeMol(String outputType, IAtomContainerSet mols, String filepath) throws IOException, IllegalArgumentException, CDKException {
        Writer out;
        if (filepath.equals("--")) {
            Writer outWriter = argumentHandler.getOutputWriter();
            if (outWriter == null) {
                out = new PrintWriter(System.out);
            } else {
                out = outWriter;
            }
        } else {
            out = new FileWriter(filepath);
        }
        if (outputType.equals("MOL")) {
            writeMolsToMolfile(mols, out);
        } else if (outputType.equals("SMI")) {
            writeMolsToSmiles(mols, out);
        }
    }
}
