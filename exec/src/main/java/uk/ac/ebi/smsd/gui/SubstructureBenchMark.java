/* Copyright (C) 2009-2018  Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>
 *
 * Contact: Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>
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
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA.
 */
package uk.ac.ebi.smsd.gui;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.io.iterator.IIteratingChemObjectReader;
import org.openscience.cdk.io.iterator.IteratingSDFReader;
import org.openscience.cdk.isomorphism.UniversalIsomorphismTester;
import org.openscience.smsd.Substructure;

/**
 *
 *  java1.8+
 *
 * 
 * 
 *
 * @author Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>
 *
 */
public class SubstructureBenchMark {

    /**
     * @param args the command line arguments
     * @throws FileNotFoundException
     * @throws IOException
     * @throws CDKException
     */
    public static void main(String[] args) throws FileNotFoundException, IOException, CDKException {
        String queryFilePath = (args.length > 0) ? args[0] : "data/actives.sdf";
        String targetFilePath = (args.length > 1) ? args[1] : "data/all.sdf";
        File qFile = new File(queryFilePath);
        File tFile = new File(targetFilePath);

        IIteratingChemObjectReader qFileReader = read(qFile);

        if (qFileReader == null) {
            throw new IOException("Unknown input type ");
        }
        List<IAtomContainer> targets = new ArrayList<>();
        IIteratingChemObjectReader tFileReader = read(tFile);
        while (tFileReader.hasNext()) {
            targets.add((IAtomContainer) tFileReader.next());
        }

        int counter = 0;
        while (qFileReader.hasNext()) {
            IAtomContainer query = (IAtomContainer) qFileReader.next();
            int smsdSolutionCount = 0;
            int uitSolutionCount = 0;
            long t0 = System.currentTimeMillis();
            for (IAtomContainer target : targets) {
                smsdSolutionCount += getSMSDSolutionCount(query, target);
            }
            long timeNow = System.currentTimeMillis();
            long smsdTime = (timeNow - t0);

//            tFileReader = read(tFile);
            long tUIT0 = System.currentTimeMillis();
            for (IAtomContainer target : targets) {
                uitSolutionCount += getUITSolutionCount(query, target);
            }
            timeNow = System.currentTimeMillis();
            long uitTime = timeNow - tUIT0;

            String out = String.format("%d SMSDt %d SMSDs %d UITt %d UITs %d ",
                    counter, smsdTime, smsdSolutionCount, uitTime, uitSolutionCount);
            System.out.println(out);
            counter++;
        }

    }

    private static int getSMSDSolutionCount(IAtomContainer queryMol, IAtomContainer target) throws CDKException {

        Substructure substructure = new Substructure(queryMol, target, true, true, true, false);

        if (substructure.isSubgraph()) {
            return 1;
        } else {
            return 0;
        }
    }

    private static int getUITSolutionCount(IAtomContainer query, IAtomContainer target) throws CDKException {
        UniversalIsomorphismTester uti = new UniversalIsomorphismTester();
        if (uti.isSubgraph(target, query)) {
            return 1;
        } else {
            return 0;
        }
    }

    /**
     *
     * @param file
     * @return IteratingChemObjectReader
     * @throws FileNotFoundException
     */
    public static IIteratingChemObjectReader read(File file) throws FileNotFoundException {

        FileReader in = new FileReader(file);
        return new IteratingSDFReader(
                in, DefaultChemObjectBuilder.getInstance());

    }
}
