/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
package com.bioinception.smsd.cli;

import static org.junit.jupiter.api.Assertions.assertEquals;

import com.bioinception.smsd.cli.SMSDcli.MolIO;
import com.bioinception.smsd.core.MolGraph;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.CsvSource;
import org.openscience.cdk.interfaces.IAtomContainer;

/** MOL and SDF valence and radical information supplies missing hydrogen counts. */
class MolReaderHydrogenRegressionTest {
  @TempDir Path temporary;

  private static String methylcarbene(int valence, int radical) {
    String atom = "    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0";
    atom = atom.substring(0, 48) + String.format("%3d", valence) + atom.substring(51);
    return "Methylcarbene\n\n\n  2  1  0  0  0  0            999 V2000\n" + atom + "\n"
        + "    1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
        + "  1  2  1  0\n" + (radical == 0 ? "" : "M  RAD  1   1   " + radical + "\n")
        + "M  END\n";
  }

  @ParameterizedTest
  @CsvSource({"2,3,1", "2,2,1", "2,0,1", "1,3,0", "15,3,0", "0,3,1"})
  void moleculeAndSdfReadersRetainRadicalValence(int valence, int radical, int hydrogens) throws Exception {
    Path path = Files.write(temporary.resolve("molecule.mol"),
        methylcarbene(valence, radical).getBytes(StandardCharsets.UTF_8));
    IAtomContainer molecule = MolIO.loadTarget("MOL", path.toString());
    IAtomContainer query = MolIO.loadQuery("MOL", path.toString()).container();
    Path sdf = Files.write(temporary.resolve("molecule.sdf"),
        (methylcarbene(valence, radical) + "$$$$\n").getBytes(StandardCharsets.UTF_8));
    IAtomContainer record = MolIO.loadTargetsSDF(sdf.toString()).get(0);
    for (IAtomContainer input : new IAtomContainer[] {molecule, query, record}) {
      assertEquals(hydrogens, input.getAtom(0).getImplicitHydrogenCount());
      assertEquals(hydrogens, new MolGraph(input).hydrogenCount(0));
      assertEquals(radical == 3 ? 2 : radical == 2 ? 1 : 0,
          input.getConnectedSingleElectronsCount(input.getAtom(0)));
    }
  }

  @Test
  void v3000SdfRetainsRadicalValence() throws Exception {
    String block = "Methylcarbene\n\n\n  0  0  0  0  0  0            999 V3000\n"
        + "M  V30 BEGIN CTAB\nM  V30 COUNTS 2 1 0 0 0\nM  V30 BEGIN ATOM\n"
        + "M  V30 1 C 0 0 0 0 VAL=2 RAD=3\nM  V30 2 C 1.5 0 0 0\nM  V30 END ATOM\n"
        + "M  V30 BEGIN BOND\nM  V30 1 1 1 2\nM  V30 END BOND\nM  V30 END CTAB\nM  END\n$$$$\n";
    Path path = Files.write(temporary.resolve("molecule.sdf"), block.getBytes(StandardCharsets.UTF_8));
    IAtomContainer molecule = MolIO.loadTargetsSDF(path.toString()).get(0);
    assertEquals(1, molecule.getAtom(0).getImplicitHydrogenCount());
    assertEquals(1, new MolGraph(molecule).hydrogenCount(0));
    assertEquals(2, molecule.getConnectedSingleElectronsCount(molecule.getAtom(0)));
  }
}
