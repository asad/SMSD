/* SPDX-License-Identifier: Apache-2.0 */
package com.bioinception.smsd.core;

import java.util.Arrays;
import java.util.Map;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.CsvSource;
import org.junit.jupiter.params.provider.ValueSource;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;
import static org.junit.jupiter.api.Assertions.*;

class FingerprintEdgeRegressionTest {
  private static IAtomContainer molecule(String smiles) throws Exception {
    return new SmilesParser(SilentChemObjectBuilder.getInstance()).parseSmiles(smiles);
  }

  private static MolGraph graph(String smiles) throws Exception {
    return FingerprintEngine.toMolGraph(molecule(smiles));
  }

  @Test void aromaticTorsionsDistinguishBenzeneFromCyclohexane() throws Exception {
    MolGraph benzene = graph("c1ccccc1"), cyclohexane = graph("C1CCCCC1");
    assertFalse(Arrays.equals(FingerprintEngine.topologicalTorsion(benzene, 2048),
        FingerprintEngine.topologicalTorsion(cyclohexane, 2048)));
    assertFalse(Arrays.equals(FingerprintEngine.topologicalTorsionCounts(benzene, 2048),
        FingerprintEngine.topologicalTorsionCounts(cyclohexane, 2048)));
    assertArrayEquals(new int[]{1337}, FingerprintEngine.toBitSet(FingerprintEngine.topologicalTorsion(benzene, 2048)).stream().toArray());
    assertArrayEquals(new int[]{1177}, FingerprintEngine.toBitSet(FingerprintEngine.topologicalTorsion(cyclohexane, 2048)).stream().toArray());
    assertArrayEquals(FingerprintEngine.topologicalTorsion(benzene, 2048),
        FingerprintEngine.topologicalTorsion(molecule("C1=CC=CC=C1"), 2048));
  }

  @Test void pathLimitsBeyondTheMoleculeDoNotAllocateUnreachableDepths() throws Exception {
    IAtomContainer ethanol = molecule("CCO");
    assertArrayEquals(FingerprintEngine.pathFingerprint(ethanol, 2, 1024),
        FingerprintEngine.pathFingerprint(ethanol, Integer.MAX_VALUE, 1024));
    MolGraph graph = graph("CCO");
    assertArrayEquals(FingerprintEngine.mcsFingerprint(graph, 2, 1024),
        FingerprintEngine.mcsFingerprint(graph, Integer.MAX_VALUE, 1024));
  }

  @ParameterizedTest @CsvSource({"C,0", "CCO,2", "CC(C)C,2", "C1CCCCC1,3", "CC.CCC,2", "C.CC,1"})
  void unlimitedRadiusStopsAtTheLargestFiniteGraphDistance(String smiles, int extent) throws Exception {
    MolGraph graph = graph(smiles);
    assertArrayEquals(FingerprintEngine.ecfp(graph, extent, 2048), FingerprintEngine.ecfp(graph, -1, 2048));
    assertArrayEquals(FingerprintEngine.fcfp(graph, extent, 2048), FingerprintEngine.fcfp(graph, -1, 2048));
    assertArrayEquals(FingerprintEngine.ecfpCounts(graph, extent, 2048), FingerprintEngine.ecfpCounts(graph, -1, 2048));
    assertArrayEquals(FingerprintEngine.fcfpCounts(graph, extent, 2048), FingerprintEngine.fcfpCounts(graph, -1, 2048));
    int expandedAtoms = 0;
    for (int[] neighbours : graph.neighbors) if (neighbours.length > 0) expandedAtoms++;
    assertEquals(graph.n + extent * expandedAtoms, Arrays.stream(FingerprintEngine.ecfpCounts(graph, -1, 1)).sum());
    assertEquals(graph.n + extent * expandedAtoms, Arrays.stream(FingerprintEngine.fcfpCounts(graph, -1, 1)).sum());
  }

  @Test void zeroLengthPathsContainOnlyTheSelectedAtomTypes() throws Exception {
    IAtomContainer ethanol = molecule("CCO"), atoms = molecule("C.C.O");
    assertArrayEquals(FingerprintEngine.pathFingerprint(atoms, 0, 1024),
        FingerprintEngine.pathFingerprint(ethanol, 0, 1024));
    assertFalse(Arrays.equals(FingerprintEngine.pathFingerprint(ethanol, 0, 1024),
        FingerprintEngine.pathFingerprint(ethanol, 1, 1024)));
    MolGraph graph = graph("CCO");
    long[] expected = FingerprintEngine.mcsFingerprint(graph, 0, 1024);
    int atomTypes = Arrays.stream(expected).mapToInt(Long::bitCount).sum();
    assertTrue(atomTypes > 0 && atomTypes <= graph.n);
    long[] withBonds = FingerprintEngine.mcsFingerprint(graph, 1, 1024);
    assertTrue(FingerprintEngine.fingerprintSubset(expected, withBonds));
    assertFalse(Arrays.equals(expected, withBonds));
    assertArrayEquals(expected, FingerprintEngine.mcsFingerprint(ethanol, new ChemOptions(), 0, 1024));
  }

  @ParameterizedTest @ValueSource(ints = {1, 8, 16, 65, 1000, 1024})
  void qualityUsesOnlyNominalBitsAndAcceptsEveryPositiveSize(int size) {
    long[] full = new long[(size + 63) / 64];
    Arrays.fill(full, -1L);
    Map<String, Double> quality = FingerprintEngine.analyzeFingerprintQuality(full, size);
    assertEquals(size, quality.get("bitsFilled"));
    assertEquals(1, quality.get("fillRate"));
    assertTrue(Double.isFinite(quality.get("chiSquared")));
    Map<String, Object> simple = FingerprintEngine.analyzeFpQuality(full, size);
    assertEquals(size, simple.get("set_bits"));
    assertEquals(1.0, simple.get("density"));
  }

  @Test void invalidDimensionsHaveExplicitArgumentErrors() throws Exception {
    IAtomContainer molecule = molecule("CC");
    MolGraph graph = graph("CC");
    assertThrows(IllegalArgumentException.class, () -> FingerprintEngine.mcsFingerprint(graph, -1, 1024));
    assertThrows(IllegalArgumentException.class, () -> FingerprintEngine.mcsFingerprint(graph, 1, 0));
    for (int size : new int[]{0, -1, 65}) {
      assertThrows(IllegalArgumentException.class, () -> FingerprintEngine.analyzeFingerprintQuality(new long[1], size));
      assertThrows(IllegalArgumentException.class, () -> FingerprintEngine.analyzeFpQuality(new long[1], size));
    }
    assertThrows(NullPointerException.class, () -> FingerprintEngine.analyzeFingerprintQuality(null, 1));
    assertThrows(NullPointerException.class, () -> FingerprintEngine.analyzeFpQuality(null, 1));
    assertThrows(IllegalArgumentException.class, () -> FingerprintEngine.pathFingerprint(molecule, -1, 1024));
  }
}
