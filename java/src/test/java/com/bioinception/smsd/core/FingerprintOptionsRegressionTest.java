/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
package com.bioinception.smsd.core;

import static com.bioinception.smsd.TestSupport.rawMolecule;
import static org.junit.jupiter.api.Assertions.*;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.List;
import org.junit.jupiter.api.Test;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IStereoElement;
import org.openscience.cdk.interfaces.ITetrahedralChirality;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.stereo.TetrahedralChirality;

/** Configurable fingerprints: classic configurations keep the published bits, switches extend them. */
class FingerprintOptionsRegressionTest {
  private static final String[] MOLECULES = {
    "C", "CC", "CCO", "CC(=O)O", "c1ccccc1", "c1ccncc1", "Oc1ccccc1", "Cc1ccccc1", "CC(=O)Oc1ccccc1C(=O)O",
    "CC(=O)Nc1ccc(O)cc1", "Cn1cnc2c1c(=O)n(C)c(=O)n2C", "CC(C)Cc1ccc(C(C)C(=O)O)cc1", "C[C@H](N)C(=O)O",
    "C[C@@H](N)C(=O)O", "C/C=C/C(=O)O", "C/C=C\\C(=O)O", "[13CH4]", "[NH4+]", "[O-]c1ccccc1", "[Na+].CC(=O)[O-]",
    "CCO.CCO", "c1ccccc1.CCO", "C1CC1", "C1CCCCCCC1", "C12C3C4C1C5C2C3C45", "C1C2CC3CC1CC(C2)C3",
    "c1ccc2ccccc2c1", "C1CCC2(C1)CCCCC2", "O=C1CCCCC1", "c1ccccc1CCc1ccccc1", "NS(=O)(=O)c1ccc(Cl)cc1",
    "OC(=O)CCC(=O)O", "CC(=O)CC(C)=O", "[H]OC", "[H]N([H])[H]", "BrC(Br)Br", "FC(F)(F)c1ccccc1", "CC#N", "C=N",
  };

  private static IAtomContainer container(String smiles) throws Exception {
    return smiles.isEmpty() ? SilentChemObjectBuilder.getInstance().newAtomContainer() : rawMolecule(smiles);
  }

  private static MolGraph graph(String smiles) throws Exception {
    return FingerprintEngine.toMolGraph(container(smiles));
  }

  private static long[] bits(String smiles, FingerprintOptions options) throws Exception {
    return FingerprintEngine.fingerprint(graph(smiles), options);
  }

  /** The test molecules and, when the benchmark corpus is present, its first molecules. */
  private static List<String> corpus() throws Exception {
    List<String> smiles = new ArrayList<>(Arrays.asList(MOLECULES));
    Path file = Paths.get("benchmarks/data/chembl_mcs_benchmark.smi");
    if (Files.isRegularFile(file)) {
      int taken = 0;
      for (String line : Files.readAllLines(file, StandardCharsets.UTF_8)) {
        if (line.trim().isEmpty() || taken >= 1500) continue;
        smiles.add(line.trim().split("\\s+")[0]);
        taken++;
      }
    }
    return smiles;
  }

  @Test
  void classicConfigurationsEqualTheFixedMethods() throws Exception {
    int checked = 0;
    for (String smiles : corpus()) {
      IAtomContainer container;
      try { container = container(smiles); } catch (Exception unparsable) { continue; }
      MolGraph graph = FingerprintEngine.toMolGraph(container);
      for (int radius : new int[] {0, 2, -1})
        for (int size : new int[] {64, 1000, 2048}) {
          assertArrayEquals(FingerprintEngine.ecfp(graph, radius, size),
              FingerprintEngine.fingerprint(graph, FingerprintOptions.ecfp(radius, size)), smiles + " ecfp");
          assertArrayEquals(FingerprintEngine.ecfpCounts(graph, radius, size),
              FingerprintEngine.fingerprintCounts(graph, FingerprintOptions.ecfp(radius, size)), smiles + " ecfp counts");
          assertArrayEquals(FingerprintEngine.fcfp(graph, radius, size),
              FingerprintEngine.fingerprint(graph, FingerprintOptions.fcfp(radius, size)), smiles + " fcfp");
          assertArrayEquals(FingerprintEngine.fcfpCounts(graph, radius, size),
              FingerprintEngine.fingerprintCounts(graph, FingerprintOptions.fcfp(radius, size)), smiles + " fcfp counts");
        }
      for (int length : new int[] {1, 3, 7})
        assertArrayEquals(FingerprintEngine.pathFingerprint(container, length, 2048),
            FingerprintEngine.fingerprint(graph, FingerprintOptions.path(length, 2048)), smiles + " path");
      for (int size : new int[] {64, 1000, 2048}) {
        assertArrayEquals(FingerprintEngine.topologicalTorsion(graph, size),
            FingerprintEngine.fingerprint(graph, FingerprintOptions.torsion(size)), smiles + " torsion");
        assertArrayEquals(FingerprintEngine.topologicalTorsionCounts(graph, size),
            FingerprintEngine.fingerprintCounts(graph, FingerprintOptions.torsion(size)), smiles + " torsion counts");
      }
      checked++;
    }
    assertTrue(checked >= MOLECULES.length, "the test molecules were checked");
    assertArrayEquals(FingerprintEngine.ecfp(graph("CCO"), 2, 2048), bits("CCO", new FingerprintOptions()),
        "The default options are the classic ECFP4 configuration");
  }

  @Test
  void containerOverloadsStandardiseLikeTheGraphForms() throws Exception {
    FingerprintOptions options = FingerprintOptions.ecfp(2, 512);
    options.cipAtoms = options.ringSize = true;
    for (String smiles : MOLECULES) {
      IAtomContainer container = container(smiles);
      assertArrayEquals(FingerprintEngine.fingerprint(FingerprintEngine.toMolGraph(container), options),
          FingerprintEngine.fingerprint(container, options), smiles);
      assertArrayEquals(FingerprintEngine.fingerprintCounts(FingerprintEngine.toMolGraph(container), options),
          FingerprintEngine.fingerprintCounts(container, options), smiles);
    }
  }

  @Test
  void stereoTermsSeparateIsomers() throws Exception {
    String left = "C[C@H](N)C(=O)O", right = "C[C@@H](N)C(=O)O";
    for (FingerprintOptions.Family family : new FingerprintOptions.Family[] {
        FingerprintOptions.Family.ECFP, FingerprintOptions.Family.FCFP, FingerprintOptions.Family.PATH,
        FingerprintOptions.Family.TORSION}) {
      FingerprintOptions plain = FingerprintOptions.of(family), chiral = FingerprintOptions.of(family);
      chiral.cipAtoms = true;
      assertArrayEquals(bits(left, plain), bits(right, plain), family + " without the CIP term");
      assertFalse(Arrays.equals(bits(left, chiral), bits(right, chiral)), family + " with the CIP term");
    }
    FingerprintOptions radiusZero = FingerprintOptions.ecfp(0, 2048);
    radiusZero.cipAtoms = true;
    assertArrayEquals(bits(left, radiusZero), bits(right, radiusZero), "radius 0 does not separate enantiomers");
    FingerprintOptions chiral = FingerprintOptions.ecfp(2, 4096);
    chiral.cipAtoms = true;
    BitSet differing = BitSet.valueOf(bits(left, chiral));
    differing.xor(BitSet.valueOf(bits(right, chiral)));
    assertTrue(differing.cardinality() <= 10, "enantiomers differ only where the centre lies");
    FingerprintOptions geometric = FingerprintOptions.ecfp();
    geometric.cipBonds = true;
    assertFalse(Arrays.equals(bits("C/C=C/C(=O)O", geometric), bits("C/C=C\\C(=O)O", geometric)));
    assertArrayEquals(bits("C/C=C/C(=O)O", FingerprintOptions.ecfp()), bits("C/C=C\\C(=O)O", FingerprintOptions.ecfp()));
  }

  @Test
  void ringSizeFragmentsAndScaffold() throws Exception {
    FingerprintOptions ringSize = FingerprintOptions.ecfp();
    ringSize.ringSize = true;
    assertFalse(Arrays.equals(bits("C1CCCCCC1", ringSize), bits("C1CCCCCCC1", ringSize)), "sizes below the cap differ");
    assertArrayEquals(bits("C1CCCCCCC1", ringSize), bits("C1CCCCCCCC1", ringSize), "sizes at the cap or more share a class");
    ringSize.ringSizeCap = 12;
    assertFalse(Arrays.equals(bits("C1CCCCCCC1", ringSize), bits("C1CCCCCCCC1", ringSize)), "a higher cap separates them");

    FingerprintOptions largest = FingerprintOptions.ecfp();
    largest.fragments = FingerprintOptions.Fragments.LARGEST;
    assertArrayEquals(bits("CC(=O)O", FingerprintOptions.ecfp()), bits("CCO.CC(=O)O", largest));
    assertArrayEquals(bits("CC(=O)[O-]", FingerprintOptions.ecfp()), bits("[Na+].CC(=O)[O-]", largest));
    FingerprintOptions perFragment = FingerprintOptions.ecfp();
    perFragment.fragments = FingerprintOptions.Fragments.PER_FRAGMENT;
    assertFalse(Arrays.equals(bits("c1ccccc1.CCO", perFragment), bits("c1ccccc1.CCO", FingerprintOptions.ecfp())));
    assertArrayEquals(bits("CCO", perFragment), bits("CCO.CCO", perFragment), "identical components share one class");

    FingerprintOptions scaffold = FingerprintOptions.ecfp();
    scaffold.scaffold = FingerprintOptions.Scaffold.MURCKO;
    long[] benzene = bits("c1ccccc1", FingerprintOptions.ecfp());
    assertArrayEquals(benzene, bits("Cc1ccccc1", scaffold), "toluene reduces to benzene");
    assertArrayEquals(benzene, bits("CC(=O)Nc1ccc(O)cc1", scaffold), "side chains are removed");
    assertArrayEquals(bits("O=C1CCCCC1", FingerprintOptions.ecfp()), bits("O=C1CCCCC1", scaffold), "a carbonyl stays");
    assertArrayEquals(bits("O=c1cccc[nH]1", FingerprintOptions.ecfp()), bits("Cn1ccccc1=O", scaffold),
        "an N-methylpyridone reduces to pyridone");
    assertArrayEquals(benzene, bits("CC(=O)c1ccccc1", scaffold), "a side chain is cut at its single bond");
    assertArrayEquals(bits("[nH]1cnc2c1c(=O)[nH]c(=O)[nH]2", FingerprintOptions.ecfp()),
        bits("Cn1cnc2c1c(=O)n(C)c(=O)n2C", scaffold), "the framework of caffeine is xanthine");
    assertArrayEquals(bits("c1ccccc1CCc1ccccc1", FingerprintOptions.ecfp()), bits("c1ccccc1CCc1ccccc1", scaffold), "linkers stay");
    assertArrayEquals(bits("CCO", FingerprintOptions.ecfp()), bits("CCO", scaffold), "an acyclic molecule is unchanged");

    FingerprintOptions chiral = FingerprintOptions.ecfp(), framedChiral = FingerprintOptions.ecfp();
    chiral.cipAtoms = framedChiral.cipAtoms = true;
    framedChiral.scaffold = FingerprintOptions.Scaffold.MURCKO;
    assertArrayEquals(bits("[C@H]12CCCC[C@@H]1CCCO2", chiral), bits("C[C@]12CCCC[C@@H]1CCCO2", framedChiral),
        "a ring-fusion centre keeps its configuration when its angular methyl is capped");
    assertArrayEquals(bits("[C@@H]12CCCC[C@@H]1CCCO2", chiral), bits("C[C@@]12CCCC[C@@H]1CCCO2", framedChiral));
    FingerprintOptions geometric = FingerprintOptions.ecfp(), framedGeometric = FingerprintOptions.ecfp();
    geometric.cipBonds = framedGeometric.cipBonds = true;
    framedGeometric.scaffold = FingerprintOptions.Scaffold.MURCKO;
    assertArrayEquals(bits("c1ccccc1/C=C\\c1ccccc1", geometric), bits("C/C(c1ccccc1)=C/c1ccccc1", framedGeometric),
        "a double bond keeps its configuration when its reference substituent is capped");
    assertArrayEquals(bits("C1CC=CCC1", geometric), bits("C1C/C=C\\CC1", geometric), "no E/Z in a ring below eight atoms");
    assertFalse(Arrays.equals(bits("C1CCC/C=C/CC1", geometric), bits("C1CCC/C=C\\CC1", geometric)), "cyclooctene keeps E/Z");
  }

  @Test
  void substructureKeys() throws Exception {
    assertEquals(166, FingerprintOptions.maccsKeys().size());
    FingerprintOptions maccs = FingerprintOptions.maccs();
    assertEquals(FingerprintOptions.Family.KEYS, maccs.family);
    assertEquals(166, maccs.size);
    BitSet benzene = BitSet.valueOf(bits("c1ccccc1", maccs));
    assertTrue(benzene.get(162) && benzene.get(164) && benzene.get(161), "six-membered ring, ring and aromatic keys");
    assertFalse(benzene.get(163) || benzene.get(160), "no oxygen or nitrogen key");
    BitSet ethanol = BitSet.valueOf(bits("CCO", maccs));
    assertTrue(ethanol.get(163) && ethanol.get(156) && ethanol.get(138) && !ethanol.get(158));
    assertTrue(BitSet.valueOf(bits("CC(=O)O.O", maccs)).get(165), "two fragments set key 166");
    assertTrue(BitSet.valueOf(bits("[13CH4]", maccs)).get(0), "an isotope sets key 1");
    assertTrue(BitSet.valueOf(bits("c1ccc2ccccc2c1", maccs)).get(124), "two aromatic rings set key 125");
    assertEquals(0, BitSet.valueOf(bits("", maccs)).cardinality());

    FingerprintOptions custom = FingerprintOptions.withKeys(Arrays.asList(new FingerprintOptions.Key("c1ccccc1"),
        new FingerprintOptions.Key("[#8]", 2), new FingerprintOptions.Key("[#8]", 5)));
    MolGraph aspirin = graph("CC(=O)Oc1ccccc1C(=O)O");
    assertEquals(BitSet.valueOf(new long[] {3L}), BitSet.valueOf(FingerprintEngine.fingerprint(aspirin, custom)));
    assertArrayEquals(new int[] {1, 4, 4}, FingerprintEngine.fingerprintCounts(aspirin, custom));
  }

  @Test
  void bitsDoNotDependOnTheAtomOrder() throws Exception {
    FingerprintOptions chiral = FingerprintOptions.ecfp(), geometric = FingerprintOptions.ecfp();
    FingerprintOptions perFragment = FingerprintOptions.ecfp(), largest = FingerprintOptions.ecfp();
    FingerprintOptions framework = FingerprintOptions.ecfp(), ringSize = FingerprintOptions.ecfp();
    chiral.cipAtoms = true;
    geometric.cipBonds = true;
    perFragment.fragments = FingerprintOptions.Fragments.PER_FRAGMENT;
    largest.fragments = FingerprintOptions.Fragments.LARGEST;
    framework.scaffold = FingerprintOptions.Scaffold.MURCKO;
    ringSize.ringSize = true;
    Object[][] cases = {
      {"C[C@H](N)C(=O)O", "N[C@@H](C)C(=O)O", chiral},
      {"C/C=C/C(=O)O", "OC(=O)/C=C/C", geometric},
      {"c1ccccc1.CCO", "OCC.c1ccccc1", perFragment},
      {"CCO.CCN", "CCN.CCO", largest},
      {"Cc1ccccc1CC(C)c1ccccc1O", "Oc1ccccc1C(C)Cc1ccccc1C", framework},
      {"C1CC1c1ccccc1", "c1ccc(cc1)C1CC1", ringSize},
    };
    for (Object[] item : cases)
      assertArrayEquals(bits((String) item[0], (FingerprintOptions) item[2]), bits((String) item[1], (FingerprintOptions) item[2]),
          (String) item[0]);
    assertFalse(Arrays.equals(bits("C[C@@H](N)C(=O)O", chiral), bits("N[C@@H](C)C(=O)O", chiral)));
  }

  /** A CDK stereo element may list its ligands in any order; the reduced graphs must read the tag the same way. */
  @Test
  void containerTagsWithAnotherLigandOrderSurviveReduction() throws Exception {
    FingerprintOptions framework = FingerprintOptions.ecfp(), largest = FingerprintOptions.ecfp();
    framework.scaffold = FingerprintOptions.Scaffold.MURCKO;
    framework.cipAtoms = true;
    largest.fragments = FingerprintOptions.Fragments.LARGEST;
    largest.cipAtoms = true;
    String left = "c1ccccc1[C@H](c1ccncc1)c1ccco1", right = "c1ccccc1[C@@H](c1ccncc1)c1ccco1";
    for (String smiles : new String[] {left, right}) {
      assertArrayEquals(FingerprintEngine.fingerprint(rawMolecule(smiles), framework),
          FingerprintEngine.fingerprint(withSwappedLigands(rawMolecule(smiles)), framework), smiles);
      assertArrayEquals(FingerprintEngine.fingerprint(rawMolecule(smiles + ".O"), largest),
          FingerprintEngine.fingerprint(withSwappedLigands(rawMolecule(smiles + ".O")), largest), smiles);
    }
    assertFalse(Arrays.equals(FingerprintEngine.fingerprint(rawMolecule(left), framework),
        FingerprintEngine.fingerprint(rawMolecule(right), framework)), "the enantiomers differ after reduction");
  }

  /** The same configuration, with the first two ligands of each tetrahedral element swapped and the winding reversed. */
  private static IAtomContainer withSwappedLigands(IAtomContainer molecule) {
    List<IStereoElement> elements = new ArrayList<>();
    for (IStereoElement element : molecule.stereoElements()) {
      if (element instanceof ITetrahedralChirality) {
        ITetrahedralChirality tetrahedral = (ITetrahedralChirality) element;
        IAtom[] ligands = tetrahedral.getLigands().clone();
        IAtom swap = ligands[0];
        ligands[0] = ligands[1];
        ligands[1] = swap;
        ITetrahedralChirality.Stereo winding = tetrahedral.getStereo() == ITetrahedralChirality.Stereo.CLOCKWISE
            ? ITetrahedralChirality.Stereo.ANTI_CLOCKWISE : ITetrahedralChirality.Stereo.CLOCKWISE;
        elements.add(new TetrahedralChirality(tetrahedral.getChiralAtom(), ligands, winding));
      } else {
        elements.add(element);
      }
    }
    molecule.setStereoElements(elements);
    return molecule;
  }

  @Test
  void optionsAreValidatedAndCopied() throws Exception {
    MolGraph ethanol = graph("CCO");
    assertThrows(IllegalArgumentException.class, () -> FingerprintEngine.fingerprint(ethanol, FingerprintOptions.ecfp(2, 0)));
    assertThrows(IllegalArgumentException.class, () -> FingerprintEngine.fingerprint(ethanol, FingerprintOptions.ecfp(-2, 64)));
    assertThrows(IllegalArgumentException.class, () -> FingerprintEngine.fingerprint(ethanol, FingerprintOptions.path(0, 64)));
    FingerprintOptions cap = FingerprintOptions.ecfp();
    cap.ringSizeCap = 2;
    assertThrows(IllegalArgumentException.class, () -> FingerprintEngine.fingerprint(ethanol, cap));
    assertThrows(IllegalArgumentException.class,
        () -> FingerprintEngine.fingerprint(ethanol, FingerprintOptions.of(FingerprintOptions.Family.KEYS)));
    assertThrows(IllegalArgumentException.class, () -> FingerprintEngine.fingerprint(ethanol,
        FingerprintOptions.withKeys(Arrays.asList(new FingerprintOptions.Key("C", 0)))));
    assertThrows(IllegalArgumentException.class, () -> FingerprintEngine.fingerprint(ethanol,
        FingerprintOptions.withKeys(Arrays.asList(new FingerprintOptions.Key("C(")))));
    assertThrows(NullPointerException.class, () -> FingerprintEngine.fingerprint((MolGraph) null, FingerprintOptions.ecfp()));
    FingerprintOptions original = FingerprintOptions.maccs();
    original.cipAtoms = true;
    FingerprintOptions copy = original.copy();
    copy.cipAtoms = false;
    copy.keys.clear();
    assertTrue(original.cipAtoms);
    assertEquals(166, original.keys.size());
  }
}
