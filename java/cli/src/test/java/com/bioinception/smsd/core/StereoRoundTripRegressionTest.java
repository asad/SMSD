/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
package com.bioinception.smsd.core;

import static org.junit.jupiter.api.Assertions.*;

import com.bioinception.smsd.cli.SMSDcli.OutputUtil;
import com.bioinception.smsd.cli.SMSDcli.OutputUtil.OutType;
import java.io.ByteArrayOutputStream;
import java.io.StringReader;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import javax.vecmath.Point2d;
import javax.vecmath.Point3d;
import org.openscience.cdk.interfaces.IDoubleBondStereochemistry;
import org.openscience.cdk.io.MDLV2000Reader;
import org.openscience.cdk.silent.AtomContainer;
import java.util.Map;
import java.util.TreeMap;
import org.junit.jupiter.api.Test;
import com.bioinception.smsd.cli.SMSDcli.MolIO;
import org.openscience.cdk.interfaces.IAtom;

import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.ValueSource;
import org.junit.jupiter.params.provider.CsvSource;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IStereoElement;
import org.openscience.cdk.isomorphism.Pattern;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;

class StereoRoundTripRegressionTest {
  @ParameterizedTest
  @ValueSource(strings = {
      "C[C@H](N)C(=O)O", "C[C@@H](N)C(=O)O", "CC[C@@H](O)C",
      "[C@H](F)(Cl)Br", "F[C@](Cl)(Br)I", "N[C@@H]1CCC[C@H]1O",
      "C/C=C/C", "C/C=C\\C", "C(=C/Cl)\\Cl", "Cl/C(C)=C/C",
      "C/C=C/C=C/C", "C/C=C\\C=C/C", "C/C=C/C=C\\C", "C/C=C\\C=C\\C",
      "C/C=C(/F)\\C(Cl)=C/C", "C/C=C(\\F)/C(Cl)=C/C"
  })
  void canonicalSmilesPreservesStereo(String smiles) throws Exception {
    SmilesParser parser = new SmilesParser(SilentChemObjectBuilder.getInstance());
    IAtomContainer original = parser.parseSmiles(smiles);
    assertRoundTrip(parser, original);

    IAtomContainer reordered = original.clone();
    IBond[] bonds = new IBond[reordered.getBondCount()];
    for (int bondIndex = 0; bondIndex < bonds.length; bondIndex++)
      bonds[bondIndex] = reordered.getBond(bonds.length - bondIndex - 1);
    reordered.setBonds(bonds);
    assertRoundTrip(parser, reordered);
  }

  @ParameterizedTest
  @ValueSource(strings = {"F/C=C1\\C=NC2=CNN=NC2=N1", "F/C=C1/C=NC2=CNN=NC2=N1",
      "C1=CC2=N/C1=C\\C1=N/C(=C\\C3=NC4=C(C5=N/C(=C\\2)CC5)CC(=O)C4=C3)C=C1",
      "C1=CC2=N/C1=C/C1=N/C(=C\\C3=NC4=C(C5=N/C(=C\\2)CC5)CC(=O)C4=C3)C=C1"})
  void canonicalSmilesPreservesExternalStereoAtAromaticCarriers(String smiles) throws Exception {
    SmilesParser parser = new SmilesParser(SilentChemObjectBuilder.getInstance());
    IAtomContainer original = Standardiser.standardise(
        parser.parseSmiles(smiles), Standardiser.TautomerMode.NONE);
    int stereoElements = smiles.startsWith("F/") ? 1 : 5;
    assertEquals(stereoElements, stereoCount(original));
    MolGraph graph = new MolGraph(original);
    long originalHash = graph.getCanonicalHash();
    String expected = assertAromaticStereoRoundTrip(graph, original);
    Map<Integer, Integer> identity = new TreeMap<>();
    for (int atom = 0; atom < graph.atomCount(); atom++) identity.put(atom, atom);
    MolGraph view = MCS.fragment(graph, identity, MatchGraph.allBonds(graph, identity).keySet()).toGraph(graph);
    assertEquals(expected, assertAromaticStereoRoundTrip(view, original));
    IAtomContainer reordered = original.clone();
    IBond[] bonds = new IBond[reordered.getBondCount()];
    for (int bondIndex = 0; bondIndex < bonds.length; bondIndex++)
      bonds[bondIndex] = reordered.getBond(bonds.length - bondIndex - 1);
    reordered.setBonds(bonds);
    assertEquals(expected, assertAromaticStereoRoundTrip(new MolGraph(reordered), original));
    IAtom[] atoms = new IAtom[reordered.getAtomCount()];
    for (int atom = 0; atom < atoms.length; atom++)
      atoms[atom] = reordered.getAtom(atoms.length - atom - 1);
    reordered.setAtoms(atoms);
    assertEquals(expected, assertAromaticStereoRoundTrip(new MolGraph(reordered), original));
    assertEquals(originalHash, graph.getCanonicalHash());
    assertEquals(stereoElements, stereoCount(original));
  }

  @Test
  void rawSmilesReaderRetainsAromaticStereoAndOriginalIndices() throws Exception {
    String first = "[Na+].F/C=c1\\cnc2c(n1)N=NNC=2.[13CH3][C@H](F)Cl.F/C=C/F";
    String second = "[Na+].F/C=c1/cnc2c(n1)N=NNC=2.[13CH3][C@H](F)Cl.F/C=C/F";
    SmilesParser parser = new SmilesParser(SilentChemObjectBuilder.getInstance());
    IAtomContainer raw = parser.parseSmiles(first);
    IAtomContainer loaded = MolIO.loadTarget("SMI", first);
    assertEquals(3, stereoCount(loaded));
    assertEquals(raw.getAtomCount(), loaded.getAtomCount());
    assertEquals(raw.getBondCount(), loaded.getBondCount());
    for (int atom = 0; atom < raw.getAtomCount(); atom++) {
      assertEquals(raw.getAtom(atom).getAtomicNumber(), loaded.getAtom(atom).getAtomicNumber());
      assertEquals(raw.getAtom(atom).getFormalCharge(), loaded.getAtom(atom).getFormalCharge());
      assertEquals(raw.getAtom(atom).getMassNumber(), loaded.getAtom(atom).getMassNumber());
      assertEquals(raw.getAtom(atom).getImplicitHydrogenCount(), loaded.getAtom(atom).getImplicitHydrogenCount());
    }
    for (int bond = 0; bond < raw.getBondCount(); bond++) {
      assertEquals(raw.indexOf(raw.getBond(bond).getBegin()), loaded.indexOf(loaded.getBond(bond).getBegin()));
      assertEquals(raw.indexOf(raw.getBond(bond).getEnd()), loaded.indexOf(loaded.getBond(bond).getEnd()));
      assertEquals(raw.getBond(bond).getOrder(), loaded.getBond(bond).getOrder());
    }
    IAtomContainer opposite = MolIO.loadQuery("SMI", second).container();
    assertEquals(3, stereoCount(opposite));
    assertFalse(Pattern.findIdentical(loaded).matches(opposite));
    Substructure.Options options = new Substructure.Options();
    options.chemistry.useChirality = true;
    options.chemistry.useBondStereo = true;
    Substructure prepared = Substructure.find(loaded).withOptions(options);
    assertTrue(prepared.match(loaded).found());
    assertFalse(prepared.match(opposite).found());
    assertEquals(SearchResult.Status.TimedOut, prepared.withTimeout(0).match(loaded).status());
    options.cancellation = () -> true;
    assertEquals(SearchResult.Status.Cancelled,
        Substructure.find(loaded).withOptions(options).match(loaded).status());
    ExecutorService executor = Executors.newFixedThreadPool(2);
    try {
      Future<IAtomContainer> firstRead =
          executor.submit(() -> MolIO.loadTarget("SMI", first));
      Future<IAtomContainer> secondRead =
          executor.submit(() -> MolIO.loadTarget("SMI", second));
      assertTrue(Pattern.findIdentical(loaded).matches(firstRead.get()));
      assertTrue(Pattern.findIdentical(opposite).matches(secondRead.get()));
    } finally {
      executor.shutdownNow();
    }
  }

  @Test
  void aromaticWriterDoesNotAssignUnspecifiedDoubleBondStereo() throws Exception {
    String smiles = "C1=CC2=N/C1=C\\C1=N/C(=C\\C3=NC4=C(C5=N/C(=C\\2)CC5)CC(=O)C4=C3)C=C1";
    IAtomContainer molecule = Standardiser.standardise(
        MolIO.loadTarget("SMI", smiles), Standardiser.TautomerMode.NONE);
    List<IStereoElement> retained = new ArrayList<>();
    int doubleBond = 0;
    for (IStereoElement<?, ?> stereo : molecule.stereoElements()) {
      if (!(stereo instanceof IDoubleBondStereochemistry)
          || (doubleBond++ & 1) == 0) retained.add(stereo);
    }
    molecule.setStereoElements(retained);
    assertEquals(3, stereoCount(molecule));
    assertAromaticStereoRoundTrip(new MolGraph(molecule), molecule);
  }

  @ParameterizedTest
  @CsvSource(delimiter = '|', value = {
      "CC[C@H]1/C2=C/C3=NC4=C(C5=N/C(=C\\C6=N/C(=C\\C(=N2)[C@@H]1C)C(C(C)=O)=C6C)[C@@H](C)[C@@H]5CCC(=O)OC/C=C(\\C)CCC[C@@H](C)CCC[C@@H](C)CCCC(C)C)[C@@H](C(=O)OC)C(=O)C4=C3C|12",
      "C=CC1=C(C)C2=N/C1=C\\C1=NC(=C\\C3=NC4=C(C5=N/C(=C\\2)[C@@H](C)[C@@H]5CCC(=O)OC/C=C(\\C)CCC[C@H](C)CCC[C@H](C)CCCC(C)C)[C@@H](C(=O)OC)C(=O)C4=C3C)/C(=C\\C)[C@H]1C|12",
      "CC1(CC(=O)O)C2=CC3=N/C(=C\\C4=NC(=C(CCC(=O)O)C4CC(=O)O)/C=C4\\N=C(/C=C(/N2)C1CCC(=O)O)C(CC(=O)O)=C4CCC(=O)O)C(C)(CC(N)=O)C3CCC(=O)O|4"})
  void canonicalSmilesPreservesExplicitSmallRingImineStereo(
      String smiles, int stereoElements) throws Exception {
    IAtomContainer molecule = Standardiser.standardise(
        MolIO.loadTarget("SMI", smiles), Standardiser.TautomerMode.NONE);
    assertEquals(stereoElements, stereoCount(molecule));
    assertAromaticStereoRoundTrip(new MolGraph(molecule), molecule);
  }

  @Test
  void aromaticStereoReaderRejectsInvalidInputAndCanBeReused() throws Exception {
    assertThrows(org.openscience.cdk.exception.CDKException.class,
        () -> MolIO.loadTarget("SMI", "c1/c=c/cc1"));
    assertThrows(org.openscience.cdk.exception.InvalidSmilesException.class,
        () -> MolIO.loadTarget("SMI", "c1/c=c/cc"));
    String smiles = "C1=CC2=N/C1=C\\C1=N/C(=C\\C3=NC4=C(C5=N/C(=C\\2)CC5)CC(=O)C4=C3)C=C1";
    IAtomContainer molecule = Standardiser.standardise(
        MolIO.loadTarget("SMI", smiles), Standardiser.TautomerMode.NONE);
    String output = new MolGraph(molecule).toCanonicalSmiles();
    assertEquals(5, stereoCount(MolIO.loadTarget("SMI", output)));
    ExecutorService executor = Executors.newFixedThreadPool(2);
    try {
      Future<IAtomContainer> first =
          executor.submit(() -> MolIO.loadTarget("SMI", output));
      Future<IAtomContainer> second =
          executor.submit(() -> MolIO.loadQuery("SMI", output).container());
      assertEquals(5, stereoCount(first.get()));
      assertEquals(5, stereoCount(second.get()));
    } finally {
      executor.shutdownNow();
    }
  }

  @ParameterizedTest
  @ValueSource(strings = {
      "C1=CC2=N/C1=C\\C1=N/C(=C\\C3=NC4=C(C5=N/C(=C\\2)CC5)CC(=O)C4=C3)C=C1",
      "C1=CC2=N/C1=C/C1=N/C(=C\\C3=NC4=C(C5=N/C(=C\\2)CC5)CC(=O)C4=C3)C=C1"})
  void completeMatchPreservesAromaticDoubleBondFocusOrders(String smiles) throws Exception {
    IAtomContainer original = Standardiser.standardise(
        MolIO.loadTarget("SMI", smiles), Standardiser.TautomerMode.NONE);
    MolGraph graph = new MolGraph(original);
    Map<Integer, Integer> identity = new TreeMap<>();
    for (int atom = 0; atom < graph.atomCount(); atom++) identity.put(atom, atom);
    SearchResult.Match match = MCS.fragment(graph, identity, MatchGraph.allBonds(graph, identity).keySet());
    IAtomContainer molecule = match.toMolecule(graph);
    assertEquals(5, stereoCount(molecule));
    for (IStereoElement<?, ?> stereo : molecule.stereoElements())
      if (stereo instanceof IDoubleBondStereochemistry)
        assertEquals(IBond.Order.DOUBLE,
            ((IDoubleBondStereochemistry) stereo).getStereoBond().getOrder());
    assertAromaticStereoRoundTrip(new MolGraph(molecule), original);
    IAtomContainer smilesRoundTrip = Standardiser.standardise(
        MolIO.loadTarget("SMI", match.toSmiles(graph)), Standardiser.TautomerMode.NONE);
    assertTrue(Pattern.findIdentical(original).matches(smilesRoundTrip));
    ByteArrayOutputStream bytes = new ByteArrayOutputStream();
    OutputUtil.writeMolecule(molecule,
        OutputUtil.OutType.MOL, bytes);
    IAtomContainer molRoundTrip;
    try (MDLV2000Reader reader = new MDLV2000Reader(
        new StringReader(bytes.toString("UTF-8")))) {
      molRoundTrip = Standardiser.standardise(reader.read(
          new AtomContainer()), Standardiser.TautomerMode.NONE);
    }
    assertEquals(3, carbonDoubleBondStereoCount(molRoundTrip));
    assertTrue(Pattern.findIdentical(carbonDoubleBondStereoOnly(original))
        .matches(carbonDoubleBondStereoOnly(molRoundTrip)));
  }

  @Test
  void completeMatchReservesAromaticFocusAfterAtomPermutation() throws Exception {
    IAtomContainer molecule = MolIO.loadTarget("SMI", "C1=CC=CC=C1");
    IBond focus = molecule.getBond(molecule.getAtom(0), molecule.getAtom(1));
    molecule.addStereoElement(new org.openscience.cdk.stereo.DoubleBondStereochemistry(focus,
        new IBond[] {molecule.getBond(molecule.getAtom(0), molecule.getAtom(5)),
            molecule.getBond(molecule.getAtom(1), molecule.getAtom(2))},
        IDoubleBondStereochemistry.Conformation.TOGETHER));
    molecule = Standardiser.standardise(molecule, Standardiser.TautomerMode.NONE);
    IAtom[] atoms = new IAtom[molecule.getAtomCount()];
    for (int atom = 0; atom < atoms.length; atom++) atoms[atom] = molecule.getAtom((atom + 1) % atoms.length);
    molecule.setAtoms(atoms);
    IBond[] bonds = new IBond[molecule.getBondCount()];
    for (int bond = 0; bond < bonds.length; bond++) bonds[bond] = molecule.getBond(bonds.length - bond - 1);
    molecule.setBonds(bonds);
    MolGraph graph = new MolGraph(molecule);
    Map<Integer, Integer> identity = new TreeMap<>();
    for (int atom = 0; atom < graph.atomCount(); atom++) identity.put(atom, atom);
    IAtomContainer exported = MCS.fragment(graph, identity, MatchGraph.allBonds(graph, identity).keySet())
        .toMolecule(graph);
    for (IStereoElement<?, ?> stereo : exported.stereoElements())
      assertEquals(IBond.Order.DOUBLE,
          ((IDoubleBondStereochemistry) stereo).getStereoBond().getOrder());
    assertAromaticStereoRoundTrip(new MolGraph(exported), molecule);
  }

  @ParameterizedTest
  @ValueSource(ints = {0, 1})
  void builderMatchReservesSpecifiedAromaticFocus(int first) throws Exception {
    int[][] neighbors = {{1, 5}, {0, 2}, {1, 3}, {2, 4}, {3, 5}, {4, 0}};
    int[][] stereo = new int[6][6];
    stereo[first][first + 1] = stereo[first + 1][first] = 1;
    MolGraph graph = new MolGraph.Builder().atomCount(6)
        .atomicNumbers(new int[] {6, 6, 6, 6, 6, 6})
        .hydrogenCounts(new int[] {1, 1, 1, 1, 1, 1})
        .ringFlags(new boolean[] {true, true, true, true, true, true})
        .aromaticFlags(new boolean[] {true, true, true, true, true, true})
        .neighbors(neighbors).bondOrders(new int[][] {{4, 4}, {4, 4}, {4, 4}, {4, 4}, {4, 4}, {4, 4}})
        .bondRingFlags(new boolean[][] {{true, true}, {true, true}, {true, true},
            {true, true}, {true, true}, {true, true}})
        .bondAromaticFlags(new boolean[][] {{true, true}, {true, true}, {true, true},
            {true, true}, {true, true}, {true, true}}).doubleBondStereo(stereo).build();
    Map<Integer, Integer> identity = new TreeMap<>();
    for (int atom = 0; atom < graph.atomCount(); atom++) identity.put(atom, atom);
    IAtomContainer exported = MCS.fragment(graph, identity, MatchGraph.allBonds(graph, identity).keySet())
        .toMolecule(graph);
    assertEquals(IBond.Order.DOUBLE,
        exported.getBond(exported.getAtom(first), exported.getAtom(first + 1)).getOrder());
    assertEquals(1, stereoCount(exported));
    assertAromaticStereoRoundTrip(graph, exported);
  }

  @Test
  void completeMatchCopiesSourceCoordinates() throws Exception {
    IAtomContainer molecule = MolIO.loadTarget("SMI", "F/C=C/F");
    for (int atom = 0; atom < molecule.getAtomCount(); atom++) {
      molecule.getAtom(atom).setPoint2d(new Point2d(atom, atom * 0.5));
      molecule.getAtom(atom).setPoint3d(new Point3d(atom, atom * 0.5, atom * 0.25));
    }
    MolGraph graph = new MolGraph(molecule);
    Map<Integer, Integer> identity = new TreeMap<>();
    for (int atom = 0; atom < graph.atomCount(); atom++) identity.put(atom, atom);
    IAtomContainer exported = MCS.fragment(graph, identity, MatchGraph.allBonds(graph, identity).keySet())
        .toMolecule(graph);
    for (int atom = 0; atom < graph.atomCount(); atom++) {
      assertEquals(molecule.getAtom(atom).getPoint2d(), exported.getAtom(atom).getPoint2d());
      assertEquals(molecule.getAtom(atom).getPoint3d(), exported.getAtom(atom).getPoint3d());
      assertNotSame(graph.mol.getAtom(atom).getPoint2d(), exported.getAtom(atom).getPoint2d());
      assertNotSame(graph.mol.getAtom(atom).getPoint3d(), exported.getAtom(atom).getPoint3d());
      exported.getAtom(atom).getPoint2d().x += 1;
      exported.getAtom(atom).getPoint3d().x += 1;
      assertEquals(molecule.getAtom(atom).getPoint2d(), graph.mol.getAtom(atom).getPoint2d());
      assertEquals(molecule.getAtom(atom).getPoint3d(), graph.mol.getAtom(atom).getPoint3d());
    }
  }

  @ParameterizedTest
  @ValueSource(strings = {"F/C=C/F", "F/C=C\\F", "FC=CF"})
  void molOutputPreservesDefinedAndUnspecifiedAlkeneStereo(String smiles) throws Exception {
    IAtomContainer original = MolIO.loadTarget("SMI", smiles);
    IAtomContainer restored = molRoundTrip(original);
    assertEquals(stereoCount(original), stereoCount(restored));
    assertTrue(Pattern.findIdentical(original).matches(restored));
    for (IAtom atom : original.atoms()) assertNull(atom.getPoint2d());
  }

  @Test
  void molOutputRetainsSuppliedCoordinates() throws Exception {
    IAtomContainer original = MolIO.loadTarget("SMI", "F/C=C/F");
    double[][] coordinates = {{-1, 1}, {0, 0}, {1, 0}, {2, -1}};
    for (int atom = 0; atom < coordinates.length; atom++)
      original.getAtom(atom).setPoint2d(new Point2d(coordinates[atom]));
    IAtomContainer restored = molRoundTrip(original);
    assertEquals(1, stereoCount(restored));
    assertTrue(Pattern.findIdentical(original).matches(restored));
    for (int atom = 0; atom < coordinates.length; atom++) {
      assertEquals(new Point2d(coordinates[atom]), original.getAtom(atom).getPoint2d());
      assertEquals(original.getAtom(atom).getPoint2d(), restored.getAtom(atom).getPoint2d());
    }
  }

  @ParameterizedTest
  @ValueSource(booleans = {false, true})
  void molOutputRejectsContradictoryOrDegenerateCoordinates(boolean degenerate) throws Exception {
    IAtomContainer molecule = MolIO.loadTarget("SMI", "F/C=C/F");
    double[][] coordinates = degenerate
        ? new double[][] {{0, 0}, {0, 0}, {0, 0}, {0, 0}}
        : new double[][] {{-1, 1}, {0, 0}, {1, 0}, {2, 1}};
    for (int atom = 0; atom < coordinates.length; atom++)
      molecule.getAtom(atom).setPoint2d(new Point2d(coordinates[atom]));
    ByteArrayOutputStream bytes = new ByteArrayOutputStream();
    assertThrows(IllegalArgumentException.class,
        () -> OutputUtil.writeMolecule(molecule,
            OutputUtil.OutType.MOL, bytes));
    assertEquals(0, bytes.size());
    for (int atom = 0; atom < coordinates.length; atom++)
      assertEquals(new Point2d(coordinates[atom]), molecule.getAtom(atom).getPoint2d());
  }

  @Test
  void molOutputRejectsGeometryThatCollapsesAtFourDecimals() throws Exception {
    IAtomContainer molecule = MolIO.loadTarget("SMI", "F/C=C/F");
    double[][] coordinates = {{-1, 0.00001}, {0, 0}, {1, 0}, {2, -0.00001}};
    for (int atom = 0; atom < coordinates.length; atom++)
      molecule.getAtom(atom).setPoint2d(new Point2d(coordinates[atom]));
    ByteArrayOutputStream output = new ByteArrayOutputStream();
    assertThrows(IllegalArgumentException.class,
        () -> OutputUtil.writeMolecule(molecule, OutType.MOL, output));
    assertEquals(0, output.size());
    for (int atom = 0; atom < coordinates.length; atom++)
      assertEquals(new Point2d(coordinates[atom]), molecule.getAtom(atom).getPoint2d());
  }

  @ParameterizedTest
  @ValueSource(strings = {"F/C=C/F", "F/C=C\\F"})
  void molOutputUsesSuppliedThreeDimensionalCoordinates(String smiles) throws Exception {
    IAtomContainer source = MolIO.loadTarget("SMI", smiles);
    double[][] coordinates = {{-1, 1, 0}, {0, 0, 0}, {1, 0, 0},
        {2, smiles.contains("\\") ? 1 : -1, 0}};
    for (int atom = 0; atom < coordinates.length; atom++) {
      source.getAtom(atom).setPoint3d(new Point3d(coordinates[atom]));
      source.getAtom(atom).setPoint2d(new Point2d());
    }
    IAtomContainer restored = molRoundTrip(source);
    assertEquals(1, stereoCount(restored));
    assertTrue(Pattern.findIdentical(source).matches(restored));
    for (int atom = 0; atom < coordinates.length; atom++) {
      assertEquals(new Point3d(coordinates[atom]), source.getAtom(atom).getPoint3d());
      assertEquals(source.getAtom(atom).getPoint3d(), restored.getAtom(atom).getPoint3d());
      assertEquals(new Point2d(), source.getAtom(atom).getPoint2d());
    }
  }

  @ParameterizedTest
  @ValueSource(booleans = {false, true})
  void molOutputUsesValidatedTwoDimensionsForIncompleteThreeDimensions(boolean complete2D)
      throws Exception {
    IAtomContainer source = MolIO.loadTarget("SMI", "F/C=C/F");
    source.getAtom(0).setPoint3d(new Point3d());
    double[][] coordinates = {{-1, 1}, {0, 0}, {1, 0}, {2, -1}};
    if (complete2D)
      for (int atom = 0; atom < coordinates.length; atom++)
        source.getAtom(atom).setPoint2d(new Point2d(coordinates[atom]));
    IAtomContainer restored = molRoundTrip(source);
    assertEquals(1, stereoCount(restored));
    assertTrue(Pattern.findIdentical(source).matches(restored));
    assertEquals(new Point3d(), source.getAtom(0).getPoint3d());
    for (int atom = 0; atom < coordinates.length; atom++) {
      assertNotNull(restored.getAtom(atom).getPoint2d());
      assertNull(restored.getAtom(atom).getPoint3d());
      if (complete2D)
        assertEquals(new Point2d(coordinates[atom]), source.getAtom(atom).getPoint2d());
      else assertNull(source.getAtom(atom).getPoint2d());
    }
  }

  private static IAtomContainer molRoundTrip(IAtomContainer molecule) throws Exception {
    ByteArrayOutputStream bytes = new ByteArrayOutputStream() {
      @Override public void close() { fail("The caller's stream must remain open"); }
    };
    OutputUtil.writeMolecule(molecule,
        OutputUtil.OutType.MOL, bytes);
    try (MDLV2000Reader reader = new MDLV2000Reader(
        new StringReader(bytes.toString("UTF-8")))) {
      return reader.read(new AtomContainer());
    }
  }

  @ParameterizedTest
  @ValueSource(strings = {"F/C=C/F", "F/C=C\\F", "[13CH3][C@H](F)[NH3+]",
      "C1=CC2=N/C1=C\\C1=N/C(=C\\C3=NC4=C(C5=N/C(=C\\2)CC5)CC(=O)C4=C3)C=C1"})
  void publicSmilesOutputsPreserveStereoAndMappingLabels(String smiles) throws Exception {
    IAtomContainer source = Standardiser.standardise(
        MolIO.loadTarget("SMI", smiles), Standardiser.TautomerMode.NONE);
    ByteArrayOutputStream output = new ByteArrayOutputStream();
    OutputUtil.writeMolecule(source, OutType.SMI, output);
    IAtomContainer restored = Standardiser.standardise(
        MolIO.loadTarget("SMI", output.toString("UTF-8").trim()), Standardiser.TautomerMode.NONE);
    assertEquals(stereoCount(source), stereoCount(restored));
    assertTrue(Pattern.findIdentical(source).matches(restored));
    source.getAtom(0).setMapIdx(41);
    IAtomContainer target = source.clone();
    Map<Integer, Integer> identity = new TreeMap<>();
    for (int atom = 0; atom < source.getAtomCount(); atom++) identity.put(atom, atom);
    output.reset();
    OutputUtil.writeMappings(source, target, java.util.Collections.singletonList(identity),
        OutType.SMI, output, false);
    String[] columns = output.toString("UTF-8").trim().split("\\t");
    assertEquals(2, columns.length);
    for (String column : columns) {
      restored = Standardiser.standardise(MolIO.loadTarget("SMI", column), Standardiser.TautomerMode.NONE);
      assertEquals(stereoCount(source), stereoCount(restored));
      assertTrue(Pattern.findIdentical(source).matches(restored));
      java.util.Set<Integer> labels = new java.util.HashSet<>();
      for (IAtom atom : restored.atoms()) {
        assertTrue(atom.getMapIdx() > 0);
        assertEquals(source.getAtom(atom.getMapIdx() - 1).getAtomicNumber(), atom.getAtomicNumber());
        labels.add(atom.getMapIdx());
      }
      assertEquals(source.getAtomCount(), labels.size());
    }
    assertEquals(41, source.getAtom(0).getMapIdx());
    assertEquals(41, target.getAtom(0).getMapIdx());
    output.reset();
    OutputUtil.writeMappings(source, target,
        java.util.Collections.singletonList(java.util.Collections.singletonMap(1, 1)),
        OutType.SMI, output, false);
    for (String column : output.toString("UTF-8").trim().split("\\t")) {
      restored = Standardiser.standardise(MolIO.loadTarget("SMI", column), Standardiser.TautomerMode.NONE);
      assertEquals(stereoCount(source), stereoCount(restored));
      assertTrue(Pattern.findIdentical(source).matches(restored));
      java.util.Set<Integer> labels = new java.util.HashSet<>();
      for (IAtom atom : restored.atoms()) if (atom.getMapIdx() > 0) labels.add(atom.getMapIdx());
      assertEquals(new java.util.HashSet<>(java.util.Arrays.asList(1, 41)), labels);
    }
  }

  @Test
  void mappedCanonicalSmilesLabelsDoNotChangeGraphIdentity() throws Exception {
    MolGraph graph = new MolGraph(MolIO.loadTarget("SMI", "[13CH3][C@H](F)[NH3+]"));
    long hash = graph.getCanonicalHash();
    String original = graph.toCanonicalSmiles();
    Map<Integer, Integer> labels = java.util.Collections.singletonMap(1, 71);
    assertTrue(graph.toCanonicalSmiles(labels).contains(":71]"));
    assertEquals(original, graph.toCanonicalSmiles());
    assertEquals(hash, graph.getCanonicalHash());
    assertThrows(NullPointerException.class, () -> graph.toCanonicalSmiles(null));
    assertThrows(IllegalArgumentException.class,
        () -> graph.toCanonicalSmiles(java.util.Collections.singletonMap(-1, 1)));
    assertThrows(IllegalArgumentException.class,
        () -> graph.toCanonicalSmiles(java.util.Collections.singletonMap(graph.atomCount(), 1)));
    assertThrows(IllegalArgumentException.class,
        () -> graph.toCanonicalSmiles(java.util.Collections.singletonMap(0, 0)));
    ExecutorService executor = Executors.newFixedThreadPool(2);
    try {
      Future<String> first = executor.submit(() -> graph.toCanonicalSmiles(labels));
      Future<String> second = executor.submit(
          () -> graph.toCanonicalSmiles(java.util.Collections.singletonMap(1, 72)));
      assertTrue(first.get().contains(":71]"));
      assertTrue(second.get().contains(":72]"));
      assertEquals(original, graph.toCanonicalSmiles());
    } finally {
      executor.shutdownNow();
    }
  }

  @Test
  void legacyMappedMolPreservesRetainedStereoAndDefersMissingLigands() throws Exception {
    IAtomContainer alkene = MolIO.loadTarget("SMI", "F/C=C/F");
    assertEquals(1, stereoCount(mappedMol(alkene, 0, 1, 2, 3)));
    assertEquals(0, stereoCount(mappedMol(alkene, 0, 1, 2)));
    IAtomContainer tetrahedral = MolIO.loadTarget("SMI", "F[C@H](Cl)CCO");
    assertEquals(1, stereoCount(mappedMol(tetrahedral, 0, 1, 2, 3)));
    assertEquals(0, stereoCount(mappedMol(tetrahedral, 1, 2, 3)));
    IAtomContainer labelled = MolIO.loadTarget("SMI", "[13CH3][C@H](F)[NH3+]");
    IAtomContainer restored = mappedMol(labelled, 0, 1, 2, 3);
    assertEquals(1, stereoCount(restored));
    assertEquals(13, restored.getAtom(0).getMassNumber());
    assertEquals(1, restored.getAtom(3).getFormalCharge());
    for (IAtom atom : labelled.atoms()) assertNull(atom.getPoint2d());
  }

  private static IAtomContainer mappedMol(IAtomContainer target, int... atoms) throws Exception {
    Map<Integer, Integer> mapping = new TreeMap<>();
    for (int atom : atoms) mapping.put(atom, atom);
    ByteArrayOutputStream output = new ByteArrayOutputStream();
    OutputUtil.writeMappings(null, target, java.util.Collections.singletonList(mapping),
        OutType.MOL, output, false);
    try (MDLV2000Reader reader = new MDLV2000Reader(new StringReader(output.toString("UTF-8")))) {
      return reader.read(new AtomContainer());
    }
  }

  @Test
  void mcsCliMolExportUsesSelectedCommonEdges(@org.junit.jupiter.api.io.TempDir java.nio.file.Path directory)
      throws Exception {
    java.nio.file.Path output = directory.resolve("mapping.mol"), result = directory.resolve("result.json");
    int status = new picocli.CommandLine(new com.bioinception.smsd.cli.SMSDcli()).execute(
        "--Q", "SMI", "--q", "CCC", "--T", "SMI", "--t", "C1CC1", "--mode", "mcs",
        "--objective", "MCES", "--map-format", "MOL", "--map-out", output.toString(),
        "--output", result.toString());
    assertEquals(0, status);
    IAtomContainer molecule = MolIO.loadTarget("MOL", output.toString());
    assertEquals(3, molecule.getAtomCount());
    assertEquals(2, molecule.getBondCount());
  }

  private static int carbonDoubleBondStereoCount(IAtomContainer molecule) {
    int count = 0;
    for (IStereoElement<?, ?> stereo : molecule.stereoElements())
      if (stereo instanceof IDoubleBondStereochemistry) {
        IBond bond = ((IDoubleBondStereochemistry) stereo).getStereoBond();
        if (bond.getBegin().getAtomicNumber() == 6 && bond.getEnd().getAtomicNumber() == 6) count++;
      }
    return count;
  }

  private static IAtomContainer carbonDoubleBondStereoOnly(IAtomContainer molecule) throws Exception {
    IAtomContainer copy = molecule.clone();
    List<IStereoElement> retained = new ArrayList<>();
    for (IStereoElement<?, ?> stereo : copy.stereoElements()) {
      if (stereo instanceof IDoubleBondStereochemistry) {
        IBond bond = ((IDoubleBondStereochemistry) stereo).getStereoBond();
        if (bond.getBegin().getAtomicNumber() != 6 || bond.getEnd().getAtomicNumber() != 6) continue;
      }
      retained.add(stereo);
    }
    copy.setStereoElements(retained);
    return copy;
  }

  private static String assertAromaticStereoRoundTrip(
      MolGraph graph, IAtomContainer original) throws Exception {
    String output = graph.toCanonicalSmiles();
    IAtomContainer reparsed = Standardiser.standardise(
        MolIO.loadTarget("SMI", output), Standardiser.TautomerMode.NONE);
    assertEquals(stereoCount(original), stereoCount(reparsed), output);
    assertTrue(Pattern.findIdentical(original).matches(reparsed), output);
    assertEquals(graph.getCanonicalHash(), new MolGraph(reparsed).getCanonicalHash(), output);
    return output;
  }

  private static void assertRoundTrip(SmilesParser parser, IAtomContainer original) throws Exception {
    String output = new MolGraph(original).toCanonicalSmiles();
    IAtomContainer reparsed = parser.parseSmiles(output);
    assertEquals(stereoCount(original), stereoCount(reparsed), output);
    assertTrue(Pattern.findIdentical(original).matches(reparsed), output);
  }

  private static int stereoCount(IAtomContainer molecule) {
    int count = 0;
    for (IStereoElement<?, ?> stereo : molecule.stereoElements()) count++;
    return count;
  }
}
