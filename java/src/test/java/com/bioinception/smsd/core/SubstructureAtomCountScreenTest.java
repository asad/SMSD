/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
package com.bioinception.smsd.core;

import static com.bioinception.smsd.TestSupport.rawGraph;
import static org.junit.jupiter.api.Assertions.*;
import java.util.*;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.CsvSource;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.silent.SilentChemObjectBuilder;

/** The element and degree count screen that ends impossible molecule searches before any node. */
class SubstructureAtomCountScreenTest {

  @ParameterizedTest
  @CsvSource({
      "OCO, OCCC",           // two oxygens against one
      "CC(C)(C)C, CCCCCC",   // enough carbons, but none with four neighbours
      "CN, CCC"              // no nitrogen at all
  })
  void shortfallsCompleteWithoutASearchNode(String querySmiles, String targetSmiles) throws Exception {
    MolGraph query = rawGraph(querySmiles), target = rawGraph(targetSmiles);
    assertFalse(new Substructure.QueryTables(query, new ChemOptions()).admits(target));
    Substructure prepared = Substructure.find(query);
    for (SearchResult result : Arrays.asList(prepared.match(target), prepared.matchAll(target))) {
      assertEquals(SearchResult.Status.Completed, result.status());
      assertTrue(result.complete());
      assertTrue(result.optimal());
      assertFalse(result.found());
      assertEquals(0, result.score());
      assertEquals(0, result.nodesVisited());
    }
    // No node is spent, so a node limit cannot stop the proof.
    Substructure.Options options = new Substructure.Options();
    options.nodeLimit = 1;
    SearchResult limited = prepared.withOptions(options).match(target);
    assertEquals(SearchResult.Status.Completed, limited.status());
    assertTrue(limited.complete());
    assertEquals(0, limited.nodesVisited());
  }

  @Test
  void withoutTheElementTestEveryPairIsSearched() throws Exception {
    Substructure.Options options = new Substructure.Options();
    options.chemistry.matchAtomType = false;
    MolGraph query = rawGraph("CN"), target = rawGraph("CCC");
    assertTrue(new Substructure.QueryTables(query, options.chemistry).admits(target));
    SearchResult result = Substructure.find(query).withOptions(options).match(target);
    assertEquals(SearchResult.Status.Completed, result.status());
    assertTrue(result.found());
    assertTrue(result.nodesVisited() > 0);
  }

  @Test
  void anExpiredTimeoutStillTimesOutOnAScreenedPair() throws Exception {
    SearchResult result = Substructure.find(rawGraph("CN")).withTimeout(0).match(rawGraph("CCC"));
    assertEquals(SearchResult.Status.TimedOut, result.status());
    assertFalse(result.complete());
    assertFalse(result.optimal());
    assertEquals(0, result.nodesVisited());
  }

  @Test
  void queriesWithMoreThanSixteenElementsAreSearchedUnscreened() throws Exception {
    String sixteen = "[Li][Be]BCNOF[Na][Mg][Al][Si]PSCl[K]";
    MolGraph screened = rawGraph(sixteen + "[Ca]"), unscreened = rawGraph(sixteen + "[Ca][Sc]");
    MolGraph target = rawGraph(sixteen + "[Sr][Sc]C");
    assertFalse(new Substructure.QueryTables(screened, new ChemOptions()).admits(target));
    assertEquals(0, Substructure.find(screened).match(target).nodesVisited());
    assertTrue(new Substructure.QueryTables(unscreened, new ChemOptions()).admits(target));
    SearchResult result = Substructure.find(unscreened).match(target);
    assertEquals(SearchResult.Status.Completed, result.status());
    assertFalse(result.found());
    assertTrue(result.nodesVisited() > 0);
  }

  /** The chemistry of each benchmark profile. */
  private static Map<String, ChemOptions> profiles() {
    Map<String, ChemOptions> profiles = new LinkedHashMap<>();
    for (String name : Arrays.asList("default", "rings", "strict", "any", "tautomer", "stereo")) {
      ChemOptions options = new ChemOptions();
      options.matchBondOrder = ChemOptions.BondOrderMode.STRICT;
      if (name.equals("rings")) options.completeRingsOnly = true;
      if (name.equals("rings") || name.equals("strict")) {
        options.aromaticityMode = ChemOptions.AromaticityMode.STRICT;
        options.matchFormalCharge = options.matchIsotope = options.ringMatchesRingOnly = true;
      }
      if (name.equals("any")) options.matchBondOrder = ChemOptions.BondOrderMode.ANY;
      if (name.equals("tautomer")) {
        options.tautomerAware = true;
        options.matchBondOrder = ChemOptions.BondOrderMode.LOOSE;
      }
      if (name.equals("stereo")) options.useChirality = options.useBondStereo = true;
      profiles.put(name, options);
    }
    return profiles;
  }

  /** A random connected molecule of mostly carbon, with branches and a few ring closures. */
  private static IAtomContainer randomMolecule(Random random, int atoms) {
    String[] elements = {"C", "C", "C", "C", "N", "O", "S"};
    IAtomContainer molecule = SilentChemObjectBuilder.getInstance().newAtomContainer();
    for (int atom = 0; atom < atoms; atom++) {
      IAtom created = molecule.getBuilder().newInstance(IAtom.class, elements[random.nextInt(elements.length)]);
      created.setImplicitHydrogenCount(0);
      molecule.addAtom(created);
      if (atom > 0)
        molecule.addBond(random.nextInt(atom), atom, random.nextInt(5) == 0 ? IBond.Order.DOUBLE : IBond.Order.SINGLE);
    }
    for (int closure = random.nextInt(3); closure > 0; closure--) {
      int first = random.nextInt(atoms), second = random.nextInt(atoms);
      if (first != second && molecule.getBond(molecule.getAtom(first), molecule.getAtom(second)) == null)
        molecule.addBond(first, second, IBond.Order.SINGLE);
    }
    return molecule;
  }

  /** A connected piece of the molecule with shuffled atoms, sometimes with one element changed. */
  private static IAtomContainer randomPiece(IAtomContainer molecule, Random random, int atoms) {
    List<Integer> chosen = new ArrayList<>();
    Set<Integer> taken = new HashSet<>();
    int start = random.nextInt(molecule.getAtomCount());
    chosen.add(start);
    taken.add(start);
    for (int index = 0; index < chosen.size() && chosen.size() < atoms; index++)
      for (IAtom neighbor : molecule.getConnectedAtomsList(molecule.getAtom(chosen.get(index))))
        if (chosen.size() < atoms && taken.add(molecule.indexOf(neighbor))) chosen.add(molecule.indexOf(neighbor));
    Collections.shuffle(chosen, random);
    IAtomContainer piece = SilentChemObjectBuilder.getInstance().newAtomContainer();
    Map<Integer, Integer> position = new HashMap<>();
    for (int atom : chosen) {
      IAtom source = molecule.getAtom(atom);
      String element = random.nextInt(4) == 0 ? (random.nextBoolean() ? "N" : "O") : source.getSymbol();
      IAtom copy = piece.getBuilder().newInstance(IAtom.class, element);
      copy.setImplicitHydrogenCount(0);
      position.put(atom, piece.getAtomCount());
      piece.addAtom(copy);
    }
    for (IBond bond : molecule.bonds()) {
      Integer first = position.get(molecule.indexOf(bond.getBegin())), second = position.get(molecule.indexOf(bond.getEnd()));
      if (first != null && second != null) piece.addBond(first, second, bond.getOrder());
    }
    return piece;
  }

  @Test
  void theScreenOnlyRejectsPairsWithoutAnEmbedding() throws Exception {
    Random random = new Random(20261020L);
    int rejected = 0, found = 0;
    for (int pair = 0; pair < 400; pair++) {
      IAtomContainer targetMolecule = randomMolecule(random, 4 + random.nextInt(14));
      IAtomContainer queryMolecule = random.nextInt(3) == 0
          ? randomMolecule(random, 1 + random.nextInt(targetMolecule.getAtomCount()))
          : randomPiece(targetMolecule, random, 1 + random.nextInt(targetMolecule.getAtomCount()));
      MolGraph query = new MolGraph(queryMolecule), target = new MolGraph(targetMolecule);
      for (Map.Entry<String, ChemOptions> profile : profiles().entrySet())
        for (boolean induced : new boolean[] {false, true}) {
          ChemOptions chemistry = profile.getValue();
          Substructure.QueryTables tables = new Substructure.QueryTables(query, chemistry);
          Substructure.SubgraphSearch search = new Substructure.SubgraphSearch(
              query, tables, target, chemistry, new SearchResult.Control(60_000, null));
          search.induced = induced;
          boolean exists = search.exists();
          boolean admitted = tables.admits(target);
          String context = profile.getKey() + " induced " + induced + " pair " + pair;
          if (!admitted) {
            rejected++;
            assertFalse(exists, context);
          }
          if (exists) found++;
          Substructure.Options options = new Substructure.Options();
          options.chemistry = chemistry;
          options.induced = induced;
          SearchResult result = Substructure.find(query).withOptions(options).match(target);
          assertEquals(SearchResult.Status.Completed, result.status(), context);
          assertEquals(exists, result.found(), context);
          if (!admitted) assertEquals(0, result.nodesVisited(), context);
        }
    }
    assertTrue(rejected > 100, "rejected " + rejected);
    assertTrue(found > 100, "found " + found);
  }
}
