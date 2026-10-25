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

/** The stereo test of each new pair in the general MCIS search. */
class VertexSearchStereoTest {
  /**
   * Centres with an implicit hydrogen, four neighbours or a lone pair, partly configured
   * molecules, ring centres, and double bonds in chains, with three substituents and in rings.
   */
  static final String[] MOLECULES = {
      "C[C@H](O)CC", "C[C@@H](O)CC", "C[C@H](O)C(C)O", "C[C@@H](O)[C@H](C)O", "C[C@H](O)[C@H](C)O",
      "C[C@](F)(Cl)CO", "C[S@](=O)CC", "C[C@H]1CC[C@@H]1O", "C[C@H]1CC[C@H]1O", "C/C=C/C",
      "C/C=C\\C", "C/C(Cl)=C/C", "CC=CC", "C/C=C/[C@H](C)O", "C/C=C\\[C@@H](C)O", "O/C=C/C=C\\C",
      "C/C=C1/CC[C@H]1C"};

  /** The atom mapping and the selected bonds of {@code match}, in a form that compares as text. */
  static String key(SearchResult.Match match) {
    List<String> bonds = new ArrayList<>();
    for (SearchResult.Bond bond : match.bondMapping().keySet())
      bonds.add(bond.toString());
    Collections.sort(bonds);
    return match.atomMapping() + " " + bonds;
  }

  /**
   * Every maximum induced common subgraph that {@link MCS#validateMapping} accepts. The atom maps
   * are enumerated with only the induced-topology and bond tests, which every part of a valid
   * MCIS passes, and the stereo test is left to the validation.
   */
  static Set<String> oracle(MolGraph query, MolGraph target, MCS.Options options) {
    List<List<int[]>> bySize = new ArrayList<>();
    for (int size = 0; size <= query.n; size++)
      bySize.add(new ArrayList<>());
    int[] mapping = new int[query.n];
    Arrays.fill(mapping, -1);
    enumerate(query, target, options, 0, 0, mapping, new boolean[target.n], bySize);
    for (int size = query.n; size > 0; size--) {
      Set<String> keys = new TreeSet<>();
      int bonds = -1;
      for (int[] candidate : bySize.get(size)) {
        Map<Integer, Integer> atoms = new TreeMap<>();
        List<SearchResult.Bond> selected = new ArrayList<>();
        for (int atom = 0; atom < query.n; atom++)
          if (candidate[atom] >= 0) {
            atoms.put(atom, candidate[atom]);
            for (int neighbor : query.neighbors[atom])
              if (neighbor < atom && candidate[neighbor] >= 0)
                selected.add(new SearchResult.Bond(neighbor, atom));
          }
        SearchResult.Match match = MCS.fragment(query, atoms, selected);
        if (selected.size() < bonds
            || !MCS.validateMapping(query, target, match, options).isEmpty())
          continue;
        if (selected.size() > bonds) {
          keys.clear();
          bonds = selected.size();
        }
        keys.add(key(match));
      }
      if (!keys.isEmpty())
        return keys;
    }
    return Collections.emptySet();
  }

  private static void enumerate(MolGraph query, MolGraph target, MCS.Options options, int atom,
      int size, int[] mapping, boolean[] used, List<List<int[]>> bySize) {
    if (atom == query.n) {
      bySize.get(size).add(mapping.clone());
      return;
    }
    enumerate(query, target, options, atom + 1, size, mapping, used, bySize);
    for (int image = 0; image < target.n; image++) {
      if (used[image]
          || !MolGraph.ChemOps.atomsCompatFast(query, atom, target, image, options.chemistry))
        continue;
      boolean induced = true;
      for (int earlier = 0; earlier < atom && induced; earlier++) {
        if (mapping[earlier] < 0)
          continue;
        boolean bonded = query.hasBond(atom, earlier);
        induced = bonded == target.hasBond(image, mapping[earlier])
            && (!bonded || MolGraph.ChemOps.bondsCompatible(
                    query, atom, earlier, target, image, mapping[earlier], options.chemistry));
      }
      if (!induced)
        continue;
      mapping[atom] = image;
      used[image] = true;
      enumerate(query, target, options, atom + 1, size + 1, mapping, used, bySize);
      used[image] = false;
      mapping[atom] = -1;
    }
  }

  static MCS.Options options(boolean chirality, boolean bondStereo, MCS.Connectivity connectivity) {
    MCS.Options options = new MCS.Options();
    options.chemistry.useChirality = chirality;
    options.chemistry.useBondStereo = bondStereo;
    options.connectivity = connectivity;
    options.timeoutMillis = 60_000;
    return options;
  }

  /** matchAll() must give every mapping of the enumeration, and match() one of them. */
  private static Set<String> assertEnumerated(
      MolGraph query, MolGraph target, MCS.Options options, String label) {
    Set<String> expected = oracle(query, target, options);
    MCS prepared = MCS.find(query).withOptions(options);
    SearchResult all = prepared.matchAll(target, 100_000);
    assertEquals(SearchResult.Status.Completed, all.status(), label);
    assertTrue(all.optimal() && all.complete(), label);
    Set<String> actual = new TreeSet<>();
    for (SearchResult.Match match : all.matches())
      assertTrue(actual.add(key(match)), label);
    assertEquals(expected, actual, label);
    SearchResult one = prepared.match(target);
    assertTrue(one.optimal(), label);
    assertEquals(!expected.isEmpty(), one.found(), label);
    if (one.found())
      assertTrue(expected.contains(key(one.matches().get(0))), label);
    return actual;
  }

  @Test
  void everyStereoMcisAgreesWithExhaustiveEnumeration() throws Exception {
    boolean[][] stereo = {{true, false}, {false, true}, {true, true}};
    int searches = 0;
    for (String querySmiles : MOLECULES)
      for (String targetSmiles : MOLECULES) {
        MolGraph query = rawGraph(querySmiles), target = rawGraph(targetSmiles);
        for (boolean[] flags : stereo)
          for (MCS.Connectivity connectivity : MCS.Connectivity.values()) {
            assertEnumerated(query, target, options(flags[0], flags[1], connectivity),
                querySmiles + " " + targetSmiles + " " + Arrays.toString(flags) + " "
                    + connectivity);
            searches++;
          }
      }
    assertEquals(MOLECULES.length * MOLECULES.length * 6, searches);
  }

  /**
   * An alkene from the Builder: atom 0 doubly bonded to atom 1 in the given configuration,
   * single bonds elsewhere and no hydrogens, with the neighbours in the order given.
   */
  private static MolGraph alkene(int[] elements, int[][] neighbors, int configuration) {
    int n = elements.length;
    int[][] orders = new int[n][];
    for (int atom = 0; atom < n; atom++) {
      orders[atom] = new int[neighbors[atom].length];
      for (int position = 0; position < orders[atom].length; position++)
        orders[atom][position] = atom <= 1 && neighbors[atom][position] <= 1 ? 2 : 1;
    }
    int[][] stereo = new int[n][n];
    stereo[0][1] = stereo[1][0] = configuration;
    return new MolGraph.Builder().atomCount(n).atomicNumbers(elements).neighbors(neighbors)
        .bondOrders(orders).hydrogenCounts(new int[n]).doubleBondStereo(stereo).build();
  }

  /**
   * Double bonds with an end of two substituents in one graph and three in the other. The test
   * reads the first mapped substituent at each end, and here a later one can change its answer,
   * so these bonds are left to the test of complete mappings; testing them as each pair is made
   * loses a maximum in the first three pairs and a tied mapping in the others. The C++ tests
   * use the same graphs.
   */
  @Test
  void crowdedDoubleBondEndsAreLeftToTheMappingTest() {
    Object[][] pairs = {
        {new int[] {6, 6, 6, 7, 6, 6}, new int[][] {{1, 2, 3}, {0, 5}, {0}, {0, 4}, {3}, {1}}, 1,
            new int[] {6, 6, 6, 6, 6, 6, 7, 6},
            new int[][] {{4, 6, 2, 1}, {0, 7}, {3, 0}, {2}, {5, 0}, {4}, {0}, {1}}, 2},
        {new int[] {6, 6, 7, 6, 6, 6, 6},
            new int[][] {{1, 2, 4, 3}, {6, 0}, {0}, {0}, {5, 0}, {4}, {1}}, 1,
            new int[] {6, 6, 6, 6, 6, 6},
            new int[][] {{1, 4, 2}, {5, 0}, {0, 3}, {2}, {0}, {1}}, 1},
        {new int[] {6, 6, 7, 6, 6, 7, 6},
            new int[][] {{4, 1, 2}, {6, 0, 5}, {0, 3}, {2}, {0}, {1}, {1}}, 1,
            new int[] {6, 6, 6, 7, 6, 7, 6, 6},
            new int[][] {{1, 3, 2}, {6, 0, 4, 5}, {0}, {0}, {1}, {1}, {7, 1}, {6}}, 1},
        {new int[] {6, 6, 6, 6, 6, 6}, new int[][] {{1, 2}, {0, 4, 3}, {0}, {1}, {5, 1}, {4}}, 2,
            new int[] {6, 6, 6, 6, 6, 6, 6},
            new int[][] {{1, 2}, {0, 3, 4, 5}, {0}, {1}, {1}, {6, 1}, {5}}, 1},
        {new int[] {6, 6, 6, 6, 6, 6, 6, 6},
            new int[][] {{2, 1}, {7, 0, 5, 4}, {3, 0}, {2}, {1}, {1, 6}, {5}, {1}}, 2,
            new int[] {6, 6, 6, 6, 6, 6},
            new int[][] {{2, 3, 1}, {5, 0}, {0}, {4, 0}, {3}, {1}}, 2},
        {new int[] {6, 6, 7, 7, 6, 7, 6},
            new int[][] {{2, 3, 1}, {6, 5, 0, 4}, {0}, {0}, {1}, {1}, {1}}, 1,
            new int[] {6, 6, 6, 7, 7, 6, 7},
            new int[][] {{1, 4, 3, 2}, {6, 5, 0}, {0}, {0}, {0}, {1}, {1}}, 1}};
    for (int index = 0; index < pairs.length; index++) {
      Object[] pair = pairs[index];
      MolGraph query = alkene((int[]) pair[0], (int[][]) pair[1], (int) pair[2]);
      MolGraph target = alkene((int[]) pair[3], (int[][]) pair[4], (int) pair[5]);
      for (MCS.Connectivity connectivity : MCS.Connectivity.values())
        assertEnumerated(query, target, options(false, true, connectivity),
            "pair " + index + " " + connectivity);
    }
  }

  /** Small tetrahedral stars with supplied hydrogen counts and reordered ligand rows. */
  private static MolGraph tetrahedralStar(int ligands, boolean explicitHydrogen, boolean permuted) {
    int atoms = ligands + 1;
    int[] elements = {6, 9, 17, 35, explicitHydrogen ? 1 : 53};
    int[] labels = new int[atoms], hydrogens = new int[atoms], tags = new int[atoms];
    int[][] neighbors = new int[atoms][], orders = new int[atoms][];
    int[] image = new int[atoms], ligandOrder = new int[ligands];
    for (int atom = 0; atom < atoms; atom++) {
      image[atom] = permuted ? (atom == 0 ? atoms - 1 : atom - 1) : atom;
      labels[image[atom]] = elements[atom];
    }
    for (int ligand = 0; ligand < ligands; ligand++) ligandOrder[ligand] = ligand + 1;
    if (permuted) { ligandOrder[0] = 2; ligandOrder[1] = 1; }
    int centre = image[0], configuration = explicitHydrogen ? 2 : 1;
    hydrogens[centre] = ligands == 3 ? 1 : 0;
    tags[centre] = permuted && ligands == 4 ? 3 - configuration : configuration;
    neighbors[centre] = new int[ligands];
    orders[centre] = new int[ligands];
    for (int position = 0; position < ligands; position++) {
      int ligand = image[ligandOrder[position]];
      neighbors[centre][position] = ligand;
      orders[centre][position] = 1;
      neighbors[ligand] = new int[] {centre};
      orders[ligand] = new int[] {1};
    }
    return new MolGraph.Builder().atomCount(atoms).atomicNumbers(labels).hydrogenCounts(hydrogens)
        .tetrahedralChirality(tags).neighbors(neighbors).bondOrders(orders).build();
  }

  @Test
  void relaxedLabelsAndHydrogenPermutationsRetainEveryOptimum() {
    for (MCS.Connectivity connectivity : MCS.Connectivity.values()) {
      for (boolean stereo : new boolean[] {false, true}) {
        MCS.Options options = options(stereo, stereo, connectivity);
        options.chemistry.matchAtomType = false;
        for (int kind = 0; kind < 3; kind++) {
          int ligands = kind == 0 ? 3 : 4;
          boolean explicitHydrogen = kind == 2;
          MolGraph query = tetrahedralStar(ligands, explicitHydrogen, false);
          MolGraph target = tetrahedralStar(ligands, explicitHydrogen, true);
          int permutations = ligands == 3 ? 6 : 24;
          Set<String> actual = assertEnumerated(query, target, options, "relaxed tetrahedral permutations");
          assertEquals(stereo ? permutations / 2 : permutations, actual.size());
        }
        MolGraph query = alkene(new int[] {6, 6, 9, 17, 35, 53},
            new int[][] {{1, 2, 3}, {0, 4, 5}, {0}, {0}, {1}, {1}}, 1);
        MolGraph target = alkene(new int[] {6, 6, 9, 17, 35, 53},
            new int[][] {{3, 2, 1}, {5, 4, 0}, {0}, {0}, {1}, {1}}, 1);
        Set<String> actual = assertEnumerated(query, target, options, "relaxed alkene permutations");
        assertEquals(stereo ? 4 : 8, actual.size());
      }
      MolGraph implicit = tetrahedralStar(3, false, false), explicitH = tetrahedralStar(4, true, true);
      MCS.Options options = options(true, true, connectivity);
      assertEnumerated(implicit, explicitH, options, "implicit to explicit hydrogen");
      assertEnumerated(explicitH, implicit, options, "explicit to implicit hydrogen");
    }
  }

  @Test
  void stereoStopsRetainOnlyFeasibleWitnesses() {
    MolGraph query = tetrahedralStar(4, false, false), target = tetrahedralStar(4, false, true);
    for (boolean all : new boolean[] {false, true}) {
      for (long cap : new long[] {1, 2}) {
        MCS.Options options = options(true, true, MCS.Connectivity.Connected);
        options.chemistry.matchAtomType = false;
        options.nodeLimit = cap;
        MCS prepared = MCS.find(query).withOptions(options);
        SearchResult stopped = all ? prepared.matchAll(target, 100_000) : prepared.match(target);
        assertEquals(SearchResult.Status.LimitReached, stopped.status());
        assertFalse(stopped.optimal() || stopped.complete());
        for (SearchResult.Match match : stopped.matches())
          assertTrue(MCS.validateMapping(query, target, match, options).isEmpty());
      }
      MCS.Options options = options(true, true, MCS.Connectivity.Connected);
      options.cancellation = () -> true;
      MCS prepared = MCS.find(query).withOptions(options);
      SearchResult cancelled = all ? prepared.matchAll(target, 100_000) : prepared.match(target);
      assertEquals(SearchResult.Status.Cancelled, cancelled.status());
      assertFalse(cancelled.optimal() || cancelled.complete());
      assertFalse(cancelled.found());
    }
  }

  /**
   * Epimers from the benchmark pairs: doxorubicin and epirubicin, two fluvastatin diastereomers
   * and two cefotaxime epimers, with chirality and double-bond stereo. A pair that breaks a
   * centre is passed over before its node, which keeps these searches small. The C++ tests pin
   * the native counts of the same searches, whose results are the same.
   */
  @Test
  void pairsThatBreakACentreArePassedOver() throws Exception {
    String doxorubicin = "COc1cccc2c1C(=O)c1c([O-])c3c(c([O-])c1C2=O)C[C@@](O)(C(=O)CO)C[C@@H]3"
        + "O[C@H]1C[C@H]([NH3+])[C@@H](O)[C@H](C)O1";
    String epirubicin = "COc1cccc2c1C(=O)c1c([O-])c3c(c([O-])c1C2=O)C[C@@](O)(C(=O)CO)C[C@@H]3"
        + "O[C@H]1C[C@H]([NH3+])[C@H](O)[C@H](C)O1";
    String fluvastatin = "CC(C)n1c(/C=C/[C@H](O)C[C@H](O)CC(=O)[O-])c(-c2ccc(F)cc2)c2ccccc21";
    String diastereomer = "CC(C)n1c(/C=C/[C@@H](O)C[C@H](O)CC(=O)O)c(-c2ccc(F)cc2)c2ccccc21";
    String cefotaxime = "CO/N=C(\\C(=O)N[C@H]1C(=O)N2C(C(=O)[O-])=C(COC(C)=O)CS[C@@H]12)"
        + "c1csc(N)n1.[Na+]";
    String epimer = "CO/N=C(\\C(=O)N[C@@H]1C(=O)N2C(C(=O)[O-])=C(COC(C)=O)CS[C@H]12)c1csc(N)n1";
    // Query, target, connectivity, atoms, bonds, match() nodes, mappings, matchAll() nodes.
    Object[][] cases = {
        {doxorubicin, epirubicin, MCS.Connectivity.Connected, 38, 42, 148L, 1, 257L},
        {doxorubicin, epirubicin, MCS.Connectivity.Fragments, 38, 42, 166L, 1, 259L},
        {fluvastatin, diastereomer, MCS.Connectivity.Fragments, 29, 31, 186L, 4, 236L},
        {cefotaxime, epimer, MCS.Connectivity.Connected, 29, 29, 156L, 1, 193L}};
    for (Object[] row : cases) {
      MolGraph query = rawGraph((String) row[0]), target = rawGraph((String) row[1]);
      MCS prepared =
          MCS.find(query).withOptions(options(true, true, (MCS.Connectivity) row[2]));
      String label = row[0] + " " + row[1] + " " + row[2];
      SearchResult one = prepared.match(target), all = prepared.matchAll(target, 1000);
      assertTrue(one.optimal() && all.optimal() && all.complete(), label);
      assertEquals(row[3], one.matches().get(0).atomCount(), label);
      assertEquals(row[4], one.matches().get(0).bondCount(), label);
      assertEquals(row[5], one.nodesVisited(), label);
      assertEquals(row[6], all.matches().size(), label);
      assertEquals(row[7], all.nodesVisited(), label);
      Set<String> keys = new HashSet<>();
      for (SearchResult.Match match : all.matches())
        keys.add(key(match));
      assertTrue(keys.contains(key(one.matches().get(0))), label);
    }
  }
}
