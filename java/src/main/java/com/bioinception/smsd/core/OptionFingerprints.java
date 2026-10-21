/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
package com.bioinception.smsd.core;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.ConcurrentHashMap;

/**
 * The fingerprints of {@link FingerprintOptions}, term for term as the native core computes
 * them: the classic configuration of each family gives the bits of the fixed
 * {@link FingerprintEngine} methods, and every further switch extends the hashed terms.
 */
final class OptionFingerprints {
  private static final long SEED = FingerprintEngine.FNV1A_SEED;
  private static final long PRIME = FingerprintEngine.FNV1A_PRIME;
  /** Mixing tags of the stereo descriptors, so a descriptor never collides with a small term. */
  private static final long CIP_ATOM_TAG = 0xC1A000L, CIP_BOND_TAG = 0xE2E000L;
  /** Neighbour-term marks of the aromatic and ring bond terms of the circular families. */
  private static final long AROMATIC_BOND_MARK = 0x9E3779B97F4A7C15L, RING_BOND_MARK = 0xC2B2AE3D27D4EB4FL;
  /** Embeddings examined per key; a key whose count is not settled within them is an error, never a guess. */
  private static final int KEY_MATCH_LIMIT = 100000;
  /**
   * Time allowed to one key search, in milliseconds, so that it normally completes; a search that
   * does time out is an error rather than a fingerprint that depends on the machine.
   */
  private static final long KEY_SEARCH_MILLIS = 3600000L;

  private OptionFingerprints() {}

  private static long mix(long hash, long value) {
    hash ^= value;
    return hash * PRIME;
  }

  static void validate(FingerprintOptions options) {
    if (options == null) throw new NullPointerException("options");
    if (options.family == null || options.fragments == null || options.scaffold == null)
      throw new NullPointerException("family, fragments and scaffold must be set");
    if (options.family == FingerprintOptions.Family.KEYS) {
      if (options.keys == null || options.keys.isEmpty())
        throw new IllegalArgumentException("keys must name at least one key");
      for (FingerprintOptions.Key key : options.keys) {
        if (key == null) throw new NullPointerException("key");
        if (key.minimum < 1) throw new IllegalArgumentException("a key minimum must be positive");
        if (key.kind == FingerprintOptions.Key.Kind.SMARTS && key.smarts.isEmpty())
          throw new IllegalArgumentException("a SMARTS key must have a pattern");
        if (key.kind == FingerprintOptions.Key.Kind.SMARTS) preparedKey(key.smarts); // a pattern that does not parse fails here
      }
    } else if (options.size <= 0) {
      throw new IllegalArgumentException("size must be positive");
    }
    if (options.radius < -1) throw new IllegalArgumentException("radius must be >= -1 (use -1 for unlimited)");
    if (options.pathLength <= 0 || options.pathLength == Integer.MAX_VALUE)
      throw new IllegalArgumentException("pathLength must be positive and representable");
    if (options.ringSizeCap < 3) throw new IllegalArgumentException("ringSizeCap must be at least 3");
  }

  // ---- Scaffold and fragment handling ----

  /** Connected components in graph order; {@code count[0]} receives their number. */
  static int[] components(MolGraph g, int[] count) {
    int[] component = new int[g.n];
    Arrays.fill(component, -1);
    int[] queue = new int[g.n];
    int found = 0;
    for (int start = 0; start < g.n; start++) {
      if (component[start] >= 0) continue;
      int head = 0, tail = 0;
      queue[tail++] = start;
      component[start] = found;
      while (head < tail) {
        int atom = queue[head++];
        for (int neighbour : g.neighbors[atom])
          if (component[neighbour] < 0) { component[neighbour] = found; queue[tail++] = neighbour; }
      }
      found++;
    }
    count[0] = found;
    return component;
  }

  /**
   * The induced graph of {@code atoms}, in index order. Each kept atom receives
   * {@code extraHydrogens[atom]} hydrogens on top of its implicit count. Tetrahedral tags and
   * double-bond configurations are kept for atoms and bonds that lose no neighbour, since their
   * meaning depends on the neighbour order, which is preserved.
   */
  static MolGraph inducedGraph(MolGraph g, int[] atoms, int[] extraHydrogens) {
    int count = atoms.length;
    int[] newIndex = new int[g.n];
    Arrays.fill(newIndex, -1);
    for (int k = 0; k < count; k++) newIndex[atoms[k]] = k;
    int[] atomicNum = new int[count], formalCharge = new int[count], massNumber = new int[count];
    int[] hydrogens = new int[count];
    boolean[] ring = new boolean[count], aromatic = new boolean[count];
    int[][] neighbors = new int[count][], orders = new int[count][];
    boolean[][] ringBonds = new boolean[count][], aromaticBonds = new boolean[count][];
    int[] tags = g.tetraChirality != null ? new int[count] : null;
    for (int k = 0; k < count; k++) {
      int atom = atoms[k];
      atomicNum[k] = g.atomicNum[atom];
      formalCharge[k] = g.formalCharge[atom];
      massNumber[k] = g.massNumber[atom];
      hydrogens[k] = g.hydrogenCount(atom) + extraHydrogens[atom];
      ring[k] = g.ring[atom];
      aromatic[k] = g.aromatic[atom];
      int kept = 0;
      for (int neighbour : g.neighbors[atom]) if (newIndex[neighbour] >= 0) kept++;
      neighbors[k] = new int[kept];
      orders[k] = new int[kept];
      ringBonds[k] = new boolean[kept];
      aromaticBonds[k] = new boolean[kept];
      int at = 0;
      for (int neighbour : g.neighbors[atom]) {
        if (newIndex[neighbour] < 0) continue;
        neighbors[k][at] = newIndex[neighbour];
        orders[k][at] = g.bondOrder(atom, neighbour);
        ringBonds[k][at] = g.bondInRing(atom, neighbour);
        aromaticBonds[k][at] = g.bondAromatic(atom, neighbour);
        at++;
      }
      if (tags != null) tags[k] = reducedTag(g, atom, k, newIndex, neighbors[k]);
    }
    int[][] stereo = null;
    for (int k = 0; k < count; k++) {
      int atom = atoms[k];
      for (int neighbour : g.neighbors[atom]) {
        if (neighbour <= atom || newIndex[neighbour] < 0 || g.bondOrder(atom, neighbour) != 2) continue;
        int configuration = reducedConfiguration(g, atom, neighbour, newIndex);
        if (configuration == 0) continue;
        if (stereo == null) stereo = new int[count][];
        int other = newIndex[neighbour];
        if (stereo[k] == null) stereo[k] = new int[count];
        if (stereo[other] == null) stereo[other] = new int[count];
        stereo[k][other] = configuration;
        stereo[other][k] = configuration;
      }
    }
    MolGraph.Builder builder = new MolGraph.Builder().atomCount(count).atomicNumbers(atomicNum)
        .formalCharges(formalCharge).massNumbers(massNumber).hydrogenCounts(hydrogens)
        .ringFlags(ring).aromaticFlags(aromatic).neighbors(neighbors).bondOrders(orders)
        .bondRingFlags(ringBonds).bondAromaticFlags(aromaticBonds);
    if (tags != null) builder.tetrahedralChirality(tags);
    if (stereo != null) builder.doubleBondStereo(stereo);
    return builder.build();
  }

  /**
   * The tetrahedral tag of atom {@code atom} as the reduced graph reads it, in which it is atom
   * {@code k} with the kept neighbours {@code neighbors}. A removed neighbour becomes the implicit
   * hydrogen that caps the atom, so the ligand order the tag winds (the CDK element's order, or the
   * derived one) is re-read against the order the reduced graph derives, and the tag flips when the
   * two differ by an odd permutation. An atom that is no longer a tetrahedral centre has no tag.
   */
  private static int reducedTag(MolGraph g, int atom, int k, int[] newIndex, int[] neighbors) {
    int tag = g.tetraChirality[atom];
    int[] ligands = tag == 0 ? null : g.tetrahedralLigands(atom);
    if (ligands == null) return 0;
    int[] actual = new int[4];
    int implicitCount = 0;
    for (int slot = 0; slot < 4; slot++) {
      actual[slot] = ligands[slot] < 0 || newIndex[ligands[slot]] < 0 ? -1 : newIndex[ligands[slot]];
      if (actual[slot] < 0) implicitCount++;
    }
    int degree = neighbors.length;
    if (!((degree == 4 && implicitCount == 0) || (degree == 3 && implicitCount == 1))) return 0;
    int implicit = degree == 4 ? -1 : neighbors[0] < k ? 1 : 0;
    int[] derived = new int[4];
    for (int slot = 0, position = 0; slot < 4; slot++) derived[slot] = slot == implicit ? -1 : neighbors[position++];
    return oddPermutation(actual, derived) ? 3 - tag : tag;
  }

  /**
   * The configuration of the double bond {@code atom}-{@code neighbour} as the reduced graph reads
   * it, or 0 when it is dropped. The stored relation refers to the first other neighbour at each
   * end; when that neighbour is removed the next kept one takes its place and the relation inverts,
   * and when none is left the bond is no longer stereogenic.
   */
  private static int reducedConfiguration(MolGraph g, int atom, int neighbour, int[] newIndex) {
    int configuration = g.dbStereo(atom, neighbour);
    if (configuration == 0) return 0;
    int[][] ends = {{atom, neighbour}, {neighbour, atom}};
    for (int[] end : ends) {
      int reference = g.doubleBondReference(end[0], end[1]);
      if (reference < 0) return 0;
      if (newIndex[reference] >= 0) continue;
      int next = -1;
      for (int other : g.neighbors[end[0]]) if (other != end[1] && newIndex[other] >= 0) { next = other; break; }
      if (next < 0) return 0;
      configuration = 3 - configuration;
    }
    return configuration;
  }

  /** Whether the two orders of the same four ligands differ by an odd permutation. */
  private static boolean oddPermutation(int[] actual, int[] derived) {
    int[] position = new int[4];
    for (int i = 0; i < 4; i++) {
      position[i] = -1;
      for (int j = 0; j < 4; j++) if (derived[j] == actual[i]) { position[i] = j; break; }
      if (position[i] < 0) return false;
    }
    int inversions = 0;
    for (int i = 0; i < 4; i++) for (int j = i + 1; j < 4; j++) if (position[i] > position[j]) inversions++;
    return (inversions & 1) == 1;
  }

  /**
   * The atoms of the Bemis-Murcko framework: ring atoms, the atoms on paths between rings, and the
   * atoms joined to those by a double or triple bond, which the framework keeps as the usual
   * convention does (an exocyclic carbonyl stays; a side chain is cut at its single bond). Chain
   * atoms are pruned from the leaves inwards in linear time; a chain atom survives exactly when two
   * directions from it lead to ring atoms, that is, when it lies on a path between rings.
   */
  static int[] frameworkAtoms(MolGraph g) {
    boolean hasRing = false;
    for (int atom = 0; atom < g.n && !hasRing; atom++) hasRing = g.ring[atom];
    if (!hasRing) return new int[0];
    int[] remaining = new int[g.n], queue = new int[g.n * 2 + 1];
    boolean[] removed = new boolean[g.n];
    int head = 0, tail = 0;
    for (int atom = 0; atom < g.n; atom++) {
      remaining[atom] = g.neighbors[atom].length;
      if (!g.ring[atom] && remaining[atom] <= 1) queue[tail++] = atom;
    }
    while (head < tail) {
      int atom = queue[head++];
      if (removed[atom]) continue;
      removed[atom] = true;
      for (int neighbour : g.neighbors[atom]) {
        if (removed[neighbour]) continue;
        if (--remaining[neighbour] <= 1 && !g.ring[neighbour]) queue[tail++] = neighbour;
      }
    }
    boolean[] keep = new boolean[g.n];
    for (int atom = 0; atom < g.n; atom++) keep[atom] = !removed[atom];
    for (int atom = 0; atom < g.n; atom++) {
      if (removed[atom]) continue;
      for (int neighbour : g.neighbors[atom]) {
        int order = g.bondOrder(atom, neighbour);
        if (removed[neighbour] && (order == 2 || order == 3)) keep[neighbour] = true;
      }
    }
    int kept = 0;
    for (boolean flag : keep) if (flag) kept++;
    int[] atoms = new int[kept];
    for (int atom = 0, at = 0; atom < g.n; atom++) if (keep[atom]) atoms[at++] = atom;
    return atoms;
  }

  /**
   * The framework of a molecule with a ring, capped with hydrogen where side chains were; a
   * molecule without a ring, or one that is its own framework, is returned unchanged.
   */
  static MolGraph frameworkGraph(MolGraph g) {
    int[] kept = frameworkAtoms(g);
    if (kept.length == 0 || kept.length == g.n) return g;
    boolean[] keep = new boolean[g.n];
    for (int atom : kept) keep[atom] = true;
    int[] extraHydrogens = new int[g.n];
    for (int atom : kept)
      for (int neighbour : g.neighbors[atom])
        if (!keep[neighbour]) {
          int order = g.bondOrder(atom, neighbour);
          extraHydrogens[atom] += order >= 1 && order <= 3 ? order : 1;
        }
    return inducedGraph(g, kept, extraHydrogens);
  }

  /** A renumbering-invariant signature of an atom: the classic ECFP atom invariant without its tautomer class. */
  private static long atomSignature(MolGraph g, int atom) {
    long h = SEED;
    h = mix(h, g.atomicNum[atom]);
    h = mix(h, g.degree[atom]);
    int bondOrderSum = 0;
    for (int neighbour : g.neighbors[atom]) bondOrderSum += g.bondOrder(atom, neighbour);
    h = mix(h, bondOrderSum);
    h = mix(h, g.massNumber[atom]);
    h = mix(h, (long) g.formalCharge[atom] + 4L);
    h = mix(h, g.attachedHydrogenCount(atom));
    h = mix(h, g.ring[atom] ? 1 : 0);
    h = mix(h, g.aromatic[atom] ? 1 : 0);
    return h;
  }

  /** A renumbering-invariant signature of a set of atoms: its atom and bond counts and sorted atom signatures. */
  private static long componentSignature(MolGraph g, int[] component, int which, int size) {
    long[] terms = new long[size];
    int ends = 0, filled = 0;
    for (int atom = 0; atom < g.n; atom++) {
      if (component[atom] != which) continue;
      // Sorted as unsigned, as the native core sorts uint64_t.
      terms[filled++] = atomSignature(g, atom) ^ Long.MIN_VALUE;
      ends += g.neighbors[atom].length;
    }
    Arrays.sort(terms);
    long s = SEED;
    s = mix(s, size);
    s = mix(s, ends / 2);
    for (long term : terms) s = mix(s, term ^ Long.MIN_VALUE);
    return s;
  }

  /**
   * The component with the most heavy atoms, then the most atoms; among components alike in both,
   * the one with the greater signature, so that the choice does not depend on the atom order, and
   * among identical components the first in atom order.
   */
  static MolGraph largestComponent(MolGraph g) {
    int[] count = new int[1];
    int[] component = components(g, count);
    if (count[0] <= 1) return g;
    int[] heavy = new int[count[0]], atoms = new int[count[0]];
    for (int atom = 0; atom < g.n; atom++) {
      atoms[component[atom]]++;
      if (g.atomicNum[atom] != 1) heavy[component[atom]]++;
    }
    long[] signature = new long[count[0]];
    int best = 0;
    for (int c = 1; c < count[0]; c++) {
      if (heavy[c] != heavy[best]) { if (heavy[c] > heavy[best]) best = c; continue; }
      if (atoms[c] != atoms[best]) { if (atoms[c] > atoms[best]) best = c; continue; }
      if (signature[best] == 0) signature[best] = componentSignature(g, component, best, atoms[best]);
      signature[c] = componentSignature(g, component, c, atoms[c]);
      if (Long.compareUnsigned(signature[c], signature[best]) > 0) best = c;
    }
    int[] kept = new int[atoms[best]];
    for (int atom = 0, at = 0; atom < g.n; atom++) if (component[atom] == best) kept[at++] = atom;
    return inducedGraph(g, kept, new int[g.n]);
  }

  /** The graph a fingerprint reads: the molecule after its scaffold and fragment handling. */
  static MolGraph subject(MolGraph g, FingerprintOptions options) {
    MolGraph graph = g;
    if (options.scaffold == FingerprintOptions.Scaffold.MURCKO) graph = frameworkGraph(graph);
    if (options.fragments == FingerprintOptions.Fragments.LARGEST) graph = largestComponent(graph);
    return graph;
  }

  // ---- Terms ----

  /** The smallest SSSR ring through each atom, as a class: 0 for an acyclic atom, else the size up to the cap. */
  static int[] ringSizeClasses(MolGraph g, int cap) {
    int[] classes = new int[g.n];
    for (int[] ring : g.computeRings())
      for (int atom : ring)
        if (classes[atom] == 0 || ring.length < classes[atom]) classes[atom] = ring.length;
    for (int atom = 0; atom < g.n; atom++) classes[atom] = Math.min(classes[atom], cap);
    return classes;
  }

  /** The CIP descriptors as codes: R 1, S 2, r 3, s 4 by atom; E 1, Z 2 by double bond, sorted at each end. */
  static final class Stereo {
    final int[] atom;
    final int[][] bondsAt;
    final Map<Long, Integer> bond = new HashMap<>();

    Stereo(MolGraph g, boolean atoms, boolean bonds) {
      atom = new int[g.n];
      bondsAt = new int[g.n][];
      if (atoms)
        for (Map.Entry<Integer, Character> label : CIPAssigner.assignRS(g).entrySet()) {
          char c = label.getValue();
          atom[label.getKey()] = c == 'R' ? 1 : c == 'S' ? 2 : c == 'r' ? 3 : c == 's' ? 4 : 0;
        }
      List<List<Integer>> perAtom = new ArrayList<>(g.n);
      for (int i = 0; i < g.n; i++) perAtom.add(new ArrayList<Integer>());
      if (bonds) {
        // A double bond in a ring of fewer than eight atoms has its geometry fixed by the ring and
        // carries no configuration, whatever marks the input wrote on it.
        Set<Long> smallRingBonds = new HashSet<>();
        for (int[] ring : g.computeRings()) {
          if (ring.length >= 8) continue;
          for (int k = 0; k < ring.length; k++) {
            int a = ring[k], b = ring[(k + 1) % ring.length];
            smallRingBonds.add(MolGraph.bondKey(Math.min(a, b), Math.max(a, b)));
          }
        }
        Map<Long, Character> labels = CIPAssigner.assignEZ(g);
        for (int first = 0; first < g.n; first++)
          for (int second : g.neighbors[first]) {
            if (second <= first || g.bondOrder(first, second) != 2) continue;
            if (smallRingBonds.contains(MolGraph.bondKey(first, second))) continue;
            Character label = labels.get(MolGraph.bondKey(first, second));
            if (label == null) continue;
            int code = label == 'E' ? 1 : 2;
            bond.put(MolGraph.bondKey(first, second), code);
            perAtom.get(first).add(code);
            perAtom.get(second).add(code);
          }
      }
      for (int i = 0; i < g.n; i++) {
        List<Integer> codes = perAtom.get(i);
        bondsAt[i] = new int[codes.size()];
        for (int k = 0; k < codes.size(); k++) bondsAt[i][k] = codes.get(k);
        Arrays.sort(bondsAt[i]);
      }
    }

    int bondCode(int from, int to) {
      Integer code = bond.get(MolGraph.bondKey(Math.min(from, to), Math.max(from, to)));
      return code == null ? 0 : code;
    }
  }

  /** Everything the atom and bond terms of a hashed family read. */
  static final class Terms {
    final MolGraph g;
    final FingerprintOptions o;
    final int[] features;
    final int[] ringSize;
    final Stereo stereo;
    final int[] component;
    final int componentCount;

    Terms(MolGraph graph, FingerprintOptions options) {
      g = graph;
      o = options;
      features = o.family == FingerprintOptions.Family.FCFP ? g.getPharmacophoreFeatures() : null;
      ringSize = o.ringSize ? ringSizeClasses(g, o.ringSizeCap) : null;
      stereo = new Stereo(g, o.cipAtoms, o.cipBonds);
      if (o.fragments == FingerprintOptions.Fragments.PER_FRAGMENT) {
        int[] count = new int[1];
        component = components(g, count);
        componentCount = count[0];
      } else {
        component = null;
        componentCount = 0;
      }
    }

    int bondOrderSum(int atom) {
      int sum = 0;
      for (int neighbour : g.neighbors[atom]) sum += g.bondOrder(atom, neighbour);
      return sum;
    }

    int piElectrons(int atom) {
      int count = 0;
      for (int neighbour : g.neighbors[atom]) {
        int order = g.bondOrder(atom, neighbour);
        if (order == 2) count += 1;
        else if (order == 3) count += 2;
      }
      if (g.aromatic[atom]) count++;
      return count;
    }

    int tautomer(int atom) { return g.tautomerClass != null ? g.tautomerClass[atom] : -1; }

    /** The circular invariant of an atom before the fragment term. */
    long circularBase(int atom) {
      long h = SEED;
      if (features != null) h = mix(h, features[atom]);
      if (o.element) h = mix(h, g.atomicNum[atom]);
      if (o.degree) h = mix(h, g.degree[atom]);
      if (o.valence) h = mix(h, bondOrderSum(atom));
      if (o.isotope) h = mix(h, g.massNumber[atom]);
      if (o.formalCharge) h = mix(h, (long) g.formalCharge[atom] + 4L);
      if (o.hydrogenCount) h = mix(h, g.attachedHydrogenCount(atom));
      if (o.ringMembership) h = mix(h, g.ring[atom] ? 1 : 0);
      if (o.aromaticity) h = mix(h, g.aromatic[atom] ? 1 : 0);
      if (o.tautomerClass && tautomer(atom) >= 0) h = mix(h, 0xCAFE00L + tautomer(atom));
      if (o.ringSize) h = mix(h, ringSize[atom]);
      return h;
    }

    /** The path atom term before the fragment term: the atomic number alone in the classic form. */
    int pathBase(int atom) {
      int t = o.element ? g.atomicNum[atom] : 0;
      if (o.degree) t = t * 37 + g.degree[atom];
      if (o.valence) t = t * 37 + bondOrderSum(atom);
      if (o.isotope) t = t * 37 + g.massNumber[atom];
      if (o.formalCharge) t = t * 37 + (g.formalCharge[atom] + 4);
      if (o.hydrogenCount) t = t * 37 + (int) g.attachedHydrogenCount(atom);
      if (o.ringMembership) t = t * 37 + (g.ring[atom] ? 1 : 0);
      if (o.aromaticity) t = t * 37 + (g.aromatic[atom] ? 1 : 0);
      if (o.tautomerClass) t = t * 37 + (tautomer(atom) + 1);
      if (o.ringSize) t = t * 37 + ringSize[atom];
      if (o.cipAtoms) t = t * 37 + stereo.atom[atom];
      return t;
    }

    boolean pathHashesBonds() { return o.bondOrder || o.aromaticBond || o.ringBond || o.cipBonds; }

    int pathBond(int from, int to) {
      int b = 17;
      if (o.bondOrder) b = b * 37 + g.bondOrder(from, to);
      if (o.ringBond) b = b * 37 + (g.bondInRing(from, to) ? 1 : 0);
      if (o.aromaticBond) b = b * 37 + (g.bondAromatic(from, to) ? 1 : 0);
      if (o.cipBonds) b = b * 37 + stereo.bondCode(from, to);
      return b;
    }

    /** The torsion atom type before the fragment term: the classic packed type, extended by mixing. */
    long torsionBase(int atom) {
      int z = o.element ? g.atomicNum[atom] : 0;
      int deg = o.degree ? g.degree[atom] : 0;
      int pi = o.valence ? piElectrons(atom) : 0;
      int inRing = o.ringMembership && g.ring[atom] ? 1 : 0;
      long t = (z << 16) | (deg << 8) | (pi << 4) | inRing;
      if (o.isotope) t = mix(t, g.massNumber[atom]);
      if (o.formalCharge) t = mix(t, (long) g.formalCharge[atom] + 4L);
      if (o.hydrogenCount) t = mix(t, g.attachedHydrogenCount(atom));
      if (o.aromaticity) t = mix(t, g.aromatic[atom] ? 1 : 0);
      if (o.tautomerClass) t = mix(t, tautomer(atom) + 1);
      if (o.ringSize) t = mix(t, ringSize[atom]);
      if (o.cipAtoms && stereo.atom[atom] != 0) t = mix(t, CIP_ATOM_TAG + stereo.atom[atom]);
      if (o.cipBonds) for (int code : stereo.bondsAt[atom]) t = mix(t, CIP_BOND_TAG + code);
      return t;
    }

    /** The class of each component: a hash of its atom and bond counts and its sorted base terms. */
    long[] fragmentTerms(long[] base) {
      int[] sizes = new int[componentCount], bonds = new int[componentCount];
      for (int atom = 0; atom < g.n; atom++) {
        sizes[component[atom]]++;
        bonds[component[atom]] += g.neighbors[atom].length;
      }
      long[][] perComponent = new long[componentCount][];
      for (int c = 0; c < componentCount; c++) perComponent[c] = new long[sizes[c]];
      int[] filled = new int[componentCount];
      for (int atom = 0; atom < g.n; atom++) {
        int c = component[atom];
        // Sorted as unsigned, as the native core sorts uint64_t.
        perComponent[c][filled[c]++] = base[atom] ^ Long.MIN_VALUE;
      }
      long[] signature = new long[componentCount];
      for (int c = 0; c < componentCount; c++) {
        Arrays.sort(perComponent[c]);
        long s = SEED;
        s = mix(s, sizes[c]);
        s = mix(s, bonds[c] / 2);
        for (long term : perComponent[c]) s = mix(s, term ^ Long.MIN_VALUE);
        signature[c] = s;
      }
      long[] result = new long[g.n];
      for (int atom = 0; atom < g.n; atom++) result[atom] = signature[component[atom]];
      return result;
    }
  }

  // ---- Families ----

  /** Receives every hash of a hashed family. */
  interface Sink { void record(long hash); }

  /** Hands every hash of a circular fingerprint to {@code sink}. */
  static void circularHashes(Terms terms, Sink sink) {
    MolGraph g = terms.g;
    FingerprintOptions o = terms.o;
    int n = g.n;
    long[] prevHash = new long[n];
    for (int atom = 0; atom < n; atom++) prevHash[atom] = terms.circularBase(atom);
    if (o.fragments == FingerprintOptions.Fragments.PER_FRAGMENT) {
      long[] fragment = terms.fragmentTerms(prevHash);
      for (int atom = 0; atom < n; atom++) prevHash[atom] = mix(prevHash[atom], fragment[atom]);
    }
    for (int atom = 0; atom < n; atom++) sink.record(prevHash[atom]);
    int maxRadius = o.radius < 0 ? maxComponentExtent(g) : o.radius;
    long[] nextHash = new long[n];
    int maxDegree = 0;
    for (int[] row : g.neighbors) maxDegree = Math.max(maxDegree, row.length);
    long[] nbHashes = new long[maxDegree];
    for (int r = 1; r <= maxRadius; r++) {
      for (int atom = 0; atom < n; atom++) {
        int deg = g.neighbors[atom].length;
        if (deg == 0) { nextHash[atom] = prevHash[atom]; continue; }
        // The neighbour hashes sort as unsigned, as the native core sorts uint64_t: the sign bit
        // is flipped for the signed sort and flipped back when the hash is read.
        for (int k = 0; k < deg; k++) {
          int neighbour = g.neighbors[atom][k];
          long term = prevHash[neighbour];
          if (o.bondOrder) term ^= g.bondOrder(atom, neighbour) * 31L;
          if (o.aromaticBond && g.bondAromatic(atom, neighbour)) term ^= AROMATIC_BOND_MARK;
          if (o.ringBond && g.bondInRing(atom, neighbour)) term ^= RING_BOND_MARK;
          nbHashes[k] = term ^ Long.MIN_VALUE;
        }
        Arrays.sort(nbHashes, 0, deg);
        long h = SEED;
        h = mix(h, r);
        h = mix(h, prevHash[atom]);
        for (int k = 0; k < deg; k++) h = mix(h, nbHashes[k] ^ Long.MIN_VALUE);
        if (r == 1) {
          // The descriptors join after the first iteration, so a centre changes the bits of its
          // own environments from radius 1 and of its neighbours' from radius 2.
          if (o.cipAtoms && terms.stereo.atom[atom] != 0) h = mix(h, CIP_ATOM_TAG + terms.stereo.atom[atom]);
          if (o.cipBonds) for (int code : terms.stereo.bondsAt[atom]) h = mix(h, CIP_BOND_TAG + code);
        }
        nextHash[atom] = h;
        sink.record(h);
      }
      System.arraycopy(nextHash, 0, prevHash, 0, n);
    }
  }

  /** Largest finite shortest-path distance, including disconnected components. */
  private static int maxComponentExtent(MolGraph graph) {
    int extent = 0;
    int[] distances = new int[graph.n], queue = new int[graph.n];
    for (int start = 0; start < graph.n; start++) {
      Arrays.fill(distances, -1);
      distances[start] = 0;
      int next = 0, count = 1;
      queue[0] = start;
      while (next < count) {
        int atom = queue[next++];
        for (int neighbour : graph.neighbors[atom]) {
          if (distances[neighbour] >= 0) continue;
          distances[neighbour] = distances[atom] + 1;
          extent = Math.max(extent, distances[neighbour]);
          queue[count++] = neighbour;
        }
      }
    }
    return extent;
  }

  /** The path fingerprint of the terms, over heavy atoms, as the classic element path fingerprint. */
  static long[] pathBits(Terms terms, int fpSize) {
    MolGraph g = terms.g;
    int[] atomTerm = new int[g.n];
    for (int atom = 0; atom < g.n; atom++) atomTerm[atom] = terms.pathBase(atom);
    if (terms.o.fragments == FingerprintOptions.Fragments.PER_FRAGMENT) {
      long[] base = new long[g.n];
      for (int atom = 0; atom < g.n; atom++) base[atom] = atomTerm[atom] & 0xFFFFFFFFL;
      long[] fragment = terms.fragmentTerms(base);
      for (int atom = 0; atom < g.n; atom++)
        atomTerm[atom] = atomTerm[atom] * 37 + (int) (fragment[atom] ^ (fragment[atom] >>> 32));
    }
    long[] fp = new long[(fpSize - 1) / 64 + 1];
    boolean bonds = terms.pathHashesBonds();
    boolean[] visited = new boolean[g.n];
    int[] path = new int[Math.min(terms.o.pathLength, g.n) + 1];
    for (int start = 0; start < g.n; start++) {
      if (g.atomicNum[start] == 1) continue;
      path[0] = start;
      visited[start] = true;
      setBit(fp, pathHash(terms, atomTerm, path, 1, bonds) % fpSize);
      extendPath(terms, atomTerm, path, 1, terms.o.pathLength, bonds, visited, fp, fpSize);
      visited[start] = false;
    }
    return fp;
  }

  /** Sets a bit for every simple path of heavy atoms extending the current one. */
  private static void extendPath(Terms terms, int[] atomTerm, int[] path, int length, int maxBonds, boolean bonds,
      boolean[] visited, long[] fp, int fpSize) {
    MolGraph g = terms.g;
    for (int next : g.neighbors[path[length - 1]]) {
      if (visited[next] || g.atomicNum[next] == 1) continue;
      path[length] = next;
      setBit(fp, pathHash(terms, atomTerm, path, length + 1, bonds) % fpSize);
      if (length < maxBonds) {
        visited[next] = true;
        extendPath(terms, atomTerm, path, length + 1, maxBonds, bonds, visited, fp, fpSize);
        visited[next] = false;
      }
    }
  }

  /** The direction-independent hash of the first {@code length} atoms of a path, and their bonds when hashed. */
  private static int pathHash(Terms terms, int[] atomTerm, int[] path, int length, boolean bonds) {
    int fwd = 17, rev = 17;
    for (int i = 0; i < length; i++) {
      int fi = i, ri = length - 1 - i;
      fwd = fwd * 31 + atomTerm[path[fi]];
      rev = rev * 31 + atomTerm[path[ri]];
      if (bonds && i < length - 1) {
        fwd = fwd * 31 + terms.pathBond(path[fi], path[fi + 1]);
        rev = rev * 31 + terms.pathBond(path[ri], path[ri - 1]);
      }
    }
    return Math.min(fwd & 0x7FFFFFFF, rev & 0x7FFFFFFF);
  }

  private static void setBit(long[] fp, int bit) { fp[bit >>> 6] |= 1L << (bit & 63); }

  /** Hands the canonical hash of every four-atom path to {@code sink}. */
  static void torsionHashes(Terms terms, Sink sink) {
    MolGraph g = terms.g;
    int n = g.n;
    long[] atomType = new long[n];
    for (int atom = 0; atom < n; atom++) atomType[atom] = terms.torsionBase(atom);
    if (terms.o.fragments == FingerprintOptions.Fragments.PER_FRAGMENT) {
      long[] fragment = terms.fragmentTerms(atomType);
      for (int atom = 0; atom < n; atom++) atomType[atom] = mix(atomType[atom], fragment[atom]);
    }
    if (n < 4) return;
    for (int a = 0; a < n; a++)
      for (int b : g.neighbors[a])
        for (int c : g.neighbors[b]) {
          if (c == a) continue;
          for (int d : g.neighbors[c]) {
            if (d == b || d <= a) continue;
            long fwd = SEED;
            fwd = mix(fwd, atomType[a]);
            fwd = mix(fwd, atomType[b]);
            fwd = mix(fwd, atomType[c]);
            fwd = mix(fwd, atomType[d]);
            long rev = SEED;
            rev = mix(rev, atomType[d]);
            rev = mix(rev, atomType[c]);
            rev = mix(rev, atomType[b]);
            rev = mix(rev, atomType[a]);
            sink.record(Long.compareUnsigned(fwd, rev) <= 0 ? fwd : rev);
          }
        }
  }

  // ---- Substructure keys ----

  /** Prepared SMARTS searches by pattern, shared by every key fingerprint. */
  private static final ConcurrentHashMap<String, Substructure> PREPARED = new ConcurrentHashMap<>();

  private static Substructure preparedKey(String pattern) {
    Substructure prepared = PREPARED.get(pattern);
    if (prepared != null) return prepared;
    Substructure.Options options = new Substructure.Options();
    options.timeoutMillis = KEY_SEARCH_MILLIS;
    prepared = Substructure.smarts(pattern).withOptions(options);
    if (PREPARED.size() >= 4096) PREPARED.clear();
    PREPARED.putIfAbsent(pattern, prepared);
    return prepared;
  }

  /** A key search must run to its end: a search cut short by time would make the bits depend on the machine. */
  private static void requireSettled(SearchResult result) {
    if (result.status() == SearchResult.Status.TimedOut || result.status() == SearchResult.Status.Cancelled)
      throw new IllegalStateException("a substructure key search did not complete within its time budget");
  }

  /**
   * How many distinct atom sets a SMARTS pattern matches, or whether it matches at all when
   * {@code atMost} is one. Beyond the embedding limit the count is settled only when it has reached
   * {@code atMost}; otherwise the key is reported as an error rather than guessed.
   */
  private static int keyMatches(Substructure query, MolGraph g, int atMost) {
    if (atMost <= 1) {
      SearchResult first = query.match(g);
      requireSettled(first);
      return first.found() ? 1 : 0;
    }
    SearchResult all = query.matchAll(g, KEY_MATCH_LIMIT);
    requireSettled(all);
    Set<List<Integer>> distinct = new HashSet<>();
    for (SearchResult.Match match : all.matches()) {
      List<Integer> atoms = new ArrayList<>(match.atomMapping().values());
      java.util.Collections.sort(atoms);
      distinct.add(atoms);
    }
    if (all.status() == SearchResult.Status.LimitReached && distinct.size() < atMost)
      throw new IllegalStateException("a substructure key matched more than 100000 embeddings before its count was settled");
    return distinct.size();
  }

  private static int aromaticRingCount(MolGraph g) {
    int count = 0;
    for (int[] ring : g.computeRings()) {
      boolean aromatic = ring.length > 0;
      for (int k = 0; k < ring.length && aromatic; k++) aromatic = g.bondAromatic(ring[k], ring[(k + 1) % ring.length]);
      if (aromatic) count++;
    }
    return count;
  }

  /** The count of each key; {@code exact} asks for the full count rather than whether the minimum is met. */
  static int[] keyCounts(MolGraph g, FingerprintOptions options, boolean exact) {
    int[] counts = new int[options.keys.size()];
    int componentCount = -1;
    for (int k = 0; k < counts.length; k++) {
      FingerprintOptions.Key key = options.keys.get(k);
      switch (key.kind) {
        case SMARTS:
          counts[k] = keyMatches(preparedKey(key.smarts), g, exact ? Integer.MAX_VALUE : key.minimum);
          break;
        case ISOTOPE:
          for (int atom = 0; atom < g.n; atom++) if (g.massNumber[atom] != 0) counts[k]++;
          break;
        case AROMATIC_RINGS:
          counts[k] = aromaticRingCount(g);
          break;
        case FRAGMENTS:
          if (componentCount < 0) {
            int[] count = new int[1];
            components(g, count);
            componentCount = count[0];
          }
          counts[k] = componentCount;
          break;
        default:
          break;
      }
    }
    return counts;
  }

  // ---- Entry points ----

  /** The bit fingerprint of a molecule under {@code options}. */
  static long[] bits(MolGraph molecule, FingerprintOptions options) {
    validate(options);
    if (molecule == null) throw new NullPointerException("MolGraph must not be null");
    MolGraph g = subject(molecule, options);
    if (options.family == FingerprintOptions.Family.KEYS) {
      int size = options.keys.size();
      long[] fp = new long[(size - 1) / 64 + 1];
      int[] counts = keyCounts(g, options, false);
      for (int k = 0; k < size; k++) if (counts[k] >= options.keys.get(k).minimum) setBit(fp, k);
      return fp;
    }
    final int fpSize = options.size;
    final long[] fp = new long[(fpSize - 1) / 64 + 1];
    if (g.n == 0) return fp;
    Terms terms = new Terms(g, options);
    Sink setBit = new Sink() {
      @Override public void record(long hash) { setBit(fp, (int) Long.remainderUnsigned(hash, fpSize)); }
    };
    switch (options.family) {
      case ECFP:
      case FCFP:
        circularHashes(terms, setBit);
        return fp;
      case PATH:
        return pathBits(terms, fpSize);
      case TORSION:
        torsionHashes(terms, setBit);
        return fp;
      default:
        return fp;
    }
  }

  /** The count fingerprint of a molecule under {@code options}: one bin per bit position, or per key. */
  static int[] counts(MolGraph molecule, FingerprintOptions options) {
    validate(options);
    if (molecule == null) throw new NullPointerException("MolGraph must not be null");
    MolGraph g = subject(molecule, options);
    if (options.family == FingerprintOptions.Family.KEYS) return keyCounts(g, options, true);
    final int fpSize = options.size;
    final int[] counts = new int[fpSize];
    if (g.n == 0) return counts;
    Terms terms = new Terms(g, options);
    Sink count = new Sink() {
      @Override public void record(long hash) { counts[(int) Long.remainderUnsigned(hash, fpSize)]++; }
    };
    switch (options.family) {
      case ECFP:
      case FCFP:
        circularHashes(terms, count);
        break;
      case PATH: {
        // Path bits are set once per path, so the count form counts the set positions.
        long[] bits = pathBits(terms, fpSize);
        for (int bit = 0; bit < fpSize; bit++) if ((bits[bit >>> 6] & (1L << (bit & 63))) != 0) counts[bit] = 1;
        break;
      }
      case TORSION:
        torsionHashes(terms, count);
        break;
      default:
        break;
    }
    return counts;
  }
}
