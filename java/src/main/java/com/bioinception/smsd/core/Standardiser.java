/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
package com.bioinception.smsd.core;

import com.fasterxml.jackson.databind.ObjectMapper;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smarts.SmartsPattern;
import org.openscience.cdk.atomtype.CDKAtomTypeMatcher;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomType;
import org.openscience.cdk.tools.CDKHydrogenAdder;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import org.openscience.cdk.tools.manipulator.AtomTypeManipulator;

/**
 * One-stop preprocessing: molecule standardisation, SMARTS matching, named predicate expansion,
 * and predicate registry.
 *
 * @author Syed Asad Rahman
 */
public final class Standardiser {

  private Standardiser() {}

  /** Tautomer handling mode. */
  public enum TautomerMode {
    NONE
  }

  /**
   * Standardise a molecule (clone, atom-type, add implicit H, ring perception, aromaticity).
   * The original is never modified.
   */
  public static IAtomContainer standardise(IAtomContainer mol, TautomerMode mode) throws Exception {
    if (mol == null) throw new NullPointerException("Input molecule must not be null");
    if (mode == null) mode = TautomerMode.NONE;

    IAtomContainer m = mol.clone();
    // Explicit per-atom type matching for robust handling of exotic valence states
    // (organometallics, hypervalent S/P, charged aromatics).  The convenience
    // shortcut percieveAtomTypesAndConfigureAtoms silently skips untyped atoms;
    // explicit iteration lets us log them and still proceed safely.
    CDKAtomTypeMatcher matcher = CDKAtomTypeMatcher.getInstance(m.getBuilder());
    for (IAtom atom : m.atoms()) {
      IAtomType matched = matcher.findMatchingAtomType(m, atom);
      if (matched != null) {
        AtomTypeManipulator.configure(atom, matched);
      }
      // Unmatched atoms retain their original properties — hydrogen adder
      // will simply skip them, which is safer than the shortcut's silent failure.
    }
    CDKHydrogenAdder hydrogens = CDKHydrogenAdder.getInstance(SilentChemObjectBuilder.getInstance());
    for (IAtom atom : m.atoms())
      if (atom.getImplicitHydrogenCount() == null && atom.getAtomTypeName() != null)
        hydrogens.addImplicitHydrogens(m, atom);
    RingPerception.markRingAtomsAndBonds(m);
    MolGraph.kekuliseAromaticBondSystems(m);
    Aromaticity arom =
        new Aromaticity(Aromaticity.Model.Daylight, Cycles.all());
    arom.apply(m);

    return m;
  }

  // ========== SMARTS matching ==========

  /**
   * Match a SMARTS pattern against a target molecule, returning unique atom mappings.
   *
   * @deprecated Use {@code Substructure.smarts(smarts).matchAll(target)}, which applies the SMARTS
   *     limits and the Daylight semantics of the native core and reports a {@link SearchResult}.
   */
  @Deprecated
  public static List<int[]> matchAll(String smarts, IAtomContainer target) {
    return new ArrayList<>(Arrays.asList(
        SmartsPattern.create(smarts).matchAll(target).uniqueAtoms().toArray()));
  }

  // ========== Named predicate expansion ==========

  /**
   * Replace each $name token of a registered name, including those inside recursive SMARTS, with
   * $(definition), as the native core does. A name is ASCII letters, digits and underscores;
   * definitions are inserted as written and $(...) is kept.
   */
  public static String expandNamedPredicates(String smarts, PredicateRegistry reg) {
    StringBuilder out = new StringBuilder(smarts.length());
    int i = 0, n = smarts.length();
    while (i < n) {
      char c = smarts.charAt(i);
      if (c == '$' && i + 1 < n && smarts.charAt(i + 1) != '(') {
        int j = i + 1;
        while (j < n && nameCharacter(smarts.charAt(j))) j++;
        String name = smarts.substring(i + 1, j);
        if (reg.contains(name)) {
          out.append("$(").append(reg.get(name)).append(')');
          i = j;
          continue;
        }
      }
      out.append(c);
      i++;
    }
    return out.toString();
  }

  private static boolean nameCharacter(char c) {
    return c >= 'a' && c <= 'z' || c >= 'A' && c <= 'Z' || c >= '0' && c <= '9' || c == '_';
  }

  // ========== Predicate registry ==========

  /**
   * Named SMARTS predicate registry: $name -> SMARTS. A definition is expanded as recursive SMARTS,
   * so it describes the atom written first. Supports loading from JSON files with
   * {@link #loadJsonFile}, which needs the optional dependency
   * {@code com.fasterxml.jackson.core:jackson-databind}.
   */
  public static class PredicateRegistry {

    /**
     * Holds the JSON reader, so Jackson is loaded only when {@link #loadJsonFile} is first called
     * and the registry otherwise works without it on the class path.
     */
    private static final class Json {
      static final ObjectMapper MAPPER = new ObjectMapper();
    }

    private final Map<String, String> map = new LinkedHashMap<>();

    /** Starts with the ten built-in predicates, defined as in the native core. */
    public PredicateRegistry() {
      // N bonded to a carbonyl carbon: amides, lactams, imides, ureas and carbamates, those in
      // aromatic rings included.
      map.put("isAmideN", "[#7X3][#6]=O");
      // Carbonyl carbon bearing OH or O-: carboxylic, carbamic and carbonic acids and their anions.
      map.put("isCarboxylC", "[CX3](=O)[OX2H1,OX1-]");
      map.put("isSulfonamideN", "[NX3;$(NS(=O)=O)]");
      // Carbonyl carbon of a carboxylic ester, formates and aromatic lactones included; anhydrides
      // are excluded.
      map.put("isEster", "[#6X3;H1,$([#6][#6])](=O)[#8X2][#6;!$([#6]=O)]");
      // Carbonyl carbon bonded to two carbons, carbonyls in aromatic rings included.
      map.put("isKetone", "[#6X3](=O)([#6])[#6]");
      map.put("isHalogen", "[F,Cl,Br,I]");
      // Every positive or negative charge in the OpenSMILES range, -15 to +15.
      map.put("isPositive", "[+,+2,+3,+4,+5,+6,+7,+8,+9,+10,+11,+12,+13,+14,+15]");
      map.put("isNegative", "[-,-2,-3,-4,-5,-6,-7,-8,-9,-10,-11,-12,-13,-14,-15]");
      // N, O or S with hydrogen, aromatic NH included.
      map.put("isHBDonor", "[#7,#8,#16;!H0]");
      // N, O or S of charge 0 to -3 with a free lone pair. Amide N, aromatic N with three
      // connections (pyrrole type), N of valence 5 (uncharged nitro and N-oxide) and S of valence
      // 4 or 6 (sulfoxides, sulfones) are excluded; aniline, nitrile, thioamide, amidine and
      // sulfonamide N are acceptors.
      map.put("isHBAcceptor", "[#7&!X3&!v5,#7&A&!v5&!$(NC=O),#8,#16&!v4&!v6;+0,-,-2,-3]");
    }

    public void put(String name, String smarts) { map.put(name, smarts); }
    public boolean contains(String name) { return map.containsKey(name); }
    public String get(String name) { return map.get(name); }

    /**
     * Load predicates from a JSON file: {"name": "SMARTS", ...}. Non-String values are skipped.
     *
     * <p>Requires {@code com.fasterxml.jackson.core:jackson-databind}, an optional dependency of
     * SMSD that applications using this method add themselves.
     *
     * @throws IllegalStateException if jackson-databind is not on the class path
     */
    @SuppressWarnings("unchecked")
    public void loadJsonFile(String path) throws Exception {
      Map<String, Object> raw;
      try {
        raw = Json.MAPPER.readValue(new File(path), Map.class);
      } catch (NoClassDefFoundError missing) {
        throw new IllegalStateException("loadJsonFile requires "
            + "com.fasterxml.jackson.core:jackson-databind on the class path", missing);
      }
      for (Map.Entry<String, Object> entry : raw.entrySet()) {
        if (entry.getValue() instanceof String) {
          map.put(entry.getKey(), (String) entry.getValue());
        }
        // silently skip non-String values
      }
    }
  }
}
