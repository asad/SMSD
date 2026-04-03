/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms. */
package com.bioinception.smsd.core;

import com.fasterxml.jackson.databind.ObjectMapper;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.aromaticity.ElectronDonation;
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
    CDKHydrogenAdder.getInstance(SilentChemObjectBuilder.getInstance()).addImplicitHydrogens(m);
    Cycles.markRingAtomsAndBonds(m);
    @SuppressWarnings("deprecation")
    Aromaticity arom =
        new Aromaticity(ElectronDonation.daylight(), Cycles.all());
    arom.apply(m);

    return m;
  }

  // ========== SMARTS matching ==========

  /** Match a SMARTS pattern against a target molecule, returning unique atom mappings. */
  public static List<int[]> matchAll(String smarts, IAtomContainer target) {
    return new ArrayList<>(Arrays.asList(
        SmartsPattern.create(smarts).matchAll(target).uniqueAtoms().toArray()));
  }

  // ========== Named predicate expansion ==========

  /** Expand $name tokens to recursive SMARTS using a registry, leaving existing $() intact. */
  public static String expandNamedPredicates(String smarts, PredicateRegistry reg) {
    StringBuilder out = new StringBuilder();
    int i = 0, n = smarts.length();
    while (i < n) {
      if (smarts.charAt(i) == '[') {
        int depth = 1, j = i + 1;
        while (j < n && depth > 0) {
          char d = smarts.charAt(j);
          if (d == '[') depth++;
          else if (d == ']') depth--;
          j++;
        }
        out.append('[').append(rewriteContent(smarts.substring(i + 1, j - 1), reg)).append(']');
        i = j;
      } else {
        out.append(smarts.charAt(i++));
      }
    }
    return out.toString();
  }

  private static String rewriteContent(String content, PredicateRegistry reg) {
    StringBuilder sb = new StringBuilder();
    int i = 0, n = content.length();
    while (i < n) {
      char c = content.charAt(i);
      if (c != '$') { sb.append(c); i++; continue; }
      if (i + 1 < n && content.charAt(i + 1) == '(') {
        int depth = 1, j = i + 2;
        while (j < n && depth > 0) {
          char d = content.charAt(j);
          if (d == '(') depth++;
          else if (d == ')') depth--;
          j++;
        }
        sb.append(content, i, j);
        i = j;
      } else {
        int j = i + 1;
        while (j < n && (Character.isLetterOrDigit(content.charAt(j)) || content.charAt(j) == '_'))
          j++;
        String name = content.substring(i + 1, j);
        if (!name.isEmpty() && reg.contains(name)) {
          sb.append("$(").append(reg.get(name)).append(")");
          i = j;
        } else {
          sb.append(c);
          i++;
        }
      }
    }
    return sb.toString();
  }

  // ========== Predicate registry ==========

  /** Named SMARTS predicate registry: $name -> SMARTS. Supports loading from JSON files. */
  public static class PredicateRegistry {

    private static final ObjectMapper MAPPER = new ObjectMapper();
    private final Map<String, String> map = new LinkedHashMap<>();

    public PredicateRegistry() {
      map.put("isAmideN", "[$([NX3;H2,H1;!$(NC=O)]),$([NX3;$(NC=O)])]");
      map.put("isCarboxylC", "[CX3](=O)[OX2H1,OX1-]");
      map.put("isSulfonamideN", "[NX3;$(NS(=O)=O)]");
      map.put("isEster", "[#6][CX3](=O)[OX2][#6]");
      map.put("isKetone", "[#6][CX3](=O)[#6]");
      map.put("isHalogen", "[F,Cl,Br,I]");
      map.put("isPositive", "[*+]");
      map.put("isNegative", "[*-]");
      map.put("isHBDonor", "[N!H0,O!H0,S!H0]");
      map.put("isHBAcceptor", "[N,O,S;H0;!+0]");
    }

    public void put(String name, String smarts) { map.put(name, smarts); }
    public boolean contains(String name) { return map.containsKey(name); }
    public String get(String name) { return map.get(name); }

    /** Load predicates from a JSON file: {"name": "SMARTS", ...}. Non-String values are skipped. */
    @SuppressWarnings("unchecked")
    public void loadJsonFile(String path) throws Exception {
      Map<String, Object> raw = MAPPER.readValue(new File(path), Map.class);
      for (Map.Entry<String, Object> entry : raw.entrySet()) {
        if (entry.getValue() instanceof String) {
          map.put(entry.getKey(), (String) entry.getValue());
        }
        // silently skip non-String values
      }
    }
  }
}
