/*
 * SPDX-License-Identifier: Apache-2.0 AND BSD-3-Clause
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
package com.bioinception.smsd.core;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * A fingerprint specification shared by the circular (ECFP and FCFP), path, topological torsion
 * and substructure-key families, with the same bits in C++ and Python.
 *
 * <p>A new instance is the classic ECFP configuration; the factories give the classic
 * configuration of the other families, which reproduce the bits of the fixed
 * {@link FingerprintEngine} methods exactly. Every further switch extends the hashed terms
 * without moving the classic ones. Pass the options to
 * {@link FingerprintEngine#fingerprint(MolGraph, FingerprintOptions)} or
 * {@link FingerprintEngine#fingerprintCounts(MolGraph, FingerprintOptions)}.
 */
public final class FingerprintOptions {

  /** The hashed families and the substructure-key family. */
  public enum Family { ECFP, FCFP, PATH, TORSION, KEYS }

  /** How the connected components of a molecule are fingerprinted. */
  public enum Fragments {
    /** The whole molecule, as the fixed methods do. */
    WHOLE,
    /** The component with the most heavy atoms only. */
    LARGEST,
    /** Every atom term carries the class of its component, so bits are component-specific. */
    PER_FRAGMENT
  }

  /** Whether the molecule is reduced before it is fingerprinted. */
  public enum Scaffold {
    /** The molecule as given. */
    NONE,
    /**
     * Its Bemis-Murcko scaffold, as {@link MCS#murckoScaffold} keeps it: ring atoms, the linkers
     * between rings and the atoms double or triple bonded to either, capped with hydrogen.
     */
    MURCKO
  }

  /**
   * One substructure key: its bit is set when the molecule holds at least {@code minimum} of what it
   * counts, which for a SMARTS pattern is its distinct matched atom sets.
   */
  public static final class Key {
    /** What a key counts. */
    public enum Kind {
      /** Distinct atom sets matched by a SMARTS pattern. */
      SMARTS,
      /** Atoms with a mass number. */
      ISOTOPE,
      /** Rings of the smallest set of smallest rings whose bonds are all aromatic. */
      AROMATIC_RINGS,
      /** Connected components. */
      FRAGMENTS
    }

    public final Kind kind;
    public final String smarts;
    public final int minimum;

    public Key(String smarts) { this(smarts, 1); }

    public Key(String smarts, int minimum) { this(Kind.SMARTS, smarts, minimum); }

    private Key(Kind kind, String smarts, int minimum) {
      this.kind = kind;
      this.smarts = smarts == null ? "" : smarts;
      this.minimum = minimum;
    }

    /** A key counting something other than a SMARTS pattern. */
    public static Key special(Kind kind, int minimum) { return new Key(kind, "", minimum); }

    @Override public String toString() {
      return kind == Kind.SMARTS ? smarts + " >= " + minimum : kind + " >= " + minimum;
    }
  }

  public Family family = Family.ECFP;
  /** Circular families: the largest radius; -1 iterates to the largest component extent. */
  public int radius = 2;
  /** Path family: the largest number of bonds in a path. */
  public int pathLength = 7;
  /** Bits of a bit fingerprint or bins of a count fingerprint; the key family uses one per key. */
  public int size = 2048;
  /** Count fingerprint rather than bits, where one entry point serves both (Python). */
  public boolean counts;

  // Atom terms.
  /** Atomic number. */
  public boolean element = true;
  /** Number of explicit neighbours (heavy atoms of a prepared graph). */
  public boolean degree = true;
  /** Bond order sum; the torsion family reads its pi-electron count here. */
  public boolean valence = true;
  public boolean formalCharge = true;
  /** Mass number. */
  public boolean isotope = true;
  /** Attached hydrogens, implicit and explicit. */
  public boolean hydrogenCount = true;
  /** The tautomer class, when one is assigned. */
  public boolean tautomerClass = true;
  public boolean aromaticity = true;

  // Bond terms.
  /** Bond order; aromatic bonds have order 4. */
  public boolean bondOrder = true;
  public boolean aromaticBond;
  public boolean ringBond;

  // Ring terms.
  public boolean ringMembership = true;
  /** Smallest ring through the atom, from the SSSR, as a class up to {@link #ringSizeCap}. */
  public boolean ringSize;
  /** Rings of this size or more share one class. */
  public int ringSizeCap = 8;

  // Stereo terms, from the library's own CIP assignment.
  /** R/S (and r/s) of tetrahedral centres. */
  public boolean cipAtoms;
  /** E/Z of stereogenic double bonds. */
  public boolean cipBonds;

  public Fragments fragments = Fragments.WHOLE;
  public Scaffold scaffold = Scaffold.NONE;

  /** The keys of the key family, in bit order. */
  public List<Key> keys = new ArrayList<>();

  /** The classic ECFP configuration. */
  public FingerprintOptions() {}

  /** The classic configuration of a family. */
  public static FingerprintOptions of(Family family) {
    FingerprintOptions options = new FingerprintOptions();
    options.family = family;
    switch (family) {
      case ECFP:
        break;
      case FCFP:
        options.element = options.valence = options.formalCharge = options.isotope = false;
        options.hydrogenCount = options.tautomerClass = options.aromaticity = false;
        break;
      case PATH:
        options.degree = options.valence = options.formalCharge = options.isotope = false;
        options.hydrogenCount = options.tautomerClass = options.aromaticity = false;
        options.ringMembership = false;
        options.bondOrder = false;
        break;
      case TORSION:
        options.formalCharge = options.isotope = options.hydrogenCount = false;
        options.tautomerClass = options.aromaticity = false;
        options.bondOrder = false;
        break;
      case KEYS:
        options.size = 0;
        break;
      default:
        break;
    }
    return options;
  }

  public static FingerprintOptions ecfp() { return ecfp(2, 2048); }

  public static FingerprintOptions ecfp(int radius, int size) {
    FingerprintOptions options = of(Family.ECFP);
    options.radius = radius;
    options.size = size;
    return options;
  }

  public static FingerprintOptions fcfp() { return fcfp(2, 2048); }

  public static FingerprintOptions fcfp(int radius, int size) {
    FingerprintOptions options = of(Family.FCFP);
    options.radius = radius;
    options.size = size;
    return options;
  }

  public static FingerprintOptions path() { return path(7, 2048); }

  public static FingerprintOptions path(int pathLength, int size) {
    FingerprintOptions options = of(Family.PATH);
    options.pathLength = pathLength;
    options.size = size;
    return options;
  }

  public static FingerprintOptions torsion() { return torsion(2048); }

  public static FingerprintOptions torsion(int size) {
    FingerprintOptions options = of(Family.TORSION);
    options.size = size;
    return options;
  }

  /** A key fingerprint over the given keys, one bit each. */
  public static FingerprintOptions withKeys(List<Key> keys) {
    FingerprintOptions options = of(Family.KEYS);
    options.keys = new ArrayList<>(keys);
    options.size = options.keys.size();
    return options;
  }

  /** The 166 MACCS structural keys; key k is bit k - 1. */
  public static FingerprintOptions maccs() { return withKeys(maccsKeys()); }

  /** A copy of these options. */
  public FingerprintOptions copy() {
    FingerprintOptions options = new FingerprintOptions();
    options.family = family;
    options.radius = radius;
    options.pathLength = pathLength;
    options.size = size;
    options.counts = counts;
    options.element = element;
    options.degree = degree;
    options.valence = valence;
    options.formalCharge = formalCharge;
    options.isotope = isotope;
    options.hydrogenCount = hydrogenCount;
    options.tautomerClass = tautomerClass;
    options.aromaticity = aromaticity;
    options.bondOrder = bondOrder;
    options.aromaticBond = aromaticBond;
    options.ringBond = ringBond;
    options.ringMembership = ringMembership;
    options.ringSize = ringSize;
    options.ringSizeCap = ringSizeCap;
    options.cipAtoms = cipAtoms;
    options.cipBonds = cipBonds;
    options.fragments = fragments;
    options.scaffold = scaffold;
    options.keys = new ArrayList<>(keys);
    return options;
  }

  // ---- MACCS keys ----

  /**
   * The MACCS key definitions as SMARTS with a minimum count of distinct matches; keys 1, 125 and
   * 166 are the isotope, aromatic ring and fragment counts. The native core holds the same table.
   * The SMARTS are adapted from the RDKit file rdkit/Chem/MACCSkeys.py, Copyright (C) 2001-2011
   * Greg Landrum and Rational Discovery LLC, and are used under the BSD 3-Clause License; see
   * licenses/rdkit/LICENSE.
   */
  private static final String[] MACCS_SMARTS = {
    "",                                                           // 1 isotope
    "[#104]",                                                     // 2
    "[#32,#33,#34,#50,#51,#52,#82,#83,#84]",                      // 3
    "[Ac,Th,Pa,U,Np,Pu,Am,Cm,Bk,Cf,Es,Fm,Md,No,Lr]",              // 4
    "[Sc,Ti,Y,Zr,Hf]",                                            // 5
    "[La,Ce,Pr,Nd,Pm,Sm,Eu,Gd,Tb,Dy,Ho,Er,Tm,Yb,Lu]",             // 6
    "[V,Cr,Mn,Nb,Mo,Tc,Ta,W,Re]",                                 // 7
    "[!#6;!#1]1~*~*~*~1",                                         // 8
    "[Fe,Co,Ni,Ru,Rh,Pd,Os,Ir,Pt]",                               // 9
    "[Be,Mg,Ca,Sr,Ba,Ra]",                                        // 10
    "*1~*~*~*~1",                                                 // 11
    "[Cu,Zn,Ag,Cd,Au,Hg]",                                        // 12
    "[#8]~[#7](~[#6])~[#6]",                                      // 13
    "[#16]-[#16]",                                                // 14
    "[#8]~[#6](~[#8])~[#8]",                                      // 15
    "[!#6;!#1]1~*~*~1",                                           // 16
    "[#6]#[#6]",                                                  // 17
    "[#5,#13,#31,#49,#81]",                                       // 18
    "*1~*~*~*~*~*~*~1",                                           // 19
    "[#14]",                                                      // 20
    "[#6]=[#6](~[!#6;!#1])~[!#6;!#1]",                            // 21
    "*1~*~*~1",                                                   // 22
    "[#7]~[#6](~[#8])~[#8]",                                      // 23
    "[#7]-[#8]",                                                  // 24
    "[#7]~[#6](~[#7])~[#7]",                                      // 25
    "[#6]=;@[#6](@*)@*",                                          // 26
    "[I]",                                                        // 27
    "[!#6;!#1]~[CH2]~[!#6;!#1]",                                  // 28
    "[#15]",                                                      // 29
    "[#6]~[!#6;!#1](~[#6])(~[#6])~*",                             // 30
    "[!#6;!#1]~[F,Cl,Br,I]",                                      // 31
    "[#6]~[#16]~[#7]",                                            // 32
    "[#7]~[#16]",                                                 // 33
    "[CH2]=*",                                                    // 34
    "[Li,Na,K,Rb,Cs,Fr]",                                         // 35
    "[#16R]",                                                     // 36
    "[#7]~[#6](~[#8])~[#7]",                                      // 37
    "[#7]~[#6](~[#6])~[#7]",                                      // 38
    "[#8]~[#16](~[#8])~[#8]",                                     // 39
    "[#16]-[#8]",                                                 // 40
    "[#6]#[#7]",                                                  // 41
    "F",                                                          // 42
    "[!#6;!#1;!H0]~*~[!#6;!#1;!H0]",                              // 43
    "[!#1;!#6;!#7;!#8;!#9;!#14;!#15;!#16;!#17;!#35;!#53]",        // 44
    "[#6]=[#6]~[#7]",                                             // 45
    "Br",                                                         // 46
    "[#16]~*~[#7]",                                               // 47
    "[#8]~[!#6;!#1](~[#8])(~[#8])",                               // 48
    "[!+0]",                                                      // 49
    "[#6]=[#6](~[#6])~[#6]",                                      // 50
    "[#6]~[#16]~[#8]",                                            // 51
    "[#7]~[#7]",                                                  // 52
    "[!#6;!#1;!H0]~*~*~*~[!#6;!#1;!H0]",                          // 53
    "[!#6;!#1;!H0]~*~*~[!#6;!#1;!H0]",                            // 54
    "[#8]~[#16]~[#8]",                                            // 55
    "[#8]~[#7](~[#8])~[#6]",                                      // 56
    "[#8R]",                                                      // 57
    "[!#6;!#1]~[#16]~[!#6;!#1]",                                  // 58
    "[#16]!:*:*",                                                 // 59
    "[#16]=[#8]",                                                 // 60
    "*~[#16](~*)~*",                                              // 61
    "*@*!@*@*",                                                   // 62
    "[#7]=[#8]",                                                  // 63
    "*@*!@[#16]",                                                 // 64
    "c:n",                                                        // 65
    "[#6]~[#6](~[#6])(~[#6])~*",                                  // 66
    "[!#6;!#1]~[#16]",                                            // 67
    "[!#6;!#1;!H0]~[!#6;!#1;!H0]",                                // 68
    "[!#6;!#1]~[!#6;!#1;!H0]",                                    // 69
    "[!#6;!#1]~[#7]~[!#6;!#1]",                                   // 70
    "[#7]~[#8]",                                                  // 71
    "[#8]~*~*~[#8]",                                              // 72
    "[#16]=*",                                                    // 73
    "[CH3]~*~[CH3]",                                              // 74
    "*!@[#7]@*",                                                  // 75
    "[#6]=[#6](~*)~*",                                            // 76
    "[#7]~*~[#7]",                                                // 77
    "[#6]=[#7]",                                                  // 78
    "[#7]~*~*~[#7]",                                              // 79
    "[#7]~*~*~*~[#7]",                                            // 80
    "[#16]~*(~*)~*",                                              // 81
    "*~[CH2]~[!#6;!#1;!H0]",                                      // 82
    "[!#6;!#1]1~*~*~*~*~1",                                       // 83
    "[NH2]",                                                      // 84
    "[#6]~[#7](~[#6])~[#6]",                                      // 85
    "[C;H2,H3][!#6;!#1][C;H2,H3]",                                // 86
    "[F,Cl,Br,I]!@*@*",                                           // 87
    "[#16]",                                                      // 88
    "[#8]~*~*~*~[#8]",                                            // 89
    "[$([!#6;!#1;!H0]~*~*~[CH2]~*),$([!#6;!#1;!H0;R]1@[R]@[R]@[CH2;R]1),"
        + "$([!#6;!#1;!H0]~[R]1@[R]@[CH2;R]1)]",                  // 90
    "[$([!#6;!#1;!H0]~*~*~*~[CH2]~*),$([!#6;!#1;!H0;R]1@[R]@[R]@[R]@[CH2;R]1),"
        + "$([!#6;!#1;!H0]~[R]1@[R]@[R]@[CH2;R]1),$([!#6;!#1;!H0]~*~[R]1@[R]@[CH2;R]1)]", // 91
    "[#8]~[#6](~[#7])~[#6]",                                      // 92
    "[!#6;!#1]~[CH3]",                                            // 93
    "[!#6;!#1]~[#7]",                                             // 94
    "[#7]~*~*~[#8]",                                              // 95
    "*1~*~*~*~*~1",                                               // 96
    "[#7]~*~*~*~[#8]",                                            // 97
    "[!#6;!#1]1~*~*~*~*~*~1",                                     // 98
    "[#6]=[#6]",                                                  // 99
    "*~[CH2]~[#7]",                                               // 100
    "[$([R]@1@[R]@[R]@[R]@[R]@[R]@[R]@[R]1),$([R]@1@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]1),"
        + "$([R]@1@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]1),$([R]@1@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]1),"
        + "$([R]@1@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]1),"
        + "$([R]@1@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]1),"
        + "$([R]@1@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]1)]", // 101
    "[!#6;!#1]~[#8]",                                             // 102
    "Cl",                                                         // 103
    "[!#6;!#1;!H0]~*~[CH2]~*",                                    // 104
    "*@*(@*)@*",                                                  // 105
    "[!#6;!#1]~*(~[!#6;!#1])~[!#6;!#1]",                          // 106
    "[F,Cl,Br,I]~*(~*)~*",                                        // 107
    "[CH3]~*~*~*~[CH2]~*",                                        // 108
    "*~[CH2]~[#8]",                                               // 109
    "[#7]~[#6]~[#8]",                                             // 110
    "[#7]~*~[CH2]~*",                                             // 111
    "*~*(~*)(~*)~*",                                              // 112
    "[#8]!:*:*",                                                  // 113
    "[CH3]~[CH2]~*",                                              // 114
    "[CH3]~*~[CH2]~*",                                            // 115
    "[$([CH3]~*~*~[CH2]~*),$([CH3]~*1~*~[CH2]1)]",                // 116
    "[#7]~*~[#8]",                                                // 117
    "[$(*~[CH2]~[CH2]~*),$(*1~[CH2]~[CH2]1)]",                    // 118
    "[#7]=*",                                                     // 119
    "[!#6;R]",                                                    // 120
    "[#7;R]",                                                     // 121
    "*~[#7](~*)~*",                                               // 122
    "[#8]~[#6]~[#8]",                                             // 123
    "[!#6;!#1]~[!#6;!#1]",                                        // 124
    "",                                                           // 125 aromatic rings
    "*!@[#8]!@*",                                                 // 126
    "*@*!@[#8]",                                                  // 127
    "[$(*~[CH2]~*~*~*~[CH2]~*),$([R]1@[CH2;R]@[R]@[R]@[R]@[CH2;R]1),"
        + "$(*~[CH2]~[R]1@[R]@[R]@[CH2;R]1),$(*~[CH2]~*~[R]1@[R]@[CH2;R]1)]", // 128
    "[$(*~[CH2]~*~*~[CH2]~*),$([R]1@[CH2]@[R]@[R]@[CH2;R]1),$(*~[CH2]~[R]1@[R]@[CH2;R]1)]", // 129
    "[!#6;!#1]~[!#6;!#1]",                                        // 130
    "[!#6;!#1;!H0]",                                              // 131
    "[#8]~*~[CH2]~*",                                             // 132
    "*@*!@[#7]",                                                  // 133
    "[F,Cl,Br,I]",                                                // 134
    "[#7]!:*:*",                                                  // 135
    "[#8]=*",                                                     // 136
    "[!C;!c;R]",                                                  // 137
    "[!#6;!#1]~[CH2]~*",                                          // 138
    "[O;!H0]",                                                    // 139
    "[#8]",                                                       // 140
    "[CH3]",                                                      // 141
    "[#7]",                                                       // 142
    "*@*!@[#8]",                                                  // 143
    "*!:*:*!:*",                                                  // 144
    "*1~*~*~*~*~*~1",                                             // 145
    "[#8]",                                                       // 146
    "[$(*~[CH2]~[CH2]~*),$([R]1@[CH2;R]@[CH2;R]1)]",              // 147
    "*~[!#6;!#1](~*)~*",                                          // 148
    "[C;H3,H4]",                                                  // 149
    "*!@*@*!@*",                                                  // 150
    "[#7;!H0]",                                                   // 151
    "[#8]~[#6](~[#6])~[#6]",                                      // 152
    "[!#6;!#1]~[CH2]~*",                                          // 153
    "[#6]=[#8]",                                                  // 154
    "*!@[CH2]!@*",                                                // 155
    "[#7]~*(~*)~*",                                               // 156
    "[#6]-[#8]",                                                  // 157
    "[#6]-[#7]",                                                  // 158
    "[#8]",                                                       // 159
    "[C;H3,H4]",                                                  // 160
    "[#7]",                                                       // 161
    "a",                                                          // 162
    "*1~*~*~*~*~*~1",                                             // 163
    "[#8]",                                                       // 164
    "[R]",                                                        // 165
    "",                                                           // 166 fragments
  };

  /** The minimum count of each MACCS key, parallel to {@link #MACCS_SMARTS}. */
  private static final int[] MACCS_MINIMUM = {
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1,   //   1 -  10
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1,   //  11 -  20
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1,   //  21 -  30
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1,   //  31 -  40
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1,   //  41 -  50
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1,   //  51 -  60
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1,   //  61 -  70
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1,   //  71 -  80
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1,   //  81 -  90
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1,   //  91 - 100
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1,   // 101 - 110
    1, 1, 1, 1, 1, 1, 1, 2, 1, 2,   // 111 - 120
    1, 1, 1, 1, 2, 1, 2, 1, 1, 2,   // 121 - 130
    2, 1, 1, 1, 1, 2, 1, 2, 1, 4,   // 131 - 140
    3, 2, 1, 1, 2, 3, 1, 1, 2, 1,   // 141 - 150
    1, 1, 1, 1, 1, 1, 1, 1, 2, 1,   // 151 - 160
    1, 1, 1, 1, 1, 2,               // 161 - 166
  };

  private static final List<Key> MACCS_KEYS;

  static {
    List<Key> keys = new ArrayList<>(MACCS_SMARTS.length);
    for (int index = 0; index < MACCS_SMARTS.length; index++) {
      int minimum = MACCS_MINIMUM[index];
      if (index == 0) keys.add(Key.special(Key.Kind.ISOTOPE, minimum));
      else if (index == 124) keys.add(Key.special(Key.Kind.AROMATIC_RINGS, minimum));
      else if (index == 165) keys.add(Key.special(Key.Kind.FRAGMENTS, minimum));
      else keys.add(new Key(MACCS_SMARTS[index], minimum));
    }
    MACCS_KEYS = Collections.unmodifiableList(keys);
  }

  /** The MACCS key definitions, in key order. */
  public static List<Key> maccsKeys() { return MACCS_KEYS; }
}
