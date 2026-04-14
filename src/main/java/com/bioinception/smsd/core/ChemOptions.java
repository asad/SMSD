/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms. */
package com.bioinception.smsd.core;

/**
 * Configuration for chemical matching constraints and engine selection.
 *
 * <p>Controls how atoms and bonds are compared during substructure and MCS searches.
 * Pre-built profiles: {@code "strict"} (strict bond order + aromaticity),
 * {@code "compat-substruct"} / {@code "compat-fmcs"} (loose, flexible).
 *
 * <pre>{@code
 * // Default options (strict bond order, flexible aromaticity)
 * ChemOptions opts = new ChemOptions();
 *
 * // Pre-configured loose profile
 * ChemOptions loose = ChemOptions.profile("compat-fmcs");
 *
 * // Fluent builder style
 * ChemOptions custom = new ChemOptions()
 *     .withBondOrderMode(ChemOptions.BondOrderMode.LOOSE)
 *     .withCompleteRingsOnly(true);
 * }</pre>
 *
 * @author Syed Asad Rahman
 * @see SMSD
 */
public final class ChemOptions {

  /**
   * Bond order comparison mode.
   *
   * <ul>
   *   <li>{@link #STRICT} - bond orders must match exactly</li>
   *   <li>{@link #LOOSE} - single/double/aromatic bonds are considered compatible</li>
   *   <li>{@link #ANY} - all bond orders are compatible</li>
   * </ul>
   */
  public enum BondOrderMode { STRICT, LOOSE, ANY }

  /**
   * Aromaticity comparison mode.
   *
   * <ul>
   *   <li>{@link #STRICT} - aromatic bonds only match aromatic bonds</li>
   *   <li>{@link #FLEXIBLE} - aromatic bonds can match single or double bonds</li>
   * </ul>
   */
  public enum AromaticityMode { STRICT, FLEXIBLE }

  /**
   * Subgraph isomorphism engine selection.
   *
   * <ul>
   *   <li>{@link #VF2} - classic VF2 algorithm</li>
   *   <li>{@link #VF2PP} - VF2++ with improved pruning (default, recommended)</li>
   *   <li>{@link #VF3} - VF3 with lightweight pruning</li>
   * </ul>
   */
  public enum MatcherEngine { VF2, VF2PP, VF3 }

  /**
   * Ring fusion matching mode.
   *
   * <ul>
   *   <li>{@link #IGNORE} - ring fusion is not considered</li>
   *   <li>{@link #PERMISSIVE} - fused rings are preferred but not required</li>
   *   <li>{@link #STRICT} - fused ring systems must match exactly</li>
   * </ul>
   */
  public enum RingFusionMode { IGNORE, PERMISSIVE, STRICT }

  /**
   * Solvent environment for Tier 2 pKa-informed tautomer weight corrections.
   *
   * <ul>
   *   <li>{@link #AQUEOUS} - water (default, Tier 1 weights unchanged)</li>
   *   <li>{@link #DMSO} - dimethyl sulfoxide (aprotic, enol more prevalent)</li>
   *   <li>{@link #METHANOL} - protic, similar to aqueous with small shifts</li>
   *   <li>{@link #CHLOROFORM} - aprotic, large keto-enol equilibrium shift</li>
   *   <li>{@link #ACETONITRILE} - polar aprotic, moderate keto-enol shift</li>
   *   <li>{@link #DIETHYL_ETHER} - non-polar aprotic, significant keto-enol shift</li>
   * </ul>
   */
  public enum Solvent {
    AQUEOUS,
    DMSO,
    METHANOL,
    CHLOROFORM,
    ACETONITRILE,
    DIETHYL_ETHER
  }

  // --- Atom matching ---

  /** Whether to match atom types (element symbols). Default: {@code true}. */
  public boolean matchAtomType = true;

  /** Whether to match formal charges on atoms. Default: {@code false} (RDKit-compatible). */
  public boolean matchFormalCharge = false;

  /** Whether to use chirality (R/S) in atom matching. Default: {@code false}. */
  public boolean useChirality = false;

  /** Whether to use bond stereochemistry (E/Z) in matching. Default: {@code false}. */
  public boolean useBondStereo = false;

  /**
   * Whether ring atoms can only match ring atoms. Default: {@code false} (RDKit-compatible).
   * <p>When enabled, matched atoms and bonds must agree on ring membership.
   */
  public boolean ringMatchesRingOnly = false;

  // --- Ring completeness ---

  /**
   * Whether to require complete ring matching. Default: {@code false}.
   * <p>When enabled, if any atom in a ring is mapped, all atoms in that ring must be mapped.
   */
  public boolean completeRingsOnly = false;

  // --- Isotope matching ---

  /** Whether to match isotope labels on atoms. Default: {@code false}. */
  public boolean matchIsotope = false;

  // --- Ring fusion ---

  /** Ring fusion matching mode. Default: {@link RingFusionMode#IGNORE}. */
  public RingFusionMode ringFusionMode = RingFusionMode.IGNORE;

  // --- Tautomer matching ---

  /**
   * Whether to consider tautomeric forms when matching. Default: {@code false}.
   * <p>When enabled, tautomeric bonds (e.g., keto/enol) are considered equivalent.
   */
  public boolean tautomerAware = false;

  /**
   * pH used for pKa-informed tautomer relevance scoring. Default: {@code 7.4} (physiological).
   * <p>Only used when {@code tautomerAware=true}. Affects which tautomeric form is considered
   * dominant; the reported {@code tautomerConfidence} reflects match quality at this pH.
   * Range: 0–14.
   */
  public double pH = 7.4;

  /** Set the simulation pH for tautomer relevance scoring (fluent). Default 7.4. */
  public ChemOptions withPH(double pH) {
    this.pH = Math.max(0, Math.min(14, pH));
    return this;
  }

  /**
   * Solvent environment for Tier 2 pKa tautomer weight corrections. Default: {@link
   * Solvent#AQUEOUS} (no correction applied).
   */
  public Solvent solvent = Solvent.AQUEOUS;

  /** Set the solvent for tautomer weight corrections (fluent). */
  public ChemOptions withSolvent(Solvent s) {
    this.solvent = s != null ? s : Solvent.AQUEOUS;
    return this;
  }

  // --- Bond matching ---

  /** Bond order matching mode. Default: {@link BondOrderMode#STRICT}. */
  public BondOrderMode matchBondOrder = BondOrderMode.STRICT;

  /** Aromaticity matching mode. Default: {@link AromaticityMode#FLEXIBLE}. */
  public AromaticityMode aromaticityMode = AromaticityMode.FLEXIBLE;

  // --- Engine ---

  /** Subgraph isomorphism engine. Default: {@link MatcherEngine#VF2PP}. */
  public MatcherEngine matcherEngine = MatcherEngine.VF2PP;

  // --- Pruning ---

  /** Enable two-hop neighbourhood label frequency pruning. Default: {@code true}. */
  public boolean useTwoHopNLF = true;

  /** Enable three-hop neighbourhood label frequency pruning. Default: {@code false}.
   *  Can be enabled for very large molecules where extra pruning helps.
   *  Auto-disabled for small molecules (n &lt; 20) in VF2++. */
  public boolean useThreeHopNLF = false;

  /** Enable bit-parallel feasibility checks for faster pruning. Default: {@code true}. */
  public boolean useBitParallelFeasibility = true;

  /**
   * When {@code true}, SMILES strings passed to the SMSD(String, String, ChemOptions) constructor
   * are pre-sanitised before CDK parsing: unbalanced parentheses are corrected, and unclosed
   * ring-digit openings are removed. Default: {@code false} (strict parsing).
   */
  public boolean lenientSmiles = false;

  /**
   * Create a {@code ChemOptions} instance with all default settings.
   *
   * <pre>{@code
   * ChemOptions opts = new ChemOptions();
   * // defaults: strict bond order, flexible aromaticity, VF2PP engine
   * }</pre>
   */
  public ChemOptions() {}

  /**
   * Create a deep copy of the given {@code ChemOptions}.
   *
   * <p>Useful when a caller needs a modified copy without mutating the original
   * (e.g., reaction-aware MCS relaxes {@code matchFormalCharge}).
   *
   * @param src the source options to copy
   * @return a new independent copy with all fields duplicated
   * @since 6.5.0
   */
  public static ChemOptions copyOf(ChemOptions src) {
    if (src == null) throw new NullPointerException("src must not be null");
    ChemOptions c = new ChemOptions();
    c.matchAtomType            = src.matchAtomType;
    c.matchFormalCharge        = src.matchFormalCharge;
    c.useChirality             = src.useChirality;
    c.useBondStereo            = src.useBondStereo;
    c.ringMatchesRingOnly      = src.ringMatchesRingOnly;
    c.completeRingsOnly        = src.completeRingsOnly;
    c.matchIsotope             = src.matchIsotope;
    c.ringFusionMode           = src.ringFusionMode;
    c.tautomerAware            = src.tautomerAware;
    c.pH                       = src.pH;
    c.solvent                  = src.solvent;
    c.matchBondOrder           = src.matchBondOrder;
    c.aromaticityMode          = src.aromaticityMode;
    c.matcherEngine            = src.matcherEngine;
    c.useTwoHopNLF             = src.useTwoHopNLF;
    c.useThreeHopNLF           = src.useThreeHopNLF;
    c.useBitParallelFeasibility = src.useBitParallelFeasibility;
    c.lenientSmiles            = src.lenientSmiles;
    return c;
  }

  /**
   * Create pre-configured options from a named profile.
   *
   * <p>Available profiles:
   * <ul>
   *   <li>{@code "strict"} - strict bond order and aromaticity</li>
   *   <li>{@code "compat-fmcs"} / {@code "compat-substruct"} - loose bond order,
   *       flexible aromaticity, no ring-only constraint</li>
   * </ul>
   *
   * <pre>{@code
   * ChemOptions loose = ChemOptions.profile("compat-fmcs");
   * }</pre>
   *
   * @param name the profile name (case-insensitive)
   * @return a new {@code ChemOptions} configured according to the profile
   */
  public static ChemOptions profile(String name) {
    if (name == null) throw new NullPointerException("profile name must not be null");
    ChemOptions opts = new ChemOptions();
    switch (name.toLowerCase()) {
      case "strict":
        opts.matchBondOrder = BondOrderMode.STRICT;
        opts.aromaticityMode = AromaticityMode.STRICT;
        break;
      case "compat-fmcs":
      case "compat-substruct":
        opts.matchBondOrder = BondOrderMode.LOOSE;
        opts.aromaticityMode = AromaticityMode.FLEXIBLE;
        opts.ringMatchesRingOnly = false;
        break;
      default:
        throw new IllegalArgumentException("Unknown ChemOptions profile: " + name);
    }
    return opts;
  }

  /**
   * Set the aromaticity matching mode (fluent).
   *
   * @param m the aromaticity mode; {@code null} defaults to {@link AromaticityMode#FLEXIBLE}
   * @return this instance for chaining
   */
  public ChemOptions withAromaticityMode(AromaticityMode m) {
    this.aromaticityMode = m != null ? m : AromaticityMode.FLEXIBLE;
    return this;
  }

  /**
   * Set the bond order matching mode (fluent).
   *
   * @param m the bond order mode; {@code null} defaults to {@link BondOrderMode#STRICT}
   * @return this instance for chaining
   */
  public ChemOptions withBondOrderMode(BondOrderMode m) {
    this.matchBondOrder = m != null ? m : BondOrderMode.STRICT;
    return this;
  }

  /**
   * Enable or disable two-hop neighbourhood label frequency pruning (fluent).
   *
   * @param on {@code true} to enable
   * @return this instance for chaining
   */
  public ChemOptions withTwoHopNLF(boolean on) {
    this.useTwoHopNLF = on;
    return this;
  }

  /**
   * Enable or disable three-hop neighbourhood label frequency pruning (fluent).
   *
   * @param on {@code true} to enable
   * @return this instance for chaining
   */
  public ChemOptions withThreeHopNLF(boolean on) {
    this.useThreeHopNLF = on;
    return this;
  }

  /**
   * Enable or disable bit-parallel feasibility checks (fluent).
   *
   * @param on {@code true} to enable
   * @return this instance for chaining
   */
  public ChemOptions withBitParallelFeasibility(boolean on) {
    this.useBitParallelFeasibility = on;
    return this;
  }

  /**
   * Enable or disable complete ring matching (fluent).
   *
   * <p>When enabled, if any atom in a ring is part of the MCS, all ring atoms must be included.
   *
   * @param on {@code true} to require complete ring matching
   * @return this instance for chaining
   */
  public ChemOptions withCompleteRingsOnly(boolean on) {
    this.completeRingsOnly = on;
    return this;
  }

  /**
   * Enable or disable isotope matching (fluent).
   *
   * @param on {@code true} to require isotope labels to match
   * @return this instance for chaining
   */
  public ChemOptions withMatchIsotope(boolean on) {
    this.matchIsotope = on;
    return this;
  }

  /**
   * Set the ring fusion matching mode (fluent).
   *
   * @param m the ring fusion mode; {@code null} defaults to {@link RingFusionMode#IGNORE}
   * @return this instance for chaining
   */
  public ChemOptions withRingFusionMode(RingFusionMode m) {
    this.ringFusionMode = m != null ? m : RingFusionMode.IGNORE;
    return this;
  }

  /**
   * Enable or disable tautomer-aware matching (fluent).
   *
   * @param on {@code true} to consider tautomeric equivalence
   * @return this instance for chaining
   */
  public ChemOptions withTautomerAware(boolean on) {
    this.tautomerAware = on;
    return this;
  }

  /**
   * Enable or disable lenient SMILES pre-sanitisation (fluent).
   *
   * @param on {@code true} to pre-sanitise SMILES before CDK parsing
   * @return this instance for chaining
   */
  public ChemOptions withLenientSmiles(boolean on) {
    this.lenientSmiles = on;
    return this;
  }

  /**
   * Pre-configured profile for tautomer-aware MCS: loose bond orders, flexible aromaticity.
   *
   * <pre>{@code
   * ChemOptions tOpts = ChemOptions.tautomerProfile();
   * Map<Integer, Integer> mcs = smsd.findMCS(false, true, 10000);
   * }</pre>
   *
   * @return a new {@code ChemOptions} configured for tautomer-aware matching
   */
  public static ChemOptions tautomerProfile() {
    ChemOptions c = new ChemOptions();
    c.tautomerAware = true;
    c.matchBondOrder = BondOrderMode.LOOSE;
    c.aromaticityMode = AromaticityMode.FLEXIBLE;
    c.matchFormalCharge = false;
    c.ringMatchesRingOnly = false;
    return c;
  }

  /**
   * Pre-configured profile that relaxes topological constraints for loose FMCS-style matching.
   *
   * <p>SMSD's default mode enforces chemical rigour: ring atoms must match ring atoms, partial
   * rings are rejected, and bond orders are checked strictly. This profile lifts those constraints
   * so that SMSD will map ring atoms to chain atoms and accept partial ring fragments — the same
   * loose behaviour as standard FMCS implementations. The result may be numerically larger but is
   * less chemically meaningful. Use this profile for interoperability benchmarks or when the
   * application explicitly needs loose topology matching.
   *
   * <pre>{@code
   * ChemOptions loose = ChemOptions.fmcsProfile();
   * Map<Integer, Integer> mcs = new SMSD(mol1, mol2, loose).findMCS();
   * }</pre>
   *
   * @return a new {@code ChemOptions} configured for loose FMCS-style matching
   */
  public static ChemOptions fmcsProfile() {
    ChemOptions c = new ChemOptions();
    c.matchBondOrder      = BondOrderMode.LOOSE;
    c.ringMatchesRingOnly = false;
    c.completeRingsOnly   = false;
    c.matchFormalCharge   = false;
    return c;
  }
}
