/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms. */
package com.bioinception.smsd.core;

/**
 * Chemical matching constraints for substructure and MCS searches.
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
 * @see MCS
 * @see Substructure
 */
public final class ChemOptions {

  /**
   * Bond order comparison mode.
   *
   * <ul>
   *   <li>{@link #STRICT} - bond orders must match exactly</li>
   *   <li>{@link #LOOSE} - bond orders are compatible, subject to aromaticity and ring constraints</li>
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
   * <p>MCS retains ring atoms and bonds only when each is covered by a wholly selected
   * relevant ring in both inputs. Shared fusion atoms may belong to one retained whole ring.
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
   * pH context for the explicit four-argument MCS.computeTautomerConfidence helper.
   * Default 7.4; range 0–14. This heuristic does not alter the MCIS or MCES objective.
   */
  public double pH = 7.4;

  /** Set the pH context for explicit tautomer confidence scoring. Default 7.4. */
  public ChemOptions withPH(double pH) {
    if (!Double.isFinite(pH))
      throw new IllegalArgumentException("pH must be finite");
    this.pH = Math.max(0, Math.min(14, pH));
    return this;
  }

  /**
   * Solvent context for the explicit four-argument MCS.computeTautomerConfidence helper.
   * Default: {@link Solvent#AQUEOUS}. This heuristic does not alter the MCIS or MCES objective.
   */
  public Solvent solvent = Solvent.AQUEOUS;

  /** Set the solvent context for explicit tautomer confidence scoring. */
  public ChemOptions withSolvent(Solvent s) {
    this.solvent = s != null ? s : Solvent.AQUEOUS;
    return this;
  }

  // --- Bond matching ---

  /** Bond order matching mode. Default: {@link BondOrderMode#STRICT}. */
  public BondOrderMode matchBondOrder = BondOrderMode.STRICT;

  /** Aromaticity matching mode. Default: {@link AromaticityMode#FLEXIBLE}. */
  public AromaticityMode aromaticityMode = AromaticityMode.FLEXIBLE;

  // --- Pruning ---

  /** Enable two-hop neighbourhood label frequency pruning. Default: {@code true}. */
  public boolean useTwoHopNLF = true;

  /** Enable three-hop neighbourhood label frequency pruning. Default: {@code false}.
   *  Can be enabled for very large molecules where extra pruning helps.
   *  Disabled for small molecules (n &lt; 20). */
  public boolean useThreeHopNLF = false;

  /** Enable bit-parallel feasibility checks for faster pruning. Default: {@code true}. */
  public boolean useBitParallelFeasibility = true;


  /**
   * Create a {@code ChemOptions} instance with all default settings.
   *
   * <pre>{@code
   * ChemOptions opts = new ChemOptions();
   * // defaults: strict bond order and flexible aromaticity
   * }</pre>
   */
  public ChemOptions() {}

  /**
   * Create a deep copy of the given {@code ChemOptions}.
   *
   * <p>Useful when a caller needs a modified copy without mutating the original
   * (e.g. temporarily relaxing {@code matchFormalCharge} for one query).
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
    c.useTwoHopNLF             = src.useTwoHopNLF;
    c.useThreeHopNLF           = src.useThreeHopNLF;
    c.useBitParallelFeasibility = src.useBitParallelFeasibility;
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
    switch (name.toLowerCase(java.util.Locale.ROOT)) {
      case "default":
        break;
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
      case "pharma":
        opts.matchFormalCharge = true;
        opts.ringMatchesRingOnly = true;
        opts.completeRingsOnly = true;
        break;
      case "reaction":
        opts.matchBondOrder = BondOrderMode.LOOSE;
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
   * Pre-configured profile for tautomer-aware MCS: loose bond orders, flexible aromaticity.
   *
   * <pre>{@code
   * ChemOptions tOpts = ChemOptions.tautomerProfile();
   * MCS.Options options=new MCS.Options();
   * options.chemistry=tOpts;
   * SearchResult result=MCS.find(query).withOptions(options).withTimeout(10000).match(target);
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
   * <p>Uses loose bond orders, flexible aromaticity and permits ring-to-chain and partial-ring matches.
   *
   * <pre>{@code
   * ChemOptions loose = ChemOptions.fmcsProfile();
   * MCS.Options options=new MCS.Options();
   * options.chemistry=loose;
   * SearchResult result=MCS.find(mol1).withOptions(options).match(mol2);
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
