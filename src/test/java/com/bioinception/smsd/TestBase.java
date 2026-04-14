/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
package com.bioinception.smsd;

import com.bioinception.smsd.core.Standardiser;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;

/**
 * Shared utilities for SMSD test classes. Provides a single SmilesParser and standardised molecule
 * parsing.
 */
public class TestBase {

  protected static final SmilesParser SP = new SmilesParser(SilentChemObjectBuilder.getInstance());

  /** Parse SMILES and standardise (no tautomer normalisation). */
  protected static IAtomContainer mol(String smi) throws Exception {
    return Standardiser.standardise(SP.parseSmiles(smi), Standardiser.TautomerMode.NONE);
  }
}
