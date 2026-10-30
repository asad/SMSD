/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
package com.bioinception.smsd.core;

import java.math.BigInteger;
import java.util.Collections;
import java.util.List;

/**
 * The relevant cycles of a {@link MolGraph} grouped into families, in ascending ring size and then
 * root; see {@link RelevantCycleFamily}. {@link #truncated()} reports an incomplete enumeration.
 * When the enumeration is complete, the relevant cycles and their number do not depend on atom
 * order; the grouping into families, and so each family's record, can.
 */
public final class RelevantCycleFamilies {
  private final List<RelevantCycleFamily> families;
  private final BigInteger relevantCycleCount;
  private final boolean truncated;

  RelevantCycleFamilies(List<RelevantCycleFamily> families, BigInteger relevantCycleCount, boolean truncated) {
    this.families = Collections.unmodifiableList(families);
    this.relevantCycleCount = relevantCycleCount;
    this.truncated = truncated;
  }

  /** The families, as an unmodifiable list. */
  public List<RelevantCycleFamily> families() { return families; }

  /** The exact number of members of all families listed: every relevant cycle unless truncated. */
  public BigInteger relevantCycleCount() { return relevantCycleCount; }

  /** Whether the enumeration is incomplete, so that families may be missing. */
  public boolean truncated() { return truncated; }

  @Override public String toString() {
    return "RelevantCycleFamilies{families=" + families.size() + ", relevantCycles=" + relevantCycleCount
        + ", truncated=" + truncated + "}";
  }
}
