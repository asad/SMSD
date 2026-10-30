/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
package com.bioinception.smsd.core;

import java.math.BigInteger;
import java.util.Arrays;

/**
 * One family of relevant cycles of a {@link MolGraph}: its ring size, the highest-numbered atom
 * shared by its members, the atoms opposite it, one member and the exact number of members. Atoms
 * are graph indices. See {@link MolGraph#computeRelevantCycleFamilies()}.
 */
public final class RelevantCycleFamily {
  private final int size, root, first, second, opposite;
  private final int[] prototype;
  private final BigInteger memberCount;

  RelevantCycleFamily(int size, int root, int first, int second, int opposite, int[] prototype,
                      BigInteger memberCount) {
    this.size = size;
    this.root = root;
    this.first = first;
    this.second = second;
    this.opposite = opposite;
    this.prototype = prototype;
    this.memberCount = memberCount;
  }

  /** The number of atoms of every member. */
  public int size() { return size; }

  /** The highest-numbered atom, which every member shares. */
  public int root() { return root; }

  /**
   * The lower-numbered of the two atoms opposite the root: an end of the opposite bond for an odd
   * size, a ring neighbour of the opposite atom for an even size.
   */
  public int first() { return first; }

  /** The higher-numbered of the two atoms opposite the root, as for {@link #first()}. */
  public int second() { return second; }

  /** The atom opposite the root for an even size, otherwise -1. */
  public int opposite() { return opposite; }

  /**
   * One member, as a copy: its atoms in ring order, starting at the root and reaching
   * {@link #first()} before {@link #second()}.
   */
  public int[] prototype() { return prototype.clone(); }

  /** The exact number of members. */
  public BigInteger memberCount() { return memberCount; }

  @Override public String toString() {
    return "RelevantCycleFamily{size=" + size + ", root=" + root + ", first=" + first + ", second=" + second
        + ", opposite=" + opposite + ", members=" + memberCount + ", prototype=" + Arrays.toString(prototype) + "}";
  }
}
