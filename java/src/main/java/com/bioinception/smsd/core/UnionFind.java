/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
package com.bioinception.smsd.core;

/** Disjoint sets over item indices, stored as a parent array. */
final class UnionFind {
  private UnionFind() {}

  /** A union-find over {@code size} items, each in a set of its own: the parent of each item is itself. */
  static int[] create(int size) {
    int[] parent = new int[size];
    for (int i = 0; i < size; i++) parent[i] = i;
    return parent;
  }

  /** The root of the set of {@code item}, halving the path to it on the way. */
  static int root(int[] parent, int item) {
    while (parent[item] != item) item = parent[item] = parent[parent[item]];
    return item;
  }

  /** Joins the sets of {@code first} and {@code second}; the root of {@code first} is placed below the root of {@code second}. */
  static void join(int[] parent, int first, int second) {
    first = root(parent, first);
    second = root(parent, second);
    if (first != second) parent[first] = second;
  }
}
