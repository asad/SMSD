/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * See the NOTICE file for attribution.
 */
package com.bioinception.smsd;

import java.util.ArrayList;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.Set;

/** Immutable test fixtures and text helpers available on Java 8. */
public final class TestSupport {
  private TestSupport() {}

  public static Map<Integer, Integer> mapping(int... pairs) {
    if ((pairs.length & 1) != 0) throw new IllegalArgumentException("Incomplete mapping pair");
    Map<Integer, Integer> values = new LinkedHashMap<>();
    for (int i = 0; i < pairs.length; i += 2) {
      if (values.put(pairs[i], pairs[i + 1]) != null)
        throw new IllegalArgumentException("Duplicate mapping key: " + pairs[i]);
    }
    return Collections.unmodifiableMap(values);
  }

  @SafeVarargs
  public static <T> List<T> list(T... elements) {
    List<T> values = new ArrayList<>(elements.length);
    for (T element : elements) values.add(Objects.requireNonNull(element));
    return Collections.unmodifiableList(values);
  }

  @SafeVarargs
  public static <T> Set<T> set(T... elements) {
    Set<T> values = new LinkedHashSet<>();
    for (T element : elements) {
      if (!values.add(Objects.requireNonNull(element)))
        throw new IllegalArgumentException("Duplicate set element: " + element);
    }
    return Collections.unmodifiableSet(values);
  }

  public static String repeat(String text, int count) {
    if (count < 0) throw new IllegalArgumentException("Negative repeat count");
    StringBuilder result = new StringBuilder();
    for (int i = 0; i < count; i++) result.append(text);
    return result.toString();
  }

  public static boolean isBlank(String text) {
    return text.codePoints().allMatch(Character::isWhitespace);
  }
}
