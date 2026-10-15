/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
package com.bioinception.smsd.core;

import com.fasterxml.jackson.annotation.JsonCreator;
import com.fasterxml.jackson.annotation.JsonProperty;
import java.util.ArrayList;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.function.BooleanSupplier;

/** Mappings and completion information from a molecular search. */
public final class SearchResult {
  public enum Status { Completed, TimedOut, Cancelled, LimitReached }

  /** An unordered bond identified by its original atom indices. */
  public static final class Bond {
    private final int first, second;

    @JsonCreator
    public Bond(@JsonProperty("first") int first, @JsonProperty("second") int second) {
      if (first < 0 || second < 0 || first == second)
        throw new IllegalArgumentException("Bond endpoints must be distinct non-negative indices");
      this.first = Math.min(first, second);
      this.second = Math.max(first, second);
    }

    @JsonProperty("first")
    public int first() {
      return first;
    }
    @JsonProperty("second")
    public int second() {
      return second;
    }

    @Override
    public boolean equals(Object other) {
      if (!(other instanceof Bond))
        return false;
      Bond bond = (Bond) other;
      return first == bond.first && second == bond.second;
    }

    public static Bond valueOf(String key) {
      int separator = key.indexOf('-');
      if (separator < 1)
        throw new IllegalArgumentException("Invalid bond key");
      return new Bond(Integer.parseInt(key.substring(0, separator)),
          Integer.parseInt(key.substring(separator + 1)));
    }

    @Override
    public int hashCode() {
      return 31 * first + second;
    }
    @Override
    public String toString() {
      return first + "-" + second;
    }
  }

  /** One atom mapping, its selected common bonds and their connected fragments. */
  public static final class Match {
    private final Map<Integer, Integer> atomMapping;
    private final Map<Bond, Bond> bondMapping;
    private final List<Set<Integer>> fragments;
    private final double score;

    @JsonCreator
    Match(@JsonProperty("atomMapping") Map<Integer, Integer> atoms,
        @JsonProperty("bondMapping") Map<Bond, Bond> bonds,
        @JsonProperty("fragments") List<Set<Integer>> fragments,
        @JsonProperty("score") double score) {
      this.atomMapping = Collections.unmodifiableMap(new TreeMap<>(atoms));
      this.bondMapping = Collections.unmodifiableMap(new LinkedHashMap<>(bonds));
      List<Set<Integer>> copied = new ArrayList<>();
      for (Set<Integer> fragment : fragments)
        copied.add(Collections.unmodifiableSet(new TreeSet<>(fragment)));
      this.fragments = Collections.unmodifiableList(copied);
      this.score = score;
    }

    @JsonProperty("atomMapping")
    public Map<Integer, Integer> atomMapping() {
      return atomMapping;
    }
    @JsonProperty("bondMapping")
    public Map<Bond, Bond> bondMapping() {
      return bondMapping;
    }
    @JsonProperty("fragments")
    public List<Set<Integer>> fragments() {
      return fragments;
    }
    public int atomCount() {
      return atomMapping.size();
    }
    public int bondCount() {
      return bondMapping.size();
    }
    @JsonProperty("score")
    public double score() {
      return score;
    }

    /** Select source graph properties without capping or re-perceiving a molecular fragment. */
    public MolGraph toGraph(MolGraph query) {
      return MCS.toGraph(query, atomMapping, bondMapping.keySet());
    }

    /** Return the common fragment using the selected bonds of the query. */
    public org.openscience.cdk.interfaces.IAtomContainer toMolecule(MolGraph query) {
      return MCS.toMolecule(query, atomMapping, bondMapping.keySet());
    }

    public org.openscience.cdk.interfaces.IAtomContainer toMolecule(
        org.openscience.cdk.interfaces.IAtomContainer query) {
      return toMolecule(new MolGraph(query));
    }

    public String toSmiles(MolGraph query) {
      return new MolGraph(toMolecule(query)).toCanonicalSmiles();
    }
    public String toSmiles(org.openscience.cdk.interfaces.IAtomContainer query) {
      return toSmiles(new MolGraph(query));
    }

    @Override
    public boolean equals(Object other) {
      if (!(other instanceof Match))
        return false;
      Match match = (Match) other;
      return atomMapping.equals(match.atomMapping) && bondMapping.equals(match.bondMapping)
          && fragments.equals(match.fragments) && Double.compare(score, match.score) == 0;
    }

    @Override
    public int hashCode() {
      return Objects.hash(atomMapping, bondMapping, fragments, score);
    }
  }

  private final List<Match> matches;
  private final Status status;
  private final boolean optimal, complete;
  private final double score;
  private final long nodesVisited, elapsedMillis;

  @JsonCreator
  SearchResult(@JsonProperty("matches") List<Match> matches, @JsonProperty("status") Status status,
      @JsonProperty("optimal") boolean optimal, @JsonProperty("complete") boolean complete,
      @JsonProperty("score") double score, @JsonProperty("nodesVisited") long nodesVisited,
      @JsonProperty("elapsedMillis") long elapsedMillis) {
    this.matches = Collections.unmodifiableList(new ArrayList<>(matches));
    this.status = Objects.requireNonNull(status, "status");
    this.optimal = optimal;
    this.complete = complete;
    this.score = score;
    this.nodesVisited = nodesVisited;
    this.elapsedMillis = elapsedMillis;
  }

  @JsonProperty("matches")
  public List<Match> matches() {
    return matches;
  }
  public boolean found() {
    return !matches.isEmpty();
  }
  @JsonProperty("status")
  public Status status() {
    return status;
  }
  @JsonProperty("optimal")
  public boolean optimal() {
    return optimal;
  }
  @JsonProperty("complete")
  public boolean complete() {
    return complete;
  }
  @JsonProperty("score")
  public double score() {
    return score;
  }
  @JsonProperty("nodesVisited")
  public long nodesVisited() {
    return nodesVisited;
  }
  @JsonProperty("elapsedMillis")
  public long elapsedMillis() {
    return elapsedMillis;
  }

  @Override
  public boolean equals(Object other) {
    if (!(other instanceof SearchResult))
      return false;
    SearchResult result = (SearchResult) other;
    return matches.equals(result.matches) && status == result.status && optimal == result.optimal
        && complete == result.complete && Double.compare(score, result.score) == 0
        && nodesVisited == result.nodesVisited && elapsedMillis == result.elapsedMillis;
  }
  @Override
  public int hashCode() {
    return Objects.hash(matches, status, optimal, complete, score, nodesVisited, elapsedMillis);
  }
  @Override
  public String toString() {
    return "SearchResult[status=" + status + ", optimal=" + optimal + ", complete=" + complete
        + ", score=" + score + ", matches=" + matches.size() + "]";
  }

  /** One deadline and cancellation state for an invocation and all its stages. */
  static final class Control {
    private final long started = System.nanoTime();
    private final long deadline;
    private final BooleanSupplier cancellation;
    private Status status = Status.Completed;

    Control(long timeoutMillis, BooleanSupplier cancellation) {
      if (timeoutMillis < 0)
        throw new IllegalArgumentException("Resolve the adaptive timeout before searching");
      long nanos = timeoutMillis >= Long.MAX_VALUE / 4 / 1_000_000L ? Long.MAX_VALUE / 4
                                                                    : timeoutMillis * 1_000_000L;
      this.deadline = started + nanos;
      this.cancellation = cancellation;
    }

    static long resolveTimeout(long timeoutMillis, int queryAtoms, int targetAtoms) {
      if (timeoutMillis < -1)
        throw new IllegalArgumentException("timeoutMillis must be -1 or non-negative");
      return timeoutMillis >= 0 ? timeoutMillis
                                : Math.min(30_000L, 500L + (long) queryAtoms * targetAtoms * 2L);
    }

    boolean stopped() {
      if (status != Status.Completed)
        return true;
      if (Thread.currentThread().isInterrupted()
          || (cancellation != null && cancellation.getAsBoolean()))
        status = Status.Cancelled;
      else if (System.nanoTime() - deadline >= 0)
        status = Status.TimedOut;
      return status != Status.Completed;
    }

    long remainingMillis() {
      if (stopped())
        return 0;
      long remaining = deadline - System.nanoTime();
      return Math.max(0L, (remaining + 999_999L) / 1_000_000L);
    }
    void limitReached() {
      if (!stopped())
        status = Status.LimitReached;
    }
    Status status() {
      return status;
    }
    long elapsedMillis() {
      return Math.max(0L, (System.nanoTime() - started) / 1_000_000L);
    }
  }
}
