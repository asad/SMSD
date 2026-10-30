/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
package com.bioinception.smsd.core;

import com.fasterxml.jackson.annotation.JsonCreator;
import com.fasterxml.jackson.annotation.JsonProperty;
import java.util.AbstractMap;
import java.util.AbstractSet;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.NoSuchElementException;
import java.util.Objects;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.function.BooleanSupplier;
import java.util.function.Function;

/** Mappings and completion information from a molecular search. */
public final class SearchResult {
  /** The calling thread searches targets in turn for about this long before other threads join. */
  private static final long HEAD_START_NANOS = 100_000L;

  /**
   * The search of each target, on up to {@code numThreads} threads (zero for one per available
   * processor), in the order of the targets. As in the native core, the calling thread first
   * searches targets in turn for about 100 microseconds, since a short batch of cheap targets ends
   * sooner than other threads can start; threads from a shared pool then join it, with the larger
   * targets first in a short batch. A failed search fails the call once the others end, with the
   * failure of the lowest failing index, the one a search in turn reports.
   */
  static List<SearchResult> searchEach(List<MolGraph> targets, int numThreads,
      Function<MolGraph, SearchResult> search) {
    Objects.requireNonNull(targets, "targets");
    if (numThreads < 0) throw new IllegalArgumentException("Thread count must be non-negative");
    List<MolGraph> inputs = new ArrayList<>(targets);
    for (MolGraph target : inputs) Objects.requireNonNull(target, "target");
    int requested = numThreads == 0 ? Runtime.getRuntime().availableProcessors() : numThreads;
    int workers = Math.max(1, Math.min(requested, inputs.size()));
    SearchResult[] results = new SearchResult[inputs.size()];
    // Every earlier target has succeeded, so a failure here is the lowest failing index.
    int done = 0;
    long started = System.nanoTime();
    while (done < results.length) {
      results[done] = search.apply(inputs.get(done));
      done++;
      if (workers > 1 && System.nanoTime() - started >= HEAD_START_NANOS) break;
    }
    if (workers > 1 && Thread.currentThread().isInterrupted())
      throw new IllegalStateException("Interrupted while searching targets");
    if (done < results.length) new Batch(inputs, results, done, workers, search).run();
    return Collections.unmodifiableList(Arrays.asList(results));
  }

  /**
   * The targets left after the head start, searched by the calling thread and by helpers from a
   * shared pool. Every index below the lowest failing one is still searched, so the reported
   * failure does not depend on the order, the number of threads or their timing.
   */
  private static final class Batch {
    /** Daemon threads shared by every batch, started on first use; idle ones end after a minute. */
    private static final class Helpers {
      static final ExecutorService POOL = Executors.newCachedThreadPool(task -> {
        Thread thread = new Thread(task, "smsd-match-many");
        thread.setDaemon(true);
        return thread;
      });
    }

    private final List<MolGraph> inputs;
    private final SearchResult[] results;
    private final int workers;
    private final Function<MolGraph, SearchResult> search;
    /** The indices left, in the order they are handed out. */
    private final int[] order;
    private final AtomicInteger next = new AtomicInteger();
    private volatile int firstFailure = Integer.MAX_VALUE;
    private Throwable failure;

    Batch(List<MolGraph> inputs, SearchResult[] results, int from, int workers,
        Function<MolGraph, SearchResult> search) {
      this.inputs = inputs;
      this.results = results;
      this.workers = workers;
      this.search = search;
      int size = results.length - from;
      order = new int[size];
      if (size <= 64L * workers) {
        // In a short batch one late, large target can leave the other threads idle, so larger
        // targets, which tend to cost more, start first, ties in input order. Long batches keep
        // input order.
        long[] keys = new long[size];
        for (int position = 0; position < size; position++) {
          int index = from + position;
          keys[position] = (long) (Integer.MAX_VALUE - inputs.get(index).n) << 32 | index;
        }
        Arrays.sort(keys);
        for (int position = 0; position < size; position++) order[position] = (int) keys[position];
      } else {
        for (int position = 0; position < size; position++) order[position] = from + position;
      }
    }

    void run() {
      int helpers = Math.min(workers, order.length) - 1;
      CountDownLatch finished = new CountDownLatch(helpers);
      List<Future<?>> started = new ArrayList<>(helpers);
      boolean ended = false;
      try {
        for (int helper = 0; helper < helpers; helper++)
          started.add(Helpers.POOL.submit(() -> {
            try {
              work();
            } finally {
              finished.countDown();
            }
          }));
        work();
        finished.await();
        ended = true;
      } catch (InterruptedException interrupted) {
        Thread.currentThread().interrupt();
        throw new IllegalStateException("Interrupted while searching targets", interrupted);
      } finally {
        if (!ended) {
          // Leave nothing to claim and stop the searches still running.
          next.set(order.length);
          for (Future<?> helper : started) helper.cancel(true);
        }
      }
      Throwable lowest;
      synchronized (this) {
        lowest = failure;
      }
      if (lowest instanceof RuntimeException) throw (RuntimeException) lowest;
      if (lowest instanceof Error) throw (Error) lowest;
      if (lowest != null) throw new IllegalStateException(lowest);
    }

    private void work() {
      for (int position; (position = next.getAndIncrement()) < order.length;) {
        int index = order[position];
        // A later index cannot be the one reported.
        if (index > firstFailure) continue;
        try {
          results[index] = search.apply(inputs.get(index));
        } catch (Throwable error) {
          fail(index, error);
        }
      }
    }

    private synchronized void fail(int index, Throwable error) {
      if (index < firstFailure) {
        firstFailure = index;
        failure = error;
      }
    }
  }

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

    private Match(AtomMap atoms, BondMap bonds, List<Set<Integer>> fragments, double score) {
      this.atomMapping = Collections.unmodifiableMap(atoms);
      this.bondMapping = Collections.unmodifiableMap(bonds);
      this.fragments = fragments;
      this.score = score;
    }

    /**
     * A match that keeps what it is given without copying it: query atoms 0 to n - 1 mapped to
     * {@code targets}, which the match then owns; each query bond, in the order of {@code bonds},
     * mapped to the bond between the images of its atoms; and {@code fragments}, an unmodifiable
     * list of unmodifiable sorted sets. Its maps iterate, compare and hash as the copied ones do.
     */
    static Match of(int[] targets, QueryBonds bonds, List<Set<Integer>> fragments, double score) {
      Bond[] images = new Bond[bonds.keys.length];
      for (int bond = 0; bond < images.length; bond++)
        images[bond] = new Bond(targets[bonds.keys[bond].first()], targets[bonds.keys[bond].second()]);
      return new Match(new AtomMap(targets), new BondMap(bonds, images), fragments, score);
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
      return MatchGraph.toGraph(query, atomMapping, bondMapping.keySet());
    }

    /** Return the common fragment using the selected bonds of the query. */
    public org.openscience.cdk.interfaces.IAtomContainer toMolecule(MolGraph query) {
      return MatchGraph.toMolecule(query, atomMapping, bondMapping.keySet());
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

    /**
     * A portable SMARTS query over the selected atoms and common bonds, retaining their source
     * context. Unsupported chemistry predicates throw {@link IllegalArgumentException}; the
     * pattern describes a subgraph query and does not encode an induced MCIS objective.
     */
    public String toSmarts(MolGraph query, ChemOptions chemistry) {
      return SmartsWriter.write(query, atomMapping.keySet(), bondMapping.keySet(), chemistry);
    }

    /** The selected common graph as SMARTS, reading atom indices from the supplied molecule. */
    public String toSmarts(org.openscience.cdk.interfaces.IAtomContainer query, ChemOptions chemistry) {
      return toSmarts(new MolGraph(query), chemistry);
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

  /** The bonds of one query in the order matches report them, shared by its matches. */
  static final class QueryBonds {
    final Bond[] keys;
    private final Map<Bond, Integer> positions;

    QueryBonds(List<Bond> bonds) {
      keys = bonds.toArray(new Bond[0]);
      positions = new HashMap<>(keys.length * 2);
      for (int bond = 0; bond < keys.length; bond++) positions.put(keys[bond], bond);
    }
  }

  /** Atoms 0 to n - 1 in ascending order, each mapped to its entry of an array. */
  private static final class AtomMap extends AbstractMap<Integer, Integer> {
    private final int[] targets;

    AtomMap(int[] targets) {
      this.targets = targets;
    }

    @Override
    public int size() {
      return targets.length;
    }
    @Override
    public boolean containsKey(Object key) {
      return key instanceof Integer && (Integer) key >= 0 && (Integer) key < targets.length;
    }
    @Override
    public Integer get(Object key) {
      return containsKey(key) ? targets[(Integer) key] : null;
    }
    @Override
    public Set<Map.Entry<Integer, Integer>> entrySet() {
      return new AbstractSet<Map.Entry<Integer, Integer>>() {
        @Override
        public int size() {
          return targets.length;
        }
        @Override
        public Iterator<Map.Entry<Integer, Integer>> iterator() {
          return new Iterator<Map.Entry<Integer, Integer>>() {
            private int next;
            @Override
            public boolean hasNext() {
              return next < targets.length;
            }
            @Override
            public Map.Entry<Integer, Integer> next() {
              if (next >= targets.length) throw new NoSuchElementException();
              int atom = next++;
              return new AbstractMap.SimpleImmutableEntry<>(atom, targets[atom]);
            }
          };
        }
      };
    }
  }

  /** Query bonds in their reported order, each mapped to the image at its position. */
  private static final class BondMap extends AbstractMap<Bond, Bond> {
    private final QueryBonds bonds;
    private final Bond[] images;

    BondMap(QueryBonds bonds, Bond[] images) {
      this.bonds = bonds;
      this.images = images;
    }

    @Override
    public int size() {
      return images.length;
    }
    @Override
    public boolean containsKey(Object key) {
      return bonds.positions.containsKey(key);
    }
    @Override
    public Bond get(Object key) {
      Integer position = bonds.positions.get(key);
      return position == null ? null : images[position];
    }
    @Override
    public Set<Map.Entry<Bond, Bond>> entrySet() {
      return new AbstractSet<Map.Entry<Bond, Bond>>() {
        @Override
        public int size() {
          return images.length;
        }
        @Override
        public Iterator<Map.Entry<Bond, Bond>> iterator() {
          return new Iterator<Map.Entry<Bond, Bond>>() {
            private int next;
            @Override
            public boolean hasNext() {
              return next < images.length;
            }
            @Override
            public Map.Entry<Bond, Bond> next() {
              if (next >= images.length) throw new NoSuchElementException();
              int bond = next++;
              return new AbstractMap.SimpleImmutableEntry<>(bonds.keys[bond], images[bond]);
            }
          };
        }
      };
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
    /** Results found in time are built until an eighth of the time limit past the deadline. */
    private final long resultsDeadline;
    private final BooleanSupplier cancellation;
    private Status status = Status.Completed;
    private int polls;

    Control(long timeoutMillis, BooleanSupplier cancellation) {
      if (timeoutMillis < 0)
        throw new IllegalArgumentException("Resolve the adaptive timeout before searching");
      long nanos = timeoutMillis >= Long.MAX_VALUE / 4 / 1_000_000L ? Long.MAX_VALUE / 4
                                                                    : timeoutMillis * 1_000_000L;
      this.deadline = started + nanos;
      this.resultsDeadline = deadline + nanos / 8;
      this.cancellation = cancellation;
    }

    /**
     * Whether building the results of a finished search should stop: on cancellation, or past
     * the results deadline, when a completed search is reported as timed out.
     */
    boolean resultsStopped() {
      Status reason = Status.Completed;
      if (Thread.currentThread().isInterrupted() || (cancellation != null && cancellation.getAsBoolean()))
        reason = Status.Cancelled;
      else if (System.nanoTime() - resultsDeadline >= 0)
        reason = Status.TimedOut;
      if (reason == Status.Completed) return false;
      if (status == Status.Completed) status = reason;
      return true;
    }

    /**
     * The time limit of a search: {@code timeoutMillis} when it is not negative, or else, for any
     * negative value, the adaptive budget of the native core, which is 10 s for a substructure
     * search and, for a common-subgraph search, 500 ms plus 2 ms per atom pair up to 30 s.
     */
    static long resolveTimeout(long timeoutMillis, int queryAtoms, int targetAtoms,
        boolean commonSubgraph) {
      if (timeoutMillis >= 0)
        return timeoutMillis;
      return commonSubgraph ? Math.min(30_000L, 500L + (long) queryAtoms * targetAtoms * 2L) : 10_000L;
    }

    /** Reject a negative node limit; any negative time limit selects the adaptive budget. */
    static void checkLimits(long timeoutMillis, long nodeLimit) {
      if (nodeLimit < 0)
        throw new IllegalArgumentException("Node limit must be non-negative");
    }

    boolean stopped() {
      return stopped(true);
    }

    /** A caller that skips the clock must still read it at a bounded interval. */
    boolean stopped(boolean readClock) {
      if (status != Status.Completed)
        return true;
      if (Thread.currentThread().isInterrupted()
          || (cancellation != null && cancellation.getAsBoolean()))
        status = Status.Cancelled;
      else if (readClock && System.nanoTime() - deadline >= 0)
        status = Status.TimedOut;
      return status != Status.Completed;
    }

    /** For short search steps: cancellation is checked on each call, the clock on every sixteenth. */
    boolean stoppedSampled() {
      return stopped((++polls & 15) == 0);
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
