/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
package com.bioinception.smsd.core;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Deque;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import org.openscience.cdk.AtomRef;
import org.openscience.cdk.BondRef;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.isomorphism.matchers.Expr;
import org.openscience.cdk.isomorphism.matchers.QueryAtom;
import org.openscience.cdk.isomorphism.matchers.QueryAtomContainer;
import org.openscience.cdk.isomorphism.matchers.QueryBond;
import org.openscience.cdk.smarts.Smarts;
import org.openscience.cdk.smarts.SmartsResult;
import org.openscience.cdk.silent.SilentChemObjectBuilder;

/**
 * Parses SMARTS within the limits shared with the native core, as {@link Substructure#smarts}
 * lists them, measuring the text as the native parser does before the CDK parser reads it. The
 * CDK parser recurses once per operator, so long text is parsed on a thread with a large stack,
 * and long chains of one operator are then rebalanced so that matching recurses only
 * logarithmically.
 */
final class SmartsInput {
  /** Text up to this length parses on the calling thread. */
  private static final int INLINE_LENGTH = 1000;
  /** Stack of the parsing thread for longer text; it is reserved, and used only as deep as the text needs. */
  private static final long PARSER_STACK_BYTES = 256L << 20;
  /** Chains of one operator with more operands are rebalanced; written patterns have far fewer. */
  private static final int BALANCED_OPERANDS = 32;

  private SmartsInput() {}

  /** The compiled query, or IllegalArgumentException for invalid SMARTS or SMARTS beyond the limits. */
  static IAtomContainer parse(String expression) {
    if (expression.length() > Substructure.MAX_SMARTS_LENGTH)
      throw new IllegalArgumentException("SMARTS exceeds the maximum length of "
          + Substructure.MAX_SMARTS_LENGTH + " characters");
    // The parser stores atoms, and joins primitives written side by side, in time quadratic in
    // their number, so the text is measured first as the native parser counts it.
    if (atoms(expression) > Substructure.MAX_SMARTS_ATOMS) throw tooManyAtoms();
    if (recursionDepth(expression) > Substructure.MAX_SMARTS_RECURSION_DEPTH)
      throw new IllegalArgumentException("Recursive SMARTS nesting exceeds maximum depth of "
          + Substructure.MAX_SMARTS_RECURSION_DEPTH);
    NativeCount count = NativeCount.of(expression);
    if (count != null) {
      if (count.atoms > Substructure.MAX_SMARTS_ATOMS) throw tooManyAtoms();
      if (count.largestExpression > Substructure.MAX_SMARTS_NODES_PER_EXPRESSION) throw tooLargeExpression();
      if (count.nodes > Substructure.MAX_SMARTS_EXPRESSION_NODES) throw tooManyNodes();
    } else if (longestBracket(expression) > Substructure.MAX_SMARTS_NODES_PER_EXPRESSION / 2) {
      // Outside the native grammar, a bracket atom of more characters may hold more nodes.
      throw tooLargeExpression();
    }
    final String read = directionsAtOpenings(expression);
    final boolean counted = count != null;
    if (read.length() <= INLINE_LENGTH) return compile(read, counted);
    Object[] outcome = new Object[1];
    Thread parser = new Thread(null, () -> {
      try {
        outcome[0] = compile(read, counted);
      } catch (Throwable failure) {
        outcome[0] = failure;
      }
    }, "smsd-smarts-parser", PARSER_STACK_BYTES);
    parser.setDaemon(true);
    parser.start();
    boolean interrupted = false;
    while (true) {
      try {
        parser.join();
        break;
      } catch (InterruptedException interruption) {
        // Parsing takes time linear in the bounded text, so it is awaited.
        interrupted = true;
      }
    }
    if (interrupted) Thread.currentThread().interrupt();
    if (outcome[0] instanceof RuntimeException) throw (RuntimeException) outcome[0];
    if (outcome[0] instanceof Error) throw (Error) outcome[0];
    return (IAtomContainer) outcome[0];
  }

  /** Parses text already measured; text outside the native grammar is counted in the parsed query. */
  private static IAtomContainer compile(String expression, boolean counted) {
    IAtomContainer query = new QueryAtomContainer(SilentChemObjectBuilder.getInstance());
    SmartsResult parsed = Smarts.parseToResult(query, expression);
    if (!parsed.ok()) throw new IllegalArgumentException("Invalid SMARTS: " + parsed.getMessage());
    if (!counted) {
      long[] nodes = nodes(query);
      if (nodes[1] > Substructure.MAX_SMARTS_NODES_PER_EXPRESSION) throw tooLargeExpression();
      if (nodes[0] > Substructure.MAX_SMARTS_EXPRESSION_NODES) throw tooManyNodes();
    }
    normalise(query);
    return query;
  }

  private static IllegalArgumentException tooManyAtoms() {
    return new IllegalArgumentException("SMARTS exceeds the maximum of " + Substructure.MAX_SMARTS_ATOMS + " atoms");
  }

  private static IllegalArgumentException tooManyNodes() {
    return new IllegalArgumentException("SMARTS exceeds the maximum of "
        + Substructure.MAX_SMARTS_EXPRESSION_NODES + " expression nodes");
  }

  private static IllegalArgumentException tooLargeExpression() {
    return new IllegalArgumentException("SMARTS atom or bond expression exceeds the maximum of "
        + Substructure.MAX_SMARTS_NODES_PER_EXPRESSION + " nodes");
  }

  /** Deepest nesting of recursive SMARTS, read as the native parser reads it. */
  static int recursionDepth(String expression) {
    int parentheses = 0, deepest = 0;
    Deque<Integer> frames = new ArrayDeque<>();
    for (int at = 0; at < expression.length(); at++) {
      char c = expression.charAt(at);
      if (c == '(') {
        parentheses++;
        if (at > 0 && expression.charAt(at - 1) == '$') {
          frames.push(parentheses);
          deepest = Math.max(deepest, frames.size());
        }
      } else if (c == ')') {
        if (!frames.isEmpty() && frames.peek() == parentheses) frames.pop();
        if (parentheses > 0) parentheses--;
      }
    }
    return deepest;
  }

  /** Characters of the longest bracket atom, not counting its recursive SMARTS. */
  static int longestBracket(String expression) {
    int longest = 0, length = 0;
    Deque<int[]> levels = new ArrayDeque<>();
    int branches = 0;
    boolean inBracket = false;
    for (int at = 0; at < expression.length(); at++) {
      char c = expression.charAt(at);
      if (inBracket) {
        if (c == '$' && at + 1 < expression.length() && expression.charAt(at + 1) == '(') {
          levels.push(new int[] {branches, length});
          branches = 0;
          inBracket = false;
          at++;
        } else if (c == ']') {
          inBracket = false;
          longest = Math.max(longest, length);
        } else {
          length++;
        }
        continue;
      }
      if (c == '[') {
        inBracket = true;
        length = 0;
      } else if (c == '(') {
        branches++;
      } else if (c == ')') {
        if (branches > 0) branches--;
        else if (!levels.isEmpty()) {
          int[] level = levels.pop();
          branches = level[0];
          length = level[1] + 1;
          inBracket = true;
        }
      }
    }
    return Math.max(longest, inBracket ? length : 0);
  }

  /**
   * Atoms written in a SMARTS, at every depth of recursive SMARTS: a bracket atom, a wildcard,
   * or an element symbol outside brackets, where Cl and Br are one atom each.
   */
  static long atoms(String expression) {
    long atoms = 0;
    // Per open recursive SMARTS: its branch depth, so that the closing parenthesis is found.
    Deque<int[]> levels = new ArrayDeque<>();
    int branches = 0;
    boolean inBracket = false;
    for (int at = 0; at < expression.length(); at++) {
      char c = expression.charAt(at);
      if (inBracket) {
        if (c == '$' && at + 1 < expression.length() && expression.charAt(at + 1) == '(') {
          levels.push(new int[] {branches});
          branches = 0;
          inBracket = false;
          at++;
        } else if (c == ']') {
          inBracket = false;
        }
        continue;
      }
      if (c == '[') {
        atoms++;
        inBracket = true;
      } else if (c == '(') {
        branches++;
      } else if (c == ')') {
        if (branches > 0) branches--;
        else if (!levels.isEmpty()) {
          branches = levels.pop()[0];
          inBracket = true;
        }
      } else if (c == '*') {
        atoms++;
      } else if ((c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z')) {
        char previous = at > 0 ? expression.charAt(at - 1) : 0;
        if (!((c == 'l' && previous == 'C') || (c == 'r' && previous == 'B'))) atoms++;
      }
    }
    return atoms;
  }

  /**
   * The expression with the bond of each ring closure written with a direction (/ or \) at its
   * closing digit alone moved to its opening digit and reversed, since a direction at the closing
   * digit is read from the closing atom, which the CDK parser does not do: F/C=C1.Cl/1 is read as
   * F/C=C\1.Cl1. A ring closure, at any depth of recursive SMARTS, written with a bond at both ends
   * and a direction at either is rejected, as the native parser rejects it.
   */
  static String directionsAtOpenings(String expression) {
    // Per open recursive SMARTS: its branch depth and its open ring closures, each the bond written
    // at its opening digit and where that digit starts.
    Deque<Object[]> levels = new ArrayDeque<>();
    Map<Integer, Object[]> open = new HashMap<>();
    StringBuilder bond = new StringBuilder();
    // Each move: where the reversed bond goes, and where the closing bond starts and ends.
    List<int[]> moves = new ArrayList<>();
    int branches = 0;
    boolean inBracket = false;
    for (int at = 0; at < expression.length(); at++) {
      char c = expression.charAt(at);
      if (inBracket) {
        if (c == '$' && at + 1 < expression.length() && expression.charAt(at + 1) == '(') {
          levels.push(new Object[] {branches, open});
          branches = 0;
          open = new HashMap<>();
          bond.setLength(0);
          inBracket = false;
          at++;
        } else if (c == ']') {
          inBracket = false;
        }
        continue;
      }
      int ring = -1, digit = at;
      if (c >= '0' && c <= '9') {
        ring = c - '0';
      } else if (c == '%') {
        int end = at + 1;
        if (end < expression.length() && expression.charAt(end) == '(') {
          int close = expression.indexOf(')', end);
          if (close < 0) return expression;
          String digits = expression.substring(end + 1, close);
          if (digits.isEmpty() || digits.length() > 9 || !digits.chars().allMatch(Character::isDigit)) return expression;
          ring = Integer.parseInt(digits);
          at = close;
        } else if (end + 1 < expression.length() && Character.isDigit(expression.charAt(end))
            && Character.isDigit(expression.charAt(end + 1))) {
          ring = Integer.parseInt(expression.substring(end, end + 2));
          at = end + 1;
        } else {
          return expression;  // the parser rejects it
        }
      }
      if (ring >= 0) {
        String written = bond.toString();
        bond.setLength(0);
        Object[] opening = open.remove(ring);
        String first = opening == null ? null : (String) opening[0];
        if (first == null) {
          open.put(ring, new Object[] {written, digit});
        } else if (!first.isEmpty() && !written.isEmpty() && (directed(first) || directed(written))) {
          throw new IllegalArgumentException("Ring closure " + ring + " has a bond direction written at both ends");
        } else if (first.isEmpty() && directed(written)) {
          moves.add(new int[] {(Integer) opening[1], digit - written.length(), digit});
        }
      } else if ("-=#:~@!/\\?,;&".indexOf(c) >= 0) {
        bond.append(c);
      } else if (c == ')' && branches == 0 && !levels.isEmpty()) {
        Object[] level = levels.pop();
        branches = (Integer) level[0];
        @SuppressWarnings("unchecked")
        Map<Integer, Object[]> outer = (Map<Integer, Object[]>) level[1];
        open = outer;
        bond.setLength(0);
        inBracket = true;
      } else {
        if (c == '[') inBracket = true;
        else if (c == '(') branches++;
        else if (c == ')' && branches > 0) branches--;
        bond.setLength(0);
      }
    }
    if (moves.isEmpty()) return expression;
    Map<Integer, String> insertions = new HashMap<>();
    Map<Integer, Integer> deletions = new HashMap<>();
    for (int[] move : moves) {
      StringBuilder reversed = new StringBuilder(expression.substring(move[1], move[2]));
      for (int at = 0; at < reversed.length(); at++) {
        char c = reversed.charAt(at);
        reversed.setCharAt(at, c == '/' ? '\\' : c == '\\' ? '/' : c);
      }
      insertions.put(move[0], reversed.toString());
      deletions.put(move[1], move[2]);
    }
    StringBuilder text = new StringBuilder(expression.length());
    for (int at = 0; at < expression.length(); ) {
      String inserted = insertions.get(at);
      if (inserted != null) text.append(inserted);
      Integer end = deletions.get(at);
      if (end != null) {
        at = end;
      } else {
        text.append(expression.charAt(at++));
      }
    }
    return text.toString();
  }

  private static boolean directed(String bond) {
    return bond.indexOf('/') >= 0 || bond.indexOf('\\') >= 0;
  }

  /**
   * Expression nodes of a parsed query at every depth of recursive SMARTS, and those of its
   * largest atom or bond expression: each primitive and each logical operator counts once,
   * read without recursion.
   */
  static long[] nodes(IAtomContainer query) {
    long total = 0, largest = 0;
    Deque<IAtomContainer> containers = new ArrayDeque<>();
    containers.push(query);
    Deque<Expr> pending = new ArrayDeque<>();
    while (!containers.isEmpty()) {
      IAtomContainer container = containers.pop();
      List<Expr> roots = new ArrayList<>();
      for (IAtom atom : container.atoms()) {
        IAtom source = AtomRef.deref(atom);
        if (source instanceof QueryAtom && ((QueryAtom) source).getExpression() != null)
          roots.add(((QueryAtom) source).getExpression());
      }
      for (IBond bond : container.bonds()) {
        IBond source = BondRef.deref(bond);
        if (source instanceof QueryBond && ((QueryBond) source).getExpression() != null)
          roots.add(((QueryBond) source).getExpression());
      }
      for (Expr root : roots) {
        long nodes = 0;
        pending.push(root);
        while (!pending.isEmpty()) {
          Expr expression = pending.pop();
          nodes++;
          switch (expression.type()) {
            case AND:
            case OR:
              if (expression.right() != null) pending.push(expression.right());
              if (expression.left() != null) pending.push(expression.left());
              break;
            case NOT:
              if (expression.left() != null) pending.push(expression.left());
              break;
            case RECURSIVE:
              if (expression.subquery() != null) containers.push(expression.subquery());
              break;
            default:
              break;
          }
        }
        total += nodes;
        largest = Math.max(largest, nodes);
      }
    }
    return new long[] {total, largest};
  }

  /**
   * Rebalances every chain of one Boolean operator with more than BALANCED_OPERANDS operands,
   * keeping the operands in order, so a long chain nests logarithmically; shorter chains are
   * left as parsed. Evaluation reads the operands in the same order, so matches are unchanged.
   *
   * <p>Reads each directional bond, / or \, as a single, non-aromatic bond, as Daylight SMARTS
   * and the native parser do. CDK parses it as a bond of order one, which an aromatic bond of a
   * Kekule structure also has; the direction itself is kept on the double bond's stereo element.
   */
  private static void normalise(IAtomContainer query) {
    Deque<Object> pending = new ArrayDeque<>();
    pending.push(query);
    while (!pending.isEmpty()) {
      Object item = pending.pop();
      if (item instanceof IAtomContainer) {
        IAtomContainer container = (IAtomContainer) item;
        for (IAtom atom : container.atoms()) {
          IAtom source = AtomRef.deref(atom);
          if (source instanceof QueryAtom && ((QueryAtom) source).getExpression() != null)
            pending.push(((QueryAtom) source).getExpression());
        }
        for (IBond bond : container.bonds()) {
          IBond source = BondRef.deref(bond);
          if (source instanceof QueryBond && ((QueryBond) source).getExpression() != null)
            pending.push(((QueryBond) source).getExpression());
        }
        continue;
      }
      Expr expression = (Expr) item;
      switch (expression.type()) {
        case AND:
        case OR: {
          List<Expr> operands = operands(expression);
          if (operands.size() > BALANCED_OPERANDS) {
            Expr balanced = build(expression.type(), operands, 0, operands.size());
            expression.setLogical(expression.type(), balanced.left(), balanced.right());
          }
          for (Expr operand : operands) pending.push(operand);
          break;
        }
        case NOT:
          if (expression.left() != null) pending.push(expression.left());
          break;
        case RECURSIVE:
          if (expression.subquery() != null) pending.push(expression.subquery());
          break;
        case ORDER:
          // Only a directional bond parses to a plain order test.
          if (expression.value() == 1) expression.setPrimitive(Expr.Type.ALIPHATIC_ORDER, 1);
          break;
        default:
          break;
      }
    }
  }

  /**
   * Atoms and expression nodes of a SMARTS as the native parser counts them, by the native
   * grammar: each primitive counts one node, an operator joining k operands k - 1, and an odd
   * run of negations one; a recursive SMARTS counts one node in its expression and is counted
   * as a pattern of its own. Nested patterns are read from a work list, without recursion.
   */
  static final class NativeCount {
    /** The element symbols, "*" and the aromatic symbols, as the native grammar reads them. */
    private static final Set<String> ELEMENTS = elementSymbols();

    private static Set<String> elementSymbols() {
      Set<String> symbols = new HashSet<>(Arrays.asList("*", "b", "c", "n", "o", "p", "s", "se", "te"));
      for (int element = 1; element < SmilesWriter.ELEMENT_SYMBOL.length; element++)
        symbols.add(SmilesWriter.ELEMENT_SYMBOL[element]);
      return symbols;
    }

    long atoms, nodes, largestExpression;
    private String text;
    private int at;
    private final Deque<String> nested = new ArrayDeque<>();

    /** The count, or null when the text is outside the native grammar. */
    static NativeCount of(String smarts) {
      NativeCount count = new NativeCount();
      count.nested.push(smarts);
      try {
        while (!count.nested.isEmpty()) count.pattern(count.nested.pop());
      } catch (IllegalStateException outside) {
        return null;
      }
      return count;
    }

    private static IllegalStateException outside() {
      return new IllegalStateException();
    }

    private char peek() { return at < text.length() ? text.charAt(at) : 0; }
    private static boolean upper(char c) { return c >= 'A' && c <= 'Z'; }
    private static boolean lower(char c) { return c >= 'a' && c <= 'z'; }
    private static boolean letter(char c) { return upper(c) || lower(c); }
    private static boolean digit(char c) { return c >= '0' && c <= '9'; }
    private static boolean space(char c) { return c == ' ' || (c >= 9 && c <= 13); }
    private static boolean bondChar(char c) {
      return c == '-' || c == '=' || c == '#' || c == ':' || c == '~' || c == '@' || c == '!' || c == '/' || c == '\\';
    }

    private void expression(long size) {
      nodes += size;
      largestExpression = Math.max(largestExpression, size);
    }

    /** Reads a non-negative int, or returns -1 when no digits follow; larger values are outside the grammar. */
    private int digits() {
      int start = at;
      long value = 0;
      while (digit(peek())) {
        value = value * 10 + (peek() - '0');
        if (value > Integer.MAX_VALUE) throw outside();
        at++;
      }
      return at == start ? -1 : (int) value;
    }

    private void optionalRange() {
      if (peek() != '{') return;
      int saved = at;
      at++;
      digits();
      if (peek() == '-') {
        at++;
        digits();
      }
      if (peek() == '}') at++;
      else at = saved;
    }

    private void pattern(String smarts) {
      text = smarts;
      at = 0;
      Deque<Boolean> branches = new ArrayDeque<>();
      Map<Integer, Long> closures = new HashMap<>();
      boolean haveAtom = false, pending = false;
      long pendingNodes = 0;
      while (at < text.length()) {
        char c = peek();
        if (c == '(') {
          at++;
          branches.push(haveAtom);
          pending = false;
        } else if (c == ')') {
          at++;
          if (branches.isEmpty()) throw outside();
          haveAtom = branches.pop();
          pending = false;
        } else if (bondChar(c)) {
          pendingNodes = bondExpression();
          pending = true;
        } else if (digit(c) || c == '%') {
          int ring;
          if (c == '%') {
            at++;
            ring = digits();
            if (ring < 0) throw outside();
          } else {
            ring = text.charAt(at++) - '0';
          }
          if (!haveAtom) throw outside();
          Long open = closures.remove(ring);
          if (open != null) {
            // The closing bond takes the bond written at the closure, else the one at the opening.
            expression(pending ? pendingNodes : open);
          } else {
            closures.put(ring, pending ? pendingNodes : 1L);
          }
          pending = false;
        } else if (c == '[' || c == '*' || letter(c)) {
          if (c == '[') {
            at++;
            expression(andLow());
            if (peek() != ']') throw outside();
            at++;
          } else if (c == '*' || c == 'a' || c == 'A') {
            at++;
            expression(1);
          } else {
            at++;
            String symbol = String.valueOf(c);
            if ((c == 'C' && peek() == 'l') || (c == 'B' && peek() == 'r')) symbol += text.charAt(at++);
            if (!ELEMENTS.contains(symbol)) throw outside();
            expression(1);
          }
          atoms++;
          // A bond written before any atom is kept for the next bond, as the native parser keeps it.
          if (haveAtom) {
            expression(pending ? pendingNodes : 1);
            pending = false;
          }
          haveAtom = true;
        } else if (space(c)) {
          at++;
        } else if (c == '.') {
          at++;
          haveAtom = false;
          pending = false;
        } else {
          throw outside();
        }
      }
      if (!branches.isEmpty() || !closures.isEmpty()) throw outside();
    }

    private long bondExpression() {
      long first = bondOr();
      if (peek() != ';') return first;
      long total = first;
      while (peek() == ';') {
        at++;
        total += bondOr() + 1;
      }
      return total;
    }

    private long bondOr() {
      long first = bondAndHigh();
      if (peek() != ',') return first;
      long total = first;
      while (peek() == ',') {
        at++;
        total += bondAndHigh() + 1;
      }
      return total;
    }

    private long bondAndHigh() {
      long total = singleBond();
      while (peek() == '&' || bondChar(peek())) {
        if (peek() == '&') at++;
        total += singleBond() + 1;
      }
      return total;
    }

    private long singleBond() {
      boolean negated = false;
      while (peek() == '!') {
        at++;
        negated = !negated;
      }
      char c = peek();
      if (!(c == '-' || c == '=' || c == '#' || c == ':' || c == '~' || c == '@' || c == '/' || c == '\\'))
        throw outside();
      at++;
      return negated ? 2 : 1;
    }

    private long andLow() {
      long total = or();
      while (peek() == ';') {
        at++;
        total += or() + 1;
      }
      return total;
    }

    private long or() {
      long total = andHigh();
      while (peek() == ',') {
        at++;
        total += andHigh() + 1;
      }
      return total;
    }

    private long andHigh() {
      long total = not();
      while (at < text.length()) {
        char c = peek();
        if (c == ',' || c == ';' || c == ']') break;
        if (c == '&') at++;
        total += not() + 1;
      }
      return total;
    }

    private long not() {
      boolean negated = false;
      while (peek() == '!') {
        at++;
        negated = !negated;
      }
      primitive();
      return negated ? 2 : 1;
    }

    private boolean twoLetterElement() {
      return at + 1 < text.length() && lower(text.charAt(at + 1))
          && ELEMENTS.contains(text.substring(at, at + 2));
    }

    private void element() {
      char c = peek();
      if (!letter(c)) throw outside();
      at++;
      String symbol = String.valueOf(c);
      if (lower(peek()) && ELEMENTS.contains(symbol + peek())) symbol += text.charAt(at++);
      if (!ELEMENTS.contains(symbol)) throw outside();
    }

    private void primitive() {
      char c = peek();
      switch (c) {
        case '#':
          at++;
          if (digits() < 0) throw outside();
          return;
        case '*':
        case 'a':
          at++;
          return;
        case 'A':
          if (twoLetterElement()) element();
          else at++;
          return;
        case 'X':
        case 'Q':
        case 'D':
        case 'R':
        case 'Z':
          if (c != 'Q' && twoLetterElement()) {
            element();
            return;
          }
          at++;
          digits();
          optionalRange();
          return;
        case 'd':
        case 'v':
        case 'r':
        case 'x':
        case 'z':
          at++;
          digits();
          optionalRange();
          return;
        case '^': {
          at++;
          int hybridisation = digits();
          if (hybridisation < 0 || hybridisation > 5) throw outside();
          return;
        }
        case '+':
        case '-':
          at++;
          if (digits() < 0) while (peek() == c) at++;
          return;
        case 'H':
          if (at + 1 < text.length() && "fgos".indexOf(text.charAt(at + 1)) >= 0) {
            element();
            return;
          }
          at++;
          digits();
          optionalRange();
          return;
        case '@':
          at++;
          if (peek() == '@' || peek() == '?') at++;
          return;
        case ':':
          at++;
          digits();
          return;
        case '$': {
          at++;
          if (peek() != '(') throw outside();
          at++;
          int depth = 1, start = at;
          while (at < text.length() && depth > 0) {
            if (text.charAt(at) == '(') depth++;
            else if (text.charAt(at) == ')') depth--;
            if (depth > 0) at++;
          }
          if (depth != 0) throw outside();
          nested.push(text.substring(start, at));
          at++;
          return;
        }
        default:
          if (letter(c)) element();
          else if (digit(c)) digits();
          else throw outside();
      }
    }
  }

  /** Operands of the chain of one operator rooted at an expression, in evaluation order. */
  private static List<Expr> operands(Expr root) {
    List<Expr> operands = new ArrayList<>();
    Deque<Expr> pending = new ArrayDeque<>();
    pending.push(root);
    while (!pending.isEmpty()) {
      Expr expression = pending.pop();
      if (expression.type() == root.type() && expression.left() != null && expression.right() != null) {
        pending.push(expression.right());
        pending.push(expression.left());
      } else {
        operands.add(expression);
      }
    }
    return operands;
  }

  /** A tree of one operator over operands[from, to), of depth logarithmic in their number. */
  private static Expr build(Expr.Type type, List<Expr> operands, int from, int to) {
    if (to - from == 1) return operands.get(from);
    int middle = (from + to) >>> 1;
    return new Expr(type, build(type, operands, from, middle), build(type, operands, middle, to));
  }
}
