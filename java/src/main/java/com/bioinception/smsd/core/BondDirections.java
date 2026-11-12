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
import java.util.List;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.smarts.SmartsResult;

/**
 * The directional bonds, / and \, of a SMARTS query, read as the native core reads them. The CDK
 * parser keeps one pair of directional bonds per double bond, as a relation on that bond, and
 * none around an aromatic bond, so the expression of each directional bond is read again from the
 * text. A bond matches as the native bond test does, and a complete mapping holds when the
 * directions can be chosen consistently with every double bond and aromatic bond between them.
 */
final class BondDirections {
  private static final int ANY = 0, SINGLE = 1, DOUBLE = 2, TRIPLE = 3, AROMATIC = 4, RING = 5,
      UP = 6, DOWN = 7, UP_OR_UNSPECIFIED = 8, DOWN_OR_UNSPECIFIED = 9, NOT = 10, AND = 11, OR = 12;
  /** Results of a test that may not know a bond's direction yet. */
  private static final int FALSE = 0, TRUE = 1, UNKNOWN = 2;
  /** Bits of the directions a bond may take: up, down and none, beside a double bond without a configuration. */
  private static final int BOND_UP = 1, BOND_DOWN = 2, BOND_UNSPECIFIED = 4;

  /** A bond expression: a primitive, or a logical operator over its operands. */
  private static final class Node {
    final int type;
    final Node[] operands;

    Node(int type, Node... operands) {
      this.type = type;
      this.operands = operands;
    }
  }

  private final int[] from, to;
  /** The expression of each directional bond, and null for the other bonds. */
  private final Node[] expressions;
  /** The neighbours of each query atom and the bonds to them, in bond order. */
  private final int[][] neighbors, bondsAt;

  private BondDirections(int atoms, int[] from, int[] to, Node[] expressions) {
    this.from = from;
    this.to = to;
    this.expressions = expressions;
    int[] degree = new int[atoms];
    for (int bond = 0; bond < from.length; bond++) {
      degree[from[bond]]++;
      degree[to[bond]]++;
    }
    neighbors = new int[atoms][];
    bondsAt = new int[atoms][];
    for (int atom = 0; atom < atoms; atom++) {
      neighbors[atom] = new int[degree[atom]];
      bondsAt[atom] = new int[degree[atom]];
    }
    Arrays.fill(degree, 0);
    for (int bond = 0; bond < from.length; bond++) {
      int first = from[bond], second = to[bond];
      neighbors[first][degree[first]] = second;
      bondsAt[first][degree[first]++] = bond;
      neighbors[second][degree[second]] = first;
      bondsAt[second][degree[second]++] = bond;
    }
  }

  /**
   * The directional bonds of a parsed query, or null when it has none. {@code text} is the text the
   * parser read, and {@code parsed} gives where each bond is written in it; a bond runs from the
   * atom written first, or from the opening atom of a ring closure, as in the native parser.
   */
  static BondDirections read(IAtomContainer query, String text, SmartsResult parsed) {
    int bonds = query.getBondCount();
    int[] from = new int[bonds], to = new int[bonds];
    Node[] expressions = new Node[bonds];
    boolean directional = false;
    for (int index = 0; index < bonds; index++) {
      IBond bond = query.getBond(index);
      from[index] = query.indexOf(bond.getBegin());
      to[index] = query.indexOf(bond.getEnd());
      int at = parsed.getBondLocation(index);
      if (at < 0 || at >= text.length() || !bondCharacter(text.charAt(at))) continue;
      Node expression = new Parser(text, at).expression();
      if (!directional(expression)) continue;
      expressions[index] = expression;
      directional = true;
    }
    return directional ? new BondDirections(query.getAtomCount(), from, to, expressions) : null;
  }

  /** Whether the bond has a directional primitive. */
  boolean directional(int bond) {
    return expressions[bond] != null;
  }

  /** Whether the directional bond can match the target bond, taking either direction. */
  boolean compatible(int bond, IBond target) {
    Node expression = expressions[bond];
    return withDirection(expression, target, 1)
        || withDirection(expression, target, -1)
        || withoutDirection(expression, target) == TRUE
        || order(target) != 1 && withDirection(expression, target, 0);
  }

  /**
   * Whether the directions of a complete mapping are consistent. A pair of directional bonds
   * constrains the bond between their atoms: a double bond takes the configuration they give it,
   * as in Daylight SMARTS, and an aromatic bond, or a double bond without a configuration, matches
   * only when either bond is written with ?. Bonds are read from {@code molecule}, the target as
   * the SMARTS atoms and bonds read it, and configurations from {@code graph}, with the same atoms.
   */
  boolean holds(int[] mapping, IAtomContainer molecule, MolGraph graph) {
    int bonds = from.length;
    int[] variableForBond = new int[bonds];
    Arrays.fill(variableForBond, -1);
    int[] allowedDirections = new int[bonds];
    int variables = 0;
    for (int bond = 0; bond < bonds; bond++) {
      Node expression = expressions[bond];
      if (expression == null) continue;
      IBond target = bond(molecule, mapping[from[bond]], mapping[to[bond]]);
      if (!carrier(target)) continue;
      // Away from double bonds, an aromatic bond the expression accepts without a direction, as
      // with /,: or !/, stays free; at a double bond it carries one, like a single bond.
      if (order(target) != 1 && !atDoubleBond(from[bond], mapping, molecule)
          && !atDoubleBond(to[bond], mapping, molecule) && withDirection(expression, target, 0)) continue;
      int allowed = 0;
      if (withDirection(expression, target, 1)) allowed |= BOND_UP;
      if (withDirection(expression, target, -1)) allowed |= BOND_DOWN;
      if (withoutDirection(expression, target) == TRUE) allowed |= BOND_UNSPECIFIED;
      if (allowed == 0) return false;
      variableForBond[bond] = variables;
      allowedDirections[variables++] = allowed;
    }
    if (variables == 0) return true;
    List<List<int[]>> constraints = new ArrayList<>(variables);
    boolean[] specified = new boolean[variables];
    for (int variable = 0; variable < variables; variable++) {
      constraints.add(new ArrayList<int[]>());
      specified[variable] = (allowedDirections[variable] & BOND_UNSPECIFIED) == 0;
    }
    // A carrier of a configured double bond takes its direction from that bond, so it says nothing
    // about an aromatic bond beside it; such carriers are left out of the aromatic foci below.
    boolean[] fixedByDoubleBond = new boolean[variables];
    for (int focus = 0; focus < bonds; focus++) {
      int mappedFirst = mapping[from[focus]], mappedSecond = mapping[to[focus]];
      if (order(bond(molecule, mappedFirst, mappedSecond)) != 2 || expressions[focus] != null
          || graph.dbStereo(mappedFirst, mappedSecond) == 0) continue;
      int[] first = carriers(from[focus], to[focus], variableForBond);
      int[] second = carriers(to[focus], from[focus], variableForBond);
      if (first.length == 0 || second.length == 0) continue;
      for (int bond : first) fixedByDoubleBond[variableForBond[bond]] = true;
      for (int bond : second) fixedByDoubleBond[variableForBond[bond]] = true;
    }
    List<int[]> optionalConstraints = new ArrayList<>();
    List<int[]> unspecifiedFoci = new ArrayList<>();
    for (int focus = 0; focus < bonds; focus++) {
      int first = from[focus], second = to[focus];
      int mappedFirst = mapping[first], mappedSecond = mapping[second];
      int order = order(bond(molecule, mappedFirst, mappedSecond));
      boolean aromatic = order == 4;
      if (order != 2 && !aromatic || expressions[focus] != null) continue;
      int[][] carriers = {carriers(first, second, variableForBond), carriers(second, first, variableForBond)};
      if (aromatic)
        for (int end = 0; end < 2; end++) carriers[end] = unfixed(carriers[end], fixedByDoubleBond, variableForBond);
      int configuration = aromatic ? 0 : graph.dbStereo(mappedFirst, mappedSecond);
      boolean configured = configuration != 0 && carriers[0].length > 0 && carriers[1].length > 0;
      for (int end = 0; end < 2; end++) {
        int atom = end == 0 ? first : second;
        if (configured) for (int bond : carriers[end]) specified[variableForBond[bond]] = true;
        // Bonds on one side of a double bond point opposite ways; an aromatic bond has no sides.
        if (aromatic) continue;
        for (int previous = 0; previous < carriers[end].length; previous++) {
          int firstBond = carriers[end][previous];
          for (int next = previous + 1; next < carriers[end].length; next++) {
            int secondBond = carriers[end][next];
            int relation = -outward(firstBond, atom) * outward(secondBond, atom);
            int[] constraint = {variableForBond[firstBond], variableForBond[secondBond], relation};
            if (configured) constrain(constraints, constraint);
            else optionalConstraints.add(constraint);
          }
        }
      }
      if (carriers[0].length == 0 || carriers[1].length == 0) continue;
      if (configuration == 0) {
        int[] focusVariables = new int[carriers[0].length + carriers[1].length];
        int count = 0;
        for (int[] end : carriers) for (int bond : end) focusVariables[count++] = variableForBond[bond];
        unspecifiedFoci.add(focusVariables);
        continue;
      }
      int firstReference = graph.doubleBondReference(mappedFirst, mappedSecond);
      int secondReference = graph.doubleBondReference(mappedSecond, mappedFirst);
      for (int firstBond : carriers[0]) {
        int firstLigand = from[firstBond] == first ? to[firstBond] : from[firstBond];
        for (int secondBond : carriers[1]) {
          int secondLigand = from[secondBond] == second ? to[secondBond] : from[secondBond];
          int relation = (configuration == 1 ? 1 : -1)
              * (mapping[firstLigand] == firstReference ? 1 : -1)
              * (mapping[secondLigand] == secondReference ? 1 : -1)
              * outward(firstBond, first) * outward(secondBond, second);
          constrain(constraints, new int[] {variableForBond[firstBond], variableForBond[secondBond], relation});
        }
      }
    }
    // An optional carrier not required by a configured focus can take no direction. This only
    // removes same-end relations; every unconfigured focus must retain such a carrier.
    for (int[] focusVariables : unspecifiedFoci) {
      boolean unspecified = false;
      for (int variable : focusVariables) unspecified |= !specified[variable];
      if (!unspecified) return false;
    }
    for (int[] constraint : optionalConstraints)
      if (specified[constraint[0]] && specified[constraint[1]]) constrain(constraints, constraint);
    return consistent(constraints, allowedDirections, specified);
  }

  /** Whether the specified bonds can take directions that satisfy every relation between them. */
  private static boolean consistent(List<List<int[]>> constraints, int[] allowedDirections, boolean[] specified) {
    int[] relativeDirections = new int[specified.length];
    Deque<Integer> pending = new ArrayDeque<>();
    for (int start = 0; start < specified.length; start++) {
      if (relativeDirections[start] != 0 || !specified[start]) continue;
      relativeDirections[start] = 1;
      // The directions the start bond may take that every bond reached from it allows.
      int allowedRoot = BOND_UP | BOND_DOWN;
      pending.add(start);
      while (!pending.isEmpty()) {
        int current = pending.poll();
        int allowed = allowedDirections[current] & (BOND_UP | BOND_DOWN);
        allowedRoot &= relativeDirections[current] == 1 ? allowed
            : (allowed & BOND_UP) << 1 | (allowed & BOND_DOWN) >> 1;
        if (allowedRoot == 0) return false;
        for (int[] constraint : constraints.get(current)) {
          int other = constraint[0], required = relativeDirections[current] * constraint[1];
          if (relativeDirections[other] == 0) {
            relativeDirections[other] = required;
            pending.add(other);
          } else if (relativeDirections[other] != required) {
            return false;
          }
        }
      }
    }
    return true;
  }

  /** Record a relation, {@code {first, second, sign}}, at both of its bonds. */
  private static void constrain(List<List<int[]>> constraints, int[] constraint) {
    constraints.get(constraint[0]).add(new int[] {constraint[1], constraint[2]});
    constraints.get(constraint[1]).add(new int[] {constraint[0], constraint[2]});
  }

  /** The directional bonds with a variable at {@code atom}, other than the one to {@code other}. */
  private int[] carriers(int atom, int other, int[] variableForBond) {
    int[] found = new int[bondsAt[atom].length];
    int count = 0;
    for (int position = 0; position < bondsAt[atom].length; position++) {
      int bond = bondsAt[atom][position];
      if (neighbors[atom][position] != other && variableForBond[bond] >= 0) found[count++] = bond;
    }
    return Arrays.copyOf(found, count);
  }

  /** The carriers whose direction no configured double bond fixes. */
  private static int[] unfixed(int[] carriers, boolean[] fixedByDoubleBond, int[] variableForBond) {
    int[] kept = new int[carriers.length];
    int count = 0;
    for (int bond : carriers) if (!fixedByDoubleBond[variableForBond[bond]]) kept[count++] = bond;
    return Arrays.copyOf(kept, count);
  }

  /** 1 when the bond is written from {@code atom}, else -1. */
  private int outward(int bond, int atom) {
    return from[bond] == atom ? 1 : -1;
  }

  /** Whether a query atom has a plain bond mapped to a double bond, which gives its other bonds a direction. */
  private boolean atDoubleBond(int atom, int[] mapping, IAtomContainer molecule) {
    for (int position = 0; position < neighbors[atom].length; position++)
      if (expressions[bondsAt[atom][position]] == null
          && order(bond(molecule, mapping[atom], mapping[neighbors[atom][position]])) == 2) return true;
    return false;
  }

  private static IBond bond(IAtomContainer molecule, int first, int second) {
    return molecule.getBond(molecule.getAtom(first), molecule.getAtom(second));
  }

  /** The order of a bond, 4 when it is aromatic, as the native graph stores it. */
  private static int order(IBond bond) {
    if (bond.isAromatic()) return 4;
    IBond.Order order = bond.getOrder();
    return order == IBond.Order.SINGLE ? 1 : order == IBond.Order.DOUBLE ? 2 : order == IBond.Order.TRIPLE ? 3 : 0;
  }

  /** Whether a bond can carry / or \: Daylight reads directional bonds as single or aromatic. */
  private static boolean carrier(IBond bond) {
    int order = order(bond);
    return order == 1 || order == 4;
  }

  private static boolean primitive(int type, IBond bond) {
    switch (type) {
      case ANY:
        return true;
      case SINGLE:
        return order(bond) == 1;
      case DOUBLE:
        return order(bond) == 2;
      case TRIPLE:
        return order(bond) == 3;
      case AROMATIC:
        return bond.isAromatic();
      case RING:
        return bond.isInRing();
      default:
        return carrier(bond);
    }
  }

  private static boolean directionalPrimitive(int type) {
    return type == UP || type == DOWN || type == UP_OR_UNSPECIFIED || type == DOWN_OR_UNSPECIFIED;
  }

  private static boolean directional(Node expression) {
    if (expression.operands.length == 0) return directionalPrimitive(expression.type);
    for (Node operand : expression.operands) if (directional(operand)) return true;
    return false;
  }

  /**
   * The expression with the bond taking one direction, 1 up or -1 down; 0 reads every
   * directional primitive as false.
   */
  private static boolean withDirection(Node expression, IBond target, int direction) {
    switch (expression.type) {
      case NOT:
        return !withDirection(expression.operands[0], target, direction);
      case AND:
        for (Node operand : expression.operands)
          if (!withDirection(operand, target, direction)) return false;
        return true;
      case OR:
        for (Node operand : expression.operands)
          if (withDirection(operand, target, direction)) return true;
        return false;
      default:
        if (!directionalPrimitive(expression.type)) return primitive(expression.type, target);
        boolean up = expression.type == UP || expression.type == UP_OR_UNSPECIFIED;
        return carrier(target) && direction == (up ? 1 : -1);
    }
  }

  /**
   * Whether the expression accepts a bond that carries no direction, as one next to a double bond
   * without a configuration does: /? and \? do, / and \ leave it open, and other primitives read
   * the bond as usual. An expression left open does not accept it.
   */
  private static int withoutDirection(Node expression, IBond target) {
    switch (expression.type) {
      case NOT: {
        int result = withoutDirection(expression.operands[0], target);
        return result == UNKNOWN ? UNKNOWN : result == TRUE ? FALSE : TRUE;
      }
      case AND:
      case OR: {
        boolean unknown = false;
        for (Node operand : expression.operands) {
          int result = withoutDirection(operand, target);
          if (expression.type == AND && result == FALSE) return FALSE;
          if (expression.type == OR && result == TRUE) return TRUE;
          if (result == UNKNOWN) unknown = true;
        }
        if (unknown) return UNKNOWN;
        return expression.type == AND ? TRUE : FALSE;
      }
      default:
        if (!directionalPrimitive(expression.type))
          return primitive(expression.type, target) ? TRUE : FALSE;
        if (!carrier(target)) return FALSE;
        return expression.type == UP_OR_UNSPECIFIED || expression.type == DOWN_OR_UNSPECIFIED ? TRUE : UNKNOWN;
    }
  }

  private static boolean bondCharacter(char c) {
    return c == '-' || c == '=' || c == '#' || c == ':' || c == '~' || c == '@' || c == '!' || c == '/' || c == '\\';
  }

  /** A bond expression in the native precedence: negation, high conjunction, disjunction, low conjunction. */
  private static final class Parser {
    private final String text;
    private int at;

    Parser(String text, int at) {
      this.text = text;
      this.at = at;
    }

    private char peek() {
      return at < text.length() ? text.charAt(at) : '\0';
    }

    Node expression() {
      Node first = disjunction();
      if (peek() != ';') return first;
      List<Node> operands = new ArrayList<>();
      operands.add(first);
      while (peek() == ';') {
        at++;
        operands.add(disjunction());
      }
      return new Node(AND, operands.toArray(new Node[0]));
    }

    private Node disjunction() {
      Node first = conjunction();
      if (peek() != ',') return first;
      List<Node> operands = new ArrayList<>();
      operands.add(first);
      while (peek() == ',') {
        at++;
        operands.add(conjunction());
      }
      return new Node(OR, operands.toArray(new Node[0]));
    }

    private Node conjunction() {
      Node first = primitive();
      if (peek() != '&' && !bondCharacter(peek())) return first;
      List<Node> operands = new ArrayList<>();
      operands.add(first);
      while (peek() == '&' || bondCharacter(peek())) {
        if (peek() == '&') at++;
        operands.add(primitive());
      }
      return new Node(AND, operands.toArray(new Node[0]));
    }

    private Node primitive() {
      boolean negated = false;
      while (peek() == '!') {
        at++;
        negated = !negated;
      }
      int type;
      switch (peek()) {
        case '-': type = SINGLE; break;
        case '=': type = DOUBLE; break;
        case '#': type = TRIPLE; break;
        case ':': type = AROMATIC; break;
        case '~': type = ANY; break;
        case '@': type = RING; break;
        case '/': type = UP; break;
        case '\\': type = DOWN; break;
        default: throw new IllegalArgumentException("Invalid SMARTS: expected a bond primitive at position " + at);
      }
      at++;
      // /? and \? also accept a double bond without a configuration.
      if ((type == UP || type == DOWN) && peek() == '?') {
        at++;
        type = type == UP ? UP_OR_UNSPECIFIED : DOWN_OR_UNSPECIFIED;
      }
      Node primitive = new Node(type);
      return negated ? new Node(NOT, primitive) : primitive;
    }
  }
}
