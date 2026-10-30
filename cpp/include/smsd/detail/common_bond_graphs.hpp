/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 *
 * Common-bond (MCES) search: each molecule as an int graph (Side), and both
 * molecules with the atom and bond pairs that may be mapped onto each other
 * (CommonBondGraphs).
 *
 * Part of common_bond_search.hpp, which describes the search as a whole.
 */
#pragma once

#include "smsd/mol_graph.hpp"
#include "smsd/detail/common_bond_classes.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

namespace smsd {
namespace common_bond_detail {

/**
 * One molecule of the search as an int graph, with the atoms numbered as in the molecule and
 * the bonds numbered from the lower atom up, in the order of its neighbours.
 */
struct Side {
    int atomCount = 0, bondCount = 0;
    /** The neighbours of an atom, and the bonds to them, are at places start[atom] to start[atom + 1] - 1. */
    std::vector<int> start, neighbour, neighbourBond;
    /** The lower and the higher atom of each bond. */
    std::vector<int> bondFirst, bondSecond;
    /** A hash of the element of each atom, for the search order and seeds. */
    std::vector<int> symbolHash;
    /** The kind and the class of each atom and of each bond; a class of -1 matches nothing. */
    std::vector<int> atomKind, bondKind, atomClass, bondClass;
    /** The label of each bond, or -1 if the other molecule lacks it. */
    std::vector<int> label;
    /** The fine label of each bond that has a label, or -1. */
    std::vector<int> fineLabel;

    /** The graph of a well-formed molecule. */
    explicit Side(const MolGraph& graph) : atomCount(graph.n), start(size_t(graph.n) + 1, 0) {
        for (int atom = 0; atom < atomCount; ++atom)
            start[atom + 1] = start[atom] + static_cast<int>(graph.neighbors[atom].size());
        neighbour.assign(start[atomCount], 0);
        neighbourBond.assign(start[atomCount], 0);
        bondCount = start[atomCount] / 2;
        bondFirst.assign(bondCount, 0);
        bondSecond.assign(bondCount, 0);
        int bonds = 0;
        for (int atom = 0; atom < atomCount; ++atom) {
            const auto& row = graph.neighbors[atom];
            for (size_t i = 0; i < row.size(); ++i) {
                neighbour[start[atom] + i] = row[i];
                if (atom < row[i]) {
                    bondFirst[bonds] = atom;
                    bondSecond[bonds] = row[i];
                    neighbourBond[start[atom] + i] = bonds++;
                }
            }
        }
        // A place from an atom to a lower neighbour takes the bond that the neighbour numbered.
        // The bonds are gathered by their higher atom, so that this takes time linear in the
        // bonds even next to a hub.
        std::vector<int> higherStart(size_t(atomCount) + 1, 0), byHigher(static_cast<size_t>(bondCount));
        std::vector<int> bondTo(size_t(atomCount), -1);
        for (int bond = 0; bond < bondCount; ++bond) ++higherStart[bondSecond[bond] + 1];
        for (int atom = 0; atom < atomCount; ++atom) higherStart[atom + 1] += higherStart[atom];
        std::vector<int> next(higherStart.begin(), higherStart.end() - 1);
        for (int bond = 0; bond < bondCount; ++bond) byHigher[next[bondSecond[bond]]++] = bond;
        for (int atom = 0; atom < atomCount; ++atom) {
            for (int i = higherStart[atom]; i < higherStart[atom + 1]; ++i) bondTo[bondFirst[byHigher[i]]] = byHigher[i];
            for (int place = start[atom]; place < start[atom + 1]; ++place)
                if (neighbour[place] < atom) neighbourBond[place] = bondTo[neighbour[place]];
        }
        symbolHash.assign(atomCount, 0);
        for (int atom = 0; atom < atomCount; ++atom) symbolHash[atom] = 1 + graph.atomicNum[atom];
        label.assign(bondCount, -1);
        fineLabel.assign(bondCount, -1);
    }

    /** Sets the kinds and classes of the atoms and bonds from those of their kinds. */
    void setKinds(const std::vector<int>& atomKinds, const std::vector<int>& atomKindClass,
                  const std::vector<int>& bondKinds, const std::vector<int>& bondKindClass) {
        atomKind = atomKinds;
        bondKind = bondKinds;
        atomClass.assign(atomCount, -1);
        bondClass.assign(bondCount, -1);
        for (int atom = 0; atom < atomCount; ++atom) atomClass[atom] = atomKindClass[atomKind[atom]];
        for (int bond = 0; bond < bondCount; ++bond) bondClass[bond] = bondKindClass[bondKind[bond]];
    }

    /**
     * Puts the label key of each bond into `keys` from `offset`: the classes of its ends, the
     * lower first, and its own class, or -1 if any of them is missing.
     */
    void labelKeys(int atomClasses, int bondClasses, std::vector<int64_t>& keys, int offset) const {
        for (int bond = 0; bond < bondCount; ++bond) {
            const int first = atomClass[bondFirst[bond]], second = atomClass[bondSecond[bond]];
            keys[offset + bond] = first < 0 || second < 0 || bondClass[bond] < 0 ? -1
                : (int64_t(std::min(first, second)) * atomClasses + std::max(first, second)) * bondClasses
                    + bondClass[bond];
        }
    }

    /**
     * Numbers the fine labels of the bonds that have a label: the kinds of the ends, the lower
     * first, and the kind of the bond. Returns their number.
     */
    int fineLabels(int atomKinds, int bondKinds) {
        std::vector<int64_t> keys(bondCount);
        for (int bond = 0; bond < bondCount; ++bond) {
            const int first = atomKind[bondFirst[bond]], second = atomKind[bondSecond[bond]];
            keys[bond] = label[bond] < 0 ? -1
                : (int64_t(std::min(first, second)) * atomKinds + std::max(first, second)) * bondKinds + bondKind[bond];
        }
        std::vector<int64_t> sorted;
        for (int64_t key : keys) if (key >= 0) sorted.push_back(key);
        std::sort(sorted.begin(), sorted.end());
        sorted.erase(std::unique(sorted.begin(), sorted.end()), sorted.end());
        for (int bond = 0; bond < bondCount; ++bond)
            fineLabel[bond] = keys[bond] < 0 ? -1
                : static_cast<int>(std::lower_bound(sorted.begin(), sorted.end(), keys[bond]) - sorted.begin());
        return static_cast<int>(sorted.size());
    }

    /**
     * Sets, for each fine label, its label in `group` and in `ends` the kinds of its ends, the
     * lower first, and of its bond.
     */
    void describeFineLabels(std::vector<int>& group, std::vector<std::array<int, 3>>& ends,
                            std::vector<uint8_t>& described) const {
        for (int bond = 0; bond < bondCount; ++bond) {
            const int fine = fineLabel[bond];
            if (fine < 0 || described[fine]) continue;
            const int first = atomKind[bondFirst[bond]], second = atomKind[bondSecond[bond]];
            group[fine] = label[bond];
            ends[fine] = {std::min(first, second), std::max(first, second), bondKind[bond]};
            described[fine] = 1;
        }
    }

    /** The number of items of each fine class, from the fine class of each item, -1 for none. */
    static std::vector<int> counts(const std::vector<int>& fine, int classes) {
        std::vector<int> result(classes, 0);
        for (int item : fine)
            if (item >= 0) ++result[item];
        return result;
    }

    /** Puts the compatibility key of each bond into `keys` from `offset`. */
    template <class BondKey>
    void bondKeys(const MolGraph& graph, BondKey& bondKey, std::vector<int64_t>& keys, int offset) const {
        for (int bond = 0; bond < bondCount; ++bond)
            keys[offset + bond] = bondKey(graph, bondFirst[bond], bondSecond[bond]);
    }

    /** The number of bonds at an atom. */
    int degree(int atom) const { return start[atom + 1] - start[atom]; }

    /** The most bonds at any atom. */
    int maxDegree() const {
        int most = 0;
        for (int atom = 0; atom < atomCount; ++atom) most = std::max(most, degree(atom));
        return most;
    }

    /** The bond between two atoms, or -1 if there is none. The shorter list of neighbours is read. */
    int bondBetween(int atom, int other) const {
        if (degree(atom) > degree(other)) std::swap(atom, other);
        for (int i = start[atom]; i < start[atom + 1]; ++i)
            if (neighbour[i] == other) return neighbourBond[i];
        return -1;
    }

    /** The dense ranks of the atoms by their number of choices, then their degree, then their element. */
    std::vector<int> invariants(const std::vector<int>& choices) const {
        std::vector<int64_t> keys(atomCount);
        for (int atom = 0; atom < atomCount; ++atom) keys[atom] = symbolHash[atom];
        const std::vector<int> symbol = denseRanks(keys);
        const int64_t degrees = int64_t(maxDegree()) + 1;
        for (int atom = 0; atom < atomCount; ++atom)
            keys[atom] = (int64_t(choices[atom]) * degrees + degree(atom)) * atomCount + symbol[atom];
        return denseRanks(keys);
    }

    /**
     * Atom signatures by radius, radius * atomCount + atom: element and degree at 0, then each
     * mixed with the sum over its neighbours.
     */
    std::vector<uint64_t> signatures(Clock& clock) const {
        std::vector<uint64_t> signature(size_t(RADII) * size_t(atomCount));
        for (int atom = 0; atom < atomCount; ++atom)
            signature[atom] = mix(uint64_t(int64_t(symbolHash[atom])) * 31U + uint64_t(degree(atom)));
        for (int radius = 1; radius < RADII; ++radius) {
            const size_t previous = size_t(radius - 1) * atomCount, current = size_t(radius) * atomCount;
            for (int atom = 0; atom < atomCount; ++atom) {
                clock.tick(1 + degree(atom));
                uint64_t sum = 0;
                for (int i = start[atom]; i < start[atom + 1]; ++i) sum += mix(signature[previous + neighbour[i]]);
                signature[current + atom] = mix(signature[previous + atom] * GOLDEN_GAMMA + sum);
            }
        }
        return signature;
    }

    /** The most atoms and bonds of a component once the bonds whose label the other molecule lacks are removed. */
    int64_t componentBound() const {
        std::vector<uint8_t> seen(atomCount, 0);
        std::vector<int> queue(atomCount);
        int mostAtoms = 0, mostBonds = 0;
        for (int root = 0; root < atomCount; ++root) {
            if (seen[root] || atomClass[root] < 0) continue;
            int size = 1, links = 0;
            queue[0] = root;
            seen[root] = 1;
            for (int head = 0; head < size; ++head) {
                const int atom = queue[head];
                for (int place = start[atom]; place < start[atom + 1]; ++place) {
                    const int other = neighbour[place];
                    if (label[neighbourBond[place]] < 0) continue;
                    ++links;
                    if (!seen[other]) {
                        seen[other] = 1;
                        queue[size++] = other;
                    }
                }
            }
            mostAtoms = std::max(mostAtoms, size);
            // Each bond is met from both ends.
            mostBonds = std::max(mostBonds, links / 2);
        }
        return score(mostAtoms, mostBonds);
    }
};

/**
 * The query and target of a CommonBondSearch as int graphs, with the atom and bond pairs that
 * may be mapped onto each other, their kinds and classes, the bond labels shared by both
 * molecules and the groups for the b-matching bounds. The molecules are only read.
 */
class CommonBondGraphs {
public:
    /** The graphs of the query and of the target. */
    Side query, target;
    /** The number of atom classes and of bond classes, each shared by both molecules. */
    int atomClasses = 0, bondClasses = 0;
    /**
     * The number of bond labels shared by both molecules. The label of a bond is the classes of
     * its ends and its own class, so a common bond and its image share a label.
     */
    int labels = 0;
    /** The lower and the higher atom class at the ends of each label. */
    std::vector<int> labelLow, labelHigh;
    /** The atom kinds of both molecules, grouped by atom class. */
    Groups atomGroups;
    /**
     * The fine labels of both molecules, the kinds of the ends and of the bond, grouped by
     * label. A query and a target fine label match when their bond kinds do and their ends do
     * either way round.
     */
    Groups labelGroups;

    /**
     * Prepares the search of `queryGraph` in `targetGraph`. `atomMatches(queryAtom, targetAtom)`
     * decides the atom pairs and `bondMatches(queryFirst, querySecond, targetFirst, targetSecond)`
     * the bond pairs; `bondKey(graph, first, second)` gives equal keys to bonds that the bond
     * predicate cannot tell apart. The target atoms in `excludedTargetAtoms` match nothing. Use
     * supports() first.
     */
    template <class AtomMatch, class BondMatch, class BondKey>
    CommonBondGraphs(const MolGraph& queryGraph, const MolGraph& targetGraph, AtomMatch& atomMatches,
                     BondMatch& bondMatches, BondKey& bondKey, const std::vector<int>& excludedTargetAtoms, Clock& clock)
        : query(queryGraph), target(targetGraph) {
        bondWords_ = words(target.bondCount);
        atomWords_ = words(target.atomCount);
        std::vector<uint8_t> excluded(target.atomCount, 0);
        for (int atom : excludedTargetAtoms)
            if (atom >= 0 && atom < target.atomCount) excluded[atom] = 1;
        atomPairs_.assign(size_t(query.atomCount) * size_t(atomWords_), 0);
        for (int queryAtom = 0; queryAtom < query.atomCount; ++queryAtom) {
            clock.tick(1 + target.atomCount);
            for (int targetAtom = 0; targetAtom < target.atomCount; ++targetAtom)
                if (!excluded[targetAtom] && atomMatches(queryAtom, targetAtom))
                    setBit(atomPairs_, atomWords_, queryAtom, targetAtom);
        }
        computeBondPairs(queryGraph, targetGraph, bondMatches, bondKey, clock);
        const Kinds atomKinds(atomPairs_, query.atomCount, target.atomCount, clock);
        const Kinds bondKinds(bondPairs_, query.bondCount, target.bondCount, clock);
        atomClasses = atomKinds.classes;
        bondClasses = bondKinds.classes;
        atomGroups = atomKinds.groups();
        query.setKinds(atomKinds.queryKind, atomKinds.queryClass, bondKinds.queryKind, bondKinds.queryClass);
        target.setKinds(atomKinds.targetKind, atomKinds.targetClass, bondKinds.targetKind, bondKinds.targetClass);
        // The labels, of classes and of kinds, of every bond whose ends and itself have classes.
        std::vector<int64_t> classKeys(size_t(query.bondCount) + target.bondCount);
        query.labelKeys(atomClasses, bondClasses, classKeys, 0);
        target.labelKeys(atomClasses, bondClasses, classKeys, query.bondCount);
        labels = shared(classKeys, query.bondCount);
        labelLow.assign(labels, 0);
        labelHigh.assign(labels, 0);
        for (int bond = 0; bond < static_cast<int>(classKeys.size()); ++bond) {
            Side& side = bond < query.bondCount ? query : target;
            const int sideBond = bond < query.bondCount ? bond : bond - query.bondCount;
            const int bondLabel = static_cast<int>(classKeys[bond]);
            side.label[sideBond] = bondLabel;
            if (bondLabel >= 0) {
                const int first = side.atomClass[side.bondFirst[sideBond]];
                const int second = side.atomClass[side.bondSecond[sideBond]];
                labelLow[bondLabel] = std::min(first, second);
                labelHigh[bondLabel] = std::max(first, second);
            }
        }
        const int queryFine = query.fineLabels(atomKinds.queryKinds, bondKinds.queryKinds);
        const int targetFine = target.fineLabels(atomKinds.targetKinds, bondKinds.targetKinds);
        labelGroups = fineLabelGroups(queryFine, targetFine, atomKinds, bondKinds);
    }

    /**
     * Whether the search can take the two graphs: neither has a bond from an atom to itself or
     * a second bond between the same atoms, every bond is listed at both its atoms, the atom
     * pairs fit an int and the bond pairs fit MAX_BOND_PAIR_WORDS.
     */
    static bool supports(const MolGraph& queryGraph, const MolGraph& targetGraph) {
        if (int64_t(queryGraph.n) * targetGraph.n > MAX_PAIRS) return false;
        if (!wellFormed(queryGraph) || !wellFormed(targetGraph)) return false;
        const int64_t queryBonds = bondCount(queryGraph), targetWords = (bondCount(targetGraph) + 63) >> 6;
        // Divided rather than multiplied, so that no product of bond counts can overflow.
        return targetWords == 0 || queryBonds <= MAX_BOND_PAIR_WORDS / targetWords;
    }

    /** The compatible atom pairs, a row of target atoms for each query atom. */
    const std::vector<uint64_t>& compatiblePairs() const { return atomPairs_; }

    /** Whether a query bond may be mapped onto a target bond. */
    bool bondsMatch(int queryBond, int targetBond) const {
        return hasBit(bondPairs_, bondWords_, queryBond, targetBond);
    }

    /**
     * The common bonds of a mapping, each as {query atom, query atom, target atom, target atom},
     * the lower query atom first, in the order of the query bonds. The mapping may be a record of
     * the search, which then leaves some common bonds out, see CommonBondConditions::select().
     */
    std::vector<std::array<int, 4>> commonBonds(const std::vector<int>& mapping) const {
        std::vector<std::array<int, 4>> bonds;
        // A record of the search may leave common bonds out, see CommonBondConditions::select().
        size_t skipped = size_t(query.atomCount) + 1;
        const size_t skippedEnd = mapping.size() > size_t(query.atomCount) ? skipped + size_t(mapping[query.atomCount]) : skipped;
        for (int queryBond = 0; queryBond < query.bondCount; ++queryBond) {
            const int first = mapping[query.bondFirst[queryBond]], second = mapping[query.bondSecond[queryBond]];
            if (first < 0 || second < 0) continue;
            if (skipped < skippedEnd && mapping[skipped] == queryBond) {
                ++skipped;
                continue;
            }
            const int targetBond = target.bondBetween(first, second);
            if (targetBond >= 0 && bondsMatch(queryBond, targetBond))
                bonds.push_back({query.bondFirst[queryBond], query.bondSecond[queryBond], first, second});
        }
        return bonds;
    }

    /**
     * The label bound, as a score. No mapping has more atoms than the most atoms of matching
     * kinds that can be paired off within each atom class, nor more common bonds than the most
     * bonds of matching fine labels that can be paired off within each label. A connected
     * mapping is also capped by the largest component of either molecule without the bonds
     * whose label the other lacks.
     */
    int64_t labelBound(bool connected) {
        int64_t mostAtoms = atomGroups.bound(Side::counts(query.atomKind, static_cast<int>(atomGroups.queryGroup.size())),
                                             Side::counts(target.atomKind, static_cast<int>(atomGroups.targetGroup.size())));
        int64_t mostBonds = labelGroups.bound(Side::counts(query.fineLabel, static_cast<int>(labelGroups.queryGroup.size())),
                                              Side::counts(target.fineLabel, static_cast<int>(labelGroups.targetGroup.size())));
        if (connected) {
            const int64_t queryPart = query.componentBound(), targetPart = target.componentBound();
            mostAtoms = std::min<int64_t>(mostAtoms, std::min(atomsOf(queryPart), atomsOf(targetPart)));
            mostBonds = std::min<int64_t>(mostBonds, std::min(bondsOf(queryPart), bondsOf(targetPart)));
        }
        return score(static_cast<int>(mostAtoms), static_cast<int>(mostBonds));
    }

private:
    int bondWords_ = 0, atomWords_ = 0;
    std::vector<uint64_t> atomPairs_, bondPairs_;

    /**
     * The compatible bond pairs. The bond predicate reads only what the bond key holds, so it
     * is asked once for each pair of bond kinds that the molecules have, on a bond of each kind.
     * A pair must match either way round, as the rows do not say which end maps onto which.
     */
    template <class BondMatch, class BondKey>
    void computeBondPairs(const MolGraph& queryGraph, const MolGraph& targetGraph, BondMatch& bondMatches,
                          BondKey& bondKey, Clock& clock) {
        std::vector<int64_t> keys(size_t(query.bondCount) + target.bondCount);
        query.bondKeys(queryGraph, bondKey, keys, 0);
        target.bondKeys(targetGraph, bondKey, keys, query.bondCount);
        const std::vector<int> kind = denseRanks(keys);
        int kinds = 0;
        for (int bondKind : kind) kinds = std::max(kinds, bondKind + 1);
        std::vector<int> queryExample(kinds, -1), targetExample(kinds, -1);
        for (int bond = query.bondCount - 1; bond >= 0; --bond) queryExample[kind[bond]] = bond;
        for (int bond = target.bondCount - 1; bond >= 0; --bond) targetExample[kind[query.bondCount + bond]] = bond;
        std::vector<uint64_t> byKind(size_t(kinds) * size_t(bondWords_), 0);
        for (int targetBond = 0; targetBond < target.bondCount; ++targetBond)
            setBit(byKind, bondWords_, kind[query.bondCount + targetBond], targetBond);
        std::vector<uint8_t> compatible(size_t(kinds) * size_t(kinds), 0);
        for (int queryKind = 0; queryKind < kinds; ++queryKind) {
            const int queryBond = queryExample[queryKind];
            for (int targetKind = 0; queryBond >= 0 && targetKind < kinds; ++targetKind) {
                const int targetBond = targetExample[targetKind];
                clock.tick();
                if (targetBond < 0) continue;
                const int queryFirst = query.bondFirst[queryBond], querySecond = query.bondSecond[queryBond];
                const int targetFirst = target.bondFirst[targetBond], targetSecond = target.bondSecond[targetBond];
                compatible[size_t(queryKind) * kinds + targetKind] =
                    bondMatches(queryFirst, querySecond, targetFirst, targetSecond)
                    && bondMatches(queryFirst, querySecond, targetSecond, targetFirst) ? 1 : 0;
            }
        }
        bondPairs_.assign(size_t(query.bondCount) * size_t(bondWords_), 0);
        for (int queryBond = 0; queryBond < query.bondCount; ++queryBond) {
            const int queryKind = kind[queryBond];
            for (int targetKind = 0; targetKind < kinds; ++targetKind) {
                clock.tick(1 + bondWords_);
                if (!compatible[size_t(queryKind) * kinds + targetKind]) continue;
                for (int word = 0; word < bondWords_; ++word)
                    bondPairs_[size_t(queryBond) * bondWords_ + word] |= byKind[size_t(targetKind) * bondWords_ + word];
            }
        }
    }

    /**
     * Replaces each key that is not negative by the place among the distinct keys of both
     * molecules that both have, the query keys first, or by -1 if only one has it. Returns the
     * number of shared keys.
     */
    static int shared(std::vector<int64_t>& keys, int queryItems) {
        const std::vector<int> rank = denseRanks(keys);
        int ranks = 0;
        for (int itemRank : rank) ranks = std::max(ranks, itemRank + 1);
        std::vector<uint8_t> inQuery(ranks, 0), inTarget(ranks, 0);
        for (int item = 0; item < static_cast<int>(keys.size()); ++item)
            if (keys[item] >= 0) (item < queryItems ? inQuery : inTarget)[rank[item]] = 1;
        std::vector<int> place(ranks, -1);
        int count = 0;
        for (int keyRank = 0; keyRank < ranks; ++keyRank) place[keyRank] = inQuery[keyRank] && inTarget[keyRank] ? count++ : -1;
        for (int item = 0; item < static_cast<int>(keys.size()); ++item) keys[item] = keys[item] < 0 ? -1 : place[rank[item]];
        return count;
    }

    /** The fine labels of both molecules grouped by label, with the pairs of them that match. */
    Groups fineLabelGroups(int queryFine, int targetFine, const Kinds& atomKinds, const Kinds& bondKinds) const {
        std::vector<int> queryGroup(queryFine, -1), targetGroup(targetFine, -1);
        std::vector<std::array<int, 3>> queryEnds(queryFine), targetEnds(targetFine);
        std::vector<uint8_t> queryDescribed(queryFine, 0), targetDescribed(targetFine, 0);
        query.describeFineLabels(queryGroup, queryEnds, queryDescribed);
        target.describeFineLabels(targetGroup, targetEnds, targetDescribed);
        return Groups(labels, std::move(queryGroup), std::move(targetGroup),
            [&](int queryLabel, int targetLabel) {
                const auto& q = queryEnds[queryLabel];
                const auto& t = targetEnds[targetLabel];
                return bondKinds.matches(q[2], t[2])
                    && ((atomKinds.matches(q[0], t[0]) && atomKinds.matches(q[1], t[1]))
                        || (atomKinds.matches(q[0], t[1]) && atomKinds.matches(q[1], t[0])));
            });
    }

    /**
     * Whether every neighbour list is free of the atom itself and of repeats, and every bond is
     * listed at both its atoms. The atoms that list each atom are gathered first, so that the
     * check takes time linear in the bonds even next to a hub.
     */
    static bool wellFormed(const MolGraph& graph) {
        if (graph.n < 0 || graph.neighbors.size() != size_t(graph.n) || graph.atomicNum.size() < size_t(graph.n))
            return false;
        std::vector<int> seen(graph.n, 0), listedStart(size_t(graph.n) + 1, 0);
        for (int atom = 0; atom < graph.n; ++atom)
            for (int other : graph.neighbors[atom]) {
                if (other < 0 || other >= graph.n || other == atom || seen[other] == atom + 1) return false;
                if (graph.bondOrder(atom, other) != graph.bondOrder(other, atom)
                    || graph.bondAromatic(atom, other) != graph.bondAromatic(other, atom)
                    || graph.bondInRing(atom, other) != graph.bondInRing(other, atom)) return false;
                seen[other] = atom + 1;
                ++listedStart[other + 1];
            }
        // Each atom is listed by as many atoms as it lists, and then by exactly those it lists.
        for (int atom = 0; atom < graph.n; ++atom) {
            if (listedStart[atom + 1] != static_cast<int>(graph.neighbors[atom].size())) return false;
            listedStart[atom + 1] += listedStart[atom];
        }
        std::vector<int> listedBy(static_cast<size_t>(listedStart[graph.n])), next(listedStart.begin(), listedStart.end() - 1);
        for (int atom = 0; atom < graph.n; ++atom)
            for (int other : graph.neighbors[atom]) listedBy[next[other]++] = atom;
        std::fill(seen.begin(), seen.end(), 0);
        for (int atom = 0; atom < graph.n; ++atom) {
            for (int other : graph.neighbors[atom]) seen[other] = atom + 1;
            for (int place = listedStart[atom]; place < listedStart[atom + 1]; ++place)
                if (seen[listedBy[place]] != atom + 1) return false;
        }
        return true;
    }

    /** The number of bonds of a well-formed graph. */
    static int64_t bondCount(const MolGraph& graph) {
        int64_t ends = 0;
        for (const auto& row : graph.neighbors) ends += static_cast<int64_t>(row.size());
        return ends / 2;
    }
};

} // namespace common_bond_detail
} // namespace smsd
