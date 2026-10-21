/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 *
 * Common-bond (MCES) search: the conditions a mapping must meet beyond its common
 * bonds (required pairs, fragment limits, stereochemistry and whole rings), and the
 * selection of common bonds that meets them.
 *
 * Part of common_bond_search.hpp, which describes the search as a whole.
 */
#pragma once

#include "smsd/stereo.hpp"
#include "smsd/detail/common_bond_graphs.hpp"
#include "smsd/detail/union_find.hpp"

#include <algorithm>
#include <array>
#include <climits>
#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

namespace smsd {
namespace common_bond_detail {

/**
 * What a mapping of CommonBondSearch must satisfy beyond its common bonds: the pairs and bonds of
 * a required mapping, the size and number of its fragments, and, where the options ask for them,
 * mapped stereochemistry and whole rings, all as MCS::validateMapping() checks them. The search
 * maps atoms; a mapping here selects some of the common bonds of its atoms, and its score is the
 * number of bonds it selects, then the number of atoms. Without stereochemistry and rings a
 * mapping selects all its common bonds or none. Stereochemistry and rings are read from the
 * selected bonds alone, so a mapping may leave a common bond out to keep a centre whose ligands
 * disagree unresolved, or to keep only whole rings. The common bonds of the atoms bound every
 * selection of them, so the bounds of the search hold.
 *
 * The rings are met by the selections in which each selected bond on a ring of either molecule
 * lies on a wholly selected relevant ring of that molecule. Such selections are closed under
 * union, so each set of bonds holds a largest one, found by taking off bonds on no whole ring
 * until none is left (keepWholeRings()). A stereo centre is met once at most two of its ligands
 * are selected, or its ligands agree; a stereo double bond once it is left out, or one end keeps
 * no selected ligand, or its ligands agree. Taking bonds away never breaks either, and every
 * other condition holds of a selection once it holds of a smaller one. So each best selection is
 * the largest whole-ring part of the common bonds less a least set of bonds that resolves the
 * conflicts, and select() tries these sets depth first, one conflict at a time, each resolved in
 * every least way.
 */
class CommonBondConditions {
    /** A stereo centre or double bond whose selected ligands disagree. */
    struct Conflict {
        bool centre;
        int atom, other, bond;
        /** The selected bonds to the ligands of a centre, or those giving each end of a double bond a ligand. */
        std::vector<int> ligands, otherLigands;
    };
    /** The ligands joined to a centre or a double bond by selected bonds, as the stereo checks read them. */
    struct SelectedLigands {
        const std::vector<int>& mapping;
        const Side& query;
        const std::vector<uint8_t>& selection;
        int first, second;
        int operator[](int atom) const {
            if (atom == first || atom == second) return mapping[atom];
            int bond = query.bondBetween(first, atom);
            if (bond >= 0 && selection[bond]) return mapping[atom];
            if (second >= 0) {
                bond = query.bondBetween(second, atom);
                if (bond >= 0 && selection[bond]) return mapping[atom];
            }
            return -1;
        }
    };

    const CommonBondGraphs& graph_;
    const Side& query_;
    const Side& target_;
    const MolGraph& queryGraph_;
    const MolGraph& targetGraph_;
    Clock& clock_;
    const bool connected_;
    const int minFragmentAtoms_, maxFragments_;
    /** The target of each query atom that a required mapping fixes, -1 for none; empty without one. */
    std::vector<int> requiredTarget_;
    /** The query bonds a required mapping selects; empty without one. */
    std::vector<uint8_t> requiredBond_;
    const bool chirality_, bondStereo_, rings_;
    /** The bonds of each relevant ring of each molecule, by the bond numbers of its side. */
    std::vector<std::vector<int>> queryCycles_, targetCycles_;
    /** The rings through each bond and through each atom of each molecule. */
    std::vector<std::vector<int>> queryBondCycles_, targetBondCycles_, queryAtomCycles_, targetAtomCycles_;
    /** Whether each bond and atom of each molecule is on a ring, as the ring check reads it. */
    std::vector<uint8_t> queryRingBond_, targetRingBond_, queryRingAtom_, targetRingAtom_;
    /** The bonds of each ring left out of the selection being checked. */
    std::vector<int> queryMissing_, targetMissing_;
    /** The mapping being checked, and its number of atoms. */
    const std::vector<int>* mapping_ = nullptr;
    int atoms_ = 0;
    /** The target bond of each common bond of the mapping, and the selected query bond on each target bond, or -1. */
    std::vector<int> image_, preimage_;
    /** The common bonds, the first commonCount_, in query bond order. */
    std::vector<int> commonBonds_;
    int commonCount_ = 0;
    /** The selection under way, and the selection of a leaf once its rings are whole. */
    std::vector<uint8_t> selected_, kept_;
    /** Bonds still to check while rings are made whole. */
    std::vector<int> queue_;
    /** Union-find of the query atoms over the selected bonds, and the size of each set. */
    std::vector<int> parent_, setSize_;
    std::vector<Conflict> conflicts_;
    /** The fewest bonds a selection needs, the most selected so far, and the selected bonds before any conflict is resolved. */
    int leastBonds_ = 0, bestBonds_ = -1, upperBonds_ = 0;
    /** The most selections wanted; ties with the best are kept until there are as many. */
    size_t room_ = 0;
    /** The stamp of each query bond counted by fewestToLeaveOut() or persistentLoss(), and the stamp under way. */
    std::vector<int> bondMark_;
    int markStamp_ = 0;
    /** The query atoms with a tetrahedral tag where chirality is checked, and none otherwise. */
    std::vector<int> centres_;
    /** The best selections so far, as records. */
    std::vector<std::vector<int>> best_;

public:
    /**
     * The conditions of a search of `queryGraph` in `targetGraph`, whose fragments have at least
     * `minFragmentAtoms` atoms and are at most `maxFragments`, or a single one when `connected`,
     * that keep the target of each query atom in `requiredTarget`, -1 for none, unless it is
     * empty, and the query bonds in `requiredBonds`. The flags say which of chirality,
     * double-bond stereochemistry and whole rings are checked.
     */
    CommonBondConditions(const MolGraph& queryGraph, const MolGraph& targetGraph, const CommonBondGraphs& graph,
                         bool connected, int minFragmentAtoms, int maxFragments, const std::vector<int>& requiredTarget,
                         const std::vector<std::array<int, 2>>& requiredBonds, bool chirality, bool bondStereo, bool rings,
                         Clock& clock)
        : graph_(graph), query_(graph.query), target_(graph.target), queryGraph_(queryGraph), targetGraph_(targetGraph),
          clock_(clock), connected_(connected), minFragmentAtoms_(minFragmentAtoms), maxFragments_(maxFragments),
          requiredTarget_(requiredTarget), chirality_(chirality), bondStereo_(bondStereo), rings_(rings),
          image_(size_t(query_.bondCount), -1), preimage_(size_t(target_.bondCount), -1),
          commonBonds_(size_t(query_.bondCount), 0), selected_(size_t(query_.bondCount), 0),
          kept_(size_t(query_.bondCount), 0), parent_(size_t(query_.atomCount), 0), setSize_(size_t(query_.atomCount), 0),
          bondMark_(size_t(query_.bondCount), 0) {
        if (chirality_)
            for (int atom = 0; atom < static_cast<int>(queryGraph_.tetraChirality.size()); ++atom)
                if (queryGraph_.tetraChirality[atom] != 0) centres_.push_back(atom);
        if (!requiredTarget_.empty()) {
            requiredBond_.assign(size_t(query_.bondCount), 0);
            for (const auto& bond : requiredBonds) {
                const int queryBond = query_.bondBetween(bond[0], bond[1]);
                if (queryBond >= 0) requiredBond_[queryBond] = 1;
            }
        }
        if (rings_) {
            queryCycles_ = cycles(query_, queryGraph_);
            targetCycles_ = cycles(target_, targetGraph_);
            queryBondCycles_ = through(queryCycles_, query_.bondCount, nullptr);
            targetBondCycles_ = through(targetCycles_, target_.bondCount, nullptr);
            queryAtomCycles_ = through(queryCycles_, query_.atomCount, &query_);
            targetAtomCycles_ = through(targetCycles_, target_.atomCount, &target_);
            queryRingBond_ = ringBonds(query_, queryGraph_);
            targetRingBond_ = ringBonds(target_, targetGraph_);
            queryRingAtom_ = ringAtoms(queryGraph_);
            targetRingAtom_ = ringAtoms(targetGraph_);
            queryMissing_.assign(queryCycles_.size(), 0);
            targetMissing_.assign(targetCycles_.size(), 0);
        }
    }

    /** Whether some atom of `graph` has a tetrahedral tag. */
    static bool hasCentres(const MolGraph& graph) {
        for (int tag : graph.tetraChirality) if (tag != 0) return true;
        return false;
    }
    /** Whether some bond of `graph` has a double-bond configuration. */
    static bool hasStereoBonds(const MolGraph& graph) {
        if (!graph.hasDbStereo) return false;
        for (int atom = 0; atom < graph.n; ++atom)
            for (int other : graph.neighbors[atom])
                if (other > atom && graph.dbStereo(atom, other) != 0) return true;
        return false;
    }
    /** Whether some atom or bond of `graph` is on a ring. */
    static bool hasRings(const MolGraph& graph) {
        for (int atom = 0; atom < graph.n; ++atom) {
            if (graph.ring[atom]) return true;
            for (int other : graph.neighbors[atom]) if (graph.bondInRing(atom, other)) return true;
        }
        return !graph.computeRelevantCycles().empty();
    }

    /**
     * A least number of the common bonds of `mapping` that the best selection of every mapping
     * extending it leaves out for its tetrahedral centres, so that the search can bound the score
     * below it. A mapped centre with at least three ligands on common bonds that disagree with
     * its image disagrees however the mapping grows, as any three ligands fix the place of the
     * fourth, so every selection leaves out all but two of these bonds. Centres that share none
     * of them leave out different bonds.
     */
    int persistentLoss(const std::vector<int>& mapping) {
        if (centres_.empty()) return 0;
        clock_.tick(static_cast<int64_t>(centres_.size()));
        if (++markStamp_ == INT_MAX) {
            std::fill(bondMark_.begin(), bondMark_.end(), 0);
            markStamp_ = 1;
        }
        int total = 0;
        for (int atom : centres_) {
            if (mapping[atom] < 0) continue;
            int ligands = 0;
            bool shared = false;
            for (int ligand : queryGraph_.neighbors[atom]) {
                const int bond = commonBond(atom, ligand, mapping);
                if (bond < 0) continue;
                ++ligands;
                shared |= bondMark_[bond] == markStamp_;
            }
            if (ligands < 3 || shared || !centreDisagreesOnCommonBonds(atom, mapping)) continue;
            total += ligands - 2;
            for (int ligand : queryGraph_.neighbors[atom]) {
                const int bond = commonBond(atom, ligand, mapping);
                if (bond >= 0) bondMark_[bond] = markStamp_;
            }
        }
        return total;
    }

    /** Whether a mapping may select fewer than all its common bonds. */
    bool dropsBonds() const { return chirality_ || bondStereo_ || rings_; }
    /**
     * Whether exchanging twin terminal atoms takes each mapping onto one that meets the conditions
     * as well, for twins that carry no stereo tag that the conditions read and whose shared atom
     * carries none either, see untagged().
     */
    bool keepsTwins() const { return requiredTarget_.empty(); }
    /**
     * Whether `atom`, of the query when `onQuery` and otherwise of the target, carries no
     * tetrahedral tag that the chirality check reads and no stereo double bond that the
     * double-bond check reads. Only the tags of an atom and of its bonds tie the order of its
     * neighbours to the conditions, so exchanging two terminal neighbours of such an atom,
     * neither tagged, leaves every check as it was.
     */
    bool untagged(bool onQuery, int atom) const {
        const MolGraph& molecule = onQuery ? queryGraph_ : targetGraph_;
        if (chirality_ && !molecule.tetraChirality.empty() && molecule.tetraChirality[atom] != 0) return false;
        if (bondStereo_)
            for (int neighbor : molecule.neighbors[atom])
                if (molecule.dbStereo(atom, neighbor) != 0) return false;
        return true;
    }
    /** Whether adding a pendant pair to a mapping that meets the conditions gives one that meets them. */
    bool keepsPendantPairs() const { return !chirality_ && !bondStereo_; }
    /** The fewest atoms of a fragment, and the most fragments. */
    int minFragmentAtoms() const { return minFragmentAtoms_; }
    int maxFragments() const { return maxFragments_; }

    /**
     * The best selections of the common bonds of `mapping` that meet the conditions, if their
     * score is at least `least`: at most `room` of them are added to `out`, as records, and their
     * score is returned. A record is the mapping, followed, when it leaves common bonds out, by
     * their number and the query bonds left out in ascending order. Returns -1, adding nothing,
     * when no selection scores `least`.
     */
    int64_t select(const std::vector<int>& mapping, int64_t least, size_t room, std::vector<std::vector<int>>& out) {
        mapping_ = &mapping;
        atoms_ = 0;
        for (int queryAtom = 0; queryAtom < query_.atomCount; ++queryAtom) atoms_ += mapping[queryAtom] >= 0;
        if (!requiredTarget_.empty())
            for (int queryAtom = 0; queryAtom < query_.atomCount; ++queryAtom)
                if (requiredTarget_[queryAtom] >= 0 && mapping[queryAtom] != requiredTarget_[queryAtom]) return -1;
        commonCount_ = 0;
        clock_.tick(query_.bondCount);
        for (int queryBond = 0; queryBond < query_.bondCount; ++queryBond) {
            const int first = mapping[query_.bondFirst[queryBond]], second = mapping[query_.bondSecond[queryBond]];
            const int targetBond = first < 0 || second < 0 ? -1 : target_.bondBetween(first, second);
            const bool common = targetBond >= 0 && graph_.bondsMatch(queryBond, targetBond);
            image_[queryBond] = common ? targetBond : -1;
            selected_[queryBond] = common;
            if (common) commonBonds_[commonCount_++] = queryBond;
        }
        leastBonds_ = atoms_ >= atomsOf(least) ? bondsOf(least) : bondsOf(least) + 1;
        if (commonCount_ < leastBonds_) return -1;
        if (rings_) keepWholeRings(selected_);
        upperBonds_ = 0;
        for (int index = 0; index < commonCount_; ++index) upperBonds_ += selected_[commonBonds_[index]];
        if (upperBonds_ < leastBonds_) return -1;
        best_.clear();
        bestBonds_ = -1;
        room_ = room;
        if (chirality_ || bondStereo_) {
            findConflicts();
            resolve(0, 0);
        } else {
            leaf();
        }
        if (best_.empty()) return -1;
        for (size_t index = 0; index < best_.size() && index < room; ++index) out.push_back(std::move(best_[index]));
        return score(atoms_, bestBonds_);
    }

private:
    const std::vector<int>& mapping() const { return *mapping_; }

    /** The conflicts of the selection: the centres first, by query atom, then the double bonds. */
    void findConflicts() {
        conflicts_.clear();
        const auto& mapped = mapping();
        if (chirality_)
            for (int atom = 0; atom < query_.atomCount; ++atom)
                if (mapped[atom] >= 0 && centreDisagrees(atom, selected_)) {
                    std::array<int, 4> ligands{};
                    stereo_detail::tetrahedralLigands(queryGraph_, atom, ligands);
                    Conflict conflict{true, atom, -1, -1, {}, {}};
                    for (int ligand : ligands) {
                        if (ligand < 0) continue;
                        const int bond = query_.bondBetween(atom, ligand);
                        if (bond >= 0 && selected_[bond]) conflict.ligands.push_back(bond);
                    }
                    std::sort(conflict.ligands.begin(), conflict.ligands.end());
                    conflicts_.push_back(std::move(conflict));
                }
        if (bondStereo_)
            for (int index = 0; index < commonCount_; ++index) {
                const int bond = commonBonds_[index];
                const int first = query_.bondFirst[bond], second = query_.bondSecond[bond];
                if (!selected_[bond] || !doubleBondDisagrees(first, second, selected_)) continue;
                conflicts_.push_back(Conflict{false, first, second, bond, endLigands(first, second), endLigands(second, first)});
            }
        clock_.tick(query_.atomCount + static_cast<int64_t>(conflicts_.size()));
    }

    /** The selected bonds that give end `end` of the double bond to `other` a ligand: to its other neighbours, or from `other` to them. */
    std::vector<int> endLigands(int end, int other) const {
        std::vector<int> bonds;
        for (int neighbour : queryGraph_.neighbors[end]) {
            if (neighbour == other) continue;
            for (int from : {end, other}) {
                const int bond = query_.bondBetween(from, neighbour);
                if (bond >= 0 && selected_[bond] && std::find(bonds.begin(), bonds.end(), bond) == bonds.end())
                    bonds.push_back(bond);
            }
        }
        std::sort(bonds.begin(), bonds.end());
        return bonds;
    }

    /**
     * Resolves the conflicts from `from` on, with `dropped` bonds left out so far: the first
     * conflict not yet resolved is resolved in each least way, and each way is followed on. Once
     * none is left, the selection is a leaf.
     */
    void resolve(size_t from, int dropped) {
        clock_.tick();
        if (upperBonds_ - dropped - fewestToLeaveOut(from) < leastBonds_ || !fragmentsHold(selected_)) return;
        size_t open = from;
        while (open < conflicts_.size() && resolved(conflicts_[open])) ++open;
        if (open == conflicts_.size()) {
            leaf();
            return;
        }
        const Conflict& conflict = conflicts_[open];
        if (conflict.centre) {
            // A centre with at most two selected ligands is unresolved.
            int selectedLigands = 0;
            std::vector<int> choices;
            for (int bond : conflict.ligands) {
                if (!selected_[bond]) continue;
                ++selectedLigands;
                if (!required(bond)) choices.push_back(bond);
            }
            const int need = selectedLigands - 2;
            if (need <= 0 || need > static_cast<int>(choices.size())) return;
            leaveOut(choices, need, 0, open, dropped);
        } else {
            // The double bond itself, or every selected ligand of one end.
            const int bond = conflict.bond;
            if (selected_[bond] && !required(bond)) {
                selected_[bond] = 0;
                resolve(open, dropped + 1);
                selected_[bond] = 1;
            }
            dropAll(conflicts_[open].ligands, open, dropped);
            dropAll(conflicts_[open].otherLigands, open, dropped);
        }
    }

    /**
     * A least number of further bonds that the conflicts from `from` on leave out: a centre keeps
     * at most two selected ligands and a double bond loses at least one bond, and conflicts that
     * share no selected bond lose different bonds.
     */
    int fewestToLeaveOut(size_t from) {
        if (++markStamp_ == INT_MAX) {
            std::fill(bondMark_.begin(), bondMark_.end(), 0);
            markStamp_ = 1;
        }
        int total = 0;
        for (size_t index = from; index < conflicts_.size(); ++index) {
            const Conflict& conflict = conflicts_[index];
            int need;
            bool shared = false;
            if (conflict.centre) {
                int selectedLigands = 0;
                for (int bond : conflict.ligands)
                    if (selected_[bond]) {
                        ++selectedLigands;
                        shared |= bondMark_[bond] == markStamp_;
                    }
                need = selectedLigands - 2;
            } else {
                need = selected_[conflict.bond] && anySelected(conflict.ligands) && anySelected(conflict.otherLigands) ? 1 : 0;
                shared = need > 0 && (bondMark_[conflict.bond] == markStamp_ || anyMarked(conflict.ligands)
                                      || anyMarked(conflict.otherLigands));
            }
            if (need <= 0 || shared) continue;
            total += need;
            if (!conflict.centre) bondMark_[conflict.bond] = markStamp_;
            mark(conflict.ligands);
            if (!conflict.centre) mark(conflict.otherLigands);
        }
        return total;
    }
    /** Whether any of `bonds` is selected. */
    bool anySelected(const std::vector<int>& bonds) const {
        for (int bond : bonds) if (selected_[bond]) return true;
        return false;
    }
    /** Whether any selected bond of `bonds` is marked for the bound under way. */
    bool anyMarked(const std::vector<int>& bonds) const {
        for (int bond : bonds) if (selected_[bond] && bondMark_[bond] == markStamp_) return true;
        return false;
    }
    /** Marks the selected bonds of `bonds` for the bound under way. */
    void mark(const std::vector<int>& bonds) {
        for (int bond : bonds) if (selected_[bond]) bondMark_[bond] = markStamp_;
    }

    /** Leaves out each choice of `need` of `choices`, from `start` on. */
    void leaveOut(const std::vector<int>& choices, int need, size_t start, size_t open, int dropped) {
        if (need == 0) {
            resolve(open, dropped);
            return;
        }
        for (size_t index = start; index + need <= choices.size(); ++index) {
            selected_[choices[index]] = 0;
            leaveOut(choices, need - 1, index + 1, open, dropped + 1);
            selected_[choices[index]] = 1;
        }
    }

    /** Leaves out every selected bond of `bonds`, unless one is required or none is selected. */
    void dropAll(const std::vector<int>& bonds, size_t open, int dropped) {
        std::vector<int> dropping;
        for (int bond : bonds) {
            if (!selected_[bond]) continue;
            if (required(bond)) return;
            dropping.push_back(bond);
        }
        if (dropping.empty()) return;
        for (int bond : dropping) selected_[bond] = 0;
        resolve(open, dropped + static_cast<int>(dropping.size()));
        for (int bond : dropping) selected_[bond] = 1;
    }

    /** Whether a required mapping selects query bond `bond`. */
    bool required(int bond) const { return !requiredBond_.empty() && requiredBond_[bond]; }

    /** Whether `conflict` is resolved by the selection under way. */
    bool resolved(const Conflict& conflict) const {
        if (conflict.centre) return !centreDisagrees(conflict.atom, selected_);
        return !selected_[conflict.bond] || !doubleBondDisagrees(conflict.atom, conflict.other, selected_);
    }

    /**
     * Checks the selection under way once its conflicts are resolved: with its rings made whole,
     * it is kept if it meets every condition and is no worse than the best so far.
     */
    void leaf() {
        clock_.tick();
        kept_ = selected_;
        if (rings_) keepWholeRings(kept_);
        int count = 0;
        for (int index = 0; index < commonCount_; ++index) count += kept_[commonBonds_[index]];
        if (count < leastBonds_ || !valid(kept_)) return;
        std::vector<int> record = recordOf(kept_, count);
        if (count > bestBonds_) {
            best_.clear();
            bestBonds_ = count;
        } else {
            // A tie, while fewer than wanted are held, may repeat one held.
            clock_.tick(static_cast<int64_t>(best_.size()));
            for (const auto& other : best_) if (other == record) return;
        }
        best_.push_back(std::move(record));
        // While there is room for ties no smaller selection is of use; once full, only a larger one.
        leastBonds_ = best_.size() < room_ ? count : count + 1;
    }

    /** The record of the mapping with the selection `selection` of `count` common bonds. */
    std::vector<int> recordOf(const std::vector<uint8_t>& selection, int count) const {
        const int left = commonCount_ - count;
        std::vector<int> record(mapping().begin(), mapping().end());
        if (left > 0) {
            record.push_back(left);
            for (int index = 0; index < commonCount_; ++index)
                if (!selection[commonBonds_[index]]) record.push_back(commonBonds_[index]);
        }
        return record;
    }

    /**
     * Whether the selection meets every condition other than the bonds on rings, which it meets
     * already: the required bonds are selected, every mapped atom has a selected bond, the
     * fragments are as many and as large as asked, the mapped ring atoms lie on whole rings, and
     * the mapped stereochemistry agrees.
     */
    bool valid(const std::vector<uint8_t>& selection) {
        if (!fragmentsHold(selection)) return false;
        const auto& mapped = mapping();
        if (rings_)
            for (int atom = 0; atom < query_.atomCount; ++atom) {
                if (mapped[atom] < 0) continue;
                if ((queryRingAtom_[atom] && !anyWhole(queryAtomCycles_[atom], queryMissing_))
                    || (targetRingAtom_[mapped[atom]] && !anyWhole(targetAtomCycles_[mapped[atom]], targetMissing_)))
                    return false;
            }
        if (chirality_)
            for (int atom = 0; atom < query_.atomCount; ++atom)
                if (mapped[atom] >= 0 && centreDisagrees(atom, selection)) return false;
        if (bondStereo_)
            for (int index = 0; index < commonCount_; ++index) {
                const int bond = commonBonds_[index];
                if (selection[bond] && doubleBondDisagrees(query_.bondFirst[bond], query_.bondSecond[bond], selection))
                    return false;
            }
        return true;
    }

    /**
     * Whether the selection keeps every required bond and gives every mapped atom a selected bond,
     * and its fragments are as many and as large as asked. Once one of these fails, leaving more
     * bonds out does not mend it, so a selection that fails here leads to no valid one.
     */
    bool fragmentsHold(const std::vector<uint8_t>& selection) {
        const auto& mapped = mapping();
        clock_.tick(query_.atomCount + commonCount_);
        for (int atom = 0; atom < query_.atomCount; ++atom) {
            parent_[atom] = atom;
            setSize_[atom] = 0;
        }
        for (int index = 0; index < commonCount_; ++index) {
            const int bond = commonBonds_[index];
            if (!selection[bond]) {
                if (required(bond)) return false;
                continue;
            }
            const int first = query_.bondFirst[bond], second = query_.bondSecond[bond];
            setSize_[first] = setSize_[second] = 1;
            const int firstRoot = find(first), secondRoot = find(second);
            if (firstRoot != secondRoot) parent_[firstRoot] = secondRoot;
        }
        for (int atom = 0; atom < query_.atomCount; ++atom)
            if (mapped[atom] >= 0 && setSize_[atom] == 0) return false;
        int fragments = 0;
        std::fill(setSize_.begin(), setSize_.end(), 0);
        for (int atom = 0; atom < query_.atomCount; ++atom)
            if (mapped[atom] >= 0 && setSize_[find(atom)]++ == 0) ++fragments;
        if ((connected_ && fragments != 1) || fragments > maxFragments_) return false;
        for (int atom = 0; atom < query_.atomCount; ++atom)
            if (mapped[atom] >= 0 && find(atom) == atom && setSize_[atom] < minFragmentAtoms_) return false;
        return true;
    }

    /** The root of the set of `atom`. */
    int find(int atom) { return detail::unionFindRoot(parent_, atom); }

    /**
     * Takes off the selected bonds of `selection` that lie on a ring of either molecule but on no
     * wholly selected relevant ring of it, until none is left. What is left is the largest
     * selection within it whose ring bonds all lie on whole rings, as such selections are closed
     * under union. The counts of bonds missing from each ring are left as they are for it.
     */
    void keepWholeRings(std::vector<uint8_t>& selection) {
        clock_.tick(commonCount_ + static_cast<int64_t>(queryCycles_.size() + targetCycles_.size()));
        for (int index = 0; index < commonCount_; ++index) {
            const int bond = commonBonds_[index];
            preimage_[image_[bond]] = selection[bond] ? bond : -1;
        }
        for (size_t cycle = 0; cycle < queryCycles_.size(); ++cycle) {
            int missing = 0;
            for (int bond : queryCycles_[cycle]) missing += !selection[bond];
            queryMissing_[cycle] = missing;
        }
        for (size_t cycle = 0; cycle < targetCycles_.size(); ++cycle) {
            int missing = 0;
            for (int bond : targetCycles_[cycle]) missing += preimage_[bond] < 0;
            targetMissing_[cycle] = missing;
        }
        queue_.clear();
        for (int index = 0; index < commonCount_; ++index)
            if (selection[commonBonds_[index]]) queue_.push_back(commonBonds_[index]);
        while (!queue_.empty()) {
            const int bond = queue_.back();
            queue_.pop_back();
            clock_.tick();
            if (!selection[bond] || onWholeRings(bond)) continue;
            selection[bond] = 0;
            preimage_[image_[bond]] = -1;
            for (int cycle : queryBondCycles_[bond])
                if (queryMissing_[cycle]++ == 0)
                    for (int other : queryCycles_[cycle]) if (selection[other]) queue_.push_back(other);
            for (int cycle : targetBondCycles_[image_[bond]])
                if (targetMissing_[cycle]++ == 0)
                    for (int other : targetCycles_[cycle]) if (preimage_[other] >= 0) queue_.push_back(preimage_[other]);
        }
        for (int index = 0; index < commonCount_; ++index) preimage_[image_[commonBonds_[index]]] = -1;
    }

    /** Whether selected query bond `bond` and its image lie on whole rings wherever they lie on rings. */
    bool onWholeRings(int bond) const {
        return (!queryRingBond_[bond] || anyWhole(queryBondCycles_[bond], queryMissing_))
            && (!targetRingBond_[image_[bond]] || anyWhole(targetBondCycles_[image_[bond]], targetMissing_));
    }

    /** Whether any of `cycles` misses no bond. */
    static bool anyWhole(const std::vector<int>& cycles, const std::vector<int>& missing) {
        for (int cycle : cycles) if (missing[cycle] == 0) return true;
        return false;
    }

    /**
     * Whether the tetrahedral centre at mapped query atom `atom` and its image disagree on the
     * ligands joined to it by selected bonds, as the stereo check of the search reads them.
     */
    bool centreDisagrees(int atom, const std::vector<uint8_t>& selection) const {
        const auto& mapped = mapping();
        const int image = mapped[atom];
        const int queryTag = atom < static_cast<int>(queryGraph_.tetraChirality.size()) ? queryGraph_.tetraChirality[atom] : 0;
        const int targetTag = image < static_cast<int>(targetGraph_.tetraChirality.size()) ? targetGraph_.tetraChirality[image] : 0;
        if (queryTag == 0 || targetTag == 0) return false;
        std::array<int, 4> queryLigands{}, targetLigands{};
        if (!stereo_detail::tetrahedralLigands(queryGraph_, atom, queryLigands)
            || !stereo_detail::tetrahedralLigands(targetGraph_, image, targetLigands)) return false;
        return !stereo_detail::tetrahedralParityMatches(queryTag, queryLigands, targetTag, targetLigands,
                                                        SelectedLigands{mapped, query_, selection, atom, -1});
    }

    /** The query bond between `atom` and `ligand` if it is a common bond of `mapping`, and -1 otherwise. */
    int commonBond(int atom, int ligand, const std::vector<int>& mapping) const {
        if (mapping[ligand] < 0) return -1;
        const int bond = query_.bondBetween(atom, ligand);
        if (bond < 0) return -1;
        const int targetBond = target_.bondBetween(mapping[atom], mapping[ligand]);
        return targetBond >= 0 && graph_.bondsMatch(bond, targetBond) ? bond : -1;
    }

    /** The images of a centre and of its ligands on common bonds of a mapping, -1 for other atoms, as the stereo checks read them. */
    struct CommonLigands {
        const CommonBondConditions& conditions;
        const std::vector<int>& mapping;
        int centre;
        int operator[](int atom) const {
            if (atom == centre) return mapping[atom];
            return conditions.commonBond(centre, atom, mapping) >= 0 ? mapping[atom] : -1;
        }
    };

    /** Whether the centre at `atom` and its image under `mapping` disagree on the ligands joined to it by common bonds. */
    bool centreDisagreesOnCommonBonds(int atom, const std::vector<int>& mapping) const {
        const int image = mapping[atom];
        const int queryTag = atom < static_cast<int>(queryGraph_.tetraChirality.size()) ? queryGraph_.tetraChirality[atom] : 0;
        const int targetTag = image < static_cast<int>(targetGraph_.tetraChirality.size()) ? targetGraph_.tetraChirality[image] : 0;
        if (queryTag == 0 || targetTag == 0) return false;
        std::array<int, 4> queryLigands{}, targetLigands{};
        if (!stereo_detail::tetrahedralLigands(queryGraph_, atom, queryLigands)
            || !stereo_detail::tetrahedralLigands(targetGraph_, image, targetLigands)) return false;
        return !stereo_detail::tetrahedralParityMatches(queryTag, queryLigands, targetTag, targetLigands,
                                                        CommonLigands{*this, mapping, atom});
    }

    /**
     * Whether the double bond between mapped query atoms `first` and `second` and its image
     * disagree on the ligands joined to its ends by selected bonds, as the stereo check of the
     * search reads them.
     */
    bool doubleBondDisagrees(int first, int second, const std::vector<uint8_t>& selection) const {
        return !stereo_detail::doubleBondParityMatches(queryGraph_, first, second, targetGraph_,
                                                       SelectedLigands{mapping(), query_, selection, first, second});
    }

    /** The relevant rings of `graph` as bonds of `side`. */
    static std::vector<std::vector<int>> cycles(const Side& side, const MolGraph& graph) {
        std::vector<std::vector<int>> result;
        for (const auto& ring : graph.computeRelevantCycles()) {
            std::vector<int> bonds(ring.size());
            bool whole = true;
            for (size_t index = 0; index < ring.size() && whole; ++index) {
                bonds[index] = side.bondBetween(ring[index], ring[(index + 1) % ring.size()]);
                whole = bonds[index] >= 0;
            }
            if (whole) result.push_back(std::move(bonds));
        }
        return result;
    }

    /** The rings through each item, in ascending order: through each bond when `side` is null, otherwise through each atom of `side`. */
    static std::vector<std::vector<int>> through(const std::vector<std::vector<int>>& cycles, int items, const Side* side) {
        std::vector<std::vector<int>> result(static_cast<size_t>(items));
        for (size_t cycle = 0; cycle < cycles.size(); ++cycle)
            for (int bond : cycles[cycle]) {
                if (!side) {
                    result[bond].push_back(static_cast<int>(cycle));
                    continue;
                }
                // Each atom of a ring is the first or second atom of two of its bonds; list the ring once.
                for (int atom : {side->bondFirst[bond], side->bondSecond[bond]}) {
                    auto& list = result[atom];
                    if (list.empty() || list.back() != static_cast<int>(cycle)) list.push_back(static_cast<int>(cycle));
                }
            }
        return result;
    }

    /** Whether each bond of `side` is on a ring, as the ring check reads it. */
    static std::vector<uint8_t> ringBonds(const Side& side, const MolGraph& graph) {
        std::vector<uint8_t> result(static_cast<size_t>(side.bondCount), 0);
        for (int bond = 0; bond < side.bondCount; ++bond) result[bond] = graph.bondInRing(side.bondFirst[bond], side.bondSecond[bond]);
        return result;
    }

    /** Whether each atom of `graph` is on a ring, as the ring check reads it. */
    static std::vector<uint8_t> ringAtoms(const MolGraph& graph) {
        std::vector<uint8_t> result(static_cast<size_t>(graph.n), 0);
        for (int atom = 0; atom < graph.n; ++atom) result[atom] = graph.ring[atom];
        return result;
    }
};

} // namespace common_bond_detail
} // namespace smsd
