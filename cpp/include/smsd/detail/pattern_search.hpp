/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 *
 * SMARTS pattern search: the pattern plan, recursive SMARTS evaluation with its
 * per-target memo, and the pattern search itself.
 *
 * Part of search.hpp.
 */
#pragma once

#include "smsd/smarts_parser.hpp"
#include "smsd/detail/search_results.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <limits>
#include <memory>
#include <numeric>
#include <optional>
#include <set>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace smsd {
namespace search_detail {

inline thread_local Budget* recursiveBudget = nullptr;
inline bool recursiveMatch(const SmartsQuery& pattern, const MolGraph& target, int anchor);

class RecursiveScope {
    Budget* previousBudget_ = recursiveBudget;
    smarts_detail::RecursiveMatcher previousMatcher_ = smarts_detail::recursiveMatcher;
public:
    explicit RecursiveScope(Budget& budget) {
        recursiveBudget = &budget;
        smarts_detail::recursiveMatcher = recursiveMatch;
    }
    ~RecursiveScope() {
        recursiveBudget = previousBudget_;
        smarts_detail::recursiveMatcher = previousMatcher_;
    }
};

/** Bond lookup for a SMARTS pattern; a repeated atom pair keeps its first bond. */
class PatternBonds {
    static constexpr int DENSE_ATOMS = 128;
    int atoms_ = 0;
    bool dense_ = false;
    std::vector<int> table_;
    std::unordered_map<int64_t, int> sparse_;
public:
    /** The bytes of the table a pattern of so many atoms would have, or zero when it would have none. */
    static size_t tableBytes(int atoms) {
        return atoms <= DENSE_ATOMS ? static_cast<size_t>(atoms) * static_cast<size_t>(atoms) * sizeof(int) : 0;
    }
    /** Small patterns have a table of atom pairs unless `table` is false; others, and those, use a map. */
    explicit PatternBonds(const SmartsQuery& pattern, bool table = true)
        : atoms_(pattern.atomCount()), dense_(table && pattern.atomCount() <= DENSE_ATOMS) {
        if (dense_) table_.assign(static_cast<size_t>(atoms_) * atoms_, -1);
        for (size_t bond = 0; bond < pattern.bonds.size(); ++bond) {
            const int first = pattern.bonds[bond].from, second = pattern.bonds[bond].to;
            if (!dense_) {
                sparse_.emplace(MolGraph::bondKey(first, second), static_cast<int>(bond));
                continue;
            }
            int& entry = table_[static_cast<size_t>(first) * atoms_ + second];
            if (entry < 0) entry = table_[static_cast<size_t>(second) * atoms_ + first] = static_cast<int>(bond);
        }
    }
    int find(int first, int second) const {
        if (dense_) return table_[static_cast<size_t>(first) * atoms_ + second];
        const auto bond = sparse_.find(MolGraph::bondKey(first, second));
        return bond == sparse_.end() ? -1 : bond->second;
    }
};

/** Chirality is resolved after a complete mapping, as in QueryView::atomMatches(). */
inline bool patternAtomMatches(const SmartsQuery& pattern, int atom, const MolGraph& target, int targetAtom) {
    return smarts_detail::evalAtomExprWithStereo(pattern.atoms[atom].expr, target, targetAtom,
                                                 [](const smarts_detail::AtomExpr&) {
        return smarts_detail::StereoTruth::UNKNOWN;
    }) != smarts_detail::StereoTruth::FALSE_VALUE;
}

/**
 * Structural equality; equal expressions select the same target atoms. Atom
 * classes always match, so their numbers are not compared.
 */
inline bool sameAtomExpr(const smarts_detail::AtomExpr& first, const smarts_detail::AtomExpr& second) {
    using smarts_detail::AtomPrimType;
    if (first.op != second.op || first.children.size() != second.children.size()) return false;
    if (first.op == smarts_detail::ExprOp::PRIM) {
        if (first.primType != second.primType) return false;
        if (first.primType != AtomPrimType::ATOM_CLASS
            && (first.intVal != second.intVal || first.hasValue != second.hasValue
                || first.rangeQuery != second.rangeQuery || first.rangeMin != second.rangeMin
                || first.rangeMax != second.rangeMax || first.aromaticity != second.aromaticity
                || first.orUnspecified != second.orUnspecified
                || first.recursiveSmarts != second.recursiveSmarts)) return false;
    }
    for (size_t child = 0; child < first.children.size(); ++child)
        if (!sameAtomExpr(first.children[child], second.children[child])) return false;
    return true;
}

/** Consistent with sameAtomExpr(). */
inline uint64_t atomExprHash(const smarts_detail::AtomExpr& expression) {
    auto mix = [](uint64_t hash, uint64_t value) {
        hash ^= value + 0x9E3779B97F4A7C15ull + (hash << 6) + (hash >> 2);
        return hash;
    };
    uint64_t hash = mix(static_cast<uint64_t>(expression.op), expression.children.size());
    if (expression.op == smarts_detail::ExprOp::PRIM) {
        hash = mix(hash, static_cast<uint64_t>(expression.primType));
        if (expression.primType != smarts_detail::AtomPrimType::ATOM_CLASS) {
            hash = mix(hash, static_cast<uint32_t>(expression.intVal));
            hash = mix(hash, (expression.hasValue ? 1u : 0u) | (expression.rangeQuery ? 2u : 0u)
                                 | (expression.orUnspecified ? 4u : 0u));
            hash = mix(hash, static_cast<uint32_t>(expression.rangeMin));
            hash = mix(hash, static_cast<uint32_t>(expression.rangeMax));
            hash = mix(hash, static_cast<uint32_t>(expression.aromaticity + 1));
            if (!expression.recursiveSmarts.empty()) hash = mix(hash, std::hash<std::string>{}(expression.recursiveSmarts));
        }
    }
    for (const auto& child : expression.children) hash = mix(hash, atomExprHash(child));
    return hash;
}

/** Element (-1 when open) and aromaticity (-1 when open) shared by every atom an expression accepts. */
struct AtomRequirement {
    int element = -1;
    int aromatic = -1;
};

inline AtomRequirement atomRequirement(const smarts_detail::AtomExpr& expression) {
    using smarts_detail::AtomPrimType;
    using smarts_detail::ExprOp;
    constexpr int ELEMENTS = 256;
    AtomRequirement result;
    switch (expression.op) {
        case ExprOp::PRIM:
            if ((expression.primType == AtomPrimType::ELEMENT || expression.primType == AtomPrimType::ATOMIC_NUM)
                && expression.intVal >= 0 && expression.intVal < ELEMENTS) {
                result.element = expression.intVal;
                if (expression.primType == AtomPrimType::ELEMENT) result.aromatic = expression.aromaticity;
            }
            else if (expression.primType == AtomPrimType::AROMATIC) result.aromatic = 1;
            else if (expression.primType == AtomPrimType::ALIPHATIC) result.aromatic = 0;
            return result;
        case ExprOp::AND: {
            // A contradiction accepts nothing; leaving it open is still a valid bound.
            bool elementConflict = false, aromaticConflict = false;
            for (const auto& child : expression.children) {
                const auto part = atomRequirement(child);
                if (part.element >= 0) {
                    if (result.element >= 0 && result.element != part.element) elementConflict = true;
                    result.element = part.element;
                }
                if (part.aromatic >= 0) {
                    if (result.aromatic >= 0 && result.aromatic != part.aromatic) aromaticConflict = true;
                    result.aromatic = part.aromatic;
                }
            }
            if (elementConflict) result.element = -1;
            if (aromaticConflict) result.aromatic = -1;
            return result;
        }
        case ExprOp::OR:
            for (size_t child = 0; child < expression.children.size(); ++child) {
                const auto part = atomRequirement(expression.children[child]);
                if (child == 0) { result = part; continue; }
                if (result.element != part.element) result.element = -1;
                if (result.aromatic != part.aromatic) result.aromatic = -1;
            }
            return result;
        case ExprOp::NOT:
            return result;
    }
    return result;
}

/** Expression nodes of the bonds of a pattern atom; a bond without a recorded count is one node. */
inline size_t bondNodesAround(const SmartsQuery& pattern, int atom) {
    if (static_cast<size_t>(atom) >= pattern.bondIndices.size()) return pattern.neighbors[atom].size();
    size_t nodes = 0;
    for (int bond : pattern.bondIndices[atom]) nodes += pattern.bondNodesOf(bond);
    return nodes;
}

/**
 * Data derived once from a compiled SMARTS pattern. It is immutable and may be
 * shared between threads; the pattern must outlive it and must not change.
 */
class PatternPlan {
public:
    /** Bond lookup, breadth-first order from atom 0 with tree parents, and stereo flag. */
    struct Layout {
        PatternBonds bonds;
        std::vector<int> order, parent;
        bool stereo;
        /** The layout reads only the bonds, the neighbour lists and the stereo flag of the pattern. */
        explicit Layout(const SmartsQuery& pattern, bool bondTable = true)
            : bonds(pattern, bondTable), parent(pattern.atomCount(), -1), stereo(pattern.checksMapping()) {
            const int count = pattern.atomCount();
            std::vector<uint8_t> seen(count, 0);
            order.reserve(count);
            for (int root = 0; root < count; ++root) {
                if (seen[root]) continue;
                seen[root] = 1;
                size_t head = order.size();
                order.push_back(root);
                for (; head < order.size(); ++head) {
                    const int atom = order[head];
                    for (int neighbor : pattern.neighbors[atom])
                        if (!seen[neighbor]) {
                            seen[neighbor] = 1;
                            parent[neighbor] = atom;
                            order.push_back(neighbor);
                        }
                }
            }
        }
    };

    explicit PatternPlan(const SmartsQuery& pattern)
        : layout_(pattern), shape_(pattern.atomCount(), pattern.neighbors) {
        const int count = pattern.atomCount();
        // Atoms with equal expressions and degrees share one domain.
        sameAs_.assign(count, -1);
        std::vector<std::pair<uint64_t, int>> keys(count);
        for (int atom = 0; atom < count; ++atom)
            keys[atom] = {atomExprHash(pattern.atoms[atom].expr) * 31 + pattern.neighbors[atom].size(), atom};
        std::sort(keys.begin(), keys.end());
        for (size_t first = 0; first < keys.size();) {
            size_t last = first;
            while (last < keys.size() && keys[last].first == keys[first].first) ++last;
            for (size_t index = first + 1; index < last; ++index) {
                const int atom = keys[index].second;
                for (size_t earlier = first; earlier < index; ++earlier) {
                    const int representative = keys[earlier].second;
                    if (sameAs_[representative] < 0
                        && pattern.neighbors[representative].size() == pattern.neighbors[atom].size()
                        && sameAtomExpr(pattern.atoms[representative].expr, pattern.atoms[atom].expr)) {
                        sameAs_[atom] = representative;
                        break;
                    }
                }
            }
            first = last;
        }
        // Element counts bound every injective embedding; requirements also order domain work.
        std::vector<int> elements(count, -1);
        elementSlot_.fill(-1);
        auto require = [&](int element, int aromatic) {
            if (element >= 0 && elementSlot_[element] < 0) {
                if (slots_ == MAX_SLOTS) return;
                elementSlot_[element] = static_cast<int8_t>(slots_++);
            }
            for (auto& known : counts_)
                if (known.element == element && known.aromatic == aromatic) { ++known.atoms; return; }
            counts_.push_back({element, aromatic, 1});
        };
        atomSlot_.assign(count, -1);
        atomAromatic_.assign(count, -1);
        // Expression nodes weigh atom tests and the bond tests of a candidate.
        atomNodes_.assign(count, 1);
        bondWork_.assign(count, 0);
        for (int atom = 0; atom < count; ++atom) {
            atomNodes_[atom] = pattern.atomNodesOf(atom);
            bondWork_[atom] = bondNodesAround(pattern, atom);
        }
        for (int atom = 0; atom < count; ++atom) {
            const auto requirement = atomRequirement(pattern.atoms[atom].expr);
            elements[atom] = requirement.element;
            atomAromatic_[atom] = static_cast<int8_t>(requirement.aromatic);
            if (requirement.element >= 0) {
                require(requirement.element, -1);
                atomSlot_[atom] = elementSlot_[requirement.element];
                if (requirement.aromatic >= 0) require(requirement.element, requirement.aromatic);
            }
            if (requirement.aromatic >= 0) require(-1, requirement.aromatic);
        }
        auto score = [&](int atom) { return elements[atom] < 0 ? 0 : elements[atom] == 6 ? 1 : 2; };
        domainOrder_.reserve(count);
        for (int atom = 0; atom < count; ++atom) if (sameAs_[atom] < 0) domainOrder_.push_back(atom);
        std::stable_sort(domainOrder_.begin(), domainOrder_.end(), [&](int first, int second) {
            return score(first) > score(second);
        });
        for (int atom = 0; atom < count; ++atom) if (sameAs_[atom] >= 0) domainOrder_.push_back(atom);
        addSubpatterns(pattern);
    }
    PatternPlan(const PatternPlan&) = delete;
    PatternPlan& operator=(const PatternPlan&) = delete;

    const Layout& layout() const { return layout_; }
    const std::vector<int>& domainOrder() const { return domainOrder_; }
    const std::vector<int>& sameAs() const { return sameAs_; }
    /** Selected bonds of a complete mapping in Match order, and its fragments as makeMatch() lists them. */
    const CompleteMappingShape& shape() const { return shape_; }
    /** The number of recursive subpatterns compiled into this pattern, at every depth. */
    size_t subpatternCount() const { return subpatterns_.size(); }
    /** The index, below subpatternCount(), of a recursive subpattern compiled into this pattern, or -1. */
    int subpatternIndex(const SmartsQuery& pattern) const {
        const auto found = std::lower_bound(subpatterns_.begin(), subpatterns_.end(), &pattern,
            [](const Subpattern& item, const SmartsQuery* key) { return std::less<const SmartsQuery*>()(item.first, key); });
        return found != subpatterns_.end() && found->first == &pattern ? static_cast<int>(found - subpatterns_.begin()) : -1;
    }
    /** The layout of the subpattern at an index from subpatternIndex(). */
    const Layout& subpatternLayout(int index) const { return *layouts_[subpatterns_[static_cast<size_t>(index)].second]; }
    /** The layout of a recursive subpattern compiled into this pattern, if any. */
    const Layout* subpattern(const SmartsQuery& pattern) const {
        const int index = subpatternIndex(pattern);
        return index < 0 ? nullptr : &subpatternLayout(index);
    }
    /** Index of the required-element group of an atom, or -1 when its element is open. */
    int atomSlot(int atom) const { return atomSlot_[atom]; }
    /** Required aromaticity of an atom: 1 aromatic, 0 aliphatic, -1 open. */
    int atomAromatic(int atom) const { return atomAromatic_[atom]; }
    /** Expression nodes of an atom, which a test of it evaluates at most. */
    size_t atomNodes(int atom) const { return atomNodes_[atom]; }
    /** Expression nodes of the bonds of an atom, which the bond tests of a candidate evaluate at most. */
    size_t bondWork(int atom) const { return bondWork_[atom]; }
    /**
     * Lists target atoms by required-element group (ascending within each group,
     * group g in atoms[start[g], start[g + 1])). False when the target has too few
     * atoms of a required element or aromaticity for any embedding.
     */
    bool groupTarget(const MolGraph& target, std::vector<int>& start, std::vector<int>& atoms) const {
        int total[2] = {0, 0};
        int elements[MAX_SLOTS][2] = {};
        for (int atom = 0; atom < target.n; ++atom) {
            const int aromatic = target.aromatic[atom] != 0 ? 1 : 0;
            ++total[aromatic];
            const int slot = slotOf(target.atomicNum[atom]);
            if (slot >= 0) ++elements[slot][aromatic];
        }
        for (const auto& count : counts_) {
            int available;
            if (count.element < 0) available = total[count.aromatic];
            else {
                const auto& slot = elements[elementSlot_[count.element]];
                available = count.aromatic < 0 ? slot[0] + slot[1] : slot[count.aromatic];
            }
            if (available < count.atoms) return false;
        }
        start.assign(static_cast<size_t>(slots_) + 1, 0);
        for (int slot = 0; slot < slots_; ++slot) start[slot + 1] = start[slot] + elements[slot][0] + elements[slot][1];
        atoms.resize(static_cast<size_t>(start[slots_]));
        int next[MAX_SLOTS];
        std::copy(start.begin(), start.begin() + slots_, next);
        for (int atom = 0; atom < target.n; ++atom) {
            const int slot = slotOf(target.atomicNum[atom]);
            if (slot >= 0) atoms[next[slot]++] = atom;
        }
        return true;
    }

private:
    static constexpr int MAX_SLOTS = 16;
    /** The most bytes the bond tables of the subpattern layouts take in all; later layouts use a map. */
    static constexpr size_t SUBPATTERN_TABLE_BYTES = size_t(1) << 24;
    struct Count { int element, aromatic, atoms; };
    /** A hash of what a layout reads from a pattern, see sameShape(). */
    static uint64_t shapeHash(const SmartsQuery& pattern, bool stereo) {
        uint64_t hash = (static_cast<uint64_t>(pattern.atomCount()) << 1) | (stereo ? 1u : 0u);
        auto mix = [&hash](uint64_t value) { hash = (hash ^ value) * 0x100000001B3ull; };
        for (const auto& bond : pattern.bonds) mix((static_cast<uint64_t>(uint32_t(bond.from)) << 32) | uint32_t(bond.to));
        return hash;
    }
    /** Whether two patterns give equal layouts: equal atom counts, bonds in the same order and neighbour lists. */
    static bool sameShape(const SmartsQuery& first, const SmartsQuery& second) {
        if (first.atomCount() != second.atomCount() || first.bonds.size() != second.bonds.size()
            || first.neighbors != second.neighbors) return false;
        for (size_t bond = 0; bond < first.bonds.size(); ++bond)
            if (first.bonds[bond].from != second.bonds[bond].from || first.bonds[bond].to != second.bonds[bond].to) return false;
        return true;
    }
    int slotOf(int element) const {
        return element >= 0 && element < static_cast<int>(elementSlot_.size()) ? elementSlot_[element] : -1;
    }
    /**
     * Compiles every distinct recursive subpattern once, then sorts them by address for lookup.
     * Subpatterns of one shape, such as the repeats of one subpattern in a long pattern, share
     * one layout, and the bond tables of all the layouts take at most SUBPATTERN_TABLE_BYTES.
     */
    void addSubpatterns(const SmartsQuery& pattern) {
        std::vector<const smarts_detail::AtomExpr*> pending;
        std::unordered_set<const SmartsQuery*> compiled;
        // The first subpattern of each layout, by the hash of its shape.
        std::unordered_multimap<uint64_t, std::pair<const SmartsQuery*, size_t>> shapes;
        size_t tableBytes = 0;
        for (const auto& atom : pattern.atoms) pending.push_back(&atom.expr);
        while (!pending.empty()) {
            const auto* expression = pending.back();
            pending.pop_back();
            for (const auto& child : expression->children) pending.push_back(&child);
            const auto* nested = expression->cachedRecursive.get();
            if (expression->primType != smarts_detail::AtomPrimType::RECURSIVE || !nested || !compiled.insert(nested).second)
                continue;
            const bool stereo = nested->checksMapping();
            const uint64_t hash = shapeHash(*nested, stereo);
            size_t layout = layouts_.size();
            const auto alike = shapes.equal_range(hash);
            for (auto shape = alike.first; shape != alike.second; ++shape)
                if (layouts_[shape->second.second]->stereo == stereo && sameShape(*shape->second.first, *nested)) {
                    layout = shape->second.second;
                    break;
                }
            if (layout == layouts_.size()) {
                const size_t bytes = PatternBonds::tableBytes(nested->atomCount());
                const bool table = bytes <= SUBPATTERN_TABLE_BYTES - tableBytes;
                if (table) tableBytes += bytes;
                layouts_.push_back(std::make_unique<const Layout>(*nested, table));
                shapes.emplace(hash, std::make_pair(nested, layout));
            }
            subpatterns_.emplace_back(nested, layout);
            for (const auto& atom : nested->atoms) pending.push_back(&atom.expr);
        }
        std::sort(subpatterns_.begin(), subpatterns_.end(), [](const Subpattern& first, const Subpattern& second) {
            return std::less<const SmartsQuery*>()(first.first, second.first);
        });
    }
    Layout layout_;
    CompleteMappingShape shape_;
    std::vector<int> sameAs_, domainOrder_;
    std::array<int8_t, 256> elementSlot_{};
    std::vector<int> atomSlot_;
    std::vector<int8_t> atomAromatic_;
    std::vector<size_t> atomNodes_, bondWork_;
    int slots_ = 0;
    std::vector<Count> counts_;
    /** A compiled subpattern and the index of its layout. */
    using Subpattern = std::pair<const SmartsQuery*, size_t>;
    std::vector<Subpattern> subpatterns_;
    std::vector<std::unique_ptr<const Layout>> layouts_;
};

struct PatternFrame {
    int atom;
    size_t begin, end, next;
    int active;
};

/** Search buffers reused by successive pattern searches on one thread. */
struct PatternWorkspace {
    std::vector<uint64_t> domains;
    std::vector<int> domainCounts, mapping, candidates, found, groupStart, groupAtoms;
    /** The pattern atom mapped to each used target atom of an induced search; read only where used is set. */
    std::vector<int> owners;
    std::vector<uint8_t> used;
    std::vector<PatternFrame> frames;
    bool busy = false;
    /** Marks the buffers as in use and releases unusually large ones afterwards. */
    class Hold {
        PatternWorkspace& workspace_;
        /** Frees a buffer of more than 2^16 items; swapping with an empty one allocates nothing. */
        template<class T>
        static void trim(std::vector<T>& buffer) noexcept {
            if (buffer.capacity() > (size_t(1) << 16)) std::vector<T>().swap(buffer);
        }
    public:
        explicit Hold(PatternWorkspace& workspace) : workspace_(workspace) { workspace_.busy = true; }
        ~Hold() {
            trim(workspace_.domains);
            trim(workspace_.used);
            for (auto* buffer : {&workspace_.domainCounts, &workspace_.mapping, &workspace_.candidates, &workspace_.found,
                                 &workspace_.groupStart, &workspace_.groupAtoms, &workspace_.owners})
                trim(*buffer);
            trim(workspace_.frames);
            workspace_.busy = false;
        }
        Hold(const Hold&) = delete;
        Hold& operator=(const Hold&) = delete;
    };
};

/**
 * Recursive SMARTS state for one target during one search call. Atom tests and
 * anchored results depend only on the target atom, so they are shared by every
 * anchor and every query atom that refers to the same compiled subpattern.
 *
 * Buffers indexed by target atom come from a byte allowance: a subpattern's
 * anchored results and mapped-atom marks when it is first met, and the tests of
 * each of its atoms when that atom is first tested. What the allowance does not
 * cover is evaluated without a cache, which gives the same answers, so many large
 * subpatterns on a large target cannot take memory in proportion to both.
 */
class RecursiveMemo {
public:
    /** The bytes of per-target-atom buffers that one search call, or one uncached evaluation, may take. */
    static constexpr size_t allowanceBytes = size_t(1) << 24;
    /** The allowance of the searches on this thread; tests lower it to compare uncached with cached answers. */
    static size_t& threadAllowance() noexcept {
        thread_local size_t allowance = allowanceBytes;
        return allowance;
    }
    /** Takes bytes from an allowance; false, with the allowance unchanged, when it holds fewer. */
    static bool take(size_t& allowance, size_t bytes) noexcept {
        if (bytes > allowance) return false;
        allowance -= bytes;
        return true;
    }
    struct Entry {
        const SmartsQuery* pattern = nullptr;
        const PatternPlan::Layout* layout = nullptr;
        std::optional<PatternPlan::Layout> ownLayout;
        int targetAtoms = 0;
        /** Per target atom: -1 unknown, 0 absent, 1 present. Empty when anchored results are not kept. */
        std::vector<int8_t> anchored;
        /** Per pattern atom: empty until the atom is first tested, then per target atom as above. */
        std::vector<std::vector<int8_t>> atomTests;
        std::vector<int> mapping, cursor;
        /** Per target atom, whether it is mapped. Empty when mapped atoms are looked up in mapping. */
        std::vector<uint8_t> used;
        Entry() = default;
        Entry(const Entry&) = delete;
        Entry& operator=(const Entry&) = delete;
        /**
         * Prepares the entry for a pattern and target size, taking its target-sized buffers
         * from allowance; anchored results are kept only when memoised. Buffers keep their
         * capacity.
         */
        void reset(const SmartsQuery& query, int atoms, const PatternPlan::Layout* shared, bool memoised, size_t& allowance) {
            const int count = query.atomCount();
            const size_t bytes = static_cast<size_t>(std::max(atoms, 0));
            pattern = &query;
            ownLayout.reset();
            if (!shared) shared = &ownLayout.emplace(query);
            layout = shared;
            targetAtoms = atoms;
            anchored.clear();
            if (memoised && take(allowance, bytes)) anchored.assign(bytes, -1);
            used.clear();
            if (take(allowance, bytes)) used.assign(bytes, 0);
            if (atomTests.size() < static_cast<size_t>(count)) atomTests.resize(static_cast<size_t>(count));
            for (auto& row : atomTests) row.clear();
            mapping.assign(count, -1);
            cursor.assign(count, 0);
        }
        /** The test results of a pattern atom, made from allowance when first needed, or null when it is spent. */
        int8_t* tests(int atom, size_t& allowance) {
            auto& row = atomTests[static_cast<size_t>(atom)];
            if (row.empty()) {
                if (!take(allowance, static_cast<size_t>(targetAtoms))) return nullptr;
                row.assign(static_cast<size_t>(targetAtoms), -1);
            }
            return row.data();
        }
        /** Whether a target atom is the image of a mapped pattern atom. */
        bool isUsed(int targetAtom) const {
            if (!used.empty()) return used[static_cast<size_t>(targetAtom)] != 0;
            return std::find(mapping.begin(), mapping.end(), targetAtom) != mapping.end();
        }
        /** Marks a target atom as mapped or free, where marks are kept. */
        void mark(int targetAtom, bool mapped) {
            if (!used.empty()) used[static_cast<size_t>(targetAtom)] = mapped ? 1 : 0;
        }
        size_t footprint() const noexcept {
            size_t bytes = anchored.capacity() + used.capacity() + atomTests.capacity() * sizeof(std::vector<int8_t>);
            for (const auto& row : atomTests) bytes += row.capacity();
            return bytes;
        }
        /** Frees every buffer; swapping with empty buffers allocates nothing. */
        void release() noexcept {
            pattern = nullptr;
            layout = nullptr;
            ownLayout.reset();
            targetAtoms = 0;
            std::vector<int8_t>().swap(anchored);
            std::vector<std::vector<int8_t>>().swap(atomTests);
            std::vector<int>().swap(mapping);
            std::vector<int>().swap(cursor);
            std::vector<uint8_t>().swap(used);
        }
    };
    /** Entries reused by successive calls on one thread, and the entry of each plan subpattern. */
    struct Pool {
        std::vector<std::unique_ptr<Entry>> entries;
        std::vector<int> slots;
        bool busy = false;
    };
    RecursiveMemo(const MolGraph& target, Pool& pool, const PatternPlan* plan)
        : target_(&target), pool_(pool), plan_(plan) { pool_.busy = true; }
    ~RecursiveMemo() noexcept {
        // Release buffers sized for unusually large targets or patterns, without allocating.
        constexpr size_t retained = size_t(1) << 16;
        for (auto& entry : pool_.entries)
            if (entry->footprint() > retained) entry->release();
        if (pool_.entries.size() > 32) pool_.entries.erase(pool_.entries.begin() + 32, pool_.entries.end());
        if (pool_.slots.capacity() > retained) std::vector<int>().swap(pool_.slots);
        pool_.busy = false;
    }
    RecursiveMemo(const RecursiveMemo&) = delete;
    RecursiveMemo& operator=(const RecursiveMemo&) = delete;
    const MolGraph* target() const { return target_; }
    /** What remains of this call's allowance. */
    size_t& allowance() { return allowance_; }
    /** The entry of a subpattern; plan subpatterns are found by index, any other by a scan of the entries. */
    Entry& entry(const SmartsQuery& pattern) {
        const int index = plan_ ? plan_->subpatternIndex(pattern) : -1;
        if (index >= 0) {
            if (!slotsReady_) {
                pool_.slots.assign(plan_->subpatternCount(), -1);
                slotsReady_ = true;
            }
            const int slot = pool_.slots[static_cast<size_t>(index)];
            if (slot >= 0) return *pool_.entries[static_cast<size_t>(slot)];
        } else {
            for (size_t slot = 0; slot < used_; ++slot)
                if (pool_.entries[slot]->pattern == &pattern) return *pool_.entries[slot];
        }
        if (used_ == pool_.entries.size()) pool_.entries.push_back(std::make_unique<Entry>());
        if (index >= 0) pool_.slots[static_cast<size_t>(index)] = static_cast<int>(used_);
        Entry& entry = *pool_.entries[used_++];
        entry.reset(pattern, target_->n, index >= 0 ? &plan_->subpatternLayout(index) : nullptr, true, allowance_);
        return entry;
    }
private:
    const MolGraph* target_;
    Pool& pool_;
    const PatternPlan* plan_;
    size_t used_ = 0;
    bool slotsReady_ = false;
    size_t allowance_ = threadAllowance();
};

inline thread_local RecursiveMemo* recursiveMemo = nullptr;

class RecursiveMemoScope {
    RecursiveMemo* previous_ = recursiveMemo;
public:
    explicit RecursiveMemoScope(RecursiveMemo& memo) { recursiveMemo = &memo; }
    ~RecursiveMemoScope() { recursiveMemo = previous_; }
    RecursiveMemoScope(const RecursiveMemoScope&) = delete;
    RecursiveMemoScope& operator=(const RecursiveMemoScope&) = delete;
};

/**
 * Whether a non-induced embedding maps pattern atom 0 to the anchor. Only
 * existence is reported, so the visiting order does not affect the result.
 * Atom tests are kept while allowance covers them.
 */
inline bool anchoredPatternExists(RecursiveMemo::Entry& entry, const MolGraph& target, int anchor, Budget& budget,
                                  size_t& allowance) {
    const SmartsQuery& pattern = *entry.pattern;
    const auto& layout = *entry.layout;
    const int count = pattern.atomCount(), targetAtoms = entry.targetAtoms;
    if (count == 0) return true;
    if (count > targetAtoms) return false;
    auto test = [&](int atom, int targetAtom) {
        int8_t* const known = entry.tests(atom, allowance);
        if (known && known[targetAtom] >= 0) return known[targetAtom] != 0;
        const bool matched = pattern.neighbors[atom].size() <= target.neighbors[targetAtom].size()
            && patternAtomMatches(pattern, atom, target, targetAtom);
        if (known && budget.status() == SearchStatus::Completed) known[targetAtom] = matched;
        return matched;
    };
    auto feasible = [&](int atom, int targetAtom) {
        if (entry.isUsed(targetAtom) || !test(atom, targetAtom)) return false;
        for (int neighbor : pattern.neighbors[atom]) {
            const int image = entry.mapping[neighbor];
            if (image < 0) continue;
            if (!target.hasBond(targetAtom, image)
                || !pattern.bondCompat(layout.bonds.find(atom, neighbor), target, targetAtom, image)) return false;
        }
        return true;
    };
    if (!feasible(0, anchor)) return false;
    bool found = false, stopped = false;
    int depth = 0;
    entry.cursor[0] = 0;
    // A candidate test evaluates at most the nodes of the atom's expression and bonds; the
    // clock is read after about 4096 of them, however long the expressions are.
    size_t work = 0;
    while (depth >= 0) {
        const int atom = layout.order[depth];
        int& image = entry.mapping[atom];
        if (image >= 0) { entry.mark(image, false); image = -1; }
        const int parent = layout.parent[atom];
        const int available = depth == 0 ? 1
            : parent >= 0 ? static_cast<int>(target.neighbors[entry.mapping[parent]].size()) : targetAtoms;
        const size_t testWork = 1 + pattern.atomNodesOf(atom) + bondNodesAround(pattern, atom);
        int candidate = -1;
        while (entry.cursor[depth] < available) {
            if ((work += testWork) >= 4096) {
                work = 0;
                if (budget.stopped()) { stopped = true; break; }
            }
            const int index = entry.cursor[depth]++;
            const int next = depth == 0 ? anchor : parent >= 0 ? target.neighbors[entry.mapping[parent]][index] : index;
            if (feasible(atom, next)) { candidate = next; break; }
        }
        if (stopped) break;
        if (candidate < 0) { --depth; continue; }
        if (!budget.visitSampled()) break;
        image = candidate;
        entry.mark(candidate, true);
        if (depth + 1 < count) {
            entry.cursor[++depth] = 0;
            continue;
        }
        if (!layout.stereo || pattern.mappingCompatible(target, entry.mapping)) { found = true; break; }
    }
    for (int atom = 0; atom < count; ++atom)
        if (entry.mapping[atom] >= 0) { entry.mark(entry.mapping[atom], false); entry.mapping[atom] = -1; }
    return found;
}

inline bool recursiveMatch(const SmartsQuery& pattern, const MolGraph& target, int anchor) {
    if (!recursiveBudget || recursiveBudget->stoppedSampled()) return false;
    Budget& budget = *recursiveBudget;
    if (recursiveMemo && recursiveMemo->target() == &target && smarts_detail::transientRecursionDepth == 0) {
        auto& entry = recursiveMemo->entry(pattern);
        size_t& allowance = recursiveMemo->allowance();
        if (entry.anchored.empty()) return anchoredPatternExists(entry, target, anchor, budget, allowance);
        int8_t& known = entry.anchored[anchor];
        if (known >= 0) return known != 0;
        const bool found = anchoredPatternExists(entry, target, anchor, budget, allowance);
        if (budget.status() == SearchStatus::Completed) known = found;
        return found;
    }
    size_t allowance = RecursiveMemo::threadAllowance();
    RecursiveMemo::Entry local;
    local.reset(pattern, target.n, nullptr, false, allowance);
    return anchoredPatternExists(local, target, anchor, budget, allowance);
}

/**
 * Substructure search for a SMARTS pattern. Atoms and candidates are visited in
 * VertexSearch order, so first matches and truncated match sets are those VertexSearch gives.
 */
class PatternSearch {
    using Frame = PatternFrame;
    const SmartsQuery& pattern_;
    const PatternPlan& plan_;
    const MolGraph& target_;
    const ChemOptions& chemistry_;
    bool induced_, all_;
    Budget& budget_;
    const PatternBonds& bonds_;
    RingConstraint targetRings_;
    int queryAtoms_, targetAtoms_;
    size_t words_;
    bool stereo_, emptyDomain_ = false;
    std::vector<uint64_t>& domains_;
    std::vector<int>& domainCounts_;
    std::vector<int>& mapping_;
    std::vector<int>& owners_;
    std::vector<uint8_t>& used_;
    int mapped_ = 0;
    std::vector<Frame>& frames_;
    std::vector<int>& candidates_;
    /** Accepted mappings in visiting order. */
    CompleteMappingRows rows_;
    /** Whether the budget stopped a frame before all its candidates were tested. */
    bool frameStopped_ = false;

    bool domainHas(int atom, int targetAtom) const {
        return (domains_[static_cast<size_t>(atom) * words_ + (targetAtom >> 6)] & (uint64_t(1) << (targetAtom & 63))) != 0;
    }
    bool bondMatches(int first, int second, int mappedFirst, int mappedSecond) const {
        const int bond = bonds_.find(first, second);
        return bond >= 0 && pattern_.bondCompat(bond, target_, mappedFirst, mappedSecond);
    }
    /**
     * The atoms an induced test reads beyond the pattern neighbours: the mapped atoms or the
     * target atom's neighbours, whichever are fewer.
     */
    size_t inducedReads(int targetAtom) const {
        return std::min(static_cast<size_t>(mapped_), target_.neighbors[targetAtom].size());
    }
    bool feasible(int atom, int targetAtom) const {
        if (used_[targetAtom] || !domainHas(atom, targetAtom)) return false;
        // Each mapped pattern neighbour needs a matching target bond.
        for (int neighbor : pattern_.neighbors[atom]) {
            const int image = mapping_[neighbor];
            if (image >= 0 && (!target_.hasBond(targetAtom, image) || !bondMatches(atom, neighbor, targetAtom, image)))
                return false;
        }
        if (!induced_) return true;
        // An induced embedding also has no target bond between the candidate and the image of a
        // mapped atom that is not a pattern neighbour; the shorter of the two lists is read.
        const auto& adjacent = target_.neighbors[targetAtom];
        if (static_cast<size_t>(mapped_) < adjacent.size()) {
            for (const Frame& frame : frames_)
                if (frame.active >= 0 && bonds_.find(atom, frame.atom) < 0 && target_.hasBond(targetAtom, frame.active))
                    return false;
            return true;
        }
        for (int neighbor : adjacent)
            if (used_[neighbor] && bonds_.find(atom, owners_[neighbor]) < 0) return false;
        return true;
    }
    std::vector<std::array<int, 4>> mappedBonds(const int* mapping) const {
        std::vector<std::array<int, 4>> bonds;
        for (int atom = 0; atom < queryAtoms_; ++atom) {
            if (mapping[atom] < 0) continue;
            for (int neighbor : pattern_.neighbors[atom])
                if (atom < neighbor && mapping[neighbor] >= 0
                    && bondMatches(atom, neighbor, mapping[atom], mapping[neighbor]))
                    bonds.push_back({atom, neighbor, mapping[atom], mapping[neighbor]});
        }
        return bonds;
    }
    bool record() {
        if (mapped_ != queryAtoms_) return false;
        if (stereo_ && !pattern_.mappingCompatible(target_, mapping_)) return false;
        if (chemistry_.completeRingsOnly) {
            std::set<int> atoms(mapping_.begin(), mapping_.end());
            std::set<int64_t> bonds;
            for (const auto& bond : mappedBonds(mapping_.data())) bonds.insert(MolGraph::bondKey(bond[2], bond[3]));
            if (!targetRings_.accepts(atoms, bonds)) return false;
        }
        rows_.append(mapping_);
        return true;
    }
    Frame frame() {
        // Choosing reads each pattern atom and its neighbours, and a candidate test evaluates at
        // most the expression nodes of the pattern atom's bonds and, when induced, reads the fewer
        // of the mapped atoms and the candidate's neighbours; the clock is read after about 4096
        // such units. Once the budget stops, the frame keeps the candidates found so far and the
        // search ends at its next node.
        size_t work = 0;
        auto stopAfter = [&](size_t units) {
            if (frameStopped_) return true;
            if ((work += units) < 4096) return false;
            work = 0;
            frameStopped_ = budget_.stopped();
            return frameStopped_;
        };
        int chosen = -1, bestFrontier = -1;
        for (int atom = 0; atom < queryAtoms_; ++atom) {
            if (stopAfter(1 + pattern_.neighbors[atom].size())) break;
            if (mapping_[atom] != -2) continue;
            int frontier = 0;
            for (int neighbor : pattern_.neighbors[atom]) frontier += mapping_[neighbor] >= 0;
            if (chosen < 0 || (frontier != 0 && bestFrontier == 0)
                || ((frontier != 0) == (bestFrontier != 0)
                    && (domainCounts_[atom] < domainCounts_[chosen]
                        || (domainCounts_[atom] == domainCounts_[chosen]
                            && pattern_.neighbors[atom].size() > pattern_.neighbors[chosen].size())))) {
                chosen = atom;
                bestFrontier = frontier;
            }
        }
        Frame result{chosen, candidates_.size(), candidates_.size(), candidates_.size(), -2};
        if (chosen < 0) return result;
        int neighborAnchor = -1;
        for (int neighbor : pattern_.neighbors[chosen])
            if (mapping_[neighbor] >= 0 && (neighborAnchor < 0
                || target_.neighbors[mapping_[neighbor]].size() < target_.neighbors[mapping_[neighborAnchor]].size()))
                neighborAnchor = neighbor;
        const size_t testCost = 1 + plan_.bondWork(chosen);
        auto test = [&](int candidate) {
            if (stopAfter(induced_ ? testCost + inducedReads(candidate) : testCost)) return false;
            if (feasible(chosen, candidate)) candidates_.push_back(candidate);
            return true;
        };
        if (neighborAnchor >= 0) {
            for (int candidate : target_.neighbors[mapping_[neighborAnchor]])
                if (!test(candidate)) break;
        } else {
            for (size_t word = 0; word < words_ && !frameStopped_; ++word) {
                uint64_t remaining = domains_[static_cast<size_t>(chosen) * words_ + word];
                while (remaining != 0) {
                    const int candidate = static_cast<int>(word * 64) + smsd::ctz64(remaining);
                    remaining &= remaining - 1;
                    if (!test(candidate)) break;
                }
            }
        }
        result.end = candidates_.size();
        return result;
    }
    /** Matches in Collector order; a single match is always built. */
    std::vector<Match> matches() const {
        std::vector<Match> result;
        if (rows_.count() <= 1) {
            if (rows_.count() == 1) {
                const std::vector<int> mapping(rows_.row(0), rows_.row(0) + queryAtoms_);
                result.push_back(makeMatch(mapping, mappedBonds(mapping.data())));
            }
            return result;
        }
        // A recorded mapping is complete, and each of its pairs of bonded atoms passed the bond
        // test when the later atom was mapped, so every bond is selected and bonds and fragments
        // follow the plan without testing the bonds again.
        return rows_.build(plan_.shape(), &budget_);
    }
public:
    PatternSearch(const SmartsQuery& pattern, const PatternPlan& plan, const MolGraph& target,
                  const ChemOptions& chemistry, bool induced, bool all, size_t limit, Budget& budget,
                  PatternWorkspace& workspace)
        : pattern_(pattern), plan_(plan), target_(target), chemistry_(chemistry), induced_(induced), all_(all),
          budget_(budget), bonds_(plan.layout().bonds), targetRings_(target, chemistry.completeRingsOnly),
          queryAtoms_(pattern.atomCount()), targetAtoms_(target.n),
          words_((static_cast<size_t>(target.n) + 63) >> 6), stereo_(plan.layout().stereo),
          domains_(workspace.domains), domainCounts_(workspace.domainCounts), mapping_(workspace.mapping),
          owners_(workspace.owners), used_(workspace.used), frames_(workspace.frames), candidates_(workspace.candidates),
          rows_(workspace.found, static_cast<size_t>(pattern.atomCount()), limit, all) {
        mapping_.clear();
        frames_.clear();
        candidates_.clear();
        if (queryAtoms_ == 0 || queryAtoms_ > targetAtoms_) return;
        // Too few atoms of a required element or aromaticity leaves no embedding.
        auto& groupStart = workspace.groupStart;
        auto& groupAtoms = workspace.groupAtoms;
        if (!plan.groupTarget(target, groupStart, groupAtoms)) { emptyDomain_ = true; return; }
        if (words_ != 0 && static_cast<size_t>(queryAtoms_) > std::numeric_limits<size_t>::max() / words_)
            throw std::length_error("Domain matrix is too large");
        domains_.assign(static_cast<size_t>(queryAtoms_) * words_, 0);
        domainCounts_.assign(queryAtoms_, 0);
        mapping_.assign(queryAtoms_, -2);
        used_.assign(targetAtoms_, 0);
        if (induced) owners_.resize(targetAtoms_);
        // Domains do not depend on the order in which they are filled; heteroatoms go
        // first because an empty domain ends the search.
        for (int atom : plan.domainOrder()) {
            if (budget.stoppedSampled()) break;
            uint64_t* domain = domains_.data() + static_cast<size_t>(atom) * words_;
            const int same = plan.sameAs()[atom];
            if (same >= 0) {
                std::copy_n(domains_.data() + static_cast<size_t>(same) * words_, words_, domain);
                domainCounts_[atom] = domainCounts_[same];
            } else {
                // An atom with a required element can only match target atoms of that element.
                const size_t degree = pattern.neighbors[atom].size();
                const int slot = plan.atomSlot(atom), aromatic = plan.atomAromatic(atom);
                const int first = slot < 0 ? 0 : groupStart[slot];
                const int last = slot < 0 ? targetAtoms_ : groupStart[slot + 1];
                // A candidate counts one unit and its test the nodes of the expression; the clock
                // is read after about 512 units, so a long expression reads it more often.
                const size_t testWork = plan.atomNodes(atom);
                size_t work = 0;
                for (int index = first; index < last; ++index) {
                    if (work >= 512) {
                        work = 0;
                        if (budget.stopped()) break;
                    }
                    ++work;
                    const int candidate = slot < 0 ? index : groupAtoms[index];
                    if (degree > target.neighbors[candidate].size()) continue;
                    if (aromatic >= 0 && (target.aromatic[candidate] != 0) != (aromatic != 0)) continue;
                    work += testWork;
                    if (patternAtomMatches(pattern, atom, target, candidate)) {
                        domain[candidate >> 6] |= uint64_t(1) << (candidate & 63);
                        ++domainCounts_[atom];
                    }
                }
            }
            if (domainCounts_[atom] == 0) { emptyDomain_ = true; break; }
        }
    }
    std::pair<bool, bool> run() {
        if (budget_.stopped()) return {false, false};
        if (queryAtoms_ == 0) {
            record();
            return {true, true};
        }
        if (queryAtoms_ > targetAtoms_ || emptyDomain_) return {true, true};
        frames_.reserve(queryAtoms_);
        frames_.push_back(frame());
        while (!frames_.empty() && !budget_.stoppedSampled()) {
            Frame& current = frames_.back();
            if (current.active >= 0) {
                used_[current.active] = 0;
                --mapped_;
                mapping_[current.atom] = -2;
                current.active = -2;
            }
            if (current.next < current.end) {
                if (!budget_.visitSampled()) break;
                const int candidate = candidates_[current.next++];
                mapping_[current.atom] = candidate;
                if (induced_) owners_[candidate] = current.atom;
                used_[candidate] = 1;
                ++mapped_;
                current.active = candidate;
                const bool accepted = record();
                if (accepted && !all_ && mapped_ == queryAtoms_) return {true, true};
                if (rows_.truncated) return {true, false};
                const Frame next = frame();
                if (next.atom >= 0) frames_.push_back(next);
            } else {
                candidates_.resize(current.begin);
                frames_.pop_back();
            }
        }
        // A frame the budget stopped was not searched through, so the stack running empty proves nothing then.
        const bool exhausted = frames_.empty() && !frameStopped_;
        return {exhausted || rows_.count() != 0, exhausted};
    }
    SearchResult finish(bool proven, bool exhausted) {
        SearchResult result;
        result.matches = matches();
        result.status = budget_.status();
        if (result.status == SearchStatus::Completed && rows_.truncated) result.status = SearchStatus::LimitReached;
        result.optimal = proven;
        result.complete = exhausted && !rows_.truncated && result.status == SearchStatus::Completed;
        result.score = rows_.count() != 0 ? double(queryAtoms_) : 0.0;
        result.nodesVisited = budget_.nodesVisited();
        result.elapsedMillis = budget_.elapsedMillis();
        return result;
    }
};

} // namespace search_detail
} // namespace smsd
