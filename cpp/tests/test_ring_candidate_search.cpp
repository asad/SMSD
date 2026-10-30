/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 *
 * The ring candidate search read from one tree per root against the frozen pairwise search:
 * every pair's paths, validity and counted work, every subtree's closing bond, every ring
 * result at and around the work limit, the internal switches, and the relevant-cycle families.
 *
 *   smsd_ring_candidate_search_tests                  the fixtures; results and families in two
 *                                                     neighbour orders, so that a Debug build stays short
 *   smsd_ring_candidate_search_tests full             the fixtures in every order with every option set
 *   smsd_ring_candidate_search_tests corpus FILE...   also each SMILES file (SMILES, then an id)
 *   smsd_ring_candidate_search_tests check-golden FILE  checks a golden file as the shared one is checked;
 *                                                     the Java ring tests write the shared file
 */
#include "smsd/smiles_parser.hpp"
#include "smsd/ring_finder.hpp"
#include "smsd/relevant_cycle_families.hpp"
#include "support/ring_finder_reference.hpp"

#include <algorithm>
#include <cctype>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <map>
#include <random>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

using smsd::MolGraph;
using smsd::RingResult;
using smsd::detail_rf::RingCandidateSearch;
using smsd::detail_rf::RingSearchNotes;
using smsd::detail_rf::RingSearchOptions;
using Kind = RingCandidateSearch::PairKind;
using Adjacency = std::vector<std::vector<int>>;

int failures = 0;
/** Every neighbour order and option set in the result and family controls, for the full run. */
bool exhaustive = false;

void require(bool value, const std::string& message) { if (!value) throw std::runtime_error(message); }

template<class Operation> void control(const std::string& name, Operation operation) {
    try { operation(); std::cout << "PASS\t" << name << '\n'; }
    catch (const std::exception& error) { ++failures; std::cout << "FAIL\t" << name << '\t' << error.what() << '\n'; }
}

MolGraph fromNeighbours(Adjacency neighbours) {
    const int n = static_cast<int>(neighbours.size());
    return MolGraph::Builder().atomCount(n).atomicNumbers(std::vector<int>(n, 6))
        .setNeighbors(std::move(neighbours)).build(false);
}

/** A portable seeded generator: the standard engine with explicit reduction. */
struct Seeded {
    std::mt19937 engine;
    explicit Seeded(uint32_t seed) : engine(seed) {}
    size_t below(size_t bound) { return static_cast<size_t>(engine() % static_cast<uint32_t>(bound)); }
};

Adjacency shuffledNeighbours(const MolGraph& g, uint32_t seed) {
    Seeded random(seed);
    Adjacency neighbours = g.neighbors;
    for (auto& row : neighbours)
        for (size_t i = row.size(); i > 1; --i) std::swap(row[i - 1], row[random.below(i)]);
    return neighbours;
}

Adjacency relabelled(const MolGraph& g, uint32_t seed) {
    Seeded random(seed);
    std::vector<int> label(static_cast<size_t>(g.n));
    for (int i = 0; i < g.n; ++i) label[i] = i;
    for (size_t i = label.size(); i > 1; --i) std::swap(label[i - 1], label[random.below(i)]);
    Adjacency neighbours(static_cast<size_t>(g.n));
    for (int u = 0; u < g.n; ++u)
        for (int w : g.neighbors[u]) neighbours[label[u]].push_back(label[w]);
    return neighbours;
}

std::string cyclodextrin(int units) {
    std::string smiles = "C%101C(O)C(O)C(OC1CO)O";
    for (int unit = 1; unit < units; ++unit) smiles += "C1C(O)C(O)C(OC1CO)O";
    return smiles + "%10";
}

/** Hexagonal rings stacked into prisms, each joined to the next by six bonds. */
Adjacency prismStack(int layers) {
    Adjacency neighbours(static_cast<size_t>(6 * layers));
    auto bond = [&](int a, int b) { neighbours[a].push_back(b); neighbours[b].push_back(a); };
    for (int layer = 0; layer < layers; ++layer)
        for (int i = 0; i < 6; ++i) {
            bond(6 * layer + i, 6 * layer + (i + 1) % 6);
            if (layer + 1 < layers) bond(6 * layer + i, 6 * (layer + 1) + i);
        }
    return neighbours;
}

/** Four-member rings linked around one larger ring. */
Adjacency fourRingChain(int rings) {
    Adjacency neighbours(static_cast<size_t>(1 + 3 * rings));
    auto bond = [&](int a, int b) { neighbours[a].push_back(b); neighbours[b].push_back(a); };
    int previous = 0;
    for (int ring = 0; ring < rings; ++ring) {
        const int firstSide = 3 * ring + 1, secondSide = firstSide + 1, next = firstSide + 2;
        bond(previous, firstSide); bond(previous, secondSide);
        bond(firstSide, next); bond(secondSide, next);
        previous = next;
    }
    bond(previous, 0);
    return neighbours;
}

/** The graph whose root 0 and bond (0, 2) read the child from the level below it. */
Adjacency belowChildFixture() { return {{1, 2, 3}, {0, 3, 4}, {0, 4, 5}, {0, 1}, {1, 2, 5}, {2, 4}}; }

/** A random graph: cycles with chords, degree at most four, pendant chains and a spare component. */
Adjacency randomGraph(uint32_t seed) {
    Seeded random(seed);
    const int ringAtoms = 6 + static_cast<int>(random.below(55));
    Adjacency neighbours(static_cast<size_t>(ringAtoms));
    auto has = [&](int a, int b) { return std::find(neighbours[a].begin(), neighbours[a].end(), b) != neighbours[a].end(); };
    auto bond = [&](int a, int b) { neighbours[a].push_back(b); neighbours[b].push_back(a); };
    for (int i = 0; i < ringAtoms; ++i) bond(i, (i + 1) % ringAtoms);
    const int chords = 1 + static_cast<int>(random.below(static_cast<size_t>(ringAtoms / 2 + 1)));
    for (int attempt = 0, added = 0; attempt < 4 * chords && added < chords; ++attempt) {
        const int a = static_cast<int>(random.below(static_cast<size_t>(ringAtoms)));
        const int b = static_cast<int>(random.below(static_cast<size_t>(ringAtoms)));
        if (a == b || has(a, b) || neighbours[a].size() >= 4 || neighbours[b].size() >= 4) continue;
        bond(a, b);
        ++added;
    }
    const int pendants = static_cast<int>(random.below(4));
    for (int p = 0; p < pendants; ++p) {
        const int at = static_cast<int>(random.below(neighbours.size()));
        if (neighbours[at].size() >= 4) continue;
        const int atom = static_cast<int>(neighbours.size());
        neighbours.emplace_back();
        bond(at, atom);
    }
    if (random.below(2) == 0) {   // a separate triangle with a tail
        const int base = static_cast<int>(neighbours.size());
        neighbours.resize(neighbours.size() + 4);
        bond(base, base + 1); bond(base + 1, base + 2); bond(base + 2, base); bond(base + 2, base + 3);
    }
    for (auto& row : neighbours)
        for (size_t i = row.size(); i > 1; --i) std::swap(row[i - 1], row[random.below(i)]);
    return neighbours;
}

struct Fixture { std::string id; MolGraph graph; bool large = false; };

std::vector<Fixture> fixtures() {
    std::vector<Fixture> list;
    const std::pair<const char*, const char*> smiles[] = {
        {"naphthalene", "c1ccc2ccccc2c1"}, {"adamantane", "C1C2CC3CC1CC(C2)C3"},
        {"cubane", "C12C3C4C1C5C2C3C45"}, {"norbornane", "C1CC2CCC1C2"}, {"bicyclooctane", "C1CC2CCC1CC2"},
        {"coronene", "c1cc2ccc3ccc4ccc5ccc6ccc1c7c2c3c4c5c67"}, {"pentaprismane", "C12C3C4C5C1C6C2C3C4C56"},
        {"hexaprismane", "C12C3C4C5C6C1C7C2C3C4C5C67"}, {"cage", "C12C3C4C1C5C3C4C25"},
        {"spiro", "C1CCC2(CC1)CCCC2"}, {"bicyclononane", "C1CC2CCCC(C1)CC2"}, {"triptycene", "c1ccc2c(c1)C1c3ccccc3C2c2ccccc12"},
        {"paracyclophane", "c1cc2ccc1CCc1ccc(cc1)CC2"}, {"porphine", "c1cc2cc3ccc(cc4ccc(cc5ccc(cc1n2)[nH]5)n4)[nH]3"},
        {"dodecahedrane", "C12C3C4C5C1C6C7C2C8C3C9C4C%10C5C6C%11C7C8C9C%10%11"},
        {"cryptand", "C1COCCOCCN2CCOCCOCCN1CCOCCOCC2"}, {"phenol-chain", "Oc1ccccc1CCc1ccncc1"},
        {"cyclopentane", "C1CCCC1"}, {"benzene", "c1ccccc1"}, {"steroid", "CC12CCC3C(CCC4CC(O)CCC34C)C1CCC2O"},
    };
    for (const auto& [id, text] : smiles) list.push_back({id, smsd::parseSMILES(text), false});
    for (int units : {6, 7, 8}) list.push_back({"cd" + std::to_string(units), smsd::parseSMILES(cyclodextrin(units)), false});
    list.push_back({"cd20", smsd::parseSMILES(cyclodextrin(20)), true});
    for (int layers = 2; layers <= 8; ++layers) list.push_back({"stack" + std::to_string(layers), fromNeighbours(prismStack(layers)), false});
    list.push_back({"stack32", fromNeighbours(prismStack(32)), true});
    for (int rings : {1, 2, 5, 12, 13}) list.push_back({"cyclobutane-chain" + std::to_string(rings), fromNeighbours(fourRingChain(rings)), false});
    list.push_back({"below-child", fromNeighbours(belowChildFixture()), false});
    for (uint32_t seed = 1; seed <= 40; ++seed) list.push_back({"random" + std::to_string(seed), fromNeighbours(randomGraph(seed)), false});
    return list;
}

/**
 * The graph as given and with sorted neighbours, then, unless the graph is large or only those two
 * orders are asked for, three neighbour shuffles and three relabellings.
 */
std::vector<std::pair<std::string, MolGraph>> variants(const Fixture& fixture, bool everyOrder = true) {
    std::vector<std::pair<std::string, MolGraph>> list;
    list.push_back({"given", fixture.graph});
    Adjacency sorted = fixture.graph.neighbors;
    for (auto& row : sorted) std::sort(row.begin(), row.end());
    list.push_back({"sorted", fromNeighbours(sorted)});
    if (fixture.large || !everyOrder) return list;
    for (uint32_t seed : {11u, 12u, 13u}) list.push_back({"shuffle" + std::to_string(seed), fromNeighbours(shuffledNeighbours(fixture.graph, seed))});
    for (uint32_t seed : {21u, 22u, 23u}) list.push_back({"relabel" + std::to_string(seed), fromNeighbours(relabelled(fixture.graph, seed))});
    return list;
}

bool sameResult(const RingResult& a, const RingResult& b) {
    return a.sssr == b.sssr && a.rcb == b.rcb && a.legacyRcb == b.legacyRcb && a.truncated == b.truncated
        && a.basisTruncated == b.basisTruncated && a.relevantTruncated == b.relevantTruncated;
}

size_t weight(const std::vector<std::vector<int>>& cycles) {
    size_t total = 0;
    for (const auto& cycle : cycles) total += cycle.size();
    return total;
}

/** Statistics of the per-pair comparison. */
struct PairStatistics {
    uint64_t graphs = 0, pairs = 0, nonTree = 0, a = 0, b1 = 0, b1Plus = 0, b2 = 0, simple = 0, covers = 0;
    void print(const std::string& label) const {
        std::cout << "PAIRS\t" << label << "\tgraphs=" << graphs << "\tpairs=" << pairs << "\tnon-tree=" << nonTree
                  << "\tA=" << a << "\tB1=" << b1 << "\tB1+=" << b1Plus << "\tB2=" << b2 << "\tsimple=" << simple
                  << "\tcovers=" << covers << '\n';
    }
};

/** Every pair against the pairwise search, and every cover against the least closing half-bond. */
void comparePairs(const std::string& id, const MolGraph& g, PairStatistics& statistics) {
    const auto systems = g.cycleSystems();
    if (systems.ofAtom.empty() || smsd::detail::cycleRank(g) <= 0) return;
    const auto edgeIndex = smsd::detail_rf::buildEdgeIndex(g);
    const auto edgeList = smsd::detail_rf::buildEdgeList(g);
    RingSearchNotes notes;
    RingCandidateSearch search(g, systems, edgeIndex, edgeList, smsd::RING_CANDIDATE_WORK_LIMIT, RingSearchOptions(), notes);
    require(search.applicable(), id + ": the tree search applies");
    const auto records = smsd::ring_reference::replayPairs(g);
    require(notes.pairs == records.size(), id + ": the pair count");
    ++statistics.graphs;
    int built = -1;
    std::vector<int> first, second;
    for (const auto& record : records) {
        ++statistics.pairs;
        const int system = systems.ofAtom[record.root];
        if (systems.simpleRing[system]) { ++statistics.simple; continue; }
        if (built != record.root) {
            search.buildTreeForTesting(record.root);
            built = record.root;
            // Each atom's cover is the least closing half-bond, by length then tree order,
            // from its subtree to outside it.
            const auto closing = search.closingForTesting();
            for (int c = 0; c < g.n; ++c) {
                if (systems.ofAtom[c] != system || c == record.root) continue;
                std::pair<int, int> best{-1, -1};
                int bestLength = 0;
                for (const auto& half : closing) {
                    if (!search.inSubtreeForTesting(c, half.second) || search.inSubtreeForTesting(c, half.first)) continue;
                    const int length = search.depthForTesting(half.first) + search.depthForTesting(half.second) + 1;
                    if (best.first < 0 || length < bestLength) { best = half; bestLength = length; }
                }
                require(search.coverForTesting(c) == best, id + ": cover of atom " + std::to_string(c));
                ++statistics.covers;
            }
        }
        const auto reading = search.readPairForTesting(record.edge);
        search.pathsForTesting(reading, first, second);
        const std::string pair = id + ": root " + std::to_string(record.root) + " bond " + std::to_string(record.edge);
        require(first == record.first && second == record.second, pair + " paths");
        require(reading.valid == record.valid, pair + " validity");
        require(reading.pathTerm == record.first.size() + record.second.size(), pair + " path atoms");
        switch (reading.kind) {
            case Kind::NonTree: ++statistics.nonTree; break;
            case Kind::ClassA: ++statistics.a; break;
            case Kind::ClassB1: ++statistics.b1; break;
            case Kind::ClassB1Plus: ++statistics.b1Plus; break;
            case Kind::ClassB2: ++statistics.b2; break;
            case Kind::SimpleRing: break;
        }
        if (reading.kind == Kind::ClassB2)
            require(reading.searchLow <= record.search && record.search <= reading.searchHigh, pair + " bracket");
        else
            require(reading.searchLow == record.search && reading.searchHigh == record.search, pair + " exact work");
    }
}

/** Results against the pairwise search at the limit, at the exact work before the last pair and below it. */
void compareResults(const std::string& id, const MolGraph& g, bool boundary, bool allOptions) {
    const auto records = smsd::ring_reference::replayPairs(g);
    const auto totals = smsd::ring_reference::workTotals(records);
    std::vector<uint64_t> limits = {smsd::RING_CANDIDATE_WORK_LIMIT};
    if (boundary) {
        limits.push_back(totals.beforeLast);
        if (totals.beforeLast > 0) limits.push_back(totals.beforeLast - 1);
        limits.push_back(uint64_t(1) << 20);
        limits.push_back(0);
    }
    std::vector<RingSearchOptions> optionSets(1);
    if (allOptions) {
        optionSets.push_back(RingSearchOptions());
        optionSets.back().loosenBracketsForTesting = true;
        optionSets.push_back(RingSearchOptions());
        optionSets.back().replayAllForTesting = true;
        optionSets.push_back(optionSets.back());
        optionSets.back().loosenBracketsForTesting = true;
    }
    for (uint64_t limit : limits)
        for (bool relevant : {false, true}) {
            const auto base = smsd::ring_reference::perceiveRingSets(g, relevant, limit);
            for (size_t o = 0; o < optionSets.size(); ++o) {
                RingSearchNotes notes;
                const auto result = smsd::detail_rf::perceiveRingSets(g, relevant, limit, optionSets[o], &notes);
                const std::string run = id + " limit " + std::to_string(limit) + (relevant ? " relevant" : " basis")
                    + " options " + std::to_string(o);
                require(sameResult(result, base), run + ": result");
                require(!notes.lastPairSummed, run + ": the last pair is never summed");
                if (notes.pairwise) continue;
                require(notes.pairs == totals.pairs, run + ": pairs");
                require(notes.fellBack == (totals.pairs > 1 && totals.beforeLast > limit), run + ": fallback decision");
                if (!notes.fellBack || notes.baseWorkExact)
                    require(notes.baseWorkLow <= totals.beforeLast
                            && (notes.baseWorkHighUnbounded || totals.beforeLast <= notes.baseWorkHigh), run + ": bounds");
                if (optionSets[o].replayAllForTesting)
                    require(notes.baseWorkExact && notes.baseWorkLow == totals.beforeLast, run + ": replayed work is exact");
            }
        }
}

std::set<std::vector<int64_t>> keysOf(const std::vector<std::vector<int>>& cycles) {
    std::set<std::vector<int64_t>> keys;
    for (const auto& cycle : cycles) keys.insert(smsd::detail_rf::canonicalCycleKey(cycle));
    return keys;
}

/** Every member of every family, unless there are more than limit. */
bool familyMembers(const MolGraph& g, const smsd::RelevantCycleFamilies& families, size_t limit,
                   std::vector<std::vector<int>>& members) {
    for (const auto& family : families.families) {
        if (family.memberCount.size() > 9 || std::stoull(family.memberCount) > limit) return false;
        const auto list = smsd::relevantCycleFamilyMembers(g, family, limit);
        require(std::to_string(list.size()) == family.memberCount, "every member is enumerated");
        require(!list.empty() && list.front() == family.prototype, "members start with the prototype");
        for (const auto& member : list) {
            require(smsd::relevantCycleFamilyOf(g, families, member) == static_cast<int>(&family - families.families.data()),
                    "each member belongs to its own family");
            members.push_back(member);
        }
    }
    return true;
}

/** Families against the relevant cycles of the pairwise search, wherever that list is complete. */
void compareFamilies(const std::string& id, const MolGraph& g) {
    const auto families = smsd::detail_rf::computeRelevantCycleFamiliesUncached(g);
    const auto base = smsd::ring_reference::perceiveRingSets(g, true);
    if (families.truncated || base.relevantTruncated) return;
    std::vector<std::vector<int>> members;
    if (!familyMembers(g, families, 20000, members)) return;
    require(keysOf(members).size() == members.size(), id + ": members are distinct");
    require(keysOf(members) == keysOf(base.rcb), id + ": the family members are the relevant cycles");
    require(families.relevantCycleCount == std::to_string(base.rcb.size()), id + ": the exact count");
}

// ----- Golden file: a JSON value, only what the file uses -----

struct Json {
    enum class Kind { Null, Bool, Number, String, Array, Object } kind = Kind::Null;
    bool boolean = false;
    std::string text;   // a number's or string's text
    std::vector<Json> items;
    std::map<std::string, Json> members;
    const Json& operator[](const std::string& key) const {
        const auto found = members.find(key);
        if (found == members.end()) throw std::runtime_error("golden: missing key " + key);
        return found->second;
    }
    uint64_t unsignedValue() const { return std::stoull(text); }
    int integer() const { return std::stoi(text); }
    std::vector<int> integers() const { std::vector<int> list; for (const auto& item : items) list.push_back(item.integer()); return list; }
};

class JsonReader {
    const std::string& source_;
    size_t at_ = 0;
    void space() { while (at_ < source_.size() && std::isspace(static_cast<unsigned char>(source_[at_]))) ++at_; }
    char peek() { space(); if (at_ >= source_.size()) throw std::runtime_error("golden: unexpected end"); return source_[at_]; }
    void expect(char c) { if (peek() != c) throw std::runtime_error(std::string("golden: expected ") + c); ++at_; }
    std::string string() {
        expect('"');
        std::string value;
        while (source_.at(at_) != '"') value += source_[at_++];
        ++at_;
        return value;
    }
public:
    explicit JsonReader(const std::string& source) : source_(source) {}
    Json value() {
        Json json;
        const char c = peek();
        if (c == '{') {
            json.kind = Json::Kind::Object;
            ++at_;
            if (peek() == '}') { ++at_; return json; }
            for (;;) {
                const std::string key = string();
                expect(':');
                json.members[key] = value();
                if (peek() == ',') { ++at_; continue; }
                expect('}');
                return json;
            }
        }
        if (c == '[') {
            json.kind = Json::Kind::Array;
            ++at_;
            if (peek() == ']') { ++at_; return json; }
            for (;;) {
                json.items.push_back(value());
                if (peek() == ',') { ++at_; continue; }
                expect(']');
                return json;
            }
        }
        if (c == '"') { json.kind = Json::Kind::String; json.text = string(); return json; }
        if (source_.compare(at_, 4, "true") == 0) { json.kind = Json::Kind::Bool; json.boolean = true; at_ += 4; return json; }
        if (source_.compare(at_, 5, "false") == 0) { json.kind = Json::Kind::Bool; at_ += 5; return json; }
        if (source_.compare(at_, 4, "null") == 0) { at_ += 4; return json; }
        json.kind = Json::Kind::Number;
        while (at_ < source_.size() && (std::isdigit(static_cast<unsigned char>(source_[at_])) || source_[at_] == '-'))
            json.text += source_[at_++];
        if (json.text.empty()) throw std::runtime_error("golden: bad value");
        return json;
    }
};

std::string list(const std::vector<int>& values) {
    std::string text = "[";
    for (size_t i = 0; i < values.size(); ++i) text += (i ? "," : "") + std::to_string(values[i]);
    return text + "]";
}

std::string lists(const std::vector<std::vector<int>>& values) {
    std::string text = "[";
    for (size_t i = 0; i < values.size(); ++i) text += (i ? "," : "") + list(values[i]);
    return text + "]";
}

const char* flag(bool value) { return value ? "true" : "false"; }

/** Whether two JSON values are equal, object members in any order. */
bool sameJson(const Json& a, const Json& b) {
    if (a.kind != b.kind) return false;
    switch (a.kind) {
        case Json::Kind::Null: return true;
        case Json::Kind::Bool: return a.boolean == b.boolean;
        case Json::Kind::Number:
        case Json::Kind::String: return a.text == b.text;
        case Json::Kind::Array:
            if (a.items.size() != b.items.size()) return false;
            for (size_t i = 0; i < a.items.size(); ++i)
                if (!sameJson(a.items[i], b.items[i])) return false;
            return true;
        case Json::Kind::Object:
            if (a.members.size() != b.members.size()) return false;
            for (const auto& [key, value] : a.members) {
                const auto found = b.members.find(key);
                if (found == b.members.end() || !sameJson(value, found->second)) return false;
            }
            return true;
    }
    return false;
}

/**
 * The values the shared golden file records for one graph, under its field names: the pairwise
 * search's pairs and work, the fallback decision at each given limit, then, at the default
 * limit, the pair classes, work bounds, counted work, candidates, basis and flags of the basis
 * search, the replays and counted work with loosened brackets, the completion of a refused
 * basis at half the work before the last pair, and the relevant-cycle families with their
 * counted work.
 */
std::string goldenEntry(const MolGraph& g, const std::vector<uint64_t>& limits) {
    const auto totals = smsd::ring_reference::workTotals(smsd::ring_reference::replayPairs(g));
    std::ostringstream out;
    out << "{\"pairs\":" << totals.pairs << ",\"workBeforeLastPair\":" << totals.beforeLast
        << ",\"lastPairCharge\":" << totals.last << ",\"fallback\":[";
    for (size_t i = 0; i < limits.size(); ++i) {
        RingSearchNotes at;
        (void)smsd::detail_rf::perceiveRingSets(g, false, limits[i], RingSearchOptions(), &at);
        out << (i ? "," : "") << '[' << limits[i] << ',' << flag(at.fellBack) << ']';
    }
    RingSearchNotes notes;
    const auto basis = smsd::detail_rf::perceiveRingSets(g, false, smsd::RING_CANDIDATE_WORK_LIMIT, RingSearchOptions(), &notes);
    out << "],\"classes\":[" << notes.nonTree << ',' << notes.classA << ',' << notes.classB1 << ',' << notes.classB1Plus
        << ',' << notes.classB2 << ',' << notes.simpleRing << "],\"baseWorkLow\":" << notes.baseWorkLow
        << ",\"baseWorkHigh\":" << notes.baseWorkHigh << ",\"replayedPairs\":" << notes.tier2Replays
        << ",\"searchWork\":" << notes.searchWork << ",\"candidates\":" << notes.candidates
        << ",\"basis\":" << lists(basis.sssr) << ",\"legacyRelevant\":" << basis.legacyRcb.size()
        << ",\"truncated\":" << flag(basis.truncated) << ",\"basisTruncated\":" << flag(basis.basisTruncated)
        << ",\"relevantTruncated\":" << flag(basis.relevantTruncated);
    RingSearchOptions loosen;
    loosen.loosenBracketsForTesting = true;
    RingSearchNotes loosened;
    (void)smsd::detail_rf::perceiveRingSets(g, false, smsd::RING_CANDIDATE_WORK_LIMIT, loosen, &loosened);
    out << ",\"loosenedReplayedPairs\":" << loosened.tier2Replays << ",\"loosenedSearchWork\":" << loosened.searchWork;
    RingSearchOptions liftOptions;
    liftOptions.liftRefusedBasis = true;
    RingSearchNotes lift;
    const auto lifted = smsd::detail_rf::perceiveRingSets(g, false, totals.beforeLast / 2, liftOptions, &lift);
    out << ",\"liftedAtHalf\":" << flag(lift.basisLifted) << ",\"liftSearchWork\":" << lift.searchWork
        << ",\"liftBasis\":" << lists(lifted.sssr);
    uint64_t familyWork = 0;
    const auto families = smsd::detail_rf::computeRelevantCycleFamiliesUncached(
        g, smsd::detail_rf::RELEVANT_CYCLE_FAMILY_WORK_LIMIT, &familyWork);
    out << ",\"families\":[";
    for (size_t i = 0; i < families.families.size(); ++i) {
        const auto& f = families.families[i];
        out << (i ? "," : "") << "{\"size\":" << f.size << ",\"root\":" << f.root << ",\"first\":" << f.first
            << ",\"second\":" << f.second << ",\"opposite\":" << f.opposite << ",\"prototype\":" << list(f.prototype)
            << ",\"memberCount\":\"" << f.memberCount << "\"}";
    }
    out << "],\"relevantCycleCount\":\"" << families.relevantCycleCount << "\",\"familiesTruncated\":"
        << flag(families.truncated) << ",\"familyWork\":" << familyWork << "}";
    return out.str();
}

/**
 * The shared golden file, which the Java ring tests write: the key hash of each listed case and
 * every recorded value of each graph, read from its neighbour lists, must be reproduced here.
 */
void checkGolden(const std::string& path) {
    std::ifstream in(path);
    require(static_cast<bool>(in), "golden file " + path);
    std::stringstream buffer;
    buffer << in.rdbuf();
    const std::string text = buffer.str();
    const Json golden = JsonReader(text).value();
    for (const auto& item : golden["hashes"]["cases"].items) {
        std::vector<uint64_t> words;
        for (const auto& word : item["words"].items) words.push_back(std::stoull(word.text));
        require(std::to_string(smsd::detail_rf::ringCandidateKeyHash(words.data(), words.size())) == item["hash"].text,
                "golden: key hash of " + std::to_string(words.size()) + " words");
    }
    size_t checked = 0;
    for (const auto& expected : golden["graphs"].items) {
        Adjacency neighbours;
        for (const auto& row : expected["neighbors"].items) neighbours.push_back(row.integers());
        const MolGraph g = fromNeighbours(neighbours);
        std::vector<uint64_t> limits;
        for (const auto& decision : expected["fallback"].items) {
            require(decision.items.size() == 2 && decision.items[0].text.find('-') == std::string::npos,
                    "golden: a fallback limit");
            limits.push_back(decision.items[0].unsignedValue());
        }
        const std::string entry = goldenEntry(g, limits);
        const Json actual = JsonReader(entry).value();
        const std::string label = expected["name"].text + " " + expected["order"].text;
        for (const auto& [key, value] : expected.members) {
            if (key == "name" || key == "source" || key == "order" || key == "neighbors") continue;
            require(actual.members.count(key) != 0 && sameJson(actual[key], value), label + ": " + key);
        }
        ++checked;
    }
    require(checked > 0, "golden: no graphs");
    std::cout << "GOLDEN\tgraphs=" << checked << '\n';
}

std::vector<std::pair<std::string, MolGraph>> readCorpus(const std::string& path) {
    std::vector<std::pair<std::string, MolGraph>> list;
    std::ifstream in(path);
    std::string line;
    while (std::getline(in, line)) {
        std::istringstream fields(line);
        std::string smiles, id;
        fields >> smiles >> id;
        if (smiles.empty()) continue;
        try {
            auto g = smsd::parseSMILES(smiles);
            if (smsd::detail::cycleRank(g) > 0) list.push_back({id.empty() ? smiles : id, std::move(g)});
        } catch (const std::exception&) {}
    }
    return list;
}

} // namespace

int main(int argc, char** argv) {
    const std::vector<std::string> arguments(argv + 1, argv + argc);
    exhaustive = !arguments.empty() && arguments[0] == "full";
    const auto all = fixtures();

    control("below-child-fixture", [&] {
        const MolGraph g = fromNeighbours(belowChildFixture());
        const auto systems = g.cycleSystems();
        const auto edgeIndex = smsd::detail_rf::buildEdgeIndex(g);
        const auto edgeList = smsd::detail_rf::buildEdgeList(g);
        RingSearchNotes notes;
        RingCandidateSearch search(g, systems, edgeIndex, edgeList, smsd::RING_CANDIDATE_WORK_LIMIT, RingSearchOptions(), notes);
        search.buildTreeForTesting(0);
        const int edge = edgeIndex.at(smsd::detail_rf::edgeKey(0, 2));
        const auto reading = search.readPairForTesting(edge);
        std::vector<int> first, second;
        search.pathsForTesting(reading, first, second);
        require(reading.kind == Kind::ClassB1Plus && reading.coverFrom == 4 && reading.coverTo == 2, "class and cover 4->2");
        require(first == std::vector<int>({0}) && second == std::vector<int>({0, 1, 4, 2}), "paths [0] and [0,1,4,2]");
        require(reading.searchLow == 15 && reading.searchHigh == 15, "work 19 - 4 = 15");
        for (const auto& record : smsd::ring_reference::replayPairs(g))
            if (record.root == 0 && record.edge == edge) require(record.search == 15, "the pairwise search reads 15");
        const auto rings = smsd::detail_rf::perceiveRingSets(g, false);
        bool found = false;
        for (const auto& cycle : rings.legacyRcb) found = found || cycle == std::vector<int>({0, 2, 4, 1});
        require(found, "candidate atoms [0,2,4,1]");
    });

    PairStatistics fixtureStatistics;
    control("pairs-match-the-pairwise-search", [&] {
        for (const auto& fixture : all)
            for (const auto& [variant, graph] : variants(fixture))
                comparePairs(fixture.id + "/" + variant, graph, fixtureStatistics);
    });
    fixtureStatistics.print("fixtures");

    control("results-match-the-pairwise-search", [&] {
        // By default the given and sorted orders, at and around the limit unless the graph is large,
        // with every option set for the given order; the full run reads every order with every set.
        for (const auto& fixture : all)
            for (const auto& [variant, graph] : variants(fixture, exhaustive))
                compareResults(fixture.id + "/" + variant, graph, exhaustive || !fixture.large,
                               !fixture.large && (exhaustive ? variant != "relabel23" : variant == "given"));
    });

    control("graphs-without-ring-systems-use-the-pairwise-search", [&] {
        // A neighbour listed twice leaves no ring systems.
        MolGraph g = fromNeighbours(fourRingChain(2));
        g.neighbors[0].push_back(g.neighbors[0].front());
        g.neighbors[g.neighbors[0].front()].push_back(0);
        for (uint64_t limit : {smsd::RING_CANDIDATE_WORK_LIMIT, uint64_t(40), uint64_t(0)})
            for (bool relevant : {false, true}) {
                RingSearchNotes notes;
                const auto result = smsd::detail_rf::perceiveRingSets(g, relevant, limit, RingSearchOptions(), &notes);
                require(notes.pairwise && sameResult(result, smsd::ring_reference::perceiveRingSets(g, relevant, limit)),
                        "repeated neighbour");
            }
    });

    control("existing-work-hooks", [&] {
        const uint64_t ample = uint64_t(1) << 34;
        for (const char* smiles : {"c1ccc2ccccc2c1", "C1C2CC3CC1CC(C2)C3", "C12C3C4C1C5C2C3C45",
                                   "C12C3C4C5C6C1C7C2C3C4C5C67", "c1cc2ccc3ccc4ccc5ccc6ccc1c7c2c3c4c5c67"}) {
            const auto g = smsd::parseSMILES(smiles);
            for (uint64_t limit : {uint64_t(0), ample})
                require(sameResult(smsd::detail_rf::perceiveRingSets(g, false, limit),
                                   smsd::ring_reference::perceiveRingSets(g, false, limit)), std::string(smiles));
        }
        const auto cd34 = smsd::parseSMILES(cyclodextrin(34));
        const auto rings = smsd::detail_rf::perceiveRingSets(cd34, false, uint64_t(1) << 20);
        require(rings.truncated && !rings.basisTruncated && rings.sssr.size() == 35 && weight(rings.sssr) == 374,
                "CD34 at 2^20 keeps its certified fallback basis of 35 rings and weight 374");
        require(sameResult(rings, smsd::ring_reference::perceiveRingSets(cd34, false, uint64_t(1) << 20)), "CD34 at 2^20");
        require(sameResult(smsd::detail_rf::perceiveRingSets(cd34, true), smsd::ring_reference::perceiveRingSets(cd34, true)),
                "CD34 at the default limit");
    });

    control("lift-refused-basis-changes-only-the-basis", [&] {
        RingSearchOptions lift;
        lift.liftRefusedBasis = true;
        int lifted = 0, refused = 0;
        for (const auto& fixture : all) {
            if (fixture.large) continue;
            const MolGraph& g = fixture.graph;
            const auto totals = smsd::ring_reference::workTotals(smsd::ring_reference::replayPairs(g));
            if (totals.pairs < 2) continue;
            const auto unlimited = smsd::ring_reference::perceiveRingSets(g, false, uint64_t(1) << 40);
            for (uint64_t limit : {totals.beforeLast - 1, totals.beforeLast / 2, totals.beforeLast / 8})
                for (bool relevant : {false, true}) {
                    const auto base = smsd::ring_reference::perceiveRingSets(g, relevant, limit);
                    RingSearchNotes notes;
                    const auto result = smsd::detail_rf::perceiveRingSets(g, relevant, limit, lift, &notes);
                    require(result.rcb == base.rcb && result.legacyRcb == base.legacyRcb && result.truncated == base.truncated
                            && result.relevantTruncated == base.relevantTruncated, fixture.id + ": only the basis changes");
                    if (!notes.basisLifted) {
                        require(sameResult(result, base), fixture.id + ": no lift, no change");
                        refused += base.basisTruncated;
                        continue;
                    }
                    ++lifted;
                    require(base.basisTruncated && !result.basisTruncated && notes.baseBasisTruncated, fixture.id + ": a refusal lifted");
                    require(result.sssr == unlimited.sssr, fixture.id + ": the lifted basis is the main-path basis");
                    require(smsd::detail_rf::certifiedMinimumBasis(g, result.sssr, smsd::detail::cycleRank(g), uint64_t(1) << 40),
                            fixture.id + ": the lifted basis is certified");
                }
        }
        std::cout << "LIFT\tlifted=" << lifted << "\tstill-refused=" << refused << '\n';
        require(lifted > 0, "some refused basis is lifted");
    });

    control("relevant-proof-never-claims-a-short-list", [&] {
        RingSearchOptions prove;
        prove.proveRelevantByFamilies = true;
        int proved = 0, open = 0;
        for (const auto& fixture : all) {
            const auto base = smsd::ring_reference::perceiveRingSets(fixture.graph, true);
            RingSearchNotes notes;
            const auto result = smsd::detail_rf::perceiveRingSets(fixture.graph, true, smsd::RING_CANDIDATE_WORK_LIMIT, prove, &notes);
            require(result.sssr == base.sssr && result.rcb == base.rcb && result.basisTruncated == base.basisTruncated,
                    fixture.id + ": the lists are unchanged");
            if (result.relevantTruncated != base.relevantTruncated) {
                require(notes.relevantProvedByFamilies && base.relevantTruncated && !base.basisTruncated, fixture.id + ": proof");
                const auto families = smsd::detail_rf::computeRelevantCycleFamiliesUncached(fixture.graph);
                std::vector<std::vector<int>> members;
                if (familyMembers(fixture.graph, families, 20000, members))
                    require(keysOf(members) == keysOf(base.rcb), fixture.id + ": a proved list is every relevant cycle");
                std::cout << "PROVED\t" << fixture.id << "\trelevant=" << base.rcb.size() << '\n';
                ++proved;
            } else if (base.relevantTruncated) {
                ++open;
            }
        }
        std::cout << "PROOF\tproved=" << proved << "\tstill-incomplete=" << open << '\n';
    });

    control("families-are-the-relevant-cycles", [&] {
        for (const auto& fixture : all)
            for (const auto& [variant, graph] : variants(fixture, exhaustive))
                if (!fixture.large) compareFamilies(fixture.id + "/" + variant, graph);
    });

    control("family-records", [&] {
        auto familiesOf = [](const std::string& smiles) {
            return smsd::computeRelevantCycleFamilies(smsd::parseSMILES(smiles));
        };
        const auto cubane = familiesOf("C12C3C4C1C5C2C3C45");
        require(cubane.families.size() == 6 && cubane.relevantCycleCount == "6" && !cubane.truncated, "cubane has six faces");
        for (const auto& family : cubane.families) require(family.memberCount == "1" && family.size == 4, "one member each");
        require(familiesOf("C1C2CC3CC1CC(C2)C3").relevantCycleCount == "4", "adamantane");
        require(familiesOf("C1CC2CCC1C2").relevantCycleCount == "2", "norbornane");
        require(familiesOf("c1ccc2ccccc2c1").relevantCycleCount == "2", "naphthalene");
        // C1CCCC1: atoms 0-1-2-3-4; root 4, opposite bond (1, 2). c1ccccc1: root 5, opposite atom 2.
        const auto five = familiesOf("C1CCCC1").families, six = familiesOf("c1ccccc1").families;
        require(five.size() == 1 && five[0].root == 4 && five[0].first == 1 && five[0].second == 2 && five[0].opposite == -1
                && five[0].prototype.front() == 4, "an odd simple ring");
        require(six.size() == 1 && six[0].root == 5 && six[0].first == 1 && six[0].second == 3 && six[0].opposite == 2
                && six[0].prototype.front() == 5, "an even simple ring");
        for (const auto& [units, total] : std::vector<std::pair<int, std::string>>{
                 {6, "70"}, {7, "135"}, {34, "17179869218"}, {40, "1099511627816"},
                 {97, "158456325028528675187087900769"}}) {
            const auto g = smsd::parseSMILES(cyclodextrin(units));
            const auto families = smsd::computeRelevantCycleFamilies(g);
            require(!families.truncated && families.relevantCycleCount == total, "CD" + std::to_string(units) + " total");
            // The glucose rings are one family each; the 2^n macrocycles split by their highest atom.
            size_t glucose = 0, largest = 0;
            for (size_t f = 0; f < families.families.size(); ++f) {
                const auto& family = families.families[f];
                require(family.size == 6 || family.size == 5 * units, "ring sizes");
                if (family.size == 6) { ++glucose; continue; }
                require((family.opposite < 0) == (units % 2 == 1), "the macrocycle's parity");
                const auto& best = families.families[largest];
                if (best.size != family.size || best.memberCount.size() < family.memberCount.size()
                    || (best.memberCount.size() == family.memberCount.size() && best.memberCount < family.memberCount))
                    largest = f;
            }
            require(glucose == static_cast<size_t>(units), "one family per glucose ring");
            const auto& macrocycle = families.families[largest];
            const auto members = smsd::relevantCycleFamilyMembers(g, macrocycle, 100);
            const size_t expected = macrocycle.memberCount.size() > 3 ? 100 : std::min<size_t>(100, std::stoul(macrocycle.memberCount));
            require(members.size() == expected, "a member limit cuts the fixed order");
            require(members.front() == macrocycle.prototype && members == smsd::relevantCycleFamilyMembers(g, macrocycle, 100),
                    "members come in a fixed order");
            for (const auto& member : members)
                require(smsd::relevantCycleFamilyOf(g, families, member) == static_cast<int>(largest), "a macrocycle member");
        }
        const auto cd8 = smsd::parseSMILES(cyclodextrin(8));
        uint64_t work = 0;
        const auto full = smsd::detail_rf::computeRelevantCycleFamiliesUncached(cd8, smsd::detail_rf::RELEVANT_CYCLE_FAMILY_WORK_LIMIT, &work);
        require(!full.truncated && smsd::detail_rf::computeRelevantCycleFamiliesUncached(cd8, work).families.size() == full.families.size(),
                "the counted work is enough");
        for (uint64_t limit : {uint64_t(0), work / 100, work / 4, work / 2, work - 1}) {
            const auto cut = smsd::detail_rf::computeRelevantCycleFamiliesUncached(cd8, limit);
            require(cut.truncated, "a lowered limit truncates");
            for (const auto& family : cut.families) {
                bool listed = false;
                for (const auto& other : full.families)
                    listed = listed || (other.root == family.root && other.first == family.first && other.second == family.second
                                        && other.opposite == family.opposite && other.memberCount == family.memberCount);
                require(listed, "a truncated list keeps only relevant families");
            }
        }
    });

    control("membership-is-exact", [&] {
        // r = 7 and a1, a2, a3, b1, b2, b3, y = 0 .. 6; r-a1-b1-a2-b2-y-b3-a3 is a ring whose arc
        // from r does not advance one level per atom (a2 is next to r), so it is no member.
        Adjacency neighbours(8);
        auto bond = [&](int a, int b) { neighbours[a].push_back(b); neighbours[b].push_back(a); };
        for (const auto& [a, b] : std::vector<std::pair<int, int>>{{7, 0}, {7, 1}, {7, 2}, {0, 3}, {1, 3}, {1, 4},
                                                                  {2, 5}, {3, 6}, {4, 6}, {5, 6}})
            bond(a, b);
        const MolGraph g = fromNeighbours(neighbours);
        const auto families = smsd::computeRelevantCycleFamilies(g);
        require(smsd::relevantCycleFamilyOf(g, families, {7, 0, 3, 1, 4, 6, 5, 2}) == -1, "the eight-ring is not a member");
        std::vector<std::vector<int>> members;
        require(familyMembers(g, families, 1000, members), "small families");
        for (const auto& member : members) require(smsd::relevantCycleFamilyOf(g, families, member) >= 0, "members are members");
        require(smsd::relevantCycleFamilyOf(g, families, {7, 0, 3}) == -1, "not a cycle");
        require(smsd::relevantCycleFamilyOf(g, families, {7, 0, 3, 1, 7}) == -1, "a repeated atom");
    });

    control("family-cache-follows-bond-changes", [&] {
        auto graph = smsd::parseSMILES("CCCCCC");
        require(smsd::computeRelevantCycleFamilies(graph).families.empty(), "a chain has none");
        graph.setBondOrder(0, 5, 1);
        require(smsd::computeRelevantCycleFamilies(graph).relevantCycleCount == "1", "closing the chain makes a ring");
    });

#ifdef SMSD_RING_SEARCH_GOLDEN_JSON
    control("golden", [&] { checkGolden(SMSD_RING_SEARCH_GOLDEN_JSON); });
#endif
    if (arguments.size() == 2 && arguments[0] == "check-golden") control("golden", [&] { checkGolden(arguments[1]); });

    if (!arguments.empty() && arguments[0] == "corpus") {
        for (size_t i = 1; i < arguments.size(); ++i) {
            PairStatistics statistics;
            const auto corpus = readCorpus(arguments[i]);
            control("corpus " + arguments[i], [&] {
                for (const auto& [id, graph] : corpus) {
                    Adjacency sorted = graph.neighbors;
                    for (auto& row : sorted) std::sort(row.begin(), row.end());
                    for (const MolGraph& g : {graph, fromNeighbours(sorted)}) {
                        comparePairs(id, g, statistics);
                        compareResults(id, g, true, false);
                        compareFamilies(id, g);
                    }
                }
            });
            statistics.print(arguments[i]);
        }
    }
    std::cout << "SUMMARY\t" << failures << '\n';
    return failures == 0 ? 0 : 1;
}
