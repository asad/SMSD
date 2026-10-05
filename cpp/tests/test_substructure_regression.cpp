/* SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 */
#include "smsd/vf2pp.hpp"
#include "smsd/smiles_parser.hpp"
#include <functional>
#include <iostream>
#include <set>
#include <stdexcept>
#include <string>

using MappingSet = std::set<std::vector<int>>;

static void require(bool ok, const std::string& message) {
    if (!ok) throw std::runtime_error(message);
}

static smsd::MolGraph graph(int n, unsigned edges) {
    std::vector<std::vector<int>> neighbors(n), orders(n);
    unsigned bit = 0;
    for (int i = 0; i < n; ++i)
        for (int j = i + 1; j < n; ++j, ++bit)
            if (edges & (1u << bit)) {
                neighbors[i].push_back(j); neighbors[j].push_back(i);
                orders[i].push_back(1); orders[j].push_back(1);
            }
    return smsd::MolGraph::Builder().atomCount(n)
        .atomicNumbers(std::vector<int>(n, 6))
        .setNeighbors(neighbors).setBondOrders(orders).build(false);
}

// Independent injective-map oracle: no search-engine compatibility/pruning calls.
static MappingSet oracle(const smsd::MolGraph& q, const smsd::MolGraph& t, bool induced) {
    MappingSet out;
    std::vector<int> map(q.n, -1);
    std::vector<bool> used(t.n, false);
    std::function<void(int)> visit = [&](int i) {
        if (i == q.n) { out.insert(map); return; }
        for (int j = 0; j < t.n; ++j) {
            if (used[j]) continue;
            bool ok = true;
            for (int k = 0; k < i && ok; ++k) {
                const bool qe = q.hasBond(i, k), te = t.hasBond(j, map[k]);
                ok = (!qe || te) && (!induced || qe == te);
            }
            if (!ok) continue;
            map[i] = j; used[j] = true;
            visit(i + 1);
            used[j] = false; map[i] = -1;
        }
    };
    visit(0);
    return out;
}

static MappingSet mappings(const smsd::MolGraph& q, const smsd::MolGraph& t,
                           const smsd::ChemOptions& c) {
    const auto all = smsd::findAllSubstructures(q, t, c);
    MappingSet out;
    for (const auto& pairs : all) {
        std::vector<int> map(q.n, -1);
        for (const auto& [qi, ti] : pairs) {
            require(qi >= 0 && qi < q.n && ti >= 0 && ti < t.n, "mapping index bounds");
            require(map[qi] == -1, "duplicate query index");
            map[qi] = ti;
        }
        require(pairs.size() == static_cast<size_t>(q.n), "incomplete mapping");
        out.insert(map);
    }
    require(out.size() == all.size(), "duplicate mappings");
    return out;
}

static void selfEnumeration() {
    smsd::ChemOptions c;
    for (const auto engine : {smsd::ChemOptions::MatcherEngine::VF2,
                              smsd::ChemOptions::MatcherEngine::VF2PP}) {
        c.matcherEngine = engine;
        auto chain = smsd::parseSMILES("CCC");
        auto copy = chain;
        require(mappings(chain, chain, c).size() == 2, "self chain needs both orientations");
        require(mappings(chain, copy, c).size() == 2, "equal copy needs both orientations");
        auto benzene = smsd::parseSMILES("c1ccccc1");
        auto reordered = smsd::parseSMILES("c1ccccc1");
        require(mappings(benzene, reordered, c).size() == 12, "benzene needs all 12 atom mappings");
    }
}

static void disconnectedPath() {
    auto mixed = smsd::parseSMILES("C1CC1.CC");
    auto copy = smsd::parseSMILES("C1CC1.CC");
    require(!smsd::detail::isSimplePathGraph(mixed), "cycle plus edge is not a path");
    smsd::ChemOptions c;
    require(smsd::isSubstructure(mixed, copy, c), "disconnected self match");
    require(smsd::findSubstructure(mixed, copy, c).size() == 5, "complete disconnected mapping");
    require(mappings(mixed, copy, c) == oracle(mixed, copy, false), "disconnected automorphisms");
}

static smsd::MolGraph aromaticPath(int n, bool aromaticBonds) {
    std::vector<std::vector<int>> neighbors(n), orders(n);
    std::vector<std::vector<bool>> arom(n);
    for (int i = 1; i < n; ++i) {
        neighbors[i - 1].push_back(i); neighbors[i].push_back(i - 1);
        orders[i - 1].push_back(1); orders[i].push_back(1);
        arom[i - 1].push_back(aromaticBonds); arom[i].push_back(aromaticBonds);
    }
    return smsd::MolGraph::Builder().atomCount(n)
        .atomicNumbers(std::vector<int>(n, 6)).aromaticFlags(std::vector<uint8_t>(n, 1))
        .setNeighbors(neighbors).setBondOrders(orders).bondAromaticFlags(arom).build(false);
}

static void strictAromaticBonds() {
    auto q = aromaticPath(3, true), t = aromaticPath(500, false);
    smsd::ChemOptions c;
    c.matchBondOrder = smsd::ChemOptions::BondOrderMode::ANY;
    c.aromaticityMode = smsd::ChemOptions::AromaticityMode::STRICT;
    for (const auto engine : {smsd::ChemOptions::MatcherEngine::VF2,
                              smsd::ChemOptions::MatcherEngine::VF2PP}) {
        c.matcherEngine = engine;
        require(!smsd::isSubstructure(q, t, c), "ANY bond order must respect strict aromatic bonds");
        require(smsd::findSubstructure(q, t, c).empty(), "strict aromaticity mapping consistency");
    }
}

static void stickyDeadline() {
    smsd::detail::TimeBudget budget(10000);
    budget.deadline = smsd::detail::TimeBudget::Clock::now() - std::chrono::milliseconds(1);
    require(budget.expired(), "the first budget check must detect an expired deadline");
    for (int i = 0; i < 1024; ++i)
        require(budget.expired(), "expired throttled checks must never resume search");
    smsd::detail::TimeBudget direct(10000);
    direct.deadline = budget.deadline;
    require(direct.expiredNow(), "direct budget check detects expiry");
    require(direct.expired(), "direct expiry must latch throttled checks");
}

static void exhaustiveGraphs() {
    int cases = 0;
    for (int qn = 0; qn <= 4; ++qn)
        for (int tn = 0; tn <= 4; ++tn)
            for (unsigned qm = 0; qm < (1u << (qn * (qn - 1) / 2)); ++qm)
                for (unsigned tm = 0; tm < (1u << (tn * (tn - 1) / 2)); ++tm) {
                    auto q = graph(qn, qm), t = graph(tn, tm);
                    for (bool induced : {false, true}) {
                        const auto expected = oracle(q, t, induced);
                        for (const auto engine : {smsd::ChemOptions::MatcherEngine::VF2,
                                                  smsd::ChemOptions::MatcherEngine::VF2PP}) {
                            smsd::ChemOptions c;
                            c.induced = induced; c.matcherEngine = engine;
                            const auto context = "oracle q=" + std::to_string(qn) + "/" + std::to_string(qm)
                                + " t=" + std::to_string(tn) + "/" + std::to_string(tm);
                            require(mappings(q, t, c) == expected, context);
                            require(smsd::isSubstructure(q, t, c) == !expected.empty(), context + " exists");
                            ++cases;
                        }
                    }
                }
    std::cout << cases << " exhaustive topology/engine/policy cases passed\n";
}

int main(int argc, char** argv) {
    try {
        const std::string group = argc > 1 ? argv[1] : "all";
        if (group == "all" || group == "self") selfEnumeration();
        if (group == "all" || group == "path") disconnectedPath();
        if (group == "all" || group == "aromatic") strictAromaticBonds();
        if (group == "all" || group == "timeout") stickyDeadline();
        if (group == "all" || group == "oracle") exhaustiveGraphs();
        std::cout << "Substructure regression checks passed\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
