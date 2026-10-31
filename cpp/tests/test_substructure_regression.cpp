/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
#include "test_support.hpp"
#include "timing_allowance.hpp"
#include "smsd/substructure.hpp"
#include "smsd/smiles_parser.hpp"
#include <chrono>
#include <functional>
#include <iostream>
#include <map>
#include <numeric>
#include <random>
#include <set>
#include <stdexcept>
#include <string>

using MappingSet = std::set<std::vector<int>>;

using smsd_test::require;

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
                           const smsd::ChemOptions& c, bool induced = false) {
    auto options = smsd_test::substructureOptions(c);
    options.induced = induced;
    const auto result = smsd::Substructure::find(q).withOptions(options).matchAll(t);
    std::vector<std::vector<std::pair<int, int>>> all;
    for (const auto& match : result.matches) all.emplace_back(match.atomMapping.begin(), match.atomMapping.end());
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
    {

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
    require(mixed.neighbors[0].size() == 2 && mixed.neighbors[3].size() == 1, "cycle plus edge fixture");
    smsd::ChemOptions c;
    require(smsd_test::substructureExists(mixed, copy, c), "disconnected self match");
    require(smsd_test::substructureAtoms(mixed, copy, c).size() == 5, "complete disconnected mapping");
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
    {

        require(!smsd_test::substructureExists(q, t, c), "ANY bond order must respect strict aromatic bonds");
        require(smsd_test::substructureAtoms(q, t, c).empty(), "strict aromaticity mapping consistency");
    }
}

static void stickyDeadline() {
    smsd::SearchOptions options;
    options.timeoutMs = 0;
    smsd::search_detail::Budget budget(options, 1, 1, false);
    require(budget.stopped(), "the first budget check must detect an expired deadline");
    for (int i = 0; i < 1024; ++i)
        require(budget.stopped(), "expired checks must never resume search");
    require(budget.status() == smsd::SearchStatus::TimedOut, "expiry records its reason");
}

// Element and degree counts that rule out an embedding end the search before it starts, and
// equivalent target atoms still leave every embedding reachable.
static void countScreens() {
    const auto options = smsd_test::substructureOptions(smsd::ChemOptions{});
    for (const auto& [query, target] : std::vector<std::pair<std::string, std::string>>{
             {"OC(=O)CC(=O)O", "OC(=O)CCCCCC"}, {"CC(C)(C)C", "CCCCCCCC"}, {"NCC(N)CN", "NCCCC(N)CCCCO"}}) {
        const auto result = smsd::Substructure::find(query).withOptions(options).match(target);
        require(result.matches.empty() && result.status == smsd::SearchStatus::Completed
            && result.complete && result.optimal && result.nodesVisited == 0, "count screen: " + query + " in " + target);
    }
    for (const auto& [query, target] : std::vector<std::pair<std::string, std::string>>{
             {"C.C", "C.C.C"}, {"CC(C)C", "CC(C)(C)C"}, {"FC(F)F", "FC(F)(F)C(F)(F)F"}, {"OC(=O)CC(=O)O", "OC(=O)CC(C(=O)O)C(=O)O"}}) {
        const auto q = smsd::parseSMILES(query), t = smsd::parseSMILES(target);
        const auto first = smsd::Substructure::find(q).withOptions(options).match(t);
        require(!first.matches.empty() && first.complete, "equivalent targets keep a match: " + query + " in " + target);
        std::vector<int> map(q.n, -1);
        for (const auto& [atom, image] : first.matches.front().atomMapping) map[atom] = image;
        require(mappings(q, t, smsd::ChemOptions{}).count(map) == 1, "first match is an embedding: " + query);
    }
}

// Target symmetry classes cost little for many isolated atoms, and a hub with many neighbours,
// which still makes them dear, cannot hold the search past its deadline.
static void twinScans() {
    constexpr int atoms = 16000;
#if !defined(SMSD_TEST_SANITISED) && (defined(__OPTIMIZE__) || (defined(_MSC_VER) && !defined(_DEBUG)))
    constexpr int64_t timeout = 30, margin = 120;
#else
    // An unoptimised or sanitised build scans the isolated atoms about twenty times more slowly; a quadratic
    // scan of 16,000 atoms would still take several seconds, so the bound keeps its meaning.
    constexpr int64_t timeout = 1000, margin = 1000;
#endif
    {
        std::vector<int> classes(atoms);
        std::iota(classes.begin(), classes.end(), 1);
        const auto waters = smsd::MolGraph::Builder().atomCount(atoms).atomicNumbers(std::vector<int>(atoms, 8))
            .hydrogenCounts(std::vector<int>(atoms, 2)).atomClasses(classes)
            .setNeighbors(std::vector<std::vector<int>>(atoms)).build(false);
        const auto result = smsd::Substructure::find("O").withTimeout(timeout).match(waters);
        require(result.status == smsd::SearchStatus::Completed && result.matches.size() == 1
            && result.elapsedMillis < timeout + margin, "isolated atoms must not make target classes quadratic");
    }
    {
        std::vector<int> elements(atoms + 1, 8), classes(atoms + 1);
        std::iota(classes.begin(), classes.end(), 0);
        elements[0] = 6;
        std::vector<std::vector<int>> neighbors(atoms + 1), orders(atoms + 1);
        for (int leaf = 1; leaf <= atoms; ++leaf) {
            neighbors[0].push_back(leaf); orders[0].push_back(1);
            neighbors[leaf].push_back(0); orders[leaf].push_back(1);
        }
        const auto hub = smsd::MolGraph::Builder().atomCount(atoms + 1).atomicNumbers(elements).atomClasses(classes)
            .setNeighbors(neighbors).setBondOrders(orders).build(false);
        const auto result = smsd::Substructure::find("CO").withTimeout(timeout).match(hub);
        require(result.elapsedMillis < timeout + margin
            && (result.status == smsd::SearchStatus::TimedOut || !result.matches.empty()),
            "target classes must stop at the deadline");
    }
    // Equal isolated atoms share a class, so the first match still takes the lowest atoms.
    const auto first = smsd::Substructure::find("O.O").match("O.O.C.O.O");
    require(first.matches.size() == 1 && first.matches.front().atomMapping == std::map<int, int>{{0, 0}, {1, 1}},
        "isolated twins keep the first embedding");
}

// atomCompatible reads exactly the properties in atomCompatibilityTerms, and only the element
// under the options for which atomCompatibleByElementOnly holds.
static void compatibilityTermsParity() {
    std::mt19937 random(20261019);
    auto pick = [&](std::initializer_list<int> values) { return *(values.begin() + random() % values.size()); };
    // Atoms copy one of a few term templates; the properties outside the terms vary freely.
    auto molecule = [&](int atoms) {
        std::vector<std::array<int, 8>> templates(3);
        for (auto& terms : templates)
            terms = {pick({6, 7}), pick({0, 1}), pick({0, 1}), pick({0, 1}), pick({0, 12, 13}), pick({0, 1, 2}),
                     pick({2, 3, 4}), pick({0, 1, 2})};
        std::vector<int> elements(atoms), charges(atoms), masses(atoms), hydrogens(atoms), classes(atoms), chirality(atoms);
        std::vector<uint8_t> aromatic(atoms), ring(atoms);
        std::vector<std::array<int, 8>> chosen(atoms);
        for (int atom = 0; atom < atoms; ++atom) {
            const auto& terms = chosen[atom] = templates[random() % templates.size()];
            elements[atom] = terms[0]; charges[atom] = terms[1]; aromatic[atom] = static_cast<uint8_t>(terms[2]);
            ring[atom] = static_cast<uint8_t>(terms[3]); masses[atom] = terms[4]; chirality[atom] = terms[5];
            hydrogens[atom] = pick({0, 1, 2}); classes[atom] = pick({0, 1, 2});
        }
        auto graph = smsd::MolGraph::Builder().atomCount(atoms).atomicNumbers(elements).formalCharges(charges)
            .massNumbers(masses).hydrogenCounts(hydrogens).atomClasses(classes).aromaticFlags(aromatic).ringFlags(ring)
            .tetrahedralChirality(chirality).setNeighbors(std::vector<std::vector<int>>(atoms)).build(false);
        for (int atom = 0; atom < atoms; ++atom) {
            graph.atomicNum[atom] = elements[atom]; graph.formalCharge[atom] = charges[atom];
            graph.aromatic[atom] = aromatic[atom]; graph.ring[atom] = ring[atom];
            graph.massNumber[atom] = masses[atom]; graph.tetraChirality[atom] = chirality[atom];
            graph.degree[atom] = chosen[atom][6];
            graph.tautomerClass[atom] = pick({-1, 0, 1});
        }
        graph.ringCount.assign(atoms, 0);
        for (int atom = 0; atom < atoms; ++atom) graph.ringCount[atom] = chosen[atom][7];
        graph.ringCountsComputed_ = true;
        return graph;
    };
    int combinations = 0, elementOnly = 0;
    for (int trial = 0; trial < 6; ++trial) {
        const auto query = molecule(14), target = molecule(14);
        for (int bits = 0; bits < 3 * 128; ++bits) {
            smsd::ChemOptions chemistry;
            chemistry.matchAtomType = bits & 1;
            chemistry.matchFormalCharge = bits & 2;
            chemistry.aromaticityMode = bits & 4 ? smsd::ChemOptions::AromaticityMode::STRICT : smsd::ChemOptions::AromaticityMode::FLEXIBLE;
            chemistry.ringMatchesRingOnly = bits & 8;
            chemistry.matchIsotope = bits & 16;
            chemistry.useChirality = bits & 32;
            chemistry.tautomerAware = bits & 64;
            chemistry.ringFusionMode = static_cast<smsd::ChemOptions::RingFusionMode>(bits / 128);
            const bool byElement = smsd::search_detail::atomCompatibleByElementOnly(chemistry);
            elementOnly += byElement;
            ++combinations;
            for (int first = 0; first < query.n; ++first)
                for (int second = 0; second < target.n; ++second) {
                    const bool compatible = smsd::search_detail::atomCompatible(query, first, target, second, chemistry);
                    if (byElement) require(compatible == (query.atomicNum[first] == target.atomicNum[second]),
                        "element-only options must compare elements alone");
                    // Atoms with equal terms on either side must agree with every atom of the other.
                    for (int other = first + 1; other < query.n; ++other)
                        if (smsd::search_detail::atomCompatibilityTerms(query, first) == smsd::search_detail::atomCompatibilityTerms(query, other))
                            require(compatible == smsd::search_detail::atomCompatible(query, other, target, second, chemistry),
                                "query atoms with equal terms must match the same target atoms");
                    for (int other = second + 1; other < target.n; ++other)
                        if (smsd::search_detail::atomCompatibilityTerms(target, second) == smsd::search_detail::atomCompatibilityTerms(target, other))
                            require(compatible == smsd::search_detail::atomCompatible(query, first, target, other, chemistry),
                                "target atoms with equal terms must match the same query atoms");
                }
        }
    }
    require(combinations == 6 * 384 && elementOnly > 0, "every option combination must be checked");
}

// Recursive SMARTS answers do not depend on how much of the per-target cache the allowance covers.
static void recursiveAllowance() {
    const std::vector<std::string> patterns = {
        "[C;$(CO)]", "[$(C=O)]~[$(C=O)]", "[C;$(C(=O)O)]~[$([#8;$(O[H]),$(O=*)])]", "[!$(C=O);#6]~[#8]",
        "[$(*~[$(N~C)])]", "[$(c1ccccc1),$(C(=O)[OH])]~*", "[$([CH2]~[CH2]~[CH2])]"};
    const std::vector<std::string> targets = {"CC(=O)O", "OC(=O)CC(=O)O", "NCC(=O)Nc1ccccc1", "CCCCCCO", "c1ccc2ccccc2c1C(=O)O"};
    auto answers = [&](size_t allowance) {
        smsd::search_detail::RecursiveMemo::threadAllowance() = allowance;
        std::vector<std::pair<std::set<std::map<int, int>>, std::map<int, int>>> found;
        for (const auto& pattern : patterns) {
            const auto query = smsd::Substructure::smarts(pattern);
            for (const auto& smiles : targets) {
                const auto target = smsd::parseSMILES(smiles);
                std::set<std::map<int, int>> all;
                for (const auto& match : query.matchAll(target).matches) all.insert(match.atomMapping);
                const auto first = query.match(target);
                found.emplace_back(std::move(all), first.matches.empty() ? std::map<int, int>{} : first.matches.front().atomMapping);
            }
        }
        smsd::search_detail::RecursiveMemo::threadAllowance() = smsd::search_detail::RecursiveMemo::allowanceBytes;
        return found;
    };
    const auto cached = answers(smsd::search_detail::RecursiveMemo::allowanceBytes);
    for (size_t allowance : {size_t(0), size_t(12), size_t(40), size_t(100)})
        require(answers(allowance) == cached, "uncached recursive SMARTS answers differ, allowance " + std::to_string(allowance));
}

// Plans and recursive caches find each of many subpatterns directly. The input limits allow
// 5000 one-atom subpatterns, half of MAX_ATOMS; more are rejected.
static void manySubpatterns() {
    std::string chain;
    for (int atom = 0; atom < 5000; ++atom) chain += "[$(C)]";
    const auto pattern = smsd::parseSMARTS(chain);
    const auto started = std::chrono::steady_clock::now();
    const auto prepared = smsd::Substructure::smarts(pattern);
    const auto planned = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - started).count();
    require(planned < 400, "planning 5000 subpatterns took " + std::to_string(planned) + " ms");
    bool rejected = false;
    try { smsd::parseSMARTS(chain + "[$(C)]"); } catch (const std::invalid_argument&) { rejected = true; }
    require(rejected, "a pattern beyond the atom limit must be rejected");
    std::string atom = "[";
    for (int index = 0; index < 2400; ++index) atom += (index ? ";$(" : "$(") + std::string(1 + index % 7, 'C') + ")";
    atom += "]";
    const auto result = smsd::Substructure::smarts(atom).withTimeout(5000).match("CCCCCCCCCC");
    require(result.status == smsd::SearchStatus::Completed && result.matches.size() == 1,
        "an atom with many subpatterns must finish well within its timeout");
}


// An induced search compares each candidate with the mapped atoms only, and reads the clock
// while it lists them, so a large query of isolated atoms keeps its time limit.
static void sparseInducedTimeout() {
    const int atoms = 32000;
    const auto graph = smsd::MolGraph::Builder().atomCount(atoms).atomicNumbers(std::vector<int>(atoms, 8))
        .setNeighbors(std::vector<std::vector<int>>(atoms)).setBondOrders(std::vector<std::vector<int>>(atoms)).build(false);
    auto options = smsd_test::substructureOptions(smsd::ChemOptions{}, 100);
    options.induced = true;
    const auto result = smsd::Substructure::find(graph).withOptions(options).match(graph);
    require(result.status == smsd::SearchStatus::TimedOut && !result.complete,
        "an induced search of 32000 isolated atoms must stop at its time limit");
    // Comparing each candidate with every query atom would take about 500 ms here.
    require(result.elapsedMillis < smsd_test::timingAllowance(350), "the induced search took " + std::to_string(result.elapsedMillis) + " ms for a 100 ms limit");
}

// An induced SMARTS search compares each candidate with the fewer of the mapped atoms and its
// own neighbours, and reads the clock while it lists candidates.
static void sparseInducedSmartsTimeout() {
    auto isolated = [](int atoms) {
        return smsd::MolGraph::Builder().atomCount(atoms).atomicNumbers(std::vector<int>(atoms, 8))
            .setNeighbors(std::vector<std::vector<int>>(atoms)).setBondOrders(std::vector<std::vector<int>>(atoms)).build(false);
    };
    auto star = [](int leaves) {
        std::vector<std::vector<int>> neighbors(leaves + 1), orders(leaves + 1);
        for (int leaf = 1; leaf <= leaves; ++leaf) {
            neighbors[0].push_back(leaf); orders[0].push_back(1);
            neighbors[leaf].push_back(0); orders[leaf].push_back(1);
        }
        return smsd::MolGraph::Builder().atomCount(leaves + 1).atomicNumbers(std::vector<int>(leaves + 1, 6))
            .setNeighbors(neighbors).setBondOrders(orders).build(false);
    };
    auto starSmarts = [](int leaves) {
        std::string smarts = "C";
        for (int leaf = 1; leaf < leaves; ++leaf) smarts += "(C)";
        return smarts + "C";
    };
    smsd::SubstructureOptions options;
    options.induced = true;
    {
        // The largest pattern MAX_ATOMS allows.
        const int atoms = smsd::MAX_ATOMS;
        std::string smarts;
        for (int atom = 0; atom < atoms; ++atom) smarts += atom ? ".[#8]" : "[#8]";
        options.timeoutMs = 100;
        const auto result = smsd::Substructure::smarts(smarts).withOptions(options).match(isolated(atoms + 6000));
        require(result.status == smsd::SearchStatus::TimedOut && !result.complete,
            "an induced SMARTS search of 10000 isolated atoms must stop at its time limit");
        // Comparing each candidate with every pattern atom would take about 650 ms here.
        require(result.elapsedMillis < smsd_test::timingAllowance(350), "the induced SMARTS search took " + std::to_string(result.elapsedMillis)
            + " ms for a 100 ms limit");
    }
    {
        options.timeoutMs = 30;
        const auto result = smsd::Substructure::smarts(starSmarts(9999)).withOptions(options).match(star(9999));
        require(result.elapsedMillis < smsd_test::timingAllowance(150) && (result.status == smsd::SearchStatus::TimedOut || !result.matches.empty()),
            "an induced star SMARTS took " + std::to_string(result.elapsedMillis) + " ms for a 30 ms limit");
    }
    {
        // Leaves read only their own neighbour, so a large star embeds well within its limit.
        options.timeoutMs = smsd_test::timingAllowance(10000);
        const auto result = smsd::Substructure::smarts(starSmarts(3000)).withOptions(options).match(star(3000));
        require(result.status == smsd::SearchStatus::Completed && result.matches.size() == 1
            && result.matches.front().atomMapping.size() == 3001 && result.matches.front().bondMapping.size() == 3000,
            "an induced star SMARTS must embed in an equal star");
    }
}

// Induced and plain SMARTS embeddings of small graphs equal the independent oracle; the pattern
// keeps the graph's atom order because its bonds are ring closures between dot-separated atoms.
static void smartsOracleGraphs() {
    auto pattern = [](const smsd::MolGraph& graph) {
        std::vector<std::string> closures(graph.n);
        int label = 10;
        for (int atom = 0; atom < graph.n; ++atom)
            for (int neighbor : graph.neighbors[atom])
                if (atom < neighbor) {
                    closures[atom] += "%" + std::to_string(label);
                    closures[neighbor] += "%" + std::to_string(label);
                    ++label;
                }
        std::string smarts;
        for (int atom = 0; atom < graph.n; ++atom) smarts += (atom ? "." : "") + std::string("[#6]") + closures[atom];
        return smarts;
    };
    int cases = 0;
    for (int qn = 1; qn <= 4; ++qn)
        for (int tn = qn; tn <= 5; ++tn)
            for (unsigned qm = 0; qm < (1u << (qn * (qn - 1) / 2)); ++qm)
                for (unsigned tm = 0; tm < (1u << (tn * (tn - 1) / 2)); tm += tn == 5 ? 7 : 1) {
                    const auto q = graph(qn, qm), t = graph(tn, tm);
                    const auto query = smsd::Substructure::smarts(pattern(q));
                    for (bool induced : {false, true}) {
                        smsd::SubstructureOptions options;
                        options.induced = induced;
                        MappingSet found;
                        for (const auto& match : query.withOptions(options).matchAll(t, 10000).matches) {
                            std::vector<int> map(q.n, -1);
                            for (const auto& [atom, image] : match.atomMapping) map[atom] = image;
                            found.insert(map);
                        }
                        require(found == oracle(q, t, induced), "SMARTS oracle q=" + std::to_string(qn) + "/" + std::to_string(qm)
                            + " t=" + std::to_string(tn) + "/" + std::to_string(tm) + (induced ? " induced" : ""));
                        ++cases;
                    }
                }
    std::cout << cases << " SMARTS topology cases passed\n";
}

// Recursive subpatterns of one shape share one layout, and a layout without a bond table finds
// the same bonds as one with it.
static void sharedSubpatternLayouts() {
    const auto pattern = smsd::parseSMARTS("[$(CCO)]~[$(NCS)]~[$(CCO)]~[$(C(C)O)]~[$([C;$(CCO)]=O)]");
    const smsd::search_detail::PatternPlan plan(pattern);
    std::vector<const smsd::SmartsQuery*> nested;
    std::function<void(const smsd::smarts_detail::AtomExpr&)> collect = [&](const smsd::smarts_detail::AtomExpr& expression) {
        if (expression.cachedRecursive) {
            nested.push_back(expression.cachedRecursive.get());
            for (const auto& atom : expression.cachedRecursive->atoms) collect(atom.expr);
        }
        for (const auto& child : expression.children) collect(child);
    };
    for (const auto& atom : pattern.atoms) collect(atom.expr);
    require(nested.size() == 6 && plan.subpatternCount() == 6, "every recursive subpattern must be compiled");
    auto layout = [&](size_t index) { return plan.subpattern(*nested[index]); };
    for (size_t index = 0; index < nested.size(); ++index) require(layout(index) != nullptr, "a subpattern lacks its layout");
    // CCO, NCS, the second CCO and the nested CCO are chains of three; C(C)O branches at its first atom.
    require(layout(0) == layout(1) && layout(0) == layout(2) && layout(0) != layout(3),
        "subpatterns of one shape must share a layout");
    for (const auto* subpattern : nested) {
        const smsd::search_detail::PatternBonds table(*subpattern), map(*subpattern, false);
        for (int first = 0; first < subpattern->atomCount(); ++first)
            for (int second = 0; second < subpattern->atomCount(); ++second)
                require(table.find(first, second) == map.find(first, second), "bond lookups must agree");
    }
    const auto query = smsd::Substructure::smarts(pattern);
    for (const char* smiles : {"OCCNCSCCOC(C)OC(CO)=O", "OCC(NCS)C(CO)C(C)OC(=O)CCO", "CCO"}) {
        const auto target = smsd::parseSMILES(smiles);
        std::set<std::map<int, int>> all;
        for (const auto& match : query.matchAll(target).matches) all.insert(match.atomMapping);
        auto transient = pattern;
        std::function<void(smsd::smarts_detail::AtomExpr&)> clear = [&](smsd::smarts_detail::AtomExpr& expression) {
            expression.cachedRecursive.reset();
            for (auto& child : expression.children) clear(child);
        };
        for (auto& atom : transient.atoms) clear(atom.expr);
        std::set<std::map<int, int>> uncompiled;
        for (const auto& match : smsd::Substructure::smarts(transient).matchAll(target).matches) uncompiled.insert(match.atomMapping);
        require(all == uncompiled, std::string("shared layouts change the matches of ") + smiles);
    }
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
                        {
                            smsd::ChemOptions c;
                            const auto context = "oracle q=" + std::to_string(qn) + "/" + std::to_string(qm)
                                + " t=" + std::to_string(tn) + "/" + std::to_string(tm);
                            require(mappings(q, t, c, induced) == expected, context);
                            auto options = smsd_test::substructureOptions(c);
                            options.induced = induced;
                            require(!smsd::Substructure::find(q).withOptions(options).match(t).matches.empty() == !expected.empty(), context + " exists");
                            ++cases;
                        }
                    }
                }
    std::cout << cases << " exhaustive topology/policy cases passed\n";
}

// Every field but the elapsed time: matches with their bonds and fragments, status, proof flags, score and nodes.
static void requireSameResult(const smsd::SearchResult& expected, const smsd::SearchResult& actual, const std::string& context) {
    require(expected.status == actual.status && expected.optimal == actual.optimal && expected.complete == actual.complete
            && expected.score == actual.score && expected.nodesVisited == actual.nodesVisited, context + ": summary differs");
    require(expected.matches.size() == actual.matches.size(), context + ": match count differs");
    for (size_t index = 0; index < expected.matches.size(); ++index) {
        const auto& first = expected.matches[index];
        const auto& second = actual.matches[index];
        require(first.atomMapping == second.atomMapping && first.bondMapping == second.bondMapping
                && first.fragments == second.fragments, context + ": match " + std::to_string(index) + " differs");
    }
}

// A molecule substructure search as substructure() runs it after its screens, keeping complete
// mappings as rows or, with `rows` false, as matches in the collector throughout.
static smsd::SearchResult vertexSubstructure(const smsd::MolGraph& query, const smsd::MolGraph& target,
                                             const smsd::SubstructureOptions& options, size_t limit, bool all, bool rows) {
    using namespace smsd::search_detail;
    const QueryView view(query);
    Budget budget(options, query.n, target.n, false);
    RecursiveScope scope(budget);
    Collector collector(nullptr, limit, all);
    if (budget.stopped()) return collector.finish(budget, false, false);
    prepare(&query, target, options.chemistry);
    VertexSearch search(view, target, options.chemistry, nullptr, options.induced, all, budget, collector,
                        -1, nullptr, true, 0, rows);
    const auto proof = search.run();
    return collector.finish(budget, proof.first, proof.second);
}

/** Both ways of keeping complete mappings for every result mode; the number of searches compared. */
static int compareRows(const smsd::MolGraph& query, const smsd::MolGraph& target,
                       const smsd::SubstructureOptions& options, const std::string& context) {
    int compared = 0;
    for (size_t limit : {size_t(1), size_t(2), size_t(7), size_t(1000)})
        for (bool all : {false, true}) {
            if (!all && limit != 1) continue;
            const auto kept = vertexSubstructure(query, target, options, limit, all, false);
            const auto rows = vertexSubstructure(query, target, options, limit, all, true);
            requireSameResult(kept, rows, context + (all ? " all " : " first ") + std::to_string(limit));
            ++compared;
        }
    return compared;
}

// Rows give the results of matches kept in the collector throughout, under every chemistry option.
static void completeMappingRows() {
    const std::vector<std::string> queries = {
        "", "C", "c1ccccc1", "c1ccc2ccccc2c1", "CC(C)(C)C", "C.C", "CC(=O)O", "C[C@H](N)C(=O)O", "C/C=C/C",
        "C1CCCCC1", "CC(=O)N", "c1ccncc1", "OC=C", "Oc1ccccc1", "CC(C)C"};
    const std::vector<std::string> targets = {
        "CC(C)Cc1ccc(cc1)[C@@H](C)C(=O)O", "CC(=O)Oc1ccccc1C(=O)O", "c1ccc2cc3ccccc3cc2c1", "CC(C)(C)c1ccc(O)cc1",
        "C/C=C/C=C/C(=O)O", "C/C=C\\C", "N[C@@H](C)C(=O)O", "N[C@H](C)C(=O)O", "CC(C)(C)C(C)(C)C", "C1CCC2CCCCC2C1",
        "O=C1NC(=O)c2ccccc21", "OC1=CC=CC=C1", "CC.CC.C", "C1CC1C1CC1", "c1ccc(cc1)-c1ccccc1",
        "CN1C=NC2=C1C(=O)N(C(=O)N2C)C", "C12CC3CC(CC(C3)C1)C2", "CC(O)=CC(C)=O", "OC(=O)CC(N)C(=O)O"};
    std::vector<std::pair<std::string, smsd::ChemOptions>> chemistries;
    chemistries.emplace_back("default", smsd::ChemOptions{});
    auto add = [&](const std::string& name, const std::function<void(smsd::ChemOptions&)>& change) {
        smsd::ChemOptions chemistry;
        change(chemistry);
        chemistries.emplace_back(name, chemistry);
    };
    add("chirality", [](smsd::ChemOptions& c) { c.useChirality = true; });
    add("bond-stereo", [](smsd::ChemOptions& c) { c.useBondStereo = true; });
    add("tautomers", [](smsd::ChemOptions& c) { c.tautomerAware = true; });
    add("complete-rings", [](smsd::ChemOptions& c) { c.completeRingsOnly = true; });
    add("ring-only", [](smsd::ChemOptions& c) { c.ringMatchesRingOnly = true; });
    add("loose", [](smsd::ChemOptions& c) { c.matchBondOrder = smsd::ChemOptions::BondOrderMode::LOOSE; });
    add("any-order", [](smsd::ChemOptions& c) { c.matchBondOrder = smsd::ChemOptions::BondOrderMode::ANY; });
    add("strict-aromatic", [](smsd::ChemOptions& c) { c.aromaticityMode = smsd::ChemOptions::AromaticityMode::STRICT; });
    add("any-element", [](smsd::ChemOptions& c) { c.matchAtomType = false; });
    add("ring-fusion", [](smsd::ChemOptions& c) { c.ringFusionMode = smsd::ChemOptions::RingFusionMode::STRICT; });
    add("everything", [](smsd::ChemOptions& c) {
        c.useChirality = c.useBondStereo = c.tautomerAware = c.completeRingsOnly = c.ringMatchesRingOnly = true;
        c.matchFormalCharge = true;
        c.aromaticityMode = smsd::ChemOptions::AromaticityMode::STRICT;
    });
    int compared = 0;
    for (const auto& querySmiles : queries) {
        auto query = querySmiles.empty()
            ? smsd::MolGraph::Builder().atomCount(0).atomicNumbers({}).setNeighbors({}).build(false)
            : smsd::parseSMILES(querySmiles);
        // Acetic acid with a wildcard for its methyl carbon.
        if (querySmiles == "CC(=O)O") query.atomicNum[0] = 0;
        for (const auto& targetSmiles : targets) {
            const auto target = smsd::parseSMILES(targetSmiles);
            for (const auto& chemistry : chemistries)
                for (bool induced : {false, true}) {
                    auto options = smsd_test::substructureOptions(chemistry.second);
                    options.induced = induced;
                    const std::string context = "rows " + querySmiles + " in " + targetSmiles + " " + chemistry.first
                        + (induced ? " induced" : "");
                    compared += compareRows(query, target, options, context);
                }
            // Budgets that stop the search: node limits and an immediate time limit.
            for (uint64_t nodes : {uint64_t(1), uint64_t(3), uint64_t(40)}) {
                auto options = smsd_test::substructureOptions({});
                options.nodeLimit = nodes;
                compared += compareRows(query, target, options, "rows " + querySmiles + " in " + targetSmiles
                                        + " nodes " + std::to_string(nodes));
            }
            compared += compareRows(query, target, smsd_test::substructureOptions({}, 0),
                                    "rows " + querySmiles + " in " + targetSmiles + " no time");
        }
    }
    std::cout << compared << " row and collector result pairs agree\n";
}

// A neighbour listed twice, or an atom listed as its own neighbour, repeats complete mappings or
// bonds, and a bond stored differently at its two ends tests differently from each; such graphs
// keep their matches in the collector, which drops the repeats and bonds that fail.
static void repeatedNeighbours() {
    const auto query = smsd::parseSMILES("CC");
    auto repeated = smsd::parseSMILES("CCC");
    repeated.neighbors[1].push_back(0);
    auto looped = smsd::parseSMILES("CCC");
    looped.neighbors[2].push_back(2);
    auto repeatedQuery = smsd::parseSMILES("CC");
    repeatedQuery.neighbors[0].push_back(1);
    const auto chain = smsd::parseSMILES("CCC");
    // A double bond written as single at one end only, through the dense bond matrix.
    const auto doubleBond = smsd::parseSMILES("C=C");
    auto oneWay = smsd::parseSMILES("CC=CC");
    require(oneWay.useDense, "a small graph stores its bonds densely");
    oneWay.bondPacked[2][1] = static_cast<uint8_t>((oneWay.bondPacked[2][1] & ~smsd::MolGraph::BOND_ORDER_MASK) | 1);
    auto oneWayQuery = smsd::parseSMILES("CC=C");
    oneWayQuery.bondPacked[2][1] = static_cast<uint8_t>((oneWayQuery.bondPacked[2][1] & ~smsd::MolGraph::BOND_ORDER_MASK) | 1);
    for (bool induced : {false, true}) {
        auto options = smsd_test::substructureOptions({});
        options.induced = induced;
        compareRows(query, repeated, options, "repeated target neighbour");
        compareRows(query, looped, options, "target self-loop");
        compareRows(repeatedQuery, chain, options, "repeated query neighbour");
        compareRows(doubleBond, oneWay, options, "target bond stored one way");
        compareRows(oneWayQuery, smsd::parseSMILES("CCC=CC"), options, "query bond stored one way");
    }
    // Ethane maps onto propane in four ways, however often a neighbour is listed.
    const auto all = smsd::Substructure::find(query).matchAll(repeated);
    require(all.matches.size() == 4 && all.complete, "a repeated target neighbour must not repeat matches");
    std::cout << "repeated neighbour graphs keep distinct matches\n";
}

// A prepared query compared by ring fusion counts its rings once and searches as a query each search counts.
static void preparedRingCounts() {
    smsd::SubstructureOptions options;
    options.chemistry.ringFusionMode = smsd::ChemOptions::RingFusionMode::STRICT;
    for (const char* querySmiles : {"C1CCCCC1", "c1ccc2ccccc2c1", "C1CC1", "CC"})
        for (const char* targetSmiles : {"C1CCC2CCCCC2C1", "c1ccc2cc3ccccc3cc2c1", "C12CC3CC(CC(C3)C1)C2", "C1CC1C1CC1", "CCO"}) {
            const auto query = smsd::parseSMILES(querySmiles);
            const auto target = smsd::parseSMILES(targetSmiles);
            const auto prepared = smsd::Substructure::find(query).withOptions(options);
            const std::string context = std::string("ring fusion ") + querySmiles + " in " + targetSmiles;
            requireSameResult(smsd::search_detail::substructure(&query, nullptr, target, options, 1000, true),
                              prepared.matchAll(target), context + " all");
            requireSameResult(smsd::search_detail::substructure(&query, nullptr, target, options, 1, false),
                              prepared.match(target), context);
        }
}

int main(int argc, char** argv) {
    return smsd_test::runGroups(argc, argv, {
        {"self", selfEnumeration}, {"path", disconnectedPath}, {"aromatic", strictAromaticBonds},
        {"timeout", stickyDeadline}, {"screen", countScreens}, {"twins", twinScans},
        {"terms", compatibilityTermsParity}, {"recursion", recursiveAllowance},
        {"subpatterns", manySubpatterns}, {"layouts", sharedSubpatternLayouts},
        {"sparse", sparseInducedTimeout}, {"sparse-smarts", sparseInducedSmartsTimeout},
        {"oracle", exhaustiveGraphs}, {"smarts-oracle", smartsOracleGraphs},
        {"fusion", preparedRingCounts}, {"rows", completeMappingRows}, {"repeats", repeatedNeighbours}});
}
