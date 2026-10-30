/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
#include "smsd/mcs.hpp"
#include "smsd/substructure.hpp"
#include "smsd/fingerprint_options.hpp"
#include <map>
#include <set>
#include <algorithm>
#include <iterator>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

static int failures = 0;
static void require(bool value, const char* message) { if (!value) throw std::runtime_error(message); }
template<class Operation> static void control(const std::string& name, Operation operation) {
    try { operation(); std::cout << "PASS\t" << name << '\n'; }
    catch (const std::exception& error) { ++failures; std::cout << "FAIL\t" << name << '\t' << error.what() << '\n'; }
}
static std::string hoop(int units) {
    std::string text = "c%101ccc(cc1)";
    for (int unit = 1; unit + 1 < units; ++unit) text += "c1ccc(cc1)";
    return text + "c1ccc%10cc1";
}
static void extraControls();

int main() {
    struct Prism { const char* id; const char* query; const char* target; int score; };
    for (const auto& test : {Prism{"pentaprismane", "C1CCCC1", "C12C3C4C5C1C6C2C3C4C56", 5},
                            Prism{"hexaprismane", "C1CCCCC1", "C12C3C4C5C6C1C7C2C3C4C5C67", 6},
                            Prism{"bridged-cage", "C1CCCC1", "C12C3C4C1C5C3C4C25", 5}}) {
        for (auto objective : {smsd::MCSObjective::MCIS, smsd::MCSObjective::MCES}) {
            const std::string name = std::string(test.id) + (objective == smsd::MCSObjective::MCIS ? "/MCIS" : "/MCES");
            control(name, [&] {
                smsd::MCSOptions options;
                options.timeoutMs = 10000;
                options.objective = objective;
                options.chemistry.completeRingsOnly = true;
                const auto result = smsd::MCS::find(test.query).withOptions(options).match(test.target);
                std::cout << "RESULT\t" << name << '\t' << result.score << '\t' << int(result.status)
                          << '\t' << result.optimal << '\t' << result.complete << '\n';
                require(result.status == smsd::SearchStatus::Completed && result.optimal && result.complete
                        && result.score == test.score, "complete-ring prism/cage optimum changed");
            });
        }
    }
    const auto macrocycle = smsd::parseSMILES(hoop(10));
    control("cpp10/numeric-R2", [&] {
        const auto result = smsd::Substructure::smarts("[R2]").withTimeout(10000).match(macrocycle);
        require(result.status == smsd::SearchStatus::Completed && result.optimal && !result.matches.empty(),
                "full minimum basis must support numeric ring membership");
    });
    control("cpp10/MACCS", [&] {
        const auto bits = smsd::batch::computeFingerprint(macrocycle, smsd::batch::FingerprintOptions::maccs());
        require(!bits.empty(), "MACCS fingerprint must remain available");
    });
    control("cpp10/ringSize", [&] {
        auto options = smsd::batch::FingerprintOptions::ecfp();
        options.ringSize = true;
        require(!smsd::batch::computeFingerprint(macrocycle, options).empty(), "ring-size fingerprint must remain available");
    });
    control("cpp10/legacy-counts", [&] {
        macrocycle.ensureRingCounts();
        std::cout << "COUNTS\tcpp10";
        for (int count : macrocycle.ringCount) std::cout << '\t' << count;
        std::cout << '\n';
    });
    extraControls();
    std::cout << "SUMMARY\t" << failures << '\n';
    return failures == 0 ? 0 : 1;
}

static void extraControls() {
    control("ring-result-aggregate-source-compatibility", [&] {
        const std::vector<std::vector<int>> smallest{{0, 1, 2}}, relevant{{0, 1, 2}};
        const smsd::RingResult legacy{smallest, relevant, true};
        require(legacy.sssr == smallest && legacy.rcb == relevant && legacy.truncated,
                "the established three-field aggregate initializer retains its meaning");
    });
    control("ring-count-storage-after-mutation", [&] {
        auto graph = smsd::parseSMILES("CCC");
        graph.ensureRingCounts();
        require(graph.ringCount.size() == 3, "ring counts cover the original chain");
        graph.setBondOrder(0, 2, 1);
        graph.ensureRingCounts();
        require(graph.ringCount.size() == 3, "recomputed ring counts retain logical atom coverage");
        require(graph.ringCount == std::vector<int>({1, 1, 1}), "closing the chain creates one ring through each atom");
        smsd::MCSOptions options;
        options.timeoutMs = 10000;
        options.chemistry.ringFusionMode = smsd::ChemOptions::RingFusionMode::STRICT;
        auto result = smsd::MCS::find("C1CC1").withOptions(options).match(graph);
        require(result.optimal && result.complete && result.score == 3, "strict fusion reads recomputed counts");
        graph.setBondOrder(0, 2, 0);
        graph.ensureRingCounts();
        require(graph.ringCount == std::vector<int>({0, 0, 0}), "reopening the ring clears all memberships");
        result = smsd::MCS::find("CCC").withOptions(options).match(graph);
        require(result.optimal && result.complete && result.score == 3, "strict fusion reads reopened-chain counts");
    });
    control("basis-and-relevant-flags", [&] {
        const auto graph = smsd::parseSMILES(hoop(10));
        const auto& basis = smsd::detail_rf::completeBasis(graph);
        const auto saved = basis;
        require(!smsd::ringsTruncated(graph), "the macrocycle has a complete minimum basis");
        (void)smsd::computeRelevantCycles(graph);
        require(basis == saved, "relevant enumeration retains the previously returned basis reference");
        require(!smsd::ringsTruncated(graph) && smsd::relevantCyclesTruncated(graph),
                "relevant enumeration limits do not truncate a complete minimum basis");
        require(graph.ringCountsTruncated() == false, "legacy fusion counts have their own complete source");
    });
    control("numeric-ring-admission-and-boolean-context", [&] {
        const auto graph = smsd::parseSMILES("C12C3C4C5C6C1C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C1C6C5C4C3C21");
        require(smsd::ringsTruncated(graph), "the 80-layer stack exhausts basis perception");
        for (const char* pattern : {"[R1]", "[!R1]", "[R1,#6]", "[$([R1])]", "[!$([R1])]"}) {
            const auto result = smsd::Substructure::smarts(pattern).withTimeout(10000).match(graph);
            require(result.status == smsd::SearchStatus::LimitReached && !result.optimal && !result.complete
                    && result.matches.empty(), "numeric basis predicates refuse before boolean evaluation");
            const auto stopped = smsd::Substructure::smarts(pattern).withTimeout(0).match(graph);
            require(stopped.status == smsd::SearchStatus::TimedOut && stopped.matches.empty(),
                    "the existing immediate budget stop has priority over numeric basis admission");
        }
        for (const char* pattern : {"[R]", "[R0]", "[!R0]", "[$([R0])]", "[#6]"}) {
            const auto result = smsd::Substructure::smarts(pattern).withTimeout(10000).match(graph);
            require(result.status == smsd::SearchStatus::Completed && result.optimal && result.complete,
                    "topology-only ring and element predicates retain their exact results");
        }
    });
    control("minimum-basis-certificate", [&] {
        using Cycles = std::vector<std::vector<int>>;
        const auto weight = [](const Cycles& cycles) { size_t w = 0; for (const auto& c : cycles) w += c.size(); return w; };
        const uint64_t ample = uint64_t(1) << 34;
        std::string cd8 = "C%101C(O)C(O)C(OC1CO)O";
        for (int unit = 1; unit < 8; ++unit) cd8 += "C1C(O)C(O)C(OC1CO)O";
        cd8 += "%10";
        for (const std::string smiles : {std::string("c1ccc2ccccc2c1"), std::string("C1C2CC3CC1CC(C2)C3"),
                                         std::string("C12C3C4C1C5C2C3C45"), std::string("C12C3C4C5C6C1C7C2C3C4C5C67"),
                                         std::string("c1cc2ccc3ccc4ccc5ccc6ccc1c7c2c3c4c5c67"), hoop(8), cd8}) {
            const auto graph = smsd::parseSMILES(smiles);
            const int rank = smsd::detail::cycleRank(graph);
            const auto exact = smsd::detail_rf::perceiveRingSets(graph, false, ample);
            require(!exact.basisTruncated && smsd::detail_rf::certifiedMinimumBasis(graph, exact.sssr, rank, ample),
                    "an exhaustive minimum basis is certified");
            const auto bounded = smsd::detail_rf::perceiveRingSets(graph, false, 0);
            if (bounded.truncated) {
                require(bounded.basisTruncated, "a certificate without work proves nothing");
                const bool minimum = static_cast<int>(bounded.sssr.size()) == rank && weight(bounded.sssr) == weight(exact.sssr);
                require(smsd::detail_rf::certifiedMinimumBasis(graph, bounded.sssr, rank, ample) == minimum,
                        "a fallback basis is certified exactly when it is a minimum basis");
            }
            // Replacing a member by its sum with a neighbouring member keeps a basis, which stays
            // minimum only when the sum is as short as the member.
            for (size_t a = 0; a < exact.sssr.size(); ++a)
                for (size_t b = 0; b < exact.sssr.size(); ++b) {
                    if (a == b) continue;
                    std::set<int64_t> first, second, sum;
                    for (size_t i = 0; i < exact.sssr[a].size(); ++i)
                        first.insert(smsd::detail_rf::edgeKey(exact.sssr[a][i], exact.sssr[a][(i + 1) % exact.sssr[a].size()]));
                    for (size_t i = 0; i < exact.sssr[b].size(); ++i)
                        second.insert(smsd::detail_rf::edgeKey(exact.sssr[b][i], exact.sssr[b][(i + 1) % exact.sssr[b].size()]));
                    std::set_symmetric_difference(first.begin(), first.end(), second.begin(), second.end(),
                                                  std::inserter(sum, sum.begin()));
                    if (sum.size() + 1 >= first.size() + second.size()) continue;  // no shared bond
                    std::map<int, std::vector<int>> adjacent;
                    for (auto key : sum) {
                        const int u = static_cast<int>(key >> 32), v = static_cast<int>(static_cast<uint32_t>(key));
                        adjacent[u].push_back(v); adjacent[v].push_back(u);
                    }
                    bool simple = true;
                    for (const auto& entry : adjacent) simple = simple && entry.second.size() == 2;
                    if (!simple) continue;
                    std::vector<int> cycle;
                    for (int previous = -1, atom = adjacent.begin()->first; cycle.size() < sum.size(); ) {
                        cycle.push_back(atom);
                        const int next = adjacent[atom][0] == previous ? adjacent[atom][1] : adjacent[atom][0];
                        previous = atom; atom = next;
                        if (atom == cycle.front()) break;
                    }
                    if (cycle.size() != sum.size()) continue;
                    Cycles swapped = exact.sssr;
                    swapped[a] = cycle;
                    require(smsd::detail_rf::certifiedMinimumBasis(graph, swapped, rank, ample) == (cycle.size() == exact.sssr[a].size()),
                            "a heavier basis is never certified");
                }
            if (rank > 1) {
                Cycles dependent = exact.sssr;
                dependent.back() = dependent.front();
                require(!smsd::detail_rf::certifiedMinimumBasis(graph, dependent, rank, ample), "a dependent set is never certified");
            }
        }
    });
    control("large-cyclodextrin-fallback-is-exact", [&] {
        std::string smiles = "C%101C(O)C(O)C(OC1CO)O";
        for (int unit = 1; unit < 34; ++unit) smiles += "C1C(O)C(O)C(OC1CO)O";
        smiles += "%10";
        const auto graph = smsd::parseSMILES(smiles);
        const auto rings = smsd::detail_rf::perceiveRingSets(graph, false, uint64_t(1) << 20);
        size_t weight = 0;
        for (const auto& cycle : rings.sssr) weight += cycle.size();
        require(rings.truncated && !rings.basisTruncated && rings.sssr.size() == 35 && weight == 34 * 11,
                "a certified fallback basis of 34 glucose rings and the macrocycle is exact");
        require(!smsd::ringsTruncated(graph) && smsd::detail_rf::completeBasis(graph).size() == 35,
                "the large cyclodextrin has an exact basis at the default work limit");
        require(smsd::RING_CANDIDATE_WORK_LIMIT == uint64_t(1) << 29, "the documented ring-perception work limit");
    });
}
