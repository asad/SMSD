/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
#include "smsd/mcs.hpp"
#include "smsd/fingerprint_options.hpp"

#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

// An 80-layer hexagonal prism stack: a synthetic cage whose ring system exceeds the
// ring-perception work limit, so it has no exact ring basis.
const char* stackSmiles = "C12C3C4C5C6C1C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C1C6C5C4C3C21";

void require(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}

/** Four-member rings linked around one larger ring. */
smsd::MolGraph fourRingChain() {
    constexpr int ringCount = 12, atomCount = 1 + 3 * ringCount;
    std::vector<std::vector<int>> neighbors(atomCount), orders(atomCount);
    const auto addBond = [&](int first, int second) {
        neighbors[first].push_back(second);
        neighbors[second].push_back(first);
        orders[first].push_back(1);
        orders[second].push_back(1);
    };
    int previous = 0;
    for (int ring = 0; ring < ringCount; ++ring) {
        const int firstSide = 3 * ring + 1, secondSide = firstSide + 1, next = firstSide + 2;
        addBond(previous, firstSide); addBond(previous, secondSide);
        addBond(firstSide, next); addBond(secondSide, next);
        previous = next;
    }
    addBond(previous, 0);
    return smsd::MolGraph::Builder().atomCount(atomCount)
        .atomicNumbers(std::vector<int>(atomCount, 6)).setNeighbors(std::move(neighbors))
        .setBondOrders(std::move(orders)).build(false);
}

/** The existing topology beside independent atoms, for necessary-element count checks. */
smsd::MolGraph withIsolatedAtoms(const smsd::MolGraph& graph, const std::vector<int>& elements) {
    auto neighbors = graph.neighbors;
    auto atomicNumbers = graph.atomicNum;
    atomicNumbers.insert(atomicNumbers.end(), elements.begin(), elements.end());
    neighbors.resize(atomicNumbers.size());
    std::vector<std::vector<int>> orders(neighbors.size());
    for (size_t atom = 0; atom < neighbors.size(); ++atom) orders[atom].assign(neighbors[atom].size(), 1);
    return smsd::MolGraph::Builder().atomCount(static_cast<int>(atomicNumbers.size()))
        .atomicNumbers(std::move(atomicNumbers)).setNeighbors(std::move(neighbors))
        .setBondOrders(std::move(orders)).build(false);
}

void status(const std::string& id, const smsd::SearchResult& result, bool incomplete) {
    const char* actualStatus = result.status == smsd::SearchStatus::Completed ? "Completed"
        : result.status == smsd::SearchStatus::TimedOut ? "TimedOut"
        : result.status == smsd::SearchStatus::Cancelled ? "Cancelled" : "LimitReached";
    std::cout << "RING_SOURCE_ACTUAL\t" << id << '\t' << actualStatus
              << '\t' << (result.optimal ? "true" : "false")
              << '\t' << (result.complete ? "true" : "false") << std::endl;
    std::cout << "RING_SOURCE_MATCH_ACTUAL\t" << id << '\t' << (!result.matches.empty() ? "true" : "false")
              << '\t' << result.matches.size() << '\t' << result.nodesVisited << std::endl;
    require(result.status == (incomplete ? smsd::SearchStatus::LimitReached : smsd::SearchStatus::Completed), id);
    require(result.optimal == !incomplete && result.complete == !incomplete, id + "/proof");
    std::cout << "RING_SOURCE_PARITY\t" << id << '\t'
              << (incomplete ? "LimitReached" : "Completed") << '\n';
}

void checkFixture(const std::string& name, const smsd::MolGraph& graph, bool basisIncomplete) {
    require(smsd::ringsTruncated(graph) == basisIncomplete, name + "/basis");
    require(smsd::relevantCyclesTruncated(graph), name + "/relevant");
    auto orderedNeighbors = graph.neighbors;
    for (auto& row : orderedNeighbors) std::sort(row.begin(), row.end());
    const auto countSource = smsd::MolGraph::Builder().atomCount(graph.n).atomicNumbers(graph.atomicNum)
        .setNeighbors(std::move(orderedNeighbors)).build(false);
    const bool countsIncomplete = smsd::ringsTruncated(countSource);
    require(graph.ringCountsTruncated() == countsIncomplete, name + "/count-source");
    std::cout << "RING_SOURCE_PARITY\t" << name << "/count-source\t"
              << (countsIncomplete ? "Incomplete" : "Complete") << '\n';
    const auto carbon = smsd::parseSMILES("C");
    for (bool strict : {false, true}) {
        smsd::MCSOptions options;
        options.timeoutMs = 60000;
        options.chemistry.completeRingsOnly = !strict;
        if (strict) options.chemistry.ringFusionMode = smsd::ChemOptions::RingFusionMode::STRICT;
        const std::string id = name + (strict ? "/strict" : "/complete");
        const bool incomplete = strict ? countsIncomplete : smsd::relevantCyclesTruncated(graph);
        status(id + "/find", smsd::MCS::find(graph).withOptions(options).match(carbon), incomplete);
        // With one carbon beside it the consensus fragment stays writable: the cage alone needs
        // more simultaneously open ring closures than SMILES output allows.
        const auto consensus = smsd::MCS::consensus({graph, carbon}, options);
        require(consensus.status == (incomplete ? smsd::SearchStatus::LimitReached : smsd::SearchStatus::Completed), id + "/consensus");
        require(!consensus.optimal && consensus.complete == !incomplete, id + "/consensus-proof");
        std::cout << "RING_SOURCE_PARITY\t" << id << "/consensus\t"
                  << (incomplete ? "LimitReached" : "Completed") << '\n';
    }
    smsd::SubstructureOptions options;
    options.timeoutMs = 60000;
    options.chemistry.completeRingsOnly = true;
    struct SmartsCase { const char* pattern; bool independentAbsence; bool noMatch; };
    const SmartsCase patterns[] = {
        {"[#7;R2]", true, true}, {"[#6;D99]", true, true},
        {"[R2]", false, false}, {"[#6;R2]", false, false},
        {"[#7,#6;R2]", false, false}, {"[!#7;R2]", false, false},
        {"[$([#7]);R2]", true, true}, {"[#7,#7;R2]", true, true},
        {"[#7;#6;R2]", false, true}
    };
    for (const auto& item : patterns) {
        const auto result = smsd::Substructure::smarts(item.pattern).withOptions(options).match(graph);
        const std::string pattern = item.pattern;
        const bool needsRelevantProof = pattern == "[R2]" || pattern == "[#6;R2]"
            || pattern == "[#7,#6;R2]" || pattern == "[!#7;R2]";
        const bool incomplete = !item.independentAbsence
            && (basisIncomplete || (needsRelevantProof && smsd::relevantCyclesTruncated(graph)));
        status(name + "/smarts/" + item.pattern, result, incomplete);
        if (item.noMatch) require(result.matches.empty(), name + "/smarts/no-match");
        if (std::string(item.pattern) == "[#7;R2]" || std::string(item.pattern) == "[#7,#7;R2]"
            || std::string(item.pattern) == "[$([#7]);R2]")
            require(result.nodesVisited == 0, name + "/smarts/element-proof-nodes");
    }
    const auto oneNitrogen = withIsolatedAtoms(graph, {7});
    const auto cardinality = smsd::Substructure::smarts("[#7;R2].[#7]").withOptions(options).match(oneNitrogen);
    status(name + "/element-cardinality", cardinality, false);
    require(cardinality.matches.empty() && cardinality.nodesVisited == 0, name + "/cardinality-proof");
    if (basisIncomplete) {
        std::vector<int> firstSixteen;
        std::string pattern;
        for (int element = 1; element <= 17; ++element) {
            if (element <= 16) firstSixteen.push_back(element);
            if (!pattern.empty()) pattern += '.';
            pattern += "[#" + std::to_string(element) + (element == 17 ? ";R2]" : "]");
        }
        const auto firstSixteenPresent = withIsolatedAtoms(graph, firstSixteen);
        const auto seventeenth = smsd::Substructure::smarts(pattern).withOptions(options).match(firstSixteenPresent);
        status(name + "/seventeenth-element", seventeenth, false);
        require(seventeenth.matches.empty() && seventeenth.nodesVisited == 0, name + "/seventeenth-proof");
    }
    for (int term = 0; term < 3; ++term) {
        auto fingerprint = smsd::batch::FingerprintOptions::ecfp();
        if (term == 0) fingerprint.ringSize = true;
        if (term == 1) fingerprint.cipBonds = true;
        if (term == 2) fingerprint = smsd::batch::FingerprintOptions::withKeys({
            smsd::batch::FingerprintKey::special(smsd::batch::KeyKind::AROMATIC_RINGS)});
        bool refused = false;
        try {
            const auto words = smsd::batch::computeFingerprint(graph, fingerprint);
            if (!basisIncomplete) {
                std::cout << "RING_SOURCE_PARITY\t" << name << "/fingerprint/" << term << "/words\t";
                for (auto word : words) std::cout << word << ',';
                std::cout << '\n';
            }
        } catch (const std::runtime_error& error) {
            refused = true;
            require(std::string(error.what()).find("Exact ring basis unavailable") != std::string::npos,
                    name + "/fingerprint/error");
        }
        require(refused == basisIncomplete, name + "/fingerprint/refusal");
        std::cout << "RING_SOURCE_PARITY\t" << name << "/fingerprint/" << term << '\t'
                  << (refused ? "Refused" : "Completed") << '\n';
    }
    for (int minimum : {1, 2}) {
        const auto keys = smsd::batch::FingerprintOptions::withKeys({smsd::batch::FingerprintKey("[R2]", minimum)});
        bool refused = false;
        try { (void)smsd::batch::computeFingerprint(graph, keys); }
        catch (const std::runtime_error& error) {
            refused = true;
            require(std::string(error.what()).find("Exact ring basis unavailable") != std::string::npos,
                    name + "/key/error");
        }
        require(refused == basisIncomplete, name + "/key/refusal");
        std::cout << "RING_SOURCE_PARITY\t" << name << "/key/" << minimum << '\t'
                  << (refused ? "Refused" : "Completed") << '\n';
        refused = false;
        try { (void)smsd::batch::computeFingerprintCounts(graph, keys); }
        catch (const std::runtime_error& error) {
            refused = true;
            require(std::string(error.what()).find("Exact ring basis unavailable") != std::string::npos,
                    name + "/key-count/error");
        }
        require(refused == basisIncomplete, name + "/key-count/refusal");
        std::cout << "RING_SOURCE_PARITY\t" << name << "/key-count/" << minimum << '\t'
                  << (refused ? "Refused" : "Completed") << '\n';
    }
}

/** Positive recursive anchors supply one necessary root element, without counting inner atoms. */
void checkRecursiveAnchors(const smsd::MolGraph& graph) {
    smsd::SubstructureOptions options;
    options.timeoutMs = 60000;
    options.chemistry.completeRingsOnly = true;
    struct RecursiveCase { const char* pattern; bool absent; };
    const RecursiveCase patterns[] = {
        {"[$([#7]-[#6]);R2]", true}, {"[$([$([#7])]);R2]", true},
        {"[$([#7]),$([#7]-[#6]);R2]", true}, {"[$([#7]),$([#6]);R2]", false},
        {"[!$([#7]);R2]", false}, {"[$([#6]-[#7]);R2]", false},
        {"[$([#6]-[#6]);R2]", false}, {"[$([#7]);$([#6]);R2]", false}
    };
    for (const auto& item : patterns) {
        const auto result = smsd::Substructure::smarts(item.pattern).withOptions(options).match(graph);
        status(std::string("stack/recursive-anchor/") + item.pattern, result, !item.absent);
        require(result.matches.empty() && result.nodesVisited == 0, "recursive-anchor/proof");
    }
    const auto oneNitrogen = withIsolatedAtoms(graph, {7});
    const auto cardinality = smsd::Substructure::smarts("[$([#7]-[#6]);R2].[$([#7]-[#6])]")
        .withOptions(options).match(oneNitrogen);
    status("stack/recursive-anchor/root-cardinality", cardinality, false);
    require(cardinality.matches.empty() && cardinality.nodesVisited == 0, "recursive-anchor/cardinality");
    const auto sharedAnchor = smsd::Substructure::smarts("[$([#7]-[#6]);$([#7]-[#6]);R2]")
        .withOptions(options).match(oneNitrogen);
    status("stack/recursive-anchor/shared-anchor", sharedAnchor, true);
    require(sharedAnchor.matches.empty() && sharedAnchor.nodesVisited == 0, "recursive-anchor/shared-anchor");
}

} // namespace

int main() {
    try {
        const auto stack = smsd::parseSMILES(stackSmiles);
        checkFixture("stack", stack, true);
        checkFixture("relevant-only", fourRingChain(), false);
        checkRecursiveAnchors(stack);
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
    return 0;
}
