/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
#include "smsd/substructure.hpp"

#include <array>
#include <future>
#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

static void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

static void completed(const smsd::SearchResult& result, bool found, int atoms) {
    require(result.status == smsd::SearchStatus::Completed && result.optimal && result.complete,
            "The query must finish with proof");
    require((!result.matches.empty()) == found, "The complete query predicates must constrain the result");
    require(result.score == (found ? atoms : 0), "The full-query score must match its atom count");
}

static smsd::Substructure query(const std::string& pattern, smsd::SubstructureOptions options = {}) {
    options.timeoutMs = 10000;
    return smsd::Substructure::smarts(pattern).withOptions(options);
}

static smsd::smarts_detail::BondExpr directionDomain(int domain) {
    using smsd::smarts_detail::BondExpr;
    using smsd::smarts_detail::BondPrimType;
    std::vector<BondExpr> alternatives;
    if (domain & 1) alternatives.push_back(BondExpr::makePrim(BondPrimType::WEDGE));
    if (domain & 2) alternatives.push_back(BondExpr::makePrim(BondPrimType::DASH));
    if (domain & 4) alternatives.push_back(BondExpr::makeAnd({
        BondExpr::makePrim(BondPrimType::WEDGE_OR_UNSPECIFIED),
        BondExpr::makePrim(BondPrimType::DASH_OR_UNSPECIFIED)}));
    return alternatives.size() == 1 ? alternatives.front() : BondExpr::makeOr(std::move(alternatives));
}

static bool containsDirection(int domain, int direction) {
    return (domain & (direction == 0 ? 4 : direction > 0 ? 1 : 2)) != 0;
}

/** Enumerate all three-valued assignments rather than repeating the parity solver. */
static bool branchedAssignment(const std::array<int, 3>& domains, int configuration) {
    for (int first : {-1, 0, 1}) for (int second : {-1, 0, 1}) for (int third : {-1, 0, 1}) {
        if (!containsDirection(domains[0], first) || !containsDirection(domains[1], second)
            || !containsDirection(domains[2], third)) continue;
        if (first && second && first != second) continue;
        if (configuration == 0) {
            if (!first || !second || !third) return true;
        } else if (first && second && third && first * third == (configuration == 1 ? -1 : 1))
            return true;
    }
    return false;
}

static bool sharedAssignment(const std::array<int, 3>& domains, int firstConfiguration, int secondConfiguration) {
    for (int first : {-1, 0, 1}) for (int shared : {-1, 0, 1}) for (int last : {-1, 0, 1}) {
        if (!containsDirection(domains[0], first) || !containsDirection(domains[1], shared)
            || !containsDirection(domains[2], last)) continue;
        const bool firstMatches = firstConfiguration == 0 ? (!first || !shared)
            : first && shared && first * shared == (firstConfiguration == 1 ? -1 : 1);
        const bool secondMatches = secondConfiguration == 0 ? (!shared || !last)
            : shared && last && shared * last == (secondConfiguration == 1 ? -1 : 1);
        if (firstMatches && secondMatches) return true;
    }
    return false;
}

static void directionDomains() {
    auto branched = smsd::parseSMARTS("F/C(/Cl)=C/Br");
    auto branchedTarget = smsd::parseSMILES("FC(Cl)=CBr");
    auto shared = smsd::parseSMARTS("F/C=C/C=C/Cl");
    auto sharedTarget = smsd::parseSMILES("FC=CC=CCl");
    const std::vector<int> branchedMapping{0, 1, 2, 3, 4}, sharedMapping{0, 1, 2, 3, 4, 5};
    for (int first = 1; first <= 7; ++first) for (int second = 1; second <= 7; ++second)
        for (int third = 1; third <= 7; ++third) {
            const std::array<int, 3> domains{first, second, third};
            branched.bonds[0].expr = directionDomain(first);
            branched.bonds[1].expr = directionDomain(second);
            branched.bonds[3].expr = directionDomain(third);
            shared.bonds[0].expr = directionDomain(first);
            shared.bonds[2].expr = directionDomain(second);
            shared.bonds[4].expr = directionDomain(third);
            for (int configuration = 0; configuration <= 2; ++configuration) {
                branchedTarget.hasDbStereo = configuration != 0;
                branchedTarget.dbStereoConf[1][3] = branchedTarget.dbStereoConf[3][1] = configuration;
                require(branched.stereoCompatible(branchedTarget, branchedMapping)
                        == branchedAssignment(domains, configuration),
                        "Same-end optional directions must agree with the three-valued assignments");
                for (int next = 0; next <= 2; ++next) {
                    sharedTarget.hasDbStereo = configuration != 0 || next != 0;
                    sharedTarget.dbStereoConf[1][2] = sharedTarget.dbStereoConf[2][1] = configuration;
                    sharedTarget.dbStereoConf[3][4] = sharedTarget.dbStereoConf[4][3] = next;
                    require(shared.stereoCompatible(sharedTarget, sharedMapping)
                            == sharedAssignment(domains, configuration, next),
                            "A shared carrier must satisfy both configured or unspecified focuses");
                }
            }
        }
}

int main() {
    try {
        const auto triangle = smsd::parseSMILES("C1=CC1");
        completed(query("[$(C=CC)]").match(triangle), true, 1);
        std::cout << "Ordinary recursive control passed\n" << std::flush;

        completed(query("C1=C-1").match(triangle), false, 2);
        completed(query("[$(C1=C-1)]").match(triangle), false, 1);
        completed(query("[!$(C1=C-1)]").match(triangle), true, 1);
        completed(query("[#6,$(C1=C-1)]").match(triangle), true, 1);
        completed(query("[$(C1=C-1),#6]").match(triangle), true, 1);

        const auto chain = smsd::parseSMILES("CC");
        const auto duplicate = query("C1-C-1").matchAll(chain, 100);
        completed(duplicate, true, 2);
        require(duplicate.matches.size() == 2, "Both compatible parallel-bond embeddings must remain");
        std::set<std::map<int, int>> mappings;
        for (const auto& match : duplicate.matches) {
            require(match.atomMapping.size() == 2 && match.bondMapping.size() == 1
                    && match.fragments == std::vector<std::vector<int>>{{0, 1}},
                    "A parallel predicate must not duplicate the selected edge or fragments");
            require(mappings.insert(match.atomMapping).second, "Embeddings must be distinct");
        }
        const auto capped = query("C1-C-1").matchAll(chain, 1);
        require(capped.matches.size() == 1 && capped.optimal && !capped.complete
                && capped.status == smsd::SearchStatus::LimitReached,
                "A capped list must preserve its known objective without claiming complete enumeration");

        for (auto mode : {smsd::ChemOptions::BondOrderMode::STRICT,
                          smsd::ChemOptions::BondOrderMode::ANY,
                          smsd::ChemOptions::BondOrderMode::LOOSE}) {
            smsd::SubstructureOptions options;
            options.chemistry.matchAtomType = false;
            options.chemistry.matchBondOrder = mode;
            options.chemistry.tautomerAware = mode == smsd::ChemOptions::BondOrderMode::LOOSE;
            completed(query("C1=C-1", options).match(triangle), false, 2);
            completed(query("C1-C-1", options).match(chain), true, 2);
            completed(query("C11", options).match(triangle), false, 1);
        }
        completed(query("[$(C11)]").match(triangle), false, 1);
        completed(query("[!$(C11)]").match(triangle), true, 1);
        completed(query("[#6,$(C11)]").match(triangle), true, 1);

        auto manual = smsd::parseSMARTS("CC");
        auto extra = manual.bonds.front();
        extra.from = 1;
        extra.to = 0;
        extra.expr = smsd::smarts_detail::BondExpr::makePrim(smsd::smarts_detail::BondPrimType::DOUBLE);
        manual.bonds.push_back(extra);
        manual.buildAdjacency();
        smsd::SubstructureOptions options;
        options.timeoutMs = 10000;
        completed(smsd::Substructure::smarts(manual).withOptions(options).match(triangle), false, 2);
        auto outer = smsd::parseSMARTS("[$(CC)]");
        outer.atoms.front().expr.cachedRecursive = std::make_shared<smsd::SmartsQuery>(manual);
        completed(smsd::Substructure::smarts(outer).withOptions(options).match(triangle), false, 1);
        completed(query("[$(C1-C-1),$(C1=C-1)]").match(chain), true, 1);
        completed(query("[$(C1=C-1)].[$(C1-C-1)]").match(triangle), false, 2);
        auto uncached = smsd::parseSMARTS("[$(C1=C-1)]");
        uncached.atoms.front().expr.cachedRecursive.reset();
        completed(smsd::Substructure::smarts(uncached).withOptions(options).match(triangle), false, 1);
        const auto reusable = query("[$(C1-C-1),$(C1=C-1)]");
        auto firstCall = std::async(std::launch::async, [&] { return reusable.match(chain); });
        auto secondCall = std::async(std::launch::async, [&] { return reusable.match(triangle); });
        completed(firstCall.get(), true, 1);
        completed(secondCall.get(), true, 1);
        completed(reusable.match(chain), true, 1);

        completed(query(R"(F1/C\1=C/F)").match("F/C=C/F"), true, 4);
        completed(query(R"(F1/C\1=C/F)").match(R"(F/C=C\F)"), false, 4);
        completed(query("F1/C/1=C/F").match("F/C=C/F"), false, 4);
        completed(query(R"(F1/?C\1=C/F)").match("FC=CF"), false, 4);
        completed(query(R"(F1/?C\?1=C/F)").match("FC=CF"), true, 4);
        completed(query(R"(F1/?C/?1=C/F)").match("F/C=C/F"), false, 4);
        completed(query(R"(F1/?C/?1=C/F)").match("FC=CF"), true, 4);
        completed(query(R"(F/?&\?C=C/F)").match("FC=CF"), true, 4);
        completed(query(R"(F/?&\?C=C/F)").match("F/C=C/F"), false, 4);
        completed(query("F/C=C/F").match("FC=CF"), false, 4);
        completed(query(R"(F/?&\?C(/Cl)=C/Br)").match("FC(Cl)=CBr"), true, 5);
        completed(query(R"(F/?&\?C(/Cl)=C/Br)").match("F/C(/Cl)=C/Br"), false, 5);
        directionDomains();
        completed(query("C[C@H](F)Cl.C1-C-1").match("C[C@H](F)Cl.CC"), true, 6);
        completed(query("C[C@H](F)Cl.C1-C-1").match("C[C@@H](F)Cl.CC"), false, 6);

        auto cancellation = std::make_shared<smsd::CancellationToken>();
        cancellation->cancel();
        options.cancellation = cancellation;
        const auto stopped = query("C11", options).match(triangle);
        require(stopped.status == smsd::SearchStatus::Cancelled && !stopped.optimal && !stopped.complete,
                "A pre-cancelled self-bond query must not manufacture a proof");
        std::cout << "SMARTS parallel-predicate regressions passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
