/* SPDX-License-Identifier: Apache-2.0 */
#include "smsd/mcs.hpp"
#include <future>
#include <iostream>
#include <limits>
#include <stdexcept>

static void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

static smsd::Match identitySelection(const smsd::MolGraph& query) {
    smsd::Match result;
    for (int atom = 0; atom < query.n; ++atom) {
        result.atomMapping.emplace(atom, atom);
        for (int neighbor : query.neighbors[atom])
            if (atom < neighbor) result.bondMapping.push_back({atom, neighbor, atom, neighbor});
    }
    return result;
}

int main() {
    try {
        const auto query = smsd::parseSMILES("CC(=O)C"), target = smsd::parseSMILES("CC(O)=C");
        const auto match = identitySelection(query);
        const auto queryWeights = query.tautomerWeight, targetWeights = target.tautomerWeight;
        const auto aqueous = smsd::ChemOptions::tautomerProfile();
        auto chloroform = aqueous; chloroform.withSolvent(smsd::ChemOptions::Solvent::CHLOROFORM);
        auto acidic = aqueous; acidic.withPH(2);
        const double baseline = smsd::MCS::computeTautomerConfidence(query, target, match, aqueous);
        const double solvent = smsd::MCS::computeTautomerConfidence(query, target, match, chloroform);
        const double lowPH = smsd::MCS::computeTautomerConfidence(query, target, match, acidic);
        require(baseline > 0 && baseline < 1, "Selected keto/enol differences must affect confidence");
        require(solvent < baseline && lowPH < baseline, "Solvent and pH must affect the explicit confidence");
        require(std::abs(baseline - 0.98) < 1e-6 && std::abs(solvent - 0.73) < 1e-6
                && std::abs(lowPH - 0.5) < 1e-6, "Context formula differs from its baseline solvent and pH contract");
        for (int repeat = 0; repeat < 4; ++repeat) {
            require(smsd::MCS::computeTautomerConfidence(query, target, match, chloroform) == solvent,
                    "Repeated solvent confidence must not compound");
            require(smsd::MCS::computeTautomerConfidence(query, target, match, aqueous) == baseline,
                    "Context order must not alter the aqueous confidence");
        }
        for (auto context : {aqueous, chloroform, acidic})
            require(smsd::MCS::computeTautomerConfidence(query, query, match, context) == 1,
                    "Identical selected chemistry must have confidence one");
        auto unchanged = match;
        unchanged.bondMapping.erase(std::remove_if(unchanged.bondMapping.begin(), unchanged.bondMapping.end(),
            [&](const auto& bond) { return query.bondOrder(bond[0], bond[1]) != target.bondOrder(bond[2], bond[3]); }),
            unchanged.bondMapping.end());
        require(smsd::MCS::computeTautomerConfidence(query, target, unchanged, aqueous) == 1,
                "Unselected bond differences must not penalise confidence");
        std::vector<std::future<double>> workers;
        for (int worker = 0; worker < 8; ++worker)
            workers.push_back(std::async(std::launch::async, [&] {
                return smsd::MCS::computeTautomerConfidence(query, target, match, chloroform);
            }));
        for (auto& worker : workers) require(worker.get() == solvent, "Concurrent context reads must agree");
        require(query.tautomerWeight == queryWeights && target.tautomerWeight == targetWeights,
                "Context confidence must not mutate molecule weights");
        auto contextualQuery = query;
        contextualQuery.applySolventCorrection(smsd::ChemOptions::Solvent::CHLOROFORM);
        require(smsd::MCS::computeTautomerConfidence(contextualQuery, target, match, chloroform) == solvent,
                "Existing graph context must not be applied twice");
        std::unordered_map<int, int> legacy(match.atomMapping.begin(), match.atomMapping.end());
        require(smsd::computeTautomerConfidence(query, target, legacy) == 1,
                "Historical element-only confidence must retain its contract");
        for (double value : {std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::infinity(),
                             -std::numeric_limits<double>::infinity()}) {
            auto invalid = aqueous; invalid.pH = value;
            bool rejected = false;
            try { smsd::MCS::computeTautomerConfidence(query, target, match, invalid); }
            catch (const std::invalid_argument&) { rejected = true; }
            require(rejected, "Non-finite confidence pH must be rejected");
            rejected = false;
            try { invalid.withPH(value); } catch (const std::invalid_argument&) { rejected = true; }
            require(rejected, "Non-finite fluent pH must be rejected");
        }
        auto malformed = match; malformed.atomMapping[query.n] = 0;
        bool rejected = false;
        try { smsd::MCS::computeTautomerConfidence(query, target, malformed, aqueous); }
        catch (const std::out_of_range&) { rejected = true; }
        require(rejected, "Out-of-range selected atoms must be rejected");
        legacy[query.n] = 0; rejected = false;
        try { smsd::computeTautomerConfidence(query, target, legacy); }
        catch (const std::out_of_range&) { rejected = true; }
        require(rejected, "Historical confidence must reject out-of-range atoms");
        std::cout << "PASS selected chemistry, context, repeatability, concurrency and bounds\n";
    } catch (const std::exception& error) {
        std::cerr << "FAIL " << error.what() << '\n';
        return 1;
    }
}
