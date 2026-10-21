/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 *
 * Public search API measurements; validation and JSON formatting are outside timers.
 * --threads N runs batch rows on N threads (0 for all hardware threads); --batch-targets FILE
 * fills each batch with distinct molecules read from FILE instead of copies of the case target.
 */
#include "smsd/mcs.hpp"
#include "smsd/substructure.hpp"
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <optional>
#include <sstream>
#include <variant>

namespace {
using Clock = std::chrono::steady_clock;

/** A call shorter than one clock tick reads 0 ns: record 1 ns and flag it; it is never a speed result. */
struct Elapsed { int64_t ns; bool belowTimerResolution; };

Elapsed elapsedSince(Clock::time_point started) {
    const int64_t ns = std::chrono::duration_cast<std::chrono::nanoseconds>(Clock::now() - started).count();
    return ns == 0 ? Elapsed{1, true} : Elapsed{ns, false};
}

void elapsedJson(std::ostream& output, const Elapsed& elapsed) {
    output << ",\"elapsed_ns\":" << elapsed.ns;
    if (elapsed.belowTimerResolution) output << ",\"below_timer_resolution\":true";
}
using Prepared = std::variant<smsd::MCS, smsd::Substructure>;
struct Case {
    std::string id, operation, queryText, targetText;
    smsd::MCSOptions common;
    smsd::SubstructureOptions substructure;
    size_t limit = 0, batchSize = 1;
    std::vector<std::string> phases;
};

std::vector<std::string> split(const std::string& text, char separator) {
    std::vector<std::string> values;
    size_t start = 0;
    for (;;) {
        const size_t end = text.find(separator, start);
        values.push_back(text.substr(start, end == std::string::npos ? end : end - start));
        if (end == std::string::npos) return values;
        start = end + 1;
    }
}

int64_t integer(const std::string& text) {
    size_t end = 0;
    const auto value = std::stoll(text, &end);
    if (end != text.size()) throw std::invalid_argument("Invalid integer: " + text);
    return value;
}

smsd::ChemOptions chemistry(const std::string& profile) {
    smsd::ChemOptions options;
    options.matchAtomType = true;
    options.matchBondOrder = smsd::ChemOptions::BondOrderMode::STRICT;
    options.aromaticityMode = smsd::ChemOptions::AromaticityMode::FLEXIBLE;
    options.matchFormalCharge = options.matchIsotope = options.useChirality = options.useBondStereo = false;
    options.ringMatchesRingOnly = options.completeRingsOnly = options.tautomerAware = false;
    options.ringFusionMode = smsd::ChemOptions::RingFusionMode::IGNORE;
    if (profile == "strict" || profile == "rings") {
        options.aromaticityMode = smsd::ChemOptions::AromaticityMode::STRICT;
        options.matchFormalCharge = options.matchIsotope = options.ringMatchesRingOnly = true;
        options.completeRingsOnly = profile == "rings";
    } else if (profile == "any") options.matchBondOrder = smsd::ChemOptions::BondOrderMode::ANY;
    else if (profile == "tautomer") {
        options.tautomerAware = true;
        options.matchBondOrder = smsd::ChemOptions::BondOrderMode::LOOSE;
    } else if (profile == "stereo") options.useChirality = options.useBondStereo = true;
    else if (profile != "default") throw std::invalid_argument("Unknown chemistry profile: " + profile);
    return options;
}

std::vector<Case> readCases(const std::string& filename) {
    std::ifstream input(filename);
    if (!input) throw std::runtime_error("Cannot read cases file");
    std::vector<Case> cases;
    std::string line;
    bool header = false;
    while (std::getline(input, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty() || line.front() == '#') continue;
        const auto fields = split(line, '\t');
        if (!header) {
            const auto expected = split("case_id\toperation\tobjective\tconnectivity\tprofile\tinduced\ttimeout_ms\tlimit\tbatch_size\tquery\ttarget\tweights\tphases", '\t');
            if (fields != expected) throw std::invalid_argument("Unexpected cases header");
            header = true;
            continue;
        }
        if (fields.size() != 13) throw std::invalid_argument("Expected 13 case fields");
        Case item; item.id = fields[0]; item.operation = fields[1];
        try {
            if (item.operation != "mcs" && item.operation != "sub" && item.operation != "smarts")
                throw std::invalid_argument("Unknown operation");
            if (fields[2] != "MCIS" && fields[2] != "MCES") throw std::invalid_argument("Unknown objective");
            item.common.objective = fields[2] == "MCIS" ? smsd::MCSObjective::MCIS : smsd::MCSObjective::MCES;
            if (fields[3] != "Connected" && fields[3] != "Fragments") throw std::invalid_argument("Unknown connectivity");
            item.common.connectivity = fields[3] == "Connected" ? smsd::MCSConnectivity::Connected : smsd::MCSConnectivity::Fragments;
            item.common.chemistry = item.substructure.chemistry = chemistry(fields[4]);
            if (fields[5] != "0" && fields[5] != "1") throw std::invalid_argument("Invalid induced flag");
            item.substructure.induced = fields[5] == "1";
            item.common.timeoutMs = item.substructure.timeoutMs = integer(fields[6]);
            const auto limit = integer(fields[7]), batchSize = integer(fields[8]);
            if (limit < 0 || batchSize < 1) throw std::invalid_argument("Invalid result or batch limit");
            item.limit = static_cast<size_t>(limit); item.batchSize = static_cast<size_t>(batchSize);
            item.queryText = fields[9]; item.targetText = fields[10];
            if (!fields[11].empty()) for (const auto& text : split(fields[11], ',')) {
                size_t end = 0;
                const double weight = std::stod(text, &end);
                if (end != text.size() || !std::isfinite(weight)) throw std::invalid_argument("Invalid atom weight");
                item.common.atomWeights.push_back(weight);
            }
            item.phases = split(fields[12], ',');
            for (const auto& phase : item.phases) {
                if (phase != "prepare" && phase != "prepared" && phase != "end_to_end" && phase != "batch")
                    throw std::invalid_argument("Unknown phase");
                if (phase == "batch" && item.limit != 0) throw std::invalid_argument("Batch rows require limit zero");
            }
            cases.push_back(std::move(item));
        } catch (const std::exception& error) { throw std::invalid_argument(item.id + ": " + error.what()); }
    }
    if (!header || cases.empty()) throw std::invalid_argument("No benchmark cases");
    return cases;
}

Prepared prepare(const Case& item, const smsd::MolGraph& query) {
    if (item.operation == "mcs") return smsd::MCS::find(query).withOptions(item.common);
    if (item.operation == "sub") return smsd::Substructure::find(query).withOptions(item.substructure);
    return smsd::Substructure::smarts(item.queryText).withOptions(item.substructure);
}

smsd::SearchResult match(const Case& item, const Prepared& query, const smsd::MolGraph& target) {
    return std::visit([&](const auto& prepared) {
        return item.limit == 0 ? prepared.match(target) : prepared.matchAll(target, item.limit);
    }, query);
}

void quoted(std::ostream& output, const std::string& value) {
    output << '"';
    constexpr char hex[] = "0123456789abcdef";
    for (unsigned char character : value) {
        if (character == '"' || character == '\\') output << '\\' << character;
        else if (character < 32) output << "\\u00" << hex[character >> 4] << hex[character & 15];
        else output << character;
    }
    output << '"';
}

void graphJson(std::ostream& output, const smsd::MolGraph& graph) {
    output << "{\"atoms\":[";
    for (int atom = 0; atom < graph.n; ++atom) {
        if (atom != 0) output << ',';
        output << "{\"element\":" << graph.atomicNum[atom] << ",\"charge\":" << graph.formalCharge[atom]
            << ",\"isotope\":" << graph.massNumber[atom] << ",\"hydrogens\":" << graph.hydrogenCount[atom]
            << ",\"aromatic\":" << bool(graph.aromatic[atom]) << ",\"ring\":" << bool(graph.ring[atom]) << '}';
    }
    output << "],\"bonds\":[";
    bool first = true;
    for (int atom = 0; atom < graph.n; ++atom) {
        auto adjacent = graph.neighbors[atom]; std::sort(adjacent.begin(), adjacent.end());
        for (int neighbor : adjacent) if (atom < neighbor) {
            if (!first) output << ','; first = false;
            output << '[' << atom << ',' << neighbor << ',' << graph.bondOrder(atom, neighbor) << ','
                << graph.bondAromatic(atom, neighbor) << ',' << graph.bondInRing(atom, neighbor) << ']';
        }
    }
    output << "]}";
}

const char* statusName(smsd::SearchStatus status) {
    switch (status) {
        case smsd::SearchStatus::Completed: return "Completed";
        case smsd::SearchStatus::TimedOut: return "TimedOut";
        case smsd::SearchStatus::Cancelled: return "Cancelled";
        case smsd::SearchStatus::LimitReached: return "LimitReached";
    }
    throw std::logic_error("Unknown search status");
}

void resultJson(std::ostream& output, const smsd::SearchResult& result) {
    output << "{\"status\":"; quoted(output, statusName(result.status));
    output << ",\"optimal\":" << result.optimal << ",\"complete\":" << result.complete
        << ",\"score\":" << result.score << ",\"nodes_visited\":" << result.nodesVisited << ",\"matches\":[";
    for (size_t index = 0; index < result.matches.size(); ++index) {
        if (index != 0) output << ',';
        const auto& selected = result.matches[index];
        output << "{\"atoms\":[";
        bool first = true;
        for (const auto& atom : selected.atomMapping) {
            if (!first) output << ','; first = false;
            output << '[' << atom.first << ',' << atom.second << ']';
        }
        output << "],\"bonds\":[";
        for (size_t bond = 0; bond < selected.bondMapping.size(); ++bond) {
            if (bond != 0) output << ',';
            const auto& pair = selected.bondMapping[bond];
            output << '[' << pair[0] << ',' << pair[1] << ',' << pair[2] << ',' << pair[3] << ']';
        }
        output << "],\"fragments\":[";
        for (size_t fragment = 0; fragment < selected.fragments.size(); ++fragment) {
            if (fragment != 0) output << ',';
            output << '[';
            for (size_t atom = 0; atom < selected.fragments[fragment].size(); ++atom) {
                if (atom != 0) output << ',';
                output << selected.fragments[fragment][atom];
            }
            output << ']';
        }
        output << "]}";
    }
    output << "]}";
}

void validate(const Case& item, const smsd::MolGraph& query, const smsd::MolGraph& target,
              const smsd::SearchResult& result, const std::optional<smsd::SmartsQuery>& pattern) {
    if (result.matches.size() > (item.limit == 0 ? 1 : item.limit)) throw std::runtime_error("Result limit exceeded");
    if (result.complete && result.status != smsd::SearchStatus::Completed) throw std::runtime_error("Stopped result claims completeness");
    smsd::SearchOptions validationOptions; validationOptions.timeoutMs = 60000;
    smsd::search_detail::Budget budget(validationOptions, pattern ? pattern->atomCount() : query.n, target.n, false);
    smsd::search_detail::RecursiveScope recursiveScope(budget);
    const auto view = pattern ? smsd::search_detail::QueryView(*pattern) : smsd::search_detail::QueryView(query);
    for (const auto& selected : result.matches) {
        if (item.operation == "mcs") {
            const auto errors = smsd::MCS::validateMapping(query, target, selected, item.common);
            if (!errors.empty()) throw std::runtime_error("Invalid MCS mapping: " + errors.front());
            const auto score = smsd::search_detail::scoreMatch(selected, item.common);
            const double primary = item.common.atomWeights.empty() ? double(score.primary) : score.primary / 1000.0;
            if (primary != result.score) throw std::runtime_error("Mapping and result scores differ");
        } else {
            if (selected.atomMapping.size() != static_cast<size_t>(view.size())) throw std::runtime_error("Partial substructure mapping");
            std::vector<int> mapping(view.size(), -1);
            std::set<int> used;
            for (const auto& atom : selected.atomMapping) {
                if (atom.first < 0 || atom.first >= view.size() || atom.second < 0 || atom.second >= target.n
                    || !used.insert(atom.second).second || !view.atomMatches(atom.first, target, atom.second, item.substructure.chemistry))
                    throw std::runtime_error("Invalid substructure atom mapping");
                mapping[atom.first] = atom.second;
            }
            std::vector<std::array<int, 4>> expected;
            for (int atom = 0; atom < view.size(); ++atom)
                for (int other = atom + 1; other < view.size(); ++other) {
                    if (view.adjacent(atom, other)) {
                        if (!view.bondMatches(atom, other, target, mapping[atom], mapping[other], item.substructure.chemistry))
                            throw std::runtime_error("Invalid substructure bond mapping");
                        expected.push_back({atom, other, mapping[atom], mapping[other]});
                    } else if (item.substructure.induced && target.hasBond(mapping[atom], mapping[other]))
                        throw std::runtime_error("Invalid induced substructure mapping");
                }
            auto actual = selected.bondMapping; std::sort(actual.begin(), actual.end());
            if (expected != actual) throw std::runtime_error("Substructure selected bonds differ from query");
            if (pattern ? !pattern->stereoCompatible(target, mapping)
                        : !smsd::search_detail::selectedStereoCompatible(query, target, item.substructure.chemistry, selected))
                throw std::runtime_error("Invalid substructure stereochemistry");
        }
        std::vector<int> mapping(view.size(), -1);
        for (const auto& atom : selected.atomMapping) mapping[atom.first] = atom.second;
        auto expectedFragments = smsd::search_detail::makeMatch(mapping, selected.bondMapping).fragments;
        auto actualFragments = selected.fragments;
        for (auto& fragment : actualFragments) std::sort(fragment.begin(), fragment.end());
        std::sort(actualFragments.begin(), actualFragments.end());
        std::sort(expectedFragments.begin(), expectedFragments.end());
        if (actualFragments != expectedFragments) throw std::runtime_error("Invalid selected-edge fragments");
    }
    if (budget.stopped()) throw std::runtime_error("Validation budget expired");
}

/** One SMILES per line, optionally after an identifier and a tab; blank lines and # comments are skipped. */
std::vector<std::string> readBatchTargets(const std::string& filename) {
    std::ifstream input(filename);
    if (!input) throw std::runtime_error("Cannot read batch targets file");
    std::vector<std::string> targets;
    std::string line;
    while (std::getline(input, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty() || line.front() == '#') continue;
        targets.push_back(split(line, '\t').back());
    }
    if (targets.empty()) throw std::invalid_argument("No batch targets");
    return targets;
}

/** How batch rows run: the thread count and, when given, the distinct targets that fill a batch in turn. */
struct BatchSettings {
    int threads = 1;
    std::vector<std::string> distinctTargets;
};

void run(const Case& item, const std::vector<std::string>& phases, int warmup, int repeats, const BatchSettings& batch,
         std::ostream& output) {
    const auto query = item.operation == "smarts" ? smsd::MolGraph{} : smsd::parseSMILES(item.queryText);
    const auto target = smsd::parseSMILES(item.targetText);
    const bool distinct = !batch.distinctTargets.empty();
    const auto pattern = item.operation == "smarts" ? std::optional<smsd::SmartsQuery>(smsd::parseSMARTS(item.queryText)) : std::nullopt;
    const std::string graphPhase = std::find(phases.begin(), phases.end(), "prepared") != phases.end() ? "prepared" : phases.front();
    for (const auto& phase : phases) {
        std::optional<Prepared> prepared;
        if (phase == "prepared" || phase == "batch") prepared = prepare(item, query);
        std::vector<smsd::MolGraph> targets;
        if (phase == "batch") {
            if (!distinct) targets.assign(item.batchSize, target);
            else for (size_t index = 0; index < item.batchSize; ++index)
                targets.push_back(smsd::parseSMILES(batch.distinctTargets[index % batch.distinctTargets.size()]));
        }
        for (int sample = -warmup; sample < repeats; ++sample) {
            std::vector<smsd::SearchResult> results;
            smsd::SearchResult scalar;
            const auto started = Clock::now();
            if (phase == "prepare") {
                auto constructed = prepare(item, query);
                const auto elapsed = elapsedSince(started);
                prepared = std::move(constructed);
                if (sample < 0) continue;
                output << "{\"case_id\":"; quoted(output, item.id);
                output << ",\"language\":\"cpp\",\"phase\":\"prepare\",\"sample\":" << sample
                    << ",\"threads\":1";
                elapsedJson(output, elapsed);
                output << ",\"results\":[]";
            } else {
                if (phase == "end_to_end") {
                    const auto parsedQuery = item.operation == "smarts" ? smsd::MolGraph{} : smsd::parseSMILES(item.queryText);
                    const auto parsedTarget = smsd::parseSMILES(item.targetText);
                    const auto constructed = prepare(item, parsedQuery);
                    scalar = match(item, constructed, parsedTarget);
                } else if (phase == "batch") {
                    results = std::visit([&](const auto& constructed) { return constructed.matchMany(targets, batch.threads); },
                                         *prepared);
                } else scalar = match(item, *prepared, target);
                const auto elapsed = elapsedSince(started);
                if (phase != "batch") results.push_back(std::move(scalar));
                for (size_t index = 0; index < results.size(); ++index)
                    validate(item, query, phase == "batch" ? targets[index] : target, results[index], pattern);
                if (sample < 0) continue;
                output << "{\"case_id\":"; quoted(output, item.id);
                output << ",\"language\":\"cpp\",\"phase\":"; quoted(output, phase);
                output << ",\"sample\":" << sample << ",\"threads\":" << (phase == "batch" ? batch.threads : 1);
                if (phase == "batch") output << ",\"batch_targets\":\"" << (distinct ? "distinct" : "copies") << '"';
                elapsedJson(output, elapsed);
                output << ",\"results\":[";
                for (size_t index = 0; index < results.size(); ++index) {
                    if (index != 0) output << ',';
                    resultJson(output, results[index]);
                }
                output << ']';
                if (phase == "batch" && distinct && sample == 0) {
                    output << ",\"batch_target_graphs\":[";
                    for (size_t index = 0; index < targets.size(); ++index) {
                        if (index != 0) output << ',';
                        graphJson(output, targets[index]);
                    }
                    output << ']';
                }
            }
            output << ",\"query_graph\":";
            if (sample == 0 && phase == graphPhase && item.operation != "smarts") graphJson(output, query); else output << "null";
            output << ",\"target_graph\":";
            if (sample == 0 && phase == graphPhase) graphJson(output, target); else output << "null";
            output << "}\n";
        }
        output.flush();
    }
}
} // namespace

int main(int argc, char** argv) {
    try {
        std::string casesFile, outputFile, selectedPhase, selectedCase, batchTargetsFile;
        int warmup = 0, repeats = 1;
        BatchSettings batch;
        for (int index = 1; index < argc; index += 2) {
            if (index + 1 >= argc) throw std::invalid_argument("Expected option value");
            const std::string option = argv[index], value = argv[index + 1];
            if (option == "--cases") casesFile = value;
            else if (option == "--output") outputFile = value;
            else if (option == "--phase") selectedPhase = value;
            else if (option == "--case") selectedCase = value;
            else if (option == "--batch-targets") batchTargetsFile = value;
            else if (option == "--threads") {
                const auto count = integer(value);
                if (count < 0 || count > 4096) throw std::invalid_argument("Invalid thread count");
                batch.threads = static_cast<int>(count);
            }
            else if (option == "--warmup" || option == "--repeats") {
                const auto count = integer(value);
                if (count < 0 || count > INT_MAX || (option == "--repeats" && count == 0)) throw std::invalid_argument("Invalid sample count");
                (option == "--warmup" ? warmup : repeats) = static_cast<int>(count);
            } else throw std::invalid_argument("Unknown option: " + option);
        }
        if (casesFile.empty() || outputFile.empty()) throw std::invalid_argument("--cases and --output are required");
        const auto cases = readCases(casesFile);
        if (!batchTargetsFile.empty()) batch.distinctTargets = readBatchTargets(batchTargetsFile);
        std::ofstream output(outputFile);
        if (!output) throw std::runtime_error("Cannot write output file");
        output << std::boolalpha << std::setprecision(17);
        bool found = false;
        for (const auto& item : cases) {
            if (!selectedCase.empty() && item.id != selectedCase) continue;
            std::vector<std::string> phases;
            for (const auto& phase : item.phases) if (selectedPhase.empty() || phase == selectedPhase) phases.push_back(phase);
            if (phases.empty()) continue;
            found = true;
            try { run(item, phases, warmup, repeats, batch, output); }
            catch (const std::exception& error) { throw std::runtime_error(item.id + ": " + error.what()); }
        }
        if (!found) throw std::invalid_argument("No cases selected");
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
