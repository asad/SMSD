/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 *
 * Shared chemistry-rule fixture, also read by the Java and Python tests.
 */
#include "smsd/smsd.hpp"
#include "smsd/cip.hpp"
#include "smsd/similarity_screen.hpp"
#include "smsd/stereo.hpp"

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#ifndef SMSD_CHEMISTRY_RULES_JSON
#error "SMSD_CHEMISTRY_RULES_JSON must name the shared fixture"
#endif

namespace {

/** A JSON value; only what the fixture uses. */
struct Json {
    enum class Kind { Null, Bool, Number, String, Array, Object } kind = Kind::Null;
    bool boolean = false;
    double number = 0;
    std::string text;
    std::vector<Json> items;
    std::map<std::string, Json> members;

    const Json& operator[](const std::string& key) const {
        const auto found = members.find(key);
        if (kind != Kind::Object || found == members.end()) throw std::runtime_error("Missing key " + key);
        return found->second;
    }
    bool has(const std::string& key) const { return kind == Kind::Object && members.count(key) != 0; }
    int integer() const { return static_cast<int>(number); }
    std::vector<int> integers() const {
        std::vector<int> values;
        for (const auto& item : items) values.push_back(item.integer());
        return values;
    }
};

class JsonReader {
    const std::string& source_;
    size_t at_ = 0;

    void space() { while (at_ < source_.size() && std::isspace(static_cast<unsigned char>(source_[at_]))) ++at_; }
    char peek() { space(); if (at_ >= source_.size()) throw std::runtime_error("Unexpected end of JSON"); return source_[at_]; }
    void expect(char c) { if (peek() != c) throw std::runtime_error(std::string("Expected ") + c); ++at_; }
    std::string string() {
        expect('"');
        std::string value;
        while (source_.at(at_) != '"') {
            if (source_[at_] == '\\') {
                const char escaped = source_.at(++at_);
                value += escaped == 'n' ? '\n' : escaped == 't' ? '\t' : escaped;
            } else {
                value += source_[at_];
            }
            ++at_;
        }
        ++at_;
        return value;
    }

public:
    explicit JsonReader(const std::string& source) : source_(source) {}

    Json value() {
        Json result;
        const char c = peek();
        if (c == '{') {
            result.kind = Json::Kind::Object;
            ++at_;
            if (peek() == '}') { ++at_; return result; }
            for (;;) {
                const std::string key = string();
                expect(':');
                result.members[key] = value();
                if (peek() == ',') { ++at_; continue; }
                expect('}');
                return result;
            }
        }
        if (c == '[') {
            result.kind = Json::Kind::Array;
            ++at_;
            if (peek() == ']') { ++at_; return result; }
            for (;;) {
                result.items.push_back(value());
                if (peek() == ',') { ++at_; continue; }
                expect(']');
                return result;
            }
        }
        if (c == '"') { result.kind = Json::Kind::String; result.text = string(); return result; }
        if (source_.compare(at_, 4, "true") == 0) { at_ += 4; result.kind = Json::Kind::Bool; result.boolean = true; return result; }
        if (source_.compare(at_, 5, "false") == 0) { at_ += 5; result.kind = Json::Kind::Bool; return result; }
        if (source_.compare(at_, 4, "null") == 0) { at_ += 4; return result; }
        size_t used = 0;
        result.kind = Json::Kind::Number;
        result.number = std::stod(source_.substr(at_, 32), &used);
        at_ += used;
        return result;
    }
};

int failures = 0;

void check(bool condition, const std::string& message) {
    if (!condition) {
        ++failures;
        std::cerr << "FAIL " << message << '\n';
    }
}

void near(double actual, double expected, double tolerance, const std::string& message) {
    std::ostringstream text;
    text.precision(17);
    text << message << ": expected " << expected << ", got " << actual;
    check(std::fabs(actual - expected) <= tolerance, text.str());
}

smsd::MolGraph molecule(const std::string& smiles) {
    return smiles.empty() ? smsd::MolGraph() : smsd::parseSMILES(smiles);
}

/** The same graph built without hydrogen counts, so every count follows the valence model. */
smsd::MolGraph withoutHydrogenCounts(const smsd::MolGraph& graph) {
    std::vector<std::vector<int>> orders(graph.n);
    std::vector<uint8_t> ring(graph.n), aromatic(graph.n);
    for (int atom = 0; atom < graph.n; ++atom) {
        ring[atom] = graph.ring[atom];
        aromatic[atom] = graph.aromatic[atom];
        for (int neighbor : graph.neighbors[atom]) orders[atom].push_back(graph.bondOrder(atom, neighbor));
    }
    return smsd::MolGraph::Builder().atomCount(graph.n).atomicNumbers(graph.atomicNum)
        .formalCharges(graph.formalCharge).ringFlags(ring).aromaticFlags(aromatic)
        .setNeighbors(graph.neighbors).setBondOrders(orders).build();
}

std::vector<int> matchedAtoms(const std::string& smarts, const smsd::MolGraph& target) {
    std::vector<int> atoms;
    for (const auto& match : smsd::Substructure::smarts(smarts).matchAll(target).matches)
        atoms.push_back(match.atomMapping.at(0));
    std::sort(atoms.begin(), atoms.end());
    return atoms;
}

/** The R/S labels of a graph by atom, as the fixture writes them. */
std::map<std::string, std::string> cipLabels(const smsd::MolGraph& graph) {
    const auto labels = smsd::cip::assignAll(graph).rsLabels;
    std::map<std::string, std::string> found;
    for (int atom = 0; atom < graph.n; ++atom)
        if (labels[atom] == smsd::cip::RSLabel::R || labels[atom] == smsd::cip::RSLabel::S)
            found[std::to_string(atom)] = labels[atom] == smsd::cip::RSLabel::R ? "R" : "S";
    return found;
}

/** All four descriptors of a graph by atom, the pseudoasymmetric r and s included. */
std::map<std::string, std::string> allCipLabels(const smsd::MolGraph& graph) {
    const auto labels = smsd::cip::assignAll(graph).rsLabels;
    std::map<std::string, std::string> found;
    for (int atom = 0; atom < graph.n; ++atom)
        if (labels[atom] != smsd::cip::RSLabel::NONE) found[std::to_string(atom)] = smsd::cip::to_string(labels[atom]);
    return found;
}

std::map<std::string, std::string> expectedLabels(const Json& item) {
    std::map<std::string, std::string> expected;
    for (const auto& label : item["labels"].members) expected[label.first] = label.second.text;
    return expected;
}

std::string describe(const std::vector<int>& values) {
    std::ostringstream text;
    text << '[';
    for (size_t i = 0; i < values.size(); ++i) text << (i ? ", " : "") << values[i];
    text << ']';
    return text.str();
}

void similarity(const Json& section) {
    for (const auto& item : section["cases"].items) {
        const std::string name = "similarity " + item["name"].text;
        const int size = item["fp_size"].integer();
        const auto first = smsd::bitWords(item["first"].integers(), size);
        const auto second = smsd::bitWords(item["second"].integers(), size);
        // An empty fingerprint is the same set whether it has zero words or none.
        std::vector<std::pair<std::vector<uint64_t>, std::vector<uint64_t>>> forms = {{first, second}};
        if (item["first"].items.empty()) forms.push_back({{}, second});
        if (item["second"].items.empty()) forms.push_back({first, {}});
        for (const auto& form : forms) {
            near(smsd::batch::fingerprintTanimoto(form.first, form.second), item["tanimoto"].number, 1e-12, name + " Tanimoto");
            near(smsd::batch::fingerprintOverlap(form.first, form.second), item["overlap"].number, 1e-12, name + " overlap");
            near(smsd::batch::fingerprintDice(form.first, form.second), item["dice"].number, 1e-12, name + " Dice");
            near(smsd::batch::fingerprintCosine(form.first, form.second), item["cosine"].number, 1e-12, name + " cosine");
            near(smsd::batch::fingerprintSoergel(form.first, form.second), item["soergel"].number, 1e-12, name + " Soergel");
        }
        const auto counts = smsd::batch::bitCounts(item["first"].integers(), item["second"].integers());
        near(smsd::batch::tanimoto(counts), item["tanimoto"].number, 1e-12, name + " Tanimoto of positions");
        near(smsd::batch::overlap(counts), item["overlap"].number, 1e-12, name + " overlap of positions");
    }
}

void circularSimilarity(const Json& section) {
    for (const auto& item : section["cases"].items) {
        const auto query = smsd::batch::detail::computeCircularFingerprintECFP(molecule(item["query"].text), 2, 2048);
        for (size_t index = 0; index < item["targets"].items.size(); ++index) {
            const auto& target = item["targets"].items[index].text;
            const auto fingerprint = smsd::batch::detail::computeCircularFingerprintECFP(molecule(target), 2, 2048);
            near(smsd::batch::fingerprintTanimoto(query, fingerprint), item["tanimoto"].items[index].number, 1e-12,
                 "circular similarity '" + item["query"].text + "' to '" + target + "'");
        }
    }
}

// The deprecated helpers stay under test, so their warnings are silenced here alone.
#if defined(_MSC_VER)
#pragma warning(push)
#pragma warning(disable : 4996)
#else
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#endif
int legacyImplicitHydrogens(int atomicNumber, int bondOrderSum, int charge) {
    return smsd::implicitHydrogens(atomicNumber, bondOrderSum, charge);
}

double legacyConfidence(const smsd::MolGraph& query, const smsd::MolGraph& target,
                        const std::unordered_map<int, int>& atoms) {
    return smsd::computeTautomerConfidence(query, target, atoms);
}
#if defined(_MSC_VER)
#pragma warning(pop)
#else
#pragma GCC diagnostic pop
#endif

void implicitHydrogens(const Json& section) {
    for (const auto& item : section["cases"].items) {
        const int hydrogens = smsd::detail::computeImplicitH(item["atomic_number"].integer(), item["aromatic"].boolean,
                                                             item["bond_order_sum"].integer(), item["charge"].integer());
        check(hydrogens == item["hydrogens"].integer(),
              "implicit hydrogens " + item["note"].text + ": got " + std::to_string(hydrogens));
        if (!item["aromatic"].boolean)
            check(legacyImplicitHydrogens(item["atomic_number"].integer(), item["bond_order_sum"].integer(),
                                          item["charge"].integer()) == item["hydrogens"].integer(),
                  "periodic-table implicit hydrogens " + item["note"].text);
    }
}

void smarts(const Json& section, bool omitCounts) {
    for (const auto& item : section["cases"].items) {
        const auto parsed = molecule(item["smiles"].text);
        const auto atoms = matchedAtoms(item["smarts"].text, omitCounts ? withoutHydrogenCounts(parsed) : parsed);
        check(atoms == item["atoms"].integers(), std::string(omitCounts ? "SMARTS without counts " : "SMARTS ")
              + item["name"].text + " " + item["smarts"].text + ": got " + describe(atoms));
    }
}

smsd::ChemOptions::Solvent solvent(const std::string& name) {
    using Solvent = smsd::ChemOptions::Solvent;
    static const std::map<std::string, Solvent> solvents = {
        {"AQUEOUS", Solvent::AQUEOUS}, {"DMSO", Solvent::DMSO}, {"METHANOL", Solvent::METHANOL},
        {"CHLOROFORM", Solvent::CHLOROFORM}, {"ACETONITRILE", Solvent::ACETONITRILE},
        {"DIETHYL_ETHER", Solvent::DIETHYL_ETHER}};
    return solvents.at(name);
}

void tautomerConfidence(const Json& section) {
    const double tolerance = section["tolerance"].number;
    std::map<std::string, double> confidences;
    for (const auto& item : section["cases"].items) {
        const std::string name = "tautomer confidence " + item["name"].text;
        const auto query = molecule(item["query"].text), target = molecule(item["target"].text);
        smsd::Match match;
        for (const auto& pair : item["atoms"].items) match.atomMapping[pair.items[0].integer()] = pair.items[1].integer();
        std::unordered_map<int, int> atoms(match.atomMapping.begin(), match.atomMapping.end());
        auto chemistry = smsd::ChemOptions::tautomerProfile();
        chemistry.solvent = solvent(item["solvent"].text);
        chemistry.pH = item["pH"].number;
        if (item.has("error")) {
            bool fourRejected = false, threeRejected = false;
            try { smsd::MCS::computeTautomerConfidence(query, target, match, chemistry); }
            catch (const std::out_of_range&) { fourRejected = true; }
            try { legacyConfidence(query, target, atoms); }
            catch (const std::out_of_range&) { threeRejected = true; }
            check(fourRejected && threeRejected, name + " must be rejected by both forms");
            continue;
        }
        for (const auto& atom : match.atomMapping)
            for (int neighbor : query.neighbors[atom.first]) {
                const auto other = match.atomMapping.find(neighbor);
                if (neighbor > atom.first && other != match.atomMapping.end() && target.hasBond(atom.second, other->second))
                    match.bondMapping.push_back({atom.first, neighbor, atom.second, other->second});
            }
        const double confidence = smsd::MCS::computeTautomerConfidence(query, target, match, chemistry);
        confidences[item["name"].text] = confidence;
        near(confidence, item["confidence"].number, tolerance, name);
        if (item.has("three_argument"))
            near(legacyConfidence(query, target, atoms), item["three_argument"].number, tolerance,
                 name + " (three-argument form)");
    }
    // A case with same_as differs from the named case only in the order of the atoms.
    for (const auto& item : section["cases"].items) {
        if (!item.has("same_as")) continue;
        const auto reference = confidences.find(item["same_as"].text);
        check(reference != confidences.end(), "tautomer confidence " + item["name"].text + " names no case");
        if (reference != confidences.end())
            near(confidences.at(item["name"].text), reference->second, tolerance,
                 "tautomer confidence " + item["name"].text + " against " + item["same_as"].text);
    }
}

void pathFingerprints(const Json& section) {
    for (const auto& item : section["cases"].items) {
        const std::string name = "fingerprint " + item["name"].text;
        const auto graph = molecule(item["smiles"].text);
        const int length = item["path_length"].integer(), pathSize = item["path_fp_size"].integer();
        const int batchSize = item["batch_fp_size"].integer();
        const auto path = smsd::bitPositions(smsd::batch::detail::computeElementPathFingerprint(graph, length, pathSize), pathSize);
        check(path == item["path_fingerprint"].integers(), name + " path bits " + describe(path));
        const auto batch = smsd::bitPositions(smsd::batch::batchFingerprint({graph}, length, batchSize, 1).at(0));
        check(batch == item["batch_fingerprint"].integers(), name + " batch bits " + describe(batch));
        smsd::ChemOptions chemistry;
        chemistry.tautomerAware = item["tautomer_aware"].boolean;
        const auto mcs = smsd::bitPositions(
            smsd::batch::detail::computeMCSFingerprintForChemistry(graph, length, pathSize, chemistry), pathSize);
        check(mcs == item["mcs_fingerprint"].integers(), name + " MCS bits " + describe(mcs));
    }
}

std::vector<std::vector<int>> integerRows(const Json& rows) {
    std::vector<std::vector<int>> result;
    for (const auto& row : rows.items) result.push_back(row.integers());
    return result;
}

void tetrahedralCentres(const Json& section) {
    for (const auto& item : section["graphs"].items) {
        const std::string name = "tetrahedral centre " + item["name"].text;
        const auto elements = item["atomic_numbers"].integers();
        const int count = static_cast<int>(elements.size());
        auto builder = smsd::MolGraph::Builder().atomCount(count).atomicNumbers(elements)
            .setNeighbors(integerRows(item["neighbors"])).setBondOrders(integerRows(item["orders"]));
        if (item.has("hydrogens")) builder.hydrogenCounts(item["hydrogens"].integers());
        std::vector<int> tags(count, 0);
        tags[item["atom"].integer()] = 1;
        const auto graph = builder.tetrahedralChirality(tags).build();
        std::array<int, 4> ligands{};
        const bool centre = smsd::stereo_detail::tetrahedralLigands(graph, item["atom"].integer(), ligands);
        if (item["ligands"].kind == Json::Kind::Null) {
            check(!centre, name + " is not a centre");
        } else {
            const std::vector<int> found(ligands.begin(), ligands.end());
            check(centre && found == item["ligands"].integers(), name + " ligands " + describe(found));
        }
    }
    smsd::SubstructureOptions options;
    options.chemistry.useChirality = true;
    for (const auto& item : section["matches"].items) {
        const auto query = molecule(item["query"].text), target = molecule(item["target"].text);
        const bool found = !smsd::Substructure::find(query).withOptions(options).match(target).matches.empty();
        check(found == item["found"].boolean, "chiral substructure " + item["name"].text);
    }
    for (const auto& item : section["labels"].items)
        check(cipLabels(molecule(item["smiles"].text)) == expectedLabels(item),
              "CIP labels of " + item["smiles"].text);
}

void mancudeLabels(const Json& section) {
    for (const auto& item : section["cases"].items) {
        const auto graph = molecule(item["smiles"].text);
        check(cipLabels(graph) == expectedLabels(item), "mancude CIP labels of " + item["name"].text);
        // The same graph without aromatic flags, in the Kekule structure the core assigns.
        auto kekule = graph;
        check(kekule.kekulize() && cipLabels(kekule) == expectedLabels(item),
              "mancude CIP labels of the Kekule form of " + item["name"].text);
    }
}

void auxiliaryDescriptors(const Json& section) {
    for (const auto& item : section["cases"].items)
        check(allCipLabels(molecule(item["smiles"].text)) == expectedLabels(item),
              "Rule 4 and 5 CIP labels of " + item["name"].text);
}

std::vector<std::pair<int, int>> nonZero(const std::vector<int>& counts) {
    std::vector<std::pair<int, int>> pairs;
    for (size_t i = 0; i < counts.size(); ++i)
        if (counts[i] != 0) pairs.emplace_back(static_cast<int>(i), counts[i]);
    return pairs;
}

std::vector<std::pair<int, int>> pairsOf(const Json& items) {
    std::vector<std::pair<int, int>> pairs;
    for (const auto& item : items.items) pairs.emplace_back(item.items.at(0).integer(), item.items.at(1).integer());
    return pairs;
}

std::string describe(const std::vector<std::pair<int, int>>& pairs) {
    std::ostringstream text;
    text << '[';
    for (size_t i = 0; i < pairs.size(); ++i) text << (i ? ", " : "") << '[' << pairs[i].first << ", " << pairs[i].second << ']';
    text << ']';
    return text.str();
}

/** The options of a fixture case: the classic configuration of its family, then its keyword arguments. */
smsd::batch::FingerprintOptions fingerprintOptionsOf(const Json& spec) {
    using smsd::batch::FingerprintOptions;
    const std::string family = spec["family"].text;
    FingerprintOptions options = spec.has("maccs") && spec["maccs"].boolean ? FingerprintOptions::maccs()
        : FingerprintOptions::of(family == "ecfp" ? smsd::batch::FingerprintFamily::ECFP
            : family == "fcfp" ? smsd::batch::FingerprintFamily::FCFP
            : family == "path" ? smsd::batch::FingerprintFamily::PATH
            : family == "torsion" ? smsd::batch::FingerprintFamily::TORSION
            : smsd::batch::FingerprintFamily::KEYS);
    const std::map<std::string, bool FingerprintOptions::*> switches = {
        {"element", &FingerprintOptions::element}, {"degree", &FingerprintOptions::degree},
        {"valence", &FingerprintOptions::valence}, {"formal_charge", &FingerprintOptions::formalCharge},
        {"isotope", &FingerprintOptions::isotope}, {"hydrogen_count", &FingerprintOptions::hydrogenCount},
        {"tautomer_class", &FingerprintOptions::tautomerClass}, {"aromaticity", &FingerprintOptions::aromaticity},
        {"bond_order", &FingerprintOptions::bondOrder}, {"aromatic_bond", &FingerprintOptions::aromaticBond},
        {"ring_bond", &FingerprintOptions::ringBond}, {"ring_membership", &FingerprintOptions::ringMembership},
        {"ring_size", &FingerprintOptions::ringSize}, {"cip_atoms", &FingerprintOptions::cipAtoms},
        {"cip_bonds", &FingerprintOptions::cipBonds},
    };
    for (const auto& entry : spec.members) {
        const std::string& key = entry.first;
        const Json& value = entry.second;
        if (key == "family" || key == "maccs") continue;
        const auto flag = switches.find(key);
        if (flag != switches.end()) { options.*(flag->second) = value.boolean; continue; }
        if (key == "radius") options.radius = value.integer();
        else if (key == "path_length") options.pathLength = value.integer();
        else if (key == "size") options.size = value.integer();
        else if (key == "ring_size_cap") options.ringSizeCap = value.integer();
        else if (key == "fragments")
            options.fragments = value.text == "largest" ? smsd::batch::FragmentMode::LARGEST
                : value.text == "per_fragment" ? smsd::batch::FragmentMode::PER_FRAGMENT : smsd::batch::FragmentMode::WHOLE;
        else if (key == "scaffold")
            options.scaffold = value.text == "murcko" ? smsd::batch::ScaffoldMode::MURCKO : smsd::batch::ScaffoldMode::NONE;
        else if (key == "keys") {
            options.keys.clear();
            for (const auto& item : value.items) options.keys.emplace_back(item.items.at(0).text, item.items.at(1).integer());
            options.size = static_cast<int>(options.keys.size());
        } else {
            throw std::runtime_error("Unknown fingerprint option " + key);
        }
    }
    return options;
}

void fingerprintOptions(const Json& section) {
    using namespace smsd::batch;
    for (const auto& item : section["cases"].items) {
        const std::string name = "fingerprint options " + item["name"].text;
        const auto graph = molecule(item["smiles"].text);
        const auto options = fingerprintOptionsOf(item["options"]);
        const int size = options.family == FingerprintFamily::KEYS ? static_cast<int>(options.keys.size()) : options.size;
        const auto bits = computeFingerprint(graph, options);
        const auto positions = smsd::bitPositions(bits, size);
        check(positions == item["bits"].integers(), name + " bits " + describe(positions));
        const auto counts = computeFingerprintCounts(graph, options);
        check(nonZero(counts) == pairsOf(item["counts"]), name + " counts " + describe(nonZero(counts)));
        if (!item.has("classic") || !item["classic"].boolean) continue;
        // The classic configuration of a family gives the bits of its fixed function.
        switch (options.family) {
        case FingerprintFamily::ECFP:
            check(bits == detail::computeCircularFingerprintECFP(graph, options.radius, size), name + " equals ECFP");
            check(counts == detail::computeCircularFingerprintECFPCounts(graph, options.radius, size), name + " equals ECFP counts");
            break;
        case FingerprintFamily::FCFP:
            check(bits == detail::computeCircularFingerprintFCFP(graph, options.radius, size), name + " equals FCFP");
            check(counts == detail::computeCircularFingerprintFCFPCounts(graph, options.radius, size), name + " equals FCFP counts");
            break;
        case FingerprintFamily::PATH:
            check(bits == detail::computeElementPathFingerprint(graph, options.pathLength, size), name + " equals path fingerprint");
            break;
        case FingerprintFamily::TORSION:
            check(bits == detail::computeTopologicalTorsion(graph, size), name + " equals torsion fingerprint");
            check(counts == detail::computeTopologicalTorsionCounts(graph, size), name + " equals torsion counts");
            break;
        case FingerprintFamily::KEYS:
            break;
        }
    }
}

void daylightSmarts(const Json& section) {
    for (const auto& item : section["cases"].items) {
        const auto target = molecule(item["smiles"].text);
        const auto found = smsd::Substructure::smarts(item["smarts"].text).matchAll(target).matches.size();
        check(static_cast<int>(found) == item["embeddings"].integer(),
              "Daylight SMARTS " + item["name"].text + " " + item["smarts"].text + ": got " + std::to_string(found));
    }
}

std::string enumName(smsd::ChemOptions::BondOrderMode mode) {
    switch (mode) {
        case smsd::ChemOptions::BondOrderMode::STRICT: return "STRICT";
        case smsd::ChemOptions::BondOrderMode::LOOSE: return "LOOSE";
        default: return "ANY";
    }
}
std::string enumName(smsd::ChemOptions::AromaticityMode mode) {
    return mode == smsd::ChemOptions::AromaticityMode::STRICT ? "STRICT" : "FLEXIBLE";
}
std::string enumName(smsd::ChemOptions::RingFusionMode mode) {
    switch (mode) {
        case smsd::ChemOptions::RingFusionMode::IGNORE: return "IGNORE";
        case smsd::ChemOptions::RingFusionMode::PERMISSIVE: return "PERMISSIVE";
        default: return "STRICT";
    }
}
std::string enumName(smsd::ChemOptions::Solvent value) {
    switch (value) {
        case smsd::ChemOptions::Solvent::AQUEOUS: return "AQUEOUS";
        case smsd::ChemOptions::Solvent::DMSO: return "DMSO";
        case smsd::ChemOptions::Solvent::METHANOL: return "METHANOL";
        case smsd::ChemOptions::Solvent::CHLOROFORM: return "CHLOROFORM";
        case smsd::ChemOptions::Solvent::ACETONITRILE: return "ACETONITRILE";
        default: return "DIETHYL_ETHER";
    }
}

/** Every field of a profile as the fixture names it. */
std::map<std::string, std::string> profileFields(const smsd::ChemOptions& options) {
    auto flag = [](bool value) { return std::string(value ? "true" : "false"); };
    std::ostringstream pH;
    pH << options.pH;
    return {{"match_atom_type", flag(options.matchAtomType)}, {"match_formal_charge", flag(options.matchFormalCharge)},
            {"use_chirality", flag(options.useChirality)}, {"use_bond_stereo", flag(options.useBondStereo)},
            {"ring_matches_ring_only", flag(options.ringMatchesRingOnly)},
            {"complete_rings_only", flag(options.completeRingsOnly)}, {"match_isotope", flag(options.matchIsotope)},
            {"ring_fusion_mode", enumName(options.ringFusionMode)}, {"tautomer_aware", flag(options.tautomerAware)},
            {"pH", pH.str()}, {"solvent", enumName(options.solvent)}, {"match_bond_order", enumName(options.matchBondOrder)},
            {"aromaticity_mode", enumName(options.aromaticityMode)}};
}

std::string fieldText(const Json& value) {
    if (value.kind == Json::Kind::Bool) return value.boolean ? "true" : "false";
    if (value.kind == Json::Kind::String) return value.text;
    std::ostringstream text;
    text << value.number;
    return text.str();
}

void chemistryProfiles(const Json& section) {
    for (const auto& item : section["cases"].items) {
        const std::string name = "chemistry profile " + item["name"].text;
        std::map<std::string, std::string> expected;
        for (const auto& field : section["defaults"].members) expected[field.first] = fieldText(field.second);
        for (const auto& field : item["changes"].members) expected[field.first] = fieldText(field.second);
        const smsd::ChemOptions options = item["profile"].text == "tautomer_profile"
            ? smsd::ChemOptions::tautomerProfile() : smsd::ChemOptions::profile(item["profile"].text);
        const auto fields = profileFields(options);
        check(fields.size() == expected.size(), name + ": the fixture names every field");
        for (const auto& field : expected)
            check(fields.at(field.first) == field.second, name + " " + field.first + ": got " + fields.at(field.first));
    }
}

/** The chemistry of a fixture: the defaults with the listed fields changed. */
smsd::ChemOptions chemistryOf(const Json& fields) {
    smsd::ChemOptions chemistry;
    if (fields.has("match_atom_type")) chemistry.matchAtomType = fields["match_atom_type"].boolean;
    if (fields.has("match_formal_charge")) chemistry.matchFormalCharge = fields["match_formal_charge"].boolean;
    if (fields.has("ring_matches_ring_only")) chemistry.ringMatchesRingOnly = fields["ring_matches_ring_only"].boolean;
    if (fields.has("tautomer_aware")) chemistry.tautomerAware = fields["tautomer_aware"].boolean;
    if (fields.has("aromaticity_mode"))
        chemistry.aromaticityMode = fields["aromaticity_mode"].text == "STRICT" ? smsd::ChemOptions::AromaticityMode::STRICT
                                                                               : smsd::ChemOptions::AromaticityMode::FLEXIBLE;
    return chemistry;
}

void similarityUpperBound(const Json& section) {
    const double tolerance = section["tolerance"].number;
    for (const auto& item : section["cases"].items) {
        const std::string name = "similarity upper bound " + item["name"].text;
        const auto chemistry = chemistryOf(item["chemistry"]);
        const auto query = molecule(item["query"].text), target = molecule(item["target"].text);
        near(smsd::similarityUpperBound(query, target, chemistry), item["bound"].number, tolerance, name);
        near(smsd::similarityUpperBound(target, query, chemistry), item["bound"].number, tolerance, name + " reversed");
        const std::vector<smsd::MolGraph> targets{target};
        check((smsd::screenTargets(query, targets, item["bound"].number, chemistry).size() == 1),
              name + ": screening passes the pair at its own bound");
    }
}

void aromaticBondHydrogens(const Json& section) {
    for (const auto& item : section["cases"].items) {
        const auto graph = molecule(item["smiles"].text);
        check(graph.hydrogenCount == item["hydrogens"].integers(),
              "hydrogens of " + item["smiles"].text + ": got " + describe(graph.hydrogenCount));
    }
}

smsd::ChemOptions stereoChemistry(const Json& item) {
    smsd::ChemOptions chemistry;
    chemistry.useChirality = item["chirality"].boolean;
    chemistry.useBondStereo = item["bond_stereo"].boolean;
    return chemistry;
}

void stereoEmbeddings(const Json& section) {
    for (const auto& item : section["cases"].items) {
        smsd::SubstructureOptions options;
        options.chemistry = stereoChemistry(item);
        const auto query = molecule(item["query"].text), target = molecule(item["target"].text);
        const auto prepared = smsd::Substructure::find(query).withOptions(options);
        check(!prepared.match(target).matches.empty() == item["found"].boolean, "stereo embedding " + item["name"].text);
        check(!prepared.matchAll(target).matches.empty() == item["found"].boolean, "stereo embeddings " + item["name"].text);
    }
    for (const auto& item : section["common_subgraphs"].items) {
        smsd::MCSOptions options;
        options.chemistry = stereoChemistry(item);
        const auto result = smsd::MCS::find(molecule(item["query"].text)).withOptions(options).match(molecule(item["target"].text));
        const int atoms = result.matches.empty() ? 0 : static_cast<int>(result.matches.front().atomMapping.size());
        check(atoms == item["atoms"].integer(), "stereo common subgraph " + item["name"].text + ": got " + std::to_string(atoms));
    }
}

void tautomerConsistency(const Json& section) {
    for (const auto& item : section["cases"].items) {
        std::map<int, int> mapping;
        for (const auto& pair : item["atoms"].items) mapping[pair.items[0].integer()] = pair.items[1].integer();
        const bool consistent = smsd::validateTautomerConsistency(molecule(item["query"].text), molecule(item["target"].text), mapping);
        check(consistent == item["consistent"].boolean, "tautomer consistency " + item["name"].text);
    }
}

void tautomerClasses(const Json& section) {
    for (const auto& item : section["cases"].items) {
        const auto graph = molecule(item["smiles"].text);
        check(graph.tautomerClass == item["classes"].integers(),
              "tautomer classes of " + item["name"].text + ": got " + describe(graph.tautomerClass));
    }
}

void murckoScaffolds(const Json& section) {
    for (const auto& item : section["cases"].items) {
        const auto scaffold = smsd::murckoScaffold(molecule(item["smiles"].text));
        check(scaffold.hydrogenCount == item["hydrogens"].integers(),
              "scaffold hydrogens of " + item["smiles"].text + ": got " + describe(scaffold.hydrogenCount));
        check(scaffold.toCanonicalSmiles() == item["scaffold"].text,
              "scaffold of " + item["smiles"].text + ": got " + scaffold.toCanonicalSmiles());
    }
}

void consensus(const Json& section) {
    for (const auto& item : section["cases"].items) {
        const std::string name = "consensus " + item["name"].text;
        std::vector<smsd::MolGraph> molecules;
        for (const auto& smiles : item["molecules"].items) molecules.push_back(molecule(smiles.text));
        smsd::MCS::Options options;
        if (item.has("objective"))
            options.objective = item["objective"].text == "MCES" ? smsd::MCS::Objective::MCES : smsd::MCS::Objective::MCIS;
        const auto result = smsd::MCS::consensus(molecules, options, item["threshold"].number);
        check(result.found == item["found"].boolean, name + ": found");
        check(result.complete, name + ": complete");
        if (!result.found) continue;
        check(static_cast<int>(result.fragment.atomMapping.size()) == item["atoms"].integer(), name + ": atoms");
        check(static_cast<int>(result.fragment.bondMapping.size()) == item["bonds"].integer(), name + ": bonds");
        check(result.supportCount == item["support"].integer(), name + ": support");
        check(result.sourceIndex == item["source_index"].integer(), name + ": source index");
        std::vector<int> supporting;
        for (const auto& mapping : result.mappings) supporting.push_back(mapping.first);
        check(supporting == item["mappings"].integers(), name + ": supporting inputs");
        check(result.graph.n == item["atoms"].integer(), name + ": graph atoms");
    }
}

void consensusNativeOrder(const Json& section) {
    for (const auto& item : section["cases"].items) {
        const std::string name = "consensus " + item["name"].text;
        std::vector<smsd::MolGraph> molecules;
        for (const auto& smiles : item["molecules"].items) molecules.push_back(molecule(smiles.text));
        smsd::MCSOptions options;
        options.objective = item["objective"].text == "MCES" ? smsd::MCSObjective::MCES : smsd::MCSObjective::MCIS;
        options.connectivity = item["connectivity"].text == "Fragments" ? smsd::MCSConnectivity::Fragments
                                                                          : smsd::MCSConnectivity::Connected;
        const auto result = smsd::MCS::consensus(molecules, options, item["threshold"].number);
        check(result.found == item["found"].boolean && result.sourceIndex == item["source_index"].integer(),
              name + ": found " + std::to_string(result.found) + " from input " + std::to_string(result.sourceIndex));
        std::vector<int> atoms, support;
        for (const auto& atom : result.fragment.atomMapping) atoms.push_back(atom.first);
        for (const auto& mapping : result.mappings) support.push_back(mapping.first);
        check(atoms == item["atoms"].integers(), name + " atoms " + describe(atoms));
        check(static_cast<int>(result.fragment.bondMapping.size()) == item["bonds"].integer(),
              name + " bonds " + std::to_string(result.fragment.bondMapping.size()));
        check(support == item["support"].integers() && result.supportCount == static_cast<int>(support.size()),
              name + " support " + describe(support));
        check(result.toSmiles() == item["smiles"].text, name + " SMILES " + result.toSmiles());
    }
}

void strictRingFusion(const Json& section) {
    for (const auto& item : section["ring_counts"].items) {
        const auto graph = molecule(item["smiles"].text);
        graph.ensureRingCounts();
        check(graph.ringCount == item["counts"].integers(),
              "ring counts of " + item["name"].text + ": got " + describe(graph.ringCount));
    }
    for (const auto& item : section["common_bonds"].items) {
        smsd::MCSOptions options;
        options.objective = smsd::MCSObjective::MCES;
        options.chemistry.ringFusionMode = smsd::ChemOptions::RingFusionMode::STRICT;
        const auto query = molecule(item["query"].text), target = molecule(item["target"].text);
        for (int turn = 0; turn < 2; ++turn) {
            const auto result = turn == 0 ? smsd::MCS::find(query).withOptions(options).match(target)
                                          : smsd::MCS::find(target).withOptions(options).match(query);
            const int bonds = result.matches.empty() ? 0 : static_cast<int>(result.matches.front().bondMapping.size());
            const int atoms = result.matches.empty() ? 0 : static_cast<int>(result.matches.front().atomMapping.size());
            check(result.optimal && bonds == item["bonds"].integer() && atoms == item["atoms"].integer(),
                  "strict ring fusion " + item["name"].text + (turn == 0 ? "" : " reversed") + ": got "
                      + std::to_string(bonds) + " bonds on " + std::to_string(atoms) + " atoms");
        }
    }
}

} // namespace

int main() {
    std::ifstream input(SMSD_CHEMISTRY_RULES_JSON);
    if (!input) {
        std::cerr << "FAIL cannot read " << SMSD_CHEMISTRY_RULES_JSON << '\n';
        return 1;
    }
    std::stringstream buffer;
    buffer << input.rdbuf();
    const std::string source = buffer.str();
    try {
        const Json rules = JsonReader(source).value();
        similarity(rules["similarity"]);
        circularSimilarity(rules["circular_similarity"]);
        implicitHydrogens(rules["implicit_hydrogens"]);
        smarts(rules["smarts"], false);
        smarts(rules["smarts_without_hydrogen_counts"], true);
        tautomerConfidence(rules["tautomer_confidence"]);
        pathFingerprints(rules["path_fingerprints"]);
        tetrahedralCentres(rules["tetrahedral_centres"]);
        daylightSmarts(rules["smarts_daylight"]);
        mancudeLabels(rules["cip_mancude"]);
        chemistryProfiles(rules["chemistry_profiles"]);
        similarityUpperBound(rules["similarity_upper_bound"]);
        aromaticBondHydrogens(rules["aromatic_bond_hydrogens"]);
        stereoEmbeddings(rules["stereo_embeddings"]);
        tautomerConsistency(rules["tautomer_consistency"]);
        tautomerClasses(rules["tautomer_classes"]);
        murckoScaffolds(rules["murcko_scaffolds"]);
        consensusNativeOrder(rules["consensus_native_order"]);
        consensus(rules["consensus"]);
        auxiliaryDescriptors(rules["cip_auxiliary_descriptors"]);
        fingerprintOptions(rules["fingerprint_options"]);
        strictRingFusion(rules["strict_ring_fusion"]);
    } catch (const std::exception& error) {
        std::cerr << "FAIL " << error.what() << '\n';
        return 1;
    }
    if (failures != 0) return 1;
    std::cout << "PASS shared chemistry rules\n";
    return 0;
}
