/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms. *
 * Python bindings for SMSD via pybind11.
 * Usage: import smsd; result = smsd.find_mcs(mol1, mol2)
 */

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <pybind11/functional.h>
#include "smsd/smsd.hpp"
#include "smsd/rascal.hpp"
#include "smsd/smiles_parser.hpp"
#include "smsd/smarts_parser.hpp"
#include "smsd/cip.hpp"
#include "smsd/layout.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <functional>
#include <map>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#if defined(_MSC_VER)
#include <intrin.h>
#endif

// Portable bit intrinsics (Ubuntu / macOS / Windows — no infinite recursion)
inline int smsd_popcount64(uint64_t x) {
#if defined(_MSC_VER)
    return static_cast<int>(__popcnt64(x));
#elif defined(__GNUC__) || defined(__clang__)
    return __builtin_popcountll(x);   // was: smsd_popcount64(x) — infinite recursion fix
#else
    x = x - ((x >> 1) & 0x5555555555555555ULL);
    x = (x & 0x3333333333333333ULL) + ((x >> 2) & 0x3333333333333333ULL);
    return static_cast<int>(((x + (x >> 4)) & 0x0F0F0F0F0F0F0F0FULL) * 0x0101010101010101ULL >> 56);
#endif
}

inline int smsd_ctz64(uint64_t x) {
#if defined(_MSC_VER)
    unsigned long index;
    if (_BitScanForward64(&index, x)) return static_cast<int>(index);
    return 64;
#elif defined(__GNUC__) || defined(__clang__)
    return __builtin_ctzll(x);
#else
    if (x == 0) return 64;
    int c = 0;
    if (!(x & 0xFFFFFFFFULL)) { c += 32; x >>= 32; }
    if (!(x & 0xFFFFULL)) { c += 16; x >>= 16; }
    if (!(x & 0xFFULL)) { c += 8; x >>= 8; }
    if (!(x & 0xFULL)) { c += 4; x >>= 4; }
    if (!(x & 0x3ULL)) { c += 2; x >>= 2; }
    if (!(x & 0x1ULL)) c += 1;
    return c;
#endif
}

namespace py = pybind11;

// ============================================================================
// Fingerprint implementation (ported from Java SearchEngine)
// ============================================================================

namespace {

// Hash a path of atomic numbers (path fingerprint, simpler variant).
// Computes forward and reverse hashes and takes the minimum for
// canonical (direction-independent) ordering.
int fpHashPath(const int* atomicNum, const int* path, int len) {
    int fwd = 17, rev = 17;
    for (int i = 0; i < len; ++i) {
        fwd = fwd * 31 + atomicNum[path[i]];
        rev = rev * 31 + atomicNum[path[len - 1 - i]];
    }
    return std::min(fwd & 0x7FFFFFFF, rev & 0x7FFFFFFF);
}

void fpSetBit(std::vector<uint64_t>& fp, int hash, int fpSize) {
    int bit = hash % fpSize;
    fp[bit >> 6] |= uint64_t(1) << (bit & 63);
}

void fpEnumeratePaths(const std::vector<std::vector<int>>& adj,
                      const std::vector<int>& atomicNum,
                      std::vector<uint64_t>& fp, int fpSize,
                      std::vector<bool>& visited, int* path,
                      int depth, int maxDepth) {
    int cur = path[depth - 1];
    for (int nb : adj[cur]) {
        if (visited[nb]) continue;
        path[depth] = nb;
        visited[nb] = true;
        fpSetBit(fp, fpHashPath(atomicNum.data(), path, depth + 1), fpSize);
        if (depth < maxDepth)
            fpEnumeratePaths(adj, atomicNum, fp, fpSize, visited, path,
                             depth + 1, maxDepth);
        visited[nb] = false;
    }
}

// MCS-aware atom hash: encodes element, ring, aromatic, tautomer class, degree.
int mcsAtomHash(const smsd::MolGraph& g, bool tautAware, int atom) {
    int h = 17;
    int label = g.atomicNum[atom];
    bool hasTautClass = !g.tautomerClass.empty() && g.tautomerClass[atom] >= 0;
    if (tautAware && hasTautClass) label = 999; // tautomer-invariant label
    h = h * 37 + label;
    h = h * 37 + (g.ring[atom] ? 1 : 0);
    h = h * 37 + (g.aromatic[atom] ? 1 : 0);
    int tc = (!g.tautomerClass.empty()) ? g.tautomerClass[atom] : -1;
    h = h * 37 + (tautAware ? 0 : tc);
    h = h * 37 + g.degree[atom];
    return h;
}

int mcsBondHash(const smsd::MolGraph& g, int from, int to) {
    int h = 17;
    h = h * 37 + g.bondOrder(from, to);
    h = h * 37 + (g.bondInRing(from, to) ? 1 : 0);
    h = h * 37 + (g.bondAromatic(from, to) ? 1 : 0);
    return h;
}

int mcsFpHashPath(const smsd::MolGraph& g, bool tautAware,
                  const int* path, int len) {
    int fwd = 17, rev = 17;
    for (int i = 0; i < len; ++i) {
        int fi = i, ri = len - 1 - i;
        int ahF = mcsAtomHash(g, tautAware, path[fi]);
        int ahR = mcsAtomHash(g, tautAware, path[ri]);
        fwd = fwd * 31 + ahF;
        rev = rev * 31 + ahR;
        if (i < len - 1) {
            int bondF = mcsBondHash(g, path[fi], path[fi + 1]);
            int bondR = mcsBondHash(g, path[ri], path[ri - 1]);
            fwd = fwd * 31 + bondF;
            rev = rev * 31 + bondR;
        }
    }
    return std::min(fwd & 0x7FFFFFFF, rev & 0x7FFFFFFF);
}

void mcsFpEnumeratePaths(const smsd::MolGraph& g, bool tautAware,
                         const std::vector<bool>& isHeavy,
                         std::vector<uint64_t>& fp, int fpSize,
                         std::vector<bool>& visited, int* path,
                         int depth, int maxDepth) {
    int cur = path[depth - 1];
    for (int nb : g.neighbors[cur]) {
        if (visited[nb] || !isHeavy[nb]) continue;
        path[depth] = nb;
        visited[nb] = true;
        fpSetBit(fp, mcsFpHashPath(g, tautAware, path, depth + 1), fpSize);
        if (depth < maxDepth)
            mcsFpEnumeratePaths(g, tautAware, isHeavy, fp, fpSize, visited,
                                path, depth + 1, maxDepth);
        visited[nb] = false;
    }
}

int popcount64(uint64_t x) {
#if defined(__GNUC__) || defined(__clang__)
    return smsd_popcount64(x);
#else
    x = x - ((x >> 1) & 0x5555555555555555ULL);
    x = (x & 0x3333333333333333ULL) + ((x >> 2) & 0x3333333333333333ULL);
    return (int)(((x + (x >> 4)) & 0x0F0F0F0F0F0F0F0FULL) * 0x0101010101010101ULL >> 56);
#endif
}

} // anonymous namespace

// ============================================================================
// Public fingerprint functions
// ============================================================================

// Path fingerprint: enumerates all simple paths up to pathLength bonds,
// hashes each path by atomic number sequence. Returns bit positions that
// are set.
static std::vector<int> pathFingerprint(const smsd::MolGraph& mol,
                                         int pathLength = 7,
                                         int fpSize = 2048) {
    int words = (fpSize + 63) / 64;
    std::vector<uint64_t> fp(words, 0);
    int n = mol.n;
    if (n == 0) return {};

    // Build heavy-atom adjacency (skip explicit H)
    std::vector<int> heavyMap(n, -1);
    int heavyCount = 0;
    for (int i = 0; i < n; ++i) {
        if (mol.atomicNum[i] != 1) { heavyMap[i] = heavyCount++; }
    }
    if (heavyCount == 0) return {};

    std::vector<int> heavyAtomicNum(heavyCount);
    std::vector<std::vector<int>> heavyAdj(heavyCount);
    for (int i = 0; i < n; ++i) {
        if (heavyMap[i] < 0) continue;
        int hi = heavyMap[i];
        heavyAtomicNum[hi] = mol.atomicNum[i];
        for (int nb : mol.neighbors[i]) {
            if (heavyMap[nb] >= 0) heavyAdj[hi].push_back(heavyMap[nb]);
        }
    }

    std::vector<int> pathBuf(pathLength + 1);
    std::vector<bool> visited(heavyCount, false);
    for (int start = 0; start < heavyCount; ++start) {
        pathBuf[0] = start;
        visited[start] = true;
        fpSetBit(fp, fpHashPath(heavyAtomicNum.data(), pathBuf.data(), 1), fpSize);
        fpEnumeratePaths(heavyAdj, heavyAtomicNum, fp, fpSize, visited,
                         pathBuf.data(), 1, pathLength);
        visited[start] = false;
    }

    // Collect set bit positions
    std::vector<int> bits;
    for (int w = 0; w < words; ++w) {
        uint64_t word = fp[w];
        while (word) {
            int bit = w * 64 + smsd_ctz64(word);
            if (bit < fpSize) bits.push_back(bit);
            word &= word - 1;
        }
    }
    return bits;
}

// MCS-aware fingerprint: encodes element, ring, aromatic, tautomer class,
// degree per atom; order, ring, aromatic per bond. Returns set bit positions.
static std::vector<int> mcsFingerprint(const smsd::MolGraph& mol,
                                        int pathLength = 7,
                                        int fpSize = 2048) {
    int words = (fpSize + 63) / 64;
    std::vector<uint64_t> fp(words, 0);
    int n = mol.n;
    if (n == 0) return {};

    std::vector<bool> isHeavy(n);
    int heavyCount = 0;
    for (int i = 0; i < n; ++i) {
        isHeavy[i] = (mol.atomicNum[i] != 1);
        if (isHeavy[i]) heavyCount++;
    }
    if (heavyCount == 0) return {};

    bool tautAware = !mol.tautomerClass.empty();
    std::vector<int> pathBuf(pathLength + 1);
    std::vector<bool> visited(n, false);
    for (int start = 0; start < n; ++start) {
        if (!isHeavy[start]) continue;
        pathBuf[0] = start;
        visited[start] = true;
        fpSetBit(fp, mcsFpHashPath(mol, tautAware, pathBuf.data(), 1), fpSize);
        mcsFpEnumeratePaths(mol, tautAware, isHeavy, fp, fpSize, visited,
                            pathBuf.data(), 1, pathLength);
        visited[start] = false;
    }

    std::vector<int> bits;
    for (int w = 0; w < words; ++w) {
        uint64_t word = fp[w];
        while (word) {
            int bit = w * 64 + smsd_ctz64(word);
            if (bit < fpSize) bits.push_back(bit);
            word &= word - 1;
        }
    }
    return bits;
}

// Check if fingerprint fp1 is a subset of fp2 (all bits in fp1 are set in fp2).
// Both are represented as sorted lists of set bit positions.
static bool fingerprintSubset(const std::vector<int>& fp1,
                               const std::vector<int>& fp2) {
    size_t j = 0;
    for (int bit : fp1) {
        while (j < fp2.size() && fp2[j] < bit) ++j;
        if (j >= fp2.size() || fp2[j] != bit) return false;
    }
    return true;
}

// Analyze fingerprint quality: density, number of set bits, expected
// collision rate.
static py::dict analyzeFpQuality(const std::vector<int>& fp, int fpSize = 2048) {
    int setBits = static_cast<int>(fp.size());
    double density = static_cast<double>(setBits) / fpSize;
    // Expected collision probability for random hashing
    double expectedCollisions = 1.0 - std::exp(-static_cast<double>(setBits) / fpSize);
    py::dict result;
    result["set_bits"] = setBits;
    result["fp_size"] = fpSize;
    result["density"] = density;
    result["expected_collision_rate"] = expectedCollisions;
    result["is_saturated"] = density > 0.75;
    return result;
}

// ============================================================================
// Module definition
// ============================================================================

PYBIND11_MODULE(_smsd, m) {
    m.doc() = "SMSD -- Substructure & MCS search for chemical graphs";

    // -----------------------------------------------------------------------
    // Enums
    // -----------------------------------------------------------------------
    py::enum_<smsd::ChemOptions::BondOrderMode>(m, "BondOrderMode")
        .value("STRICT", smsd::ChemOptions::BondOrderMode::STRICT)
        .value("LOOSE",  smsd::ChemOptions::BondOrderMode::LOOSE)
        .value("ANY",    smsd::ChemOptions::BondOrderMode::ANY)
        .export_values();

    py::enum_<smsd::ChemOptions::AromaticityMode>(m, "AromaticityMode")
        .value("STRICT",   smsd::ChemOptions::AromaticityMode::STRICT)
        .value("FLEXIBLE", smsd::ChemOptions::AromaticityMode::FLEXIBLE)
        .export_values();

    py::enum_<smsd::ChemOptions::RingFusionMode>(m, "RingFusionMode")
        .value("IGNORE",     smsd::ChemOptions::RingFusionMode::IGNORE)
        .value("PERMISSIVE", smsd::ChemOptions::RingFusionMode::PERMISSIVE)
        .value("STRICT",     smsd::ChemOptions::RingFusionMode::STRICT)
        .export_values();

    py::enum_<smsd::ChemOptions::Solvent>(m, "Solvent")
        .value("AQUEOUS",       smsd::ChemOptions::Solvent::AQUEOUS)
        .value("DMSO",          smsd::ChemOptions::Solvent::DMSO)
        .value("METHANOL",      smsd::ChemOptions::Solvent::METHANOL)
        .value("CHLOROFORM",    smsd::ChemOptions::Solvent::CHLOROFORM)
        .value("ACETONITRILE",  smsd::ChemOptions::Solvent::ACETONITRILE)
        .value("DIETHYL_ETHER", smsd::ChemOptions::Solvent::DIETHYL_ETHER)
        .export_values();

    // -----------------------------------------------------------------------
    // ChemOptions
    // -----------------------------------------------------------------------
    py::class_<smsd::ChemOptions>(m, "ChemOptions")
        .def(py::init<>())
        // Atom matching
        .def_readwrite("match_atom_type",       &smsd::ChemOptions::matchAtomType)
        .def_readwrite("match_formal_charge",   &smsd::ChemOptions::matchFormalCharge)
        .def_readwrite("match_isotope",         &smsd::ChemOptions::matchIsotope)
        .def_readwrite("use_chirality",         &smsd::ChemOptions::useChirality)
        .def_readwrite("use_bond_stereo",       &smsd::ChemOptions::useBondStereo)
        .def_readwrite("ring_matches_ring_only",&smsd::ChemOptions::ringMatchesRingOnly)
        // Ring completeness
        .def_readwrite("complete_rings_only",   &smsd::ChemOptions::completeRingsOnly)
        // Tautomer
        .def_readwrite("tautomer_aware",        &smsd::ChemOptions::tautomerAware)
        // Solvent
        .def_readwrite("solvent",               &smsd::ChemOptions::solvent)
        .def("with_solvent", &smsd::ChemOptions::withSolvent,
             py::arg("solvent"),
             "Set solvent for tautomer weight corrections (fluent)")
        // Bond matching
        .def_readwrite("match_bond_order",      &smsd::ChemOptions::matchBondOrder)
        .def_readwrite("aromaticity_mode",      &smsd::ChemOptions::aromaticityMode)
        // Ring fusion
        .def_readwrite("ring_fusion_mode",      &smsd::ChemOptions::ringFusionMode)
        // Named profiles
        .def_static("tautomer_profile", &smsd::ChemOptions::tautomerProfile,
                     "Create ChemOptions configured for tautomer-aware matching")
        .def_static("profile", &smsd::ChemOptions::profile,
                     py::arg("name") = "default",
                     "Create ChemOptions from named profile ('strict', 'default')")
        .def("__repr__", [](const smsd::ChemOptions& o) {
            return "<ChemOptions match_atom_type=" +
                   std::string(o.matchAtomType ? "True" : "False") +
                   " tautomer_aware=" +
                   std::string(o.tautomerAware ? "True" : "False") + ">";
        });

    // -----------------------------------------------------------------------
    // McsOptions
    // -----------------------------------------------------------------------
    py::class_<smsd::McsOptions>(m, "McsOptions")
        .def(py::init<>())
        .def_readwrite("induced",           &smsd::McsOptions::induced)
        .def_readwrite("connected_only",    &smsd::McsOptions::connectedOnly)
        .def_readwrite("disconnected_mcs",  &smsd::McsOptions::disconnectedMCS)
        .def_readwrite("maximize_bonds",    &smsd::McsOptions::maximizeBonds)
        .def_readwrite("min_fragment_size", &smsd::McsOptions::minFragmentSize)
        .def_readwrite("max_fragments",     &smsd::McsOptions::maxFragments)
        .def_readwrite("timeout_ms",        &smsd::McsOptions::timeoutMs)
        .def_readwrite("extra_seeds",       &smsd::McsOptions::extraSeeds)
        .def_readwrite("template_fuzzy_atoms", &smsd::McsOptions::templateFuzzyAtoms)
        .def_readwrite("reaction_aware",     &smsd::McsOptions::reactionAware)
        .def_readwrite("near_mcs_delta",     &smsd::McsOptions::nearMcsDelta)
        .def_readwrite("near_mcs_candidates", &smsd::McsOptions::nearMcsCandidates)
        .def_readwrite("bond_change_aware",  &smsd::McsOptions::bondChangeAware)
        .def("__repr__", [](const smsd::McsOptions& o) {
            return "<McsOptions induced=" +
                   std::string(o.induced ? "True" : "False") +
                   " connected_only=" +
                   std::string(o.connectedOnly ? "True" : "False") +
                   " timeout_ms=" + std::to_string(o.timeoutMs) + ">";
        });

    // -----------------------------------------------------------------------
    // MolGraph
    // -----------------------------------------------------------------------
    py::class_<smsd::MolGraph>(m, "MolGraph")
        .def(py::init<>())
        .def_readonly("n",          &smsd::MolGraph::n)
        .def_readonly("atomic_num", &smsd::MolGraph::atomicNum)
        .def_readonly("degree",     &smsd::MolGraph::degree)
        .def_readonly("aromatic",   &smsd::MolGraph::aromatic)
        .def_readonly("ring",       &smsd::MolGraph::ring)
        .def_readonly("hydrogen_count", &smsd::MolGraph::hydrogenCount)
        .def_readonly("atom_class", &smsd::MolGraph::atomClass)
        .def_readonly("formal_charge", &smsd::MolGraph::formalCharge)
        .def_readonly("mass_number",   &smsd::MolGraph::massNumber)
        .def_readwrite("name",         &smsd::MolGraph::name)
        .def_readwrite("program_line", &smsd::MolGraph::programLine)
        .def_readwrite("comment",      &smsd::MolGraph::comment)
        .def_readwrite("properties",   &smsd::MolGraph::properties)
        .def_readonly("label",         &smsd::MolGraph::label)
        .def_readonly("tautomer_class",&smsd::MolGraph::tautomerClass)
        .def_readonly("morgan_rank",   &smsd::MolGraph::morganRank)
        .def_readonly("canonical_label", &smsd::MolGraph::canonicalLabel)
        .def_readonly("orbit",         &smsd::MolGraph::orbit)
        .def("bond_order", [](const smsd::MolGraph& g, int i, int j) {
                 if (i < 0 || i >= g.n || j < 0 || j >= g.n)
                     throw py::index_error("atom index out of range");
                 return g.bondOrder(i, j);
             }, py::arg("i"), py::arg("j"),
             "Get bond order between atoms i and j (0 if no bond)")
        .def("bond_in_ring", [](const smsd::MolGraph& g, int i, int j) {
                 if (i < 0 || i >= g.n || j < 0 || j >= g.n)
                     throw py::index_error("atom index out of range");
                 return g.bondInRing(i, j);
             }, py::arg("i"), py::arg("j"))
        .def("bond_aromatic", [](const smsd::MolGraph& g, int i, int j) {
                 if (i < 0 || i >= g.n || j < 0 || j >= g.n)
                     throw py::index_error("atom index out of range");
                 return g.bondAromatic(i, j);
             }, py::arg("i"), py::arg("j"))
        .def("ensure_canonical", [](const smsd::MolGraph& g) {
                 g.ensureCanonical();
             }, "Compute and cache Morgan ranks, canonical labels, and orbits")
        .def("ensure_ring_counts", [](const smsd::MolGraph& g) {
                 g.ensureRingCounts();
             }, "Compute and cache per-atom ring counts")
        .def("prewarm", [](const smsd::MolGraph& g) {
                 g.ensureCanonical();
                 g.ensureRingCounts();
                 g.getNLF1();
                 g.getNLF2();
                 g.getNLF3();
                 g.getNeighborsByDegDesc();
             },
             py::call_guard<py::gil_scoped_release>(),
             "Compute and cache common lazy invariants used by matching and batch APIs")
        .def("has_bond", [](const smsd::MolGraph& g, int i, int j) {
                 if (i < 0 || i >= g.n || j < 0 || j >= g.n)
                     throw py::index_error("atom index out of range");
                 return g.hasBond(i, j);
             }, py::arg("i"), py::arg("j"))
        .def("__repr__", [](const smsd::MolGraph& g) {
            return "<MolGraph n=" + std::to_string(g.n) + ">";
        })
        .def("__len__", [](const smsd::MolGraph& g) { return g.n; });

    // -----------------------------------------------------------------------
    // MolGraph.Builder
    // -----------------------------------------------------------------------
    py::class_<smsd::MolGraph::Builder>(m, "MolGraphBuilder")
        .def(py::init<>())
        .def("atom_count", [](smsd::MolGraph::Builder& b, int n) -> smsd::MolGraph::Builder& {
            return b.atomCount(n);
        }, py::arg("n"), py::return_value_policy::reference_internal,
             "Set the number of atoms")
        .def("atomic_numbers", [](smsd::MolGraph::Builder& b, std::vector<int> v) -> smsd::MolGraph::Builder& {
            return b.atomicNumbers(std::move(v));
        }, py::arg("nums"), py::return_value_policy::reference_internal,
             "Set atomic numbers for all atoms")
        .def("formal_charges", [](smsd::MolGraph::Builder& b, std::vector<int> v) -> smsd::MolGraph::Builder& {
            return b.formalCharges(std::move(v));
        }, py::arg("charges"), py::return_value_policy::reference_internal,
             "Set formal charges")
        .def("ring_flags", [](smsd::MolGraph::Builder& b, std::vector<uint8_t> v) -> smsd::MolGraph::Builder& {
            return b.ringFlags(std::move(v));
        }, py::arg("flags"), py::return_value_policy::reference_internal,
             "Set ring membership flags")
        .def("aromatic_flags", [](smsd::MolGraph::Builder& b, std::vector<uint8_t> v) -> smsd::MolGraph::Builder& {
            return b.aromaticFlags(std::move(v));
        }, py::arg("flags"), py::return_value_policy::reference_internal,
             "Set aromaticity flags")
        .def("neighbors", [](smsd::MolGraph::Builder& b, std::vector<std::vector<int>> v) -> smsd::MolGraph::Builder& {
            return b.setNeighbors(std::move(v));
        }, py::arg("adj"), py::return_value_policy::reference_internal,
             "Set adjacency lists")
        .def("bond_orders", [](smsd::MolGraph::Builder& b, std::vector<std::vector<int>> v) -> smsd::MolGraph::Builder& {
            return b.setBondOrders(std::move(v));
        }, py::arg("orders"), py::return_value_policy::reference_internal,
             "Set bond orders (parallel to neighbors)")
        .def("atom_ids", [](smsd::MolGraph::Builder& b, std::vector<int> v) -> smsd::MolGraph::Builder& {
            return b.atomIds(std::move(v));
        }, py::arg("ids"), py::return_value_policy::reference_internal,
             "Set optional external atom IDs (1-based, gaps allowed) for index translation")
        .def("name", [](smsd::MolGraph::Builder& b, std::string v) -> smsd::MolGraph::Builder& {
            return b.name(std::move(v));
        }, py::arg("value"), py::return_value_policy::reference_internal,
             "Set molecule name/title")
        .def("program_line", [](smsd::MolGraph::Builder& b, std::string v) -> smsd::MolGraph::Builder& {
            return b.programLine(std::move(v));
        }, py::arg("value"), py::return_value_policy::reference_internal,
             "Set MDL program line/header text")
        .def("comment", [](smsd::MolGraph::Builder& b, std::string v) -> smsd::MolGraph::Builder& {
            return b.comment(std::move(v));
        }, py::arg("value"), py::return_value_policy::reference_internal,
             "Set molecule comment text")
        .def("properties", [](smsd::MolGraph::Builder& b, std::map<std::string, std::string> v) -> smsd::MolGraph::Builder& {
            return b.properties(std::move(v));
        }, py::arg("values"), py::return_value_policy::reference_internal,
             "Set SDF-style molecule properties")
        .def("build", &smsd::MolGraph::Builder::build,
             "Build the immutable MolGraph");

    // -----------------------------------------------------------------------
    // ParseOptions (for lenient SMILES parsing)
    // -----------------------------------------------------------------------
    py::class_<smsd::ParseOptions>(m, "ParseOptions")
        .def(py::init<>())
        .def_readwrite("lenient", &smsd::ParseOptions::lenient,
                       "If true, recover from malformed SMILES instead of throwing");

    // -----------------------------------------------------------------------
    // SMILES parser / writer
    // -----------------------------------------------------------------------
    m.def("parse_smiles",
          static_cast<smsd::MolGraph (*)(const std::string&)>(&smsd::parseSMILES),
          py::arg("smiles"),
          "Parse a SMILES string into a MolGraph");
    m.def("parse_smiles_lenient",
          static_cast<smsd::MolGraph (*)(const std::string&, const smsd::ParseOptions&)>(&smsd::parseSMILES),
          py::arg("smiles"), py::arg("opts"),
          "Parse a SMILES string with options (e.g. lenient mode)");

    m.def("to_smiles", &smsd::toSMILES,
          py::arg("mol"),
          "Generate canonical SMILES from a MolGraph");

    m.def("read_mol_block", &smsd::readMolBlock,
          py::arg("mol_block"),
          "Parse an MDL MOL V2000 block into a MolGraph");
    m.def("write_mol_block", &smsd::writeMolBlock,
          py::arg("mol"),
          "Write a MolGraph as an MDL MOL V2000 block");
    m.def("write_mol_block_v3000", &smsd::writeMolBlockV3000,
          py::arg("mol"),
          "Write a MolGraph as an MDL MOL V3000 block");
    m.def("write_sdf_record", &smsd::writeSDFRecord,
          py::arg("mol"),
          "Write a MolGraph as a single SDF record");

    // -----------------------------------------------------------------------
    // Substructure search
    // -----------------------------------------------------------------------
    m.def("is_substructure", &smsd::isSubstructure,
          py::arg("query"), py::arg("target"),
          py::arg("opts") = smsd::ChemOptions(),
          py::arg("timeout_ms") = 10000,
          py::call_guard<py::gil_scoped_release>(),
          "Check if query is a substructure of target");

    m.def("find_substructure", &smsd::findSubstructure,
          py::arg("query"), py::arg("target"),
          py::arg("opts") = smsd::ChemOptions(),
          py::arg("timeout_ms") = 10000,
          py::call_guard<py::gil_scoped_release>(),
          "Find one substructure mapping (query->target atom pairs)");

    // -----------------------------------------------------------------------
    // MCS search
    // -----------------------------------------------------------------------
    m.def("find_mcs",
          [](const smsd::MolGraph& g1, const smsd::MolGraph& g2,
             const smsd::ChemOptions& chem, smsd::McsOptions opts,
             int64_t timeout_ms) {
              if (timeout_ms > 0) opts.timeoutMs = timeout_ms;
              return smsd::findMCS(g1, g2, chem, opts);
          },
          py::arg("g1"), py::arg("g2"),
          py::arg("chem") = smsd::ChemOptions(),
          py::arg("opts") = smsd::McsOptions(),
          py::arg("timeout_ms") = -1,
          py::call_guard<py::gil_scoped_release>(),
          "Find Maximum Common Substructure. timeout_ms overrides opts.timeout_ms if > 0.");

    m.def("find_all_mcs",
          [](const smsd::MolGraph& g1, const smsd::MolGraph& g2,
             const smsd::ChemOptions& chem, smsd::McsOptions opts,
             int maxResults, int64_t timeout_ms) {
              if (timeout_ms > 0) opts.timeoutMs = timeout_ms;
              auto results = smsd::findAllMCS(g1, g2, chem, opts, maxResults);
              // Convert vector<map<int,int>> to vector<dict>
              std::vector<std::map<int,int>> out;
              out.reserve(results.size());
              for (auto& m : results) out.push_back(std::move(m));
              return out;
          },
          py::arg("g1"), py::arg("g2"),
          py::arg("chem") = smsd::ChemOptions(),
          py::arg("opts") = smsd::McsOptions(),
          py::arg("max_results") = 10,
          py::arg("timeout_ms") = -1,
          py::call_guard<py::gil_scoped_release>(),
          "Find multiple distinct MCS mappings of the maximum size.\n"
          "Returns a list of dicts, each mapping query atom index to target atom index.\n"
          "max_results caps the number of solutions (default 10).\n"
          "timeout_ms overrides opts.timeout_ms if > 0.");

    m.def("mcs_size",
          [](const smsd::MolGraph& g1, const smsd::MolGraph& g2,
             const smsd::ChemOptions& chem, smsd::McsOptions opts,
             int64_t timeout_ms) {
              if (timeout_ms > 0) opts.timeoutMs = timeout_ms;
              return static_cast<int>(smsd::findMCS(g1, g2, chem, opts).size());
          },
          py::arg("g1"), py::arg("g2"),
          py::arg("chem") = smsd::ChemOptions(),
          py::arg("opts") = smsd::McsOptions(),
          py::arg("timeout_ms") = 500,
          py::call_guard<py::gil_scoped_release>(),
          "Return MCS size only (faster for screening). Default timeout 500ms.");

    m.def("find_mcs_progressive",
          [](const smsd::MolGraph& g1, const smsd::MolGraph& g2,
             const smsd::ChemOptions& chem, smsd::McsOptions opts,
             py::object callback, int64_t timeout_ms) {
              if (timeout_ms > 0) opts.timeoutMs = timeout_ms;
              // Check is_none() BEFORE releasing GIL (py::object access needs GIL)
              bool hasCallback = !callback.is_none();
              if (!hasCallback) {
                  py::gil_scoped_release release;
                  return smsd::findMCS(g1, g2, chem, opts);
              }
              // Capture callback BY VALUE (copy py::object) to avoid dangling reference
              py::object cb = callback;
              smsd::McsProgressFn fn = [cb](
                  const std::map<int,int>& best, int bestSize, int64_t elapsedMs) {
                  py::gil_scoped_acquire acq;
                  cb(best, bestSize, elapsedMs);
              };
              {
                  py::gil_scoped_release release;
                  return smsd::findMCS(g1, g2, chem, opts, fn);
              }
          },
          py::arg("g1"), py::arg("g2"),
          py::arg("chem") = smsd::ChemOptions(),
          py::arg("opts") = smsd::McsOptions(),
          py::arg("callback") = py::none(),
          py::arg("timeout_ms") = -1,
          "Find MCS with optional progress callback.\n"
          "callback(best_mapping, best_size, elapsed_ms) is called after each pipeline level.");

    m.def("translate_to_atom_ids",
          [](const std::map<int,int>& mapping,
             const smsd::MolGraph& g1, const smsd::MolGraph& g2) {
              return smsd::translateToAtomIds(mapping, g1, g2);
          },
          py::arg("mapping"), py::arg("g1"), py::arg("g2"),
          "Translate MCS mapping from internal 0-based indices to external atom IDs.\n"
          "If neither graph has atomId set, the mapping is returned unchanged.");

    m.def("validate_tautomer_consistency", &smsd::validateTautomerConsistency,
          py::arg("g1"), py::arg("g2"), py::arg("mcs"),
          "Validate that an MCS mapping conserves mobilisable protons across "
          "tautomeric groups. Returns True if consistent, False if proton "
          "counts mismatch (impossible proton transfer required).");

    m.def("mcs_to_smiles", &smsd::mcsToSmiles,
          py::arg("g"), py::arg("mapping"),
          "Extract the MCS induced subgraph as a canonical SMILES string. "
          "Takes the query MolGraph and the MCS mapping (query->target). "
          "Returns \"\" if the mapping is empty.");

    m.def("find_mcs_smiles",
          [](const smsd::MolGraph& g1, const smsd::MolGraph& g2,
             const smsd::ChemOptions& chem, smsd::McsOptions opts,
             int64_t timeout_ms) {
              if (timeout_ms > 0) opts.timeoutMs = timeout_ms;
              return smsd::findMcsSmiles(g1, g2, chem, opts);
          },
          py::arg("g1"), py::arg("g2"),
          py::arg("chem") = smsd::ChemOptions(),
          py::arg("opts") = smsd::McsOptions(),
          py::arg("timeout_ms") = -1,
          py::call_guard<py::gil_scoped_release>(),
          "Compute MCS and return it as a canonical SMILES string. "
          "timeout_ms overrides opts.timeout_ms if > 0.");

    // -----------------------------------------------------------------------
    // Reaction-aware MCS (v6.4.0)
    // -----------------------------------------------------------------------
    m.def("reaction_aware_mcs",
          [](const smsd::MolGraph& g1, const smsd::MolGraph& g2,
             const smsd::ChemOptions& chem, smsd::McsOptions opts,
             int64_t timeout_ms) {
              if (timeout_ms > 0) opts.timeoutMs = timeout_ms;
              opts.reactionAware = true;
              return smsd::reactionAwareMCS(g1, g2, chem, opts);
          },
          py::arg("g1"), py::arg("g2"),
          py::arg("chem") = smsd::ChemOptions(),
          py::arg("opts") = smsd::McsOptions(),
          py::arg("timeout_ms") = 10000,
          py::call_guard<py::gil_scoped_release>(),
          "Reaction-aware MCS: find MCS candidates, generate near-MCS\n"
          "variants (K-1, K-2), and re-rank by heteroatom coverage,\n"
          "rare-element importance, and connectivity.\n"
          "Returns the best-scoring candidate mapping as a dict.");

    m.def("map_reaction_aware",
          [](const std::string& smi1, const std::string& smi2,
             smsd::ChemOptions chem, smsd::McsOptions opts,
             int64_t timeout_ms) {
              if (timeout_ms > 0) opts.timeoutMs = timeout_ms;
              return smsd::mapReactionAware(smi1, smi2, chem, opts);
          },
          py::arg("smi1"), py::arg("smi2"),
          py::arg("chem") = smsd::ChemOptions(),
          py::arg("opts") = smsd::McsOptions(),
          py::arg("timeout_ms") = 10000,
          py::call_guard<py::gil_scoped_release>(),
          "Convenience: parse two SMILES, run reaction-aware MCS,\n"
          "return the best mapping as a dict {reactant_idx: product_idx}.");

    // -----------------------------------------------------------------------
    // Bond-change scoring (v6.6.0)
    // -----------------------------------------------------------------------
    m.def("bond_change_score",
          [](const std::map<int,int>& mapping,
             const smsd::MolGraph& g1, const smsd::MolGraph& g2) {
              return smsd::bondChangeScore(mapping, g1, g2);
          },
          py::arg("mapping"), py::arg("g1"), py::arg("g2"),
          "Score an MCS mapping by bond-change plausibility.\n"
          "Lower score = more chemically plausible.\n"
          "C-C break=3.0, C-N/O=1.5, heteroatom=0.5.");

    // -----------------------------------------------------------------------
    // Batch MCS with non-overlap constraints (v6.6.0)
    // -----------------------------------------------------------------------
    m.def("batch_mcs_constrained",
          [](const std::vector<smsd::MolGraph>& queries,
             const std::vector<smsd::MolGraph>& targets,
             const smsd::ChemOptions& chem,
             const smsd::McsOptions& opts) {
              return smsd::batchMcsConstrained(queries, targets, chem, opts);
          },
          py::arg("queries"), py::arg("targets"),
          py::arg("chem") = smsd::ChemOptions(),
          py::arg("opts") = smsd::McsOptions(),
          py::call_guard<py::gil_scoped_release>(),
          "Find MCS across multiple molecule pairs with non-overlap\n"
          "atom exclusion constraints. Larger pairs processed first.\n"
          "Returns list of mappings (same order as input).");

    // -----------------------------------------------------------------------
    // RASCAL screening
    // -----------------------------------------------------------------------
    m.def("similarity_upper_bound", &smsd::similarityUpperBound,
          py::arg("g1"), py::arg("g2"),
          "RASCAL O(V+E) similarity upper bound");

    m.def("screen_targets", &smsd::screenTargets,
          py::arg("query"), py::arg("targets"), py::arg("threshold"),
          py::call_guard<py::gil_scoped_release>(),
          "Batch RASCAL screening: returns indices of targets above threshold");

    // -----------------------------------------------------------------------
    // Fingerprints
    // -----------------------------------------------------------------------
    m.def("path_fingerprint", &pathFingerprint,
          py::arg("mol"), py::arg("path_length") = 7, py::arg("fp_size") = 2048,
          "Compute path-based fingerprint. Returns list of set bit positions.");

    m.def("mcs_fingerprint",
          [](const smsd::MolGraph& mol, int pathLength, int fpSize) {
              auto fp = smsd::batch::detail::computeMcsFingerprint(mol, pathLength, fpSize);
              std::vector<int> bits;
              for (size_t w = 0; w < fp.size(); ++w)
                  for (int b = 0; b < 64; ++b)
                      if (fp[w] & (1ULL << b)) {
                          int pos = static_cast<int>(w * 64 + b);
                          if (pos < fpSize) bits.push_back(pos);
                      }
              return bits;
          },
          py::arg("mol"), py::arg("path_length") = 7, py::arg("fp_size") = 2048,
          "Compute MCS-aware path fingerprint. Returns list of set bit positions.");

    m.def("fingerprint_subset", &fingerprintSubset,
          py::arg("fp1"), py::arg("fp2"),
          "Check if fingerprint fp1 is a subset of fp2 (all bits in fp1 are set in fp2)");

    m.def("analyze_fp_quality", &analyzeFpQuality,
          py::arg("fp"), py::arg("fp_size") = 2048,
          "Analyze fingerprint quality: density, set bits, collision rate");

    // -----------------------------------------------------------------------
    // SMARTS writer
    // -----------------------------------------------------------------------
    m.def("to_smarts",
          [](const smsd::MolGraph& g,
             bool inclArom, bool inclCharge,
             bool inclRing, bool inclH, bool inclIso) {
              smsd::SmartsWriteOptions opts;
              opts.includeAromaticity = inclArom;
              opts.includeCharge      = inclCharge;
              opts.includeRingMember  = inclRing;
              opts.includeHCount      = inclH;
              opts.includeIsotope     = inclIso;
              return smsd::toSMARTS(g, opts);
          },
          py::arg("mol"),
          py::arg("include_aromaticity") = true,
          py::arg("include_charge")      = true,
          py::arg("include_ring_member") = false,
          py::arg("include_h_count")     = false,
          py::arg("include_isotope")     = false,
          "Generate canonical SMARTS with configurable predicates");

    // -----------------------------------------------------------------------
    // GPU / compute backend
    // -----------------------------------------------------------------------
    m.def("gpu_is_available",
          []() { return smsd::gpu::isAvailable(); },
          "True when a GPU is available at runtime (CUDA on Linux/Windows, "
          "Metal on macOS/Apple Silicon).");

    m.def("gpu_device_info",
          []() { return smsd::gpu::deviceInfo(); },
          "Human-readable compute backend string.\n"
          "Examples:\n"
          "  'GPU: Tesla T4 (compute 7.5)  [OpenMP 4.5, 8 threads]'\n"
          "  'Metal GPU: Apple M2 Pro  [OpenMP 5.0, 10 threads]'\n"
          "  'CPU: OpenMP 4.5, 16 threads'");

    // -----------------------------------------------------------------------
    // Parallel batch operations (OpenMP on CPU, CUDA on GPU when available)
    // -----------------------------------------------------------------------
    m.def("batch_substructure",
          [](const smsd::MolGraph& query,
             const std::vector<smsd::MolGraph>& targets,
             const smsd::ChemOptions& opts,
             int numThreads) {
              return smsd::batch::batchSubstructure(query, targets, opts, numThreads);
          },
          py::arg("query"), py::arg("targets"),
          py::arg("opts")        = smsd::ChemOptions(),
          py::arg("num_threads") = 0,
          py::call_guard<py::gil_scoped_release>(),
          "Parallel 1-query-vs-N substructure check. Returns list[bool]. "
          "num_threads=0 uses all available processors.");

    m.def("batch_mcs",
          [](const smsd::MolGraph& query,
             const std::vector<smsd::MolGraph>& targets,
             const smsd::ChemOptions& chem,
             const smsd::McsOptions& opts,
             int numThreads) {
              return smsd::batch::batchMCS(query, targets, chem, opts, numThreads);
          },
          py::arg("query"), py::arg("targets"),
          py::arg("chem")        = smsd::ChemOptions(),
          py::arg("opts")        = smsd::McsOptions(),
          py::arg("num_threads") = 0,
          py::call_guard<py::gil_scoped_release>(),
          "Parallel 1-query-vs-N MCS. Returns list[dict] of atom mappings. "
          "num_threads=0 uses all available processors.");

    m.def("batch_mcs_size",
          [](const smsd::MolGraph& query,
             const std::vector<smsd::MolGraph>& targets,
             const smsd::ChemOptions& chem,
             const smsd::McsOptions& opts,
             int numThreads) {
              return smsd::batch::batchMCSSize(query, targets, chem, opts, numThreads);
          },
          py::arg("query"), py::arg("targets"),
          py::arg("chem")        = smsd::ChemOptions(),
          py::arg("opts")        = smsd::McsOptions(),
          py::arg("num_threads") = 0,
          py::call_guard<py::gil_scoped_release>(),
          "Parallel 1-query-vs-N MCS size only. Returns list[int]. "
          "num_threads=0 uses all available processors.");

    // -----------------------------------------------------------------------
    // Batch screening + fingerprint
    // -----------------------------------------------------------------------
    m.def("screen_and_match",
          [](const smsd::MolGraph& query,
             const std::vector<smsd::MolGraph>& targets,
             double threshold,
             const smsd::ChemOptions& chem,
             const smsd::McsOptions& opts,
             int numThreads) {
              return smsd::batch::screenAndMatch(query, targets, chem, opts, threshold, numThreads);
          },
          py::arg("query"), py::arg("targets"), py::arg("threshold"),
          py::arg("chem")        = smsd::ChemOptions(),
          py::arg("opts")        = smsd::McsOptions(),
          py::arg("num_threads") = 0,
          py::call_guard<py::gil_scoped_release>(),
          "RASCAL pre-screen + exact MCS on hits. Returns list[(index, mapping)].");

    m.def("screen_and_mcs_size",
          [](const smsd::MolGraph& query,
             const std::vector<smsd::MolGraph>& targets,
             double threshold,
             const smsd::ChemOptions& chem,
             const smsd::McsOptions& opts,
             int numThreads) {
              return smsd::batch::screenAndMCSSize(query, targets, chem, opts, threshold, numThreads);
          },
          py::arg("query"), py::arg("targets"), py::arg("threshold"),
          py::arg("chem")        = smsd::ChemOptions(),
          py::arg("opts")        = smsd::McsOptions(),
          py::arg("num_threads") = 0,
          py::call_guard<py::gil_scoped_release>(),
          "RASCAL pre-screen + exact MCS size on hits. Returns list[(index, size)].");

    m.def("batch_fingerprint",
          [](const std::vector<smsd::MolGraph>& mols,
             int pathLength, int fpSize, int numThreads) {
              return smsd::batch::batchFingerprint(mols, pathLength, fpSize, numThreads);
          },
          py::arg("mols"), py::arg("path_length") = 7,
          py::arg("fp_size") = 2048, py::arg("num_threads") = 0,
          py::call_guard<py::gil_scoped_release>(),
          "Parallel path fingerprint generation for N molecules.");

    m.def("batch_fingerprint_screen",
          [](const std::vector<uint64_t>& queryFP,
             const std::vector<std::vector<uint64_t>>& targetFPs,
             int numThreads) {
              return smsd::batch::batchFingerprintScreen(queryFP, targetFPs, numThreads);
          },
          py::arg("query_fp"), py::arg("target_fps"),
          py::arg("num_threads") = 0,
          py::call_guard<py::gil_scoped_release>(),
          "Parallel fingerprint subset screening: returns indices where query is subset of target.");

    // -----------------------------------------------------------------------
    // Circular Fingerprints (ECFP / FCFP)
    // -----------------------------------------------------------------------
    m.def("circular_fingerprint",
          [](const smsd::MolGraph& mol, int radius, int fpSize,
             const std::string& mode) {
              auto fp = (mode == "fcfp")
                  ? smsd::batch::detail::computeCircularFingerprintFCFP(mol, radius, fpSize)
                  : smsd::batch::detail::computeCircularFingerprintECFP(mol, radius, fpSize);
              // Convert uint64_t words to list of set-bit positions for Python
              std::vector<int> bits;
              for (size_t w = 0; w < fp.size(); ++w)
                  for (int b = 0; b < 64; ++b)
                      if (fp[w] & (1ULL << b)) {
                          int pos = static_cast<int>(w * 64 + b);
                          if (pos < fpSize) bits.push_back(pos);
                      }
              return bits;
          },
          py::arg("mol"), py::arg("radius") = 2,
          py::arg("fp_size") = 2048, py::arg("mode") = "ecfp",
          "Circular fingerprint (Morgan). mode='ecfp' (structural) or 'fcfp' "
          "(pharmacophoric). Returns list of set bit positions.");

    // -----------------------------------------------------------------------
    // Count-based Circular Fingerprints (ECFP/FCFP counts) — v6.0.1
    // -----------------------------------------------------------------------
    m.def("circular_fingerprint_counts",
          [](const smsd::MolGraph& mol, int radius, int fpSize,
             const std::string& mode) {
              auto counts = (mode == "fcfp")
                  ? smsd::batch::detail::computeCircularFingerprintFCFPCounts(mol, radius, fpSize)
                  : smsd::batch::detail::computeCircularFingerprintECFPCounts(mol, radius, fpSize);
              // Return list of (bit_position, count) tuples for non-zero bins
              std::vector<std::pair<int, int>> result;
              for (int i = 0; i < fpSize; ++i) {
                  if (counts[i] > 0) {
                      result.emplace_back(i, counts[i]);
                  }
              }
              return result;
          },
          py::arg("mol"), py::arg("radius") = 2,
          py::arg("fp_size") = 2048, py::arg("mode") = "ecfp",
          "Count-based circular fingerprint (Morgan). mode='ecfp' (structural) "
          "or 'fcfp' (pharmacophoric). Returns list of (bit_position, count) "
          "tuples for non-zero bins. Preserves multiplicity information lost by "
          "binary fingerprints.");

    // -----------------------------------------------------------------------
    // Topological Torsion Fingerprint (Nilakantan et al. 1987) — v6.0.1
    // -----------------------------------------------------------------------
    m.def("topological_torsion",
          [](const smsd::MolGraph& mol, int fpSize) {
              auto fp = smsd::batch::detail::computeTopologicalTorsion(mol, fpSize);
              // Convert uint64_t words to list of set-bit positions
              std::vector<int> bits;
              for (size_t w = 0; w < fp.size(); ++w)
                  for (int b = 0; b < 64; ++b)
                      if (fp[w] & (1ULL << b)) {
                          int pos = static_cast<int>(w * 64 + b);
                          if (pos < fpSize) bits.push_back(pos);
                      }
              return bits;
          },
          py::arg("mol"), py::arg("fp_size") = 2048,
          "Topological torsion fingerprint (Nilakantan 1987). Enumerates all "
          "4-atom linear paths and hashes atom types (atomicNum, heavyDeg, "
          "piElectrons, ring). Returns list of set bit positions.");

    m.def("topological_torsion_counts",
          [](const smsd::MolGraph& mol, int fpSize) {
              auto counts = smsd::batch::detail::computeTopologicalTorsionCounts(mol, fpSize);
              // Return list of (bit_position, count) tuples for non-zero bins
              std::vector<std::pair<int, int>> result;
              for (int i = 0; i < fpSize; ++i) {
                  if (counts[i] > 0) {
                      result.emplace_back(i, counts[i]);
                  }
              }
              return result;
          },
          py::arg("mol"), py::arg("fp_size") = 2048,
          "Count-based topological torsion fingerprint. Returns list of "
          "(bit_position, count) tuples for non-zero bins.");

    // -----------------------------------------------------------------------
    // Tanimoto similarity
    // -----------------------------------------------------------------------
    m.def("tanimoto",
          [](const std::vector<int>& fp1, const std::vector<int>& fp2,
             int fpSize) {
              // Convert set-bit lists to uint64_t arrays, compute Tanimoto
              int numWords = (fpSize + 63) / 64;
              std::vector<uint64_t> a(numWords, 0ULL), b(numWords, 0ULL);
              for (int bit : fp1) if (bit >= 0 && bit < fpSize)
                  a[bit / 64] |= (1ULL << (bit % 64));
              for (int bit : fp2) if (bit >= 0 && bit < fpSize)
                  b[bit / 64] |= (1ULL << (bit % 64));
              return smsd::batch::fingerprintTanimoto(a, b);
          },
          py::arg("fp1"), py::arg("fp2"), py::arg("fp_size") = 2048,
          "Tanimoto similarity between two fingerprints (lists of set bit positions).");

    // -----------------------------------------------------------------------
    // Dice similarity
    // -----------------------------------------------------------------------
    m.def("dice",
          [](const std::vector<int>& fp1, const std::vector<int>& fp2,
             int fpSize) {
              int numWords = (fpSize + 63) / 64;
              std::vector<uint64_t> a(numWords, 0ULL), b(numWords, 0ULL);
              for (int bit : fp1) if (bit >= 0 && bit < fpSize)
                  a[bit / 64] |= (1ULL << (bit % 64));
              for (int bit : fp2) if (bit >= 0 && bit < fpSize)
                  b[bit / 64] |= (1ULL << (bit % 64));
              return smsd::batch::fingerprintDice(a, b);
          },
          py::arg("fp1"), py::arg("fp2"), py::arg("fp_size") = 2048,
          "Dice similarity between two fingerprints (lists of set bit positions). "
          "Dice = 2*|A&B| / (|A| + |B|). Weights shared features more heavily "
          "than Tanimoto.");

    // -----------------------------------------------------------------------
    // Cosine similarity
    // -----------------------------------------------------------------------
    m.def("cosine",
          [](const std::vector<int>& fp1, const std::vector<int>& fp2,
             int fpSize) {
              int numWords = (fpSize + 63) / 64;
              std::vector<uint64_t> a(numWords, 0ULL), b(numWords, 0ULL);
              for (int bit : fp1) if (bit >= 0 && bit < fpSize)
                  a[bit / 64] |= (1ULL << (bit % 64));
              for (int bit : fp2) if (bit >= 0 && bit < fpSize)
                  b[bit / 64] |= (1ULL << (bit % 64));
              return smsd::batch::fingerprintCosine(a, b);
          },
          py::arg("fp1"), py::arg("fp2"), py::arg("fp_size") = 2048,
          "Cosine similarity between two fingerprints (lists of set bit positions). "
          "Cosine = |A&B| / sqrt(|A| * |B|). Normalises for fingerprint density.");

    // -----------------------------------------------------------------------
    // Soergel distance
    // -----------------------------------------------------------------------
    m.def("soergel",
          [](const std::vector<int>& fp1, const std::vector<int>& fp2,
             int fpSize) {
              int numWords = (fpSize + 63) / 64;
              std::vector<uint64_t> a(numWords, 0ULL), b(numWords, 0ULL);
              for (int bit : fp1) if (bit >= 0 && bit < fpSize)
                  a[bit / 64] |= (1ULL << (bit % 64));
              for (int bit : fp2) if (bit >= 0 && bit < fpSize)
                  b[bit / 64] |= (1ULL << (bit % 64));
              return smsd::batch::fingerprintSoergel(a, b);
          },
          py::arg("fp1"), py::arg("fp2"), py::arg("fp_size") = 2048,
          "Soergel distance between two fingerprints (lists of set bit positions). "
          "Soergel = 1 - Tanimoto. A proper metric distance for clustering.");

    // -----------------------------------------------------------------------
    // Count-vector Tanimoto
    // -----------------------------------------------------------------------
    m.def("count_tanimoto",
          [](const std::vector<int>& fp1, const std::vector<int>& fp2) {
              return smsd::batch::countTanimoto(fp1, fp2);
          },
          py::arg("fp1"), py::arg("fp2"),
          "Count-vector Tanimoto: sum(min(a,b)) / sum(max(a,b)). "
          "Generalisation of Tanimoto for integer count vectors (e.g. ECFP counts).");

    // -----------------------------------------------------------------------
    // Count-vector Dice
    // -----------------------------------------------------------------------
    m.def("count_dice",
          [](const std::vector<int>& fp1, const std::vector<int>& fp2) {
              return smsd::batch::countDice(fp1, fp2);
          },
          py::arg("fp1"), py::arg("fp2"),
          "Count-vector Dice: 2*sum(min(a,b)) / (sum(a) + sum(b)). "
          "Generalisation of Dice for integer count vectors.");

    // -----------------------------------------------------------------------
    // Count-vector Cosine
    // -----------------------------------------------------------------------
    m.def("count_cosine",
          [](const std::vector<int>& fp1, const std::vector<int>& fp2) {
              return smsd::batch::countCosine(fp1, fp2);
          },
          py::arg("fp1"), py::arg("fp2"),
          "Count-vector Cosine: dot(a,b) / (|a| * |b|). "
          "Generalisation of cosine similarity for integer count vectors.");

    // -----------------------------------------------------------------------
    // Fingerprint Format Conversions
    // -----------------------------------------------------------------------
    m.def("to_hex",
          [](const std::vector<int>& fp, int fpSize) {
              int numWords = (fpSize + 63) / 64;
              std::vector<uint64_t> words(numWords, 0ULL);
              for (int bit : fp) if (bit >= 0 && bit < fpSize)
                  words[bit / 64] |= (1ULL << (bit % 64));
              return smsd::batch::toHex(words);
          },
          py::arg("fp"), py::arg("fp_size") = 2048,
          "Convert fingerprint (set bit positions) to hexadecimal string.");

    m.def("from_hex",
          [](const std::string& hex) {
              auto words = smsd::batch::fromHex(hex);
              std::vector<int> bits;
              for (size_t w = 0; w < words.size(); ++w)
                  for (int b = 0; b < 64; ++b)
                      if (words[w] & (1ULL << b))
                          bits.push_back(static_cast<int>(w * 64 + b));
              return bits;
          },
          py::arg("hex"),
          "Parse hexadecimal string back to fingerprint (set bit positions).");

    m.def("to_binary_string",
          [](const std::vector<int>& fp, int fpSize) {
              int numWords = (fpSize + 63) / 64;
              std::vector<uint64_t> words(numWords, 0ULL);
              for (int bit : fp) if (bit >= 0 && bit < fpSize)
                  words[bit / 64] |= (1ULL << (bit % 64));
              return smsd::batch::toBinaryString(words, fpSize);
          },
          py::arg("fp"), py::arg("fp_size") = 2048,
          "Convert fingerprint to binary string of 0s and 1s.");

    // -----------------------------------------------------------------------
    // CIP stereo descriptors
    // -----------------------------------------------------------------------
    py::enum_<smsd::cip::RSLabel>(m, "RSLabel",
        "CIP tetrahedral stereo descriptor (R, S, or NONE)")
        .value("NONE", smsd::cip::RSLabel::NONE)
        .value("R",    smsd::cip::RSLabel::R)
        .value("S",    smsd::cip::RSLabel::S)
        .export_values();

    py::enum_<smsd::cip::EZLabel>(m, "EZLabel",
        "CIP double-bond stereo descriptor (E, Z, or NONE)")
        .value("NONE", smsd::cip::EZLabel::NONE)
        .value("E",    smsd::cip::EZLabel::E)
        .value("Z",    smsd::cip::EZLabel::Z)
        .export_values();

    py::class_<smsd::cip::CIPDescriptors>(m, "CIPDescriptors",
        "CIP stereo descriptors for all centres and double bonds")
        .def_readonly("rs_labels", &smsd::cip::CIPDescriptors::rsLabels,
            "Per-atom R/S labels (NONE if not a stereocentre)")
        .def_readonly("ez_bonds",  &smsd::cip::CIPDescriptors::ezBonds,
            "List of (atom1, atom2, EZLabel) for stereogenic double bonds")
        .def("__repr__", [](const smsd::cip::CIPDescriptors& d) {
            int nRS = 0, nEZ = 0;
            for (auto l : d.rsLabels)
                if (l != smsd::cip::RSLabel::NONE) nRS++;
            nEZ = static_cast<int>(d.ezBonds.size());
            return "<CIPDescriptors stereocentres=" + std::to_string(nRS) +
                   " stereo_bonds=" + std::to_string(nEZ) + ">";
        });

    m.def("assign_rs",
          [](const smsd::MolGraph& g, int centre) {
              return smsd::cip::assignRS(g, centre);
          },
          py::arg("mol"), py::arg("centre"),
          "Assign CIP R/S descriptor to a tetrahedral stereocentre.\n"
          "Uses two-pass approach: Rules 1-2, then 4b-c + Rule 5.\n"
          "Returns RSLabel (NONE, R, S, r, s).");

    m.def("assign_ez", &smsd::cip::assignEZ,
          py::arg("mol"), py::arg("a1"), py::arg("a2"),
          "Assign CIP E/Z descriptor to a double bond between atoms a1 and a2");

    m.def("assign_cip", &smsd::cip::assignAll,
          py::arg("mol"),
          "Compute CIP descriptors for all stereocentres and double bonds");

    m.def("assign_cip_from_smiles", &smsd::cip::assignFromSMILES,
          py::arg("smiles"),
          "Parse SMILES and compute CIP descriptors");

    // -----------------------------------------------------------------------
    // SMARTS-based MCS
    // -----------------------------------------------------------------------
    m.def("find_mcs_smarts",
          [](const std::string& smartsStr, const smsd::MolGraph& target,
             int maxMatches) {
              return smsd::findMcsSmarts(smartsStr, target, maxMatches);
          },
          py::arg("smarts"), py::arg("target"),
          py::arg("max_matches") = 1000,
          "Find the largest SMARTS substructure match in a target molecule.\n"
          "Returns a dict mapping SMARTS atom indices to target atom indices.\n"
          "If no match, returns an empty dict.");

    m.def("smarts_match",
          [](const std::string& smartsStr, const smsd::MolGraph& target) {
              auto q = smsd::parseSMARTS(smartsStr);
              return q.matches(target);
          },
          py::arg("smarts"), py::arg("target"),
          "Check if a SMARTS pattern matches a target molecule. Returns bool.");

    m.def("smarts_find_all",
          [](const std::string& smartsStr, const smsd::MolGraph& target,
             int maxMatches) {
              auto q = smsd::parseSMARTS(smartsStr);
              return q.findAll(target, maxMatches);
          },
          py::arg("smarts"), py::arg("target"),
          py::arg("max_matches") = 1000,
          "Find all SMARTS substructure matches in a target molecule.\n"
          "Returns list of dicts (SMARTS atom idx -> target atom idx).");

    py::class_<smsd::SmartsQuery>(m, "SmartsQuery")
        .def("atom_count", &smsd::SmartsQuery::atomCount,
             "Number of atoms in the compiled SMARTS query")
        .def("matching_order", &smsd::SmartsQuery::matchingOrder,
             "Return the historical/debug BFS matching order for the query")
        .def("matches", &smsd::SmartsQuery::matches,
             py::arg("target"),
             py::call_guard<py::gil_scoped_release>(),
             "Check if this compiled SMARTS query matches a target molecule")
        .def("find_all", &smsd::SmartsQuery::findAll,
             py::arg("target"),
             py::arg("max_matches") = 1000,
             py::call_guard<py::gil_scoped_release>(),
             "Find all matches of this compiled SMARTS query in a target molecule")
        .def("matches_many",
             [](const smsd::SmartsQuery& query,
                const std::vector<smsd::MolGraph>& targets) {
                 std::vector<bool> results;
                 results.reserve(targets.size());
                 for (const auto& target : targets) {
                     results.push_back(query.matches(target));
                 }
                 return results;
             },
             py::arg("targets"),
             py::call_guard<py::gil_scoped_release>(),
             "Match one compiled SMARTS query against many targets")
        .def("__len__", &smsd::SmartsQuery::atomCount)
        .def("__repr__", [](const smsd::SmartsQuery& q) {
            return "<SmartsQuery atoms=" + std::to_string(q.atomCount()) + ">";
        });

    m.def("compile_smarts",
          [](const std::string& smartsStr, int maxRecursionDepth) {
              return smsd::parseSMARTS(smartsStr, maxRecursionDepth);
          },
          py::arg("smarts"),
          py::arg("max_recursion_depth") = 20,
          py::call_guard<py::gil_scoped_release>(),
          "Parse SMARTS once and return a reusable compiled query object.");

    // -----------------------------------------------------------------------
    // Ring layout utilities
    // -----------------------------------------------------------------------

    m.def("compute_sssr",
          [](const smsd::MolGraph& g) {
              return smsd::computeSSSR(g);
          },
          py::arg("mol"),
          "Return the Smallest Set of Smallest Rings (minimum cycle basis),\n"
          "pre-sorted ascending by ring size.");

    m.def("layout_sssr",
          [](const smsd::MolGraph& g) {
              return smsd::layoutSSSR(g);
          },
          py::arg("mol"),
          "Return SSSR rings optimized for 2D coordinate generation:\n"
          "largest ring system first, fused rings by adjacency, smallest first.");

    // -----------------------------------------------------------------------
    // Point2D + reduceCrossings
    // -----------------------------------------------------------------------

    py::class_<smsd::Point2D>(m, "Point2D")
        .def(py::init<>())
        .def(py::init([](double x, double y) {
            smsd::Point2D p; p.x = x; p.y = y; return p;
        }), py::arg("x"), py::arg("y"))
        .def_readwrite("x", &smsd::Point2D::x)
        .def_readwrite("y", &smsd::Point2D::y)
        .def("__repr__", [](const smsd::Point2D& p) {
            return "Point2D(" + std::to_string(p.x) + ", " + std::to_string(p.y) + ")";
        });

    m.def("reduce_crossings",
          [](const smsd::MolGraph& g, std::vector<std::vector<double>>& coordsList,
             int maxIter) {
              if (coordsList.size() < static_cast<size_t>(g.n))
                  coordsList.resize(static_cast<size_t>(g.n), std::vector<double>(2, 0.0));
              for (auto& coord : coordsList) {
                  if (coord.size() < 2) coord.resize(2, 0.0);
              }
              // Convert list of [x, y] pairs to vector<Point2D>
              std::vector<smsd::Point2D> coords(coordsList.size());
              for (size_t i = 0; i < coordsList.size(); i++) {
                  if (coordsList[i].size() >= 2) {
                      coords[i].x = coordsList[i][0];
                      coords[i].y = coordsList[i][1];
                  }
              }
              int result = smsd::reduceCrossings(g, coords, maxIter);
              // Copy back
              for (size_t i = 0; i < coords.size(); i++) {
                  coordsList[i][0] = coords[i].x;
                  coordsList[i][1] = coords[i].y;
              }
              return py::make_tuple(result, coordsList);
          },
          py::arg("mol"), py::arg("coords"), py::arg("max_iter") = 1000,
          "Reduce bond crossings in a 2D layout by optimizing ring orientations.\n"
          "coords: list of [x, y] pairs for each atom (modified in place).\n"
          "Returns (crossings_remaining, updated_coords).");

    // -----------------------------------------------------------------------
    // Force-directed layout
    // -----------------------------------------------------------------------

    m.def("force_directed_layout",
          [](const smsd::MolGraph& g, std::vector<std::vector<double>>& coordsList,
             int maxIter, double targetBondLength) {
              if (coordsList.size() < static_cast<size_t>(g.n))
                  coordsList.resize(static_cast<size_t>(g.n), std::vector<double>(2, 0.0));
              for (auto& coord : coordsList) {
                  if (coord.size() < 2) coord.resize(2, 0.0);
              }
              std::vector<smsd::Point2D> coords(coordsList.size());
              for (size_t i = 0; i < coordsList.size(); i++) {
                  if (coordsList[i].size() >= 2) {
                      coords[i].x = coordsList[i][0];
                      coords[i].y = coordsList[i][1];
                  }
              }
              double stress = smsd::forceDirectedLayout(g, coords, maxIter, targetBondLength);
              for (size_t i = 0; i < coords.size(); i++) {
                  coordsList[i][0] = coords[i].x;
                  coordsList[i][1] = coords[i].y;
              }
              return py::make_tuple(stress, coordsList);
          },
          py::arg("mol"), py::arg("coords"),
          py::arg("max_iter") = 500, py::arg("target_bond_length") = 1.5,
          "Force-directed 2D layout minimisation.\n"
          "Iteratively moves atoms to minimise bond-length stress, non-bonded\n"
          "repulsion, and crossing penalties.\n"
          "coords: list of [x, y] pairs for each atom (modified in place).\n"
          "Returns (final_stress, updated_coords).");

    // -----------------------------------------------------------------------
    // Stress majorisation (SMACOF)
    // -----------------------------------------------------------------------

    m.def("stress_majorisation",
          [](const smsd::MolGraph& g, std::vector<std::vector<double>>& coordsList,
             int maxIter, double targetBondLength) {
              if (coordsList.size() < static_cast<size_t>(g.n))
                  coordsList.resize(static_cast<size_t>(g.n), std::vector<double>(2, 0.0));
              for (auto& coord : coordsList) {
                  if (coord.size() < 2) coord.resize(2, 0.0);
              }
              std::vector<smsd::Point2D> coords(coordsList.size());
              for (size_t i = 0; i < coordsList.size(); i++) {
                  if (coordsList[i].size() >= 2) {
                      coords[i].x = coordsList[i][0];
                      coords[i].y = coordsList[i][1];
                  }
              }
              double stress = smsd::stressMajorisation(g, coords, maxIter, targetBondLength);
              for (size_t i = 0; i < coords.size(); i++) {
                  coordsList[i][0] = coords[i].x;
                  coordsList[i][1] = coords[i].y;
              }
              return py::make_tuple(stress, coordsList);
          },
          py::arg("mol"), py::arg("coords"),
          py::arg("max_iter") = 300, py::arg("target_bond_length") = 1.5,
          "Stress majorisation (SMACOF algorithm) for 2D layout.\n"
          "Minimises weighted stress to produce graph-distance-proportional layouts.\n"
          "coords: list of [x, y] pairs for each atom (modified in place).\n"
          "Returns (final_stress, updated_coords).");

    // -----------------------------------------------------------------------
    // Template matching
    // -----------------------------------------------------------------------

    m.def("match_template",
          [](const smsd::MolGraph& g, double targetBondLength) {
              auto result = smsd::matchTemplate(g, targetBondLength);
              if (result.empty()) return py::list();
              py::list out;
              for (auto& p : result) {
                  out.append(py::make_tuple(p.x, p.y));
              }
              return out;
          },
          py::arg("mol"), py::arg("target_bond_length") = 1.5,
          "Check if a molecule matches a known scaffold template.\n"
          "Returns list of (x, y) coordinate tuples if matched, empty list otherwise.\n"
          "Supports: benzene, naphthalene, indole, purine, steroid, morphinan,\n"
          "quinoline, biphenyl, cyclohexane, piperidine.");
}
