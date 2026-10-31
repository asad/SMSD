/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 *
 * Python bindings for SMSD via pybind11.
 * Usage: import smsd; result = smsd.MCS.find(query).match(target)
 */

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "smsd/smsd.hpp"
#include "smsd/screen_accelerator.hpp"
#include "smsd/similarity_screen.hpp"
#include "smsd/smiles_parser.hpp"
#include "smsd/smarts_parser.hpp"
#include "smsd/cip.hpp"
#include "smsd/layout.hpp"
#include "smsd/depict.hpp"
#include "smsd/substructure.hpp"

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstring>
#include <functional>
#include <map>
#include <optional>
#include <limits>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <vector>

namespace py = pybind11;

// The build defines the package version; without it the package rejects the extension on import.
#ifndef SMSD_VERSION_INFO
#define SMSD_VERSION_INFO "unknown"
#endif

// ============================================================================
// Binding helpers
// ============================================================================

namespace {

/** Keep metadata mutations attached to the owning molecule. */
struct PropertyView {
    py::object owner;
    std::map<std::string, std::string>& values() const {
        return owner.cast<smsd::MolGraph&>().properties;
    }
};

/** The words of an fp_size-bit fingerprint given as set-bit positions; other positions are ignored. */
std::vector<uint64_t> fingerprintWords(const std::vector<int>& positions, int fpSize) {
    if (fpSize <= 0) throw std::invalid_argument("fp_size must be positive");
    return smsd::bitWords(positions, fpSize);
}

/** Fail with IndexError unless both atoms are in the graph. */
void requireAtomPair(const smsd::MolGraph& g, int i, int j) {
    if (i < 0 || i >= g.n || j < 0 || j >= g.n) throw py::index_error("atom index out of range");
}

/** The (bin, count) pairs of the non-zero bins of a dense count vector. */
std::vector<std::pair<int, int>> nonZeroBins(const std::vector<int>& counts) {
    std::vector<std::pair<int, int>> bins;
    for (size_t i = 0; i < counts.size(); ++i)
        if (counts[i] > 0) bins.emplace_back(static_cast<int>(i), counts[i]);
    return bins;
}

/** A lower-case copy of a Python string, for option names given as text. */
std::string lowered(const py::object& value) {
    std::string text = py::cast<std::string>(value);
    for (char& c : text) c = smsd::parser_detail::asciiLower(c);
    return text;
}

/** A FingerprintFamily from its enum member or its name as text. */
smsd::batch::FingerprintFamily fingerprintFamily(const py::object& value) {
    if (py::isinstance<smsd::batch::FingerprintFamily>(value)) return value.cast<smsd::batch::FingerprintFamily>();
    if (!py::isinstance<py::str>(value)) throw py::type_error("family must be a FingerprintFamily or its name");
    const std::string name = lowered(value);
    if (name == "ecfp" || name == "circular") return smsd::batch::FingerprintFamily::ECFP;
    if (name == "fcfp") return smsd::batch::FingerprintFamily::FCFP;
    if (name == "path") return smsd::batch::FingerprintFamily::PATH;
    if (name == "torsion") return smsd::batch::FingerprintFamily::TORSION;
    if (name == "keys") return smsd::batch::FingerprintFamily::KEYS;
    throw std::invalid_argument("family must be 'ecfp', 'fcfp', 'path', 'torsion' or 'keys', got '" + name + "'");
}

smsd::batch::FragmentMode fragmentMode(const py::object& value) {
    if (py::isinstance<smsd::batch::FragmentMode>(value)) return value.cast<smsd::batch::FragmentMode>();
    if (!py::isinstance<py::str>(value)) throw py::type_error("fragments must be a FragmentMode or its name");
    const std::string name = lowered(value);
    if (name == "whole") return smsd::batch::FragmentMode::WHOLE;
    if (name == "largest") return smsd::batch::FragmentMode::LARGEST;
    if (name == "per_fragment") return smsd::batch::FragmentMode::PER_FRAGMENT;
    throw std::invalid_argument("fragments must be 'whole', 'largest' or 'per_fragment', got '" + name + "'");
}

smsd::batch::ScaffoldMode scaffoldMode(const py::object& value) {
    if (py::isinstance<smsd::batch::ScaffoldMode>(value)) return value.cast<smsd::batch::ScaffoldMode>();
    if (!py::isinstance<py::str>(value)) throw py::type_error("scaffold must be a ScaffoldMode or its name");
    const std::string name = lowered(value);
    if (name == "none") return smsd::batch::ScaffoldMode::NONE;
    if (name == "murcko") return smsd::batch::ScaffoldMode::MURCKO;
    throw std::invalid_argument("scaffold must be 'none' or 'murcko', got '" + name + "'");
}

int integerOption(const py::object& value, const char* name);

/** Keys given as FingerprintKey objects, SMARTS strings or (smarts, minimum) pairs. */
std::vector<smsd::batch::FingerprintKey> fingerprintKeys(const py::object& value) {
    std::vector<smsd::batch::FingerprintKey> keys;
    if (py::isinstance<py::str>(value)) throw py::type_error("keys must be a sequence of keys, not one string");
    for (py::handle item : value) {
        if (py::isinstance<smsd::batch::FingerprintKey>(item)) {
            keys.push_back(item.cast<smsd::batch::FingerprintKey>());
        } else if (py::isinstance<py::str>(item)) {
            keys.emplace_back(item.cast<std::string>(), 1);
        } else if (py::isinstance<py::tuple>(item) || py::isinstance<py::list>(item)) {
            const auto pair = py::reinterpret_borrow<py::sequence>(item);
            if (pair.size() != 2 || !py::isinstance<py::str>(pair[0])) throw py::type_error("a key pair is (smarts, minimum)");
            keys.emplace_back(pair[0].cast<std::string>(),
                              integerOption(py::reinterpret_borrow<py::object>(pair[1]), "a key minimum"));
        } else {
            throw py::type_error("a key is a FingerprintKey, a SMARTS string or a (smarts, minimum) pair");
        }
    }
    return keys;
}

/** An integer option that is not a bool and fits a C int. */
int integerOption(const py::object& value, const char* name) {
    if (py::isinstance<py::bool_>(value) || !py::isinstance<py::int_>(value))
        throw py::type_error(std::string(name) + " must be an integer");
    long long wide = 0;
    try {
        wide = value.cast<long long>();
    } catch (const py::cast_error&) {
        throw py::value_error(std::string(name) + " is out of range");
    }
    if (wide < std::numeric_limits<int>::min() || wide > std::numeric_limits<int>::max())
        throw py::value_error(std::string(name) + " is out of range");
    return static_cast<int>(wide);
}

bool switchOption(const py::object& value, const char* name) {
    if (!py::isinstance<py::bool_>(value)) throw py::type_error(std::string(name) + " must be True or False");
    return value.cast<bool>();
}

/** Applies one keyword argument of FingerprintOptions(...). */
void applyFingerprintOption(smsd::batch::FingerprintOptions& options, const std::string& name, const py::object& value) {
    using smsd::batch::FingerprintOptions;
    static const std::map<std::string, bool FingerprintOptions::*> switches = {
        {"counts", &FingerprintOptions::counts}, {"element", &FingerprintOptions::element},
        {"degree", &FingerprintOptions::degree}, {"valence", &FingerprintOptions::valence},
        {"formal_charge", &FingerprintOptions::formalCharge}, {"isotope", &FingerprintOptions::isotope},
        {"hydrogen_count", &FingerprintOptions::hydrogenCount}, {"tautomer_class", &FingerprintOptions::tautomerClass},
        {"aromaticity", &FingerprintOptions::aromaticity}, {"bond_order", &FingerprintOptions::bondOrder},
        {"aromatic_bond", &FingerprintOptions::aromaticBond}, {"ring_bond", &FingerprintOptions::ringBond},
        {"ring_membership", &FingerprintOptions::ringMembership}, {"ring_size", &FingerprintOptions::ringSize},
        {"cip_atoms", &FingerprintOptions::cipAtoms}, {"cip_bonds", &FingerprintOptions::cipBonds},
    };
    static const std::map<std::string, int FingerprintOptions::*> integers = {
        {"radius", &FingerprintOptions::radius}, {"path_length", &FingerprintOptions::pathLength},
        {"size", &FingerprintOptions::size}, {"ring_size_cap", &FingerprintOptions::ringSizeCap},
    };
    const auto flag = switches.find(name);
    if (flag != switches.end()) { options.*(flag->second) = switchOption(value, name.c_str()); return; }
    const auto number = integers.find(name);
    if (number != integers.end()) { options.*(number->second) = integerOption(value, name.c_str()); return; }
    if (name == "fragments") { options.fragments = fragmentMode(value); return; }
    if (name == "scaffold") { options.scaffold = scaffoldMode(value); return; }
    if (name == "keys") {
        options.keys = fingerprintKeys(value);
        if (options.family == smsd::batch::FingerprintFamily::KEYS) options.size = static_cast<int>(options.keys.size());
        return;
    }
    throw py::type_error("FingerprintOptions() got an unexpected keyword argument '" + name + "'");
}

std::string describeFingerprintOptions(const smsd::batch::FingerprintOptions& options) {
    using smsd::batch::FingerprintFamily;
    using smsd::batch::FragmentMode;
    using smsd::batch::ScaffoldMode;
    const char* family = options.family == FingerprintFamily::ECFP ? "ecfp" : options.family == FingerprintFamily::FCFP ? "fcfp"
                       : options.family == FingerprintFamily::PATH ? "path" : options.family == FingerprintFamily::TORSION ? "torsion" : "keys";
    std::string text = std::string("FingerprintOptions(family='") + family + "'";
    if (options.family == FingerprintFamily::KEYS) text += ", keys=" + std::to_string(options.keys.size());
    else text += ", size=" + std::to_string(options.size);
    if (options.family == FingerprintFamily::ECFP || options.family == FingerprintFamily::FCFP)
        text += ", radius=" + std::to_string(options.radius);
    if (options.family == FingerprintFamily::PATH) text += ", path_length=" + std::to_string(options.pathLength);
    if (options.counts) text += ", counts=True";
    const std::pair<const char*, bool> switches[] = {
        {"element", options.element}, {"degree", options.degree}, {"valence", options.valence},
        {"formal_charge", options.formalCharge}, {"isotope", options.isotope}, {"hydrogen_count", options.hydrogenCount},
        {"tautomer_class", options.tautomerClass}, {"aromaticity", options.aromaticity}, {"bond_order", options.bondOrder},
        {"aromatic_bond", options.aromaticBond}, {"ring_bond", options.ringBond}, {"ring_membership", options.ringMembership},
        {"ring_size", options.ringSize}, {"cip_atoms", options.cipAtoms}, {"cip_bonds", options.cipBonds},
    };
    if (options.family != FingerprintFamily::KEYS) {
        std::string on;
        for (const auto& entry : switches) if (entry.second) on += (on.empty() ? "" : ", ") + std::string(entry.first);
        text += ", terms=[" + on + "]";
        if (options.ringSize) text += ", ring_size_cap=" + std::to_string(options.ringSizeCap);
    }
    if (options.fragments != FragmentMode::WHOLE)
        text += std::string(", fragments='") + (options.fragments == FragmentMode::LARGEST ? "largest" : "per_fragment") + "'";
    if (options.scaffold != ScaffoldMode::NONE) text += ", scaffold='murcko'";
    return text + ")";
}

/** Whether a circular fingerprint mode names the functional-class variant; only "ecfp" and "fcfp" are known. */
bool circularMode(const std::string& mode) {
    if (mode == "fcfp") return true;
    if (mode == "ecfp") return false;
    throw std::invalid_argument("mode must be 'ecfp' or 'fcfp', got '" + mode + "'");
}

/** How a coordinate pair with fewer than two values is read. */
enum class ShortPairs { ZeroFill, Reject };

/**
 * The (x, y) pairs of `coords` as points. The result holds at least `minimumSize` points, and exactly
 * that many when `truncate` is set; points without a pair are (0, 0). A pair shorter than two is (0, 0)
 * under ZeroFill and raises IndexError under Reject, as indexing it would.
 */
std::vector<smsd::Point2D> points2D(const py::sequence& coords, size_t minimumSize = 0, bool truncate = false,
                                    ShortPairs shortPairs = ShortPairs::ZeroFill) {
    const size_t given = static_cast<size_t>(coords.size());
    const size_t size = truncate ? minimumSize : std::max(given, minimumSize);
    std::vector<smsd::Point2D> points(size, smsd::Point2D{0.0, 0.0});
    for (size_t i = 0; i < given && i < size; ++i) {
        const py::sequence pair = coords[i].cast<py::sequence>();
        if (shortPairs == ShortPairs::ZeroFill && pair.size() < 2) continue;
        points[i].x = pair[0].cast<double>();
        points[i].y = pair[1].cast<double>();
    }
    return points;
}

/** The points as a list of (x, y) tuples. */
py::list pointList(const std::vector<smsd::Point2D>& points) {
    py::list out(points.size());
    for (size_t i = 0; i < points.size(); ++i)
        out[i] = py::make_tuple(points[i].x, points[i].y);
    return out;
}

struct PythonGraphView {
    smsd::batch::detail::MolGraphView graphs;
    std::vector<py::object> owners;
};

PythonGraphView graphView(const py::sequence& targets) {
    PythonGraphView view;
    view.graphs.reserve(targets.size());
    view.owners.reserve(targets.size());
    // Each item is held by an owned reference, as a sequence may create its items afresh.
    for (py::object target : targets) {
        view.graphs.push_back(target.cast<const smsd::MolGraph&>());
        view.owners.push_back(std::move(target));
    }
    return view;
}

/** Capture a fixed sequence without copying a repeated molecule at every slot. */
void addCorpusTargets(smsd::batch::TargetCorpus& corpus, const py::iterable& molecules) {
    if (PyList_CheckExact(molecules.ptr()) || PyTuple_CheckExact(molecules.ptr())) {
        smsd::batch::detail::MolGraphView references;
        references.reserve(py::len(molecules));
        for (py::handle molecule : molecules) references.push_back(molecule.cast<const smsd::MolGraph&>());
        auto captured = smsd::batch::detail::TargetStorage::capture(references);
        py::gil_scoped_release release;
        smsd::batch::detail::CorpusAccess::addCapturedTargets(corpus, std::move(captured));
        return;
    }
    // Iterators may mutate and yield the same object again: capture each yielded value.
    std::vector<smsd::MolGraph> graphs;
    for (py::handle molecule : molecules) {
        const auto& graph = molecule.cast<const smsd::MolGraph&>();
        const std::lock_guard<std::recursive_mutex> lock(graph.lazyCacheMutex());
        graphs.push_back(graph);
    }
    py::gil_scoped_release release;
    corpus.addTargets(std::move(graphs));
}

/** Private search snapshots own their graphs; repeated slots share only these detached copies. */
py::list corpusSearchTargets(const smsd::batch::TargetCorpus& corpus) {
    std::vector<smsd::MolGraph> graphs;
    std::vector<size_t> slots;
    {
        py::gil_scoped_release release;
        const auto targets = smsd::batch::detail::CorpusAccess::targetView(corpus);
        graphs.reserve(targets.storedSize());
        slots.reserve(targets.size());
        std::unordered_map<const smsd::MolGraph*, size_t> copied;
        for (size_t index = 0; index < targets.size(); ++index) {
            const auto& graph = targets[index];
            auto found = copied.find(&graph);
            if (found == copied.end()) {
                const std::lock_guard<std::recursive_mutex> lock(graph.lazyCacheMutex());
                graphs.push_back(graph);
                found = copied.emplace(&graph, graphs.size() - 1).first;
            }
            slots.push_back(found->second);
        }
    }
    std::vector<py::object> values;
    values.reserve(graphs.size());
    for (auto& graph : graphs) values.push_back(py::cast(std::move(graph)));
    py::list result(slots.size());
    for (size_t index = 0; index < slots.size(); ++index) result[index] = values[slots[index]];
    return result;
}

/**
 * The targets of a similarity-bound screen or a fingerprint batch, read from a Python iterable
 * while the GIL is held. MolGraph items are referenced rather than copied, and SMILES strings are
 * parsed by the batch's workers with the parser parse_smiles uses; a parsed graph is dropped once
 * it is scored or fingerprinted. Any other item goes through convert, the package's molecule
 * conversion, in input order. Reading stops at the first item that conversion or iteration
 * rejects, and that error is raised unless a string before it fails to parse, so the error raised
 * is that of the first item that fails, as when the items are converted one by one.
 */
class ScreenTargets {
    struct Item {
        const smsd::MolGraph* graph;   // null for a SMILES string
        const char* text;              // the string's UTF-8 text, held by its Python object
        size_t length;
    };

    std::vector<Item> items_;
    std::vector<py::object> owners_;
    std::optional<py::error_already_set> stopped_;

    /** Hold the graph convert gives for an item; false when it raises, which stops reading. */
    bool convertItem(const py::object& convert, const py::object& item) {
        if (convert.is_none()) {
            if (!PyUnicode_Check(item.ptr()))
                PyErr_Format(PyExc_TypeError, "Expected MolGraph or SMILES string, got %s", Py_TYPE(item.ptr())->tp_name);
            else if (PyUnicode_GetLength(item.ptr()) == 0)
                PyErr_SetString(PyExc_ValueError, "SMILES string must not be empty");
            else
                PyErr_SetString(PyExc_TypeError, "SMILES string must be valid UTF-8");
            stopped_.emplace();
            return false;
        }
        try {
            py::object graph = convert(item);
            items_.push_back({&graph.cast<const smsd::MolGraph&>(), nullptr, 0});
            owners_.push_back(std::move(graph));
            return true;
        } catch (py::error_already_set& error) {
            stopped_.emplace(std::move(error));
            return false;
        } catch (const py::cast_error&) {
            PyErr_SetString(PyExc_TypeError, "convert must return a MolGraph");
            stopped_.emplace();
            return false;
        }
    }

public:
    /** A target as the screen reads it: a held graph, or one parsed for this read. */
    struct Target {
        const smsd::MolGraph* held;
        std::optional<smsd::MolGraph> parsed;
        operator const smsd::MolGraph&() const { return parsed ? *parsed : *held; }
    };

    ScreenTargets(const py::handle& targets, const py::object& convert) {
        auto* graphType = reinterpret_cast<PyTypeObject*>(py::type::of<smsd::MolGraph>().ptr());
        if (PyList_Check(targets.ptr()) || PyTuple_Check(targets.ptr())) {
            const auto count = static_cast<size_t>(Py_SIZE(targets.ptr()));
            items_.reserve(count);
            owners_.reserve(count);
        }
        const auto iterator = py::reinterpret_steal<py::object>(PyObject_GetIter(targets.ptr()));
        if (!iterator) throw py::error_already_set();
        while (PyObject* next = PyIter_Next(iterator.ptr())) {
            py::object item = py::reinterpret_steal<py::object>(next);
            Py_ssize_t length = 0;
            const char* text = nullptr;
            if (PyObject_TypeCheck(next, graphType)) {
                items_.push_back({&item.cast<const smsd::MolGraph&>(), nullptr, 0});
            } else if (PyUnicode_Check(next) && (text = PyUnicode_AsUTF8AndSize(next, &length)) != nullptr
                       && length > 0) {
                items_.push_back({nullptr, text, static_cast<size_t>(length)});
            } else {
                // An empty string, or one that is not valid UTF-8, is left to the conversion, which rejects it.
                if (PyErr_Occurred()) PyErr_Clear();
                if (!convertItem(convert, item)) return;
            }
            owners_.push_back(std::move(item));
        }
        if (PyErr_Occurred()) stopped_.emplace();
    }

    size_t size() const { return items_.size(); }

    /** The target at index; a SMILES string that fails to parse throws as parse_smiles does. */
    Target operator[](size_t index) const {
        const Item& item = items_[index];
        if (item.graph != nullptr) return {item.graph, std::nullopt};
        return {nullptr, smsd::parseSMILES(std::string(item.text, item.length))};
    }

    /** Raise the error that stopped reading, if any; needs the GIL. */
    void raiseStop() {
        if (!stopped_) return;
        py::error_already_set error = std::move(*stopped_);
        stopped_.reset();
        throw error;
    }
};

/** TypeError for a screen argument the binding cannot take, raised after the targets' errors. */
[[noreturn]] void rejectScreenArgument(const char* name, const char* expected, const py::handle& value) {
    throw py::type_error(std::string("screen_targets(): incompatible function arguments; ") + name
                         + " must be " + expected + ", not " + Py_TYPE(value.ptr())->tp_name);
}

/**
 * Screen Python targets by the similarity bound on up to numThreads threads (0 = the OpenMP thread limit),
 * with the GIL released while the targets are parsed and scored. The threshold and chemistry are
 * taken after the targets, as when the targets were converted before the native call, and the
 * chemistry is copied before the GIL is released.
 */
std::vector<int> screenPythonTargets(const smsd::MolGraph& query, const py::object& targets,
                                     const py::object& thresholdValue, const py::object& chemistryValue,
                                     int numThreads, const py::object& convert) {
    if (numThreads < 0) throw std::invalid_argument("Thread count must be non-negative");
    ScreenTargets view(targets, convert);
    std::optional<double> threshold;
    try {
        threshold = thresholdValue.cast<double>();
    } catch (const py::cast_error&) {
        // Raised after the targets' errors, as below.
    }
    const bool chemistryTaken = py::isinstance<smsd::ChemOptions>(chemistryValue);
    const smsd::ChemOptions chemistry = chemistryTaken ? chemistryValue.cast<const smsd::ChemOptions&>()
                                                       : smsd::ChemOptions();
    std::vector<double> scores;
    {
        // A SMILES string that fails to parse propagates from here, before any later error.
        py::gil_scoped_release release;
        scores = smsd::screening_detail::similarityScores(query, view, chemistry, numThreads);
    }
    view.raiseStop();
    if (!threshold) rejectScreenArgument("threshold", "a float", thresholdValue);
    if (!chemistryTaken) rejectScreenArgument("chem", "a ChemOptions", chemistryValue);
    return smsd::screening_detail::indicesReaching(scores, *threshold);
}

/**
 * The threads a Python batch asks for: num_threads, or the OpenMP worker count (which follows
 * OMP_NUM_THREADS) when it is 0. A build without OpenMP searches on the calling thread.
 */
int batchThreads(int numThreads) {
#ifdef _OPENMP
    return numThreads == 0 ? omp_get_max_threads() : numThreads;
#else
    (void)numThreads;
    return 1;
#endif
}

/**
 * Search every target with the native batch dispatcher, whose threads end with the call. The
 * failure of the lowest failing target propagates once the GIL is held again.
 */
template<class Query>
std::vector<smsd::SearchResult> matchView(const Query& query, const PythonGraphView& view, int numThreads,
                                        const smsd::batch::TargetCorpus* corpus = nullptr) {
    const int threads = batchThreads(numThreads);
    py::gil_scoped_release release;
    if constexpr (std::is_same<Query, smsd::Substructure>::value) {
        if (corpus != nullptr) {
            // An addition after the snapshot changes its size; the captured graphs still use the CPU.
            if (auto survivors = query.corpusCandidates(*corpus, threads))
                if (survivors->size() == view.graphs.size())
                    return query.matchManyScreened(view.graphs, *survivors, threads);
        }
    }
    return query.matchMany(view.graphs, threads);
}

template<class Query>
std::vector<smsd::SearchResult> matchMany(const Query& query, const py::sequence& targets, int numThreads) {
    if (numThreads < 0) throw std::invalid_argument("Thread count must be non-negative");
    return matchView(query, graphView(targets), numThreads);
}

/** The set positions of fingerprint(mol, kind=kind), computed as its single-molecule binding does. */
std::function<std::vector<int>(const smsd::MolGraph&)> kindFingerprint(const std::string& kind, int radius,
                                                                       int pathLength, int fpSize) {
    namespace fp = smsd::batch::detail;
    if (kind == "path")
        return [=](const smsd::MolGraph& mol) {
            return smsd::bitPositions(fp::computeElementPathFingerprint(mol, pathLength, fpSize), fpSize);
        };
    if (kind == "mcs")
        return [=](const smsd::MolGraph& mol) {
            return smsd::bitPositions(
                fp::computeMCSFingerprintForChemistry(mol, pathLength, fpSize, smsd::ChemOptions()), fpSize);
        };
    if (kind == "ecfp" || kind == "circular" || kind == "fcfp") {
        const bool functional = kind == "fcfp";
        return [=](const smsd::MolGraph& mol) {
            return smsd::bitPositions(functional ? fp::computeCircularFingerprintFCFP(mol, radius, fpSize)
                                                 : fp::computeCircularFingerprintECFP(mol, radius, fpSize), fpSize);
        };
    }
    if (kind == "torsion")
        return [=](const smsd::MolGraph& mol) {
            return smsd::bitPositions(fp::computeTopologicalTorsion(mol, fpSize), fpSize);
        };
    throw std::invalid_argument("Unknown fingerprint kind: '" + kind
                                + "'. Use 'path', 'mcs', 'ecfp', 'circular', 'fcfp', or 'torsion'.");
}

/**
 * The fingerprint of every molecule, in input order, as fingerprint() gives each one under the
 * FingerprintOptions or the kind that spec holds, on up to numThreads threads (0 = the OpenMP
 * thread limit) with the GIL released. Molecules are read as screen_targets reads its targets,
 * so the error raised is that of the first molecule that fails.
 */
py::object pythonFingerprints(const py::object& molecules, const py::object& spec, int radius, int pathLength,
                              int fpSize, int numThreads, const py::object& convert) {
    using smsd::batch::FingerprintOptions;
    if (numThreads < 0) throw std::invalid_argument("Thread count must be non-negative");
    const bool optionsGiven = py::isinstance<FingerprintOptions>(spec);
    // The options are copied while the GIL is held, as for the searches.
    const FingerprintOptions options = optionsGiven ? spec.cast<const FingerprintOptions&>() : FingerprintOptions();
    std::function<std::vector<int>(const smsd::MolGraph&)> bits;
    if (!optionsGiven) bits = kindFingerprint(spec.cast<std::string>(), radius, pathLength, fpSize);
    ScreenTargets view(molecules, convert);
    const int threads = batchThreads(numThreads);
    py::object fingerprints;
    if (optionsGiven && options.counts) {
        std::vector<std::vector<std::pair<int, int>>> counts;
        {
            py::gil_scoped_release release;
            counts = smsd::batch::detail::fingerprintEach(view, threads, [&options](const smsd::MolGraph& mol) {
                return nonZeroBins(smsd::batch::computeFingerprintCounts(mol, options));
            });
        }
        fingerprints = py::cast(counts);
    } else {
        const int size = options.family == smsd::batch::FingerprintFamily::KEYS
                       ? static_cast<int>(options.keys.size()) : options.size;
        std::vector<std::vector<int>> positions;
        {
            py::gil_scoped_release release;
            positions = optionsGiven
                ? smsd::batch::detail::fingerprintEach(view, threads, [&options, size](const smsd::MolGraph& mol) {
                      return smsd::bitPositions(smsd::batch::computeFingerprint(mol, options), size);
                  })
                : smsd::batch::detail::fingerprintEach(view, threads, bits);
        }
        fingerprints = py::cast(positions);
    }
    view.raiseStop();
    return fingerprints;
}

// ----------------------------------------------------------------------------
// Direct construction of the Python search results
// ----------------------------------------------------------------------------
//
// smsd.search registers its immutable Match and SearchResult classes once.
// Searches then build those instances here in a single call, so a prepared
// search does not materialise intermediate wrapper objects, lists and dicts
// that Python would only copy again.  Module initialisation starts the cached
// references below afresh: any held before belong to an interpreter that has
// since been finalised, so they are dropped without being released, and no
// Python object is released after its interpreter has finalised.

enum ResultField {
    AtomMappingField, BondMappingField, FragmentsField,
    MatchesField, StatusField, OptimalField, CompleteField, ScoreField,
    NodesVisitedField, ElapsedField, ResultFieldCount
};

/** The attribute of each field; those of Match come first, then those of SearchResult. */
const char* const resultFieldNames[ResultFieldCount] = {
    "atom_mapping", "bond_mapping", "fragments",
    "matches", "status", "optimal", "complete", "score", "nodes_visited", "elapsed_ms"};

struct PythonResultTypes {
    PyObject* match = nullptr;
    PyObject* result = nullptr;
    PyObject* status[4] = {nullptr, nullptr, nullptr, nullptr};
    PyObject* names[ResultFieldCount] = {};
};

PythonResultTypes& resultTypes() {
    static PythonResultTypes types;
    return types;
}

/** Own a new reference, or raise the pending Python error. */
py::object owned(PyObject* object) {
    if (object == nullptr) throw py::error_already_set();
    return py::reinterpret_steal<py::object>(object);
}

/** Store a new reference in a cached slot before releasing the one it held. */
void replaceReference(PyObject*& slot, PyObject* value) {
    PyObject* previous = slot;
    slot = value;
    Py_XDECREF(previous);
}

/**
 * Start the cached references afresh for the interpreter that imports the module: the field
 * names and the members of the search status enumeration, which results share rather than
 * allocating a copy per search.
 */
void initialiseResultTypes(const py::object& searchStatus) {
    auto& types = resultTypes();
    types = PythonResultTypes{};
    for (int field = 0; field < ResultFieldCount; ++field)
        types.names[field] = owned(PyUnicode_InternFromString(resultFieldNames[field])).release().ptr();
    const std::pair<smsd::SearchStatus, const char*> members[] = {
        {smsd::SearchStatus::Completed, "COMPLETED"}, {smsd::SearchStatus::TimedOut, "TIMED_OUT"},
        {smsd::SearchStatus::Cancelled, "CANCELLED"}, {smsd::SearchStatus::LimitReached, "LIMIT_REACHED"}};
    for (const auto& member : members) {
        const auto index = static_cast<size_t>(member.first);
        if (index < 4) types.status[index] = py::object(searchStatus.attr(member.second)).release().ptr();
    }
}

/**
 * Map native atom indices to the caller's numbering, or keep them as they are. Every number
 * is checked when the numbering is taken, before any result is built: only an exact int is
 * known to hash and compare without running Python code.
 */
class AtomNumbering {
    py::object indices_;
public:
    /** Keep native indices. */
    AtomNumbering() = default;

    explicit AtomNumbering(py::handle indices) {
        if (indices.is_none()) return;
        indices_ = PyTuple_CheckExact(indices.ptr()) ? py::reinterpret_borrow<py::object>(indices)
                                                     : owned(PySequence_Tuple(indices.ptr()));
        for (Py_ssize_t index = 0, size = PyTuple_GET_SIZE(indices_.ptr()); index < size; ++index)
            if (!PyLong_CheckExact(PyTuple_GET_ITEM(indices_.ptr(), index)))
                throw py::type_error("Atom numbering must hold only int values");
    }

    /** Return a new reference to the caller's number for a native atom. */
    PyObject* number(int index) const {
        if (!indices_) {
            PyObject* value = PyLong_FromLong(index);
            if (value == nullptr) throw py::error_already_set();
            return value;
        }
        if (index < 0 || index >= PyTuple_GET_SIZE(indices_.ptr()))
            throw py::index_error("Atom index out of range");
        PyObject* value = PyTuple_GET_ITEM(indices_.ptr(), index);
        Py_INCREF(value);
        return value;
    }
};

/** Create an instance without running __init__, as object.__new__(type) would. */
py::object newInstance(PyObject* type) {
    const py::object noArguments = owned(PyTuple_New(0));
    return owned(PyBaseObject_Type.tp_new(reinterpret_cast<PyTypeObject*>(type), noArguments.ptr(), nullptr));
}

/** Store a field directly; the frozen class rejects assignment through __setattr__. */
void setField(const py::object& instance, ResultField field, PyObject* value) {
    if (PyObject_GenericSetAttr(instance.ptr(), resultTypes().names[field], value) < 0)
        throw py::error_already_set();
}

/** A tuple of the objects, built only once all of them exist. */
py::object tupleOf(std::vector<py::object>& items) {
    py::object tuple = owned(PyTuple_New(static_cast<Py_ssize_t>(items.size())));
    for (size_t index = 0; index < items.size(); ++index)
        PyTuple_SET_ITEM(tuple.ptr(), static_cast<Py_ssize_t>(index), items[index].release().ptr());
    return tuple;
}

py::object pythonMatch(const smsd::Match& match, const AtomNumbering& query, const AtomNumbering& target) {
    py::object atoms = owned(PyDict_New());
    for (const auto& pair : match.atomMapping) {
        py::object first = py::reinterpret_steal<py::object>(query.number(pair.first));
        py::object second = py::reinterpret_steal<py::object>(target.number(pair.second));
        if (PyDict_SetItem(atoms.ptr(), first.ptr(), second.ptr()) < 0) throw py::error_already_set();
    }
    py::object bonds = owned(PyTuple_New(static_cast<Py_ssize_t>(match.bondMapping.size())));
    for (size_t bond = 0; bond < match.bondMapping.size(); ++bond) {
        PyObject* entry = PyTuple_New(4);
        if (entry == nullptr) throw py::error_already_set();
        PyTuple_SET_ITEM(bonds.ptr(), static_cast<Py_ssize_t>(bond), entry);
        const auto& atomsOfBond = match.bondMapping[bond];
        for (int position = 0; position < 4; ++position)
            PyTuple_SET_ITEM(entry, position, (position < 2 ? query : target).number(atomsOfBond[position]));
    }
    py::object fragments = owned(PyTuple_New(static_cast<Py_ssize_t>(match.fragments.size())));
    for (size_t fragment = 0; fragment < match.fragments.size(); ++fragment) {
        const auto& members = match.fragments[fragment];
        PyObject* entry = PyTuple_New(static_cast<Py_ssize_t>(members.size()));
        if (entry == nullptr) throw py::error_already_set();
        PyTuple_SET_ITEM(fragments.ptr(), static_cast<Py_ssize_t>(fragment), entry);
        for (size_t member = 0; member < members.size(); ++member)
            PyTuple_SET_ITEM(entry, static_cast<Py_ssize_t>(member), query.number(members[member]));
    }
    py::object instance = newInstance(resultTypes().match);
    setField(instance, AtomMappingField, owned(PyDictProxy_New(atoms.ptr())).ptr());
    setField(instance, BondMappingField, bonds.ptr());
    setField(instance, FragmentsField, fragments.ptr());
    return instance;
}

const PythonResultTypes& registeredResultTypes() {
    const auto& types = resultTypes();
    if (types.match == nullptr || types.result == nullptr)
        throw std::runtime_error("Search result types have not been registered");
    return types;
}

/**
 * When turning a search's matches into Python objects stops, as native matchAll() stops
 * building them (see Budget::resultsStopped()): once the cancellation token, if any, is
 * cancelled, or past the time from the search's start that native matchAll() takes to build
 * them. A search left without either builds every match.
 */
struct ResultsDeadline {
    std::chrono::steady_clock::time_point at;
    bool bounded = false;
    std::shared_ptr<const smsd::CancellationToken> cancellation;

    /** Why building should stop now: cancellation first, then the deadline; Completed while neither applies. */
    smsd::SearchStatus stopReason() const {
        if (cancellation && cancellation->isCancelled()) return smsd::SearchStatus::Cancelled;
        if (bounded && std::chrono::steady_clock::now() >= at) return smsd::SearchStatus::TimedOut;
        return smsd::SearchStatus::Completed;
    }

    /** No deadline, with `cancellation`. */
    static ResultsDeadline cancelledBy(std::shared_ptr<const smsd::CancellationToken> cancellation) {
        ResultsDeadline deadline;
        deadline.cancellation = std::move(cancellation);
        return deadline;
    }

    /** The deadline `milliseconds` from `started`, unbounded when negative or too far away to count. */
    static ResultsDeadline after(int64_t milliseconds, std::chrono::steady_clock::time_point started,
                                 std::shared_ptr<const smsd::CancellationToken> cancellation) {
        ResultsDeadline deadline = cancelledBy(std::move(cancellation));
        // Compared in milliseconds, so a limit near the largest count does not overflow.
        const auto remaining = std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::steady_clock::time_point::max() - started).count();
        if (milliseconds >= 0 && milliseconds < remaining) {
            deadline.at = started + std::chrono::milliseconds(milliseconds);
            deadline.bounded = true;
        }
        return deadline;
    }

    template<class Query>
    static ResultsDeadline of(const Query& query, const smsd::MolGraph& target,
                              std::chrono::steady_clock::time_point started) {
        return after(query.resultsTimeMillis(target), started, query.cancellation());
    }
};

py::object pythonResult(const smsd::SearchResult& result, const AtomNumbering& query, const AtomNumbering& target,
                        const ResultsDeadline& deadline = ResultsDeadline{}) {
    const auto& types = registeredResultTypes();
    std::vector<py::object> built;
    built.reserve(result.matches.size());
    // As native matchAll() builds its matches, conversion counts towards the time limit and
    // honours cancellation: after about 4096 atoms and bonds the token and the clock are read,
    // and once the token is cancelled or the results deadline has passed the matches built so
    // far are returned, as a cancelled or timed-out result.
    auto stopped = smsd::SearchStatus::Completed;
    size_t work = 0;
    for (const auto& match : result.matches) {
        if (!built.empty() && work >= 4096) {
            work = 0;
            stopped = deadline.stopReason();
            if (stopped != smsd::SearchStatus::Completed) break;
        }
        built.push_back(pythonMatch(match, query, target));
        work += 1 + match.atomMapping.size() + match.bondMapping.size();
    }
    py::object matches = tupleOf(built);
    // A search that stopped by itself keeps its reason; one that completed or reached its limit
    // takes the reason its conversion stopped.
    const bool cut = stopped != smsd::SearchStatus::Completed;
    const auto resultStatus = cut && (result.status == smsd::SearchStatus::Completed
                                      || result.status == smsd::SearchStatus::LimitReached)
        ? stopped : result.status;
    const auto statusIndex = static_cast<size_t>(resultStatus);
    py::object status = statusIndex < 4 && types.status[statusIndex] != nullptr
        ? py::reinterpret_borrow<py::object>(types.status[statusIndex]) : py::cast(resultStatus);
    py::object instance = newInstance(types.result);
    setField(instance, MatchesField, matches.ptr());
    setField(instance, StatusField, status.ptr());
    setField(instance, OptimalField, result.optimal ? Py_True : Py_False);
    setField(instance, CompleteField, result.complete && !cut ? Py_True : Py_False);
    setField(instance, ScoreField, owned(PyFloat_FromDouble(result.score)).ptr());
    setField(instance, NodesVisitedField, owned(PyLong_FromUnsignedLongLong(result.nodesVisited)).ptr());
    setField(instance, ElapsedField, owned(PyLong_FromLongLong(result.elapsedMillis)).ptr());
    return instance;
}

py::object pythonResult(const smsd::SearchResult& result, py::handle queryIndices, py::handle targetIndices,
                        const ResultsDeadline& deadline = ResultsDeadline{}) {
    const AtomNumbering query(queryIndices), target(targetIndices);
    return pythonResult(result, query, target, deadline);
}

/**
 * Fail unless `type` is a class whose instances object.__new__ creates and whose dataclass
 * fields are exactly the names from `first` up to `last`, which the native results set.
 */
void requireResultClass(py::handle type, int first, int last, const char* role) {
    const std::string name(role);
    if (!PyType_Check(type.ptr())) throw py::type_error(name + " must be a class");
    if (reinterpret_cast<PyTypeObject*>(type.ptr())->tp_new != PyBaseObject_Type.tp_new)
        throw py::type_error(name + " must be created by object.__new__");
    const py::object fields = py::getattr(type, "__dataclass_fields__", py::none());
    if (!PyDict_Check(fields.ptr())) throw py::type_error(name + " must be a dataclass");
    bool same = PyDict_GET_SIZE(fields.ptr()) == last - first;
    for (int field = first; same && field < last; ++field) {
        const int contains = PyDict_Contains(fields.ptr(), resultTypes().names[field]);
        if (contains < 0) throw py::error_already_set();
        same = contains == 1;
    }
    if (!same) throw py::type_error(name + " fields differ from those the native results set");
}

void registerResultTypes(py::handle match, py::handle result) {
    requireResultClass(match, AtomMappingField, MatchesField, "Match");
    requireResultClass(result, MatchesField, ResultFieldCount, "SearchResult");
    auto& types = resultTypes();
    replaceReference(types.match, py::reinterpret_borrow<py::object>(match).release().ptr());
    replaceReference(types.result, py::reinterpret_borrow<py::object>(result).release().ptr());
}

/** Search one target and return the registered Python result directly. */
template<class Query>
py::object matchPython(const Query& query, const AtomNumbering& queryNumbering, const smsd::MolGraph& target,
                       const AtomNumbering& targetNumbering) {
    registeredResultTypes();
    const auto started = std::chrono::steady_clock::now();
    smsd::SearchResult result;
    {
        py::gil_scoped_release release;
        result = query.match(target);
    }
    return pythonResult(result, queryNumbering, targetNumbering, ResultsDeadline::of(query, target, started));
}

template<class Query>
py::object matchAllPython(const Query& query, const AtomNumbering& queryNumbering, const smsd::MolGraph& target,
                          int limit, const AtomNumbering& targetNumbering) {
    if (limit < 1) throw std::invalid_argument("Match limit must be positive");
    registeredResultTypes();
    const auto started = std::chrono::steady_clock::now();
    smsd::SearchResult result;
    {
        py::gil_scoped_release release;
        result = query.matchAll(target, static_cast<size_t>(limit));
    }
    return pythonResult(result, queryNumbering, targetNumbering, ResultsDeadline::of(query, target, started));
}

/** Search (graph, atom numbering) pairs and return one Python result per pair. */
template<class Query>
py::list matchManyPython(const Query& query, const AtomNumbering& queryNumbering, const py::sequence& targets,
                         int numThreads, const smsd::batch::TargetCorpus* corpus = nullptr) {
    if (numThreads < 0) throw std::invalid_argument("Thread count must be non-negative");
    registeredResultTypes();
    PythonGraphView view;
    std::vector<AtomNumbering> numbering;
    const auto count = static_cast<size_t>(targets.size());
    view.graphs.reserve(count);
    view.owners.reserve(count);
    numbering.reserve(count);
    // Each pair is held by an owned reference, as a sequence may create its items afresh.
    for (py::object target : targets) {
        if (!PyTuple_Check(target.ptr()) || PyTuple_GET_SIZE(target.ptr()) != 2)
            throw py::type_error("Each target must be a (graph, atom numbering) pair");
        view.graphs.push_back(py::handle(PyTuple_GET_ITEM(target.ptr(), 0)).cast<const smsd::MolGraph&>());
        numbering.emplace_back(PyTuple_GET_ITEM(target.ptr(), 1));
        view.owners.push_back(std::move(target));
    }
    const auto results = matchView(query, view, numThreads, corpus);
    std::vector<py::object> built;
    built.reserve(results.size());
    // Each search took its own time, so only cancellation stops turning the results over.
    const auto deadline = ResultsDeadline::cancelledBy(query.cancellation());
    for (size_t index = 0; index < results.size(); ++index)
        built.push_back(pythonResult(results[index], queryNumbering, numbering[index], deadline));
    py::list converted(built.size());
    for (size_t index = 0; index < built.size(); ++index)
        PyList_SET_ITEM(converted.ptr(), static_cast<Py_ssize_t>(index), built[index].release().ptr());
    return converted;
}

/** Corpus decisions with the same atom-numbering conversion as ordinary batch searches. */
template<class Query>
py::list matchCorpusPython(const Query& query, const smsd::batch::TargetCorpus& corpus,
    const py::sequence& targets, int minimumAtoms, bool gpuPrefilter, int numThreads, py::handle queryIndices) {
    if (numThreads < 0) throw std::invalid_argument("Thread count must be non-negative");
    registeredResultTypes();
    const AtomNumbering queryNumbering(queryIndices);
    std::vector<AtomNumbering> numbering;
    std::vector<py::object> owners;
    numbering.reserve(static_cast<size_t>(targets.size()));
    for (py::object target : targets) {
        if (!PyTuple_Check(target.ptr()) || PyTuple_GET_SIZE(target.ptr()) != 2)
            throw py::type_error("Each target must be a (graph, atom numbering) pair");
        numbering.emplace_back(PyTuple_GET_ITEM(target.ptr(), 1));
        owners.push_back(std::move(target));
    }
    std::vector<smsd::CorpusSearchResult> outcomes;
    {
        py::gil_scoped_release release;
        const int threads = batchThreads(numThreads);
        if constexpr (std::is_same<Query, smsd::MCS>::value)
            outcomes = query.matchCorpus(corpus, minimumAtoms, gpuPrefilter, threads);
        else
            outcomes = query.matchCorpus(corpus, gpuPrefilter, threads);
    }
    if (outcomes.size() != numbering.size())
        throw std::runtime_error("Corpus changed while its atom numbering was captured");
    const auto deadline = ResultsDeadline::cancelledBy(query.cancellation());
    py::list converted;
    for (const auto& outcome : outcomes) {
        py::object result = outcome.result
            ? pythonResult(*outcome.result, queryNumbering, numbering[outcome.targetIndex], deadline) : py::none();
        py::object bound = outcome.atomUpperBound < 0 ? py::none() : py::cast(outcome.atomUpperBound);
        converted.append(py::make_tuple(outcome.targetIndex, bound, outcome.provenBelowThreshold,
                                        outcome.thresholdReached, result));
    }
    return converted;
}

/** A prepared query with the caller's numbering of its atoms, checked once when it is made. */
template<class Query>
struct NumberedQuery {
    Query query;
    AtomNumbering numbering;
};

/**
 * The private search entry points of a prepared query class: those that take the query
 * numbering on every call, and a numbered query whose calls take only their targets.
 */
template<class Query>
void bindPythonSearches(py::class_<Query>& queries) {
    using Numbered = NumberedQuery<Query>;
    py::class_<Numbered>(queries, "_NumberedQuery")
        .def("_match", [](const Numbered& numbered, const smsd::MolGraph& target) {
            return matchPython(numbered.query, numbered.numbering, target, AtomNumbering());
        }, py::arg("target"), "Search a graph whose atoms keep their native indices.")
        .def("_match_all", [](const Numbered& numbered, const smsd::MolGraph& target, int limit) {
            return matchAllPython(numbered.query, numbered.numbering, target, limit, AtomNumbering());
        }, py::arg("target"), py::arg("limit"), "Enumerate matches in a graph whose atoms keep their native indices.")
        .def("_match_many", [](const Numbered& numbered, const py::sequence& targets, int numThreads) {
            return matchManyPython(numbered.query, numbered.numbering, targets, numThreads);
        }, py::arg("targets"), py::arg("num_threads"), "Search (graph, atom numbering) pairs, one result per pair.");
    queries
        .def("_with_numbering", [](const Query& query, py::handle queryIndices) {
            return Numbered{query, AtomNumbering(queryIndices)};
        }, py::arg("query_indices"), "This query with its atom numbering, or None for native indices.")
        .def("_match", [](const Query& query, const smsd::MolGraph& target, py::handle queryIndices,
                          py::handle targetIndices) {
            const AtomNumbering queryNumbering(queryIndices), targetNumbering(targetIndices);
            return matchPython(query, queryNumbering, target, targetNumbering);
        })
        .def("_match_all", [](const Query& query, const smsd::MolGraph& target, int limit, py::handle queryIndices,
                              py::handle targetIndices) {
            const AtomNumbering queryNumbering(queryIndices), targetNumbering(targetIndices);
            return matchAllPython(query, queryNumbering, target, limit, targetNumbering);
        })
        .def("_match_many", [](const Query& query, const py::sequence& targets, int numThreads,
                               py::handle queryIndices) {
            return matchManyPython(query, AtomNumbering(queryIndices), targets, numThreads);
        });
}

} // anonymous namespace

// ============================================================================
// Public fingerprint functions
// ============================================================================

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

/** Issue a DeprecationWarning for a deprecated attribute, as the Python package does. */
static void warnDeprecated(const char* message) {
    if (PyErr_WarnEx(PyExc_DeprecationWarning, message, 1) < 0) throw py::error_already_set();
}

// The deprecated field stays readable and writable from Python, with a warning on each use.
static constexpr const char* kAromaticityModelDeprecated =
    "ChemOptions.aromaticity_model is deprecated: no search reads it; pass the model to "
    "MolGraphBuilder.build() or perceive_aromaticity()";
static smsd::AromaticityModel deprecatedAromaticityModel(const smsd::ChemOptions& options) {
    warnDeprecated(kAromaticityModelDeprecated);
    return options.aromaticityModel;
}
static void setDeprecatedAromaticityModel(smsd::ChemOptions& options, smsd::AromaticityModel model) {
    warnDeprecated(kAromaticityModelDeprecated);
    options.aromaticityModel = model;
}

// smsd::batch::fingerprintQuality under the keys of the Java analyzeFpQuality.
static py::dict analyzeFpQuality(const std::vector<int>& fp, int fpSize = 2048) {
    const auto quality = smsd::batch::fingerprintQuality(fp, fpSize);
    py::dict result;
    result["set_bits"] = quality.setBits;
    result["fp_size"] = quality.fpSize;
    result["density"] = quality.density;
    result["expected_collision_rate"] = quality.expectedCollisionRate;
    result["is_saturated"] = quality.saturated;
    return result;
}

// ============================================================================
// Module definition
// ============================================================================

PYBIND11_MODULE(_smsd, m) {
    m.doc() = "SMSD -- Substructure & MCS search for chemical graphs";
    // The package compares this with its own version when it is imported.
    m.attr("__version__") = SMSD_VERSION_INFO;

    // Parser messages quote the input byte at fault, which need not be whole UTF-8; such bytes are
    // shown as escapes, so the error stays a ValueError with the parser's message.
    py::register_local_exception_translator([](std::exception_ptr pointer) {
        try {
            if (pointer) std::rethrow_exception(pointer);
        } catch (const std::invalid_argument& error) {
            const char* text = error.what();
            PyObject* message = PyUnicode_DecodeUTF8(text, static_cast<Py_ssize_t>(std::strlen(text)), "backslashreplace");
            if (message == nullptr) return;
            PyErr_SetObject(PyExc_ValueError, message);
            Py_DECREF(message);
        }
    });
    // An unavailable exact ring basis is a RuntimeError of its own type, as in C++ and Java; a
    // batch failure keeps the single-call failure as its cause.
    static PyObject* ringBasisUnavailable = nullptr;
    ringBasisUnavailable = py::register_exception<smsd::detail_rf::RingBasisUnavailable>(
        m, "RingBasisUnavailable", PyExc_RuntimeError).ptr();
    py::register_local_exception_translator([](std::exception_ptr pointer) {
        try {
            if (pointer) std::rethrow_exception(pointer);
        } catch (const smsd::detail_rf::RingBasisUnavailable& error) {
            py::detail::handle_nested_exception(error, pointer);
            if (PyErr_Occurred()) py::raise_from(ringBasisUnavailable, error.what());
            else PyErr_SetString(ringBasisUnavailable, error.what());
        }
    });

    py::enum_<smsd::ChemOptions::BondOrderMode>(m, "BondOrderMode",
          R"doc(Bond-order matching mode.


Example::

    >>> import smsd
    >>> chemistry = smsd.ChemOptions().with_bond_order_mode(smsd.BondOrderMode.STRICT)
    >>> assert chemistry.match_bond_order == smsd.BondOrderMode.STRICT
)doc")
        .value("STRICT", smsd::ChemOptions::BondOrderMode::STRICT)
        .value("LOOSE",  smsd::ChemOptions::BondOrderMode::LOOSE)
        .value("ANY",    smsd::ChemOptions::BondOrderMode::ANY)
        .export_values();

    py::enum_<smsd::ChemOptions::AromaticityMode>(m, "AromaticityMode",
          R"doc(Aromatic bond compatibility mode during matching.


Example::

    >>> import smsd
    >>> chemistry = smsd.ChemOptions().with_aromaticity_mode(smsd.AromaticityMode.STRICT)
    >>> assert chemistry.aromaticity_mode == smsd.AromaticityMode.STRICT
)doc")
        .value("STRICT",   smsd::ChemOptions::AromaticityMode::STRICT)
        .value("FLEXIBLE", smsd::ChemOptions::AromaticityMode::FLEXIBLE)
        .export_values();

    py::enum_<smsd::AromaticityModel>(m, "AromaticityModel",
          R"doc(Model used when perceiving graph aromaticity.


Example::

    >>> import smsd
    >>> graph = smsd.perceive_aromaticity("C1=CC=CC=C1", model=smsd.AromaticityModel.DAYLIGHT_LIKE)
    >>> assert all(graph.aromatic)
)doc")
        .value("DAYLIGHT_LIKE", smsd::AromaticityModel::DAYLIGHT_LIKE)
        .export_values();

    py::enum_<smsd::ChemOptions::RingFusionMode>(m, "RingFusionMode",
          R"doc(Ring-fusion matching restriction.


Example::

    >>> import smsd
    >>> chemistry = smsd.ChemOptions().with_ring_fusion_mode(smsd.RingFusionMode.STRICT)
    >>> assert chemistry.ring_fusion_mode == smsd.RingFusionMode.STRICT
)doc")
        .value("IGNORE",     smsd::ChemOptions::RingFusionMode::IGNORE)
        .value("PERMISSIVE", smsd::ChemOptions::RingFusionMode::PERMISSIVE)
        .value("STRICT",     smsd::ChemOptions::RingFusionMode::STRICT)
        .export_values();

    py::enum_<smsd::ChemOptions::Solvent>(m, "Solvent",
          R"doc(Solvent used by the optional tautomer-confidence heuristic.


Example::

    >>> import smsd
    >>> chemistry = smsd.ChemOptions().with_solvent(smsd.Solvent.AQUEOUS)
    >>> assert chemistry.solvent == smsd.Solvent.AQUEOUS
)doc")
        .value("AQUEOUS",       smsd::ChemOptions::Solvent::AQUEOUS)
        .value("DMSO",          smsd::ChemOptions::Solvent::DMSO)
        .value("METHANOL",      smsd::ChemOptions::Solvent::METHANOL)
        .value("CHLOROFORM",    smsd::ChemOptions::Solvent::CHLOROFORM)
        .value("ACETONITRILE",  smsd::ChemOptions::Solvent::ACETONITRILE)
        .value("DIETHYL_ETHER", smsd::ChemOptions::Solvent::DIETHYL_ETHER)
        .export_values();

    py::class_<smsd::ChemOptions>(m, "ChemOptions",
          R"doc(Atom, bond, ring and stereo matching settings.


Example::

    >>> import smsd
    >>> chemistry = smsd.ChemOptions.profile("strict")
    >>> options = smsd.MCS.Options()
    >>> options.chemistry = chemistry
    >>> assert smsd.MCS.find("CCO").with_options(options).match("CCO").best.atom_count == 3
)doc")
        .def(py::init<>())
        .def(py::init<const smsd::ChemOptions&>(), py::arg("options"), "An independent copy of options")
        .def("__copy__", [](const smsd::ChemOptions& o) { return smsd::ChemOptions(o); })
        .def("__deepcopy__", [](const smsd::ChemOptions& o, const py::dict&) { return smsd::ChemOptions(o); },
             py::arg("memo"))
        // Atom matching
        .def_readwrite("match_atom_type",       &smsd::ChemOptions::matchAtomType)
        .def_readwrite("match_formal_charge",   &smsd::ChemOptions::matchFormalCharge)
        .def_readwrite("match_isotope",         &smsd::ChemOptions::matchIsotope)
        .def("with_match_isotope", &smsd::ChemOptions::withMatchIsotope, py::arg("on"),
             py::return_value_policy::reference_internal, "Set isotope matching and return these options"
          R"doc(

Example::

    >>> import smsd
    >>> options = smsd.ChemOptions()
    >>> assert options.with_match_isotope(True) is options
    >>> assert options.match_isotope == True
)doc")
        .def_readwrite("use_chirality",         &smsd::ChemOptions::useChirality)
        .def_readwrite("use_bond_stereo",       &smsd::ChemOptions::useBondStereo)
        .def_readwrite("ring_matches_ring_only",&smsd::ChemOptions::ringMatchesRingOnly)
        // Ring completeness
        .def_readwrite("complete_rings_only",   &smsd::ChemOptions::completeRingsOnly)
        .def("with_complete_rings_only", &smsd::ChemOptions::withCompleteRingsOnly, py::arg("on"),
             py::return_value_policy::reference_internal, "Set complete-ring matching and return these options"
          R"doc(

Example::

    >>> import smsd
    >>> options = smsd.ChemOptions()
    >>> assert options.with_complete_rings_only(True) is options
    >>> assert options.complete_rings_only == True
)doc")
        // Tautomer
        .def_readwrite("tautomer_aware",        &smsd::ChemOptions::tautomerAware)
        .def("with_tautomer_aware", &smsd::ChemOptions::withTautomerAware, py::arg("on"),
             py::return_value_policy::reference_internal, "Set tautomer-aware matching and return these options"
          R"doc(

Example::

    >>> import smsd
    >>> options = smsd.ChemOptions()
    >>> assert options.with_tautomer_aware(True) is options
    >>> assert options.tautomer_aware == True
)doc")
        // Solvent
        .def_readwrite("solvent",               &smsd::ChemOptions::solvent)
        .def("with_solvent", &smsd::ChemOptions::withSolvent, py::arg("solvent"),
             py::return_value_policy::reference_internal,
             "Set the solvent for tautomer weight corrections and return these options"
          R"doc(

Example::

    >>> import smsd
    >>> options = smsd.ChemOptions()
    >>> assert options.with_solvent(smsd.Solvent.METHANOL) is options
    >>> assert options.solvent == smsd.Solvent.METHANOL
)doc")
        // Bond matching
        .def_readwrite("match_bond_order",      &smsd::ChemOptions::matchBondOrder)
        .def("with_bond_order_mode", &smsd::ChemOptions::withBondOrderMode, py::arg("mode"),
             py::return_value_policy::reference_internal, "Set the bond order mode and return these options"
          R"doc(

Example::

    >>> import smsd
    >>> options = smsd.ChemOptions()
    >>> assert options.with_bond_order_mode(smsd.BondOrderMode.ANY) is options
    >>> assert options.match_bond_order == smsd.BondOrderMode.ANY
)doc")
        .def_readwrite("aromaticity_mode",      &smsd::ChemOptions::aromaticityMode)
        .def("with_aromaticity_mode", &smsd::ChemOptions::withAromaticityMode, py::arg("mode"),
             py::return_value_policy::reference_internal, "Set the aromaticity mode and return these options"
          R"doc(

Example::

    >>> import smsd
    >>> options = smsd.ChemOptions()
    >>> assert options.with_aromaticity_mode(smsd.AromaticityMode.STRICT) is options
    >>> assert options.aromaticity_mode == smsd.AromaticityMode.STRICT
)doc")
        .def_property("aromaticity_model", &deprecatedAromaticityModel, &setDeprecatedAromaticityModel,
             "Deprecated: no search reads this field. Pass the model to MolGraphBuilder.build() or "
             "perceive_aromaticity(); reading or writing it warns with DeprecationWarning.")
        // pH
        .def_readwrite("pH",                    &smsd::ChemOptions::pH)
        .def("with_ph", &smsd::ChemOptions::withPH, py::arg("pH"),
             py::return_value_policy::reference_internal,
             "Set a finite pH for the tautomer confidence helper, clamped to [0, 14], and return these options"
          R"doc(

Example::

    >>> import smsd
    >>> options = smsd.ChemOptions()
    >>> assert options.with_ph(7.4) is options
    >>> assert options.pH == 7.4
)doc")
        // Ring fusion
        .def_readwrite("ring_fusion_mode",      &smsd::ChemOptions::ringFusionMode)
        .def("with_ring_fusion_mode", &smsd::ChemOptions::withRingFusionMode, py::arg("mode"),
             py::return_value_policy::reference_internal, "Set the ring fusion mode and return these options"
          R"doc(

Example::

    >>> import smsd
    >>> options = smsd.ChemOptions()
    >>> assert options.with_ring_fusion_mode(smsd.RingFusionMode.STRICT) is options
    >>> assert options.ring_fusion_mode == smsd.RingFusionMode.STRICT
)doc")
        // Named profiles
        .def_static("tautomer_profile", &smsd::ChemOptions::tautomerProfile,
                     "ChemOptions for tautomer-aware matching with loose bond orders; otherwise the defaults"
          R"doc(

Example::

    >>> import smsd
    >>> assert smsd.ChemOptions.tautomer_profile().tautomer_aware
)doc")
        .def_static("profile", &smsd::ChemOptions::profile,
                     py::arg("name") = "default",
                     "ChemOptions from a named profile, ignoring case: 'default', 'strict', 'pharma', "
                     "'compat-fmcs', 'compat-substruct' or 'reaction'"
          R"doc(

Example::

    >>> import smsd
    >>> options = smsd.ChemOptions.profile("default")
    >>> assert options.match_atom_type and not options.match_formal_charge
)doc")
        .def("__repr__", [](const smsd::ChemOptions& o) {
            auto flag = [](bool value) { return std::string(value ? "True" : "False"); };
            auto bondOrder = [](smsd::ChemOptions::BondOrderMode mode) {
                return mode == smsd::ChemOptions::BondOrderMode::STRICT ? "STRICT"
                     : mode == smsd::ChemOptions::BondOrderMode::LOOSE ? "LOOSE" : "ANY";
            };
            auto ringFusion = [](smsd::ChemOptions::RingFusionMode mode) {
                return mode == smsd::ChemOptions::RingFusionMode::IGNORE ? "IGNORE"
                     : mode == smsd::ChemOptions::RingFusionMode::PERMISSIVE ? "PERMISSIVE" : "STRICT";
            };
            auto solvent = [](smsd::ChemOptions::Solvent value) {
                switch (value) {
                    case smsd::ChemOptions::Solvent::AQUEOUS: return "AQUEOUS";
                    case smsd::ChemOptions::Solvent::DMSO: return "DMSO";
                    case smsd::ChemOptions::Solvent::METHANOL: return "METHANOL";
                    case smsd::ChemOptions::Solvent::CHLOROFORM: return "CHLOROFORM";
                    case smsd::ChemOptions::Solvent::ACETONITRILE: return "ACETONITRILE";
                    default: return "DIETHYL_ETHER";
                }
            };
            return "<ChemOptions match_atom_type=" + flag(o.matchAtomType)
                 + " match_formal_charge=" + flag(o.matchFormalCharge)
                 + " match_isotope=" + flag(o.matchIsotope)
                 + " use_chirality=" + flag(o.useChirality)
                 + " use_bond_stereo=" + flag(o.useBondStereo)
                 + " ring_matches_ring_only=" + flag(o.ringMatchesRingOnly)
                 + " complete_rings_only=" + flag(o.completeRingsOnly)
                 + " ring_fusion_mode=" + ringFusion(o.ringFusionMode)
                 + " match_bond_order=" + bondOrder(o.matchBondOrder)
                 + " aromaticity_mode=" + (o.aromaticityMode == smsd::ChemOptions::AromaticityMode::STRICT ? "STRICT" : "FLEXIBLE")
                 + " tautomer_aware=" + flag(o.tautomerAware)
                 + " pH=" + py::str(py::float_(o.pH)).cast<std::string>()
                 + " solvent=" + solvent(o.solvent) + ">";
        });

    py::class_<PropertyView>(m, "_Properties")
        .def("__len__", [](const PropertyView& view) { return view.values().size(); })
        .def("__contains__", [](const PropertyView& view, const std::string& key) {
            return view.values().count(key) != 0;
        })
        .def("__getitem__", [](const PropertyView& view, const std::string& key) {
            auto entry = view.values().find(key);
            if (entry == view.values().end()) throw py::key_error(key);
            return entry->second;
        })
        .def("__setitem__", [](PropertyView& view, const std::string& key, const std::string& value) {
            view.values()[key] = value;
        })
        .def("__delitem__", [](PropertyView& view, const std::string& key) {
            if (view.values().erase(key) == 0) throw py::key_error(key);
        })
        .def("__iter__", [](const PropertyView& view) { return py::iter(py::cast(view.values())); })
        .def("keys", [](const PropertyView& view) { return py::cast(view.values()).attr("keys")(); })
        .def("items", [](const PropertyView& view) { return py::cast(view.values()).attr("items")(); })
        .def("values", [](const PropertyView& view) { return py::cast(view.values()).attr("values")(); })
        .def("get", [](const PropertyView& view, const std::string& key, py::object fallback) -> py::object {
            auto entry = view.values().find(key);
            return entry == view.values().end() ? fallback : py::cast(entry->second);
        }, py::arg("key"), py::arg("default") = py::none())
        .def("update", [](PropertyView& view, const py::dict& values) {
            const auto entries = values.cast<std::map<std::string, std::string>>();
            for (const auto& entry : entries) view.values()[entry.first] = entry.second;
        }, py::arg("values"))
        .def("clear", [](PropertyView& view) { view.values().clear(); })
        .def("copy", [](const PropertyView& view) { return py::cast(view.values()); })
        .def("__repr__", [](const PropertyView& view) { return py::repr(py::cast(view.values())); })
        .def("__eq__", [](const PropertyView& view, py::object other) {
            py::object values = py::cast(view.values());
            return values.equal(py::hasattr(other, "keys") ? py::object(py::dict(other)) : other);
        }, py::is_operator());

    py::class_<smsd::MolGraph>(m, "MolGraph",
          R"doc(Native molecular graph with represented atom and bond properties.


Example::

    >>> import smsd
    >>> graph = smsd.parse_smiles("CCO")
    >>> assert graph.n == len(graph) == 3 and graph.atomic_num == [6, 6, 8]
    >>> assert graph.has_bond(0, 1) and graph.bond_order(1, 2) == 1
)doc")
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
        .def_property("properties", [](py::object owner) { return PropertyView{owner}; },
            [](smsd::MolGraph& graph, const std::map<std::string, std::string>& values) {
                graph.properties = values;
            })
        .def_readonly("label",         &smsd::MolGraph::label)
        .def_readonly("tautomer_class",&smsd::MolGraph::tautomerClass)
        .def_property_readonly("morgan_rank", [](const smsd::MolGraph& graph) {
            const std::lock_guard<std::recursive_mutex> cacheLock(graph.lazyCacheMutex());
            return graph.morganRank;
        })
        // Computed on first access, as the C++ getters compute them; morgan_rank is filled by the same step.
        .def_property_readonly("canonical_label", [](const smsd::MolGraph& graph) {
            const std::lock_guard<std::recursive_mutex> cacheLock(graph.lazyCacheMutex());
            return graph.getCanonicalLabeling();
        })
        .def_property_readonly("orbit", [](const smsd::MolGraph& graph) {
            const std::lock_guard<std::recursive_mutex> cacheLock(graph.lazyCacheMutex());
            return graph.getOrbit();
        })
        .def_readonly("neighbors",     &smsd::MolGraph::neighbors)
        .def("bonds", [](const smsd::MolGraph& g) {
                 std::vector<std::tuple<int,int,int>> result;
                 for (int i = 0; i < g.n; ++i)
                     for (int j : g.neighbors[i])
                         if (j > i) result.emplace_back(i, j, g.bondOrder(i, j));
                 return result;
             }, "Return all bonds as (atom_i, atom_j, order) tuples, 0-indexed"
          R"doc(

Example::

    >>> import smsd
    >>> graph = smsd.parse_smiles("CCO")
    >>> assert graph.bonds() == [(0, 1, 1), (1, 2, 1)]
)doc")
        .def("bond_order", [](const smsd::MolGraph& g, int i, int j) {
                 requireAtomPair(g, i, j);
                 return g.bondOrder(i, j);
             }, py::arg("i"), py::arg("j"),
             "Get bond order between atoms i and j (0 if no bond)"
          R"doc(

Example::

    >>> import smsd
    >>> graph = smsd.parse_smiles("CC=O")
    >>> assert graph.bond_order(1, 2) == 2 and graph.bond_order(0, 2) == 0
)doc")
        .def("bond_in_ring", [](const smsd::MolGraph& g, int i, int j) {
                 requireAtomPair(g, i, j);
                 return g.bondInRing(i, j);
             }, py::arg("i"), py::arg("j"),
          R"doc(Whether a represented bond belongs to a ring.


Example::

    >>> import smsd
    >>> graph = smsd.parse_smiles("C1CC1")
    >>> assert graph.bond_in_ring(0, 1)
)doc")
        .def("bond_aromatic", [](const smsd::MolGraph& g, int i, int j) {
                 requireAtomPair(g, i, j);
                 return g.bondAromatic(i, j);
             }, py::arg("i"), py::arg("j"),
          R"doc(Whether a represented bond is marked aromatic.


Example::

    >>> import smsd
    >>> graph = smsd.parse_smiles("c1ccccc1")
    >>> assert graph.bond_aromatic(0, 1)
)doc")
        .def("ensure_canonical", [](const smsd::MolGraph& g) {
                 g.ensureCanonical();
             }, "Compute and cache Morgan ranks, canonical labels, and orbits"
          R"doc(

Example::

    >>> import smsd
    >>> graph = smsd.parse_smiles("CCO")
    >>> graph.ensure_canonical()
    >>> assert sorted(graph.canonical_label) == [0, 1, 2]
)doc")
        .def("get_automorphism_generators", [](const smsd::MolGraph& g) {
                 return g.getAutomorphismGenerators();
             }, "Return automorphism generators as list of permutation lists.\n"
                "Each generator is a list of length n where gen[i] is the image\n"
                "of atom i. The full automorphism group is the closure of these\n"
                "generators."
          R"doc(

Example::

    >>> import smsd
    >>> graph = smsd.parse_smiles("CC")
    >>> generators = graph.get_automorphism_generators()
    >>> assert generators and all(sorted(generator) == [0, 1] for generator in generators)
)doc")
        .def("automorphism_generators_truncated", [](const smsd::MolGraph& g) {
                 return g.automorphismGeneratorsTruncated();
             }, "True if generators were capped (group may be incomplete)"
          R"doc(

Example::

    >>> import smsd
    >>> graph = smsd.parse_smiles("CC")
    >>> _ = graph.get_automorphism_generators()
    >>> assert not graph.automorphism_generators_truncated()
)doc")
        .def("get_canonical_hash", [](const smsd::MolGraph& g) -> uint64_t {
                 return g.getCanonicalHash();
             }, "Return the canonical hash of this molecule. "
                "Covers connectivity, bond order, aromaticity, ring membership, "
                "formal charge, mass number and stereo (tetrahedral chirality and "
                "E/Z double-bond configuration). The hash is a compact screening key, "
                "not a proof of identity or difference: use a validated graph mapping "
                "when identity matters."
          R"doc(

Example::

    >>> import smsd
    >>> graph = smsd.parse_smiles("CCO")
    >>> assert graph.get_canonical_hash() == smsd.parse_smiles("OCC").get_canonical_hash()
)doc")
        .def("ensure_ring_counts", [](const smsd::MolGraph& g) {
                 g.ensureRingCounts();
             }, "Compute and cache legacy bounded Horton-candidate counts per atom, "
                "using sorted neighbour order"
          R"doc(

Example::

    >>> import smsd
    >>> graph = smsd.parse_smiles("c1ccccc1")
    >>> graph.ensure_ring_counts()
    >>> assert graph.num_rings() == 1
)doc")
        .def("ring_system_of", [](const smsd::MolGraph& g, int atom) {
                 return g.ringSystemOf(atom);
             }, py::arg("atom"), "Ring system ID for an atom (-1 if acyclic)"
          R"doc(

Example::

    >>> import smsd
    >>> graph = smsd.parse_smiles("CCc1ccccc1")
    >>> assert graph.ring_system_of(0) == -1 and graph.ring_system_of(2) >= 0
)doc")
        .def("ring_system_sig", [](const smsd::MolGraph& g, int rsId) {
                 return g.ringSystemSig(rsId);
             }, py::arg("rs_id"), "Signature of a ring system; identical ring systems share it"
          R"doc(

Example::

    >>> import smsd
    >>> graph = smsd.parse_smiles("c1ccccc1")
    >>> system = graph.ring_system_of(0)
    >>> assert graph.ring_system_sig(system) == graph.ring_system_sig(graph.ring_system_of(1))
)doc")
        .def("ring_system_count", [](const smsd::MolGraph& g) {
                 return g.ringSystemCount();
             }, "Number of distinct ring systems"
          R"doc(

Example::

    >>> import smsd
    >>> assert smsd.parse_smiles("c1ccccc1.C1CC1").ring_system_count() == 2
)doc")
        .def("prewarm", [](const smsd::MolGraph& g) {
                 smsd::batch::detail::prewarmGraph(g);
             },
             py::call_guard<py::gil_scoped_release>(),
             "Compute and cache common lazy invariants used by matching and batch APIs"
          R"doc(

Example::

    >>> import smsd
    >>> graph = smsd.parse_smiles("CCO")
    >>> graph.prewarm()
    >>> assert smsd.MCS.find(graph).match("CCO").best.atom_count == 3
)doc")
        .def("num_rings", [](const smsd::MolGraph& g) {
                 g.ensureRingCounts();
                 auto sssr = smsd::computeSSSR(g);
                 return static_cast<int>(sssr.size());
             },
             "Return the SSSR ring count for this molecule"
          R"doc(

Example::

    >>> import smsd
    >>> assert smsd.parse_smiles("c1ccccc1").num_rings() == 1
)doc")
        .def("perceive_aromaticity", [](smsd::MolGraph& g, smsd::AromaticityModel model) -> smsd::MolGraph& {
                 g.perceiveAromaticity(model);
                 return g;
             },
             py::return_value_policy::reference_internal,
             py::arg("model") = smsd::AromaticityModel::DAYLIGHT_LIKE,
             py::call_guard<py::gil_scoped_release>(),
             "Recompute ring membership and aromaticity using the native model"
          R"doc(

Example::

    >>> import smsd
    >>> graph = smsd.parse_smiles("C1=CC=CC=C1")
    >>> assert graph.perceive_aromaticity() is graph and all(graph.aromatic)
)doc")
        .def("kekulize", [](smsd::MolGraph& g) -> smsd::MolGraph& {
                 if (!g.kekulize())
                     throw py::value_error("kekulize could not assign a valid Kekule form for this aromatic system");
                 return g;
             },
             py::return_value_policy::reference_internal,
             py::call_guard<py::gil_scoped_release>(),
             "Convert aromatic atoms and bonds into an explicit Kekule assignment"
          R"doc(

Example::

    >>> import smsd
    >>> graph = smsd.parse_smiles("c1ccccc1")
    >>> assert graph.kekulize() is graph and not any(graph.aromatic)
)doc")
        .def("dearomatize", [](smsd::MolGraph& g) -> smsd::MolGraph& {
                 if (!g.dearomatize())
                     throw py::value_error("dearomatize could not assign a valid non-aromatic Kekule form for this aromatic system");
                 return g;
             },
             py::return_value_policy::reference_internal,
             py::call_guard<py::gil_scoped_release>(),
             "Alias of kekulize(): remove aromatic flags via an explicit Kekule assignment"
          R"doc(

Example::

    >>> import smsd
    >>> graph = smsd.parse_smiles("c1ccccc1")
    >>> assert graph.dearomatize() is graph and not any(graph.aromatic)
)doc")
        .def("has_bond", [](const smsd::MolGraph& g, int i, int j) {
                 requireAtomPair(g, i, j);
                 return g.hasBond(i, j);
             }, py::arg("i"), py::arg("j"),
          R"doc(Whether two valid atom indices are adjacent.


Example::

    >>> import smsd
    >>> graph = smsd.parse_smiles("CCO")
    >>> assert graph.has_bond(0, 1) and not graph.has_bond(0, 2)
)doc")
        .def("__repr__", [](const smsd::MolGraph& g) {
            return "<MolGraph n=" + std::to_string(g.n) + ">";
        })
        .def("__len__", [](const smsd::MolGraph& g) { return g.n; });

    py::class_<smsd::MolGraph::Builder>(m, "MolGraphBuilder",
          R"doc(Construct a molecular graph from explicit atom and adjacency data.


Example::

    >>> import smsd
    >>> graph = (smsd.MolGraphBuilder().atom_count(2).atomic_numbers([6, 8])
    ...          .neighbors([[1], [0]]).bond_orders([[1], [1]]).hydrogen_counts([3, 1]).build())
    >>> assert graph.n == 2 and graph.atomic_num == [6, 8]
)doc")
        .def(py::init<>())
        .def("atom_count", [](smsd::MolGraph::Builder& b, int n) -> smsd::MolGraph::Builder& {
            return b.atomCount(n);
        }, py::arg("n"), py::return_value_policy::reference_internal,
             "Set the number of atoms"
          R"doc(

Example::

    >>> import smsd
    >>> builder = smsd.MolGraphBuilder().atom_count(2).atomic_numbers([6, 8]).neighbors([[1], [0]]).bond_orders([[1], [1]])
    >>> _ = builder.atom_count(2)
    >>> graph = builder.build(perceive_aromaticity=False)
    >>> assert graph.n == 2
)doc")
        .def("atomic_numbers", [](smsd::MolGraph::Builder& b, std::vector<int> v) -> smsd::MolGraph::Builder& {
            return b.atomicNumbers(std::move(v));
        }, py::arg("nums"), py::return_value_policy::reference_internal,
             "Set atomic numbers for all atoms"
          R"doc(

Example::

    >>> import smsd
    >>> builder = smsd.MolGraphBuilder().atom_count(2).atomic_numbers([6, 8]).neighbors([[1], [0]]).bond_orders([[1], [1]])
    >>> _ = builder.atomic_numbers([6, 8])
    >>> graph = builder.build(perceive_aromaticity=False)
    >>> assert graph.atomic_num == [6, 8]
)doc")
        .def("formal_charges", [](smsd::MolGraph::Builder& b, std::vector<int> v) -> smsd::MolGraph::Builder& {
            return b.formalCharges(std::move(v));
        }, py::arg("charges"), py::return_value_policy::reference_internal,
             "Set formal charges"
          R"doc(

Example::

    >>> import smsd
    >>> builder = smsd.MolGraphBuilder().atom_count(2).atomic_numbers([6, 8]).neighbors([[1], [0]]).bond_orders([[1], [1]])
    >>> _ = builder.formal_charges([0, 0])
    >>> graph = builder.build(perceive_aromaticity=False)
    >>> assert graph.formal_charge == [0, 0]
)doc")
        .def("mass_numbers", [](smsd::MolGraph::Builder& b, std::vector<int> v) -> smsd::MolGraph::Builder& {
            return b.massNumbers(std::move(v));
        }, py::arg("masses"), py::return_value_policy::reference_internal,
          R"doc(Set isotope mass numbers, with zero meaning unspecified.


Example::

    >>> import smsd
    >>> builder = smsd.MolGraphBuilder().atom_count(2).atomic_numbers([6, 8]).neighbors([[1], [0]]).bond_orders([[1], [1]])
    >>> _ = builder.mass_numbers([13, 0])
    >>> graph = builder.build(perceive_aromaticity=False)
    >>> assert graph.mass_number == [13, 0]
)doc")
        .def("hydrogen_counts", [](smsd::MolGraph::Builder& b, std::vector<int> v) -> smsd::MolGraph::Builder& {
            return b.hydrogenCounts(std::move(v));
        }, py::arg("counts"), py::return_value_policy::reference_internal,
          R"doc(Set per-atom hydrogen counts.


Example::

    >>> import smsd
    >>> builder = smsd.MolGraphBuilder().atom_count(2).atomic_numbers([6, 8]).neighbors([[1], [0]]).bond_orders([[1], [1]])
    >>> _ = builder.hydrogen_counts([3, 1])
    >>> graph = builder.build(perceive_aromaticity=False)
    >>> assert graph.hydrogen_count == [3, 1]
)doc")
        .def("ring_flags", [](smsd::MolGraph::Builder& b, std::vector<uint8_t> v) -> smsd::MolGraph::Builder& {
            return b.ringFlags(std::move(v));
        }, py::arg("flags"), py::return_value_policy::reference_internal,
             "Set ring membership flags"
          R"doc(

Example::

    >>> import smsd
    >>> builder = smsd.MolGraphBuilder().atom_count(2).atomic_numbers([6, 8]).neighbors([[1], [0]]).bond_orders([[1], [1]])
    >>> _ = builder.ring_flags([False, False])
    >>> graph = builder.build(perceive_aromaticity=False)
    >>> assert not any(graph.ring)
)doc")
        .def("aromatic_flags", [](smsd::MolGraph::Builder& b, std::vector<uint8_t> v) -> smsd::MolGraph::Builder& {
            return b.aromaticFlags(std::move(v));
        }, py::arg("flags"), py::return_value_policy::reference_internal,
             "Set aromaticity flags"
          R"doc(

Example::

    >>> import smsd
    >>> builder = smsd.MolGraphBuilder().atom_count(2).atomic_numbers([6, 8]).neighbors([[1], [0]]).bond_orders([[1], [1]])
    >>> _ = builder.aromatic_flags([False, False])
    >>> graph = builder.build(perceive_aromaticity=False)
    >>> assert not any(graph.aromatic)
)doc")
        .def("neighbors", [](smsd::MolGraph::Builder& b, std::vector<std::vector<int>> v) -> smsd::MolGraph::Builder& {
            return b.setNeighbors(std::move(v));
        }, py::arg("adj"), py::return_value_policy::reference_internal,
             "Set adjacency lists"
          R"doc(

Example::

    >>> import smsd
    >>> builder = smsd.MolGraphBuilder().atom_count(2).atomic_numbers([6, 8]).neighbors([[1], [0]]).bond_orders([[1], [1]])
    >>> _ = builder.neighbors([[1], [0]])
    >>> graph = builder.build(perceive_aromaticity=False)
    >>> assert graph.neighbors == [[1], [0]]
)doc")
        .def("bond_orders", [](smsd::MolGraph::Builder& b, std::vector<std::vector<int>> v) -> smsd::MolGraph::Builder& {
            return b.setBondOrders(std::move(v));
        }, py::arg("orders"), py::return_value_policy::reference_internal,
             "Set bond orders (parallel to neighbors)"
          R"doc(

Example::

    >>> import smsd
    >>> builder = smsd.MolGraphBuilder().atom_count(2).atomic_numbers([6, 8]).neighbors([[1], [0]]).bond_orders([[1], [1]])
    >>> _ = builder.bond_orders([[1], [1]])
    >>> graph = builder.build(perceive_aromaticity=False)
    >>> assert graph.bond_order(0, 1) == 1
)doc")
        .def("bond_ring_flags", [](smsd::MolGraph::Builder& b, std::vector<std::vector<bool>> v) -> smsd::MolGraph::Builder& {
            return b.bondRingFlags(std::move(v));
        }, py::arg("flags"), py::return_value_policy::reference_internal,
          R"doc(Set bond ring flags parallel to the adjacency lists.


Example::

    >>> import smsd
    >>> builder = smsd.MolGraphBuilder().atom_count(2).atomic_numbers([6, 8]).neighbors([[1], [0]]).bond_orders([[1], [1]])
    >>> _ = builder.bond_ring_flags([[False], [False]])
    >>> graph = builder.build(perceive_aromaticity=False)
    >>> assert not graph.bond_in_ring(0, 1)
)doc")
        .def("bond_aromatic_flags", [](smsd::MolGraph::Builder& b, std::vector<std::vector<bool>> v) -> smsd::MolGraph::Builder& {
            return b.bondAromaticFlags(std::move(v));
        }, py::arg("flags"), py::return_value_policy::reference_internal,
          R"doc(Set bond aromaticity flags parallel to the adjacency lists.


Example::

    >>> import smsd
    >>> builder = smsd.MolGraphBuilder().atom_count(2).atomic_numbers([6, 8]).neighbors([[1], [0]]).bond_orders([[1], [1]])
    >>> _ = builder.bond_aromatic_flags([[False], [False]])
    >>> graph = builder.build(perceive_aromaticity=False)
    >>> assert not graph.bond_aromatic(0, 1)
)doc")
        .def("atom_ids", [](smsd::MolGraph::Builder& b, std::vector<int> v) -> smsd::MolGraph::Builder& {
            return b.atomIds(std::move(v));
        }, py::arg("ids"), py::return_value_policy::reference_internal,
             "Set optional external atom IDs (1-based, gaps allowed) for index translation"
          R"doc(

Example::

    >>> import smsd
    >>> builder = smsd.MolGraphBuilder().atom_count(2).atomic_numbers([6, 8]).neighbors([[1], [0]]).bond_orders([[1], [1]])
    >>> _ = builder.atom_ids([101, 205])
    >>> graph = builder.build(perceive_aromaticity=False)
    >>> assert smsd.translate_to_atom_ids({0: 1}, graph, graph) == {101: 205}
)doc")
        .def("name", [](smsd::MolGraph::Builder& b, std::string v) -> smsd::MolGraph::Builder& {
            return b.name(std::move(v));
        }, py::arg("value"), py::return_value_policy::reference_internal,
             "Set molecule name/title"
          R"doc(

Example::

    >>> import smsd
    >>> builder = smsd.MolGraphBuilder().atom_count(2).atomic_numbers([6, 8]).neighbors([[1], [0]]).bond_orders([[1], [1]])
    >>> _ = builder.name("Methanol")
    >>> graph = builder.build(perceive_aromaticity=False)
    >>> assert graph.name == "Methanol"
)doc")
        .def("program_line", [](smsd::MolGraph::Builder& b, std::string v) -> smsd::MolGraph::Builder& {
            return b.programLine(std::move(v));
        }, py::arg("value"), py::return_value_policy::reference_internal,
             "Set MDL program line/header text"
          R"doc(

Example::

    >>> import smsd
    >>> builder = smsd.MolGraphBuilder().atom_count(2).atomic_numbers([6, 8]).neighbors([[1], [0]]).bond_orders([[1], [1]])
    >>> _ = builder.program_line("Molecule example")
    >>> graph = builder.build(perceive_aromaticity=False)
    >>> assert graph.program_line == "Molecule example"
)doc")
        .def("comment", [](smsd::MolGraph::Builder& b, std::string v) -> smsd::MolGraph::Builder& {
            return b.comment(std::move(v));
        }, py::arg("value"), py::return_value_policy::reference_internal,
             "Set molecule comment text"
          R"doc(

Example::

    >>> import smsd
    >>> builder = smsd.MolGraphBuilder().atom_count(2).atomic_numbers([6, 8]).neighbors([[1], [0]]).bond_orders([[1], [1]])
    >>> _ = builder.comment("Two heavy atoms")
    >>> graph = builder.build(perceive_aromaticity=False)
    >>> assert graph.comment == "Two heavy atoms"
)doc")
        .def("properties", [](smsd::MolGraph::Builder& b, std::map<std::string, std::string> v) -> smsd::MolGraph::Builder& {
            return b.properties(std::move(v));
        }, py::arg("values"), py::return_value_policy::reference_internal,
             "Set SDF-style molecule properties"
          R"doc(

Example::

    >>> import smsd
    >>> builder = smsd.MolGraphBuilder().atom_count(2).atomic_numbers([6, 8]).neighbors([[1], [0]]).bond_orders([[1], [1]])
    >>> _ = builder.properties({"sample": "methanol"})
    >>> graph = builder.build(perceive_aromaticity=False)
    >>> assert graph.properties["sample"] == "methanol"
)doc")
        .def("build", [](const smsd::MolGraph::Builder& b,
                         bool perceive_aromaticity,
                         smsd::AromaticityModel aromaticity_model) {
            return b.build(perceive_aromaticity, aromaticity_model);
        },
             py::arg("perceive_aromaticity") = true,
             py::arg("aromaticity_model") = smsd::AromaticityModel::DAYLIGHT_LIKE,
             "Build the immutable MolGraph, optionally running aromaticity perception"
          R"doc(

Example::

    >>> import smsd
    >>> builder = smsd.MolGraphBuilder().atom_count(2).atomic_numbers([6, 8]).neighbors([[1], [0]]).bond_orders([[1], [1]])
    >>> graph = builder.hydrogen_counts([3, 1]).build()
    >>> assert graph.n == 2 and graph.bond_order(0, 1) == 1
)doc");

    py::class_<smsd::ParseOptions>(m, "ParseOptions")
        .def(py::init<>())
        .def_readwrite("lenient", &smsd::ParseOptions::lenient,
                       "If true, recover from malformed SMILES instead of throwing");

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
          "Parse an MDL MOL V2000 block into a MolGraph"
          R"doc(

Example::

    >>> import smsd
    >>> restored = smsd.read_mol_block(smsd.write_mol_block("CCO"))
    >>> assert restored.n == 3
)doc");
    m.def("write_mol_block", &smsd::writeMolBlock,
          py::arg("mol"),
          "Write a MolGraph as an MDL MOL V2000 block"
          R"doc(

Example::

    >>> import smsd
    >>> text = smsd.write_mol_block("CCO")
    >>> assert "V2000" in text and smsd.read_mol_block(text).n == 3
)doc");
    m.def("write_mol_block_v3000", &smsd::writeMolBlockV3000,
          py::arg("mol"),
          "Write a MolGraph as an MDL MOL V3000 block"
          R"doc(

Example::

    >>> import smsd
    >>> text = smsd.write_mol_block_v3000("CCO")
    >>> assert "V3000" in text and smsd.read_mol_block(text).n == 3
)doc");
    m.def("write_sdf_record", &smsd::writeSDFRecord,
          py::arg("mol"),
          "Write a MolGraph as a single SDF record"
          R"doc(

Example::

    >>> import smsd
    >>> text = smsd.write_sdf_record("CCO")
    >>> assert text.rstrip().endswith("$$$$") and "V2000" in text
)doc");

    auto searchStatus = py::enum_<smsd::SearchStatus>(m, "SearchStatus",
          R"doc(Completion or stopping reason from the native search.


Example::

    >>> import smsd
    >>> result = smsd.MCS.find("CCO").with_timeout(0).match("CCN")
    >>> assert result.status == smsd.SearchStatus.TIMED_OUT and not result.optimal
)doc")
        .value("COMPLETED", smsd::SearchStatus::Completed)
        .value("TIMED_OUT", smsd::SearchStatus::TimedOut)
        .value("CANCELLED", smsd::SearchStatus::Cancelled)
        .value("LIMIT_REACHED", smsd::SearchStatus::LimitReached);
    initialiseResultTypes(searchStatus);

    py::class_<smsd::CancellationToken, std::shared_ptr<smsd::CancellationToken>>(m, "CancellationToken",
          R"doc(A cancellation signal shared with search options.


Example::

    >>> import smsd
    >>> token = smsd.CancellationToken()
    >>> assert not token.is_cancelled
    >>> token.cancel()
    >>> assert token.is_cancelled
)doc")
        .def(py::init<>())
        .def("cancel", &smsd::CancellationToken::cancel,
          R"doc(Request cancellation; the token cannot be reset.


Example::

    >>> import smsd
    >>> token = smsd.CancellationToken()
    >>> token.cancel()
    >>> assert token.is_cancelled
)doc")
        .def_property_readonly("is_cancelled", &smsd::CancellationToken::isCancelled,
          R"doc(Whether cancellation has been requested.


Example::

    >>> import smsd
    >>> token = smsd.CancellationToken()
    >>> assert not token.is_cancelled
    >>> token.cancel()
    >>> assert token.is_cancelled
)doc");

    py::class_<smsd::SearchOptions>(m, "_SearchOptions")
        .def(py::init<>())
        .def_readwrite("chemistry", &smsd::SearchOptions::chemistry)
        .def_readwrite("timeout_ms", &smsd::SearchOptions::timeoutMs)
        .def_readwrite("node_limit", &smsd::SearchOptions::nodeLimit)
        .def_readwrite("cancellation", &smsd::SearchOptions::cancellation);

    py::class_<smsd::Match>(m, "Match")
        .def(py::init<>())
        .def_readwrite("atom_mapping", &smsd::Match::atomMapping)
        .def_readwrite("bond_mapping", &smsd::Match::bondMapping)
        .def_readwrite("fragments", &smsd::Match::fragments)
        .def("to_graph", &smsd::Match::toGraph, py::arg("query"))
        .def("to_smiles", &smsd::Match::toSmiles, py::arg("query"))
        .def("to_smarts", &smsd::Match::toSmarts, py::arg("query"), py::arg("chemistry") = smsd::ChemOptions());

    py::class_<smsd::SearchResult>(m, "SearchResult")
        .def_readonly("matches", &smsd::SearchResult::matches)
        .def_readonly("status", &smsd::SearchResult::status)
        .def_readonly("optimal", &smsd::SearchResult::optimal)
        .def_readonly("complete", &smsd::SearchResult::complete)
        .def_readonly("score", &smsd::SearchResult::score)
        .def_readonly("nodes_visited", &smsd::SearchResult::nodesVisited)
        .def_readonly("elapsed_ms", &smsd::SearchResult::elapsedMillis);

    m.def("_register_search_types", &registerResultTypes, py::arg("match"), py::arg("result"));
    // Converting stops once `cancellation` is cancelled or, unless `time_limit_ms` is None,
    // that many milliseconds from now have passed, as for the results of a search.
    m.def("_search_result",
          [](const smsd::SearchResult& result, py::handle queryIndices, py::handle targetIndices,
             std::shared_ptr<smsd::CancellationToken> cancellation, py::object timeLimitMillis) {
              const auto deadline = timeLimitMillis.is_none()
                  ? ResultsDeadline::cancelledBy(std::move(cancellation))
                  : ResultsDeadline::after(timeLimitMillis.cast<int64_t>(), std::chrono::steady_clock::now(),
                                           std::move(cancellation));
              return pythonResult(result, queryIndices, targetIndices, deadline);
          },
          py::arg("result"), py::arg("query_indices"), py::arg("target_indices"),
          py::arg("cancellation") = nullptr, py::arg("time_limit_ms") = py::none());
    m.def("_python_match",
          [](const smsd::Match& match, py::handle queryIndices, py::handle targetIndices) {
              registeredResultTypes();
              return pythonMatch(match, AtomNumbering(queryIndices), AtomNumbering(targetIndices));
          });

    auto mcs = py::class_<smsd::MCS>(m, "MCS");
    py::enum_<smsd::MCS::Objective>(mcs, "Objective",
          R"doc(Optimise common atoms (MCIS) or selected common bonds (MCES).


Example::

    >>> import smsd
    >>> result = smsd.MCS.find("CC").with_objective(smsd.MCS.Objective.MCES).match("CCC")
    >>> assert result.score == 1 and result.best.bond_count == 1
)doc")
        .value("MCIS", smsd::MCS::Objective::MCIS)
        .value("MCES", smsd::MCS::Objective::MCES);
    py::enum_<smsd::MCS::Connectivity>(mcs, "Connectivity",
          R"doc(Require a connected common graph or allow fragments.


Example::

    >>> import smsd
    >>> options = smsd.MCS.Options()
    >>> options.connectivity = smsd.MCS.Connectivity.FRAGMENTS
    >>> assert smsd.MCS.find("C.O").with_options(options).match("C.O").best.atom_count == 2
)doc")
        .value("CONNECTED", smsd::MCS::Connectivity::Connected)
        .value("FRAGMENTS", smsd::MCS::Connectivity::Fragments);
    py::class_<smsd::MCS::Options, smsd::SearchOptions>(mcs, "Options",
          R"doc(MCS objective, connectivity, chemistry and stopping limits.


Example::

    >>> import smsd
    >>> options = smsd.MCS.Options()
    >>> options.timeout_ms, options.node_limit = 1000, 100000
    >>> options.objective = smsd.MCS.Objective.MCIS
    >>> assert smsd.MCS.find("CCO").with_options(options).match("CCN").found
)doc")
        .def(py::init<>())
        .def(py::init<const smsd::MCS::Options&>())
        .def_readwrite("objective", &smsd::MCS::Options::objective)
        .def_readwrite("connectivity", &smsd::MCS::Options::connectivity)
        .def_readwrite("atom_weights", &smsd::MCS::Options::atomWeights)
        .def_readwrite("min_fragment_atoms", &smsd::MCS::Options::minFragmentAtoms)
        .def_readwrite("max_fragments", &smsd::MCS::Options::maxFragments)
        .def_readwrite("excluded_target_atoms", &smsd::MCS::Options::excludedTargetAtoms);
    mcs.def_static("find", py::overload_cast<const smsd::MolGraph&>(&smsd::MCS::find), py::arg("query"))
        .def("with_options", &smsd::MCS::withOptions, py::arg("options"))
        .def("with_timeout", &smsd::MCS::withTimeout, py::arg("milliseconds"))
        .def("match", py::overload_cast<const smsd::MolGraph&>(&smsd::MCS::match, py::const_),
             py::arg("target"), py::call_guard<py::gil_scoped_release>())
        .def("match_all", [](const smsd::MCS& query, const smsd::MolGraph& target, int limit) {
            if (limit < 1) throw std::invalid_argument("Match limit must be positive");
            py::gil_scoped_release release;
            return query.matchAll(target, static_cast<size_t>(limit));
        }, py::arg("target"), py::arg("limit") = 1000)
        .def("match_many", &matchMany<smsd::MCS>, py::arg("targets"), py::arg("num_threads") = 0);
    mcs.def("_match_corpus_prefiltered", [](const smsd::MCS& query, const smsd::batch::TargetCorpus& corpus,
        const py::sequence& targets, int minimumAtoms, bool gpuPrefilter, int numThreads, py::handle queryIndices) {
        return matchCorpusPython(query, corpus, targets, minimumAtoms, gpuPrefilter, numThreads, queryIndices);
    }, py::arg("corpus"), py::arg("targets"), py::arg("minimum_atoms"), py::arg("gpu_prefilter"),
       py::arg("num_threads"), py::arg("query_indices"));
    bindPythonSearches(mcs);

    py::class_<smsd::MCS::ConsensusResult>(mcs, "ConsensusResult")
        .def_readonly("fragment", &smsd::MCS::ConsensusResult::fragment)
        .def_readonly("source_index", &smsd::MCS::ConsensusResult::sourceIndex)
        .def_readonly("graph", &smsd::MCS::ConsensusResult::graph)
        .def_readonly("mappings", &smsd::MCS::ConsensusResult::mappings)
        .def_readonly("support_count", &smsd::MCS::ConsensusResult::supportCount)
        .def_readonly("found", &smsd::MCS::ConsensusResult::found)
        .def_readonly("status", &smsd::MCS::ConsensusResult::status)
        .def_readonly("optimal", &smsd::MCS::ConsensusResult::optimal)
        .def_readonly("complete", &smsd::MCS::ConsensusResult::complete)
        .def_readonly("nodes_visited", &smsd::MCS::ConsensusResult::nodesVisited)
        .def_readonly("elapsed_ms", &smsd::MCS::ConsensusResult::elapsedMillis)
        .def("to_smiles", &smsd::MCS::ConsensusResult::toSmiles);
    // The options are copied while the GIL is held, so Python threads that change them, such
    // as clearing their cancellation token, cannot race with the search.
    mcs.def_static("consensus",
        [](const std::vector<smsd::MolGraph>& molecules, const smsd::MCS::Options& options, double threshold) {
            const smsd::MCS::Options copied = options;
            py::gil_scoped_release release;
            return smsd::MCS::consensus(molecules, copied, threshold);
        },
        py::arg("molecules"), py::arg("options") = smsd::MCS::Options(), py::arg("threshold") = 1.0);
    mcs.def_static("match_constrained",
        [](const std::vector<smsd::MolGraph>& queries, const smsd::MolGraph& target, const smsd::MCS::Options& options) {
            const smsd::MCS::Options copied = options;
            py::gil_scoped_release release;
            return smsd::MCS::matchConstrained(queries, target, copied);
        },
        py::arg("queries"), py::arg("target"), py::arg("options") = smsd::MCS::Options());
    mcs.def_static("compute_tautomer_confidence", &smsd::MCS::computeTautomerConfidence,
        py::arg("query"), py::arg("target"), py::arg("match"), py::arg("chemistry"),
        py::call_guard<py::gil_scoped_release>());

    auto substructure = py::class_<smsd::Substructure>(m, "Substructure");
    py::class_<smsd::Substructure::Options, smsd::SearchOptions>(substructure, "Options",
          R"doc(Molecule and SMARTS options, including induced topology.


Example::

    >>> import smsd
    >>> options = smsd.Substructure.Options()
    >>> options.induced = True
    >>> assert not smsd.Substructure.find("CCC").with_options(options).match("C1CC1").found
)doc")
        .def(py::init<>())
        .def(py::init<const smsd::Substructure::Options&>())
        .def_readwrite("induced", &smsd::Substructure::Options::induced);
    substructure.def_static("find", py::overload_cast<const smsd::MolGraph&>(&smsd::Substructure::find), py::arg("query"))
        .def_static("smarts", py::overload_cast<const std::string&>(&smsd::Substructure::smarts), py::arg("query"))
        .def("with_options", &smsd::Substructure::withOptions, py::arg("options"))
        .def("with_timeout", &smsd::Substructure::withTimeout, py::arg("milliseconds"))
        .def("match", py::overload_cast<const smsd::MolGraph&>(&smsd::Substructure::match, py::const_),
             py::arg("target"), py::call_guard<py::gil_scoped_release>())
        .def("match_all", [](const smsd::Substructure& query, const smsd::MolGraph& target, int limit) {
            if (limit < 1) throw std::invalid_argument("Match limit must be positive");
            py::gil_scoped_release release;
            return query.matchAll(target, static_cast<size_t>(limit));
        }, py::arg("target"), py::arg("limit") = 1000)
        .def("match_many", &matchMany<smsd::Substructure>, py::arg("targets"), py::arg("num_threads") = 0)
        .def("_match_corpus", [](const smsd::Substructure& query, const smsd::batch::TargetCorpus& corpus,
                                const py::sequence& targets, int numThreads, py::handle queryIndices) {
            return matchManyPython(query, AtomNumbering(queryIndices), targets, numThreads, &corpus);
        }, py::arg("corpus"), py::arg("targets"), py::arg("num_threads"), py::arg("query_indices"));
    substructure.def("_match_corpus_prefiltered", [](const smsd::Substructure& query, const smsd::batch::TargetCorpus& corpus,
        const py::sequence& targets, bool gpuPrefilter, int numThreads, py::handle queryIndices) {
        return matchCorpusPython(query, corpus, targets, 1, gpuPrefilter, numThreads, queryIndices);
    }, py::arg("corpus"), py::arg("targets"), py::arg("gpu_prefilter"), py::arg("num_threads"), py::arg("query_indices"));
    bindPythonSearches(substructure);

    mcs.def_static("canonicalize_mapping",
          [](const smsd::MolGraph& g1, const smsd::MolGraph& g2,
             const std::map<int,int>& mapping) {
              return smsd::MCS::canonicalizeMapping(g1, g2, mapping);
          },
          py::arg("g1"), py::arg("g2"), py::arg("mapping"),
          py::call_guard<py::gil_scoped_release>(),
          "Canonicalise an atom-atom mapping under the automorphism groups of\n"
          "both molecules. Returns the lexicographically smallest equivalent\n"
          "mapping. Two mappings that differ only by a symmetry of g1 or g2\n"
          "give the same canonical result.");

    m.def("translate_to_atom_ids",
          [](const std::map<int,int>& mapping,
             const smsd::MolGraph& g1, const smsd::MolGraph& g2) {
              return smsd::translateToAtomIds(mapping, g1, g2);
          },
          py::arg("mapping"), py::arg("g1"), py::arg("g2"),
          "Translate MCS mapping from internal 0-based indices to external atom IDs.\n"
          "If neither graph has atomId set, the mapping is returned unchanged."
          R"doc(

Example::

    >>> import smsd
    >>> graph = smsd.parse_smiles("CCO")
    >>> assert smsd.translate_to_atom_ids({0: 1}, graph, graph) == {0: 1}
)doc");

    m.def("validate_tautomer_consistency", &smsd::validateTautomerConsistency,
          py::arg("g1"), py::arg("g2"), py::arg("mcs"),
          "Validate that an MCS mapping conserves mobilisable protons across "
          "tautomeric groups. Returns True if consistent, False if proton "
          "counts mismatch (impossible proton transfer required).");

    m.def("similarity_upper_bound",
          [](const smsd::MolGraph& g1, const smsd::MolGraph& g2, const smsd::ChemOptions& chemistry) {
              return smsd::similarityUpperBound(g1, g2, chemistry);
          },
          py::arg("g1"), py::arg("g2"), py::arg("chem") = smsd::ChemOptions(),
          "Linear-time similarity upper bound");

    m.def("screen_targets", &screenPythonTargets,
          py::arg("query"), py::arg("targets"), py::arg("threshold"), py::arg("chem") = smsd::ChemOptions(),
          py::arg("num_threads") = 0, py::arg("convert") = py::none(),
          "Batch similarity-bound screening: the indices of the targets whose bound reaches the "
          "threshold, in order. Targets are MolGraph objects, SMILES strings parsed natively, or "
          "items that convert(item) turns into a MolGraph; num_threads=0 uses all available processors.");

    m.def("path_fingerprint",
          [](const smsd::MolGraph& mol, int pathLength, int fpSize) {
              return smsd::bitPositions(
                  smsd::batch::detail::computeElementPathFingerprint(mol, pathLength, fpSize), fpSize);
          },
          py::arg("mol"), py::arg("path_length") = 7, py::arg("fp_size") = 2048,
          "Element path fingerprint over heavy atoms: each path of up to path_length bonds is "
          "hashed by its atomic numbers. Returns list of set bit positions. Distinct from "
          "batch_fingerprint, which also hashes atom labels, bond orders and hydrogens."
          R"doc(

Example::

    >>> import smsd
    >>> bits = smsd.path_fingerprint(smsd.parse_smiles("CCO"), fp_size=128)
    >>> assert bits and all(0 <= bit < 128 for bit in bits)
)doc");

    m.def("mcs_fingerprint",
          [](const smsd::MolGraph& mol, int pathLength, int fpSize, bool tautomerAware) {
              smsd::ChemOptions chemistry;
              chemistry.tautomerAware = tautomerAware;
              return smsd::bitPositions(
                  smsd::batch::detail::computeMCSFingerprintForChemistry(mol, pathLength, fpSize, chemistry),
                  fpSize);
          },
          py::arg("mol"), py::arg("path_length") = 7, py::arg("fp_size") = 2048,
          py::arg("tautomer_aware") = false,
          "Compute MCS-aware path fingerprint. Returns list of set bit positions."
          R"doc(

Example::

    >>> import smsd
    >>> bits = smsd.mcs_fingerprint(smsd.parse_smiles("CCO"), fp_size=128)
    >>> assert bits and all(0 <= bit < 128 for bit in bits)
)doc");

    m.def("fingerprint_subset", &fingerprintSubset,
          py::arg("fp1"), py::arg("fp2"),
          "Check if fingerprint fp1 is a subset of fp2 (all bits in fp1 are set in fp2)"
          R"doc(

Example::

    >>> import smsd
    >>> assert smsd.fingerprint_subset([1, 3], [1, 3, 8])
    >>> assert not smsd.fingerprint_subset([1, 5], [1, 3, 8])
)doc");

    m.def("analyze_fp_quality", &analyzeFpQuality,
          py::arg("fp"), py::arg("fp_size") = 2048,
          "Analyse fingerprint quality: density, set bits, collision rate"
          R"doc(

Example::

    >>> import smsd
    >>> quality = smsd.analyze_fp_quality([1, 3], fp_size=64)
    >>> assert quality["set_bits"] == 2 and quality["density"] == 2 / 64
)doc");

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
          "Generate canonical SMARTS with configurable predicates"
          R"doc(

Example::

    >>> import smsd
    >>> pattern = smsd.to_smarts("CCO")
    >>> assert smsd.Substructure.smarts(pattern).match("CCO").found
)doc");

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

    // The screening accelerator of a separately installed package, such as smsd-cuda, which the
    // smsd package attaches when it finds one; TargetCorpus screens then compute their label
    // overlaps there and give the same hits.
    m.def("_attach_screen_accelerator",
          [](const py::capsule& capsule) {
              const char* name = capsule.name();
              if (name == nullptr || std::strcmp(name, SMSD_SCREEN_ACCELERATOR_CAPSULE) != 0)
                  throw py::value_error("the capsule does not hold a screening accelerator");
              const std::string reason = smsd::screen_accelerator::attach(
                  capsule.get_pointer<const smsd_screen_accelerator_v1>(), SMSD_VERSION_INFO);
              if (!reason.empty()) throw py::value_error(reason);
          },
          py::arg("capsule"),
          "Attach the screening accelerator the capsule holds; its owner must stay loaded.");
    m.def("_detach_screen_accelerator", &smsd::screen_accelerator::detach,
          "Detach the screening accelerator, so that screens run on the CPU.");
    m.def("_screen_accelerator",
          []() -> py::object {
              const smsd_screen_accelerator_v1* table = smsd::screen_accelerator::active();
              if (table == nullptr) return py::none();
              return py::make_tuple(std::string(table->name), smsd::screen_accelerator::deviceInfo());
          },
          "The active screening accelerator as (name, device), or None.");
    m.def("_screen_accelerator_failures",
          []() {
              return py::make_tuple(smsd::screen_accelerator::failureCount(),
                                    smsd::screen_accelerator::latestFailure());
          },
          "The number of accelerated screens that ran on the CPU instead, and the latest reason.");

    m.def("batch_fingerprint",
          [](const std::vector<smsd::MolGraph>& mols,
             int pathLength, int fpSize, int numThreads) {
              return smsd::batch::batchFingerprint(mols, pathLength, fpSize, numThreads);
          },
          py::arg("mols"), py::arg("path_length") = 7,
          py::arg("fp_size") = 2048, py::arg("num_threads") = 0,
          py::call_guard<py::gil_scoped_release>(),
          "Parallel path fingerprint generation for N molecules over all atoms, hashing atom "
          "labels and bond orders; returns fingerprint words. Distinct from path_fingerprint.");

    py::class_<smsd::batch::TargetCorpus>(m, "TargetCorpus",
        "A fingerprinted target collection for repeated queries. Call prewarm() once after "
        "loading targets, then run substructure, MCS and screening queries many times without "
        "recomputing graph invariants.")
        .def(py::init<>())
        .def("add_target", &smsd::batch::TargetCorpus::addTarget, py::arg("mol"),
             py::call_guard<py::gil_scoped_release>(),
             "Add a single target molecule.")
        .def("add_targets", &addCorpusTargets,
             py::arg("mols"),
             "Add multiple target molecules.")
        .def("prewarm",
             [](smsd::batch::TargetCorpus& self, int numThreads) {
                 self.prewarm(numThreads);
             },
             py::arg("num_threads") = 0,
             py::call_guard<py::gil_scoped_release>(),
             "Pre-compute graph invariants and ECFP fingerprints for all targets. "
             "num_threads=0 uses all available processors.")
        .def("__len__", &smsd::batch::TargetCorpus::size,
             py::call_guard<py::gil_scoped_release>(),
             "Number of targets in the corpus.")
        .def("targets", &smsd::batch::TargetCorpus::targets, py::call_guard<py::gil_scoped_release>())
        .def("_search_targets", &corpusSearchTargets)
        .def_property_readonly("is_prewarmed", [](const smsd::batch::TargetCorpus& self) {
                 py::gil_scoped_release release;
                 return self.isPrewarmed();
             },
             "True if prewarm() has been called since last target addition.")
        .def("screen",
             [](const smsd::batch::TargetCorpus& self,
                const smsd::MolGraph& query, double threshold, int nThreads) {
                 return self.screen(query, threshold, nThreads);
             },
             py::arg("query"),
             py::arg("threshold"),
             py::arg("num_threads") = 0,
             py::call_guard<py::gil_scoped_release>(),
             "Tanimoto fingerprint screening. Returns indices of targets above threshold. "
             "Requires prewarm() to have been called.")
        .def("screen_many",
             [](const smsd::batch::TargetCorpus& self, const py::iterable& queries,
                double threshold, int numThreads) {
                 if (numThreads < 0) throw std::invalid_argument("Thread count must be non-negative");
                 std::vector<smsd::MolGraph> graphs;
                 for (py::handle value : queries) {
                     const auto item = py::reinterpret_borrow<py::object>(value);
                     const auto& graph = item.cast<const smsd::MolGraph&>();
                     const std::lock_guard<std::recursive_mutex> graphLock(graph.lazyCacheMutex());
                     graphs.push_back(graph);
                 }
                 py::gil_scoped_release release;
                 return self.screenMany(graphs, threshold, numThreads);
             },
             py::arg("queries"), py::arg("threshold"), py::arg("num_threads") = 0,
             "ECFP Tanimoto screens in query order. Requires prewarm().")
        .def("screen_targets",
             [](const smsd::batch::TargetCorpus& self, const smsd::MolGraph& query, double threshold,
                const smsd::ChemOptions& given, int numThreads) {
                 if (numThreads < 0) throw std::invalid_argument("Thread count must be non-negative");
                 // The options are copied while the GIL is held, as for the searches.
                 const smsd::ChemOptions chemistry = given;
                 py::gil_scoped_release release;
                 return self.screenTargets(query, threshold, chemistry, numThreads);
             },
             py::arg("query"), py::arg("threshold"), py::arg("chem") = smsd::ChemOptions(),
             py::arg("num_threads") = 0,
             "Similarity-bound screening of the targets, as screen_targets gives for them; needs no "
             "prewarm(). num_threads=0 uses all available processors.")
        .def("screen_targets_many",
             [](const smsd::batch::TargetCorpus& self, const py::iterable& queries,
                double threshold, const smsd::ChemOptions& given, int numThreads) {
                 if (numThreads < 0) throw std::invalid_argument("Thread count must be non-negative");
                 std::vector<smsd::MolGraph> graphs;
                 for (py::handle value : queries) {
                     const auto item = py::reinterpret_borrow<py::object>(value);
                     const auto& graph = item.cast<const smsd::MolGraph&>();
                     const std::lock_guard<std::recursive_mutex> graphLock(graph.lazyCacheMutex());
                     graphs.push_back(graph);
                 }
                 const smsd::ChemOptions chemistry = given;
                 py::gil_scoped_release release;
                 return self.screenTargetsMany(graphs, threshold, chemistry, numThreads);
             },
             py::arg("queries"), py::arg("threshold"), py::arg("chem") = smsd::ChemOptions(),
             py::arg("num_threads") = 0,
             "Similarity-bound screens in query order; needs no prewarm().");

    m.def("circular_fingerprint",
          [](const smsd::MolGraph& mol, int radius, int fpSize,
             const std::string& mode) {
              auto fp = circularMode(mode)
                  ? smsd::batch::detail::computeCircularFingerprintFCFP(mol, radius, fpSize)
                  : smsd::batch::detail::computeCircularFingerprintECFP(mol, radius, fpSize);
              return smsd::bitPositions(fp, fpSize);
          },
          py::arg("mol"), py::arg("radius") = 2,
          py::arg("fp_size") = 2048, py::arg("mode") = "ecfp",
          "Circular fingerprint (Morgan). mode='ecfp' (structural) or 'fcfp' "
          "(pharmacophoric); any other mode raises ValueError. Returns list of set bit positions."
          R"doc(

Example::

    >>> import smsd
    >>> bits = smsd.circular_fingerprint(smsd.parse_smiles("CCO"), fp_size=128)
    >>> assert bits and all(0 <= bit < 128 for bit in bits)
)doc");

    m.def("circular_fingerprint_counts",
          [](const smsd::MolGraph& mol, int radius, int fpSize,
             const std::string& mode) {
              return nonZeroBins(circularMode(mode)
                  ? smsd::batch::detail::computeCircularFingerprintFCFPCounts(mol, radius, fpSize)
                  : smsd::batch::detail::computeCircularFingerprintECFPCounts(mol, radius, fpSize));
          },
          py::arg("mol"), py::arg("radius") = 2,
          py::arg("fp_size") = 2048, py::arg("mode") = "ecfp",
          "Count-based circular fingerprint (Morgan). mode='ecfp' (structural) "
          "or 'fcfp' (pharmacophoric); any other mode raises ValueError. Returns list of "
          "(bit_position, count) tuples for non-zero bins, preserving the multiplicity "
          "that binary fingerprints lose."
          R"doc(

Example::

    >>> import smsd
    >>> counts = smsd.circular_fingerprint_counts(smsd.parse_smiles("CCCCCC"), fp_size=128)
    >>> assert counts and any(count > 1 for position, count in counts)
)doc");

    m.def("topological_torsion",
          [](const smsd::MolGraph& mol, int fpSize) {
              auto fp = smsd::batch::detail::computeTopologicalTorsion(mol, fpSize);
              return smsd::bitPositions(fp, fpSize);
          },
          py::arg("mol"), py::arg("fp_size") = 2048,
          "Topological torsion fingerprint (Nilakantan 1987). Enumerates all "
          "4-atom linear paths and hashes atom types (atomicNum, heavyDeg, "
          "piElectrons, ring). Returns list of set bit positions.");

    m.def("topological_torsion_counts",
          [](const smsd::MolGraph& mol, int fpSize) {
              return nonZeroBins(smsd::batch::detail::computeTopologicalTorsionCounts(mol, fpSize));
          },
          py::arg("mol"), py::arg("fp_size") = 2048,
          "Count-based topological torsion fingerprint. Returns list of "
          "(bit_position, count) tuples for non-zero bins."
          R"doc(

Example::

    >>> import smsd
    >>> counts = smsd.topological_torsion_counts(smsd.parse_smiles("CCCCC"), fp_size=128)
    >>> assert counts and all(0 <= position < 128 and count > 0 for position, count in counts)
)doc");

    // -----------------------------------------------------------------------
    // Configurable fingerprints
    // -----------------------------------------------------------------------
    using smsd::batch::FingerprintFamily;
    using smsd::batch::FingerprintKey;
    using smsd::batch::FingerprintOptions;
    using smsd::batch::FragmentMode;
    using smsd::batch::ScaffoldMode;

    py::enum_<FingerprintFamily>(m, "FingerprintFamily",
        "A fingerprint family: the circular ECFP and FCFP, PATH and TORSION families, or substructure KEYS."
          R"doc(

Example::

    >>> import smsd
    >>> assert smsd.FingerprintOptions(family="ecfp").family == smsd.FingerprintFamily.ECFP
)doc")
        .value("ECFP", FingerprintFamily::ECFP)
        .value("FCFP", FingerprintFamily::FCFP)
        .value("PATH", FingerprintFamily::PATH)
        .value("TORSION", FingerprintFamily::TORSION)
        .value("KEYS", FingerprintFamily::KEYS);

    py::enum_<FragmentMode>(m, "FragmentMode",
        "How the connected components of a molecule are fingerprinted: the WHOLE molecule, the "
        "LARGEST component only, or PER_FRAGMENT, where every atom term carries its component's class."
          R"doc(

Example::

    >>> import smsd
    >>> options = smsd.FingerprintOptions(fragments="largest")
    >>> assert options.fragments == smsd.FragmentMode.LARGEST
)doc")
        .value("WHOLE", FragmentMode::WHOLE)
        .value("LARGEST", FragmentMode::LARGEST)
        .value("PER_FRAGMENT", FragmentMode::PER_FRAGMENT);

    py::enum_<ScaffoldMode>(m, "ScaffoldMode",
        "Whether the molecule is reduced before it is fingerprinted: NONE, or its MURCKO framework."
          R"doc(

Example::

    >>> import smsd
    >>> options = smsd.FingerprintOptions(scaffold="murcko")
    >>> assert options.scaffold == smsd.ScaffoldMode.MURCKO
)doc")
        .value("NONE", ScaffoldMode::NONE)
        .value("MURCKO", ScaffoldMode::MURCKO);

    py::class_<FingerprintKey>(m, "FingerprintKey",
        "A substructure key: its bit is set when a molecule holds at least `minimum` distinct "
        "matches of its SMARTS pattern. The MACCS set also counts isotopes, aromatic rings and "
        "fragments, which `kind` names."
          R"doc(

Example::

    >>> import smsd
    >>> key = smsd.FingerprintKey("[#8]", minimum=1)
    >>> options = smsd.FingerprintOptions(family="keys", keys=[key])
    >>> assert smsd.fingerprint("CCO", options=options) == [0]
    >>> assert smsd.fingerprint("CCN", options=options) == []
)doc")
        .def(py::init<std::string, int>(), py::arg("smarts"), py::arg("minimum") = 1)
        .def_readonly("smarts", &FingerprintKey::smarts)
        .def_readonly("minimum", &FingerprintKey::minimum)
        .def_property_readonly("kind", [](const FingerprintKey& key) {
            switch (key.kind) {
            case smsd::batch::KeyKind::SMARTS: return "smarts";
            case smsd::batch::KeyKind::ISOTOPE: return "isotope";
            case smsd::batch::KeyKind::AROMATIC_RINGS: return "aromatic_rings";
            case smsd::batch::KeyKind::FRAGMENTS: return "fragments";
            }
            return "smarts";
        })
        .def("__repr__", [](const FingerprintKey& key) {
            if (key.kind != smsd::batch::KeyKind::SMARTS) {
                const char* kind = key.kind == smsd::batch::KeyKind::ISOTOPE ? "isotope"
                                 : key.kind == smsd::batch::KeyKind::AROMATIC_RINGS ? "aromatic_rings" : "fragments";
                return "FingerprintKey(<" + std::string(kind) + ">, minimum=" + std::to_string(key.minimum) + ")";
            }
            return "FingerprintKey(" + py::repr(py::str(key.smarts)).cast<std::string>()
                 + ", minimum=" + std::to_string(key.minimum) + ")";
        });

    py::class_<FingerprintOptions>(m, "FingerprintOptions",
        "A fingerprint specification shared with C++ and Java. FingerprintOptions(family=...) starts "
        "from the classic configuration of the family, which gives the bits of the fixed fingerprint "
        "functions, and applies the keyword arguments; every further switch extends the hashed "
        "terms. Assigning `family` later keeps the current switches. Pass it as the `options` of "
        "fingerprint() or fingerprint_from_smiles()."
          R"doc(

Example::

    >>> import smsd
    >>> options = smsd.FingerprintOptions.ecfp(radius=2, size=128)
    >>> bits = smsd.fingerprint("CCO", options=options)
    >>> assert bits and all(0 <= bit < 128 for bit in bits)
)doc")
        .def(py::init([](py::object family, py::kwargs keywords) {
            FingerprintOptions options = FingerprintOptions::of(fingerprintFamily(family));
            for (auto item : keywords)
                applyFingerprintOption(options, py::cast<std::string>(item.first), py::reinterpret_borrow<py::object>(item.second));
            return options;
        }), py::arg("family") = "ecfp")
        .def_property("family",
            [](const FingerprintOptions& options) { return options.family; },
            [](FingerprintOptions& options, py::object value) { options.family = fingerprintFamily(value); })
        .def_readwrite("radius", &FingerprintOptions::radius)
        .def_readwrite("path_length", &FingerprintOptions::pathLength)
        .def_readwrite("size", &FingerprintOptions::size)
        .def_readwrite("counts", &FingerprintOptions::counts)
        .def_readwrite("element", &FingerprintOptions::element)
        .def_readwrite("degree", &FingerprintOptions::degree)
        .def_readwrite("valence", &FingerprintOptions::valence)
        .def_readwrite("formal_charge", &FingerprintOptions::formalCharge)
        .def_readwrite("isotope", &FingerprintOptions::isotope)
        .def_readwrite("hydrogen_count", &FingerprintOptions::hydrogenCount)
        .def_readwrite("tautomer_class", &FingerprintOptions::tautomerClass)
        .def_readwrite("aromaticity", &FingerprintOptions::aromaticity)
        .def_readwrite("bond_order", &FingerprintOptions::bondOrder)
        .def_readwrite("aromatic_bond", &FingerprintOptions::aromaticBond)
        .def_readwrite("ring_bond", &FingerprintOptions::ringBond)
        .def_readwrite("ring_membership", &FingerprintOptions::ringMembership)
        .def_readwrite("ring_size", &FingerprintOptions::ringSize)
        .def_readwrite("ring_size_cap", &FingerprintOptions::ringSizeCap)
        .def_readwrite("cip_atoms", &FingerprintOptions::cipAtoms)
        .def_readwrite("cip_bonds", &FingerprintOptions::cipBonds)
        .def_property("fragments",
            [](const FingerprintOptions& options) { return options.fragments; },
            [](FingerprintOptions& options, py::object value) { options.fragments = fragmentMode(value); })
        .def_property("scaffold",
            [](const FingerprintOptions& options) { return options.scaffold; },
            [](FingerprintOptions& options, py::object value) { options.scaffold = scaffoldMode(value); })
        .def_property("keys",
            [](const FingerprintOptions& options) { return options.keys; },
            [](FingerprintOptions& options, py::object value) {
                options.keys = fingerprintKeys(value);
                if (options.family == FingerprintFamily::KEYS) options.size = static_cast<int>(options.keys.size());
            })
        .def("copy", [](const FingerprintOptions& options) { return options; },
             "A copy of these options."
          R"doc(

Example::

    >>> import smsd
    >>> options = smsd.FingerprintOptions.ecfp(size=128)
    >>> copied = options.copy()
    >>> copied.size = 256
    >>> assert options.size == 128 and copied.size == 256
)doc")
        .def_static("ecfp", &FingerprintOptions::ecfp, py::arg("radius") = 2, py::arg("size") = 2048,
                    "The classic ECFP configuration; radius 2 is ECFP4."
          R"doc(

Example::

    >>> import smsd
    >>> options = smsd.FingerprintOptions.ecfp(radius=2, size=128)
    >>> assert options.family == smsd.FingerprintFamily.ECFP and smsd.fingerprint("CCO", options=options)
)doc")
        .def_static("fcfp", &FingerprintOptions::fcfp, py::arg("radius") = 2, py::arg("size") = 2048,
                    "The classic FCFP configuration; radius 2 is FCFP4."
          R"doc(

Example::

    >>> import smsd
    >>> options = smsd.FingerprintOptions.fcfp(radius=2, size=128)
    >>> assert options.family == smsd.FingerprintFamily.FCFP and smsd.fingerprint("CCO", options=options)
)doc")
        .def_static("path", &FingerprintOptions::path, py::arg("path_length") = 7, py::arg("size") = 2048,
                    "The classic path configuration: heavy-atom paths hashed by atomic number, as path_fingerprint."
          R"doc(

Example::

    >>> import smsd
    >>> options = smsd.FingerprintOptions.path(path_length=4, size=128)
    >>> assert options.family == smsd.FingerprintFamily.PATH and smsd.fingerprint("CCO", options=options)
)doc")
        .def_static("torsion", &FingerprintOptions::torsion, py::arg("size") = 2048,
                    "The classic topological torsion configuration, as topological_torsion."
          R"doc(

Example::

    >>> import smsd
    >>> options = smsd.FingerprintOptions.torsion(size=128)
    >>> assert options.family == smsd.FingerprintFamily.TORSION and smsd.fingerprint("CCCC", options=options)
)doc")
        .def_static("maccs", &FingerprintOptions::maccs,
                    "The 166 MACCS structural keys; key k is position k - 1."
          R"doc(

Example::

    >>> import smsd
    >>> options = smsd.FingerprintOptions.maccs()
    >>> assert options.size == 166 and smsd.fingerprint("CCO", options=options)
)doc")
        .def_static("maccs_keys", [] { return FingerprintOptions::maccsKeys(); },
                    "The MACCS key definitions, in key order."
          R"doc(

Example::

    >>> import smsd
    >>> keys = smsd.FingerprintOptions.maccs_keys()
    >>> assert len(keys) == 166 and all(isinstance(key, smsd.FingerprintKey) for key in keys)
)doc")
        .def("__repr__", &describeFingerprintOptions);

    m.def("fingerprint_options",
          [](const smsd::MolGraph& mol, const FingerprintOptions& given) -> py::object {
              // The options are copied while the GIL is held, as for the searches.
              const FingerprintOptions options = given;
              if (options.counts) {
                  std::vector<int> counts;
                  { py::gil_scoped_release release;
                    counts = smsd::batch::computeFingerprintCounts(mol, options); }
                  return py::cast(nonZeroBins(counts));
              }
              std::vector<uint64_t> fp;
              { py::gil_scoped_release release;
                fp = smsd::batch::computeFingerprint(mol, options); }
              const int size = options.family == FingerprintFamily::KEYS
                             ? static_cast<int>(options.keys.size()) : options.size;
              return py::cast(smsd::bitPositions(fp, size));
          },
          py::arg("mol"), py::arg("options"),
          "The fingerprint of a graph under FingerprintOptions: the set positions, or the "
          "(position, count) pairs of the non-zero bins when options.counts is set.");

    m.def("fingerprints", &pythonFingerprints,
          py::arg("mols"), py::arg("spec"), py::arg("radius") = 2, py::arg("path_length") = 7,
          py::arg("fp_size") = 2048, py::arg("num_threads") = 0, py::arg("convert") = py::none(),
          "The fingerprint of every molecule, in input order, under the FingerprintOptions or the "
          "fingerprint() kind that spec gives. Molecules are MolGraph objects, SMILES strings parsed "
          "natively, or items that convert(item) turns into a MolGraph; num_threads=0 uses all "
          "available processors.");

    // -----------------------------------------------------------------------
    m.def("tanimoto",
          [](const std::vector<int>& fp1, const std::vector<int>& fp2, int fpSize) {
              return smsd::batch::fingerprintTanimoto(fingerprintWords(fp1, fpSize),
                                                      fingerprintWords(fp2, fpSize));
          },
          py::arg("fp1"), py::arg("fp2"), py::arg("fp_size") = 2048,
          "Tanimoto coefficient popcount(A & B) / popcount(A | B) between two fingerprints (lists of set bit positions).");
    // Set forms used by the package-level helpers: every position counts once, whatever its value.
    m.def("_set_tanimoto",
          [](const std::vector<int>& fp1, const std::vector<int>& fp2) {
              return smsd::batch::tanimoto(smsd::batch::bitCounts(fp1, fp2));
          },
          py::arg("fp1"), py::arg("fp2"),
          "Tanimoto coefficient of two sets of bit positions; two empty sets score 1.0.");
    m.def("_set_overlap",
          [](const std::vector<int>& fp1, const std::vector<int>& fp2) {
              return smsd::batch::overlap(smsd::batch::bitCounts(fp1, fp2));
          },
          py::arg("fp1"), py::arg("fp2"),
          "Overlap coefficient |A&B| / min(|A|, |B|) of two sets of bit positions; two empty sets score 1.0.");

    m.def("dice",
          [](const std::vector<int>& fp1, const std::vector<int>& fp2,
             int fpSize) {
              return smsd::batch::fingerprintDice(fingerprintWords(fp1, fpSize),
                                                  fingerprintWords(fp2, fpSize));
          },
          py::arg("fp1"), py::arg("fp2"), py::arg("fp_size") = 2048,
          "Dice similarity between two fingerprints (lists of set bit positions). "
          "Dice = 2*|A&B| / (|A| + |B|). Weights shared features more heavily "
          "than Tanimoto."
          R"doc(

Example::

    >>> import smsd
    >>> assert smsd.dice([1, 3], [3, 5], fp_size=64) == 0.5
)doc");

    m.def("cosine",
          [](const std::vector<int>& fp1, const std::vector<int>& fp2,
             int fpSize) {
              return smsd::batch::fingerprintCosine(fingerprintWords(fp1, fpSize),
                                                  fingerprintWords(fp2, fpSize));
          },
          py::arg("fp1"), py::arg("fp2"), py::arg("fp_size") = 2048,
          "Cosine similarity between two fingerprints (lists of set bit positions). "
          "Cosine = |A&B| / sqrt(|A| * |B|). Normalises for fingerprint density."
          R"doc(

Example::

    >>> import smsd
    >>> assert smsd.cosine([1, 3], [3, 5], fp_size=64) == 0.5
)doc");

    m.def("soergel",
          [](const std::vector<int>& fp1, const std::vector<int>& fp2,
             int fpSize) {
              return smsd::batch::fingerprintSoergel(fingerprintWords(fp1, fpSize),
                                                  fingerprintWords(fp2, fpSize));
          },
          py::arg("fp1"), py::arg("fp2"), py::arg("fp_size") = 2048,
          "Soergel distance between two fingerprints (lists of set bit positions). "
          "Soergel = 1 - Tanimoto. A proper metric distance for clustering."
          R"doc(

Example::

    >>> import smsd
    >>> assert smsd.soergel([1, 3], [1, 3], fp_size=64) == 0.0
)doc");

    // -----------------------------------------------------------------------
    m.def("count_tanimoto",
          [](const std::vector<int>& fp1, const std::vector<int>& fp2) {
              return smsd::batch::countTanimoto(fp1, fp2);
          },
          py::arg("fp1"), py::arg("fp2"),
          "Count-vector Tanimoto: sum(min(a,b)) / sum(max(a,b)) for integer count vectors.");

    m.def("count_dice",
          [](const std::vector<int>& fp1, const std::vector<int>& fp2) {
              return smsd::batch::countDice(fp1, fp2);
          },
          py::arg("fp1"), py::arg("fp2"),
          "Count-vector Dice: 2*sum(min(a,b)) / (sum(a) + sum(b)). "
          "Generalisation of Dice for integer count vectors.");

    m.def("count_cosine",
          [](const std::vector<int>& fp1, const std::vector<int>& fp2) {
              return smsd::batch::countCosine(fp1, fp2);
          },
          py::arg("fp1"), py::arg("fp2"),
          "Count-vector Cosine: dot(a,b) / (|a| * |b|). "
          "Generalisation of cosine similarity for integer count vectors.");

    m.def("to_hex",
          [](const std::vector<int>& fp, int fpSize) {
              return smsd::batch::toHex(fingerprintWords(fp, fpSize));
          },
          py::arg("fp"), py::arg("fp_size") = 2048,
          "Convert fingerprint (set bit positions) to hexadecimal string.");

    m.def("from_hex",
          [](const std::string& hex) {
              return smsd::bitPositions(smsd::batch::fromHex(hex));
          },
          py::arg("hex"),
          "Parse hexadecimal string back to fingerprint (set bit positions).");

    m.def("to_binary_string",
          [](const std::vector<int>& fp, int fpSize) {
              return smsd::batch::toBinaryString(fingerprintWords(fp, fpSize), fpSize);
          },
          py::arg("fp"), py::arg("fp_size") = 2048,
          "Convert fingerprint to binary string of 0s and 1s.");

    py::enum_<smsd::cip::RSLabel>(m, "RSLabel",
        "CIP tetrahedral stereo descriptor (R, S, pseudoasymmetric r, s, or NONE)"
          R"doc(

Example::

    >>> import smsd
    >>> assert smsd.assign_cip("N[C@@H](C)C(=O)O").rs_labels[1] == smsd.RSLabel.S
)doc")
        .value("NONE", smsd::cip::RSLabel::NONE)
        .value("R",    smsd::cip::RSLabel::R)
        .value("S",    smsd::cip::RSLabel::S)
        .value("r",    smsd::cip::RSLabel::r)
        .value("s",    smsd::cip::RSLabel::s)
        .export_values();

    py::enum_<smsd::cip::EZLabel>(m, "EZLabel",
        "CIP double-bond stereo descriptor (E, Z, or NONE)"
          R"doc(

Example::

    >>> import smsd
    >>> assert smsd.assign_cip("C/C=C/C").ez_bonds[0][2] == smsd.EZLabel.E
)doc")
        .value("NONE", smsd::cip::EZLabel::NONE)
        .value("E",    smsd::cip::EZLabel::E)
        .value("Z",    smsd::cip::EZLabel::Z)
        .export_values();

    py::class_<smsd::cip::CIPDescriptors>(m, "CIPDescriptors",
        "CIP stereo descriptors for all centres and double bonds"
          R"doc(

Example::

    >>> import smsd
    >>> descriptors = smsd.assign_cip("N[C@@H](C)C(=O)O")
    >>> assert descriptors.rs_labels[1] == smsd.RSLabel.S
)doc")
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
          "Ranks the ligands by Rule 1a, then those it leaves tied by every sequence rule.\n"
          "Returns RSLabel (NONE, R, S, r, s).");

    m.def("assign_ez", &smsd::cip::assignEZ,
          py::arg("mol"), py::arg("a1"), py::arg("a2"),
          "Assign CIP E/Z descriptor to a double bond between atoms a1 and a2");

    m.def("assign_cip", &smsd::cip::assignAll,
          py::arg("mol"),
          "Compute the CIP descriptors of every stereocentre and stereogenic double bond.\n"
          "Returns CIPDescriptors with rs_labels (one R/S label per atom) and ez_bonds\n"
          "(a list of (atom1, atom2, EZLabel) tuples).\n"
          "Example: assign_cip('N[C@@H](C)C(=O)O').rs_labels  # [NONE, S, NONE, ...]"
          R"doc(

Example::

    >>> import smsd
    >>> descriptors = smsd.assign_cip("N[C@@H](C)C(=O)O")
    >>> assert len(descriptors.rs_labels) == 6 and descriptors.rs_labels[1] == smsd.RSLabel.S
)doc");

    m.def("assign_cip_from_smiles", &smsd::cip::assignFromSMILES,
          py::arg("smiles"),
          "Parse SMILES and compute CIP descriptors");

    m.def("compute_sssr",
          [](const smsd::MolGraph& g) {
              return smsd::computeSSSR(g);
          },
          py::arg("mol"),
          "Return the Smallest Set of Smallest Rings (minimum cycle basis) as lists of\n"
          "atom indices, sorted by ring size.\n"
          "Example: compute_sssr('c1ccc2ccccc2c1')  # [[0, 1, 2, 3, 8, 9], [3, 4, 5, 6, 7, 8]]"
          R"doc(

Example::

    >>> import smsd
    >>> rings = smsd.compute_sssr("c1ccccc1")
    >>> assert len(rings) == 1 and set(rings[0]) == set(range(6))
)doc");

    m.def("rings_truncated", &smsd::ringsTruncated,
          py::arg("mol"),
          "Whether the bounded Horton/minimum-basis perception stopped early or failed to reach\n"
          "the graph's cycle rank. This does not report a separate relevant-cycle enumeration cut.\n"
          "The rings are the same whether or not this is asked.\n"
          "Example: rings_truncated('c1ccc2ccccc2c1')  # False"
          R"doc(

Example::

    >>> import smsd
    >>> assert not smsd.rings_truncated("c1ccccc1")
)doc");

    m.def("relevant_cycles_truncated", &smsd::relevantCyclesTruncated,
          py::arg("mol"),
          "Whether full relevant-cycle enumeration was incomplete, including its path or work\n"
          "bounds. This flag is independent of rings_truncated for the minimum-basis perception."
          R"doc(

Example::

    >>> import smsd
    >>> assert not smsd.relevant_cycles_truncated("c1ccccc1")
)doc");

    m.def("layout_sssr",
          [](const smsd::MolGraph& g) {
              return smsd::layoutSSSR(g);
          },
          py::arg("mol"),
          "Return the SSSR rings in the order used for 2D coordinate generation:\n"
          "largest ring system first, fused rings by shared-edge adjacency within a\n"
          "system, and smaller rings first within a system.\n"
          "Example: layout_sssr('c1ccc2ccccc2c1')"
          R"doc(

Example::

    >>> import smsd
    >>> rings = smsd.layout_sssr("c1ccccc1")
    >>> assert len(rings) == 1 and len(rings[0]) == 6
)doc");

    // Exact member counts may exceed 64 bits; they arrive in decimal and leave as Python int.
    const auto exactInteger = [](const std::string& decimal) {
        PyObject* value = PyLong_FromString(decimal.c_str(), nullptr, 10);
        if (value == nullptr) throw py::error_already_set();
        return py::reinterpret_steal<py::int_>(value);
    };

    py::class_<smsd::RelevantCycleFamily>(m, "RelevantCycleFamily",
          R"doc(A family of relevant cycles: its ring size, the highest-numbered atom shared by
its members (root), the atoms opposite it, one member and the exact number of members.
For an odd size, first and second are the atoms of the bond opposite root and opposite is
-1; for an even size, opposite is the atom opposite root and first and second are its ring
neighbours.

Example::

    >>> import smsd
    >>> family = smsd.compute_relevant_cycle_families("c1ccccc1").families[0]
    >>> assert (family.size, family.root, family.opposite, family.member_count) == (6, 5, 2, 1)
)doc")
        .def_readonly("size", &smsd::RelevantCycleFamily::size, "Ring size of every member.")
        .def_readonly("root", &smsd::RelevantCycleFamily::root, "The highest-numbered atom shared by every member.")
        .def_readonly("first", &smsd::RelevantCycleFamily::first,
                      "For an odd size, the lower-numbered atom of the bond opposite root; for an even size, "
                      "the lower-numbered ring neighbour of the opposite atom.")
        .def_readonly("second", &smsd::RelevantCycleFamily::second, "As first, the higher-numbered atom.")
        .def_readonly("opposite", &smsd::RelevantCycleFamily::opposite,
                      "The atom opposite root for an even ring size, otherwise -1.")
        .def_readonly("prototype", &smsd::RelevantCycleFamily::prototype,
                      "One member as atom indices, starting at root and reaching first before second.")
        .def_property_readonly("member_count",
                               [exactInteger](const smsd::RelevantCycleFamily& family) { return exactInteger(family.memberCount); },
                               "The exact number of members.")
        .def("__repr__", [](const smsd::RelevantCycleFamily& family) {
            return "RelevantCycleFamily(size=" + std::to_string(family.size) + ", root=" + std::to_string(family.root)
                + ", member_count=" + family.memberCount + ")";
        });

    py::class_<smsd::RelevantCycleFamilies>(m, "RelevantCycleFamilies",
          R"doc(Relevant cycles grouped into families, in ascending ring size, then root.
relevant_cycle_count is the exact number of members of the families listed; truncated
reports an incomplete enumeration, whose families are relevant but may not be all of them.
When the enumeration is complete, the relevant cycles and their number do not depend on
atom order; the grouping into families, and so each family's record, can.

Example::

    >>> import smsd
    >>> result = smsd.compute_relevant_cycle_families("c1ccc2ccccc2c1")
    >>> assert len(result.families) == 2 and result.relevant_cycle_count == 2 and not result.truncated
)doc")
        .def_readonly("families", &smsd::RelevantCycleFamilies::families, "The families, as a list.")
        .def_property_readonly("relevant_cycle_count",
                               [exactInteger](const smsd::RelevantCycleFamilies& families) {
                                   return exactInteger(families.relevantCycleCount);
                               },
                               "The exact number of members of the families listed.")
        .def_readonly("truncated", &smsd::RelevantCycleFamilies::truncated,
                      "Whether the enumeration was incomplete, so that families may be missing.")
        .def("__repr__", [](const smsd::RelevantCycleFamilies& families) {
            return "RelevantCycleFamilies(families=" + std::to_string(families.families.size())
                + ", relevant_cycle_count=" + families.relevantCycleCount
                + ", truncated=" + (families.truncated ? "True" : "False") + ")";
        });

    m.def("compute_relevant_cycle_families", &smsd::computeRelevantCycleFamilies,
          py::arg("mol"),
          "Return the relevant cycles grouped into families (RelevantCycleFamilies): each family\n"
          "gives its ring size, the highest-numbered atom shared by its members, the atoms\n"
          "opposite it, one member and the exact number of members; truncated reports an\n"
          "incomplete enumeration."
          R"doc(

Example::

    >>> import smsd
    >>> result = smsd.compute_relevant_cycle_families("C12C3C4C1C5C2C3C45")
    >>> assert len(result.families) == 6 and result.relevant_cycle_count == 6
)doc");

    m.def("relevant_cycle_family_of", &smsd::relevantCycleFamilyOf,
          py::arg("mol"), py::arg("families"), py::arg("ring"),
          "Return the index in families.families of the family that has ring (a list of atom\n"
          "indices in ring order) as a member, or -1 when ring is not a member of a listed family."
          R"doc(

Example::

    >>> import smsd
    >>> families = smsd.compute_relevant_cycle_families("c1ccccc1")
    >>> assert smsd.relevant_cycle_family_of("c1ccccc1", families, [0, 1, 2, 3, 4, 5]) == 0
)doc");

    m.def("relevant_cycle_family_members", &smsd::relevantCycleFamilyMembers,
          py::arg("mol"), py::arg("family"), py::arg("limit"),
          "Return up to limit members of family as lists of atom indices, in a fixed order that\n"
          "starts with the family's prototype."
          R"doc(

Example::

    >>> import smsd
    >>> family = smsd.compute_relevant_cycle_families("c1ccccc1").families[0]
    >>> assert smsd.relevant_cycle_family_members("c1ccccc1", family, 10) == [family.prototype]
)doc");

    py::class_<smsd::Point2D>(m, "Point2D",
          R"doc(A Cartesian two-dimensional coordinate.


Example::

    >>> import smsd
    >>> point = smsd.Point2D(1.0, 2.0)
    >>> assert (point.x, point.y) == (1.0, 2.0)
)doc")
        .def(py::init<>())
        .def(py::init([](double x, double y) {
            smsd::Point2D p; p.x = x; p.y = y; return p;
        }), py::arg("x"), py::arg("y"))
        .def_readwrite("x", &smsd::Point2D::x)
        .def_readwrite("y", &smsd::Point2D::y)
        .def("__repr__", [](const smsd::Point2D& p) {
            return "Point2D(" + std::to_string(p.x) + ", " + std::to_string(p.y) + ")";
        });

    // The refiners read one (x, y) pair per atom; a shorter list is padded with (0, 0).
    m.def("reduce_crossings",
          [](const smsd::MolGraph& g, const py::sequence& coordsList, int maxIter) {
              auto coords = points2D(coordsList, static_cast<size_t>(g.n));
              int result;
              {
                  py::gil_scoped_release release;
                  result = smsd::reduceCrossings(g, coords, maxIter);
              }
              return py::make_tuple(result, pointList(coords));
          },
          py::arg("mol"), py::arg("coords"), py::arg("max_iter") = 1000,
          "Reduce bond crossings in a 2D layout by optimising ring orientations.\n"
          "coords: a sequence of (x, y) pairs, one per atom; the input is not changed.\n"
          "max_iter: maximum number of refinement iterations (default 1000).\n"
          "Returns (crossings_remaining, updated_coords).\n"
          "Example: crossings, coords = reduce_crossings(mol, generate_coords_2d(mol))"
          R"doc(

Example::

    >>> import smsd
    >>> graph = smsd.parse_smiles("c1ccccc1")
    >>> crossings, coordinates = smsd.reduce_crossings(graph, smsd.generate_coords_2d(graph), max_iter=10)
    >>> assert crossings == 0 and len(coordinates) == 6
)doc");

    m.def("force_directed_layout",
          [](const smsd::MolGraph& g, const py::sequence& coordsList, int maxIter, double targetBondLength) {
              auto coords = points2D(coordsList, static_cast<size_t>(g.n));
              double stress;
              {
                  py::gil_scoped_release release;
                  stress = smsd::forceDirectedLayout(g, coords, maxIter, targetBondLength);
              }
              return py::make_tuple(stress, pointList(coords));
          },
          py::arg("mol"), py::arg("coords"),
          py::arg("max_iter") = 500, py::arg("target_bond_length") = 1.5,
          "Force-directed 2D layout minimisation.\n"
          "Iteratively moves atoms to minimise bond-length stress, non-bonded\n"
          "repulsion and crossing penalties; suitable for refining an existing layout\n"
          "or for producing crossing-free coordinates for polycyclic systems.\n"
          "coords: a sequence of (x, y) pairs, one per atom; the input is not changed.\n"
          "max_iter: maximum number of force iterations (default 500).\n"
          "target_bond_length: desired bond length in coordinate units (default 1.5).\n"
          "Returns (final_stress, updated_coords)."
          R"doc(

Example::

    >>> import smsd, math
    >>> graph = smsd.parse_smiles("CCO")
    >>> stress, coordinates = smsd.force_directed_layout(graph, smsd.generate_coords_2d(graph), max_iter=10)
    >>> assert math.isfinite(stress) and len(coordinates) == 3
)doc");

    m.def("stress_majorisation",
          [](const smsd::MolGraph& g, const py::sequence& coordsList, int maxIter, double targetBondLength) {
              auto coords = points2D(coordsList, static_cast<size_t>(g.n));
              double stress;
              {
                  py::gil_scoped_release release;
                  stress = smsd::stressMajorisation(g, coords, maxIter, targetBondLength);
              }
              return py::make_tuple(stress, pointList(coords));
          },
          py::arg("mol"), py::arg("coords"),
          py::arg("max_iter") = 300, py::arg("target_bond_length") = 1.5,
          "Stress majorisation for 2D layout.\n"
          "Minimises the weighted stress sum_ij w_ij * (d_ij - D_ij)^2, where D_ij is the\n"
          "graph distance scaled by target_bond_length, converging monotonically to a\n"
          "local minimum.\n"
          "coords: a sequence of (x, y) pairs, one per atom; the input is not changed.\n"
          "max_iter: maximum number of majorisation iterations (default 300).\n"
          "target_bond_length: desired bond length in coordinate units (default 1.5).\n"
          "Returns (final_stress, updated_coords)."
          R"doc(

Example::

    >>> import smsd, math
    >>> graph = smsd.parse_smiles("CCO")
    >>> stress, coordinates = smsd.stress_majorisation(graph, smsd.generate_coords_2d(graph), max_iter=10)
    >>> assert math.isfinite(stress) and len(coordinates) == 3
)doc");

    m.def("match_template",
          [](const smsd::MolGraph& g, double targetBondLength) {
              return pointList(smsd::matchTemplate(g, targetBondLength));
          },
          py::arg("mol"), py::arg("target_bond_length") = 1.5,
          "Check whether a molecule matches a known scaffold template.\n"
          "Returns the template's (x, y) coordinate tuples scaled to target_bond_length\n"
          "if it does, and an empty list otherwise.\n"
          "Templates: benzene, naphthalene, indole, purine, steroid, morphinan,\n"
          "quinoline, biphenyl, cyclohexane, piperidine.\n"
          "Example: match_template('c1ccccc1')  # six coordinate pairs"
          R"doc(

Example::

    >>> import smsd
    >>> assert len(smsd.match_template("c1ccccc1")) == 6
)doc");

    // The mapping checks run bounded searches, so they release the GIL as the searches do; the
    // match and options are copied while it is held, so Python threads that change them cannot
    // race with the check.
    mcs.def_static("validate_mapping",
          [](const smsd::MolGraph& query, const smsd::MolGraph& target, const smsd::Match& match,
             const smsd::MCS::Options& options) {
              const smsd::Match copiedMatch = match;
              const smsd::MCS::Options copied = options;
              py::gil_scoped_release release;
              return smsd::MCS::validateMapping(query, target, copiedMatch, copied);
          },
          py::arg("query"), py::arg("target"), py::arg("match"), py::arg("options") = smsd::MCS::Options());
    mcs.def_static("validate_mapping",
          [](const smsd::MolGraph& query, const smsd::MolGraph& target, const std::map<int, int>& mapping,
             const smsd::MCS::Options& options) {
              const smsd::MCS::Options copied = options;
              py::gil_scoped_release release;
              return smsd::MCS::validateMapping(query, target, mapping, copied);
          },
          py::arg("query"), py::arg("target"), py::arg("match"), py::arg("options") = smsd::MCS::Options());
    mcs.def_static("is_mapping_maximal",
          [](const smsd::MolGraph& query, const smsd::MolGraph& target, const smsd::Match& match,
             const smsd::MCS::Options& options) {
              const smsd::Match copiedMatch = match;
              const smsd::MCS::Options copied = options;
              py::gil_scoped_release release;
              return smsd::MCS::isMappingMaximal(query, target, copiedMatch, copied);
          },
          py::arg("query"), py::arg("target"), py::arg("match"), py::arg("options") = smsd::MCS::Options());
    mcs.def_static("is_mapping_maximal",
          [](const smsd::MolGraph& query, const smsd::MolGraph& target, const std::map<int, int>& mapping,
             const smsd::MCS::Options& options) {
              const smsd::MCS::Options copied = options;
              py::gil_scoped_release release;
              return smsd::MCS::isMappingMaximal(query, target, mapping, copied);
          },
          py::arg("query"), py::arg("target"), py::arg("match"), py::arg("options") = smsd::MCS::Options());

    substructure.def_static("decompose_r_groups",
          [](const smsd::MolGraph& core, const std::vector<smsd::MolGraph>& molecules,
             const smsd::Substructure::Options& options) {
              std::vector<smsd::RGroupResult> results;
              { py::gil_scoped_release release;
                results = smsd::Substructure::decomposeRGroups(core, molecules, options); }
              py::list out;
              for (auto& result : results) {
                  py::dict item;
                  item["core"] = result.core;
                  py::dict groups;
                  for (auto& group : result.rgroups) groups[py::cast(group.first)] = group.second;
                  item["rgroups"] = groups;
                  item["search"] = result.search;
                  out.append(item);
              }
              return out;
          }, py::arg("core"), py::arg("targets"), py::arg("options") = smsd::Substructure::Options());

    m.def("extract_subgraph",
          &smsd::extractSubgraph,
          py::arg("mol"), py::arg("atom_indices"),
          "Extract a subgraph containing only the specified atoms.\n"
          "Returns a new MolGraph with atoms re-indexed from 0.\n"
          "Preserves all atom and bond properties, a stereocentre whose neighbours all\n"
          "stay and the E/Z configuration of a double bond whose ends keep the neighbours\n"
          "that define it."
          R"doc(

Example::

    >>> import smsd
    >>> fragment = smsd.extract_subgraph("CCO", [0, 1])
    >>> assert fragment.n == 2 and fragment.atomic_num == [6, 6]
)doc");

    m.def("murcko_scaffold",
          &smsd::murckoScaffold,
          py::arg("mol"),
          "Compute the Bemis-Murcko scaffold (ring systems and linkers).\n"
          "Keeps the ring atoms, the atoms on shortest paths between rings\n"
          "and every atom joined to these by a double or triple bond;\n"
          "each stripped substituent is replaced by hydrogen. A stereocentre\n"
          "keeps its configuration when all its neighbours stay, and a\n"
          "double bond when each end keeps a ligand."
          R"doc(

Example::

    >>> import smsd
    >>> scaffold = smsd.murcko_scaffold("CCc1ccccc1")
    >>> assert scaffold.n == 6 and all(scaffold.ring)
)doc");

    m.def("count_components",
          static_cast<int (*)(const smsd::MolGraph&)>(&smsd::search_detail::countComponents),
          py::arg("mol"),
          "Count the number of connected components in a molecule.\n"
          "Returns 1 for normal connected molecules, >1 for salts/mixtures."
          R"doc(

Example::

    >>> import smsd
    >>> assert smsd.count_components("CCO.[Na+]") == 2
)doc");

    m.def("split_components",
          &smsd::search_detail::splitComponents,
          py::arg("mol"),
          "Split a molecule into its connected components.\n"
          "Returns a list of MolGraphs, one per connected component.\n"
          "Useful for separating salts, counterions, or mixture components."
          R"doc(

Example::

    >>> import smsd
    >>> assert sorted(graph.n for graph in smsd.split_components("CCO.[Na+]")) == [1, 3]
)doc");

    m.def("same_canonical_graph",
          &smsd::search_detail::sameCanonicalGraph,
          py::arg("mol1"), py::arg("mol2"),
          "Check if two molecules are canonically identical.\n"
          "Compares canonical SMILES — fast graph isomorphism test.\n"
          "Does not consider stereochemistry unless encoded in SMILES."
          R"doc(

Example::

    >>> import smsd
    >>> assert smsd.same_canonical_graph(smsd.parse_smiles("CCO"), smsd.parse_smiles("OCC"))
)doc");

    m.def("read_mol_file",
          &smsd::readMolFile,
          py::arg("filename"),
          py::call_guard<py::gil_scoped_release>(),
          "Read a MOL file from disk into a MolGraph.\n"
          "Supports V2000 and V3000 formats.\n"
          "Raises ValueError if the file cannot be opened or read.");

    m.def("read_sdf",
          &smsd::readSDF,
          py::arg("filename"),
          py::call_guard<py::gil_scoped_release>(),
          "Read all molecules from an SDF file.\n"
          "Returns a list of MolGraphs. Malformed entries become empty graphs.\n"
          "The whole file is loaded into memory; for very large files, read the\n"
          "records one at a time with read_mol_block().");

    m.def("write_sdf",
          &smsd::writeSDF,
          py::arg("molecules"), py::arg("filename"),
          py::call_guard<py::gil_scoped_release>(),
          "Write a list of MolGraphs to an SDF file.\n"
          "Each molecule is written as a V2000 MOL block followed by $$$$.");

    m.def("classify_pharmacophore",
          &smsd::batch::detail::classifyPharmacophore,
          py::arg("mol"), py::arg("atom_index"),
          "Classify an atom into pharmacophoric feature classes (FCFP invariant).\n"
          "Returns a bitmask: bit 0=H-bond donor, 1=H-bond acceptor,\n"
          "2=positive ionisable, 3=negative ionisable, 4=aromatic, 5=hydrophobic.\n"
          "Based on Rogers & Hahn 2010 ECFP/FCFP definitions."
          R"doc(

Example::

    >>> import smsd
    >>> oxygen = smsd.classify_pharmacophore(smsd.parse_smiles("CCO"), 2)
    >>> assert oxygen & 1  # hydroxyl donor
)doc");

    m.def("implicit_h",
          [](int atomicNum, int bondOrderSum, int formalCharge, bool aromatic) {
              return smsd::detail::computeImplicitH(atomicNum, aromatic, bondOrderSum, formalCharge);
          },
          py::arg("atomic_num"), py::arg("bond_order_sum"), py::arg("formal_charge"),
          py::arg("aromatic") = false,
          "Compute implicit hydrogen count for an atom with the SMILES valence model.\n"
          "Uses element- and charge-dependent default valences; an aromatic atom counts\n"
          "one extra bond order. Returns 0 if no default valence fits."
          R"doc(

Example::

    >>> import smsd
    >>> assert smsd.implicit_h(6, 1, 0) == 3  # methyl carbon
)doc");

    // The quality measures read exactly one (x, y) pair per atom, ignoring any beyond the last atom.
    m.def("count_crossings",
          [](const smsd::MolGraph& g, const py::sequence& coordsList) {
              return smsd::detail_layout::countCrossings(g, points2D(coordsList, static_cast<size_t>(g.n), true));
          },
          py::arg("mol"), py::arg("coords"),
          "Count the number of bond crossings in a 2D layout.\n"
          "coords: a sequence of (x, y) pairs, one per atom.\n"
          "Returns the number of pairs of bonds that cross each other; 0 is optimal\n"
          "for a planar graph."
          R"doc(

Example::

    >>> import smsd
    >>> graph = smsd.parse_smiles("c1ccccc1")
    >>> assert smsd.count_crossings(graph, smsd.generate_coords_2d(graph)) == 0
)doc");

    m.def("is_degenerate_layout",
          [](const py::sequence& coordsList) {
              return smsd::detail_layout::isDegenerateLayout(points2D(coordsList));
          },
          py::arg("coords"),
          "Check if a 2D layout is degenerate (all points collinear or stacked).\n"
          "Returns True if the layout has zero width or height.\n"
          "Use before reduce_crossings() to detect bad initial layouts."
          R"doc(

Example::

    >>> import smsd
    >>> assert smsd.is_degenerate_layout([(0, 0), (1, 0), (2, 0)])
    >>> assert not smsd.is_degenerate_layout([(0, 0), (1, 0), (0, 1)])
)doc");

    m.def("counts_to_array",
          &smsd::batch::countsToArray,
          py::arg("counts"), py::arg("fp_size") = 2048,
          "Convert a count fingerprint (dict) to a dense integer array.\n"
          "Returns a list of length fp_size where index i = count of feature i.\n"
          "Useful for ML pipelines that require fixed-width feature vectors.");

    py::class_<smsd::Point3D>(m, "Point3D",
          R"doc(A Cartesian three-dimensional coordinate.


Example::

    >>> import smsd
    >>> point = smsd.Point3D(1.0, 2.0, 3.0)
    >>> assert (point.x, point.y, point.z) == (1.0, 2.0, 3.0)
)doc")
        .def(py::init<>())
        .def(py::init([](double x, double y, double z) {
            return smsd::Point3D{x, y, z};
        }), py::arg("x"), py::arg("y"), py::arg("z"))
        .def_readwrite("x", &smsd::Point3D::x)
        .def_readwrite("y", &smsd::Point3D::y)
        .def_readwrite("z", &smsd::Point3D::z)
        .def("__repr__", [](const smsd::Point3D& p) {
            return "Point3D(" + std::to_string(p.x) + ", "
                   + std::to_string(p.y) + ", " + std::to_string(p.z) + ")";
        });

    m.def("generate_coords_2d",
          [](const smsd::MolGraph& g, double targetBondLength) {
              std::vector<smsd::Point2D> coords;
              { py::gil_scoped_release release;
                coords = smsd::generateCoords2D(g, targetBondLength); }
              return pointList(coords);
          },
          py::arg("mol"), py::arg("target_bond_length") = 1.5,
          "Generate 2D coordinates for a molecule.\n"
          "Returns a list of (x, y) tuples, one per atom."
          R"doc(

Example::

    >>> import smsd, math
    >>> coordinates = smsd.generate_coords_2d("CCO")
    >>> assert len(coordinates) == 3 and all(math.isfinite(value) for point in coordinates for value in point)
)doc");

    m.def("generate_coords_3d",
          [](const smsd::MolGraph& g, double targetBondLength) {
              std::vector<smsd::Point3D> coords;
              { py::gil_scoped_release release;
                coords = smsd::generateCoords3D(g, targetBondLength); }
              py::list out(coords.size());
              for (size_t i = 0; i < coords.size(); i++)
                  out[i] = py::make_tuple(coords[i].x, coords[i].y, coords[i].z);
              return out;
          },
          py::arg("mol"), py::arg("target_bond_length") = 1.5,
          "Generate an approximate 3D geometry for a molecule.\n"
          "Returns a list of (x, y, z) tuples, one per atom. For production conformers,\n"
          "dedicated tools such as ETKDG or MMFF are recommended; this gives a reasonable\n"
          "starting geometry."
          R"doc(

Example::

    >>> import smsd, math
    >>> coordinates = smsd.generate_coords_3d("CCO")
    >>> assert len(coordinates) == 3 and all(len(point) == 3 for point in coordinates)
    >>> assert all(math.isfinite(value) for point in coordinates for value in point)
)doc");

    m.def("layout_quality",
          [](const smsd::MolGraph& g, const py::sequence& coordsList, double targetBondLength) {
              return smsd::layoutQuality(g, points2D(coordsList, static_cast<size_t>(g.n), true), targetBondLength);
          },
          py::arg("mol"), py::arg("coords"), py::arg("target_bond_length") = 1.5,
          "Compute a layout quality score (0.0 = perfect).\n"
          "Combines bond length uniformity, overlap count and crossing count, so that\n"
          "different layout strategies can be compared."
          R"doc(

Example::

    >>> import smsd, math
    >>> graph = smsd.parse_smiles("CCO")
    >>> quality = smsd.layout_quality(graph, smsd.generate_coords_2d(graph))
    >>> assert math.isfinite(quality) and quality >= 0
)doc");

    m.def("resolve_overlaps",
          [](const py::sequence& coordsList, double threshold, int maxIter) {
              auto coords = points2D(coordsList);
              int remaining = smsd::resolveOverlaps(coords, threshold, maxIter);
              return py::make_tuple(remaining, pointList(coords));
          },
          py::arg("coords"), py::arg("threshold") = 0.3, py::arg("max_iter") = 100,
          "Resolve atom-atom overlaps by pushing the atoms apart.\n"
          "Returns (remaining_overlaps, updated_coords)."
          R"doc(

Example::

    >>> import smsd
    >>> remaining, coordinates = smsd.resolve_overlaps([(0, 0), (2, 0)], threshold=0.3)
    >>> assert remaining == 0 and len(coordinates) == 2
)doc");

    // The transforms take exactly the pairs given; a pair shorter than two raises IndexError.
    const auto transformed = [](const py::sequence& coordsList) {
        return points2D(coordsList, 0, false, ShortPairs::Reject);
    };

    m.def("translate_2d",
          [transformed](const py::sequence& coordsList, double dx, double dy) {
              auto coords = transformed(coordsList);
              smsd::transform::translate2D(coords, dx, dy);
              return pointList(coords);
          },
          py::arg("coords"), py::arg("dx"), py::arg("dy"),
          "Translate 2D coordinates by (dx, dy)."
          R"doc(

Example::

    >>> import smsd
    >>> assert smsd.translate_2d([(0, 0), (1, 2)], 3, -1) == [(3.0, -1.0), (4.0, 1.0)]
)doc");

    m.def("rotate_2d",
          [transformed](const py::sequence& coordsList, double angle) {
              auto coords = transformed(coordsList);
              smsd::transform::rotate2D(coords, angle);
              return pointList(coords);
          },
          py::arg("coords"), py::arg("angle"),
          "Rotate 2D coordinates by angle (radians) about the centroid."
          R"doc(

Example::

    >>> import smsd, math
    >>> coordinates = smsd.rotate_2d([(-1, 0), (1, 0)], math.pi / 2)
    >>> assert abs(coordinates[0][0]) < 1e-9 and abs(coordinates[0][1] + 1) < 1e-9
)doc");

    m.def("scale_2d",
          [transformed](const py::sequence& coordsList, double factor) {
              auto coords = transformed(coordsList);
              smsd::transform::scale2D(coords, factor);
              return pointList(coords);
          },
          py::arg("coords"), py::arg("factor"),
          "Scale 2D coordinates uniformly about the centroid."
          R"doc(

Example::

    >>> import smsd
    >>> assert smsd.scale_2d([(-1, 0), (1, 0)], 2) == [(-2.0, 0.0), (2.0, 0.0)]
)doc");

    m.def("mirror_x",
          [transformed](const py::sequence& coordsList) {
              auto coords = transformed(coordsList);
              smsd::transform::mirrorX(coords);
              return pointList(coords);
          },
          py::arg("coords"),
          "Mirror 2D coordinates about the X axis (flip vertically)."
          R"doc(

Example::

    >>> import smsd
    >>> assert smsd.mirror_x([(1, 2), (3, 4)]) == [(1.0, 4.0), (3.0, 2.0)]
)doc");

    m.def("mirror_y",
          [transformed](const py::sequence& coordsList) {
              auto coords = transformed(coordsList);
              smsd::transform::mirrorY(coords);
              return pointList(coords);
          },
          py::arg("coords"),
          "Mirror 2D coordinates about the Y axis (flip horizontally)."
          R"doc(

Example::

    >>> import smsd
    >>> assert smsd.mirror_y([(1, 2), (3, 4)]) == [(3.0, 2.0), (1.0, 4.0)]
)doc");

    m.def("center_2d",
          [transformed](const py::sequence& coordsList) {
              auto coords = transformed(coordsList);
              smsd::transform::center2D(coords);
              return pointList(coords);
          },
          py::arg("coords"),
          "Centre 2D coordinates on the origin."
          R"doc(

Example::

    >>> import smsd
    >>> assert smsd.center_2d([(1, 2), (3, 4)]) == [(-1.0, -1.0), (1.0, 1.0)]
)doc");

    m.def("align_2d",
          [transformed](const py::sequence& coordsList, const py::sequence& refList) {
              auto coords = transformed(coordsList);
              const auto ref = transformed(refList);
              double rmsd = smsd::transform::align2D(coords, ref);
              return py::make_tuple(rmsd, pointList(coords));
          },
          py::arg("coords"), py::arg("reference"),
          "Align 2D coordinates to a reference by the optimal rotation.\n"
          "Returns (RMSD, aligned_coords)."
          R"doc(

Example::

    >>> import smsd
    >>> rmsd, coordinates = smsd.align_2d([(1, 0), (2, 0)], [(0, 0), (1, 0)])
    >>> assert rmsd < 1e-9 and len(coordinates) == 2
)doc");

    m.def("bounding_box_2d",
          [transformed](const py::sequence& coordsList) {
              auto bb = smsd::transform::boundingBox2D(transformed(coordsList));
              return py::make_tuple(bb[0], bb[1], bb[2], bb[3]);
          },
          py::arg("coords"),
          "Compute the bounding box. Returns (minX, minY, maxX, maxY)."
          R"doc(

Example::

    >>> import smsd
    >>> assert smsd.bounding_box_2d([(1, 2), (-3, 4)]) == (-3.0, 2.0, 1.0, 4.0)
)doc");

    m.def("normalise_bond_length",
          [transformed](const smsd::MolGraph& g, const py::sequence& coordsList, double target) {
              auto coords = transformed(coordsList);
              smsd::transform::normaliseBondLength(g, coords, target);
              return pointList(coords);
          },
          py::arg("mol"), py::arg("coords"), py::arg("target") = 1.5,
          "Normalise bond lengths to the target value."
          R"doc(

Example::

    >>> import smsd, math
    >>> coordinates = smsd.normalise_bond_length("CC", [(0, 0), (4, 0)], target=1.5)
    >>> assert abs(math.dist(*coordinates) - 1.5) < 1e-9
)doc");

    m.def("canonical_orientation",
          [transformed](const smsd::MolGraph& g, const py::sequence& coordsList) {
              auto coords = transformed(coordsList);
              smsd::transform::canonicalOrientation(g, coords);
              return pointList(coords);
          },
          py::arg("mol"), py::arg("coords"),
          "Rotate to the canonical orientation (30-degree grid alignment)."
          R"doc(

Example::

    >>> import smsd
    >>> graph = smsd.parse_smiles("CCO")
    >>> assert len(smsd.canonical_orientation(graph, smsd.generate_coords_2d(graph))) == 3
)doc");

    py::class_<smsd::DepictOptions>(m, "DepictOptions",
        "SVG depiction options for bond, atom-label and highlight styling."
          R"doc(

Example::

    >>> import smsd
    >>> options = smsd.DepictOptions()
    >>> options.width, options.height = 400, 200
    >>> svg = smsd.depict_svg("CCO", opts=options)
    >>> assert 'width="400"' in svg and 'height="200"' in svg
)doc")
        .def(py::init<>())
        .def_readwrite("bond_length",       &smsd::DepictOptions::bondLength,
            "Bond length in pixels (default 30). All ACS proportions scale from this.")
        .def_readwrite("line_width",         &smsd::DepictOptions::lineWidth,
            "Bond line width (0 = auto from ACS ratio).")
        .def_readwrite("bold_width",         &smsd::DepictOptions::boldWidth,
            "Wedge bond max width (0 = auto from ACS ratio).")
        .def_readwrite("bond_spacing",       &smsd::DepictOptions::bondSpacing,
            "Double bond offset distance (0 = auto, 18% of bond_length).")
        .def_readwrite("font_size",          &smsd::DepictOptions::fontSize,
            "Atom label font size (0 = auto from ACS ratio).")
        .def_readwrite("subscript_scale",    &smsd::DepictOptions::subscriptScale,
            "Subscript/superscript size relative to font_size (default 0.70).")
        .def_readwrite("show_carbon_labels", &smsd::DepictOptions::showCarbonLabels,
            "Show 'C' labels at carbon positions (default False).")
        .def_readwrite("show_atom_indices",  &smsd::DepictOptions::showAtomIndices,
            "Show atom index numbers (default False).")
        .def_readwrite("show_map_numbers",   &smsd::DepictOptions::showMapNumbers,
            "Show atom-atom map numbers (default True).")
        .def_readwrite("highlight_radius",   &smsd::DepictOptions::highlightRadius,
            "Atom highlight circle radius (0 = auto, 40% of bond_length).")
        .def_readwrite("highlight_opacity",  &smsd::DepictOptions::highlightOpacity,
            "Highlight fill opacity (default 0.30).")
        .def_readwrite("match_width",        &smsd::DepictOptions::matchWidth,
            "Highlighted bond width (0 = auto, 2.2x line_width).")
        .def_readwrite("padding",            &smsd::DepictOptions::padding,
            "SVG padding in pixels (default 30).")
        .def_readwrite("width",              &smsd::DepictOptions::width,
            "Fixed SVG width (0 = auto-fit to molecule).")
        .def_readwrite("height",             &smsd::DepictOptions::height,
            "Fixed SVG height (0 = auto-fit to molecule).")
        .def_readwrite("font_family",        &smsd::DepictOptions::fontFamily,
            "CSS font family (default 'Arial, Helvetica, sans-serif').");

    m.def("depict",
          [](const smsd::MolGraph& g, const smsd::DepictOptions& opts) {
              return smsd::depict(g, opts);
          },
          py::arg("mol"), py::arg("opts") = smsd::DepictOptions(),
          "Render a molecule as an SVG string with automatic layout.");

    m.def("depict_with_mapping",
          [](const smsd::MolGraph& g, const std::map<int,int>& mapping,
             const smsd::DepictOptions& opts) {
              return smsd::depictWithMapping(g, mapping, opts);
          },
          py::arg("mol"), py::arg("mapping"), py::arg("opts") = smsd::DepictOptions(),
          "Render a molecule with MCS/substructure atoms highlighted and numbered.");

    m.def("depict_pair",
          [](const smsd::MolGraph& g1, const smsd::MolGraph& g2,
             const std::map<int,int>& mapping, const smsd::DepictOptions& opts) {
              return smsd::depictPair(g1, g2, mapping, opts);
          },
          py::arg("mol1"), py::arg("mol2"), py::arg("mapping"),
          py::arg("opts") = smsd::DepictOptions(),
          "Render two molecules side-by-side with MCS mapping highlighted.");

}
