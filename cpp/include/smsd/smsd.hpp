/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms. *
 * SMSD — Substructure & Maximum Common Substructure search for chemical graphs.
 * Single-include header: pulls in all SMSD components.
 */
#pragma once

#include "smsd/mol_graph.hpp"
#include "smsd/ring_finder.hpp"
#include "smsd/vf2pp.hpp"
#include "smsd/mcs.hpp"
#include "smsd/rascal.hpp"
#include "smsd/batch.hpp"
#include "smsd/smiles_parser.hpp"
#include "smsd/mol_reader.hpp"
#include "smsd/gpu.hpp"
#include "smsd/cip.hpp"
#include "smsd/layout.hpp"
#include "smsd/depict.hpp"
#include "smsd/clique_solver.hpp"
#include "smsd/hungarian.hpp"
#include "smsd/periodic_table.hpp"
#include "smsd/scaffold_library.hpp"
