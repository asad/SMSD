# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
# Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
# See the NOTICE file for attribution, trademark, and algorithm IP terms.
"""Write an MCS highlight to stdout and its proof status to stderr."""

import sys

import smsd
from smsd.highlight import mcs_svg

query, target = "c1ccccc1O", "c1ccccc1N"
options = smsd.MCS.Options()
options.timeout_ms = 1000
options.node_limit = 100000
result = smsd.MCS.find(query).with_options(options).match(target)
print(f"{result.status.name}: optimal={result.optimal}, complete={result.complete}", file=sys.stderr)
print(mcs_svg(query, target, result, width=600, height=300))
