# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
# Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
# See the NOTICE file for attribution, trademark, and algorithm IP terms.
"""Skip the release tooling tests on interpreters without the standard tomllib module."""

import sys

# The package supports Python 3.9 and later, but the release scripts read
# pyproject.toml with tomllib, which Python provides from 3.11.
if sys.version_info < (3, 11):
    collect_ignore = ["test_collect_release_wheels.py", "test_release_manifests.py",
                      "test_repair_windows_wheel.py", "test_supply_chain_pins.py"]
