# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
# Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
# See the NOTICE file for attribution, trademark, and algorithm IP terms.
"""Compute NumPy feature rows from native molecular fingerprints."""

import smsd
from smsd.adapters import fingerprints_numpy

molecules = ["CCO", "CCN", "c1ccccc1O"]
features = fingerprints_numpy(molecules, kind="ecfp", fp_size=128)
assert features.shape == (3, 128)
assert str(features.dtype) == "uint8"

options = smsd.FingerprintOptions(family="ecfp", size=128, counts=True)
counts = fingerprints_numpy(molecules, options=options)
assert counts.shape == (3, 128)
assert str(counts.dtype) == "int64"
print(features)
print(counts)
