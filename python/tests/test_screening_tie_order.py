# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
# Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
# See the NOTICE file for attribution, trademark, and algorithm IP terms.
"""Ranked fingerprint ties retain target input order."""
import smsd


def test_batch_similarity_ties_keep_input_order():
    ranked = smsd.batch_similarity("CC", ["CO", "CC", "CO", "CC", "CO"])
    assert [index for index, score in ranked] == [1, 3, 0, 2, 4]
    assert ranked[:2] == [(1, 1.0), (3, 1.0)]
    assert ranked[2][1] == ranked[3][1] == ranked[4][1] < 1.0
