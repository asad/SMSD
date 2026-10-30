# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
# Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
# See the NOTICE file for attribution, trademark, and algorithm IP terms.
"""Basis-only batch context preserves single causes and ordinary first errors."""

import pytest
import smsd


BASIS_CAUSE = ("Exact ring basis unavailable: ring perception reached its work limit before the minimum cycle"
               " basis was complete")
UNSETTLED = "a substructure key matched more than 100000 embeddings before its count was settled"
# A 32-layer hexagonal prism stack, whose exact basis the ring perception completes.
STACK_32_SMILES = "C12C3C4C5C6C1C17C28C32C43C54C61C15C46C34C23C82C71C17C28C32C43C64C51C15C46C34C23C82C71C17C28C32C43C64C51C15C46C34C23C82C71C17C28C32C43C64C51C15C46C34C23C82C71C17C28C32C43C64C51C15C46C34C23C82C71C17C28C32C43C64C51C15C46C34C23C82C71C17C28C32C43C64C51C15C46C34C23C82C71C17C28C32C43C64C51C15C46C34C23C82C71C17C28C32C43C64C51C15C46C34C23C82C71C17C28C32C43C64C51C15C46C34C23C82C71C17C28C32C43C64C51C15C46C34C23C82C71C17C28C32C43C64C51C15C46C34C23C82C71C17C28C32C43C64C51C15C46C34C23C82C71C17C28C32C43C64C51C15C46C34C23C82C71C17C28C32C43C64C51C15C46C34C23C82C71C1C2C3C4C6C51"
# An 80-layer stack: a synthetic cage whose ring system exceeds the ring-perception work limit.
CAGE_SMILES = "C12C3C4C5C6C1C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C1C6C5C4C3C21"


@pytest.fixture(scope="module")
def graphs():
    ordinary = smsd.parse_smiles("CCO")
    incomplete = smsd.parse_smiles(CAGE_SMILES)
    assert not smsd.rings_truncated(ordinary)
    assert smsd.rings_truncated(incomplete)
    return ordinary, incomplete


def error_of(call):
    try:
        call()
    except Exception as error:
        return error
    raise AssertionError("The fixture must raise an exception")


@pytest.mark.parametrize("counts", [False, True])
@pytest.mark.parametrize("threads", [1, 2, 0])
def test_basis_batch_names_first_input_and_preserves_single_cause(graphs, counts, threads):
    ordinary, incomplete = graphs
    options = smsd.FingerprintOptions("ecfp", ring_size=True, counts=counts)
    smsd.fingerprint(ordinary, options=options)
    single = error_of(lambda: smsd.fingerprint(incomplete, options=options))
    assert type(single) is smsd.RingBasisUnavailable and isinstance(single, RuntimeError)
    assert str(single) == BASIS_CAUSE
    assert single.__cause__ is None
    molecules = [ordinary] * 32
    molecules[9] = incomplete
    molecules[17] = "C1CC"
    actual = error_of(lambda: smsd.fingerprints(molecules, options=options, num_threads=threads))
    assert type(actual) is smsd.RingBasisUnavailable
    assert str(actual) == BASIS_CAUSE + " (input index 9, zero-based)"
    assert isinstance(actual.__cause__, RuntimeError)
    assert str(actual.__cause__) == str(single)


@pytest.mark.parametrize("counts", [False, True])
@pytest.mark.parametrize("threads", [1, 2, 0])
@pytest.mark.parametrize("earlier", ["C(C", None])
def test_earlier_ordinary_error_is_not_wrapped(graphs, counts, threads, earlier):
    ordinary, incomplete = graphs
    options = smsd.FingerprintOptions("ecfp", ring_size=True, counts=counts)
    expected = error_of(lambda: smsd.fingerprint(earlier, options=options))
    molecules = [ordinary] * 32
    molecules[3], molecules[9], molecules[17] = earlier, incomplete, "C1CC"
    actual = error_of(lambda: smsd.fingerprints(molecules, options=options, num_threads=threads))
    assert (type(actual), str(actual)) == (type(expected), str(expected))
    assert actual.__cause__ is None


@pytest.mark.parametrize("counts", [False, True])
@pytest.mark.parametrize("threads", [1, 2, 0])
def test_ring_counting_key_names_the_ring_cause_and_its_input(graphs, counts, threads):
    ordinary, incomplete = graphs
    options = smsd.FingerprintOptions("keys", keys=["[R2]"], counts=counts)
    expected = error_of(lambda: smsd.fingerprint(incomplete, options=options))
    assert type(expected) is smsd.RingBasisUnavailable
    assert str(expected) == BASIS_CAUSE
    molecules = [ordinary] * 32
    molecules[9], molecules[17] = incomplete, "C1CC"
    actual = error_of(lambda: smsd.fingerprints(molecules, options=options, num_threads=threads))
    assert type(actual) is smsd.RingBasisUnavailable
    assert str(actual) == BASIS_CAUSE + " (input index 9, zero-based)"
    assert str(actual.__cause__) == BASIS_CAUSE


def test_stack_of_32_layers_has_an_exact_basis():
    stack = smsd.parse_smiles(STACK_32_SMILES)
    assert not smsd.rings_truncated(stack)
    smsd.fingerprint(stack, options=smsd.FingerprintOptions("ecfp", ring_size=True))


@pytest.mark.parametrize("threads", [1, 2, 0])
def test_embedding_limit_settles_a_met_minimum_and_refuses_an_open_count(threads):
    stack = smsd.parse_smiles(STACK_32_SMILES)
    path = "C~C~C~C~C~C~C~C~C~C~C~C"
    smsd.fingerprint(stack, options=smsd.FingerprintOptions("keys", keys=[path]))
    counts = smsd.FingerprintOptions("keys", keys=[path], counts=True)
    expected = error_of(lambda: smsd.fingerprint(stack, options=counts))
    assert type(expected) is RuntimeError and str(expected) == UNSETTLED
    molecules = [smsd.parse_smiles("CCO")] * 32
    molecules[9] = stack
    actual = error_of(lambda: smsd.fingerprints(molecules, options=counts, num_threads=threads))
    assert type(actual) is RuntimeError
    assert str(actual) == UNSETTLED + " (input index 9, zero-based)"
