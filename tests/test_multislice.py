import cmath
from dataclasses import FrozenInstanceError
import math
from pathlib import Path

import numpy as np
import pytest

from qcbed.multislice import (
    MultisliceResult, legacy_slice_sequence_to_zero_based, unit_incident_wave,
    beam_indices_to_row_col, beam_indices_to_legacy_ib, extract_beam_amplitudes,
    beam_intensities, simulate_slice_sequence,
)
from qcbed.phase_grating import phase_grating_from_potential
from qcbed.propagator import PropagationParameters, propagation_kernel
from qcbed.reciprocal import (DirectLattice2D, ReciprocalLattice2D,
                             ReciprocalIndexGrid, ReciprocalVectorGrid)


def run(gratings, kernels, sequence, beams=((0, 0),), ordering="unshifted", incident=None):
    return simulate_slice_sequence(gratings, kernels, sequence, beam_indices_hk=beams,
                                   reciprocal_ordering=ordering, propagator_ordering=ordering,
                                   incident_wave=incident)


def test_legacy_sequence_conversion_and_region_blocks():
    assert legacy_slice_sequence_to_zero_based([1, 2, 1, 2, 3, 3], 3) == (0, 1, 0, 1, 2, 2)
    assert legacy_slice_sequence_to_zero_based(np.array([2, 1], dtype=np.int64), np.int64(2)) == (1, 0)


@pytest.mark.parametrize("bad", [[], [0], [-1], [3], [1, 0], [True], [np.bool_(True)],
                                 [1.0], [.5], [np.nan], [np.inf], ["1"], [[1]], None, "1"])
def test_invalid_legacy_sequences(bad):
    with pytest.raises((TypeError, ValueError), match="sequence"):
        legacy_slice_sequence_to_zero_based(bad, 2)


@pytest.mark.parametrize("count", [0, -1, True, 2., np.nan])
def test_invalid_slice_count(count):
    with pytest.raises((TypeError, ValueError), match="slice_type_count"):
        legacy_slice_sequence_to_zero_based([1], count)


@pytest.mark.parametrize("shape", [(4, 6), (3, 5), (1, 1)])
@pytest.mark.parametrize("ordering", ["centered", "unshifted"])
@pytest.mark.parametrize("sequence", [(0,), (0, 0), (0, 0, 0, 0)])
def test_unit_plane_default_explicit_and_repeated(shape, ordering, sequence):
    unit = unit_incident_wave(shape[1], shape[0])
    assert unit.dtype == np.complex128 and unit.shape == shape
    np.testing.assert_array_equal(unit, np.ones(shape))
    result = run([unit], [unit], sequence, ordering=ordering)
    explicit = run([unit], [unit], sequence, ordering=ordering, incident=unit)
    expected = np.zeros(shape, dtype=np.complex128)
    expected[(shape[0]//2, shape[1]//2) if ordering == "centered" else (0, 0)] = 1
    np.testing.assert_allclose(result.final_reciprocal_wave, expected, atol=3e-16, rtol=0)
    np.testing.assert_array_equal(result.final_reciprocal_wave, explicit.final_reciprocal_wave)
    np.testing.assert_allclose(result.amplitudes, [1], atol=3e-16, rtol=0)
    np.testing.assert_allclose(result.intensities, [1], atol=6e-16, rtol=0)
    assert result.slice_sequence_zero_based == sequence


@pytest.mark.parametrize("constant, expected", [(.6+.8j, .6+.8j), (.3+.4j, .3+.4j)])
def test_constant_phase_or_absorptive_grating(constant, expected):
    unit = np.ones((4, 6), dtype=np.complex128)
    result = run([np.full_like(unit, constant)], [unit], (0,), beams=((0, 0), (1, 0)))
    np.testing.assert_allclose(result.amplitudes, [expected, 0], atol=2e-16, rtol=0)
    np.testing.assert_allclose(result.intensities, [abs(expected)**2, 0], atol=3e-16, rtol=0)


@pytest.mark.parametrize("sequence, expected", [
    ((0, 1), [-.5+.5j, .5+.5j]),
    ((0, 1, 0, 1), [.5-.5j, -.5+.5j]),
    ((0, 1, 0, 1, 2, 2), [.125-.125j, -.125+.125j]),
])
def test_independent_two_slice_and_repeated_hand_calculations(sequence, expected):
    # Nx=2 transform: c=((u0+u1)/2,(u0-u1)/2), inverse=(c0+c1,c0-c1).
    # After A: c=(.5+.5i,.5+.5i), real=(1+i,0).
    # After B: c=(-.5+.5i,.5+.5i). Two more A/B yield (.5-.5i,-.5+.5i).
    # Region C has T=.5, P=1; two occurrences attenuate amplitudes by .25.
    gratings = [np.array([[1, 1j]]), np.array([[1j, 1]]), np.full((1, 2), .5+0j)]
    kernels = [np.array([[1, 1j]]), np.array([[1, -1j]]), np.ones((1, 2), dtype=np.complex128)]
    result = run(gratings, kernels, sequence, beams=((0, 0), (-1, 0)))
    np.testing.assert_array_equal(result.amplitudes, expected)
    np.testing.assert_array_equal(result.final_reciprocal_wave, [expected])
    np.testing.assert_array_equal(result.intensities, [.03125, .03125] if len(sequence) == 6 else [.5, .5])


def test_no_unused_final_return_transform_or_ms_plus_one(monkeypatch):
    import qcbed.multislice as module
    original = module.reciprocal_to_real_wave
    calls = []
    def traced(wave, *, ordering):
        calls.append(wave.copy())
        return original(wave, ordering=ordering)
    monkeypatch.setattr(module, "reciprocal_to_real_wave", traced)
    unit = np.ones((2, 3), dtype=np.complex128)
    for sequence in ((0,), (0, 0), (0, 0, 0, 0)):
        calls.clear()
        result = run([unit], [unit], sequence)
        assert len(calls) == len(sequence)-1
        assert result.slice_sequence_zero_based == sequence  # No sentinel appended.
        assert not hasattr(result, "real_wave")


BEAMS = ((0, 0), (2, 1), (-1, -1), (0, -2), (-3, 0), (1, -1), (-1, -1))
UNSHIFTED_POSITIONS = ((0, 0), (1, 2), (3, 5), (2, 0), (0, 3), (3, 1), (3, 5))
CENTERED_POSITIONS = ((2, 3), (3, 5), (1, 2), (0, 3), (2, 0), (1, 4), (1, 2))


@pytest.mark.parametrize("ordering, positions", [("centered", CENTERED_POSITIONS), ("unshifted", UNSHIFTED_POSITIONS)])
def test_rectangular_signed_beam_mapping_and_duplicates(ordering, positions):
    labels = (np.array([[15, 16, 17, 12, 13, 14], [21, 22, 23, 18, 19, 20],
                       [3, 4, 5, 0, 1, 2], [9, 10, 11, 6, 7, 8]]) if ordering == "centered"
              else np.arange(24).reshape(4, 6))
    wave = (labels+1j*(100+2*labels)).astype(np.complex128)
    assert beam_indices_to_row_col(BEAMS, mesh_x=6, mesh_y=4, ordering=ordering) == positions
    assert beam_indices_to_legacy_ib(BEAMS, mesh_x=6, mesh_y=4) == (1, 9, 24, 13, 4, 20, 24)
    expected = [100j, 8+116j, 23+146j, 12+124j, 3+106j, 19+138j, 23+146j]
    np.testing.assert_array_equal(extract_beam_amplitudes(wave, BEAMS, ordering=ordering), expected)
    # The IB helper is explicitly one-based; an independent C-order flatten checks it.
    flat = (np.arange(24)+1j*(100+2*np.arange(24))).astype(np.complex128)
    np.testing.assert_array_equal(flat[np.array([1, 9, 24, 13, 4, 20, 24])-1], expected)


def test_reference_input_beam_list_mapping_only():
    lines = (Path(__file__).parents[1]/"reference_data"/"FOR007.DAT").read_text().splitlines()
    # Fixture's actual READ record location; no parser/configuration extension.
    start = next(i for i, line in enumerate(lines) if line.strip() == "out.dat")+1
    count = int(lines[start])
    beams = tuple(tuple(map(int, line.split())) for line in lines[start+1:start+1+count])
    assert beams == ((0, 0), (1, 1), (1, -1), (-1, -1), (-1, 1), (2, 0), (0, -2), (-2, 0), (0, 2))
    assert beam_indices_to_legacy_ib(beams, mesh_x=16, mesh_y=16) == (1, 18, 242, 256, 32, 3, 225, 15, 33)


@pytest.mark.parametrize("h,k", [(3, 0), (-4, 0), (0, 2), (0, -3), (6, 0), (-6, 0)])
@pytest.mark.parametrize("ordering", ["centered", "unshifted"])
def test_outside_grid_and_aliases_are_rejected(h, k, ordering):
    with pytest.raises(ValueError, match="outside signed grid"):
        beam_indices_to_row_col([(h, k)], mesh_x=6, mesh_y=4, ordering=ordering)


def test_odd_grid_endpoints_and_empty_beams():
    assert beam_indices_to_row_col([(-2, -1), (2, 1)], mesh_x=5, mesh_y=3,
                                   ordering="centered") == ((0, 0), (2, 4))
    assert beam_indices_to_row_col([(-2, -1), (2, 1)], mesh_x=5, mesh_y=3,
                                   ordering="unshifted") == ((2, 3), (1, 2))
    unit = np.ones((3, 5), dtype=np.complex128)
    result = run([unit], [unit], (0,), beams=())
    assert result.amplitudes.shape == result.intensities.shape == (0,)
    assert result.amplitudes.dtype == np.complex128 and result.intensities.dtype == np.float64


def test_source_propagation_phase_sampled_and_no_rmt():
    grid = ReciprocalVectorGrid(ReciprocalLattice2D(DirectLattice2D(2, 4, 90)), ReciprocalIndexGrid(4, 4))
    kernel = propagation_kernel(grid, PropagationParameters(.025, 2), ordering="unshifted")
    # Explicit grating creates the H=+1 mode from unit illumination.
    grating = np.array([[1, -1j, -1, 1j]]*4)
    result = run([grating], [kernel], (0,), beams=((1, 0), (0, 0), (1, 0)))
    expected = .9992290362424846+.03925981571423079j
    np.testing.assert_allclose(result.amplitudes, [expected, 0, expected], atol=3e-16, rtol=0)
    np.testing.assert_allclose(result.intensities, [1, 0, 1], atol=3e-16, rtol=0)


def nested_transform(values, *, to_reciprocal):
    ny, nx = values.shape
    out = np.zeros_like(values)
    sign = 1 if to_reciprocal else -1
    scale = nx*ny if to_reciprocal else 1
    for r in range(ny):
        for c in range(nx):
            out[r, c] = sum(values[j, i]*cmath.exp(sign*2j*math.pi*(i*c/nx+j*r/ny))
                            for j in range(ny) for i in range(nx))/scale
    return out


def numerical_case(shape):
    ny, nx = shape
    incident = np.array([[complex((c+2*r+1)/9, (2*c-r)/7) for c in range(nx)] for r in range(ny)])
    gratings = [np.array([[cmath.exp(1j*(c+2*r+s)/5) for c in range(nx)] for r in range(ny)]) for s in (0, 2)]
    kernels = [np.array([[cmath.exp(1j*(c*c+r+s*c)/11) for c in range(nx)] for r in range(ny)]) for s in (1, 3)]
    return incident, gratings, kernels


def nested_sequence(incident, gratings, kernels, sequence):
    wave = incident.copy()
    for occurrence, index in enumerate(sequence):
        reciprocal = nested_transform(wave*gratings[index], to_reciprocal=True)*kernels[index]
        if occurrence != len(sequence)-1:
            wave = nested_transform(reciprocal, to_reciprocal=False)
    return reciprocal


@pytest.mark.parametrize("shape", [(2, 3), (3, 5), (4, 6)])
@pytest.mark.parametrize("ordering", ["centered", "unshifted"])
def test_repeated_sequence_against_nested_dft_and_norm(shape, ordering):
    incident, gratings, kernels = numerical_case(shape)
    sequence = (0, 1, 0, 1)
    expected = nested_sequence(incident, gratings, kernels, sequence)
    supplied = kernels if ordering == "unshifted" else [np.fft.fftshift(k) for k in kernels]
    result = run(gratings, supplied, sequence, ordering=ordering, incident=incident)
    # Independent NumPy permutation only; expected numerical transforms above are nested sums.
    actual = result.final_reciprocal_wave if ordering == "unshifted" else np.fft.ifftshift(result.final_reciprocal_wave)
    np.testing.assert_allclose(actual, expected, atol=8e-15, rtol=0)
    initial_norm = np.sum(abs(incident)**2)/incident.size
    assert abs(np.sum(abs(actual)**2)-initial_norm) < 4e-15


@pytest.mark.parametrize("ordering", ["centered", "unshifted"])
def test_absorption_and_boundary_loss_are_distinct(ordering):
    shape = (4, 6)
    unit = np.ones(shape, dtype=np.complex128)
    absorptive = phase_grating_from_potential(np.full_like(unit, 2-.5j), voltage_kv=200)
    absorbed = run([absorptive], [unit], (0, 0), ordering=ordering)
    single_amplitude = .9996345863381000-.001457140527998951j
    np.testing.assert_allclose(absorbed.amplitudes, [single_amplitude**2], atol=6e-16, rtol=0)
    assert np.sum(abs(absorbed.final_reciprocal_wave)**2) == pytest.approx(abs(single_amplitude)**4, abs=8e-16)
    assert absorbed.intensities[0] < 1
    grid = ReciprocalVectorGrid(ReciprocalLattice2D(DirectLattice2D(2, 4, 90)),
                               ReciprocalIndexGrid(6, 4, centered=ordering == "centered"))
    kernel = propagation_kernel(grid, PropagationParameters(.025, 0), ordering=ordering)
    # No absorption; initial wave comprises equal DC and x Nyquist amplitudes.
    incident = np.array([[2, 0, 2, 0, 2, 0]]*4, dtype=np.complex128)
    result = run([unit], [kernel], (0, 0), beams=((0, 0), (-3, 0)), ordering=ordering, incident=incident)
    np.testing.assert_allclose(result.amplitudes, [1, 0], atol=3e-16, rtol=0)
    assert np.sum(abs(incident)**2)/incident.size == 2
    assert np.sum(abs(result.final_reciprocal_wave)**2) == pytest.approx(1, abs=6e-16)


def test_result_immutability_and_no_input_mutation():
    incident, gratings, kernels = numerical_case((4, 6))
    incident = incident[:, ::-1]
    gratings = [g[:, ::-1] for g in gratings]
    inputs = [incident, *gratings, *kernels]
    copies = [v.copy() for v in inputs]
    for value in inputs:
        value.setflags(write=False)
    sequence, beams = [0, 1, 0, 1], [[0, 0], [1, 0]]
    result = run(gratings, kernels, sequence, beams=beams, incident=incident)
    for value, original in zip(inputs, copies):
        np.testing.assert_array_equal(value, original)
    sequence[0], beams[0][0] = 1, 1
    assert result.slice_sequence_zero_based == (0, 1, 0, 1)
    assert result.beam_indices_hk == ((0, 0), (1, 0))
    for value, dtype in ((result.final_reciprocal_wave, np.complex128),
                         (result.amplitudes, np.complex128), (result.intensities, np.float64)):
        assert value.dtype == dtype and value.flags.c_contiguous
        assert not any(np.shares_memory(value, source) for source in inputs)
        with pytest.raises(ValueError):
            value.setflags(write=True)
        with pytest.raises(ValueError):
            value.flat[0] = 0
    with pytest.raises(FrozenInstanceError):
        result.reciprocal_ordering = "centered"


@pytest.mark.parametrize("sequence", [(), (-1,), (2,), (True,), (1.,), (np.nan,), ("0",), None])
def test_invalid_zero_based_sequence(sequence):
    unit = np.ones((2, 3), dtype=np.complex128)
    with pytest.raises((TypeError, ValueError), match="slice_sequence"):
        run([unit, unit], [unit, unit], sequence)


@pytest.mark.parametrize("bad", [np.zeros((2, 3)), np.zeros((2, 3), dtype=np.complex64),
                                 np.zeros((1, 3), dtype=np.complex128), np.zeros(6, dtype=np.complex128),
                                 np.zeros((0, 3), dtype=np.complex128),
                                 np.full((2, 3), np.nan, dtype=np.complex128),
                                 np.full((2, 3), np.inf, dtype=np.complex128)])
@pytest.mark.parametrize("slot", ["grating", "propagator", "incident"])
def test_invalid_arrays_including_unused_stored_type(bad, slot):
    unit = np.ones((2, 3), dtype=np.complex128)
    gratings = [unit, bad if slot == "grating" else unit]
    kernels = [unit, bad if slot == "propagator" else unit]
    with pytest.raises((TypeError, ValueError), match="complex128|shape|finite"):
        run(gratings, kernels, (0,), incident=bad if slot == "incident" else None)


def test_collections_and_ordering_validation():
    unit = np.ones((2, 3), dtype=np.complex128)
    for gratings, kernels in (([], []), ([unit], []), ([unit, unit], [unit])):
        with pytest.raises(ValueError, match="collection lengths"):
            run(gratings, kernels, (0,))
    for bad in (None, "centered", 1):
        with pytest.raises(TypeError, match="phase_gratings"):
            run(bad, [unit], (0,))
    with pytest.raises(ValueError, match="ordering.*match"):
        simulate_slice_sequence([unit], [unit], [0], beam_indices_hk=[(0, 0)],
                                reciprocal_ordering="centered", propagator_ordering="unshifted")
    with pytest.raises(ValueError, match="ordering"):
        run([unit], [unit], [0], ordering=None)
    with pytest.raises(TypeError, match="ordering"):
        simulate_slice_sequence([unit], [unit], [0], beam_indices_hk=[])


@pytest.mark.parametrize("bad", [[(1., 0)], [(True, 0)], [(0, np.nan)], [(0, np.inf)],
                                 [("1", 0)], [(1,)], [(1, 2, 3)], [1], None, "0,0"])
def test_invalid_beam_lists(bad):
    with pytest.raises((ValueError, TypeError), match="beam"):
        beam_indices_to_row_col(bad, mesh_x=6, mesh_y=4, ordering="unshifted")


def test_beam_intensity_values_and_failures():
    np.testing.assert_array_equal(beam_intensities(np.array([3+4j, -2+1j, 0j])), [25., 5., 0.])
    for bad in ([1j], np.ones(2), np.ones(2, dtype=np.complex64),
                np.ones((1, 2), dtype=np.complex128), np.array([complex(np.nan, 0)]),
                np.array([complex(0, np.inf)]), np.array([1e308+0j])):
        with pytest.raises((TypeError, ValueError), match="complex128|shape|finite"):
            beam_intensities(bad)
    unit = np.ones((2, 3), dtype=np.complex128)
    with pytest.raises(ValueError, match="finite"):
        run([np.full_like(unit, 1e308)], [unit], (0,))
    with pytest.raises(ValueError, match="positive"):
        unit_incident_wave(0, 3)
    with pytest.raises(TypeError, match="integer"):
        unit_incident_wave(True, 3)


def test_result_defensively_copies_external_wave():
    wave = np.ones((2, 3), dtype=np.complex128)
    result = MultisliceResult(wave, ((0, 0),), (0,), 1, "unshifted")
    wave[:] = 0
    np.testing.assert_array_equal(result.final_reciprocal_wave, 1)
    np.testing.assert_array_equal(result.amplitudes, [1])
