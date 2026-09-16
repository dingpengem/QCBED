import cmath
from dataclasses import FrozenInstanceError
import math

import numpy as np
import pytest

from qcbed.wave_propagation import (
    WaveStepResult, apply_phase_grating, real_to_reciprocal_wave,
    apply_propagator, reciprocal_to_real_wave, propagate_one_step,
)
from qcbed.phase_grating import phase_grating_from_potential
from qcbed.propagator import PropagationParameters, propagation_kernel
from qcbed.reciprocal import (DirectLattice2D, ReciprocalLattice2D,
                             ReciprocalIndexGrid, ReciprocalVectorGrid)


def step(wave, grating, propagator, ordering="unshifted"):
    return propagate_one_step(wave, grating, propagator, ordering=ordering,
                              propagator_ordering=ordering)


def bin_at(h, k, shape, ordering):
    ny, nx = shape
    return ((k+ny//2) % ny, (h+nx//2) % nx) if ordering == "centered" else (k % ny, h % nx)


@pytest.mark.parametrize("shape", [(4, 6), (3, 5), (1, 1)])
@pytest.mark.parametrize("ordering", ["centered", "unshifted"])
@pytest.mark.parametrize("amplitude", [1+0j, .6+.8j, .3+.4j])
def test_uniform_phase_and_absorption(shape, ordering, amplitude):
    wave = np.ones(shape, dtype=np.complex128)
    result = step(wave, np.full(shape, amplitude), wave, ordering)
    expected_c = np.zeros(shape, dtype=np.complex128)
    expected_c[bin_at(0, 0, shape, ordering)] = amplitude
    np.testing.assert_allclose(result.real_wave, amplitude, atol=4e-16, rtol=0)
    np.testing.assert_allclose(result.reciprocal_wave, expected_c, atol=2e-16, rtol=0)
    assert np.sum(abs(result.real_wave)**2) == pytest.approx(wave.size*abs(amplitude)**2, abs=1e-14)


@pytest.mark.parametrize("h,k", [(1, 0), (-1, 0), (0, 1), (0, -1), (1, -1)])
@pytest.mark.parametrize("ordering", ["centered", "unshifted"])
def test_signed_modes_rectangular_axes_and_sampling_after_propagator(h, k, ordering):
    # x has 4 samples, y has 8: explicit fourth/eighth roots, no FFT oracle.
    roots_x = [1, -1j, -1, 1j]
    q = math.sqrt(.5)
    roots_y = [1, q-1j*q, -1j, -q-1j*q, -1, -q+1j*q, 1j, q+1j*q]
    wave = np.array([[roots_x[(h*c) % 4]*roots_y[(k*r) % 8]
                      for c in range(4)] for r in range(8)], dtype=np.complex128)
    expected = np.zeros((8, 4), dtype=np.complex128)
    target = bin_at(h, k, wave.shape, ordering)
    expected[target] = 1
    np.testing.assert_allclose(real_to_reciprocal_wave(wave, ordering=ordering), expected, atol=2e-16)
    np.testing.assert_allclose(reciprocal_to_real_wave(expected, ordering=ordering), wave, atol=2e-16)
    kernel = np.ones_like(wave)
    kernel[target] = .6+.8j
    multiplied = apply_propagator(expected, kernel, wave_ordering=ordering, propagator_ordering=ordering)
    np.testing.assert_array_equal(multiplied, expected*(.6+.8j))
    result = step(wave, np.ones_like(wave), kernel, ordering)
    np.testing.assert_allclose(result.reciprocal_wave, expected*(.6+.8j), atol=2e-16)
    np.testing.assert_allclose(result.real_wave, wave*(.6+.8j), atol=4e-16)


def nested_step(wave, grating, kernel, ordering):
    """Independent scalar DFT in each direction, including complete intermediate."""
    ny, nx = wave.shape
    sampled = np.zeros_like(wave)
    out = np.zeros_like(wave)
    hs = list(range(-(nx//2), nx-nx//2)) if ordering == "centered" else list(range(nx))
    ks = list(range(-(ny//2), ny-ny//2)) if ordering == "centered" else list(range(ny))
    for j, k in enumerate(ks):
        for i, h in enumerate(hs):
            total = 0j
            for r in range(ny):
                for c in range(nx):
                    total += wave[r, c]*grating[r, c]*cmath.exp(2j*math.pi*(h*c/nx+k*r/ny))
            sampled[j, i] = total/(nx*ny)*kernel[j, i]
    for r in range(ny):
        for c in range(nx):
            for j, k in enumerate(ks):
                for i, h in enumerate(hs):
                    out[r, c] += sampled[j, i]*cmath.exp(-2j*math.pi*(h*c/nx+k*r/ny))
    return out, sampled


def fixed_case(shape, ordering):
    ny, nx = shape
    wave = np.array([[complex((2*c-r+1)/7, (c+3*r-2)/11) for c in range(nx)] for r in range(ny)])
    # Spatially varying multipliers make moving or shifting the grating observable.
    grating = np.array([[cmath.exp(1j*(r+2*c)/5) for c in range(nx)] for r in range(ny)])
    kernel = np.empty(shape, dtype=np.complex128)
    for k in range(-(ny//2), ny-ny//2):
        for h in range(-(nx//2), nx-nx//2):
            kernel[bin_at(h, k, shape, ordering)] = cmath.exp(1j*(h*h+2*k*k+h-k)/9)
    return wave, grating, kernel


@pytest.mark.parametrize("shape", [(4, 6), (3, 5), (4, 4)])
@pytest.mark.parametrize("ordering", ["centered", "unshifted"])
def test_complete_nested_dft_and_norm(shape, ordering):
    wave, grating, kernel = fixed_case(shape, ordering)
    expected_real, expected_c = nested_step(wave, grating, kernel, ordering)
    result = step(wave, grating, kernel, ordering)
    np.testing.assert_allclose(result.reciprocal_wave, expected_c, atol=8e-15, rtol=0)
    np.testing.assert_allclose(result.real_wave, expected_real, atol=2e-14, rtol=0)
    norm = np.sum(abs(wave)**2)
    assert abs(np.sum(abs(result.real_wave)**2)-norm) < 3e-14
    assert abs(wave.size*np.sum(abs(result.reciprocal_wave)**2)-norm) < 3e-14
    opposite = "centered" if ordering == "unshifted" else "unshifted"
    other = step(*fixed_case(shape, opposite), opposite)
    np.testing.assert_allclose(other.real_wave, result.real_wave, atol=5e-16, rtol=0)


@pytest.mark.parametrize("ordering", ["centered", "unshifted"])
def test_existing_source_grating_absorbs_once(ordering):
    wave = np.ones((4, 6), dtype=np.complex128)
    grating = phase_grating_from_potential(np.full_like(wave, 2-.5j), voltage_kv=200)
    result = step(wave, grating, wave, ordering)
    # Independently checked printed-source PGMAP case from its numerical audit.
    expected = .9996345863381000-.001457140527998951j
    np.testing.assert_allclose(result.real_wave, expected, atol=4e-16, rtol=0)
    assert np.sum(abs(result.real_wave)**2) < wave.size
    assert np.sum(abs(result.real_wave)**2) == pytest.approx(wave.size*abs(expected)**2, abs=1e-14)


@pytest.mark.parametrize("ordering", ["centered", "unshifted"])
@pytest.mark.parametrize("shape,boundary", [((4, 6), "legacy"), ((3, 5), "none")])
def test_existing_propagator_boundaries_and_odd_support(shape, boundary, ordering):
    ny, nx = shape
    grid = ReciprocalVectorGrid(ReciprocalLattice2D(DirectLattice2D(2, 4, 90)),
                               ReciprocalIndexGrid(nx, ny, centered=ordering == "centered"))
    kernel = propagation_kernel(grid, PropagationParameters(.025, 0), ordering=ordering, boundary=boundary)
    uniform = np.ones(shape, dtype=np.complex128)
    np.testing.assert_allclose(step(uniform, uniform, kernel, ordering).real_wave, uniform, atol=3e-16)
    if boundary == "legacy":
        # Unit-modulus real wave entirely at H=-Nx/2. No absorptive grating.
        nyquist = np.array([[(-1.)**c for c in range(nx)] for _ in range(ny)], dtype=np.complex128)
        result = step(nyquist, uniform, kernel, ordering)
        np.testing.assert_allclose(result.real_wave, 0, atol=2e-16)
        assert np.all(result.reciprocal_wave[:, nx//2 if ordering == "unshifted" else 0] == 0)
        assert np.sum(abs(result.real_wave)**2) < 1e-28
    else:
        wave, grating, _ = fixed_case(shape, ordering)
        np.testing.assert_allclose(step(wave, grating, kernel, ordering).real_wave, wave*grating, atol=8e-16)


def test_nonzero_source_propagator_on_one_beam():
    grid = ReciprocalVectorGrid(ReciprocalLattice2D(DirectLattice2D(2, 4, 90)), ReciprocalIndexGrid(4, 4))
    kernel = propagation_kernel(grid, PropagationParameters(.025, 2), ordering="unshifted")
    wave = np.array([[1, -1j, -1, 1j]]*4, dtype=np.complex128)
    result = step(wave, np.ones_like(wave), kernel)
    expected = .9992290362424846+.03925981571423079j
    assert result.reciprocal_wave[0, 1] == pytest.approx(expected, abs=3e-16, rel=0)
    np.testing.assert_allclose(result.real_wave, expected*wave, atol=3e-16, rtol=0)


def test_immutable_snapshots_and_noncontiguous_inputs():
    inputs = tuple(v[:, ::-1] for v in fixed_case((4, 6), "unshifted"))
    originals = tuple(v.copy() for v in inputs)
    for v in inputs:
        v.setflags(write=False)
    result = step(*inputs)
    for original, v in zip(originals, inputs):
        np.testing.assert_array_equal(original, v)
    for out in (result.real_wave, result.reciprocal_wave):
        assert out.shape == (4, 6) and out.dtype == np.complex128 and out.flags.c_contiguous
        assert not any(np.shares_memory(out, v) for v in inputs)
        with pytest.raises(ValueError):
            out[0, 0] = 0
        with pytest.raises(ValueError):
            out.setflags(write=True)
    with pytest.raises(FrozenInstanceError):
        result.ordering = "centered"
    source = np.ones((2, 3), dtype=np.complex128)
    snapshot = WaveStepResult(source, source, "unshifted")
    source[:] = 7
    np.testing.assert_array_equal(snapshot.real_wave, 1)


BAD_ARRAYS = [np.ones((2, 3)), np.ones((2, 3), dtype=np.complex64), [[1j]],
              np.ones(6, dtype=np.complex128), np.ones((0, 3), dtype=np.complex128),
              np.full((2, 3), np.nan, dtype=np.complex128),
              np.full((2, 3), np.inf, dtype=np.complex128),
              np.full((2, 3), -np.inf, dtype=np.complex128)]


@pytest.mark.parametrize("bad", BAD_ARRAYS)
@pytest.mark.parametrize("slot", [0, 1, 2])
def test_invalid_step_arrays(bad, slot):
    inputs = [np.ones((2, 3), dtype=np.complex128) for _ in range(3)]
    inputs[slot] = bad
    with pytest.raises((TypeError, ValueError), match="complex128|shape|finite"):
        step(*inputs)


def test_shape_ordering_and_overflow_errors():
    good = np.ones((2, 3), dtype=np.complex128)
    wrong = np.ones((1, 3), dtype=np.complex128)
    for inputs in ((good, wrong, good), (good, good, wrong)):
        with pytest.raises(ValueError, match="shape"):
            step(*inputs)
    for bad in (None, True, [], "Centered"):
        with pytest.raises(ValueError, match="ordering"):
            step(good, good, good, bad)
    with pytest.raises(ValueError, match="ordering.*match"):
        propagate_one_step(good, good, good, ordering="centered", propagator_ordering="unshifted")
    with pytest.raises(ValueError, match="ordering.*match"):
        apply_propagator(good, good, wave_ordering="unshifted", propagator_ordering="centered")
    with pytest.raises(TypeError, match="ordering"):
        propagate_one_step(good, good, good)
    huge = np.full_like(good, 1e308)
    for operation in (lambda: apply_phase_grating(huge, huge),
                      lambda: apply_propagator(huge, huge, wave_ordering="unshifted", propagator_ordering="unshifted"),
                      lambda: reciprocal_to_real_wave(huge, ordering="unshifted")):
        with pytest.raises(ValueError, match="finite"):
            operation()


@pytest.mark.parametrize("bad", BAD_ARRAYS)
def test_standalone_transforms_validate(bad):
    for operation in (real_to_reciprocal_wave, reciprocal_to_real_wave):
        with pytest.raises((TypeError, ValueError), match="complex128|shape|finite"):
            operation(bad, ordering="unshifted")


def test_result_validation():
    good = np.ones((2, 3), dtype=np.complex128)
    with pytest.raises(ValueError, match="shape"):
        WaveStepResult(good, good.T, "unshifted")
    with pytest.raises(ValueError, match="ordering"):
        WaveStepResult(good, good, "unknown")
