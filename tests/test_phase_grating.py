import cmath
from dataclasses import replace
import math
from pathlib import Path

import numpy as np
import pytest

from qcbed import (AtomSite, Species, SliceGeometry, DirectLattice2D,
    ReciprocalLattice2D, ReciprocalIndexGrid, ReciprocalVectorGrid,
    AbsorptionCoefficients, elastic_structure_factors, absorptive_structure_factors)
from qcbed.phase_grating import (
    reciprocal_cutoff_mask, legacy_nyquist_mask, projected_potential_direct,
    electron_interaction_constant, phase_grating_from_potential,
    phase_grating_from_structure_factors,
)


def grid(mesh=(4, 4), centered=True, a=20, b=10, gamma=90):
    return ReciprocalVectorGrid(ReciprocalLattice2D(DirectLattice2D(a, b, gamma)),
                                ReciprocalIndexGrid(*mesh, centered=centered))


def coefficients(g, entries):
    result = np.zeros(g.indices.h.shape, dtype=np.complex128)
    for h, k, value in entries:
        result[(g.indices.h == h) & (g.indices.k == k)] = value
    return result


PHASE_CASES = [
    (100, .0009243913556182865, .9995362029095982-.0018479273565984993j),
    (200, .0007288360746343612, .9996345863381000-.001457140527998951j),
    (300, .0006526122090303823, .9996728956001657-.0013047982143670319j),
]


@pytest.mark.parametrize("kv,sigma,expected", PHASE_CASES)
def test_independent_interaction_constant_and_pgmap(kv, sigma, expected):
    # Hard-coded independent 65-digit Decimal evaluation of the printed
    # wavelength/SIGMA equations; Taylor cos/sin, Decimal exp for V=2-.5i.
    assert electron_interaction_constant(kv) == pytest.approx(sigma, rel=4e-16)
    result = phase_grating_from_potential(np.full((2, 4), 2-.5j, dtype=np.complex128), voltage_kv=kv)
    np.testing.assert_allclose(result, expected, rtol=0, atol=3e-16)


@pytest.mark.parametrize("mesh", [(4, 4), (6, 4)])
@pytest.mark.parametrize("centered", [True, False])
def test_dc_uniform_and_no_transform_normalization(mesh, centered):
    g = grid(mesh, centered)
    c = coefficients(g, [(0, 0, 2-.5j)])
    potential = projected_potential_direct(c, g)
    np.testing.assert_array_equal(potential, np.full((mesh[1], mesh[0]), 2-.5j))
    result = phase_grating_from_structure_factors(c.real.astype(np.complex128),
        c.imag.astype(np.complex128), g, voltage_kv=200, dsap_inverse_angstrom=0)
    np.testing.assert_allclose(result, .9996345863381000-.001457140527998951j, rtol=0, atol=3e-16)


@pytest.mark.parametrize("sign", [-1, 1])
def test_single_h_coefficient_negative_fourier_sign(sign):
    g = grid((4, 2))
    c = coefficients(g, [(sign, 0, 1)])
    expected = [[1, -1j*sign, -1, 1j*sign]] * 2
    np.testing.assert_allclose(projected_potential_direct(c, g), expected, atol=4e-16, rtol=0)


def test_single_k_coefficient_no_transpose_and_no_endpoint_duplicate():
    g = grid((6, 4))
    c = coefficients(g, [(0, 1, 1)])
    expected = np.repeat(np.array([1, -1j, -1, 1j])[:, None], 6, axis=1)
    np.testing.assert_allclose(projected_potential_direct(c, g), expected, atol=4e-16, rtol=0)
    # Last sample is y=3/4, not y=1: the repeated endpoint is omitted.
    assert expected[-1, 0] != expected[0, 0]


def test_friedel_pair_real_potential_and_unit_elastic_modulus():
    g = grid()
    c = coefficients(g, [(1, 0, 1+2j), (-1, 0, 1-2j)])
    potential = projected_potential_direct(c, g)
    # (1+2i)e^-it+(1-2i)e^it = 2*cos(t)+4*sin(t).
    np.testing.assert_allclose(potential, [[2, 4, -2, -4]]*4, atol=2e-15, rtol=0)
    assert np.max(abs(potential.imag)) < 1e-15
    np.testing.assert_allclose(abs(phase_grating_from_potential(potential, voltage_kv=200)), 1, atol=3e-16, rtol=0)


def test_absorption_direction_without_sign_clipping():
    v = np.array([[-100j, 0j, 100j]], dtype=np.complex128)
    result = phase_grating_from_potential(v, voltage_kv=200)
    assert abs(result[0, 0]) < 1
    assert result[0, 1] == 1+0j
    assert abs(result[0, 2]) > 1  # Preserve equation, do not force attenuation.
    assert np.all(result.imag == 0)


def nested_reference(entries, nx, ny):
    # Explicit mathematical sum, independent of grid helpers or production
    # array broadcasting. Signed H/K entries make the coordinate meaning clear.
    result = np.zeros((ny, nx), dtype=np.complex128)
    for y in range(ny):
        for x in range(nx):
            for h, k, value in entries:
                result[y, x] += value*cmath.exp(-2j*math.pi*(h*x/nx + k*y/ny))
    return result


@pytest.mark.parametrize("mesh", [(4, 4), (6, 4), (5, 3)])
@pytest.mark.parametrize("centered", [True, False])
def test_explicit_nested_reference_every_pixel(mesh, centered):
    entries = [(0, 0, 2-.5j), (1, -1, 1+2j), (-1, 1, -.25+3j), (0, 1, .75)]
    g = grid(mesh, centered)
    actual = projected_potential_direct(coefficients(g, entries), g)
    np.testing.assert_allclose(actual, nested_reference(entries, *mesh), rtol=0, atol=8e-15)


def test_fourier_sum_on_tiny_mesh_explicit_roots():
    g = grid((4, 2), centered=False)
    c = coefficients(g, [(0, 0, 1), (1, 0, 2), (0, -1, 3j)])
    expected = [[3+3j, 1+1j, -1+3j, 1+5j],
                [3-3j, 1-5j, -1-3j, 1-1j]]
    np.testing.assert_allclose(projected_potential_direct(c, g), expected, rtol=0, atol=1e-15)


def test_periodic_shift_is_cyclic_not_a_transpose():
    g = grid((4, 2))
    # Translation by one x pixel multiplies H=1 by exp(+2*pi*i/4)=i.
    c = coefficients(g, [(1, 0, 1)])
    shifted = coefficients(g, [(1, 0, 1j)])
    expected = np.roll(projected_potential_direct(c, g), 1, axis=1)
    np.testing.assert_allclose(projected_potential_direct(shifted, g), expected, atol=4e-16, rtol=0)


def test_cutoff_boundary_and_zero_origin():
    g = grid((6, 4))
    below = reciprocal_cutoff_mask(g, np.nextafter(.025, 0))
    equal = reciprocal_cutoff_mask(g, .025)
    above = reciprocal_cutoff_mask(g, np.nextafter(.025, 1))
    assert not below[2, 4]
    assert equal[2, 4] and above[2, 4]
    assert equal[2, 2]  # H=-1 is retained too.
    zero = reciprocal_cutoff_mask(g, 0)
    assert zero.sum() == 1 and zero[2, 3]
    assert equal.dtype == np.bool_


@pytest.mark.parametrize("centered", [True, False])
def test_nyquist_mask_is_reciprocal_not_spatial(centered):
    g = grid((6, 4), centered)
    mask = legacy_nyquist_mask(g)
    expected = np.ones((4, 6), dtype=bool)
    expected[0 if centered else 2, :] = False
    expected[:, 0 if centered else 3] = False
    np.testing.assert_array_equal(mask, expected)
    # Large boundary coefficients must not survive the composed operation.
    c = coefficients(g, [(-3, 0, 100), (0, -2, 100), (0, 0, 2-.5j)])
    result = phase_grating_from_structure_factors(c, np.zeros_like(c), g,
        voltage_kv=200, dsap_inverse_angstrom=100)
    np.testing.assert_allclose(result, .9996345863381000-.001457140527998951j, atol=3e-16, rtol=0)
    assert np.all(abs(result[0, :]) > 0) and np.all(abs(result[:, 0]) > 0)


def test_oblique_geometry_affects_cutoff_not_fractional_fourier_kernel():
    orth = grid(a=2, b=2)
    oblique = grid(a=2, b=2, gamma=60)
    assert not reciprocal_cutoff_mask(orth, .3)[3, 3]  # STL=sqrt(2)/4
    assert reciprocal_cutoff_mask(oblique, .3)[3, 3]   # STL=1/sqrt(12)
    c = coefficients(orth, [(1, 1, 2j)])
    np.testing.assert_array_equal(projected_potential_direct(c, orth), projected_potential_direct(c, oblique))


def test_existing_atomic_apis_thickness_and_slice_shifts():
    atom = (AtomSite(1, 0, 0, 1, "S"),)
    species = (Species("AL", 13),)
    s = SliceGeometry(2, cell_a_angstrom=20, cell_b_angstrom=10, lattice_gamma_degrees=90)
    outputs = []
    for geometry in (s, replace(s, thickness_angstrom=7, x_shift_fractional=.25, y_shift_fractional=-.3)):
        g = grid(a=geometry.cell_a_angstrom, b=geometry.cell_b_angstrom, gamma=geometry.lattice_gamma_degrees)
        e = elastic_structure_factors(species, atom, g,
            applied_debye_waller_angstrom_squared=(.86,), vacancy_fraction=0)
        a = absorptive_structure_factors(species, atom, g, voltage_kv=200,
            applied_debye_waller_angstrom_squared=(.86,), coefficients=(AbsorptionCoefficients(0, 1, 1),),
            vacancy_fraction=0, absorption_raw="YE")
        outputs.append(phase_grating_from_structure_factors(e, a, g, voltage_kv=200, dsap_inverse_angstrom=0))
    # Independent Al origin coefficient 1.409767975005-9.292395257301284i.
    np.testing.assert_allclose(outputs[0], .9932497254026199-.0010205542781868584j, atol=3e-16, rtol=0)
    np.testing.assert_array_equal(*outputs)


def test_zero_potential_exact_identity_and_input_immutability():
    g = grid((6, 4))
    c = np.zeros((4, 6), dtype=np.complex128)
    before = g.stl_inverse_angstrom.copy()
    c.setflags(write=False)
    v = projected_potential_direct(c, g)
    np.testing.assert_array_equal(v, c)
    t = phase_grating_from_structure_factors(c, c, g, voltage_kv=200, dsap_inverse_angstrom=1)
    np.testing.assert_array_equal(t, np.ones((4, 6), dtype=np.complex128))
    assert t.dtype == v.dtype == np.complex128
    assert t.flags.c_contiguous and not np.shares_memory(t, c)
    assert not np.any(c) and not c.flags.writeable
    np.testing.assert_array_equal(g.stl_inverse_angstrom, before)


@pytest.mark.parametrize("bad", [0, -1, float("nan"), float("inf"), -float("inf"), True, "200"])
def test_invalid_voltage(bad):
    with pytest.raises((ValueError, TypeError), match="voltage_kv"):
        phase_grating_from_potential(np.zeros((2, 2), dtype=np.complex128), voltage_kv=bad)


@pytest.mark.parametrize("bad", [0, -1, float("nan"), float("inf")])
def test_invalid_derived_wavelength_is_rejected(monkeypatch, bad):
    # There is no independent wavelength override; guard the derived quantity.
    monkeypatch.setattr("qcbed.phase_grating.electron_wavelength_angstrom", lambda _: bad)
    with pytest.raises(ValueError, match="wavelength_angstrom"):
        electron_interaction_constant(200)


@pytest.mark.parametrize("bad", [0, -1, float("nan"), float("inf")])
def test_invalid_thickness_stays_at_geometry_boundary(bad):
    # PGMAP has no thickness argument; existing slice geometry validates it.
    with pytest.raises(ValueError, match="thickness_angstrom"):
        SliceGeometry(bad)


@pytest.mark.parametrize("bad", [-1, float("nan"), float("inf"), True, "1"])
def test_invalid_cutoff(bad):
    with pytest.raises((TypeError, ValueError), match="dsap_inverse_angstrom"):
        reciprocal_cutoff_mask(grid(), bad)


@pytest.mark.parametrize("bad", [np.zeros((4, 4)), np.zeros((4, 4), dtype=np.complex64), [[1j]]])
def test_wrong_dtypes(bad):
    with pytest.raises(TypeError, match="complex128"):
        projected_potential_direct(bad, grid())
    with pytest.raises(TypeError, match="complex128"):
        phase_grating_from_potential(bad, voltage_kv=200)


@pytest.mark.parametrize("bad", [np.zeros((3, 4), dtype=np.complex128),
    np.zeros(16, dtype=np.complex128), np.zeros((0, 4), dtype=np.complex128),
    np.full((4, 4), np.nan, dtype=np.complex128), np.full((4, 4), np.inf, dtype=np.complex128)])
def test_invalid_shapes_and_nonfinite_coefficients(bad):
    with pytest.raises(ValueError, match="shape|finite"):
        projected_potential_direct(bad, grid())


def test_invalid_grid_odd_legacy_mask_and_overflow():
    with pytest.raises(TypeError, match="reciprocal_grid"):
        reciprocal_cutoff_mask(None, 1)
    with pytest.raises(ValueError, match="even"):
        legacy_nyquist_mask(grid((5, 3)))
    with pytest.raises(ValueError, match="shape"):
        phase_grating_from_structure_factors(np.zeros((2, 2), dtype=np.complex128),
            np.zeros((2, 2), dtype=np.complex128), grid(), voltage_kv=200, dsap_inverse_angstrom=1)
    with pytest.raises(ValueError, match="finite"):
        projected_potential_direct(np.full((4, 4), 1e308, dtype=np.complex128), grid())
    with pytest.raises(ValueError, match="finite"):
        phase_grating_from_potential(np.full((2, 2), 1e308j, dtype=np.complex128), voltage_kv=200)
    with pytest.raises(ValueError, match="wavelength"):
        electron_interaction_constant(1e308)


def test_saved_interaction_constant_only():
    # This is a scalar saved-output comparison, not a saved PGARRAY fixture.
    lines = (Path(__file__).resolve().parents[1]/"reference_data"/"FOR0013.DAT").read_text().splitlines()
    line = next(line for line in lines if "INTERACTION CONSTANT (SIGMA)" in line)
    recorded = float(line.split()[-1])
    assert electron_interaction_constant(200) == pytest.approx(recorded, abs=5e-9, rel=0)


@pytest.mark.parametrize("mesh", [(4, 4), (6, 4)])
@pytest.mark.parametrize("centered", [True, False])
@pytest.mark.parametrize("cutoff", [0, .025, 1])
def test_fft_complete_phase_grating_matches_direct(mesh, centered, cutoff):
    g = grid(mesh, centered)
    rng = np.random.default_rng(781)
    e = (rng.normal(size=g.indices.h.shape)+1j*rng.normal(size=g.indices.h.shape)).astype(np.complex128)
    a = (rng.normal(size=e.shape)+1j*rng.normal(size=e.shape)).astype(np.complex128)
    original_e, original_a = e.copy(), a.copy()
    options = dict(voltage_kv=200, dsap_inverse_angstrom=cutoff)
    direct = phase_grating_from_structure_factors(e, a, g, **options)
    explicit_direct = phase_grating_from_structure_factors(e, a, g, method="direct", **options)
    fft = phase_grating_from_structure_factors(e, a, g, method="fft", **options)
    np.testing.assert_array_equal(direct, explicit_direct)
    np.testing.assert_allclose(fft, direct, atol=8e-15, rtol=2e-14)
    np.testing.assert_array_equal(e, original_e)
    np.testing.assert_array_equal(a, original_a)
    # Huge coefficients at both excluded Nyquist boundaries must affect neither method.
    e[~legacy_nyquist_mask(g)] = 1e100
    np.testing.assert_allclose(phase_grating_from_structure_factors(e, a, g, method="fft", **options),
                               fft, atol=0, rtol=0)


@pytest.mark.parametrize("method", ["FFT", "auto", "", None, True, []])
def test_invalid_phase_grating_method(method):
    c = np.zeros((4, 4), dtype=np.complex128)
    with pytest.raises(ValueError, match="method"):
        phase_grating_from_structure_factors(c, c, grid(), voltage_kv=200,
                                             dsap_inverse_angstrom=1, method=method)


def test_fft_does_not_bypass_odd_mesh_mask_restriction():
    g = grid((5, 3))
    c = np.zeros((3, 5), dtype=np.complex128)
    with pytest.raises(ValueError, match="even"):
        phase_grating_from_structure_factors(c, c, g, voltage_kv=200, dsap_inverse_angstrom=1, method="fft")
