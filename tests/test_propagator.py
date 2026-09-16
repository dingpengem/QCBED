import cmath
from dataclasses import FrozenInstanceError, replace
import math

import numpy as np
import pytest

from qcbed.reciprocal import (DirectLattice2D, ReciprocalLattice2D,
    ReciprocalIndexGrid, ReciprocalVectorGrid)
from qcbed.grids import oshift
from qcbed.fft_backend import centered_to_unshifted
from qcbed.propagator import (PropagationParameters, propagation_phase_radians,
    fractional_slice_shift_phase_radians, propagation_boundary_mask, propagation_kernel)


P = PropagationParameters(.025, 2)
ATOL = 3e-12  # Allows ideal-angle vs printed geometry-pi difference.


def grid(a=2, b=4, gamma=90, mesh=(6, 4), centered=True, dilation=1):
    return ReciprocalVectorGrid(ReciprocalLattice2D(DirectLattice2D(a, b, gamma)),
                                ReciprocalIndexGrid(*mesh, centered=centered), dilation)


def value(h, k, p=P, **kwargs):
    return propagation_kernel(grid(**kwargs), p, ordering="centered", boundary="none")[k+2, h+3]


# Independent orthogonal metric: g=(H/2,K/4), t=(th/2,tk/4).
# phase=6.2831853*.025*(|g|^2+2*t.g) for DZ=2, d=1.
HAND_CASES = [
    (1, 0, 0, 0, .9992290362424846+.03925981571423079j),
    (-1, 0, 0, 0, .9992290362424846+.03925981571423079j),
    (0, 1, 0, 0, .9999518089594381+.009817319325932055j),
    (0, -1, 0, 0, .9999518089594381+.009817319325932055j),
    (1, 1, 0, 0, .9987954562079246+.049067674271395065j),
    (1, 0, .2, 0, .99848909745399+.05495017984971928j),
    (-1, 0, .2, 0, .9997224302186348+.023559764806694177j),
    (0, 1, 0, -.2, .9999826511362174+.005890452154296261j),
    (0, -1, 0, -.2, .9999055462888008+.013744035102359588j),
    (1, 1, .2, -.4, .9983792783536255+.05691060142094816j),
]


@pytest.mark.parametrize("h,k,th,tk,expected", HAND_CASES)
def test_hand_calculated_signs_and_tilts(h, k, th, tk, expected):
    p = replace(P, tilt_h_fractional=th, tilt_k_fractional=tk)
    assert value(h, k, p) == pytest.approx(expected, abs=ATOL, rel=0)


def test_zero_thickness_and_origin():
    p = replace(P, slice_thickness_angstrom=0, tilt_h_fractional=.3, tilt_k_fractional=-.7)
    g = grid()
    np.testing.assert_array_equal(propagation_phase_radians(g, p), np.zeros((4, 6)))
    np.testing.assert_array_equal(propagation_kernel(g, p, ordering="centered", boundary="none"), np.ones((4, 6)))
    p = replace(p, slice_thickness_angstrom=2, x_shift_fractional=.31, y_shift_fractional=-.17)
    assert value(0, 0, p) == 1+0j


@pytest.mark.parametrize("axis", ["x", "y"])
def test_nonzero_fractional_shifts_have_negative_phase(axis):
    p = replace(P, slice_thickness_angstrom=0, **{f"{axis}_shift_fractional": .25})
    h, k = (1, 0) if axis == "x" else (0, 1)
    assert value(h, k, p) == pytest.approx(1.7948965149208059e-9-1j, abs=3e-16, rel=0)
    assert value(-h, -k, p) == pytest.approx(1.7948965149208059e-9+1j, abs=3e-16, rel=0)
    phases = fractional_slice_shift_phase_radians(grid(), **{f"{axis}_shift_fractional": .25})
    assert phases[k+2, h+3] == -1.570796325
    assert value(k, h, p) == 1+0j  # Different reciprocal axis.


def test_shift_not_scaled_by_thickness_or_dilation():
    p = replace(P, slice_thickness_angstrom=0, x_shift_fractional=.25, y_shift_fractional=-.5)
    g = grid(dilation=3)
    actual = propagation_kernel(g, p, ordering="centered", boundary="none")
    # H=K=1 gives phase=+TWOPI/4, regardless of dilation.
    assert actual[3, 4] == pytest.approx(1.7948965149208059e-9+1j, abs=3e-16, rel=0)


def test_full_cell_shift_periodicity_is_limited_by_printed_twopi():
    g = grid()
    base = propagation_kernel(g, P, ordering="centered", boundary="none")
    for field in ("x_shift_fractional", "y_shift_fractional"):
        moved = propagation_kernel(g, replace(P, **{field: 1}), ordering="centered", boundary="none")
        # |TWOPI-2*pi|=7.179586e-9; maximum |H|=3 -> 2.154e-8.
        np.testing.assert_allclose(moved, base, atol=2.2e-8, rtol=0)
    assert abs(value(1, 0, replace(P, x_shift_fractional=1))-value(1, 0)) > 7e-9


def test_square_and_dilated_cell():
    assert value(0, 1, a=2, b=2) == pytest.approx(value(1, 0, a=2, b=2), abs=ATOL, rel=0)
    # At fixed internal tilt, metric /d^2 and effective DZ*d leave phase /d.
    actual = propagation_phase_radians(grid(dilation=2), P)
    np.testing.assert_allclose(actual, propagation_phase_radians(grid(), P)/2, rtol=0, atol=1e-16)


@pytest.mark.parametrize("gamma", [60, 120])
@pytest.mark.parametrize("centered", [True, False])
@pytest.mark.parametrize("ideal_angle", [True, False])
def test_full_grid_independent_cartesian_reference(gamma, centered, ideal_angle):
    g = grid(gamma=gamma, centered=centered, dilation=1.25)
    p = PropagationParameters(.025, 3, .2, -.3, .125, -.25)
    # Independent ideal oblique dual vectors for a=2,b=4, sin(gamma)=sqrt(3)/2.
    cosine = (.5 if gamma == 60 else -.5) if ideal_angle else math.cos(3.14159265358*gamma/180)
    sine = math.sqrt(1-cosine*cosine)
    ax, ay = .5/1.25, -cosine/(2*sine)/1.25
    bx, by = 0., 1/(4*sine)/1.25
    tx = p.tilt_h_fractional*ax+p.tilt_k_fractional*bx
    ty = p.tilt_h_fractional*ay+p.tilt_k_fractional*by
    expected = np.empty((4, 6), dtype=np.complex128)
    for row in range(4):
        for col in range(6):
            h, k = g.indices.h[row, col], g.indices.k[row, col]
            gx, gy = h*ax+k*bx, h*ay+k*by
            phase = 6.2831853*(.025/2*(gx*gx+gy*gy+2*(tx*gx+ty*gy))*3*1.25
                               - h*.125 + k*.25)
            expected[row, col] = cmath.exp(1j*phase)
    actual = propagation_kernel(g, p, ordering="centered" if centered else "unshifted", boundary="none")
    # Exact angles differ from source pi by ~1e-11 radians; compare the
    # printed-angle Cartesian identity separately at arithmetic precision.
    np.testing.assert_allclose(actual, expected, atol=1e-11 if ideal_angle else 3e-15, rtol=0)


@pytest.mark.parametrize("input_centered", [True, False])
@pytest.mark.parametrize("ordering", ["centered", "unshifted"])
def test_ordering_and_exact_boundary_zeros(input_centered, ordering):
    g = grid(centered=input_centered)
    actual = propagation_kernel(g, P, ordering=ordering)
    centered = propagation_kernel(grid(), P, ordering="centered")
    expected = centered if ordering == "centered" else oshift(centered)
    np.testing.assert_array_equal(actual, expected)
    zero_mask = np.zeros((4, 6), dtype=bool)
    zero_mask[0 if ordering == "centered" else 2, :] = True
    zero_mask[:, 0 if ordering == "centered" else 3] = True
    assert np.all(actual[zero_mask] == 0j)
    np.testing.assert_allclose(abs(actual[~zero_mask]), 1, atol=2e-16, rtol=0)
    assert actual.dtype == np.complex128 and actual.shape == (4, 6)
    assert actual.flags.c_contiguous


def test_odd_generic_kernel_and_even_only_boundary():
    g = grid(mesh=(5, 3))
    centered = propagation_kernel(g, P, ordering="centered", boundary="none")
    unshifted = propagation_kernel(g, P, ordering="unshifted", boundary="none")
    np.testing.assert_array_equal(unshifted, centered_to_unshifted(centered))
    np.testing.assert_allclose(abs(unshifted), 1, atol=2e-16, rtol=0)
    for call in (lambda: propagation_boundary_mask(g), lambda: propagation_kernel(g, P, ordering="unshifted")):
        with pytest.raises(ValueError, match="even"):
            call()


def test_no_bandwidth_cutoff_or_rmt_and_immutable_inputs():
    g = grid(a=.02, b=.04)
    before = g.indices.h.copy(), g.stl_inverse_angstrom.copy()
    p = replace(P, slice_thickness_angstrom=0)
    k = propagation_kernel(g, p, ordering="centered")
    mask = propagation_boundary_mask(g)
    np.testing.assert_array_equal(k[mask], np.ones(mask.sum()))
    assert g.stl_inverse_angstrom.max() > 6  # Not an atomic-scattering cutoff.
    np.testing.assert_array_equal(g.indices.h, before[0])
    np.testing.assert_array_equal(g.stl_inverse_angstrom, before[1])
    assert not g.indices.h.flags.writeable
    assert propagation_phase_radians(g, p).dtype == np.float64
    with pytest.raises(FrozenInstanceError):
        p.wavelength_angstrom = .1


@pytest.mark.parametrize("field", ["wavelength_angstrom", "slice_thickness_angstrom", "tilt_h_fractional",
    "tilt_k_fractional", "x_shift_fractional", "y_shift_fractional"])
@pytest.mark.parametrize("bad", [float("nan"), float("inf"), -float("inf"), True, ".1"])
def test_invalid_parameters(field, bad):
    with pytest.raises((TypeError, ValueError), match=field):
        replace(P, **{field: bad})


def test_invalid_domains_and_kernel_options():
    for wave in (0, -1):
        with pytest.raises(ValueError, match="wavelength"):
            replace(P, wavelength_angstrom=wave)
    with pytest.raises(ValueError, match="thickness"):
        replace(P, slice_thickness_angstrom=-1)
    with pytest.raises(TypeError, match="parameters"):
        propagation_kernel(grid(), None, ordering="centered")
    with pytest.raises(TypeError, match="reciprocal_grid"):
        propagation_phase_radians(None, P)
    for bad in (None, True, "FFT", []):
        with pytest.raises(ValueError, match="ordering"):
            propagation_kernel(grid(), P, ordering=bad)
        with pytest.raises(ValueError, match="boundary"):
            propagation_kernel(grid(), P, ordering="centered", boundary=bad)


def test_unrepresentable_phases_rejected():
    with pytest.raises(ValueError, match="phase.*finite"):
        propagation_phase_radians(grid(), replace(P, wavelength_angstrom=1e308, slice_thickness_angstrom=1e308))
    with pytest.raises(ValueError, match="phase.*finite"):
        fractional_slice_shift_phase_radians(grid(), x_shift_fractional=1e308)


@pytest.mark.parametrize("field", ["x_shift_fractional", "y_shift_fractional"])
def test_shift_function_validates_its_own_inputs(field):
    with pytest.raises(ValueError, match=field):
        fractional_slice_shift_phase_radians(grid(), **{field: float("nan")})
