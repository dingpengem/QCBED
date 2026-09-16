import cmath
import math

import numpy as np
import pytest

from qcbed.fft_backend import (centered_to_unshifted, unshifted_to_centered,
    reciprocal_to_real_sum, real_to_reciprocal_coefficients)
from qcbed.grids import reciprocal_index_grid, oshift
from qcbed.reciprocal import (DirectLattice2D, ReciprocalLattice2D,
    ReciprocalIndexGrid, ReciprocalVectorGrid)
from qcbed.phase_grating import projected_potential_direct, projected_potential_fft, legacy_nyquist_mask


MESHES = [(4, 4), (5, 5), (6, 4), (5, 4), (4, 5), (5, 3), (1, 1), (1, 4)]
ATOL = 8e-14  # Accumulated trigonometric/summation roundoff on <=25 points.
RTOL = 2e-14


def grid(nx, ny, centered):
    return ReciprocalVectorGrid(ReciprocalLattice2D(DirectLattice2D(20, 10, 90)),
                                ReciprocalIndexGrid(nx, ny, centered))


def nested_sum(values, h, k):
    ny, nx = values.shape
    out = np.zeros_like(values)
    for row in range(ny):
        for col in range(nx):
            for j in range(ny):
                for i in range(nx):
                    out[row, col] += values[j, i]*cmath.exp(-2j*math.pi*(h[j, i]*col/nx+k[j, i]*row/ny))
    return out


@pytest.mark.parametrize("nx,ny", MESHES)
def test_exact_order_mapping_and_no_value_loss(nx, ny):
    h, k = reciprocal_index_grid(nx, ny)
    values = (1+np.arange(nx*ny).reshape(ny, nx)+1j).astype(np.complex128)
    bins = centered_to_unshifted(values)
    for row in range(ny):
        for col in range(nx):
            assert bins[int(k[row, col]) % ny, int(h[row, col]) % nx] == values[row, col]
    np.testing.assert_array_equal(unshifted_to_centered(bins), values)
    # Establish the equivalence of NumPy's shift permutations, including odd.
    np.testing.assert_array_equal(bins, np.fft.ifftshift(values, axes=(0, 1)))
    np.testing.assert_array_equal(unshifted_to_centered(bins), np.fft.fftshift(bins, axes=(0, 1)))
    uh, uk = reciprocal_index_grid(nx, ny, centered=False)
    np.testing.assert_array_equal(centered_to_unshifted(h.astype(np.complex128)), uh)
    np.testing.assert_array_equal(centered_to_unshifted(k.astype(np.complex128)), uk)
    if nx % 2 == ny % 2 == 0:
        np.testing.assert_array_equal(bins, oshift(values))


def test_explicit_odd_centering_permutation():
    centered = np.array([[10, 20, 30, 40, 50]], dtype=np.complex128)
    np.testing.assert_array_equal(centered_to_unshifted(centered), [[30, 40, 50, 10, 20]])


@pytest.mark.parametrize("ordering", ["centered", "unshifted"])
def test_dc_unscaled_and_inverse_normalization(ordering):
    c = np.zeros((4, 6), dtype=np.complex128)
    c[(2, 3) if ordering == "centered" else (0, 0)] = 2-3j
    np.testing.assert_array_equal(reciprocal_to_real_sum(c, ordering=ordering), np.full((4, 6), 2-3j))
    np.testing.assert_array_equal(real_to_reciprocal_coefficients(np.full((4, 6), 2-3j), ordering=ordering), c)
    # Spatial delta gives every reciprocal coefficient exactly 1/(4*6).
    delta = np.zeros((4, 6), dtype=np.complex128)
    delta[0, 0] = 24
    np.testing.assert_array_equal(real_to_reciprocal_coefficients(delta, ordering=ordering), np.ones((4, 6)))


@pytest.mark.parametrize("h,k", [(1, 0), (-1, 0), (0, 1), (0, -1)])
@pytest.mark.parametrize("ordering", ["centered", "unshifted"])
def test_signed_sinusoid(h, k, ordering):
    hh, kk = reciprocal_index_grid(4, 4, centered=ordering == "centered")
    c = np.zeros((4, 4), dtype=np.complex128)
    c[(hh == h) & (kk == k)] = 1
    roots = np.array([1, -1j, -1, 1j])
    expected = np.empty((4, 4), dtype=np.complex128)
    for row in range(4):
        for col in range(4):
            expected[row, col] = roots[(h*col+k*row) % 4]
    np.testing.assert_array_equal(reciprocal_to_real_sum(c, ordering=ordering), expected)


def test_friedel_pair_cosine():
    c = np.zeros((3, 4), dtype=np.complex128)
    c[0, 1] = c[0, 3] = 1
    np.testing.assert_array_equal(reciprocal_to_real_sum(c), [[2, 0, -2, 0]]*3)


def test_negative_nyquist_preserved_until_explicit_mask():
    g = grid(6, 4, True)
    c = np.zeros((4, 6), dtype=np.complex128)
    c[0, 3], c[2, 0] = 7+2j, -3+4j
    bins = centered_to_unshifted(c)
    assert bins[2, 0] == 7+2j and bins[0, 3] == -3+4j
    assert np.any(reciprocal_to_real_sum(bins))
    np.testing.assert_array_equal(c*legacy_nyquist_mask(g), np.zeros_like(c))


@pytest.mark.parametrize("nx,ny", MESHES)
@pytest.mark.parametrize("centered", [True, False])
def test_seeded_complex_reference_and_round_trip(nx, ny, centered):
    rng = np.random.default_rng(20260916)
    c = (rng.normal(size=(ny, nx))+1j*rng.normal(size=(ny, nx))).astype(np.complex128)
    before = c.copy()
    c.setflags(write=False)
    g = grid(nx, ny, centered)
    ordering = "centered" if centered else "unshifted"
    fft = reciprocal_to_real_sum(c, ordering=ordering)
    explicit = nested_sum(c, g.indices.h, g.indices.k)
    direct = projected_potential_direct(c, g)
    np.testing.assert_allclose(fft, explicit, rtol=RTOL, atol=ATOL)
    np.testing.assert_allclose(fft, direct, rtol=RTOL, atol=ATOL)
    np.testing.assert_array_equal(projected_potential_fft(c, g), fft)
    restored = real_to_reciprocal_coefficients(fft, ordering=ordering)
    np.testing.assert_allclose(restored, c, rtol=2e-15, atol=2e-15)
    # Check the opposite round trip for an arbitrary spatial array as well.
    np.testing.assert_allclose(reciprocal_to_real_sum(real_to_reciprocal_coefficients(c, ordering=ordering),
                               ordering=ordering), c, rtol=2e-15, atol=2e-15)
    for out in (fft, restored, centered_to_unshifted(c), unshifted_to_centered(c)):
        assert out.dtype == np.complex128 and out.shape == (ny, nx)
        assert out.flags.c_contiguous and not np.shares_memory(c, out)
    np.testing.assert_array_equal(c, before)


@pytest.mark.parametrize("operation", [centered_to_unshifted, unshifted_to_centered,
    reciprocal_to_real_sum, real_to_reciprocal_coefficients])
@pytest.mark.parametrize("bad", [np.zeros((2, 3)), np.zeros((2, 3), dtype=np.complex64), [[1j]],
    np.zeros(6, dtype=np.complex128), np.zeros((0, 3), dtype=np.complex128),
    np.full((2, 3), np.nan, dtype=np.complex128), np.full((2, 3), np.inf, dtype=np.complex128)])
def test_invalid_arrays(operation, bad):
    with pytest.raises((TypeError, ValueError), match="complex128|shape|finite"):
        operation(bad)


@pytest.mark.parametrize("bad", ["fft", "Centered", "", None, True, 0, []])
def test_invalid_ordering(bad):
    for op in (reciprocal_to_real_sum, real_to_reciprocal_coefficients):
        with pytest.raises(ValueError, match="ordering"):
            op(np.zeros((2, 3), dtype=np.complex128), ordering=bad)


def test_fft_grid_shape_and_finite_result():
    with pytest.raises(ValueError, match="shape"):
        projected_potential_fft(np.zeros((3, 4), dtype=np.complex128), grid(4, 4, True))
    with pytest.raises(TypeError, match="reciprocal_grid"):
        projected_potential_fft(np.zeros((3, 4), dtype=np.complex128), None)
    with pytest.raises(ValueError, match="finite"):
        reciprocal_to_real_sum(np.full((4, 4), 1e308, dtype=np.complex128))
