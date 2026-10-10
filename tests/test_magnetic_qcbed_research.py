"""Focused physical and mathematical checks for research-only magnetic prototype."""
import sys
from pathlib import Path

import numpy as np
import pytest
from qcbed.magnetism.research import mu_0

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from qcbed.magnetism.research import (  # noqa: E402
    ElectronBeam,
    diffraction_intensity,
    finite_difference_sensitivity,
    marginal_magnetic_fisher,
    periodic_field_from_magnetization,
    scalar_magnetic_multislice,
)


def test_relativistic_beam_is_physical():
    beam = ElectronBeam(200_000)
    assert 0 < beam.velocity_m_s < 299_792_458
    assert 1e12 < beam.k0_rad_per_m < 1e13
    assert beam.interaction_constant_rad_per_v_m > 0


def test_uniform_magnetization_only_q0_is_separated():
    m = np.zeros((4, 12, 10, 3))
    m[..., 0] = 20_000.0
    fld = periodic_field_from_magnetization(m, (2e-10, 1e-10, 1e-10))
    assert np.allclose(fld.mean_magnetization_a_per_m, [20_000, 0, 0])
    assert np.max(np.abs(fld.a_t_m)) < 1e-18
    assert np.max(np.abs(fld.b_t)) < 1e-12


def test_transverse_sinusoidal_magnetization_has_analytic_Az_and_Bx():
    nz, ny, nx = 4, 32, 16
    dz, dy, dx = 2e-10, 1e-10, 1e-10
    M0 = 2e5
    q = 2 * np.pi / (ny * dy)
    y = np.arange(ny) * dy
    m = np.zeros((nz, ny, nx, 3))
    m[..., 0] = M0 * np.cos(q * y)[None, :, None]
    fld = periodic_field_from_magnetization(m, (dz, dy, dx))
    expected_az = mu_0 * M0 / q * np.sin(q * y)
    expected_bx = mu_0 * M0 * np.cos(q * y)
    assert np.allclose(fld.a_t_m[..., 2], expected_az[None, :, None], atol=1e-20)
    assert np.allclose(fld.b_t[..., 0], expected_bx[None, :, None], atol=1e-10)
    assert np.max(np.abs(fld.a_t_m[..., :2])) < 1e-22
    # A and B odd under moment reversal
    reversed_field = periodic_field_from_magnetization(-m, (dz, dy, dx))
    assert np.allclose(reversed_field.a_t_m, -fld.a_t_m)


def test_magnetic_multislice_has_correct_zero_field_limit_and_conserves_norm():
    rng = np.random.default_rng(4)
    wave = rng.standard_normal((32, 32)) + 1j * rng.standard_normal((32, 32))
    z = np.zeros((3, 32, 32))
    out = scalar_magnetic_multislice(wave, z, z, np.array([1e-9]*3), (1e-10,1e-10), ElectronBeam(200000))
    assert np.isclose(np.sum(abs(out)**2), np.sum(abs(wave)**2), rtol=1e-12)
    assert np.isclose(np.sum(diffraction_intensity(out)), np.sum(abs(wave)**2), rtol=1e-12)
    # Uniform incident wave and no scattering stays uniform.
    plane = np.ones((32,32), dtype=complex)
    result = scalar_magnetic_multislice(plane, z, z, np.array([1e-9]*3), (1e-10,1e-10), ElectronBeam(200000))
    assert np.allclose(result, plane, atol=1e-12)


def test_toy_magnetic_grating_produces_diffracted_intensity():
    nz, ny, nx = 2, 64, 32
    dz, dy, dx = 0.5e-9, 1e-10, 1e-10
    y = np.arange(ny)*dy
    m = np.zeros((nz,ny,nx,3))
    m[...,0] = 2e5*np.cos(2*np.pi*y/(ny*dy))[None,:,None]
    fld = periodic_field_from_magnetization(m, (dz,dy,dx))
    out = scalar_magnetic_multislice(
        np.ones((ny,nx), dtype=complex),
        np.zeros((nz,ny,nx)), fld.a_t_m[...,2],
        np.full(nz,dz), (dy,dx), ElectronBeam(200000)
    )
    spots = diffraction_intensity(out)
    cy, cx = ny//2, nx//2
    assert spots[cy-1,cx] > 0
    assert spots[cy+1,cx] > 0
    assert np.isclose(np.sum(spots), ny*nx, rtol=1e-12)


def test_fisher_detects_perfectly_degenerate_magnetic_and_nuisance():
    # p0 and p1 have IDENTICAL responses: no independent magnetic information.
    sim = lambda p: np.array([2*p[0]+2*p[1], 3*p[0]+3*p[1], p[0]+p[1]])
    res = finite_difference_sensitivity(sim, np.array([1.,1.]), np.array([0.001,0.001]),
                                       np.ones(3), np.ones(3, dtype=int))
    assert res.numerical_rank == 1
    profiled = marginal_magnetic_fisher(res.fisher, [0])
    assert np.abs(profiled[0,0]) < 1e-9
    with pytest.raises(ValueError):
        finite_difference_sensitivity(sim, np.array([1.,1.]), np.array([0.1,0.1]),
                                      np.ones(3), np.array([1,.5,1]))


def test_nio_toy_type_ii_zero_net_and_canting():
    from qcbed.magnetism.research import nio_type_ii_sites, gaussian_magnetization_grid, mu_B
    mu_b = mu_B
    p0, m0, lengths = nio_type_ii_sites()
    assert p0.shape == (32,3) and m0.shape == (32,3)
    assert np.linalg.norm(m0.sum(axis=0)) < 1e-15 * mu_b
    angle = 0.08
    pos, moments, cell = nio_type_ii_sites(canting_angle_rad=angle)
    cant_direction = np.array([1.,-1.,0.])/np.sqrt(2)
    assert np.allclose(moments.sum(axis=0), 32 * 1.9 * mu_b*np.sin(angle)*cant_direction, atol=1e-28)
    grid = gaussian_magnetization_grid((16,16,16), cell, pos, moments, width_sigma_m=0.7e-10)
    assert grid.shape == (16,16,16,3)
    assert np.allclose(grid.mean(axis=(0,1,2)), moments.sum(axis=0)/np.prod(cell), atol=1e-4)
