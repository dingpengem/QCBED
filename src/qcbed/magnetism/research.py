"""Research-only magnetic QCBED prototype (SI units throughout).

Physics model: scalar, high-energy paraxial multislice in Coulomb gauge,
retaining only the longitudinal periodic magnetic vector potential A_z.
The transverse A field, spin/Pauli term, uniform applied B, sample-boundary
fields and all experimental detector effects are NOT modeled here.

This module is standalone: it has NOT been integrated with the QCBED repo.
Use it to check units, mathematical identities, toy magnetic diffraction,
and finite-difference identifiability before a physics integration.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

import numpy as np
from math import pi

# SI constants; exact: c, e, Planck constant; CODATA 2022: m_e, mu_0.
c = 299_792_458.0
e = 1.602_176_634e-19
hbar = 6.626_070_15e-34 / (2 * pi)
m_e = 9.109_383_713_9e-31
mu_0 = 1.256_637_061_27e-6
mu_B = e * hbar / (2 * m_e)


@dataclass(frozen=True)
class ElectronBeam:
    accelerating_voltage_v: float

    def __post_init__(self) -> None:
        if not np.isfinite(self.accelerating_voltage_v) or self.accelerating_voltage_v <= 0:
            raise ValueError("Accelerating voltage must be a positive finite number in V")

    @property
    def total_energy_j(self) -> float:
        return m_e * c**2 + e * self.accelerating_voltage_v

    @property
    def momentum_kg_m_s(self) -> float:
        return np.sqrt(self.total_energy_j**2 - (m_e * c**2) ** 2) / c

    @property
    def velocity_m_s(self) -> float:
        return self.momentum_kg_m_s * c**2 / self.total_energy_j

    @property
    def k0_rad_per_m(self) -> float:
        return self.momentum_kg_m_s / hbar

    @property
    def interaction_constant_rad_per_v_m(self) -> float:
        """Electric scalar-potential phase coefficient +e/(hbar*v)."""
        return e / (hbar * self.velocity_m_s)


@dataclass(frozen=True)
class PeriodicMagneticField:
    """Only nonzero periodic spatial Fourier harmonics (q != 0)."""
    a_t_m: np.ndarray  # [nz,ny,nx,3], Coulomb-gauge vector potential, T*m
    b_t: np.ndarray  # [nz,ny,nx,3], periodic induction, tesla
    mean_magnetization_a_per_m: np.ndarray  # q=0 term, NOT included in A/B


def _validate_spacing(spacing_m: tuple[float, float, float]) -> tuple[float, float, float]:
    if len(spacing_m) != 3 or not np.all(np.isfinite(spacing_m)) or any(d <= 0 for d in spacing_m):
        raise ValueError("spacing_m=(dz,dy,dx) must contain finite positive meters")
    return tuple(float(v) for v in spacing_m)


def periodic_field_from_magnetization(
    magnetization_a_per_m: np.ndarray,
    spacing_m: tuple[float, float, float],
) -> PeriodicMagneticField:
    """Compute A_q = i*mu0*(q x M_q)/|q|^2, q!=0 (Coulomb gauge).

    Input magnetization M is real [nz,ny,nx,(Mx,My,Mz)] in A/m.
    Axes are orthogonal Cartesian x,y,z and FFT axes are (z,y,x).
    q=0 is set to zero IN BOTH A and B; use an external model for it.
    Periodic boundary conditions, including z, are required.
    """
    m = np.asarray(magnetization_a_per_m, dtype=np.float64)
    if m.ndim != 4 or m.shape[-1] != 3 or min(m.shape[:3]) < 2 or not np.all(np.isfinite(m)):
        raise ValueError("magnetization must be finite [nz,ny,nx,3], all dimensions >=2")
    dz, dy, dx = _validate_spacing(spacing_m)
    nz, ny, nx, _ = m.shape
    qz = 2 * np.pi * np.fft.fftfreq(nz, d=dz)[:, None, None]
    qy = 2 * np.pi * np.fft.fftfreq(ny, d=dy)[None, :, None]
    qx = 2 * np.pi * np.fft.fftfreq(nx, d=dx)[None, None, :]
    q = np.stack(np.broadcast_arrays(qx, qy, qz), axis=-1)
    q2 = np.sum(q * q, axis=-1)
    nonzero = q2 > 0
    m_hat = np.fft.fftn(m, axes=(0, 1, 2))
    a_hat = np.zeros_like(m_hat)
    a_hat[nonzero] = (1j * mu_0 * np.cross(q, m_hat)[nonzero] / q2[nonzero, None])
    # Derive B using the same discrete spectral operator as A: B_q = i q x A_q.
    b_hat = 1j * np.cross(q, a_hat)
    a = np.fft.ifftn(a_hat, axes=(0, 1, 2))
    b = np.fft.ifftn(b_hat, axes=(0, 1, 2))
    # Nyquist modes in even-sized grids can produce complex-valued Fourier
    # derivatives. Reject them if inconsistent instead of silently discarding.
    scale_a = max(float(np.max(np.abs(a.real))), 1e-30)
    scale_b = max(float(np.max(np.abs(b.real))), 1e-30)
    if np.max(np.abs(a.imag)) > 1e-9 * scale_a or np.max(np.abs(b.imag)) > 1e-9 * scale_b:
        raise ValueError("Complex field from unresolved/Nyquist Fourier mode; refine spatial sampling")
    return PeriodicMagneticField(
        a_t_m=a.real.copy(),
        b_t=b.real.copy(),
        mean_magnetization_a_per_m=np.mean(m, axis=(0, 1, 2)),
    )


def scalar_magnetic_multislice(
    incident_wave: np.ndarray,
    electrostatic_phase_rad: np.ndarray,
    az_t_m: np.ndarray,
    dz_m: np.ndarray,
    pixel_size_yx_m: tuple[float, float],
    beam: ElectronBeam,
) -> np.ndarray:
    """Strang-split paraxial propagator, periodic A_z approximation ONLY.

    electrostatic_phase_rad[s,y,x] = +(e/(hbar*v))*integral(phi dz),
    i.e. the *phase in radians*, not a voltage or QCBED E_g coefficient.
    az_t_m[s,y,x] = Coulomb-gauge A_z at slice center, in T*m.
    dz_m[s] is slice thickness in meters. Convention: incident[y,x].

    T_s=exp(+i electric_phase - i e/hbar * A_z * dz).
    P_s/2(k)=exp(-i dz |k_perp|^2/(4*k0)).
    No use of A_perp: this is NOT a fully gauge-invariant field propagator.
    """
    wave = np.asarray(incident_wave, dtype=np.complex128).copy()
    phase = np.asarray(electrostatic_phase_rad, dtype=np.float64)
    az = np.asarray(az_t_m, dtype=np.float64)
    depths = np.asarray(dz_m, dtype=np.float64)
    if wave.ndim != 2 or phase.ndim != 3 or az.shape != phase.shape or phase.shape[1:] != wave.shape:
        raise ValueError("wave[y,x] and phase/az[s,y,x] shapes must match")
    if depths.shape != (phase.shape[0],) or not np.all(np.isfinite(depths)) or np.any(depths <= 0):
        raise ValueError("dz_m must have one positive finite thickness per slice")
    if not np.all(np.isfinite(phase)) or not np.all(np.isfinite(az)):
        raise ValueError("Nonfinite phase or A_z")
    dy, dx = pixel_size_yx_m
    if dy <= 0 or dx <= 0 or not np.isfinite(dy + dx):
        raise ValueError("pixel sizes must be positive finite meters")
    ky = 2 * np.pi * np.fft.fftfreq(wave.shape[0], d=dy)[:, None]
    kx = 2 * np.pi * np.fft.fftfreq(wave.shape[1], d=dx)[None, :]
    k2 = kx * kx + ky * ky
    for s, step in enumerate(depths):
        half_p = np.exp(-1j * step * k2 / (4 * beam.k0_rad_per_m))
        wave = np.fft.ifft2(np.fft.fft2(wave) * half_p)
        wave *= np.exp(1j * (phase[s] - (e / hbar) * az[s] * step))
        wave = np.fft.ifft2(np.fft.fft2(wave) * half_p)
    return wave


def diffraction_intensity(exit_wave: np.ndarray) -> np.ndarray:
    """Unitary FFT diffraction intensity [arbitrary incident normalization]."""
    return np.abs(np.fft.fftshift(np.fft.fft2(exit_wave, norm="ortho"))) ** 2


@dataclass(frozen=True)
class SensitivityResult:
    jacobian: np.ndarray
    fisher: np.ndarray
    singular_values: np.ndarray
    numerical_rank: int
    marginal_covariance: np.ndarray


def finite_difference_sensitivity(
    simulator: Callable[[np.ndarray], np.ndarray],
    parameter_values: np.ndarray,
    steps: np.ndarray,
    sigma: np.ndarray,
    final_weight: np.ndarray,
    *,
    rank_rtol: float = 1e-10,
) -> SensitivityResult:
    """Finite-difference Fisher analysis; sigma is std. dev., NOT variance.

    Uses the existing QCBED idea: final_weight is binary 0/1 and sigma is
    the statistical denominator. Assumes locally independent Gaussian errors.
    Rank is evaluated on the whitened Jacobian, not unscaled parameter space;
    choosing meaningful parameter units/scales is essential.
    The pseudoinverse covariance is not a valid uncertainty for null modes.
    """
    p = np.asarray(parameter_values, dtype=float)
    h = np.asarray(steps, dtype=float)
    sig = np.asarray(sigma, dtype=float)
    mask = np.asarray(final_weight)
    if p.ndim != 1 or h.shape != p.shape or np.any(h <= 0) or not np.all(np.isfinite(p + h)):
        raise ValueError("finite 1-D parameters and positive finite differences required")
    if not np.all(np.isin(mask, [0, 1])):
        raise ValueError("Final weights must be binary 0 or 1")
    if mask.shape != sig.shape or not np.any(mask == 1):
        raise ValueError("sigma and mask shapes must match and contain active pixels")
    if not np.all(np.isfinite(sig[mask == 1])) or np.any(sig[mask == 1] <= 0):
        raise ValueError("sigma must be positive and finite on all active pixels")
    cols = []
    for j, step in enumerate(h):
        plus = p.copy(); plus[j] += step
        minus = p.copy(); minus[j] -= step
        y_plus = np.asarray(simulator(plus), dtype=float)
        y_minus = np.asarray(simulator(minus), dtype=float)
        if y_plus.shape != sig.shape or y_minus.shape != sig.shape:
            raise ValueError("simulator output and sigma shapes must match")
        cols.append(((y_plus - y_minus) / (2 * step))[mask == 1])
    jac = np.stack(cols, axis=1)
    if not np.all(np.isfinite(jac)):
        raise ValueError("nonfinite intensity Jacobian")
    white = jac / sig[mask == 1, None]
    fisher = white.T @ white
    svals = np.linalg.svd(white, compute_uv=False)
    rank = int(np.sum(svals > rank_rtol * svals[0])) if svals.size and svals[0] > 0 else 0
    return SensitivityResult(jac, fisher, svals, rank, np.linalg.pinv(fisher, rcond=rank_rtol))


def marginal_magnetic_fisher(
    fisher: np.ndarray,
    magnetic_indices: list[int] | tuple[int, ...],
) -> np.ndarray:
    """Schur complement, nuisance parameters profiled out (local Gaussian)."""
    F = np.asarray(fisher, dtype=float)
    inds = np.asarray(magnetic_indices, dtype=int)
    if F.ndim != 2 or F.shape[0] != F.shape[1] or inds.ndim != 1 or len(set(inds.tolist())) != inds.size or np.any(inds < 0) or np.any(inds >= F.shape[0]):
        raise ValueError("Invalid Fisher matrix / magnetic indices")
    nuis = np.asarray([j for j in range(F.shape[0]) if j not in set(inds.tolist())])
    Fmm = F[np.ix_(inds, inds)]
    if not len(nuis):
        return Fmm
    Fmn = F[np.ix_(inds, nuis)]
    Fnn = F[np.ix_(nuis, nuis)]
    # Pseudoinverse when nuisance block is singular. Singular magnetic
    # Schur complement signifies non-identifiability; interpret with caution.
    return Fmm - Fmn @ np.linalg.pinv(Fnn) @ Fmn.T


def nio_type_ii_sites(
    lattice_parameter_m: float = 4.17e-10,
    moment_mu_b_per_ni: float = 1.9,
    canting_angle_rad: float = 0.0,
    afm_axis_cubic: tuple[float, float, float] = (1.0, 1.0, -2.0),
    uniform_cant_axis_cubic: tuple[float, float, float] = (1.0, -1.0, 0.0),
) -> tuple[np.ndarray, np.ndarray, tuple[float, float, float]]:
    """Illustrative cubic NiO 2x2x2 magnetic supercell, 32 Ni moments.

    Returns (positions_xyz_m [32,3], moments_xyz_A_m2 [32,3], lengths_xyz_m).
    Type-II alternating (111) Ni layers: eta=(-1)^(x+y+z) for positions
    expressed in units of the conventional cubic lattice parameter a.
    Oxygen carries NO moment in this simple model. It omits real distortion,
    S/T domain mixtures, anisotropic magnetic form factors and orbital moments.
    This is a geometry/units fixture; not a calibrated NiO scattering model.
    """
    mu_b = mu_B
    if lattice_parameter_m <= 0 or moment_mu_b_per_ni < 0 or not np.isfinite(canting_angle_rad):
        raise ValueError("Invalid NiO lattice/moment/canting parameters")
    n = np.asarray(afm_axis_cubic, dtype=float)
    h = np.asarray(uniform_cant_axis_cubic, dtype=float)
    if n.shape != (3,) or h.shape != (3,) or np.linalg.norm(n) == 0 or np.linalg.norm(h) == 0:
        raise ValueError("AFM and canting directions must be nonzero 3D vectors")
    n /= np.linalg.norm(n)
    h /= np.linalg.norm(h)
    if abs(n @ h) > 1e-10:
        raise ValueError("This symmetric canting model requires perpendicular n and h")
    fcc = np.array([[0, 0, 0], [0, .5, .5], [.5, 0, .5], [.5, .5, 0.]])
    positions = []
    moments = []
    for iz in range(2):
        for iy in range(2):
            for ix in range(2):
                for f in fcc:
                    pos = f + np.array([ix, iy, iz], dtype=float)
                    s = int(round(float(np.sum(pos))))
                    eta = 1 if s % 2 == 0 else -1
                    moments.append(mu_b * moment_mu_b_per_ni * (
                        eta * np.cos(canting_angle_rad) * n +
                        np.sin(canting_angle_rad) * h
                    ))
                    positions.append(lattice_parameter_m * pos)
    a2 = 2 * lattice_parameter_m
    return np.asarray(positions), np.asarray(moments), (a2, a2, a2)


def gaussian_magnetization_grid(
    grid_shape_zyx: tuple[int, int, int],
    cell_lengths_xyz_m: tuple[float, float, float],
    positions_xyz_m: np.ndarray,
    moments_xyz_a_m2: np.ndarray,
    width_sigma_m: float,
) -> np.ndarray:
    """Periodic Gaussian moments (toy form factor), returns M[z,y,x,3] A/m.

    M_q = exp(-sigma^2*q^2/2)/V * sum_j mu_j exp(-i q.r_j).
    The Gaussian is for numerical fixtures, NOT an element-specific Ni form
    factor. For quantitative NiO use DFT or calibrated magnetic form factors.
    """
    ns = tuple(int(v) for v in grid_shape_zyx)
    lengths = np.asarray(cell_lengths_xyz_m, dtype=float)
    pos = np.asarray(positions_xyz_m, dtype=float)
    moments = np.asarray(moments_xyz_a_m2, dtype=float)
    if len(ns) != 3 or min(ns) < 2 or lengths.shape != (3,) or np.any(lengths <= 0):
        raise ValueError("grid and cell lengths must be three positive dimensions")
    if pos.ndim != 2 or pos.shape[1] != 3 or moments.shape != pos.shape or not np.all(np.isfinite(pos + moments)):
        raise ValueError("positions and moments must have matching [site,3] shapes")
    if width_sigma_m <= 0 or not np.isfinite(width_sigma_m):
        raise ValueError("Gaussian sigma must be positive finite meters")
    nz, ny, nx = ns
    lx, ly, lz = lengths
    qz = 2*np.pi*np.fft.fftfreq(nz, d=lz/nz)[:,None,None]
    qy = 2*np.pi*np.fft.fftfreq(ny, d=ly/ny)[None,:,None]
    qx = 2*np.pi*np.fft.fftfreq(nx, d=lx/nx)[None,None,:]
    gaussian = np.exp(-0.5 * width_sigma_m**2 * (qx*qx + qy*qy + qz*qz))
    mh = np.zeros((nz,ny,nx,3), dtype=complex)
    for r, mu in zip(pos, moments):
        factor = np.exp(-1j*(qx*r[0] + qy*r[1] + qz*r[2])) * gaussian
        mh += factor[...,None] * mu
    mh *= (nz*ny*nx)/np.prod(lengths)
    # The real part represents the periodically sampled physical density.
    # Nyquist artifacts from insufficient sampling can remain; convergence
    # with increasing grid must be established before physical predictions.
    result = np.fft.ifftn(mh, axes=(0,1,2)).real
    return result
