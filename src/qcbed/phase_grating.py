"""Direct and NumPy Fourier paths for PGMAP; no legacy FFT equivalence claim."""

from __future__ import annotations

import math
from typing import Literal

import numpy as np
from numpy.typing import NDArray

from .absorption import combine_elastic_absorptive_structure_factors
from .electron import electron_wavelength_angstrom
from .fft_backend import reciprocal_to_real_sum
from .reciprocal import ReciprocalVectorGrid
from .structure_factors import _finite


ComplexGrid = NDArray[np.complex128]


def _grid(value: ReciprocalVectorGrid) -> ReciprocalVectorGrid:
    if not isinstance(value, ReciprocalVectorGrid):
        raise TypeError("reciprocal_grid must be ReciprocalVectorGrid")
    return value


def _complex_grid(value: ComplexGrid, name: str) -> ComplexGrid:
    if not isinstance(value, np.ndarray) or value.dtype != np.complex128:
        raise TypeError(f"{name} must be a complex128 NumPy array")
    if value.ndim != 2 or 0 in value.shape:
        raise ValueError(f"{name} must have nonempty shape (mesh_y, mesh_x)")
    if not np.all(np.isfinite(value)):
        raise ValueError(f"{name} must be finite")
    return value


def reciprocal_cutoff_mask(
    reciprocal_grid: ReciprocalVectorGrid, dsap_inverse_angstrom: float,
) -> NDArray[np.bool_]:
    """Keep STL<=DSAP, including equality and g=0; no Nyquist/detector mask."""
    grid = _grid(reciprocal_grid)
    cutoff = _finite(dsap_inverse_angstrom, "dsap_inverse_angstrom", nonnegative=True)
    return grid.stl_inverse_angstrom <= cutoff


def legacy_nyquist_mask(reciprocal_grid: ReciprocalVectorGrid) -> NDArray[np.bool_]:
    """Keep SFMAP2's interior H/K rectangle on even meshes.

    Exclude H=-mesh_x/2 or K=-mesh_y/2 in either supplied ordering. This
    corresponds to the first centered row/column, not real-space pixels.
    Odd dimensions are unsupported by the legacy half-swap convention.
    """
    indices = _grid(reciprocal_grid).indices
    if indices.mesh_x % 2 or indices.mesh_y % 2:
        raise ValueError("legacy Nyquist suppression requires even mesh_x and mesh_y")
    return (indices.h != -indices.mesh_x//2) & (indices.k != -indices.mesh_y//2)


def projected_potential_direct(
    coefficients: ComplexGrid, reciprocal_grid: ReciprocalVectorGrid,
) -> ComplexGrid:
    """Unnormalized negative-sign Fourier sum on x=col/Nx, y=row/Ny.

    V[y,x] = sum_HK C[H,K]*exp(-2*pi*i*(H*x+K*y)). C is the already
    combined source-scaled elastic+i*absorptive coefficient grid. Applies
    no masks, no extra cell-area/thickness factor and no output shift.
    This serial reference follows the mathematical CDFT2D specification;
    historical executable equivalence still requires original build validation.
    """
    grid = _grid(reciprocal_grid)
    values = _complex_grid(coefficients, "coefficients")
    if values.shape != grid.indices.h.shape:
        raise ValueError("coefficients shape must match reciprocal_grid (mesh_y, mesh_x)")
    ny, nx = values.shape
    x = np.arange(nx, dtype=np.float64)[None, :] / nx
    y = np.arange(ny, dtype=np.float64)[:, None] / ny
    result = np.zeros((ny, nx), dtype=np.complex128)
    try:
        with np.errstate(over="raise", invalid="raise"):
            for row, col in np.ndindex(values.shape):
                phase = -math.tau * (grid.indices.h[row, col]*x + grid.indices.k[row, col]*y)
                result += values[row, col] * (np.cos(phase) + 1j*np.sin(phase))
    except FloatingPointError as exc:
        raise ValueError("projected potential must remain finite in float64") from exc
    return result


def projected_potential_fft(
    coefficients: ComplexGrid, reciprocal_grid: ReciprocalVectorGrid,
) -> ComplexGrid:
    """NumPy CPU equivalent of projected_potential_direct, with no masks.

    Negative exponent, unnormalized sum; input ordering comes from the grid.
    Output uses (row/Ny,col/Nx) spatial samples with no shift or transpose.
    Validated against the Python direct sum, not a compiled historical CDFT executable.
    """
    grid = _grid(reciprocal_grid)
    values = _complex_grid(coefficients, "coefficients")
    if values.shape != grid.indices.h.shape:
        raise ValueError("coefficients shape must match reciprocal_grid (mesh_y, mesh_x)")
    ordering = "centered" if grid.indices.centered else "unshifted"
    return reciprocal_to_real_sum(values, ordering=ordering)


def electron_interaction_constant(voltage_kv: float) -> float:
    """Source SIGMA in inverse (volt*angstrom), using derived wavelength.

    SIGMA=0.020886416*WAVEL*(1+1.9569341e-3*KV). The unit label follows
    the projected-potential interpretation; printed literal rounding remains
    unverified. Wavelength is derived through the existing public helper.
    """
    voltage = _finite(voltage_kv, "voltage_kv")
    if voltage <= 0:
        raise ValueError("voltage_kv must be positive")
    wavelength = _finite(electron_wavelength_angstrom(voltage), "derived wavelength_angstrom")
    if wavelength <= 0:
        raise ValueError("derived wavelength_angstrom must be positive in float64")
    sigma = _finite(0.020886416 * wavelength * (1.0 + 1.9569341e-3 * voltage),
                    "interaction constant")
    if sigma <= 0:
        raise ValueError("interaction constant must be positive in float64")
    return sigma


def phase_grating_from_potential(
    projected_potential: ComplexGrid, *, voltage_kv: float,
) -> ComplexGrid:
    """PGMAP: exp(SIGMA*Im(V))*(cos(-SIGMA*Re(V))+i*sin(-SIGMA*Re(V))).

    Negative imaginary potential attenuates amplitude. Positive imaginary
    potential amplifies it as written; no clipping or sign repair occurs.
    V is already projected: no slice-thickness multiplier is accepted.
    """
    values = _complex_grid(projected_potential, "projected_potential")
    sigma = electron_interaction_constant(voltage_kv)
    try:
        with np.errstate(over="raise", invalid="raise", under="ignore"):
            phase = -sigma * values.real
            amplitude = np.exp(sigma * values.imag)
            result = np.empty(values.shape, dtype=np.complex128)
            result.real = np.cos(phase) * amplitude
            result.imag = np.sin(phase) * amplitude
    except FloatingPointError as exc:
        raise ValueError("phase grating must remain finite in float64") from exc
    return result


def phase_grating_from_structure_factors(
    elastic: ComplexGrid,
    absorptive: ComplexGrid,
    reciprocal_grid: ReciprocalVectorGrid,
    *,
    voltage_kv: float,
    dsap_inverse_angstrom: float,
    method: Literal["direct", "fft"] = "direct",
) -> ComplexGrid:
    """Compose combination, shared DSAP/Nyquist masks, Fourier sum and PGMAP.

    Requires even meshes for the source-backed boundary mask. Supply arrays
    from the existing elastic/absorptive APIs on the same grid. This is a
    Python calculation, not the unverified legacy FFT/refinement workflow.
    method='direct' remains the default oracle; 'fft' explicitly selects
    the NumPy backend. Both methods share masks and exponentiation.
    """
    if not isinstance(method, str) or method not in ("direct", "fft"):
        raise ValueError("method must be 'direct' or 'fft'")
    combined = combine_elastic_absorptive_structure_factors(elastic, absorptive)
    physical_mask = reciprocal_cutoff_mask(reciprocal_grid, dsap_inverse_angstrom)
    boundary_mask = legacy_nyquist_mask(reciprocal_grid)
    if combined.shape != physical_mask.shape:
        raise ValueError("structure-factor shape must match reciprocal_grid")
    selected = np.where(physical_mask & boundary_mask, combined, 0.0j)
    reconstruct = projected_potential_direct if method == "direct" else projected_potential_fft
    potential = reconstruct(selected, reciprocal_grid)
    return phase_grating_from_potential(potential, voltage_kv=voltage_kv)
