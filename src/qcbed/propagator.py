"""Physical PROPG kernel, separate from Fourier normalization and propagation."""

from dataclasses import dataclass
from typing import Literal

import numpy as np
from numpy.typing import NDArray

from .fft_backend import centered_to_unshifted, unshifted_to_centered
from .phase_grating import legacy_nyquist_mask
from .reciprocal import ReciprocalVectorGrid
from .structure_factors import _finite


_SOURCE_TWOPI = 6.2831853


@dataclass(frozen=True)
class PropagationParameters:
    """Internal PROPG parameters; tilt_h/k are reciprocal-basis coefficients.

    Wavelength and undilated slice thickness are in angstroms. Tilts are
    dimensionless, after the caller's scan-basis transformation and dilation
    handling; no radians/mrad conversion is implied. Shifts are fractional
    direct-cell translations. Zero thickness is allowed for kernel validation.
    """

    wavelength_angstrom: float
    slice_thickness_angstrom: float
    tilt_h_fractional: float = 0.0
    tilt_k_fractional: float = 0.0
    x_shift_fractional: float = 0.0
    y_shift_fractional: float = 0.0

    def __post_init__(self) -> None:
        for name in ("wavelength_angstrom", "slice_thickness_angstrom", "tilt_h_fractional",
                     "tilt_k_fractional", "x_shift_fractional", "y_shift_fractional"):
            object.__setattr__(self, name, _finite(getattr(self, name), name))
        if self.wavelength_angstrom <= 0:
            raise ValueError("wavelength_angstrom must be positive")
        if self.slice_thickness_angstrom < 0:
            raise ValueError("slice_thickness_angstrom must be nonnegative")


def _grid(value: ReciprocalVectorGrid) -> ReciprocalVectorGrid:
    if not isinstance(value, ReciprocalVectorGrid):
        raise TypeError("reciprocal_grid must be ReciprocalVectorGrid")
    return value


def propagation_phase_radians(
    reciprocal_grid: ReciprocalVectorGrid, parameters: PropagationParameters,
) -> NDArray[np.float64]:
    """Unmasked propagation-only phase -6.2831853*SG*DZ*d in grid order.

    SG=-lambda/2*(|g|^2+2*t.g), t=(tilt_h*a_star+tilt_k*b_star)/d.
    Uses the explicit PROPG metric grouping. No shift, RMT, aperture or FFT.
    """
    grid = _grid(reciprocal_grid)
    if not isinstance(parameters, PropagationParameters):
        raise TypeError("parameters must be PropagationParameters")
    h, k = grid.indices.h, grid.indices.k
    a, b = grid.lattice.a_star_inverse_angstrom, grid.lattice.b_star_inverse_angstrom
    d = grid.dilation
    try:
        with np.errstate(over="raise", invalid="raise", divide="raise"):
            a2, b2 = a*a/(d*d), b*b/(d*d)
            adotb = a*b*grid.lattice.cos_gamma_star/(d*d)
            a1 = 2.0*parameters.tilt_h_fractional + h
            b1 = 2.0*parameters.tilt_k_fractional + k
            sg = -0.5*parameters.wavelength_angstrom*(a1*h*a2+b1*k*b2+(a1*k+b1*h)*adotb)
            phase = -_SOURCE_TWOPI * (sg*parameters.slice_thickness_angstrom*d)
    except FloatingPointError as exc:
        raise ValueError("propagation phase must remain finite in float64") from exc
    if not np.all(np.isfinite(phase)):
        raise ValueError("propagation phase must remain finite in float64")
    return np.asarray(phase, dtype=np.float64, order="C")


def fractional_slice_shift_phase_radians(
    reciprocal_grid: ReciprocalVectorGrid, *, x_shift_fractional: float = 0.0,
    y_shift_fractional: float = 0.0,
) -> NDArray[np.float64]:
    """Unmasked -6.2831853*(H*X0+K*Y0), in the supplied grid ordering.

    No thickness/dilation scaling or modulo reduction of shifts is applied.
    The source's printed TWOPI makes full-cell periodicity approximate.
    """
    grid = _grid(reciprocal_grid)
    x = _finite(x_shift_fractional, "x_shift_fractional")
    y = _finite(y_shift_fractional, "y_shift_fractional")
    try:
        with np.errstate(over="raise", invalid="raise"):
            return -_SOURCE_TWOPI*(grid.indices.h*x+grid.indices.k*y)
    except FloatingPointError as exc:
        raise ValueError("slice-shift phase must remain finite in float64") from exc


def propagation_boundary_mask(reciprocal_grid: ReciprocalVectorGrid) -> NDArray[np.bool_]:
    """Source PROPG boundary: suppress negative Nyquist axes on even grids.

    Same frequencies as the existing phase-grating Nyquist mask; no DSAP
    cutoff. Mask is returned in input-grid order; odd meshes are rejected.
    """
    return legacy_nyquist_mask(_grid(reciprocal_grid))


def propagation_kernel(
    reciprocal_grid: ReciprocalVectorGrid, parameters: PropagationParameters, *,
    ordering: Literal["centered", "unshifted"],
    boundary: Literal["legacy", "none"] = "legacy",
) -> NDArray[np.complex128]:
    """Return exp(i*(propagation phase+shift phase)), optionally boundary-masked.

    Explicit output ordering is independent of input-grid ordering. 'legacy'
    suppresses the even-mesh negative Nyquist axes; 'none' evaluates the generic
    physical expression on every bin and supports odd meshes. No RMT or wave
    multiplication is included. Legacy PRARRAY is in unshifted order.
    """
    grid = _grid(reciprocal_grid)
    if not isinstance(ordering, str) or ordering not in ("centered", "unshifted"):
        raise ValueError("ordering must be 'centered' or 'unshifted'")
    if not isinstance(boundary, str) or boundary not in ("legacy", "none"):
        raise ValueError("boundary must be 'legacy' or 'none'")
    mask = propagation_boundary_mask(grid) if boundary == "legacy" else None
    physical = propagation_phase_radians(grid, parameters)
    shift = fractional_slice_shift_phase_radians(grid, x_shift_fractional=parameters.x_shift_fractional,
                                                 y_shift_fractional=parameters.y_shift_fractional)
    try:
        with np.errstate(over="raise", invalid="raise"):
            phase = physical + shift
            kernel = np.empty(phase.shape, dtype=np.complex128)
            kernel.real, kernel.imag = np.cos(phase), np.sin(phase)
    except FloatingPointError as exc:
        raise ValueError("combined propagation phase must remain finite in float64") from exc
    if mask is not None:
        kernel[~mask] = 0.0j
    if grid.indices.centered and ordering == "unshifted":
        return centered_to_unshifted(kernel)
    if not grid.indices.centered and ordering == "centered":
        return unshifted_to_centered(kernel)
    return kernel
