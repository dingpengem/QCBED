"""Source FAATM and raw absorptive sums, without refinement or propagation."""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
import math
from numbers import Integral

import numpy as np
from numpy.typing import NDArray

from .legacy_io import AtomSite, Species
from .reciprocal import ReciprocalVectorGrid
from .structure_factors import _finite, _ATOMIC_PHASE_TWOPI


@dataclass(frozen=True)
class AbsorptionCoefficients:
    """Applied INELA/B/C weights, in species order when supplied as a sequence.

    Signed finite weights are permitted by the unrefined source input. No
    refinement masks, coefficient tying, or positivity projection is applied.
    """

    inela: float
    inelb: float
    inelc: float

    def __post_init__(self) -> None:
        for name in ("inela", "inelb", "inelc"):
            object.__setattr__(self, name, _finite(getattr(self, name), name))


def _atomic_domain(atomic_number: int, stl: float) -> tuple[int, float]:
    if isinstance(atomic_number, bool) or not isinstance(atomic_number, Integral):
        raise TypeError("atomic_number must be an integer, not bool")
    if not 1 <= atomic_number <= 98:
        raise ValueError("atomic_number must be in 1..98")
    stl = _finite(stl, "stl_inverse_angstrom", nonnegative=True)
    limit = 1.5 if atomic_number == 1 else 6.0
    if stl > limit:
        raise ValueError(f"stl_inverse_angstrom must be in 0..{limit}")
    return int(atomic_number), stl


def _positive(value: float, name: str) -> float:
    value = _finite(value, name)
    if value <= 0:
        raise ValueError(f"{name} must be positive; zero is unsupported")
    return value


def absorptive_atomic_scattering_factor(
    atomic_number: int,
    stl_inverse_angstrom: float,
    *,
    voltage_kv: float,
    applied_debye_waller_angstrom_squared: float,
    coefficients: AbsorptionCoefficients,
) -> float:
    """Return source FAATM as a binary64 Python float (not a scattering length).

    No area, occupancy, external exp(-B*STL^2), or phase is included. The
    domain is the joint SFMAP elastic/absorptive domain Z=1..98, STL=0..6
    inverse angstroms (H:0..1.5), positive applied B and voltage in kV.
    No empirical validity bounds for this absorption fit are given in source.
    """
    z, s = _atomic_domain(atomic_number, stl_inverse_angstrom)
    voltage = _positive(voltage_kv, "voltage_kv")
    volts = _positive(voltage * 1000.0, "voltage in volts")
    b = _positive(applied_debye_waller_angstrom_squared,
                  "applied_debye_waller_angstrom_squared")
    if not isinstance(coefficients, AbsorptionCoefficients):
        raise TypeError("coefficients must be AbsorptionCoefficients")
    try:
        energy_exponent = -0.0041225 * math.sqrt(volts)
        # Same expression as exp(t)/(1-exp(t)), with a stable denominator.
        energy = 0.5 + math.exp(energy_exponent) / (-math.expm1(energy_exponent))
        power = 0.64579 / (1.0 + 0.079647 * z**0.40196)
        b_power = b**power
        if b_power == 0.0:
            raise ValueError("applied Debye-Waller power underflows in float64")
        prefactor = (0.0012822 * energy * (1.0 - math.exp(-0.86245 * z**0.377))
                     * z**1.6191 * b_power)
        s2 = s * s
        mean_term = coefficients.inela * math.exp(-1.0e12 * s2)
        broad_term = coefficients.inelb * math.exp(-(2.7964 * s**1.83)
                                                  / (1.0 + 0.18088 * z**0.48919))
        tail_term = (coefficients.inelc
                     * ((0.27183 * z**0.42387 + 0.45264 * z**-0.2099) * s2)
                     / ((1.0 + 0.60054 * z**0.63377 + 3.5729 * s**4) * b_power))
        value = prefactor * (mean_term + broad_term + tail_term)
    except (OverflowError, ZeroDivisionError) as exc:
        raise ValueError("absorptive atomic factor is not representable in float64") from exc
    return _finite(value, "absorptive atomic factor")


def absorptive_structure_factors(
    species: Sequence[Species],
    atoms: Sequence[AtomSite],
    reciprocal_grid: ReciprocalVectorGrid,
    *,
    voltage_kv: float,
    applied_debye_waller_angstrom_squared: Sequence[float],
    coefficients: Sequence[AbsorptionCoefficients],
    vacancy_fraction: float,
    absorption_raw: str,
) -> NDArray[np.complex128]:
    """Return AFRE+i*AFIM, including the source's minus sign and area product.

    A = -sum exp(+i*TWOPI*(H*x+K*y))*area*FAATM*exp(-B*STL^2)*OCC.
    Absorption is enabled exactly when absorption_raw[:2]=='YE'. All inputs
    are validated even when disabled. Applied coefficients and B are required
    per species; the parser does not yet retain these parameters.
    """
    if not isinstance(reciprocal_grid, ReciprocalVectorGrid):
        raise TypeError("reciprocal_grid must be ReciprocalVectorGrid")
    if not isinstance(absorption_raw, str):
        raise TypeError("absorption_raw must be str")
    species, atoms = tuple(species), tuple(atoms)
    coefficients = tuple(coefficients)
    thermal = tuple(applied_debye_waller_angstrom_squared)
    if not species:
        raise ValueError("species must not be empty")
    if len(thermal) != len(species) or len(coefficients) != len(species):
        raise ValueError("coefficients and applied_debye_waller_angstrom_squared require one value per species")
    voltage = _positive(voltage_kv, "voltage_kv")
    thermal = tuple(_positive(b, f"applied_debye_waller_angstrom_squared[{i}]")
                    for i, b in enumerate(thermal))
    vacancy = _finite(vacancy_fraction, "vacancy_fraction")
    if not 0 <= vacancy <= 1:
        raise ValueError("vacancy_fraction must be in [0, 1]")
    max_stl = float(np.max(reciprocal_grid.stl_inverse_angstrom))
    for i, entry in enumerate(species):
        if not isinstance(entry, Species):
            raise TypeError(f"species[{i}] must be Species")
        absorptive_atomic_scattering_factor(entry.atomic_number, max_stl, voltage_kv=voltage,
            applied_debye_waller_angstrom_squared=thermal[i], coefficients=coefficients[i])
    sites = []
    for i, atom in enumerate(atoms):
        if not isinstance(atom, AtomSite):
            raise TypeError(f"atoms[{i}] must be AtomSite")
        if atom.record_type != "S":
            raise ValueError(f"atoms[{i}].record_type must be S; unresolved G records are unsupported")
        index = atom.species_index_zero_based(len(species))
        x = _finite(atom.x_fractional, f"atoms[{i}].x_fractional")
        y = _finite(atom.y_fractional, f"atoms[{i}].y_fractional")
        occupancy = _finite(atom.starting_occupancy, f"atoms[{i}].starting_occupancy",
                            nonnegative=True) * (1.0 - vacancy)
        sites.append((index, x, y, occupancy))
    area = _positive(reciprocal_grid.lattice.direct.area_angstrom_squared
                     * reciprocal_grid.dilation * reciprocal_grid.dilation, "dilated cell area")
    result = np.zeros(reciprocal_grid.indices.h.shape, dtype=np.complex128)
    if absorption_raw[:2] != "YE":
        return result
    for row, col in np.ndindex(result.shape):
        stl = float(reciprocal_grid.stl_inverse_angstrom[row, col])
        h = float(reciprocal_grid.indices.h[row, col])
        k = float(reciprocal_grid.indices.k[row, col])
        factors = [absorptive_atomic_scattering_factor(entry.atomic_number, stl,
            voltage_kv=voltage, applied_debye_waller_angstrom_squared=thermal[i],
            coefficients=coefficients[i]) for i, entry in enumerate(species)]
        real = imag = 0.0
        for index, x, y, occupancy in sites:
            phase = _finite(_ATOMIC_PHASE_TWOPI * (h*x + k*y), "atomic phase")
            damping = math.exp(-_finite(thermal[index] * (stl*stl), "thermal exponent"))
            real -= math.cos(phase) * area * factors[index] * damping * occupancy
            imag -= math.sin(phase) * area * factors[index] * damping * occupancy
        result[row, col] = complex(_finite(real, "absorptive real sum"),
                                   _finite(imag, "absorptive imaginary sum"))
    return result


def combine_elastic_absorptive_structure_factors(
    elastic: NDArray[np.complex128], absorptive: NDArray[np.complex128],
) -> NDArray[np.complex128]:
    """Return PMAP=(SFRE-AFIM)+i*(SFIM+AFRE), with no FFT/exponentiation.

    Inputs must be finite complex128 NumPy arrays with identical nonempty
    2D shapes and the same reciprocal ordering (the caller's responsibility).
    """
    for name, array in (("elastic", elastic), ("absorptive", absorptive)):
        if not isinstance(array, np.ndarray) or array.dtype != np.complex128:
            raise TypeError(f"{name} must be a complex128 NumPy array")
        if array.ndim != 2 or 0 in array.shape:
            raise ValueError(f"{name} must have a nonempty 2D shape (mesh_y, mesh_x)")
        if not np.all(np.isfinite(array)):
            raise ValueError(f"{name} must be finite")
    if elastic.shape != absorptive.shape:
        raise ValueError("elastic and absorptive must have identical shapes")
    result = np.empty(elastic.shape, dtype=np.complex128)
    try:
        with np.errstate(over="raise", invalid="raise"):
            result.real = elastic.real - absorptive.imag
            result.imag = elastic.imag + absorptive.real
    except FloatingPointError as exc:
        raise ValueError("combined structure factors must remain finite in float64") from exc
    return result
