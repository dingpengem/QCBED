"""Raw elastic atomic sum from SFMAP/SFMAP2, before map operations."""

from __future__ import annotations

from collections.abc import Sequence
import math
from numbers import Real

import numpy as np
from numpy.typing import NDArray

from .legacy_io import AtomSite, Species
from .reciprocal import ReciprocalVectorGrid
from .scattering import elastic_atomic_scattering_factor


_ATOMIC_PHASE_TWOPI = 6.283185307


def _finite(value: float, name: str, *, nonnegative: bool = False) -> float:
    if isinstance(value, bool) or not isinstance(value, Real):
        raise TypeError(f"{name} must be a real scalar, not bool")
    try:
        result = float(value)
    except OverflowError as exc:
        raise ValueError(f"{name} must be finite") from exc
    if not math.isfinite(result):
        raise ValueError(f"{name} must be finite")
    if nonnegative and result < 0:
        raise ValueError(f"{name} must be nonnegative")
    return result


def elastic_structure_factors(
    species: Sequence[Species],
    atoms: Sequence[AtomSite],
    reciprocal_grid: ReciprocalVectorGrid,
    *,
    applied_debye_waller_angstrom_squared: Sequence[float],
    vacancy_fraction: float,
) -> NDArray[np.complex128]:
    """Evaluate the unmasked elastic sum in source atom order.

    Returns a C-contiguous complex128 array (mesh_y, mesh_x), in the supplied
    grid ordering. Each term is fe(Z, STL, option=2)/dilated_area times
    exp(+i*6.283185307*(H*x+K*y)), exp(-B*STL**2), and
    starting_occupancy*(1-vacancy_fraction). B is APPLIEDDW, one per species;
    both it and the slice vacancy fraction must be supplied explicitly.

    Only fixed S atoms are supported. No slice-shift phase occurs in this
    sum: it belongs to PROPG. No DSAP/Nyquist mask, refinement multiplier,
    absorption, or propagation is applied. See docs/structure_factors.md.
    """
    if not isinstance(reciprocal_grid, ReciprocalVectorGrid):
        raise TypeError("reciprocal_grid must be ReciprocalVectorGrid")
    species = tuple(species)
    atoms = tuple(atoms)
    if not species:
        raise ValueError("species must not be empty")
    thermal = tuple(applied_debye_waller_angstrom_squared)
    if len(thermal) != len(species):
        raise ValueError("applied_debye_waller_angstrom_squared requires one value per species")
    thermal = tuple(_finite(b, f"applied_debye_waller_angstrom_squared[{i}]",
                            nonnegative=True) for i, b in enumerate(thermal))
    vacancy = _finite(vacancy_fraction, "vacancy_fraction")
    if not 0 <= vacancy <= 1:
        raise ValueError("vacancy_fraction must be in [0, 1]")
    for i, entry in enumerate(species):
        if not isinstance(entry, Species):
            raise TypeError(f"species[{i}] must be Species")
        # Reuse the scattering API's atomic-number validation, including bool.
        elastic_atomic_scattering_factor(entry.atomic_number, 0.0, option=2)
    validated_atoms = []
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
        validated_atoms.append((index, x, y, occupancy))

    area = _finite(reciprocal_grid.lattice.direct.area_angstrom_squared
                   * reciprocal_grid.dilation * reciprocal_grid.dilation,
                   "dilated cell area", nonnegative=True)
    if area == 0:
        raise ValueError("dilated cell area must be positive in float64")
    indices = reciprocal_grid.indices
    result = np.zeros(indices.h.shape, dtype=np.complex128)
    # K outer/H inner corresponds to Fortran's H-fastest PMAP ordering.
    for row, col in np.ndindex(result.shape):
        stl = float(reciprocal_grid.stl_inverse_angstrom[row, col])
        s2 = stl * stl
        factors = [elastic_atomic_scattering_factor(entry.atomic_number, stl, option=2)
                   / area for entry in species]
        h, k = float(indices.h[row, col]), float(indices.k[row, col])
        real = imag = 0.0
        for index, x, y, occupancy in validated_atoms:
            phase = _finite(_ATOMIC_PHASE_TWOPI * (h * x + k * y), "atomic phase")
            thermal_exponent = _finite(thermal[index] * s2, "thermal exponent")
            thermal_factor = math.exp(-thermal_exponent)
            # Separate real/imag additions preserve the visible source order.
            real += factors[index] * math.cos(phase) * thermal_factor * occupancy
            imag += factors[index] * math.sin(phase) * thermal_factor * occupancy
        result[row, col] = complex(_finite(real, "elastic real sum"),
                                   _finite(imag, "elastic imaginary sum"))
    return result
