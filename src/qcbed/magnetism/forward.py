"""R0 magnetic phase-grating adapter (periodic longitudinal A_z approximation)."""
from __future__ import annotations

from collections.abc import Sequence

import numpy as np

from .research import e, hbar


def magnetic_phase_grating(projected_az_weber: np.ndarray) -> np.ndarray:
    """Return exp(-i e/hbar integral A_z dz), for integrated A_z in weber."""
    az = np.asarray(projected_az_weber, dtype=np.float64)
    if az.ndim != 2 or 0 in az.shape or not np.all(np.isfinite(az)):
        raise ValueError("projected_az_weber must be a finite, nonempty 2D array")
    return np.exp(-1j * (e / hbar) * az).astype(np.complex128)


def add_magnetic_gratings(
    electrostatic_gratings: Sequence[np.ndarray],
    projected_az_weber: Sequence[np.ndarray],
) -> tuple[np.ndarray, ...]:
    """Combine stored electrostatic gratings and longitudinal magnetic phases."""
    gratings = tuple(electrostatic_gratings)
    fields = tuple(projected_az_weber)
    if not gratings or len(gratings) != len(fields):
        raise ValueError("Stored electrostatic and magnetic slices must match")
    combined = []
    for idx, (grating, field) in enumerate(zip(gratings, fields)):
        grating = np.asarray(grating)
        if (grating.ndim != 2 or grating.dtype != np.complex128
                or not np.all(np.isfinite(grating))):
            raise ValueError(f"electrostatic_gratings[{idx}] must be finite complex128 2D")
        magnetic = magnetic_phase_grating(field)
        if magnetic.shape != grating.shape:
            raise ValueError("magnetic and electrostatic grating grids must match")
        combined.append(grating * magnetic)
    return tuple(combined)


def project_az_midpoint(az_t_m: np.ndarray, dz_m: Sequence[float]) -> tuple[np.ndarray, ...]:
    """Approximate slice integrals of longitudinal vector potential in weber."""
    field = np.asarray(az_t_m, dtype=float)
    thickness = np.asarray(dz_m, dtype=float)
    if (field.ndim != 3 or 0 in field.shape
            or thickness.shape != (field.shape[0],)
            or not np.all(np.isfinite(field))
            or not np.all(np.isfinite(thickness))
            or np.any(thickness <= 0)):
        raise ValueError("Expect finite A_z[s,y,x] and positive slice thicknesses")
    return tuple(field[i] * thickness[i] for i in range(field.shape[0]))
