"""Electron microscope scalar relations from the legacy Fortran source."""

from __future__ import annotations

import math
from numbers import Real


def electron_wavelength_angstrom(voltage_kv: float) -> float:
    """Return relativistic electron wavelength in angstroms for voltage in kV.

    This is the exact statement in ``legacy_fortran/oneitpara.f``:
    ``WAVEL=0.3878314/DSQRT(KV*(1.0+0.97846707E-03*KV))``.
    Uses Python binary64 floats (float64-compatible) and returns a Python float.
    The input is always interpreted as kV; callers must convert volts explicitly.
    Units cannot be inferred from a bare number's magnitude.

    Exact executable equivalence is not established: default-kind Fortran
    literal rounding may depend on the original compiler flags.

    Raises TypeError for booleans or non-real scalars, and ValueError for
    non-finite or non-positive voltage.
    """
    if isinstance(voltage_kv, bool) or not isinstance(voltage_kv, Real):
        raise TypeError("voltage_kv must be a real numeric scalar in kV, not bool")
    try:
        voltage = float(voltage_kv)
    except OverflowError as exc:
        raise ValueError("voltage_kv must be finite and positive in kV") from exc
    if not math.isfinite(voltage) or voltage <= 0.0:
        raise ValueError("voltage_kv must be finite and positive in kV")
    return float(0.3878314 / math.sqrt(voltage * (1.0 + 0.97846707e-3 * voltage)))
