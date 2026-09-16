"""Normalized weighted root error from amoeba.f:772-782."""

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from .experimental import ExperimentalIntensityData
from .intensity_scaling import (
    _scalar, apply_legacy_intensity_scaling, compute_legacy_intensity_scaling,
)


@dataclass(frozen=True)
class LegacyObjectiveResult:
    """Float64 scalar result and the two source accumulation diagnostics.

    Arithmetic-generated NaN/Inf are retained, including in the diagnostics.
    """

    objective: np.float64
    scaling: np.float64
    weighted_squared_error: np.float64
    weighted_experimental_norm: np.float64

    def __post_init__(self):
        for name in self.__dataclass_fields__:
            object.__setattr__(self, name, _scalar(getattr(self, name), name, finite=False))


def _point_contributions(scaled, experimental):
    """Return numerator/normalizer terms in (tilt, beam) order, without sums."""
    x, v, w = experimental.intensities, experimental.variances, experimental.weights
    with np.errstate(all="ignore"):
        term = scaled-x
        return w*term*term/v, w*x*x/v


def compute_legacy_objective(
    calculated: NDArray[np.float64], experimental: ExperimentalIntensityData, *, iref: float,
) -> LegacyObjectiveResult:
    """Compute sqrt(sum(w*(s*c-x)*(s*c-x)/v) / sum(w*x*x/v)).

    Existing YNORM functions validate finite float64 (tilt, beam) inputs and
    compute/apply the global scale s. Point products are left-associated;
    accumulation is beam outermost, tilt innermost. No masking or singular
    arithmetic repair occurs. GRADIENT preprocessing is not implemented.
    """
    scaling = compute_legacy_intensity_scaling(calculated, experimental, iref=iref)
    scaled = apply_legacy_intensity_scaling(calculated, scaling)
    errors, norms = _point_contributions(scaled, experimental)
    error, norm = np.float64(0), np.float64(0)
    with np.errstate(all="ignore"):
        for beam in range(calculated.shape[1]):
            for tilt in range(calculated.shape[0]):
                norm = norm+norms[tilt, beam]
                error = error+errors[tilt, beam]
        objective = np.sqrt(error/norm)
    return LegacyObjectiveResult(objective, scaling, error, norm)
