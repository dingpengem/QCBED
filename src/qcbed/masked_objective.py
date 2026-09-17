"""Explicit non-legacy objective policy; no changes to legacy arithmetic."""

from dataclasses import dataclass

import numpy as np

from .experimental import ExperimentalIntensityData
from .intensity_scaling import _scalar, _values
from .objective import compute_legacy_objective


@dataclass(frozen=True)
class MaskedObjectiveResult:
    objective: np.float64
    scaling: np.float64
    weighted_squared_error: np.float64
    weighted_experimental_norm: np.float64


def compute_masked_objective(calculated, experimental, *, iref):
    """An explicit non-legacy numerical policy for datasets whose zero-weight
    masked points contain zero variance.

    Select weight != 0 before any variance division in both YNORM and error.
    Keep weights, IREF branches, and beam-major accumulation on selected points.
    Valid points require finite positive variance. Empty masks, singular scales
    and undefined normalized errors raise instead of repairing input values.
    No historical Fortran executable equivalence is claimed.
    """
    c = _values(calculated, "calculated")
    if not isinstance(experimental, ExperimentalIntensityData):
        raise TypeError("experimental must be ExperimentalIntensityData")
    x, v, w = experimental.intensities, experimental.variances, experimental.weights
    if c.shape != x.shape:
        raise ValueError("calculated and experimental arrays must have matching (tilt, beam) shape")
    mode = _scalar(iref, "iref", finite=True)
    valid = w != 0
    bad = valid & ((v <= 0) | ~np.isfinite(v))
    if np.any(bad):
        t, b = np.argwhere(bad)[0]
        raise ValueError(f"masked objective: nonzero-weight point requires finite positive variance; "
                         f"tilt {t}, beam {b} {experimental.beam_indices_hk[b]}, "
                         f"variance={v[t,b]}, weight={w[t,b]}")
    # Selecting in beam-major order also prevents arithmetic on masked values.
    selected = valid.T.ravel()
    if not np.any(selected):
        raise ValueError("masked objective requires at least one nonzero-weight point")
    cv, xv, vv, wv = (a.T.ravel()[selected] for a in (c, x, v, w))
    with np.errstate(all="ignore"):
        if mode > 0:
            scaling = mode
        else:
            numerator, denominator = np.float64(0), np.float64(0)
            for theory, observed, variance, weight in zip(cv, xv, vv, wv):
                if mode == 0:
                    numerator += weight*theory*observed/variance
                    denominator += weight*theory*theory/variance
                else:
                    numerator += weight*abs(observed)
                    denominator += weight*abs(theory)
            if not np.isfinite(denominator) or denominator == 0:
                raise ValueError("masked objective: YNORM denominator must be finite and nonzero")
            scaling = numerator/denominator
        if not np.isfinite(scaling):
            raise ValueError("masked objective: YNORM must be finite")
        error, norm = np.float64(0), np.float64(0)
        for theory, observed, variance, weight in zip(cv, xv, vv, wv):
            term = scaling*theory-observed
            error += weight*term*term/variance
            norm += weight*observed*observed/variance
        if not np.isfinite(norm) or norm <= 0:
            raise ValueError("masked objective: experimental normalizer must be finite and positive")
        if not np.isfinite(error) or error < 0:
            raise ValueError("masked objective: weighted squared error must be finite and nonnegative")
        objective = np.sqrt(error/norm)
        if not np.isfinite(objective):
            raise ValueError("masked objective must be finite")
    return MaskedObjectiveResult(objective, scaling, error, norm)


def compute_selected_objective(calculated, experimental, *, iref, objective_mode="legacy"):
    """Explicit dispatcher; omitted mode retains the untouched legacy policy."""
    if objective_mode == "legacy":
        return compute_legacy_objective(calculated, experimental, iref=iref)
    if objective_mode == "masked":
        return compute_masked_objective(calculated, experimental, iref=iref)
    raise ValueError("objective_mode must be 'legacy' or 'masked'")
