"""Global YNORM branches from amoeba.f:736-770; no objective evaluation."""

from numbers import Real

import numpy as np
from numpy.typing import NDArray

from .experimental import ExperimentalIntensityData
from .multislice import _snapshot


def _values(values, name):
    if not isinstance(values, np.ndarray) or values.dtype != np.float64:
        raise TypeError(f"{name} must be a float64 NumPy array")
    if values.ndim != 2:
        raise ValueError(f"{name} must have (tilt, beam) shape")
    if not np.all(np.isfinite(values)):
        raise ValueError(f"{name} must be finite")
    return values


def _scalar(value, name, *, finite):
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, Real):
        raise TypeError(f"{name} must be a real scalar excluding bool")
    result = np.float64(value)
    if finite and not np.isfinite(result):
        raise ValueError(f"{name} must be finite")
    return result


def compute_legacy_intensity_scaling(
    calculated: NDArray[np.float64], experimental: ExperimentalIntensityData, *, iref: float,
) -> np.float64:
    """Compute a single YNORM for all beams/tilts, without mutating arrays.

    iref>0: fixed factor; iref==0: sum(w*c*x/v)/sum(w*c*c/v);
    iref<0: sum(w*abs(x))/sum(w*abs(c)). Loops are beam then tilt,
    retaining source multiplication and accumulation order, without masking.
    Finite inputs may produce NaN/Inf under unguarded IEEE float64 arithmetic.
    In particular zero weight does not bypass division by zero variance.
    GRADIENT is not performed: calculated must be the array to be scaled.
    """
    c = _values(calculated, "calculated")
    if not isinstance(experimental, ExperimentalIntensityData):
        raise TypeError("experimental must be ExperimentalIntensityData")
    x, v, w = (_values(getattr(experimental, name), name)
               for name in ("intensities", "variances", "weights"))
    if any(values.shape != c.shape for values in (x, v, w)):
        raise ValueError("calculated and experimental arrays must have exactly matching (tilt, beam) shape")
    mode = _scalar(iref, "iref", finite=True)
    if mode > 0:  # amoeba.f:736-739
        return mode
    numerator, denominator = np.float64(0), np.float64(0)
    with np.errstate(divide="ignore", invalid="ignore", over="ignore", under="ignore"):
        for beam in range(c.shape[1]):
            for tilt in range(c.shape[0]):
                if mode == 0:  # amoeba.f:743-755, TOTALPOINTS in beam-major order
                    numerator = numerator + w[tilt, beam]*c[tilt, beam]*x[tilt, beam]/v[tilt, beam]
                    denominator = denominator + w[tilt, beam]*c[tilt, beam]*c[tilt, beam]/v[tilt, beam]
                else:  # amoeba.f:758-770, beam outermost, tilt innermost
                    numerator = numerator + w[tilt, beam]*abs(x[tilt, beam])
                    denominator = denominator + w[tilt, beam]*abs(c[tilt, beam])
        return np.float64(numerator/denominator)


def apply_legacy_intensity_scaling(
    calculated: NDArray[np.float64], scaling: float,
) -> NDArray[np.float64]:
    """Copy YNORM*THEORINTS (amoeba.f:197,777) into immutable float64 storage.

    Scaling may be NaN/Inf produced by the source's unguarded ratio. Those
    values propagate without clipping or repair; calculated inputs are finite.
    """
    values = _values(calculated, "calculated")
    factor = _scalar(scaling, "scaling", finite=False)
    with np.errstate(over="ignore", invalid="ignore", under="ignore"):
        return _snapshot(factor*values)
