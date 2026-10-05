"""Detector spatial-response metadata, not complete detector calibration.

No calibration estimation, file loading, normalization, or correction is performed.
"""
from dataclasses import dataclass
from numbers import Integral, Real
import math

import numpy as np
from numpy.typing import NDArray


def _text(value, name):
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f'{name} must be nonempty text')
    return value


def _real(value, name):
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, Real) or not math.isfinite(value):
        raise ValueError(f'{name} must be a finite real number, excluding bool')
    return float(value)


def _integer(value, name):
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, Integral):
        raise TypeError(f'{name} must be an integer, excluding bool')
    return int(value)


def _tuple(value, name):
    if not isinstance(value, (tuple, list)):
        raise TypeError(f'{name} must be a tuple or list')
    return tuple(value)


def _provenance(value):
    return tuple(_text(v, 'provenance entry') for v in _tuple(value, 'provenance'))


def _array(value, dtype, name):
    """One defensive copy, backed by immutable bytes; no dtype coercion."""
    if not isinstance(value, np.ndarray) or value.dtype != np.dtype(dtype):
        raise TypeError(f'{name} must be a {np.dtype(dtype)} NumPy array')
    if value.ndim != 2 or not all(value.shape):
        raise ValueError(f'{name} must be a nonempty two-dimensional array')
    if not np.all(np.isfinite(value)):
        raise ValueError(f'{name} must contain finite values')
    return np.frombuffer(value.tobytes(order='C'), dtype=value.dtype).reshape(value.shape)


@dataclass(frozen=True)
class PSFCalibrationMetrics:
    """Optional diagnostics; no acceptance threshold is implied.

    goodness_of_fit uses the calibration method's declared units (legacy: intensity
    squared); poisson_mismatch is dimensionless. Counts refer to accepted solutions.
    """
    goodness_of_fit: float | None = None
    poisson_mismatch: float | None = None
    accepted_measurements: int | None = None

    def __post_init__(self):
        for name in ('goodness_of_fit', 'poisson_mismatch'):
            value = getattr(self, name)
            if value is not None:
                value = _real(value, name)
                if value < 0:
                    raise ValueError(f'{name} must be nonnegative')
                object.__setattr__(self, name, value)
        if self.accepted_measurements is not None:
            value = _integer(self.accepted_measurements, 'accepted_measurements')
            if value < 0:
                raise ValueError('accepted_measurements must be nonnegative')
            object.__setattr__(self, 'accepted_measurements', value)


@dataclass(frozen=True, eq=False)
class DetectorResponseModel:
    """Optional PSF/OTF/MTF arrays or an opaque artifact reference.

    Arrays share a declared (row, column) grid: PSF/MTF float64, OTF complex128.
    Values, signs, centering and normalization are not modified. Artifact references
    are stored verbatim, never resolved or read; provenance must describe sampling,
    centering and normalization when known. Presence does not imply application to
    any experimental frame. Array-bearing objects use identity equality.
    """
    calibration_method: str
    psf: NDArray[np.float64] | None = None
    otf: NDArray[np.complex128] | None = None
    mtf: NDArray[np.float64] | None = None
    artifact_reference: str | None = None
    metrics: PSFCalibrationMetrics | None = None
    provenance: tuple[str, ...] = ()

    def __post_init__(self):
        _text(self.calibration_method, 'calibration_method')
        if self.artifact_reference is not None:
            _text(self.artifact_reference, 'artifact_reference')
        if self.metrics is not None and not isinstance(self.metrics, PSFCalibrationMetrics):
            raise TypeError('metrics must be PSFCalibrationMetrics')
        shapes = set()
        for name, dtype in (('psf', np.float64), ('otf', np.complex128), ('mtf', np.float64)):
            value = getattr(self, name)
            if value is not None:
                value = _array(value, dtype, name)
                shapes.add(value.shape)
                object.__setattr__(self, name, value)
        if len(shapes) > 1:
            raise ValueError('supplied PSF/OTF/MTF arrays must share one grid shape')
        object.__setattr__(self, 'provenance', _provenance(self.provenance))
