"""Explicit uncertainty units: sigma is standard deviation, variance is sigma².

Saved text sigma maps are a compatibility/replay frontend, not a production
requirement. Neither masks nor statistical floors are inferred here.
"""
from pathlib import Path
import numpy as np


def _sigma_map(value):
    if not isinstance(value,np.ndarray) or value.dtype.kind not in 'iuf':
        raise TypeError('sigma_map must be a real numerical NumPy array')
    if value.ndim != 2 or not all(value.shape):
        raise ValueError('sigma_map must be nonempty and two-dimensional (row, column)')
    try:
        with np.errstate(over='raise',invalid='raise'):
            sigma = np.array(value,dtype=np.float64,copy=True)
    except FloatingPointError as exc:
        raise ValueError('Non-finite sigma map') from exc
    if not np.isfinite(sigma).all():
        raise ValueError('sigma_map must contain finite values')
    if np.any(sigma < 0):
        raise ValueError('sigma_map must be nonnegative; no absolute-value correction')
    return np.frombuffer(sigma.tobytes(),dtype=np.float64).reshape(sigma.shape)


def variance_map_from_sigma_map(sigma_map):
    """Owned immutable float64 variance, same shape, exactly sigma**2.

    Reject negative/nonfinite sigma and overflow; retain zero without an epsilon.
    """
    sigma = _sigma_map(sigma_map)
    try:
        with np.errstate(over='raise',invalid='raise'):
            variance = sigma**2
    except FloatingPointError as exc:
        raise ValueError('Non-finite variance map from sigma squared') from exc
    return np.frombuffer(variance.tobytes(),dtype=np.float64).reshape(variance.shape)


def load_legacy_sigma_map(path, *, expected_shape=None):
    """Load a whitespace-delimited real 2-D saved sigma image, without squaring.

    Shape is dynamic; callers may require a known corresponding pattern shape.
    Does not infer whether an arbitrary file historically stored sigma or variance.
    """
    sigma = _sigma_map(np.loadtxt(Path(path),dtype=np.float64,ndmin=2))
    if expected_shape is not None and sigma.shape != tuple(expected_shape):
        raise ValueError(f'sigma shape {sigma.shape} does not match expected shape {tuple(expected_shape)}')
    return sigma
