from dataclasses import FrozenInstanceError
from pathlib import Path

import numpy as np
import pytest

from qcbed import (ExperimentalIntensityData, LegacyObjectiveResult,
                   compute_legacy_objective as compute, parse_experimental_intensities,
                   compute_legacy_intensity_scaling, apply_legacy_intensity_scaling)
from qcbed.objective import _point_contributions


def data(x, v=None, w=None):
    x = np.asarray(x, dtype=np.float64)
    return ExperimentalIntensityData(x.shape[0], 1, 0,
        tuple((b, 0) for b in range(x.shape[1])), x,
        np.ones_like(x) if v is None else np.asarray(v, dtype=np.float64),
        np.ones_like(x) if w is None else np.asarray(w, dtype=np.float64))


def oracle(c, d, scale):
    # Literal flat-M traversal, independently of the production point helper.
    arrays = [a.flatten(order="F") for a in (c, d.intensities, d.variances, d.weights)]
    norm, error = np.float64(0), np.float64(0)
    with np.errstate(all="ignore"):
        for theory, observed, variance, weight in zip(*arrays):
            term = scale*theory-observed
            norm = norm+weight*observed*observed/variance
            error = error+weight*term*term/variance
        return np.sqrt(error/norm), error, norm


@pytest.mark.parametrize("mode,scale,expected", [(1, 1, 1.5), (0, -6/13, np.sqrt(56/65)),
                                                (-1, 10/13, np.sqrt(24/13))])
def test_hand_case(mode, scale, expected):
    c = np.array([[1., 2.], [3., 4.]])
    d = data([[2, -4], [6, 0]], [[1, 2], [3, 4]], [[1, 2], [0, 2]])
    # Error(s)=13*s*s+12*s+20; normalizer=20.
    result = compute(c, d, iref=mode)
    assert result.scaling == scale
    assert result.objective == pytest.approx(expected, rel=2e-15)
    assert result.weighted_experimental_norm == 20
    if mode == 1:
        assert result.weighted_squared_error == 45


@pytest.mark.parametrize("shape", [(4, 1), (1, 4), (3, 2)])
@pytest.mark.parametrize("mode", [0, -1, 1.25])
def test_oracle_prescaling_shape_dtype_and_immutability(shape, mode):
    c = np.arange(np.prod(shape), dtype=np.float64).reshape(shape)
    d = data(3*c-7, c+.5, c/3)
    originals = [a.copy() for a in (c, d.intensities, d.variances, d.weights)]
    result = compute(c, d, iref=mode)
    scale = compute_legacy_intensity_scaling(c, d, iref=mode)
    expected = oracle(c, d, scale)
    assert (result.objective, result.weighted_squared_error, result.weighted_experimental_norm) == expected
    scaled = apply_legacy_intensity_scaling(c, scale)
    assert oracle(scaled, d, np.float64(1)) == expected
    points, norms = _point_contributions(scaled, d)
    assert points.shape == norms.shape == shape
    assert points.dtype == norms.dtype == np.float64
    for name in result.__dataclass_fields__:
        assert isinstance(getattr(result, name), np.float64)
    with pytest.raises(FrozenInstanceError):
        result.objective = 0
    for a, before in zip((c, d.intensities, d.variances, d.weights), originals):
        np.testing.assert_array_equal(a, before)


def test_beam_outer_accumulation():
    d = data(np.ones((2, 2)), w=[[1e16, -1e16], [1, 1]])
    result = compute(np.zeros((2, 2)), d, iref=1)
    # Beam-major (((1e16+1)-1e16)+1)=1; row-major gives 2.
    assert result.weighted_squared_error == result.weighted_experimental_norm == 1


@pytest.mark.parametrize("c,x,v,w,expected", [
    (0, 2, 1, 1, 1), (-2, -2, 1, 1, 0), (1, 0, 1, 1, np.inf),
    (0, 0, 1, 1, np.nan), (1, 2, 0, 1, np.nan), (1, 2, 0, 0, np.nan),
    (1, 2, 1, 0, np.nan), (1, 2, -1, 1, .5),
    (1e308, 1, 1, 1, np.inf), (1, 1e308, 1, 1, np.nan),
])
def test_unguarded_arithmetic(c, x, v, w, expected):
    with np.errstate(all="raise"):
        result = compute(np.array([[float(c)]]), data([[x]], [[v]], [[w]]), iref=1)
    np.testing.assert_equal(result.objective, expected)


def test_negative_ratio_and_multiplication_order():
    result = compute(np.array([[2., 0.]]), data([[1, 1]], w=[[-1, 2]]), iref=1)
    assert result.objective == 1
    result = compute(np.array([[3., 1.]]), data([[1, 1]], w=[[-1, 2]]), iref=1)
    assert np.isnan(result.objective)  # sqrt(-4/1), no clipping
    result = compute(np.array([[1e200]]), data([[1]], w=[[1e-200]]), iref=1)
    assert np.isfinite(result.weighted_squared_error)  # w*t*t, not w*(t*t)
    assert result.weighted_squared_error == 1e200


@pytest.mark.parametrize("mode", [0, -1])
def test_nonfinite_scaling_propagates(mode):
    result = compute(np.zeros((1, 1)), data([[1]]), iref=mode)
    assert not np.isfinite(result.scaling)
    assert np.isnan(result.objective)


@pytest.mark.parametrize("bad", [np.ones((1, 1), dtype=np.float32), np.ones(1),
    np.ones((2, 1)), np.array([[np.nan]]), np.array([[np.inf]]), [[1.]], True])
def test_invalid_calculated(bad):
    with pytest.raises((TypeError, ValueError), match="calculated"):
        compute(bad, data([[1]]), iref=1)


@pytest.mark.parametrize("bad", [True, "1", np.nan, np.inf, -np.inf])
def test_invalid_iref(bad):
    with pytest.raises((TypeError, ValueError), match="iref"):
        compute(np.ones((1, 1)), data([[1]]), iref=bad)


def test_invalid_experimental_and_result():
    with pytest.raises(TypeError, match="experimental"):
        compute(np.ones((1, 1)), None, iref=1)
    with pytest.raises(TypeError, match="objective"):
        LegacyObjectiveResult(True, 1, 1, 1)


@pytest.mark.parametrize("filename", ["in.dat", "startingin.dat", "bestin.dat"])
@pytest.mark.parametrize("mode", [0, -1, 1.25])
def test_reference_files(filename, mode):
    beams = ((0, 0), (1, 1), (1, -1), (-1, -1), (-1, 1), (2, 0), (0, -2), (-2, 0), (0, 2))
    d = parse_experimental_intensities(Path(__file__).parents[1]/"reference_data"/filename,
                                      beam_indices_hk=beams)
    c = np.ones((1521, 9), dtype=np.float64)
    result = compute(c, d, iref=mode)
    assert np.count_nonzero(d.variances == 0) == 117
    assert np.isnan(result.objective)  # Raw observations, not post-GRADIENT equivalence.
    expected = oracle(c, d, result.scaling)
    np.testing.assert_equal([result.objective, result.weighted_squared_error,
                             result.weighted_experimental_norm], expected)
    if mode < 0:
        assert result.scaling == pytest.approx(12975.673278983971, abs=1e-9, rel=0)
