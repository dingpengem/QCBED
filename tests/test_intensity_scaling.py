from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest

from qcbed import (ExperimentalIntensityData, parse_experimental_intensities,
                   compute_legacy_intensity_scaling as compute,
                   apply_legacy_intensity_scaling as apply)


def data(x, variance=None, weight=None):
    x = np.asarray(x, dtype=np.float64)
    return ExperimentalIntensityData(x.shape[0], 1, 0, tuple((i, 0) for i in range(x.shape[1])),
        x, np.ones_like(x) if variance is None else np.asarray(variance, dtype=np.float64),
        np.ones_like(x) if weight is None else np.asarray(weight, dtype=np.float64))


def oracle(c, d, mode):
    if mode > 0:
        return np.float64(mode)
    # Independent one-dimensional translation of legacy M=1..TOTALPOINTS.
    flat = [a.T.flatten() for a in (c, d.intensities, d.variances, d.weights)]
    total1, total2 = np.float64(0), np.float64(0)
    with np.errstate(all="ignore"):
        for theory, observed, variance, weight in zip(*flat):
            if mode == 0:
                term1 = weight*theory
                term1 *= observed
                term1 /= variance
                term2 = weight*theory
                term2 *= theory
                term2 /= variance
            else:
                term1, term2 = weight*np.abs(observed), weight*np.abs(theory)
            total1 += term1
            total2 += term2
        return total1/total2


@pytest.mark.parametrize("mode,expected", [(2.5, 2.5), (0, -6/13), (-1, 10/13)])
def test_hand_calculated_global_factor(mode, expected):
    c = np.array([[1., 2.], [3., 4.]])
    d = data([[2, -4], [6, 0]], [[1, 2], [3, 4]], [[1, 2], [0, 2]])
    # LS sums=-6,13; absolute-count sums=10,13 (zero weight contributes zero).
    factor = compute(c, d, iref=mode)
    assert isinstance(factor, np.float64)
    assert factor == expected
    np.testing.assert_array_equal(apply(c, factor), c*expected)


@pytest.mark.parametrize("shape", [(4, 1), (1, 4), (3, 2)])
@pytest.mark.parametrize("mode", [0, -1, 1.25])
def test_literal_order_oracle_recomputation_and_no_mutation(shape, mode):
    c = np.arange(1, np.prod(shape)+1, dtype=np.float64).reshape(shape)
    d = data(3*c-7, c+.5, c/3)
    before = [v.copy() for v in (c, d.intensities, d.variances, d.weights)]
    factor = compute(c, d, iref=mode)
    assert factor == oracle(c, d, mode)
    second = compute(2*c, d, iref=mode)
    assert second == (factor if mode > 0 else factor/2)
    scaled = apply(c, factor)
    np.testing.assert_array_equal(scaled, c*factor)
    for value, original in zip((c, d.intensities, d.variances, d.weights), before):
        np.testing.assert_array_equal(value, original)
    with pytest.raises(ValueError):
        scaled.setflags(write=True)


def test_signed_ls_no_clipping_and_absolute_count_branch():
    c = np.array([[1.], [2.]])
    d = data([[-2], [-4]])
    assert compute(c, d, iref=0) == -2
    assert compute(c, d, iref=-1) == 2
    assert compute(-c, d, iref=-3) == 2


@pytest.mark.parametrize("mode", [0, -1, 2])
def test_zero_variance_zero_weights_and_ignored_fields(mode):
    c = np.array([[1.], [2.]])
    d = data([[3], [4]], [[0], [1]], [[0], [1]])
    factor = compute(c, d, iref=mode)
    if mode == 0:
        assert np.isnan(factor)  # 0/0 is not skipped despite zero weight.
        assert np.isnan(apply(c, factor)).all()
    else:
        assert factor == 2
    with np.errstate(all="raise"):
        assert np.isnan(compute(c, d, iref=0))  # API makes IEEE policy explicit.


@pytest.mark.parametrize("mode", [0, -1])
def test_zero_denominator_and_zero_experiment(mode):
    d = data([[0], [0]])
    assert compute(np.ones((2, 1)), d, iref=mode) == 0
    assert np.isnan(compute(np.zeros((2, 1)), d, iref=mode))
    assert np.isnan(compute(np.ones((2, 1)), replace(d, weights=np.zeros((2, 1))), iref=mode))
    factor = compute(np.zeros((2, 1)), data([[1], [2]]), iref=mode)
    assert np.isnan(factor) if mode == 0 else np.isposinf(factor)


def test_source_accumulation_order_is_beam_major():
    c = np.ones((2, 2))
    d = data([[1e16, -1e16], [1, 1]])
    # Beam-major accumulation: ((1e16+1)-1e16)+1 = 1, not row-major 2.
    assert compute(c, d, iref=0) == .25


def test_signed_weight_cancellation_and_unsimplified_overflow():
    c = np.array([[1.], [2.]])
    d = data([[1], [0]], weight=[[4], [-1]])
    assert np.isposinf(compute(c, d, iref=0))  # Nonzero numerator / zero denominator.
    huge = np.array([[1e200]])
    assert np.isnan(compute(huge, data([[1e200]], [[1e300]]), iref=0))


def test_labels_do_not_change_factor_and_empty_beam_source_sums():
    c = np.array([[1., 2.]])
    d = data([[2, 8]])
    assert compute(c, d, iref=-1) == compute(c, replace(d, beam_indices_hk=((99, -2), (99, -2))), iref=-1)
    empty = data(np.empty((1, 0)))
    assert np.isnan(compute(np.empty((1, 0)), empty, iref=0))


@pytest.mark.parametrize("bad", [True, "1", np.nan, np.inf, -np.inf])
def test_invalid_iref(bad):
    with pytest.raises((ValueError, TypeError), match="iref"):
        compute(np.ones((1, 1)), data([[1]]), iref=bad)


@pytest.mark.parametrize("bad", [np.ones((1, 1), dtype=np.float32), np.ones(1),
    np.ones((2, 1)), np.array([[np.nan]]), np.array([[np.inf]]), [[1.]]])
def test_invalid_calculated(bad):
    with pytest.raises((ValueError, TypeError)):
        compute(bad, data([[1]]), iref=1)


def test_apply_ieee_outputs_and_invalid_types():
    result = apply(np.array([[0., 1., -1.]]), np.inf)
    assert np.isnan(result[0, 0]) and np.isposinf(result[0, 1]) and np.isneginf(result[0, 2])
    with pytest.raises(TypeError):
        apply(np.ones((1, 1)), True)
    with pytest.raises(TypeError):
        compute(np.ones((1, 1)), None, iref=0)


@pytest.mark.parametrize("filename", ["in.dat", "startingin.dat", "bestin.dat"])
def test_reference_file_regression(filename):
    beams = ((0, 0), (1, 1), (1, -1), (-1, -1), (-1, 1), (2, 0), (0, -2), (-2, 0), (0, 2))
    d = parse_experimental_intensities(Path(__file__).parents[1]/"reference_data"/filename,
                                       beam_indices_hk=beams)
    # Synthetic theory deliberately isolates scaling from unresolved GRADIENT.
    c = np.ones((1521, 9), dtype=np.float64)
    factor = compute(c, d, iref=-1)
    # Decimal arithmetic on source tokens: weighted abs sum=97953357.58305,
    # weight sum=7549. Tolerance permits literal-loop float64 accumulation.
    assert factor == pytest.approx(12975.673278983971, abs=1e-9, rel=0)
    assert factor == oracle(c, d, -1)
    assert compute(2*c, d, iref=-1) == factor/2
    assert np.isnan(compute(c, d, iref=0))  # Actual zero variances remain unguarded.
    assert compute(c, d, iref=1.25) == 1.25
