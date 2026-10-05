"""AMMPSF source-equation and synthetic regressions, not executed DM equivalence."""
import json
from dataclasses import replace

import numpy as np
import pytest

from qcbed.detector import (PSFCalibrationConfig, PSFCalibrationError,
                           calibrate_psf_from_aperture)
from qcbed.detector import psf as m


def synthetic(shape=(64, 80)):
    y, x = np.indices(shape)
    cy, cx = np.array(shape)//2
    # Union of offset rectangles: no assumption of circular aperture geometry.
    ideal = 100.*(((x > cx-18) & (x < cx+10) & (y > cy-16) & (y < cy+12)) |
                  ((x > cx-6) & (x < cx+20) & (y > cy-2) & (y < cy+18)))
    kernel = np.exp(-((x-cx)**2+(y-cy)**2)/(2*.7**2))
    kernel /= kernel.sum()
    measured = np.fft.ifft2(np.fft.fft2(ideal)*np.fft.fft2(np.fft.ifftshift(kernel))).real
    return ideal, kernel, measured


QUICK = PSFCalibrationConfig(accepted_measurements=2, max_attempts=5, max_iterations=30)


@pytest.mark.parametrize('shape', [(64, 64), (64, 80), (53, 61)])
def test_synthetic_blur_and_normalization(shape):
    ideal, known, aperture = synthetic(shape)
    original = aperture.copy()
    result = calibrate_psf_from_aperture(aperture, config=QUICK, seed=3)
    np.testing.assert_array_equal(aperture, original)
    for array in (result.psf, result.otf, result.mtf):
        assert array.shape == shape and np.isfinite(array).all()
        assert not array.flags.writeable
    assert result.psf.dtype == np.float64
    assert result.otf.dtype == np.complex128
    assert np.unravel_index(result.psf.argmax(), shape) == tuple(np.array(shape)//2)
    assert result.psf.sum() == pytest.approx(1, abs=2e-15)
    assert result.otf[0, 0] == pytest.approx(1, abs=2e-15)
    np.testing.assert_allclose(result.otf, np.fft.fft2(np.fft.ifftshift(result.psf)), rtol=0, atol=0)
    np.testing.assert_array_equal(result.mtf, abs(result.otf))
    assert result.mtf[0, shape[1]//2] < .5*result.mtf[0, 0]
    predicted = np.fft.ifft2(np.fft.fft2(ideal)*result.otf).real
    # Substantial improvement over treating the measured blur as an ideal image.
    assert np.mean((predicted-aperture)**2) < .1*np.mean((ideal-aperture)**2)
    assert np.sqrt(np.mean((predicted-aperture)**2)) < 1.  # intensity units, peak=100
    assert result.metrics.poisson_mismatch < 1
    assert result.metrics.accepted_measurements == 2
    record = json.loads(result.provenance[0])
    assert len([r for r in record['candidates'] if r['accepted']]) == 2
    assert record['metrics_scope'] == 'best accepted candidate, not ensemble fit'


@pytest.mark.parametrize('shape', [(7, 10), (8, 9)])
def test_fft_centering_and_shift_direction(shape):
    psf = np.zeros(shape)
    center = tuple(np.array(shape)//2)
    psf[center] = 1
    np.testing.assert_array_equal(m._otf_from_psf(psf), np.ones(shape))
    psf = np.roll(psf, 1, axis=1)
    expected = np.broadcast_to(np.exp(-2j*np.pi*np.arange(shape[1])/shape[1]), shape)
    np.testing.assert_allclose(m._otf_from_psf(psf), expected, atol=1e-15)
    np.testing.assert_allclose(np.fft.fftshift(m._real_inverse(expected)), psf, atol=2e-16)


def test_bilinear_interpolation_and_zero_extension():
    image = np.arange(12.).reshape(3, 4)
    out = m._sample(image, np.array([[.5, -1., 4.]]), np.array([[.5, 0., 0.]]))
    np.testing.assert_array_equal(out, [[2.5, 0., 0.]])


def test_shifted_threshold_aperture_source_formula():
    current = np.array([[0., 0., 0.], [0., 9., 0.], [0., 0., 0.]])
    sharp = m._sharp_aperture(current, 3., 1, 100.)
    assert sharp.sum() == 100 and sharp[1, 1] == 100
    # Independently enumerated bilinear samples for dim=3 at threshold 3:
    # center survives 5 shifts, each axial neighbor survives 1, corners none.
    expected = np.array([[0., 1., 0.], [1., 5., 1.], [0., 1., 0.]])*100/9
    np.testing.assert_allclose(m._sharp_aperture(current, 3., 3, 100.), expected)
    with pytest.raises(PSFCalibrationError):
        m._sharp_aperture(current, 10., 1, 100.)


def test_rotational_coverage_and_distinct_core_stages(monkeypatch):
    # A constant field stays constant after the legacy coverage correction.
    coverage = m._peak_normalize(m._rotational_sum(np.ones((15, 19))))
    constant = m._rotational_sum(np.full((15, 19), 2.))/coverage
    np.testing.assert_allclose(constant, constant[7, 9], rtol=1e-14)
    seen = []
    original = m._rotational_sum
    def record(image, start=1, stop=16):
        seen.append((start, stop))
        return original(image, start, stop)
    monkeypatch.setattr(m, '_rotational_sum', record)
    rough = np.zeros((15, 19)); rough[7, 9] = 1
    y, x = np.indices(rough.shape); radius = np.hypot(x-9, y-7)
    m._regularize_psf(rough, coverage, radius <= 2, radius <= 5, 1)
    assert seen == [(1, 2), (3, 16), (1, 16), (1, 16)]


def test_otf_positive_negative_and_path_control():
    positive = np.array([[1., .2], [.3, .4]], dtype=complex)
    unchanged, dim, flag = m._positive_otf(positive, 3)
    np.testing.assert_array_equal(unchanged, positive)
    assert dim == 5 and not flag
    negative = positive.copy(); negative[0, 1] = -1e-10 + 2e-11j
    original = negative.copy()
    corrected, dim, flag = m._positive_otf(negative, 3)
    np.testing.assert_array_equal(negative, original)
    np.testing.assert_allclose(corrected, (negative+2e-10)/(1+2e-10), rtol=0, atol=0)
    assert dim == 3 and flag and corrected.real.min() > 0


@pytest.mark.parametrize('value', [0, -1e31, np.nan, np.inf])
def test_unusable_otf(value):
    with pytest.raises(PSFCalibrationError):
        m._positive_otf(np.array([[1, value]], dtype=complex), 1)


def test_fourier_safeguard_explicit_and_finite():
    numerator = np.array([[2., 3., 4.]], dtype=complex)
    denominator = np.array([[1., 1e-20j, 0.]], dtype=complex)
    result, guarded = m._fourier_divide(numerator, denominator, 1e-12)
    np.testing.assert_allclose(result, [[2, -3e12j, 0]])
    assert guarded == 2 and np.isfinite(result).all()
    with pytest.raises(PSFCalibrationError, match='all components'):
        m._fourier_divide(numerator, np.zeros_like(denominator), 1e-12)


def test_legacy_quality_equations():
    measured = np.array([[1., 2.], [4., 6.]])
    predicted = np.array([[2., 3.], [5., 8.]])
    weights = np.array([[.5, .75], [.75, 0.]])
    gof = np.sum(weights*(predicted-measured)**2)/weights.sum()
    mismatch = np.sqrt(gof/(np.sum(weights*measured)/weights.sum()))
    assert m._quality(measured, predicted) == pytest.approx((gof, mismatch))


def test_seed_determinism_and_global_rng_isolation():
    aperture = synthetic()[2]
    global_before = np.random.get_state()
    a = calibrate_psf_from_aperture(aperture, config=QUICK, seed=7)
    b = calibrate_psf_from_aperture(aperture, config=QUICK, seed=7)
    np.testing.assert_array_equal(a.psf, b.psf)
    assert a.metrics == b.metrics and a.provenance == b.provenance
    global_after = np.random.get_state()
    for old, new in zip(global_before, global_after):
        np.testing.assert_array_equal(old, new)
    c = calibrate_psf_from_aperture(aperture, config=QUICK, rng=np.random.default_rng(7))
    np.testing.assert_array_equal(a.otf, c.otf)
    with pytest.raises(ValueError, match='seed or rng'):
        calibrate_psf_from_aperture(aperture, seed=1, rng=np.random.default_rng(1))


@pytest.mark.parametrize('bad', ['aperture.tif', np.ones(64), np.ones((64, 64, 1)),
                               np.ones((64, 64), dtype=complex), np.ones((64, 64), dtype=bool)])
def test_input_type(bad):
    with pytest.raises(TypeError):
        calibrate_psf_from_aperture(bad)


@pytest.mark.parametrize('bad', [np.full((64, 64), np.nan), np.full((64, 64), np.inf),
                               np.zeros((64, 64)), np.ones((64, 64))])
def test_input_degenerate(bad):
    with pytest.raises(PSFCalibrationError):
        calibrate_psf_from_aperture(bad)


def test_radius_derived_size_validation(monkeypatch):
    with pytest.raises(ValueError, match='core radii'):
        calibrate_psf_from_aperture(np.arange(50*80.).reshape(50, 80))
    # No hard-coded global minimum: custom radii allow a small rectangular frame.
    def candidate(measured, *args):
        p = np.zeros_like(measured); p[measured.shape[0]//2, measured.shape[1]//2] = 1
        return (p, 1., 1., 1), 0, 0
    monkeypatch.setattr(m, '_candidate', candidate)
    config = PSFCalibrationConfig(small_core_radius_px=1, large_core_radius_px=2,
                                  accepted_measurements=1)
    result = calibrate_psf_from_aperture(np.arange(35).reshape(5, 7), config=config)
    assert result.psf.shape == (5, 7)


def test_candidate_rejection_selection_and_peak_ensemble(monkeypatch):
    _, _, aperture = synthetic()
    a = np.zeros_like(aperture); a[32, 40] = .8; a[32, 41] = .2
    b = np.zeros_like(aperture); b[32, 40] = .6; b[32, 39] = .4
    candidates = iter([(a, 3., 26., 1), (a, 2., 25., 3), (b, 1., 2., 5)])
    monkeypatch.setattr(m, '_candidate', lambda *args: (next(candidates), 0, 0))
    result = calibrate_psf_from_aperture(aperture, config=QUICK, seed=1)
    expected = a/.8+b/.6
    np.testing.assert_allclose(result.psf, expected/expected.sum())
    assert result.metrics.goodness_of_fit == 1
    rows = json.loads(result.provenance[0])['candidates']
    assert [r['accepted'] for r in rows] == [False, True, True]
    assert rows[-1]['selected_iteration'] == 5


@pytest.mark.parametrize('result', [None, (np.ones((64, 80)), 1., 30., 1)])
def test_no_acceptable_candidate(monkeypatch, result):
    monkeypatch.setattr(m, '_candidate', lambda *args: (result, 0, 0))
    with pytest.raises(PSFCalibrationError, match='Accepted 0/2'):
        calibrate_psf_from_aperture(synthetic()[2], config=QUICK, seed=1)


def test_candidate_budget_failure_is_not_a_delta_fallback():
    config = replace(QUICK, max_iterations=1, max_attempts=2)
    with pytest.raises(PSFCalibrationError, match='iteration budget'):
        calibrate_psf_from_aperture(synthetic()[2], config=config, seed=3)


@pytest.mark.parametrize('updates', [dict(small_core_radius_px=0), dict(large_core_radius_px=np.inf),
    dict(small_core_radius_px=30), dict(accepted_measurements=True), dict(max_attempts=1),
    dict(max_iterations=0), dict(fourier_relative_floor=0), dict(fourier_relative_floor=1)])
def test_configuration_validation(updates):
    with pytest.raises((ValueError, TypeError)):
        PSFCalibrationConfig(**updates)


def test_legacy_default_ensemble():
    result = calibrate_psf_from_aperture(synthetic()[2], seed=3)
    assert result.metrics.accepted_measurements == 10
    record = json.loads(result.provenance[0])
    assert record['parameters']['small_core_radius_px'] == 4.99
    assert record['parameters']['large_core_radius_px'] == 24.99
    assert sum(r['accepted'] for r in record['candidates']) == 10
