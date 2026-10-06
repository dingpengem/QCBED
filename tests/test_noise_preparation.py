"""Source-equation tests and synthetic fits, not claimed DM end-to-end equivalence."""
import numpy as np
import pytest
from qcbed.data_preparation import (NoiseModel, NoiseFitConfig, characterize_noise,
                                   fit_noise_model, variance_map_from_noise_model)
from qcbed.data_preparation.noise import (_smooth, _valid_mask, _reference_rms,
                                         _blemishes, _slice_points, _rebin, _fit_points)


def pattern():
    y, x = np.indices((128, 192))
    return 200 + 10*x + np.random.default_rng(3).normal(0, 5, x.shape)


def test_fourth_order_smoothing_against_independent_fft():
    raw = np.random.default_rng(8).normal(size=(24, 32))
    psf = np.zeros_like(raw)
    for dy in (-1, 0, 1):
        for dx in (-1, 0, 1):
            psf[dy % 24, dx % 32] = (1 if dy == 0 else .5)*(1 if dx == 0 else .5)/4
    assert psf.sum() == 1
    transfer = np.fft.fft2(psf)
    expected = np.fft.ifft2(np.fft.fft2(raw)*(4*transfer-6*transfer**2+4*transfer**3-transfer**4)).real
    np.testing.assert_allclose(_smooth(raw), expected, atol=2e-15)
    np.testing.assert_array_equal(_smooth(np.full((12, 20), 7.)), np.full((12, 20), 7.))


def test_masks_gain_shift_holes_and_supplied_override():
    raw = np.ones((24, 32))
    raw[8, 8] = 0; raw[9, 9] = -1; raw[10, 10] = 40001; raw[11, 11] = 40000
    auto = _valid_mask(raw, None, 40000)
    assert not auto[:4].any() and not auto[-4:].any()
    assert not auto[:, :4].any() and not auto[:, -4:].any()
    assert not auto[8, 8] and not auto[9, 9] and not auto[10, 10] and auto[11, 11]
    supplied = np.ones(raw.shape, bool); supplied[7, 7] = False
    before = supplied.copy()
    mask = _valid_mask(raw, supplied, 40000)
    assert mask[10, 10]  # supplied branch does not apply gain threshold
    assert not mask[7, 7] and not mask[11, 11]
    np.testing.assert_array_equal(supplied, before)


def test_literal_central_fourier_crop_and_scale():
    raw = np.full((24, 32), 7.)
    mask = _valid_mask(raw, None, 40000)
    # Reduced inverse FFT quadruples the DC amplitude; source then divides by 2.
    assert _reference_rms(raw, mask) == pytest.approx(14.)
    raw += np.tile([1., -1.], (24, 16))
    assert _reference_rms(raw, mask) == pytest.approx(14.)


def test_blemish_threshold_and_impulse():
    raw = np.zeros((64, 80)); raw[30, 30] = 1000
    residual = raw - _smooth(raw)
    assert residual[30, 30] > 0
    mask = np.ones(raw.shape, bool)
    blemish, rms = _blemishes(residual, mask, 10)
    assert rms > 0 and not blemish[30, 30] and blemish[10, 10]
    permissive, _ = _blemishes(residual, mask, 100)
    assert permissive[30, 30]


def test_slicing_strict_lower_inclusive_upper_and_rms():
    smooth = np.array([[0., 1., 2., 3., 4.]])
    high = np.array([[9., 2., -4., 6., -8.]])
    mask = np.array([[False, True, True, True, True]])
    points = _slice_points(smooth, high, mask, 2, 1)
    np.testing.assert_allclose(points, [[3.5, np.sqrt(50), 2], [1.5, np.sqrt(10), 2]])


def test_rebin_population_weighting_and_omitted_remainder():
    points = np.array([[9., 3., 2.], [6., 2., 1.], [3., 1., 3.], [1., .5, 1.]])
    fit, remainder = _rebin(points)
    np.testing.assert_allclose(fit, [[8., np.sqrt(22/3), 3.], [3., 1., 3.]])
    np.testing.assert_array_equal(remainder, [1., .5, 1.])
    assert fit[:, 2].sum() + remainder[2] == points[:, 2].sum()


@pytest.mark.parametrize('dark', [None, 0., 3.])
def test_known_coefficients_and_dark_constant(dark):
    intensity = np.linspace(1, 10000, 80)
    c = 3 if dark is None else dark
    sigma = .02*intensity + 1.4*np.sqrt(intensity) + c
    points = np.column_stack((intensity, sigma, np.arange(1, 81)))
    model, gof = _fit_points(points, dark, 'counts')
    np.testing.assert_allclose([model.a, model.b, model.c], [.02, 1.4, c], atol=1e-10)
    assert gof < 1e-20
    if dark is not None:
        assert model.c == dark


def test_weighted_objective_and_nonnegative_fit():
    intensity = np.arange(1., 31.)**2
    sigma = 1 + .2*np.sqrt(intensity)
    sigma[0] += 4
    weights = np.arange(1., 31.)
    points = np.column_stack((intensity, sigma, weights))
    model, gof = _fit_points(points, None, 'counts')
    assert min(model.a, model.b, model.c) >= 0
    prediction = model.a*intensity + model.b*np.sqrt(intensity) + model.c
    assert gof == pytest.approx(np.sum(weights*(sigma-prediction)**2)/(weights.sum()*(1-3/30)))
    _, unweighted = _fit_points(np.column_stack((intensity, sigma, np.ones(30))), None, 'counts')
    assert gof != pytest.approx(unweighted)


def test_variance_abs_intensity_and_no_normalization():
    raw = np.array([[-9., 0., 9.], [16., -16., 1.]])
    before = raw.copy()
    model = NoiseModel(.2, 2., 3., 'counts')
    expected = (.2*abs(raw) + 2*np.sqrt(abs(raw)) + 3)**2
    np.testing.assert_array_equal(variance_map_from_noise_model(raw, model), expected)
    np.testing.assert_array_equal(raw, before)


def test_end_to_end_rectangular_diagnostics_and_determinism():
    raw = pattern(); before = raw.copy()
    config = NoiseFitConfig(sampling_stride=1)
    state = np.random.get_state()
    result = characterize_noise(raw, config=config)
    repeated = characterize_noise(raw, config=config)
    assert result.model == repeated.model == fit_noise_model(raw, config=config)
    np.testing.assert_array_equal(result.high_frequency, repeated.high_frequency)
    np.testing.assert_array_equal(np.random.get_state()[1], state[1])
    np.testing.assert_array_equal(raw, before)
    np.testing.assert_array_equal(result.noise_mask, result.valid_mask & result.blemish_mask & (result.smooth_estimate > 0))
    residual = raw-result.smooth_estimate
    np.testing.assert_allclose(result.high_frequency, residual*result.high_frequency_normalization)
    assert result.high_frequency_normalization > 0 and len(result.fit_points) > 3
    assert np.isfinite(variance_map_from_noise_model(raw, result.model)).all()
    for array in (result.high_frequency, result.fit_points, result.noise_mask):
        with pytest.raises(ValueError):
            array.setflags(write=True)
    fixed = characterize_noise(raw, config=config, dark_noise=7)
    assert fixed.model.c == 7


@pytest.mark.parametrize('raw', [np.ones(5), np.ones((0, 4)), np.full((16, 20), np.nan),
                                np.full((16, 20), np.inf), np.ones((15, 20)), np.ones((8, 20))])
def test_invalid_patterns(raw):
    with pytest.raises(ValueError):
        fit_noise_model(raw)


@pytest.mark.parametrize('raw', [[[1]], np.ones((16, 20), complex), np.ones((16, 20), bool)])
def test_invalid_pattern_types(raw):
    with pytest.raises(TypeError):
        fit_noise_model(raw)


def test_invalid_masks_and_empty_or_degenerate_support():
    raw = pattern()
    with pytest.raises(ValueError, match='shape'):
        fit_noise_model(raw, valid_mask=np.ones((3, 4), bool))
    with pytest.raises(TypeError, match='Boolean'):
        fit_noise_model(raw, valid_mask=np.ones(raw.shape))
    with pytest.raises(ValueError, match='No valid'):
        fit_noise_model(raw, valid_mask=np.zeros(raw.shape, bool))
    with pytest.raises(ValueError, match='Degenerate'):
        fit_noise_model(np.ones((16, 20)))
    with pytest.raises(ValueError, match='more than three'):
        _fit_points(np.ones((3, 3)), None, 'counts')
    with pytest.raises(ValueError, match='underdetermined'):
        _fit_points(np.ones((5, 3)), None, 'counts')


@pytest.mark.parametrize('kwargs', [dict(gain_normalization_limit=0), dict(blemish_sigma=-1),
    dict(intensity_slices=0), dict(sampling_stride=1.5), dict(intensity_slices=True),
    dict(blemish_sigma=np.nan)])
def test_invalid_config(kwargs):
    with pytest.raises((TypeError, ValueError)):
        NoiseFitConfig(**kwargs)


@pytest.mark.parametrize('dark', [-1, np.nan, True])
def test_invalid_dark(dark):
    with pytest.raises(ValueError):
        fit_noise_model(pattern(), dark_noise=dark)


def test_variance_overflow_and_negative_coefficients():
    with pytest.raises(ValueError, match='Non-finite'):
        variance_map_from_noise_model(np.full((3, 4), 1e300), NoiseModel(1, 0, 0, 'counts'))
    with pytest.raises(ValueError, match='nonnegative'):
        variance_map_from_noise_model(np.ones((3, 4)), NoiseModel(-1, 0, 0, 'counts'))


def test_supplied_pattern_sigma_empirical_model_consistency():
    # These files lack correction/mask/fit settings: this checks their numerical
    # association only, not recovery of the historical characterization pipeline.
    from pathlib import Path
    root = Path(__file__).resolve().parents[1] / 'tests/fixtures/legacy_reference'
    intensity = np.abs(np.loadtxt(Path(__file__).resolve().parents[1] / 'reference_data/pattern.dat')).ravel()
    sigma = np.loadtxt(root / 'sigma.dat').ravel()
    assert intensity.shape == sigma.shape
    training = np.arange(intensity.size) % 2 == 0
    points = np.column_stack((intensity[training], sigma[training], np.ones(training.sum())))
    model, _ = _fit_points(points, None, 'stored intensity units')
    held_out = intensity[~training].reshape(1, -1)
    calculated_sigma = np.sqrt(variance_map_from_noise_model(held_out, model)).ravel()
    # Stored sigma is rounded (up to 6 displayed significant digits).
    np.testing.assert_allclose(calculated_sigma, sigma[~training], rtol=0, atol=1.)


def test_source_defaults():
    config = NoiseFitConfig()
    assert (config.gain_normalization_limit, config.blemish_sigma,
            config.intensity_slices, config.sampling_stride) == (40000, 10, 100, 16)
