"""Synthetic whole-frame Fourier correction checks; not DM golden fixtures."""
import numpy as np
import pytest

from qcbed.detector import DetectorResponseModel, correct_cbed_with_otf


def response(otf):
    return DetectorResponseModel('synthetic', otf=np.asarray(otf, dtype=np.complex128))


@pytest.mark.parametrize('shape', [(8, 8), (7, 11), (8, 13)])
def test_identity_scale_sign_and_immutability(shape):
    raw = np.random.default_rng(3).normal(10, 20, shape)
    before = raw.copy()
    model = response(np.ones(shape))
    original_otf = model.otf.copy()
    corrected = correct_cbed_with_otf(raw, model)
    np.testing.assert_allclose(corrected, raw, atol=2e-14)
    np.testing.assert_array_equal(raw, before)
    np.testing.assert_array_equal(model.otf, original_otf)
    assert corrected.dtype == np.float64 and np.isfinite(corrected).all()
    assert not np.shares_memory(corrected, raw)
    np.testing.assert_allclose(correct_cbed_with_otf(raw, response(2 * model.otf)), raw / 2, atol=2e-14)


@pytest.mark.parametrize('shape', [(8, 10), (9, 13)])
def test_known_asymmetric_blur_recovery_and_complex_phase(shape):
    ideal = np.random.default_rng(5).uniform(0, 100, shape)
    psf = np.zeros(shape)
    y, x = np.array(shape) // 2
    psf[y, x] = .8
    psf[y, x + 1] = .2
    otf = np.fft.fft2(np.fft.ifftshift(psf))
    raw = np.fft.ifft2(np.fft.fft2(ideal) * otf).real
    recovered = correct_cbed_with_otf(raw, response(otf))
    np.testing.assert_allclose(recovered, ideal, atol=1e-12, rtol=1e-13)
    magnitude_only = correct_cbed_with_otf(raw, response(np.abs(otf)))
    assert np.max(np.abs(magnitude_only - ideal)) > 1


def test_shifted_delta_centering():
    raw = np.arange(77.).reshape(7, 11)
    psf = np.zeros_like(raw)
    psf[3, 6] = 1
    otf = np.fft.fft2(np.fft.ifftshift(psf))
    np.testing.assert_allclose(correct_cbed_with_otf(np.roll(raw, 1, axis=1), response(otf)), raw, atol=5e-14)


def test_small_and_zero_transfer_bins():
    raw = np.tile([1., -1., 1., -1.], (4, 1))
    otf = np.ones((4, 4), complex)
    otf[0, 2] = 1e-15
    np.testing.assert_allclose(correct_cbed_with_otf(raw, response(otf), fourier_relative_floor=1e-3), raw * 1e3)
    otf[0, 2] = 0
    np.testing.assert_array_equal(correct_cbed_with_otf(raw, response(otf)), np.zeros_like(raw))
    with pytest.raises(ValueError, match='all components are zero'):
        correct_cbed_with_otf(raw, response(np.zeros_like(otf)))


def test_imaginary_residual_policy():
    raw = np.ones((4, 6))
    with pytest.raises(ValueError, match='significant imaginary residual'):
        correct_cbed_with_otf(raw, response(np.full(raw.shape, 1j)))
    np.testing.assert_allclose(correct_cbed_with_otf(raw, response(np.full(raw.shape, 1 + 1e-14j))), raw)
    with pytest.raises(ValueError, match='imaginary residual'):
        correct_cbed_with_otf(raw, response(np.full(raw.shape, 1 + 1e-14j)),
                              imaginary_absolute_tolerance=0, imaginary_relative_tolerance=0)


@pytest.mark.parametrize('raw', [np.ones(3), np.ones((0, 3)), np.full((3, 3), np.nan), np.full((3, 3), np.inf)])
def test_invalid_arrays(raw):
    with pytest.raises(ValueError):
        correct_cbed_with_otf(raw, response(np.ones((3, 3))))


@pytest.mark.parametrize('raw', [[[1]], np.ones((3, 3), complex), np.ones((3, 3), bool), np.full((3, 3), '1')])
def test_invalid_types(raw):
    with pytest.raises(TypeError):
        correct_cbed_with_otf(raw, response(np.ones((3, 3))))


def test_missing_response_and_shape():
    raw = np.ones((3, 3))
    with pytest.raises(TypeError):
        correct_cbed_with_otf(raw, None)
    with pytest.raises(ValueError, match='provide an OTF'):
        correct_cbed_with_otf(raw, DetectorResponseModel('metadata'))
    with pytest.raises(ValueError, match='identical'):
        correct_cbed_with_otf(raw, response(np.ones((3, 4))))


@pytest.mark.parametrize('kwargs', [dict(fourier_relative_floor=0), dict(fourier_relative_floor=1),
    dict(fourier_relative_floor=np.nan), dict(imaginary_relative_tolerance=-1),
    dict(imaginary_absolute_tolerance=np.inf), dict(fourier_relative_floor=True)])
def test_invalid_settings(kwargs):
    with pytest.raises(ValueError):
        correct_cbed_with_otf(np.ones((3, 3)), response(np.ones((3, 3))), **kwargs)


def test_overflow_reported():
    with pytest.raises(ValueError, match='Unusable OTF correction'):
        correct_cbed_with_otf(np.full((4, 4), 1e308), response(np.ones((4, 4))))
