"""Profile geometry and compatibility-only legacy intensity localization."""
import numpy as np
import pytest
from qcbed.data_preparation import estimate_disk_geometry
from qcbed.data_preparation.legacy_disk_localization import (DiskGeometry,
    sum_bin_image, ellipse_integrated_intensity, locate_central_disk_intensity as locate_central_disk)


def disk(shape=(128, 192), center=(135, 43), radii=(14, 9), noise=False):
    y, x = np.indices(shape)
    r = ((x-center[0])/radii[0])**2 + ((y-center[1])/radii[1])**2
    image = 100.*(r < 1) + 10*np.exp(-r)
    if noise:
        image += np.random.default_rng(23).normal(0, 1, shape) + .001*x
    return image


def test_geometry_separate_reference_threshold_and_axes():
    image = disk(); before = image.copy()
    result = estimate_disk_geometry(image)
    first = image > image.mean()
    threshold = image[first].mean()/2
    selected = image > threshold
    assert result.initial_threshold == image.mean()
    assert result.final_threshold == threshold
    assert result.selected_pixels == selected.sum()
    assert result.selected_mean_intensity == image[selected].mean()
    assert result.geometry == DiskGeometry(13.5, 8.5)
    np.testing.assert_array_equal(image, before)


@pytest.mark.parametrize('noise', [False, True])
def test_displaced_center_multiscale_rectangular(noise):
    image = disk(noise=noise); before = image.copy()
    geometry = DiskGeometry(14, 9)
    result = locate_central_disk(image, geometry, binning_steps=3)
    np.testing.assert_allclose(result.central_disk_center_px, (135, 43), atol=1)
    assert result.mode == 'automatic' and result.geometry == geometry
    assert [level.bin_factor for level in result.levels] == [8, 4, 2, 1]
    assert result.levels[-1].radius_scale == .9
    assert result.integrated_intensity == result.levels[-1].integrated_intensity
    for previous, current in zip(result.levels, result.levels[1:]):
        assert np.max(abs(np.array(current.center_px)-2*np.array(previous.center_px))) <= 1
    np.testing.assert_array_equal(image, before)
    assert not hasattr(result, 'zone_axis_origin_px')


def test_sum_binning_conservation_and_blocks():
    image = np.arange(96.).reshape(8, 12)
    before = image.copy(); result = sum_bin_image(image, 4)
    expected = np.array([[image[y:y+4, x:x+4].sum() for x in (0, 4, 8)] for y in (0, 4)])
    np.testing.assert_array_equal(result, expected)
    assert result.sum() == image.sum()
    np.testing.assert_array_equal(image, before)
    with pytest.raises(ValueError, match='divide'):
        sum_bin_image(image, 3)


def test_ellipse_objective_strict_boundary_and_not_brightest_pixel():
    image = disk(); image[10, 10] = 500
    geometry = DiskGeometry(14, 9)
    y, x = np.indices(image.shape)
    expected = image[((x-135)/14)**2+((y-43)/9)**2 < 1].sum()
    assert ellipse_integrated_intensity(image, geometry, (135, 43)) == expected
    assert expected > ellipse_integrated_intensity(image, geometry, (10, 10))
    result = locate_central_disk(image, geometry, binning_steps=3)
    np.testing.assert_allclose(result.central_disk_center_px, (135, 43), atol=1)
    assert ellipse_integrated_intensity(np.ones((9, 11)), DiskGeometry(1, 1), (4, 3)) == 1


def test_preset_bypasses_search_and_divisibility(monkeypatch):
    import qcbed.data_preparation.legacy_disk_localization as module
    def fail(*args, **kwargs):
        raise AssertionError('automatic binning invoked')
    monkeypatch.setattr(module, 'sum_bin_image', fail)
    image = np.ones((13, 21))
    result = locate_central_disk(image, DiskGeometry(3, 2), preset_center_px=(17.25, 3.5))
    assert result.central_disk_center_px == (17.25, 3.5)
    assert result.mode == 'preset' and result.levels == ()
    assert result.integrated_intensity > 0


@pytest.mark.parametrize('center', [(-1, 3), (3, -1), (21, 3), (3, 13), (np.nan, 2), (1,), (True, 2)])
def test_invalid_preset(center):
    with pytest.raises((ValueError, TypeError)):
        locate_central_disk(np.ones((13, 21)), DiskGeometry(3, 2), preset_center_px=center)


@pytest.mark.parametrize('image', [np.ones((8, 12)), np.zeros((8, 12))])
def test_degenerate_reference(image):
    with pytest.raises(ValueError, match='threshold'):
        estimate_disk_geometry(image)


@pytest.mark.parametrize('image', [np.ones(5), np.ones((0, 5)), np.full((8, 12), np.nan), np.full((8, 12), np.inf)])
def test_invalid_images(image):
    with pytest.raises(ValueError):
        estimate_disk_geometry(image)
    with pytest.raises(ValueError):
        locate_central_disk(image, DiskGeometry(3, 2))


@pytest.mark.parametrize('steps', [0, -1, 50, True, 1.5])
def test_invalid_binning(steps):
    with pytest.raises((TypeError, ValueError)):
        locate_central_disk(disk(), DiskGeometry(14, 9), binning_steps=steps)


def test_indivisible_and_empty_search():
    with pytest.raises(ValueError, match='divide'):
        locate_central_disk(np.ones((33, 40)), DiskGeometry(3, 2), binning_steps=2)
    with pytest.raises(ValueError, match='No in-frame'):
        locate_central_disk(np.ones((32, 40)), DiskGeometry(100, 100), binning_steps=2)


def test_signed_objective_tie_determinism():
    image = -np.ones((16, 24))
    result = locate_central_disk(image, DiskGeometry(4, 4), binning_steps=2)
    assert result == locate_central_disk(image, DiskGeometry(4, 4), binning_steps=2)
    assert np.isfinite(result.integrated_intensity)


def test_excessive_binning_cannot_silently_select_empty_ellipse():
    with pytest.raises(ValueError, match='no sampled pixels'):
        locate_central_disk(disk(), DiskGeometry(14, 9), binning_steps=5)


def test_reported_level_scores_and_search_phase():
    image = disk(); geometry = DiskGeometry(14, 9)
    result = locate_central_disk(image, geometry, binning_steps=3)
    first = result.levels[0]
    assert (first.center_px[0]-14/8) % 1 == 0
    assert (first.center_px[1]-9/8) % 1 == 0
    for level in result.levels:
        scaled = DiskGeometry(14*level.radius_scale/level.bin_factor,
                              9*level.radius_scale/level.bin_factor)
        expected = ellipse_integrated_intensity(sum_bin_image(image, level.bin_factor),
                                                scaled, level.center_px)
        assert level.integrated_intensity == expected
