"""Patch A1 contracts only: no calibration, extraction or geometry algorithms."""
from dataclasses import FrozenInstanceError, replace

import numpy as np
import pytest

from qcbed.detector import DetectorResponseModel, PSFCalibrationMetrics
from qcbed.data_preparation import (BeamObservation, DiskGeometry, NoiseModel,
    PreparedExperimentalDataset, ReciprocalFrame, ReflectionROI, ScanGeometry, SignalMode)


def frame():
    return ReciprocalFrame((10., 20.), (11., 23.), (3., 0.), (0., -4.))


def scan():
    return ScanGeometry((-.508156, .407641), (.407641, .508156),
                        (-.407641, -.508156), 2, 3)


def dataset(n=2, **overrides):
    values = dict(intensities=np.arange(6*n, dtype=np.float64).reshape(6, n)-3,
        variances=np.full((6, n), 4., dtype=np.float64),
        weights=np.full((6, n), .7, dtype=np.float64),
        beams=tuple(BeamObservation(f'b{i}', (i+.25, -i), (1, 1), enabled=i != 0)
                    for i in range(n)), rois=(), reciprocal_frame=frame(),
        scan_geometry=scan(), mode=SignalMode.CONVENTIONAL, intensity_units='counts')
    values.update(overrides)
    return PreparedExperimentalDataset(**values)


@pytest.mark.parametrize('indices', [None, (1, -1), (1, -1, 2)])
def test_beam_id_assignment_and_enabled(indices):
    first = BeamObservation('first', [-.25, 12.3], indices, False, 'direct')
    second = replace(first, id='second', enabled=True)
    assert first.indices == second.indices
    assert first.center_px == (-.25, 12.3)
    assert first.id != second.id and not first.enabled and second.enabled
    with pytest.raises(FrozenInstanceError):
        first.id = 'changed'


@pytest.mark.parametrize('updates', [dict(id=''), dict(center_px=(1,)),
    dict(center_px=(1, 2, 3)), dict(center_px=(True, 2)), dict(center_px=(np.inf, 1)),
    dict(center_px=('1', 2)), dict(indices=(1.5, 2)), dict(indices=(True, 2)),
    dict(indices=(1,)), dict(enabled=1), dict(label=3)])
def test_invalid_beam(updates):
    with pytest.raises((ValueError, TypeError)):
        replace(BeamObservation('id', (1, 2)), **updates)


def test_independent_origins_and_owned_coordinates():
    center = [1., 2.]
    value = ReciprocalFrame(center, (3, 4), (1, 0), (0, 1))
    center[0] = 100
    assert value.central_disk_center_px == (1, 2)
    changed = replace(value, central_disk_center_px=(20, 30))
    assert changed.zone_axis_origin_px == value.zone_axis_origin_px == (3, 4)
    assert changed.basis_vector_1_px == (1, 0)


@pytest.mark.parametrize('field', ['central_disk_center_px', 'zone_axis_origin_px',
                                  'basis_vector_1_px', 'basis_vector_2_px'])
@pytest.mark.parametrize('value', [(1,), (1, 2, 3), (np.nan, 1), (False, 1)])
def test_frame_invalid_pairs(field, value):
    with pytest.raises((ValueError, TypeError)):
        replace(frame(), **{field: value})


@pytest.mark.parametrize('field', ['basis_vector_1_px', 'basis_vector_2_px'])
def test_frame_zero_basis(field):
    with pytest.raises(ValueError):
        replace(frame(), **{field: (0, 0)})


@pytest.mark.parametrize('value', [0, -1, np.inf, np.nan, True])
@pytest.mark.parametrize('field', ['radius_x_px', 'radius_y_px'])
def test_ellipse_radii(field, value):
    with pytest.raises(ValueError):
        replace(DiskGeometry(2, 3), **{field: value})


def test_roi_relationship_and_separate_center():
    roi = ReflectionROI('r', 'b0', (9, 7), DiskGeometry(2, 4), False)
    data = dataset(rois=[roi, replace(roi, id='another')])
    assert data.rois[0].center_px != data.beams[0].center_px
    assert data.rois[0].geometry.radius_x_px == 2
    with pytest.raises(ValueError, match='unknown beam'):
        dataset(rois=[replace(roi, beam_id='missing')])
    with pytest.raises(ValueError, match='unique stable'):
        dataset(rois=[roi, roi])
    with pytest.raises(ValueError, match='unique stable'):
        dataset(beams=[data.beams[0], data.beams[0]])


@pytest.mark.parametrize('updates', [dict(corner_a=(1,)), dict(corner_b=(1, 2, 3)),
    dict(corner_c=(np.inf, 0)), dict(scan_rows=0), dict(scan_columns=-1),
    dict(scan_rows=1.5), dict(scan_columns=True)])
def test_scan_validation(updates):
    with pytest.raises((ValueError, TypeError)):
        replace(scan(), **updates)


def test_scan_golden_values_are_stored_not_computed():
    assert scan().corner_a == (-.508156, .407641)
    assert scan().corner_b == (.407641, .508156)
    assert scan().corner_c == (-.407641, -.508156)
    assert replace(scan(), corner_a=(30, -20)).corner_a == (30, -20)


@pytest.mark.parametrize('count', [1, 2, 9, 23])
def test_dynamic_order_and_independent_weights(count):
    data = dataset(count)
    assert data.intensities.shape == (6, count)
    assert [b.id for b in data.beams] == [f'b{i}' for i in range(count)]
    assert len({b.indices for b in data.beams}) == 1
    assert not np.array_equal(data.weights, 1/data.variances)
    assert data.intensities[0, 0] == -3
    assert data.detector_response is None
    assert data.mode is SignalMode.CONVENTIONAL


@pytest.mark.parametrize('field', ['intensities', 'variances', 'weights'])
@pytest.mark.parametrize('bad', [np.zeros((2, 6)), np.zeros((6, 2, 1)),
    np.zeros((6, 2), dtype=np.float32), np.full((6, 2), np.nan)])
def test_dataset_arrays_reject_inconsistent_values(field, bad):
    with pytest.raises((ValueError, TypeError)):
        dataset(**{field: bad})


def test_arrays_are_owned_and_strongly_readonly():
    source = np.arange(12., dtype=np.float64).reshape(2, 6).T
    data = dataset(intensities=source, provenance=['source: test'])
    expected = source.copy()
    source[:] = 200
    np.testing.assert_array_equal(data.intensities, expected)
    for array in (data.intensities, data.variances, data.weights):
        with pytest.raises(ValueError):
            array.setflags(write=True)
        with pytest.raises(ValueError):
            array[0, 0] = 5
    assert data.provenance == ('source: test',)


def test_optional_metadata_and_zero_variance_weight():
    response = DetectorResponseModel('unknown method', artifact_reference='relative/calibration.npz')
    noise = NoiseModel(.2, .3, .4, 'counts')
    data = dataset(detector_response=response, noise_model=noise, mode=SignalMode.GRADIENT,
        variances=np.zeros((6, 2)), weights=np.zeros((6, 2)))
    assert data.detector_response is response
    assert response.artifact_reference == 'relative/calibration.npz'
    assert not hasattr(data, 'psf_corrected')
    assert data.noise_model == noise
    with pytest.raises(TypeError):
        dataset(mode='gradient')
    for name in ('variances', 'weights'):
        with pytest.raises(ValueError):
            dataset(**{name: -np.ones((6, 2))})


@pytest.mark.parametrize('values', [(np.nan, 1, 2), (1, np.inf, 2), (True, 1, 2)])
def test_noise_finite_coefficients(values):
    with pytest.raises(ValueError):
        NoiseModel(*values, intensity_units='counts')


def test_detector_optional_metadata_and_shapes():
    metrics = PSFCalibrationMetrics(2, 3, 4)
    psf = np.ones((3, 4)) * 2
    model = DetectorResponseModel('aperture', psf=psf,
        otf=np.full((3, 4), 1+2j, dtype=np.complex128), mtf=np.ones((3, 4)),
        metrics=metrics, provenance=['raw aperture'])
    psf[:] = 0
    assert np.all(model.psf == 2)  # no normalization or aliasing
    assert model.metrics == metrics
    assert model.provenance == ('raw aperture',)
    with pytest.raises(ValueError):
        model.otf.setflags(write=True)
    assert DetectorResponseModel('unspecified').psf is None
    with pytest.raises(ValueError, match='grid shape'):
        DetectorResponseModel('aperture', psf=np.ones((3, 4)), mtf=np.ones((4, 3)))


@pytest.mark.parametrize('updates', [dict(psf=np.ones(3)), dict(otf=np.ones((3, 4))),
    dict(mtf=np.full((3, 4), np.inf)), dict(psf=np.zeros((0, 2))),
    dict(metrics={}), dict(provenance='text'), dict(artifact_reference='')])
def test_detector_invalid(updates):
    with pytest.raises((ValueError, TypeError)):
        DetectorResponseModel('aperture', **updates)


@pytest.mark.parametrize('updates', [dict(goodness_of_fit=-1), dict(poisson_mismatch=np.nan),
    dict(accepted_measurements=-1), dict(accepted_measurements=True), dict(accepted_measurements=1.5)])
def test_metrics_invalid(updates):
    with pytest.raises((ValueError, TypeError)):
        PSFCalibrationMetrics(**updates)
