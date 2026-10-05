"""Immutable preprocessing models, independent of YAML, Qt and FOR007 structures.

Geometry points are (x, y) in detector pixels, arrays are (row, column) or explicitly
(tilt, beam). Nothing is fitted, indexed, normalized, clipped, or inferred here.
"""
from dataclasses import dataclass
from enum import Enum
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from .extraction import ExtractionConfig, WeightMode

import numpy as np
from numpy.typing import NDArray

from ..detector.models import (DetectorResponseModel, _array, _integer,
                               _provenance, _real, _text, _tuple)


def _pair(value, name):
    value = _tuple(value, name)
    if len(value) != 2:
        raise ValueError(f'{name} must have two components (x,y or basis coefficients)')
    return tuple(_real(v, name) for v in value)


def _enabled(value):
    if type(value) is not bool:
        raise TypeError('enabled must be bool')


class SignalMode(str, Enum):
    CONVENTIONAL = 'conventional'
    GRADIENT = 'gradient'


@dataclass(frozen=True)
class NoiseModel:
    """Coefficients of sigma(I)=a*I+b*sqrt(I)+c, not variance.

    If intensity units are U, sigma is U, a is dimensionless, b is sqrt(U), c is U.
    Coefficients are finite but not constrained by an assumed fit policy. No
    evaluation, fitting, negative-intensity handling or variance construction occurs.
    """
    a: float
    b: float
    c: float
    intensity_units: str

    def __post_init__(self):
        for name in ('a', 'b', 'c'):
            object.__setattr__(self, name, _real(getattr(self, name), name))
        _text(self.intensity_units, 'intensity_units')


@dataclass(frozen=True)
class DiskGeometry:
    """Axis-aligned ellipse semi-axes (radii), in detector pixels."""
    radius_x_px: float
    radius_y_px: float

    def __post_init__(self):
        for name in ('radius_x_px', 'radius_y_px'):
            value = _real(getattr(self, name), name)
            if value <= 0:
                raise ValueError(f'{name} must be positive')
            object.__setattr__(self, name, value)


@dataclass(frozen=True)
class BeamObservation:
    """Stable ID is independent of optional 2-index HK or 3-index HKL assignment.

    Duplicate assignments are valid. center_px is (x=column, y=row), may be
    fractional/negative, and is never wrapped or clipped to a frame.
    """
    id: str
    center_px: tuple[float, float]
    indices: tuple[int, int] | tuple[int, int, int] | None = None
    enabled: bool = True
    label: str | None = None

    def __post_init__(self):
        _text(self.id, 'beam id')
        object.__setattr__(self, 'center_px', _pair(self.center_px, 'center_px'))
        if self.indices is not None:
            values = _tuple(self.indices, 'indices')
            if len(values) not in (2, 3):
                raise ValueError('indices must contain HK or HKL')
            object.__setattr__(self, 'indices', tuple(_integer(v, 'index') for v in values))
        _enabled(self.enabled)
        if self.label is not None:
            _text(self.label, 'label')


@dataclass(frozen=True)
class ReflectionROI:
    """An ellipse associated by beam ID, with an explicit independent pixel center.

    This initial ROI type has no mask-generation behavior. Future shapes may be
    separate typed ROI variants rather than overloading ellipse dimensions.
    """
    id: str
    beam_id: str
    center_px: tuple[float, float]
    geometry: DiskGeometry
    enabled: bool = True

    def __post_init__(self):
        _text(self.id, 'ROI id')
        _text(self.beam_id, 'beam_id')
        object.__setattr__(self, 'center_px', _pair(self.center_px, 'center_px'))
        if not isinstance(self.geometry, DiskGeometry):
            raise TypeError('geometry must be DiskGeometry')
        _enabled(self.enabled)


@dataclass(frozen=True)
class ReciprocalFrame:
    """Independent detector pixel points and displacement vectors, all (x,y).

    basis_vector_*_px are displacements from zone_axis_origin_px, not absolute
    endpoints and not canonical reciprocal HK vectors. No origin or basis is
    derived. Nonzero vectors are required; linear independence is a later
    conversion/readiness concern, not an assumed interpretation of measurements.
    """
    central_disk_center_px: tuple[float, float]
    zone_axis_origin_px: tuple[float, float]
    basis_vector_1_px: tuple[float, float]
    basis_vector_2_px: tuple[float, float]

    def __post_init__(self):
        for name in ('central_disk_center_px', 'zone_axis_origin_px',
                     'basis_vector_1_px', 'basis_vector_2_px'):
            value = _pair(getattr(self, name), name)
            if name.startswith('basis_vector') and value == (0., 0.):
                raise ValueError(f'{name} must be nonzero')
            object.__setattr__(self, name, value)


@dataclass(frozen=True)
class ScanGeometry:
    """Dimensionless (A,B) coefficients, using canonical ScanConfig field names.

    This measured geometry has no canonical reciprocal HK basis assignment yet;
    reusing the full ScanConfig would require inventing that missing information.
    B corner advances columns; C advances rows. Counts may be singleton here;
    this is a data model, not a declaration of engine execution readiness.
    """
    corner_a: tuple[float, float]
    corner_b: tuple[float, float]
    corner_c: tuple[float, float]
    scan_rows: int
    scan_columns: int

    def __post_init__(self):
        for name in ('corner_a', 'corner_b', 'corner_c'):
            object.__setattr__(self, name, _pair(getattr(self, name), name))
        for name in ('scan_rows', 'scan_columns'):
            value = _integer(getattr(self, name), name)
            if value <= 0:
                raise ValueError(f'{name} must be positive')
            object.__setattr__(self, name, value)


@dataclass(frozen=True, eq=False)
class PreparedExperimentalDataset:
    """Owned read-only float64 arrays in (tilt, beam) order.

    tilt=row*scan_columns+column. Every ordered beam, including disabled beams,
    owns exactly one column; toggling enabled never drops or reorders data. ROI
    order is independent and associations use stable IDs (multiple ROIs per beam
    are allowed). Variances are nonnegative squared-intensity values, not sigma;
    weights are independent finite nonnegative values, including zero. Signed
    intensities are valid. No weights are inferred and no normalization occurs.

    Optional detector_response describes calibration availability only. Whether
    correction, smoothing or baseline subtraction occurred belongs in explicit
    provenance; neither its presence nor absence implies a processing policy.
    Missing scientific variance is represented by None, never inferred from support.
    Refinement readiness is checked separately by the adapter.
    """
    intensities: NDArray[np.float64]
    variances: NDArray[np.float64] | None
    weights: NDArray[np.float64]
    beams: tuple[BeamObservation, ...]
    rois: tuple[ReflectionROI, ...]
    reciprocal_frame: ReciprocalFrame
    scan_geometry: ScanGeometry
    mode: SignalMode
    intensity_units: str
    detector_response: DetectorResponseModel | None = None
    noise_model: NoiseModel | None = None
    provenance: tuple[str, ...] = ()
    physical_disk_geometry: DiskGeometry | None = None
    extraction_config: 'ExtractionConfig | None' = None
    weight_mode: 'WeightMode | None' = None
    lattice_coordinates: tuple[tuple[int, int], ...] = ()
    center_sources: tuple[str, ...] = ()
    extraction_block_indices: tuple[int, ...] = ()

    def __post_init__(self):
        if not isinstance(self.reciprocal_frame, ReciprocalFrame):
            raise TypeError('reciprocal_frame must be ReciprocalFrame')
        if not isinstance(self.scan_geometry, ScanGeometry):
            raise TypeError('scan_geometry must be ScanGeometry')
        if not isinstance(self.mode, SignalMode):
            raise TypeError('mode must be SignalMode')
        _text(self.intensity_units, 'intensity_units')
        for name, cls in (('detector_response', DetectorResponseModel), ('noise_model', NoiseModel)):
            value = getattr(self, name)
            if value is not None and not isinstance(value, cls):
                raise TypeError(f'{name} must be {cls.__name__}')
        for name, cls in (('beams', BeamObservation), ('rois', ReflectionROI)):
            entries = _tuple(getattr(self, name), name)
            if any(not isinstance(v, cls) for v in entries):
                raise TypeError(f'{name} must contain {cls.__name__} entries')
            if len({v.id for v in entries}) != len(entries):
                raise ValueError(f'{name} must have unique stable IDs')
            object.__setattr__(self, name, entries)
        if not self.beams:
            raise ValueError('prepared dataset requires at least one beam')
        ids = {v.id for v in self.beams}
        if any(v.beam_id not in ids for v in self.rois):
            raise ValueError('ROI refers to an unknown beam ID')
        shape = (self.scan_geometry.scan_rows*self.scan_geometry.scan_columns, len(self.beams))
        for name in ('intensities', 'variances', 'weights'):
            value = getattr(self, name)
            if name == 'variances' and value is None:
                continue
            if not isinstance(value, np.ndarray) or value.shape != shape:
                raise ValueError(f'{name} must have (tilt, beam) shape {shape}')
            value = _array(value, np.float64, name)
            if name != 'intensities' and np.any(value < 0):
                raise ValueError(f'{name} must be nonnegative')
            object.__setattr__(self, name, value)
        object.__setattr__(self, 'provenance', _provenance(self.provenance))

        from .extraction import ExtractionConfig, WeightMode
        for name, cls in (('physical_disk_geometry', DiskGeometry),
                          ('extraction_config', ExtractionConfig), ('weight_mode', WeightMode)):
            value = getattr(self, name)
            if value is not None and not isinstance(value, cls):
                raise TypeError(f'{name} must be {cls.__name__}')
        for name in ('lattice_coordinates', 'center_sources', 'extraction_block_indices'):
            values = _tuple(getattr(self, name), name)
            if values and len(values) != len(self.beams):
                raise ValueError(f'{name} must follow beam order with one entry per beam')
            if name == 'lattice_coordinates':
                parsed = []
                for pair in values:
                    pair = _tuple(pair, name)
                    if len(pair) != 2:
                        raise ValueError('lattice coordinates require (n1,n2)')
                    parsed.append(tuple(_integer(v, name) for v in pair))
                values = tuple(parsed)
            elif name == 'extraction_block_indices':
                values = tuple(_integer(v, name) for v in values)
                if any(v < 0 for v in values) or len(set(values)) != len(values):
                    raise ValueError('Extraction block indices must be unique and nonnegative')
            else:
                for value in values:
                    _text(value, name)
            object.__setattr__(self, name, values)
        if self.extraction_config is not None:
            if self.extraction_config.final_pattern_matching_shape_px != (
                    self.scan_geometry.scan_rows, self.scan_geometry.scan_columns):
                raise ValueError('Extraction final shape must match scan geometry')
