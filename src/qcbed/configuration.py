"""Immutable user configuration, independent of file formats and numerical state.

Lengths are angstroms, thermal factors angstrom squared, energy kV, lattice
angles degrees. Scan coordinates and slice shifts are dimensionless. Loading,
physical mapping and execution are deliberately outside these declarations.
"""
from __future__ import annotations

from dataclasses import dataclass, fields
from enum import IntEnum
from math import isfinite
from numbers import Integral, Real
from types import UnionType
from typing import Literal, get_args, get_origin, get_type_hints

from .models import MeshConfig, SliceGeometry
from .execution import ForwardExecutionConfig
from .refinement_categories import ParameterCategory


class VacancyKind(IntEnum):
    CONSTANT = 0
    LEGACY_PROFILE_1 = 1
    LEGACY_PROFILE_2 = 2


class SurfaceKind(IntEnum):
    NONE = 0
    LEGACY_1 = 1
    LEGACY_2 = 2


def _typed(value, annotation, path):
    """Validate declared scalar types and defensively freeze collection inputs."""
    origin, args = get_origin(annotation), get_args(annotation)
    if origin is UnionType:
        for option in args:
            try:
                return _typed(value, option, path)
            except (TypeError, ValueError):
                pass
        raise TypeError(f'{path}: expected {annotation}')
    if origin is tuple:
        if not isinstance(value, (tuple, list)):
            raise TypeError(f'{path}: expected a sequence')
        if len(args) == 2 and args[1] is Ellipsis:
            types = (args[0],)*len(value)
        else:
            types = args
            if len(value) != len(types):
                raise ValueError(f'{path}: expected {len(types)} entries')
        return tuple(_typed(v, t, f'{path}[{i}]') for i, (v, t) in enumerate(zip(value, types)))
    if origin is Literal:
        if value not in args or not any(type(value) is type(v) for v in args):
            raise ValueError(f'{path}: expected one of {args}')
        return value
    if annotation is float:
        if isinstance(value, bool) or not isinstance(value, Real) or not isfinite(value):
            raise ValueError(f'{path}: expected a finite real number excluding bool')
        return float(value)
    if annotation is int:
        if isinstance(value, bool) or not isinstance(value, Integral):
            raise TypeError(f'{path}: expected an integer excluding bool')
        return int(value)
    if not isinstance(value, annotation):
        raise TypeError(f'{path}: expected {annotation.__name__}')
    return value


class _Config:
    def __post_init__(self):
        hints = get_type_hints(type(self))
        for field in fields(self):
            object.__setattr__(self, field.name, _typed(
                getattr(self, field.name), hints[field.name], f'{type(self).__name__}.{field.name}'))


@dataclass(frozen=True)
class MicroscopeConfig(_Config):
    energy_kv: float
    zone_axis: tuple[int, int, int]


@dataclass(frozen=True)
class SpeciesConfig(_Config):
    id: str
    label: str
    atomic_number: int
    label_raw: str = ''


@dataclass(frozen=True)
class AtomConfig(_Config):
    id: str
    species: str
    x_fractional: float
    y_fractional: float
    starting_occupancy: float


@dataclass(frozen=True)
class StoredSliceConfig(_Config):
    id: str
    geometry: SliceGeometry
    atoms: tuple[AtomConfig, ...]
    record_type: Literal['S', 'G'] = 'S'
    record_type_raw: str = ''


@dataclass(frozen=True)
class VacancyConfig(_Config):
    kind: VacancyKind = VacancyKind.CONSTANT
    values: tuple[float, ...] = (0.,)


@dataclass(frozen=True)
class SurfaceConfig(_Config):
    kind: SurfaceKind = SurfaceKind.NONE
    coefficients: tuple[float, float, float] = (0., 0., 1.)


@dataclass(frozen=True)
class ThermalFactorConfig(_Config):
    species: str
    value_angstrom_squared: float


@dataclass(frozen=True)
class ReflectionConfig(_Config):
    id: str
    h: int
    k: int
    elastic: float
    absorptive: float


@dataclass(frozen=True)
class RegionConfig(_Config):
    id: str
    depth_angstrom: float
    slices: tuple[StoredSliceConfig, ...]
    debye_waller: tuple[ThermalFactorConfig, ...]
    reflections: tuple[ReflectionConfig, ...] = ()
    symmetry: Literal['NONE', 'INVE', 'SQUA', 'RECT', 'HEXA'] = 'NONE'
    vacancy: VacancyConfig = VacancyConfig()
    free_surface: SurfaceConfig = SurfaceConfig()


@dataclass(frozen=True)
class AbsorptionSpeciesConfig(_Config):
    species: str
    a: float
    b: float
    c: float
    c_mode: Literal['independent', 'linked_to_b'] = 'independent'


@dataclass(frozen=True)
class AbsorptionConfig(_Config):
    enabled: bool = False
    coefficients: tuple[AbsorptionSpeciesConfig, ...] = ()
    legacy_raw: str | None = None


@dataclass(frozen=True)
class SampleConfig(_Config):
    mesh: MeshConfig
    species: tuple[SpeciesConfig, ...]
    regions: tuple[RegionConfig, ...]
    dsap_inverse_angstrom: float
    dilation: float = 0.
    absorption: AbsorptionConfig = AbsorptionConfig()


@dataclass(frozen=True)
class ScanConfig(_Config):
    scan_columns: int
    scan_rows: int
    basis_1: tuple[float, float, float]
    basis_2: tuple[float, float, float]
    corner_a: tuple[float, float]
    corner_b: tuple[float, float]
    corner_c: tuple[float, float]


@dataclass(frozen=True)
class ExperimentalDataConfig(_Config):
    """External data specification; loading is separate from compilation."""
    path: str
    beam_indices_hk: tuple[tuple[int, int], ...]
    gradient_option: int
    header: Literal['absent', 'validate'] = 'absent'
    format: Literal['legacy_intensity_blocks'] = 'legacy_intensity_blocks'


@dataclass(frozen=True)
class ParameterTarget(_Config):
    """Stable IDs/components, resolved to source indices by compile_config."""
    category: ParameterCategory
    identifiers: tuple[str, ...] = ()


@dataclass(frozen=True)
class ParameterBinding(_Config):
    target: ParameterTarget
    range_factors: tuple[float, float]
    plamda: float


@dataclass(frozen=True)
class OptimizerConfig(_Config):
    method: str
    loops: int
    maxfsimp: int
    maxfpara: int
    ftol: float
    chimin: float
    maxtemp: float
    mintemp: float
    coolrate: float
    seed: int
    paraopt: int = 0
    loud: int = 0


@dataclass(frozen=True)
class TargetBound(_Config):
    target: ParameterTarget
    lower: float | None = None
    upper: float | None = None


@dataclass(frozen=True)
class CategoryBound(_Config):
    category: ParameterCategory
    lower: float | None = None
    upper: float | None = None


@dataclass(frozen=True)
class RelativeBound(_Config):
    category: ParameterCategory
    fraction: float
    absolute_fallback: float


@dataclass(frozen=True)
class BoundsConfig(_Config):
    policy: Literal['none', 'reject'] = 'none'
    specific: tuple[TargetBound, ...] = ()
    absolute: tuple[CategoryBound, ...] = ()
    relative: tuple[RelativeBound, ...] = ()


@dataclass(frozen=True)
class RefinementConfig(_Config):
    optimizer: OptimizerConfig
    parameters: tuple[ParameterBinding, ...]
    iref: float
    objective_mode: Literal['legacy', 'masked'] = 'legacy'
    invalid_trial_policy: Literal['raise', 'reject'] = 'raise'
    active_categories: tuple[ParameterCategory, ...] | None = None
    bounds: BoundsConfig = BoundsConfig()


@dataclass(frozen=True)
class RunConfig(_Config):
    output_directory: str
    output_filename: str = 'out.dat'
    snapshot_every: int | None = 100
    checkpoint_every: int = 25
    max_total_evaluations: int | None = 1000
    max_wall_time_seconds: float | None = 3600.


@dataclass(frozen=True)
class LegacyExportConfig(_Config):
    """Retained export metadata, not a request to write legacy files."""
    dsmax_inverse_angstrom: float | None = None
    filename: str | None = None
    title_raw: str | None = None
    output_filename_raw: str | None = None
    experimental_filename_raw: str | None = None


@dataclass(frozen=True)
class QcbedConfig(_Config):
    title: str
    microscope: MicroscopeConfig
    sample: SampleConfig
    scan: ScanConfig
    experiment: ExperimentalDataConfig
    refinement: RefinementConfig
    execution: ForwardExecutionConfig
    run: RunConfig
    legacy_export: LegacyExportConfig = LegacyExportConfig()
    schema_version: int = 1
