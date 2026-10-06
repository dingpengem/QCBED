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
    max_stage_iterations: int | None = None

    def __post_init__(self):
        super().__post_init__()
        if self.max_stage_iterations is not None and self.max_stage_iterations < 0:
            raise ValueError("Maximum stage iterations must be nonnegative")


@dataclass(frozen=True)
class LegacyExportConfig(_Config):
    """Retained export metadata, not a request to write legacy files."""
    dsmax_inverse_angstrom: float | None = None
    filename: str | None = None
    title_raw: str | None = None
    output_filename_raw: str | None = None
    experimental_filename_raw: str | None = None


@dataclass(frozen=True)
class StructureFileConfig(_Config):
    path: str
    format: str | None = None


@dataclass(frozen=True)
class CbedReflectionBasisConfig(_Config):
    g1_hkl: tuple[int, int, int]
    g2_hkl: tuple[int, int, int]

    def __post_init__(self):
        for values in (self.g1_hkl, self.g2_hkl):
            if isinstance(values, (tuple, list)) and len(values) == 4:
                raise ValueError("4-index Miller–Bravais notation is not supported; provide canonical 3-index input")
        super().__post_init__()


@dataclass(frozen=True)
class ExperimentalBeamConfig(_Config):
    id: str
    lattice_coordinates: tuple[int, int]


@dataclass(frozen=True)
class AtomicSiteGroup(_Config):
    """Explicit physical ownership in the ordered original crystallographic cell.

    Member indices refer to the immutable source structure, never simulation atoms.
    B is isotropic Debye–Waller B in Å², occupancy a chemical fraction, absorption
    the existing independent species (a,b,c) coefficients, not a new potential.
    """
    id: str
    element: str
    site_indices: tuple[int, ...]
    debye_waller_b: float | None
    occupancy: float
    absorption: tuple[float, float, float] | None = None

    def __post_init__(self):
        super().__post_init__()
        from .crystallography.structure import AtomicSite
        if not self.id.strip():
            raise ValueError('site-group ID must be nonempty')
        if not self.site_indices or len(set(self.site_indices)) != len(self.site_indices) or min(self.site_indices) < 0:
            raise ValueError('site-group members must be distinct nonnegative original indices')
        AtomicSite(self.element, (0., 0., 0.), self.occupancy, self.debye_waller_b)
        if self.absorption is not None:
            from .absorption import AbsorptionCoefficients
            AbsorptionCoefficients(*self.absorption)


@dataclass(frozen=True)
class ProjectedElasticConfig(_Config):
    """Physical slice interval and original hkl; a Friedel pair has one owner."""
    hkl: tuple[int, int, int]
    interval: tuple[float, float]
    mode: Literal['real', 'complex', 'magnitude'] = 'real'
    refine: bool = False
    initial: tuple[float, float] | None = None
    reference: str | None = None
    lower: tuple[float | None, float | None] = (None, None)
    upper: tuple[float | None, float | None] = (None, None)
    step_size: float = 0.01
    equality_group: str | None = None

    def __post_init__(self):
        super().__post_init__()
        if self.hkl == (0, 0, 0):
            raise ValueError('Elastic (000) is read-only: uniform phase is not an intensity target')
        if next(v for v in self.hkl if v) < 0:
            raise ValueError('Use the Friedel owner with first nonzero hkl component positive')
        if not 0 <= self.interval[0] < self.interval[1] <= 1:
            raise ValueError('Require a physical slice interval within [0, 1]')
        if self.equality_group is not None and not self.equality_group.strip():
            raise ValueError('Equality group name must not be empty')
        if self.step_size <= 0:
            raise ValueError('Projected elastic search step must be positive')


@dataclass(frozen=True)
class CrystallographyConfig(_Config):
    structure: StructureFileConfig
    zone_axis: tuple[int, int, int]
    cbed_basis: CbedReflectionBasisConfig
    mesh: MeshConfig
    depth_angstrom: float
    slice_boundaries_fractional: tuple[float, ...]
    debye_waller: tuple[ThermalFactorConfig, ...]
    dsap_inverse_angstrom: float
    beams: tuple[ExperimentalBeamConfig, ...]
    scan_basis_transform: tuple[tuple[int, int], tuple[int, int]]
    absorption: AbsorptionConfig = AbsorptionConfig()
    site_groups: tuple[AtomicSiteGroup, ...] = ()
    projected_elastic: tuple[ProjectedElasticConfig, ...] = ()


@dataclass(frozen=True)
class CrystalScanConfig(_Config):
    scan_columns: int
    scan_rows: int
    corner_a: tuple[float, float]
    corner_b: tuple[float, float]
    corner_c: tuple[float, float]


@dataclass(frozen=True)
class CrystalExperimentalDataConfig(_Config):
    path: str
    gradient_option: int
    header: Literal['absent', 'validate'] = 'absent'
    format: Literal['legacy_intensity_blocks'] = 'legacy_intensity_blocks'


@dataclass(frozen=True)
class QcbedConfig(_Config):
    title: str
    microscope: MicroscopeConfig
    sample: SampleConfig | None
    scan: ScanConfig | CrystalScanConfig
    experiment: ExperimentalDataConfig | CrystalExperimentalDataConfig
    refinement: RefinementConfig
    execution: ForwardExecutionConfig
    run: RunConfig
    legacy_export: LegacyExportConfig = LegacyExportConfig()
    schema_version: int = 1
    crystallography: CrystallographyConfig | None = None
