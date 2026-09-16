"""Prepare fixed S-atom samples and simulate one tilt in the Python convention."""

from dataclasses import dataclass, field
from typing import Literal

import numpy as np

from .absorption import AbsorptionCoefficients, _positive, absorptive_structure_factors
from .electron import electron_wavelength_angstrom
from .fft_backend import ComplexGrid, Ordering, _array, _ordering
from .legacy_io import AtomSite, LegacyInput, Species
from .models import MeshConfig, SliceGeometry
from .multislice import (
    MultisliceResult, _beams, _items, _sequence, _snapshot,
    beam_indices_to_row_col, simulate_slice_sequence,
)
from .phase_grating import phase_grating_from_structure_factors, legacy_nyquist_mask
from .propagator import PropagationParameters, propagation_kernel
from .reciprocal import DirectLattice2D, ReciprocalLattice2D, ReciprocalIndexGrid, ReciprocalVectorGrid
from .scattering import elastic_atomic_scattering_factor
from .structure_factors import _finite, elastic_structure_factors


@dataclass(frozen=True)
class SliceAppliedParameters:
    """Already-applied APPLIEDDW (one B per species), VAC and APPLIEDDILATION.

    B is in angstrom squared; vacancy and dilation are dimensionless. No
    vacancy/surface model, clipping or refinement projection is performed.
    """

    debye_waller_angstrom_squared: tuple[float, ...]
    vacancy_fraction: float
    dilation: float

    def __post_init__(self) -> None:
        thermal = tuple(_finite(value, f"debye_waller_angstrom_squared[{i}]", nonnegative=True)
                        for i, value in enumerate(_items(self.debye_waller_angstrom_squared,
                                                        "debye_waller_angstrom_squared")))
        vacancy = _finite(self.vacancy_fraction, "vacancy_fraction")
        if not 0 <= vacancy <= 1:
            raise ValueError("vacancy_fraction must be in [0, 1]")
        object.__setattr__(self, "debye_waller_angstrom_squared", thermal)
        object.__setattr__(self, "vacancy_fraction", vacancy)
        object.__setattr__(self, "dilation", _positive(self.dilation, "dilation"))


@dataclass(frozen=True)
class SamplePreparationParameters:
    """Explicit state not supplied by the current LegacyInput parser model.

    Voltage, DSAP and species absorption weights are global. Applied thermal,
    vacancy and dilation state is per stored slice type. All fields are
    required: None coefficients explicitly means absorption is disabled.
    Only unity SFFAC/ABSFAC and the established legacy boundary are supported.
    """

    applied_voltage_kv: float
    slices: tuple[SliceAppliedParameters, ...]
    dsap_inverse_angstrom: float
    absorption_coefficients: tuple[AbsorptionCoefficients, ...] | None
    slice_sequence_zero_based: tuple[int, ...]
    beam_indices_hk: tuple[tuple[int, int], ...]
    reciprocal_ordering: Ordering
    phase_grating_method: Literal["direct", "fft"]
    boundary_policy: Literal["legacy"]
    reflection_multipliers: Literal["unity"]

    def __post_init__(self) -> None:
        _ordering(self.reciprocal_ordering)
        if not isinstance(self.phase_grating_method, str) or self.phase_grating_method not in ("direct", "fft"):
            raise ValueError("phase_grating_method must be 'direct' or 'fft'")
        if not isinstance(self.boundary_policy, str) or self.boundary_policy != "legacy":
            raise ValueError("boundary_policy must be 'legacy'; preparation uses the source Nyquist mask")
        if not isinstance(self.reflection_multipliers, str) or self.reflection_multipliers != "unity":
            raise ValueError("reflection_multipliers must be 'unity'; nonunit SFFAC/ABSFAC are unsupported")
        slices = _items(self.slices, "slices")
        if not slices:
            raise ValueError("slices must contain at least one stored slice type")
        for i, value in enumerate(slices):
            if not isinstance(value, SliceAppliedParameters):
                raise TypeError(f"slices[{i}] must be SliceAppliedParameters")
        coefficients = self.absorption_coefficients
        if coefficients is not None:
            coefficients = _items(coefficients, "absorption_coefficients")
            for i, value in enumerate(coefficients):
                if not isinstance(value, AbsorptionCoefficients):
                    raise TypeError(f"absorption_coefficients[{i}] must be AbsorptionCoefficients")
        object.__setattr__(self, "slices", slices)
        object.__setattr__(self, "absorption_coefficients", coefficients)
        object.__setattr__(self, "applied_voltage_kv", _positive(self.applied_voltage_kv, "applied_voltage_kv"))
        object.__setattr__(self, "dsap_inverse_angstrom", _finite(self.dsap_inverse_angstrom,
                                                               "dsap_inverse_angstrom", nonnegative=True))
        object.__setattr__(self, "slice_sequence_zero_based",
                           _sequence(self.slice_sequence_zero_based, len(slices), one_based=False))
        object.__setattr__(self, "beam_indices_hk", _beams(self.beam_indices_hk))


def _cell(geometry: SliceGeometry) -> DirectLattice2D:
    if not isinstance(geometry, SliceGeometry):
        raise TypeError("geometry must be SliceGeometry")
    return DirectLattice2D(geometry.cell_a_angstrom, geometry.cell_b_angstrom,
                           geometry.lattice_gamma_degrees)


def _geometry_parameters(geometry: SliceGeometry, wavelength: float, **tilts) -> PropagationParameters:
    if geometry.thickness_angstrom <= 0:
        raise ValueError("slice thickness_angstrom must be positive")
    return PropagationParameters(wavelength, geometry.thickness_angstrom,
                                 x_shift_fractional=geometry.x_shift_fractional,
                                 y_shift_fractional=geometry.y_shift_fractional, **tilts)


@dataclass(frozen=True, eq=False)
class PreparedSliceType:
    """Reusable geometry and immutable real-space grating; no incident tilt."""

    geometry: SliceGeometry
    reciprocal_grid: ReciprocalVectorGrid
    phase_grating: ComplexGrid = field(repr=False)

    def __post_init__(self) -> None:
        if not isinstance(self.reciprocal_grid, ReciprocalVectorGrid):
            raise TypeError("reciprocal_grid must be ReciprocalVectorGrid")
        if _cell(self.geometry) != self.reciprocal_grid.lattice.direct:
            raise ValueError("geometry cell must match reciprocal_grid")
        values = _array(self.phase_grating, "phase_grating")
        if values.shape != self.reciprocal_grid.indices.h.shape:
            raise ValueError("phase_grating shape must match reciprocal_grid")
        object.__setattr__(self, "phase_grating", _snapshot(values))


@dataclass(frozen=True, eq=False)
class PreparedMultisliceSample:
    """Immutable cache boundary for the supported applied-state snapshot."""

    mesh: MeshConfig
    slice_types: tuple[PreparedSliceType, ...]
    parameters: SamplePreparationParameters
    absorption_raw: str
    source_title_raw: str
    wavelength_angstrom: float = field(init=False)
    beam_positions_row_col: tuple[tuple[int, int], ...] = field(init=False)

    def __post_init__(self) -> None:
        if not isinstance(self.mesh, MeshConfig):
            raise TypeError("mesh must be MeshConfig")
        if not isinstance(self.parameters, SamplePreparationParameters):
            raise TypeError("parameters must be SamplePreparationParameters")
        if not isinstance(self.absorption_raw, str) or not isinstance(self.source_title_raw, str):
            raise TypeError("absorption_raw and source_title_raw must be str")
        types = _items(self.slice_types, "slice_types")
        if len(types) != len(self.parameters.slices):
            raise ValueError("slice_types count must match parameters.slices")
        if (self.absorption_raw[:2] == "YE") != (self.parameters.absorption_coefficients is not None):
            raise ValueError("absorption_coefficients must agree with the raw absorption flag")
        wavelength = _positive(electron_wavelength_angstrom(self.parameters.applied_voltage_kv),
                               "wavelength_angstrom")
        for i, value in enumerate(types):
            try:
                if not isinstance(value, PreparedSliceType):
                    raise TypeError("must be PreparedSliceType")
                grid = value.reciprocal_grid
                if (grid.indices.mesh_x, grid.indices.mesh_y) != (self.mesh.mesh_x, self.mesh.mesh_y):
                    raise ValueError("mesh shape must match sample mesh")
                if grid.indices.centered != (self.parameters.reciprocal_ordering == "centered"):
                    raise ValueError("reciprocal ordering must match sample")
                if grid.dilation != self.parameters.slices[i].dilation:
                    raise ValueError("grid dilation must match applied slice parameters")
                if grid.lattice.direct != types[0].reciprocal_grid.lattice.direct:
                    raise ValueError("undilated lateral cell must match slice type 0")
                legacy_nyquist_mask(grid)
                _geometry_parameters(value.geometry, wavelength)
            except (ValueError, TypeError) as exc:
                raise type(exc)(f"slice type {i}: {exc}") from exc
        positions = beam_indices_to_row_col(self.parameters.beam_indices_hk, mesh_x=self.mesh.mesh_x,
                                            mesh_y=self.mesh.mesh_y, ordering=self.parameters.reciprocal_ordering)
        object.__setattr__(self, "slice_types", types)
        object.__setattr__(self, "wavelength_angstrom", wavelength)
        object.__setattr__(self, "beam_positions_row_col", positions)

    @property
    def absorption_enabled(self) -> bool:
        return self.absorption_raw[:2] == "YE"


def prepare_multislice_sample(
    model: LegacyInput, parameters: SamplePreparationParameters,
) -> PreparedMultisliceSample:
    """Prepare each stored S-atom slice once, validating all types before sums.

    The supported lateral contract is exactly equal undilated a,b,gamma across
    slices and one shared even mesh. Source per-slice dilation remains allowed.
    Applied voltage comes from parameters, not an inferred refinement state.
    """
    if not isinstance(model, LegacyInput):
        raise TypeError("model must be LegacyInput")
    if not isinstance(parameters, SamplePreparationParameters):
        raise TypeError("parameters must be SamplePreparationParameters")
    if not isinstance(model.mesh, MeshConfig):
        raise TypeError("model.mesh must be MeshConfig")
    _positive(model.energy_kv, "model.energy_kv")
    wavelength = _positive(electron_wavelength_angstrom(parameters.applied_voltage_kv), "wavelength_angstrom")
    if not isinstance(model.absorption_raw, str):
        raise TypeError("model.absorption_raw must be str")
    absorption = model.absorption_raw[:2] == "YE"
    if absorption != (parameters.absorption_coefficients is not None):
        raise ValueError("absorption_coefficients must be supplied when absorption_raw starts YE, otherwise explicitly None")
    geometries = _items(model.slices, "model.slices")
    atom_groups = tuple(_items(atoms, f"atoms_by_slice[{i}]")
                        for i, atoms in enumerate(_items(model.atoms_by_slice, "atoms_by_slice")))
    species = _items(model.species, "species")
    if not species:
        raise ValueError("species must not be empty")
    if len(geometries) != len(parameters.slices) or len(atom_groups) != len(geometries):
        raise ValueError("stored slice type counts must match geometry, atoms_by_slice and parameters.slices")
    if absorption and len(parameters.absorption_coefficients) != len(species):
        raise ValueError("absorption_coefficients requires one entry per global species")
    for i, entry in enumerate(species):
        if not isinstance(entry, Species):
            raise TypeError(f"species[{i}] must be Species")
        try:
            elastic_atomic_scattering_factor(entry.atomic_number, 0.0)
        except (ValueError, TypeError) as exc:
            raise type(exc)(f"species[{i}]: {exc}") from exc
    beam_indices_to_row_col(parameters.beam_indices_hk, mesh_x=model.mesh.mesh_x,
                            mesh_y=model.mesh.mesh_y, ordering=parameters.reciprocal_ordering)
    indices = ReciprocalIndexGrid(model.mesh.mesh_x, model.mesh.mesh_y,
                                   centered=parameters.reciprocal_ordering == "centered")
    grids = []
    # Complete cross-slice preflight before any expensive structure-factor sum.
    for i, (geometry, atoms, applied) in enumerate(zip(geometries, atom_groups, parameters.slices)):
        try:
            cell = _cell(geometry)
            if grids and cell != grids[0].lattice.direct:
                raise ValueError("undilated lateral cell must match slice type 0; resampling is unsupported")
            _geometry_parameters(geometry, wavelength)
            grid = ReciprocalVectorGrid(ReciprocalLattice2D(cell), indices, applied.dilation)
            legacy_nyquist_mask(grid)
            if len(applied.debye_waller_angstrom_squared) != len(species):
                raise ValueError("debye_waller_angstrom_squared requires one value per species")
            for j, b in enumerate(applied.debye_waller_angstrom_squared):
                if absorption and b <= 0:
                    raise ValueError(f"debye_waller_angstrom_squared[{j}] must be positive for absorption")
            for j, atom in enumerate(atoms):
                if not isinstance(atom, AtomSite):
                    raise TypeError(f"atom {j} must be AtomSite")
                if atom.record_type != "S":
                    raise ValueError(f"atom {j} record_type must be S; G randomization is unsupported")
                atom.species_index_zero_based(len(species))
                _finite(atom.x_fractional, f"atom {j} x_fractional")
                _finite(atom.y_fractional, f"atom {j} y_fractional")
                _finite(atom.starting_occupancy, f"atom {j} starting_occupancy", nonnegative=True)
            # Raw lower-level sums evaluate the whole grid before masking.
            max_stl = float(grid.stl_inverse_angstrom.max())
            for j, entry in enumerate(species):
                try:
                    elastic_atomic_scattering_factor(entry.atomic_number, max_stl)
                except (ValueError, TypeError) as exc:
                    raise type(exc)(f"species[{j}] full-grid scattering domain: {exc}") from exc
            grids.append(grid)
        except (ValueError, TypeError, OverflowError) as exc:
            raise ValueError(f"slice type {i}: {exc}") from exc
    prepared = []
    for i, (geometry, atoms, grid, applied) in enumerate(zip(geometries, atom_groups, grids, parameters.slices)):
        try:
            elastic = elastic_structure_factors(species, atoms, grid,
                applied_debye_waller_angstrom_squared=applied.debye_waller_angstrom_squared,
                vacancy_fraction=applied.vacancy_fraction)
            absorptive = (absorptive_structure_factors(species, atoms, grid,
                voltage_kv=parameters.applied_voltage_kv,
                applied_debye_waller_angstrom_squared=applied.debye_waller_angstrom_squared,
                coefficients=parameters.absorption_coefficients,
                vacancy_fraction=applied.vacancy_fraction, absorption_raw=model.absorption_raw)
                if absorption else np.zeros_like(elastic))
            grating = phase_grating_from_structure_factors(elastic, absorptive, grid,
                voltage_kv=parameters.applied_voltage_kv, dsap_inverse_angstrom=parameters.dsap_inverse_angstrom,
                method=parameters.phase_grating_method)
            prepared.append(PreparedSliceType(geometry, grid, grating))
        except (ValueError, TypeError, OverflowError) as exc:
            raise ValueError(f"slice type {i}: {exc}") from exc
    return PreparedMultisliceSample(model.mesh, tuple(prepared), parameters,
                                    model.absorption_raw, model.title_raw)


def build_tilt_propagators(
    prepared: PreparedMultisliceSample, *,
    tilt_h_undilated_fractional: float, tilt_k_undilated_fractional: float,
) -> tuple[ComplexGrid, ...]:
    """One kernel per stored type; divide common cell-basis tilt by its dilation.

    Inputs are dimensionless coefficients after scan-basis transformation but
    before per-slice dilation division, not radians or mrad. Linearity of the
    source caller's basis transform permits this ordering of the division.
    """
    if not isinstance(prepared, PreparedMultisliceSample):
        raise TypeError("prepared must be PreparedMultisliceSample")
    h = _finite(tilt_h_undilated_fractional, "tilt_h_undilated_fractional")
    k = _finite(tilt_k_undilated_fractional, "tilt_k_undilated_fractional")
    kernels = []
    for i, stored in enumerate(prepared.slice_types):
        try:
            d = stored.reciprocal_grid.dilation
            parameters = _geometry_parameters(stored.geometry, prepared.wavelength_angstrom,
                tilt_h_fractional=h/d, tilt_k_fractional=k/d)
            kernels.append(_snapshot(propagation_kernel(stored.reciprocal_grid, parameters,
                ordering=prepared.parameters.reciprocal_ordering, boundary=prepared.parameters.boundary_policy)))
        except (ValueError, TypeError, OverflowError) as exc:
            raise ValueError(f"slice type {i}: {exc}") from exc
    return tuple(kernels)


def simulate_one_tilt(
    prepared: PreparedMultisliceSample, *,
    tilt_h_undilated_fractional: float, tilt_k_undilated_fractional: float,
) -> MultisliceResult:
    """Construct current-tilt propagators and execute with unit illumination."""
    kernels = build_tilt_propagators(prepared, tilt_h_undilated_fractional=tilt_h_undilated_fractional,
                                    tilt_k_undilated_fractional=tilt_k_undilated_fractional)
    parameters = prepared.parameters
    return simulate_slice_sequence(tuple(s.phase_grating for s in prepared.slice_types), kernels,
        parameters.slice_sequence_zero_based, beam_indices_hk=parameters.beam_indices_hk,
        reciprocal_ordering=parameters.reciprocal_ordering, propagator_ordering=parameters.reciprocal_ordering)
