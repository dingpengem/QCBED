"""Supported source-backed trial assembly; no optimizer or mutable best state."""

from dataclasses import dataclass

import numpy as np

from .absorption import AbsorptionCoefficients, absorptive_structure_factors
from .experimental import ExperimentalIntensityData
from .multislice import _snapshot, simulate_slice_sequence
from .masked_objective import compute_selected_objective
from .phase_grating import phase_grating_from_structure_factors
from .refinement import RefinementState, apply_refinement_parameters
from .scan import generate_scan_tilts
from .simulation import (SamplePreparationParameters, SliceAppliedParameters,
                         prepare_multislice_sample, build_tilt_propagators)
from .structure_factors import elastic_structure_factors


def _reflection_scales(indices, reflections, symmetry):
    """SFMAP2: all exact matches before symmetry, first record wins."""
    elastic = np.ones(indices.h.shape, dtype=np.float64)
    absorptive = np.ones_like(elastic)
    for pos in np.ndindex(elastic.shape):
        h, k = int(indices.h[pos]), int(indices.k[pos])
        found = next((r for r in reflections if (r.h, r.k) == (h, k)), None)
        if found is None:
            transforms = {
                "NONE": (), "INVE": ((-h, -k),),
                "RECT": ((-h, k), (h, -k), (-h, -k)),
                "SQUA": ((k, h), (-h, k), (h, -k), (-h, -k), (-k, -h), (-k, h), (k, -h)),
                "HEXA": ((k, -(h+k)), (-(h+k), h), (k, h), (h, -(h+k)), (-(h+k), k),
                         (-h, -k), (-k, h+k), (h+k, -h), (-k, -h), (-h, h+k), (h+k, -k)),
            }[symmetry]
            found = next((r for r in reflections if (r.h, r.k) in transforms), None)
        if found is not None:
            elastic[pos], absorptive[pos] = found.elastic, found.absorptive
    return elastic, absorptive


def _prepare_trial(state, *, original_voltage):
    # In this subset the three source UPDPARS/MSEQ passes are identical:
    # VAC=0, APPLIEDDILATION=1, and no free-surface correction.
    for i, r in enumerate(state.regions):
        if r.vacancy_flag != 0 or r.vacancy_values != (0.0,) or r.free_surface[0] != 0:
            raise ValueError(f"region {i}: vacancy/free-surface derived-state hazards are unsupported for trial evaluation")
        if r.depth_angstrom <= 0:
            raise ValueError(f"region {i}: positive depth required by the supported forward sequence")
    if any(atom.record_type != "S" for row in state.model.atoms_by_slice for atom in row):
        raise ValueError("G-atom randomization is unsupported for deterministic trials")
    if state.model.absorption_enabled and state.model.energy_kv != original_voltage:
        raise ValueError("voltage refinement with absorption is unsupported: UPDPARS leaves VOLTS stale")
    applied, sequence, region_ids = [], [], []
    offset = 0
    for i, r in enumerate(state.regions):
        applied.extend(SliceAppliedParameters(r.debye_waller, 0.0, 1.0) for _ in range(r.slice_count))
        region_ids.extend([i]*r.slice_count)
        thickness, occurrence = 0.0, 0
        while thickness < r.depth_angstrom:
            index = offset + occurrence % r.slice_count
            step = state.model.slices[index].thickness_angstrom
            if step <= 0 or thickness+step == thickness:
                raise ValueError("slice thickness must advance the source sequence")
            sequence.append(index)
            thickness += step
            occurrence += 1
        offset += r.slice_count
    coefficients = tuple(AbsorptionCoefficients(*row) for row in state.absorption) if state.model.absorption_enabled else None
    # Beam order supplied by the experimental dataset at the public boundary.
    return tuple(applied), tuple(sequence), tuple(region_ids), coefficients


def _gradient(values, rows, columns):
    """gradients.f:13–48: zero-padded eight-neighbor mean minus center."""
    images = values.reshape(rows, columns, -1)
    p = np.pad(images, ((1, 1), (1, 1), (0, 0)))
    # Match literal source addition order; axes here are row, column, beam.
    total = p[1:-1, 2:] + p[2:, 1:-1]
    total = total+p[1:-1, :-2]
    total = total+p[:-2, 1:-1]
    total = total+p[2:, 2:]
    total = total+p[:-2, 2:]
    total = total+p[2:, :-2]
    total = total+p[:-2, :-2]
    return (total/8-images).reshape(values.shape)


@dataclass(frozen=True, eq=False)
class RefinementTrialResult:
    objective: np.float64
    scaling: np.float64
    physical_state: RefinementState
    calculated_intensities: np.ndarray  # before GRADIENT
    objective_intensities: np.ndarray  # after optional GRADIENT, before YNORM
    amplitudes: np.ndarray
    objective_mode: str = "legacy"

    def __post_init__(self):
        for name in ("calculated_intensities", "objective_intensities", "amplitudes"):
            object.__setattr__(self, name, _snapshot(getattr(self, name)))


def evaluate_refinement_trial(
    state: RefinementState, trial_vector: np.ndarray, experimental: ExperimentalIntensityData, *,
    coordinate_mode="physical", phase_grating_method="fft", objective_mode="legacy",
    evaluation_callback=None,
) -> RefinementTrialResult:
    """Map a fresh trial, rebuild physics, run serial tilts, compute YNORM/objective.

    Supports S atoms, zero vacancy and no surface correction; rejects unresolved
    physical updates explicitly. All source reflection symmetries are supported.
    The objective API computes/applies YNORM exactly once. No best-state updates.
    """
    if not isinstance(experimental, ExperimentalIntensityData):
        raise TypeError("experimental must be ExperimentalIntensityData")
    if objective_mode not in ("legacy", "masked"):
        raise ValueError("objective_mode must be 'legacy' or 'masked'")
    mapped = apply_refinement_parameters(state, trial_vector, coordinate_mode=coordinate_mode)
    scan = mapped.scan
    if (scan.scan_rows, scan.scan_columns) != (experimental.scan_rows, experimental.scan_columns):
        raise ValueError("experimental scan dimensions must match trial scan")
    applied, sequence, regions, coefficients = _prepare_trial(
        mapped, original_voltage=state.model.energy_kv)
    parameters = SamplePreparationParameters(mapped.model.energy_kv, applied, mapped.dsap,
        coefficients, sequence, experimental.beam_indices_hk, "centered", phase_grating_method, "legacy", "unity")
    base = prepare_multislice_sample(mapped.model, parameters)
    gratings = []
    for i, stored in enumerate(base.slice_types):
        r = regions[i]
        reflections = mapped.reflections[r]
        if all(v.elastic == 1 and v.absorptive == 1 for v in reflections):
            gratings.append(stored.phase_grating)
            continue
        grid = stored.reciprocal_grid
        elastic = elastic_structure_factors(mapped.model.species, mapped.model.atoms_by_slice[i], grid,
            applied_debye_waller_angstrom_squared=applied[i].debye_waller_angstrom_squared, vacancy_fraction=0.0)
        absorptive = (absorptive_structure_factors(mapped.model.species, mapped.model.atoms_by_slice[i], grid,
            voltage_kv=mapped.model.energy_kv, applied_debye_waller_angstrom_squared=applied[i].debye_waller_angstrom_squared,
            coefficients=coefficients, vacancy_fraction=0.0, absorption_raw=mapped.model.absorption_raw)
            if coefficients is not None else np.zeros_like(elastic))
        es, ab = _reflection_scales(grid.indices, reflections, mapped.symmetry[r])
        gratings.append(phase_grating_from_structure_factors(elastic*es, absorptive*ab, grid,
            voltage_kv=mapped.model.energy_kv, dsap_inverse_angstrom=mapped.dsap, method=phase_grating_method))
    tilts = generate_scan_tilts(scan)
    amplitudes = np.empty((len(tilts), len(experimental.beam_indices_hk)), dtype=np.complex128)
    intensities = np.empty(amplitudes.shape, dtype=np.float64)
    for i, (h, k) in enumerate(tilts):
        kernels = build_tilt_propagators(base, tilt_h_undilated_fractional=h, tilt_k_undilated_fractional=k)
        result = simulate_slice_sequence(tuple(gratings), kernels, sequence,
            beam_indices_hk=experimental.beam_indices_hk, reciprocal_ordering="centered", propagator_ordering="centered")
        amplitudes[i], intensities[i] = result.amplitudes, result.intensities
    objective_values = _gradient(intensities, scan.scan_rows, scan.scan_columns) if experimental.gradient_option == 1 else intensities
    result = compute_selected_objective(objective_values, experimental, iref=mapped.iref,
                                        objective_mode=objective_mode)
    trial = RefinementTrialResult(result.objective, result.scaling, mapped, intensities,
                                   objective_values, amplitudes, objective_mode)
    if evaluation_callback is not None:
        evaluation_callback(trial)
    return trial
