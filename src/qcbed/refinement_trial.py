"""Supported source-backed trial assembly; no optimizer or mutable best state."""

from dataclasses import dataclass, replace

import numpy as np

from .absorption import AbsorptionCoefficients, absorptive_structure_factors
from .execution import _execution
from .tilt_series import simulate_tilt_series_batched
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
        from .depth_geometry import requested_depth_to_sample_geometry
        extent = requested_depth_to_sample_geometry(r.depth_angstrom,
            (g.thickness_angstrom for g in state.model.slices[offset:offset+r.slice_count]))
        sequence.extend(offset+i for i in extent.sequence)
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


def _evaluate_refinement_trial(
    state: RefinementState, trial_vector: np.ndarray, experimental: ExperimentalIntensityData, *,
    coordinate_mode="physical", phase_grating_method="fft", objective_mode="legacy",
    evaluation_callback=None, execution=None, before_evaluation=None,
) -> RefinementTrialResult:
    """Map a fresh trial, rebuild physics, run selected tilts, compute YNORM/objective.

    Supports S atoms, zero vacancy and no surface correction; rejects unresolved
    physical updates explicitly. All source reflection symmetries are supported.
    The objective API computes/applies YNORM exactly once. No best-state updates.
    execution=None retains the scalar reference path; ForwardExecutionConfig
    selects the existing batched NumPy engine without changing objective policy.
    """
    execution = _execution(execution)
    if not isinstance(experimental, ExperimentalIntensityData):
        raise TypeError("experimental must be ExperimentalIntensityData")
    if objective_mode not in ("legacy", "masked"):
        raise ValueError("objective_mode must be 'legacy' or 'masked'")
    if before_evaluation is not None:
        before_evaluation()
    mapped = apply_refinement_parameters(state, trial_vector, coordinate_mode=coordinate_mode)
    from .refinement_bounds import _scope
    policy, bounds, _ = _scope.get()
    if policy == "reject":
        violations = bounds.violations(mapped)
        if violations:
            return BoundRejectedTrial(trial_vector, coordinate_mode, mapped,
                "", "bound_rejection", violations=violations)
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
    if execution.execution_backend == "batched":
        # Preserve the trial-specific SFMAP2 reflection multipliers. Reuse the
        # prepared geometry/parameters and replace only the final gratings.
        prepared = replace(base, slice_types=tuple(
            replace(stored, phase_grating=grating)
            for stored, grating in zip(base.slice_types, gratings)))
        series = simulate_tilt_series_batched(prepared, tilts, batch_size=execution.batch_size)
        amplitudes, intensities = series.amplitudes, series.intensities
    else:
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


@dataclass(frozen=True, eq=False)
class InvalidRefinementTrial:
    """Rejected evaluation metadata only: no fabricated amplitudes or YNORM."""
    trial_vector: np.ndarray
    coordinate_mode: str
    physical_state: RefinementState | None
    exception_type: str
    reason: str
    objective: float = float("inf")

    def __post_init__(self):
        object.__setattr__(self, "trial_vector", _snapshot(self.trial_vector))


def _invalid_policy(value):
    if value not in ("raise", "reject"):
        raise ValueError("invalid_trial_policy must be 'raise' or 'reject'")
    return value


@dataclass(frozen=True, eq=False)
class BoundRejectedTrial(InvalidRefinementTrial):
    """Constraint rejection, not a numerical exception; shared rejecting protocol."""
    violations: tuple = ()

    def __post_init__(self):
        super().__post_init__()
        object.__setattr__(self, "violations", tuple(self.violations))


def _numerical_trial_error(exc):
    # Do not hide geometry/configuration errors, unsupported physics, invalid
    # experimental variances, user callback errors, or unexpected exceptions.
    if isinstance(exc, (FloatingPointError, OverflowError)):
        return True
    message = str(exc)
    if message.startswith("masked objective: nonzero-weight point") or "experimental normalizer" in message:
        return False
    return isinstance(exc, ValueError) and any(text in message for text in (
        "must be finite", "must remain finite", "must contain finite", "nonfinite", "non-finite"))


def _evaluate_with_invalid_policy(
    state: RefinementState, trial_vector: np.ndarray, experimental: ExperimentalIntensityData, *,
    coordinate_mode="physical", phase_grating_method="fft", objective_mode="legacy",
    evaluation_callback=None, execution=None, before_evaluation=None,
    invalid_trial_policy="raise",
) -> RefinementTrialResult | InvalidRefinementTrial:
    """Evaluate unchanged physics, optionally rejecting numerical trial failures.

    reject is a Python production-run numerical safety policy; not proven
    historical Fortran behavior. raise preserves original errors and nonfinite
    low-level objective results. Callbacks/budget checks are outside the catch.
    """
    _invalid_policy(invalid_trial_policy)
    if invalid_trial_policy == "reject":
        from .refinement import build_refinement_parameter_layout
        count = len(build_refinement_parameter_layout(state).parameters)
        if (not isinstance(trial_vector, np.ndarray) or trial_vector.dtype != np.float64 or
                trial_vector.shape != (count,)):
            # Preserve the original structural/type error, outside the numerical
            # failure catch. Only a correctly shaped proposal can be rejected.
            apply_refinement_parameters(state, trial_vector, coordinate_mode=coordinate_mode)
    if before_evaluation is not None:
        before_evaluation()
    try:
        # Preserve the original arithmetic/errstate in raise mode. Reject mode
        # turns overflow into an exception at its first arithmetic occurrence;
        # no finite validation is disabled and no array is repaired.
        from contextlib import nullcontext
        context = np.errstate(over="raise", invalid="raise", divide="raise") if invalid_trial_policy == "reject" else nullcontext()
        with context:
            result = _evaluate_refinement_trial(state, trial_vector, experimental,
                coordinate_mode=coordinate_mode, phase_grating_method=phase_grating_method,
                objective_mode=objective_mode, execution=execution)
        if not isinstance(result, BoundRejectedTrial) and invalid_trial_policy == "reject" and (not np.isfinite(result.objective) or not np.isfinite(result.scaling)):
            raise FloatingPointError("YNORM/objective must be finite")
    except (ValueError, FloatingPointError, OverflowError) as exc:
        if invalid_trial_policy != "reject" or not _numerical_trial_error(exc):
            raise
        try:
            mapped = apply_refinement_parameters(state, trial_vector, coordinate_mode=coordinate_mode)
        except (ValueError, FloatingPointError, OverflowError):
            mapped = None
        result = InvalidRefinementTrial(trial_vector, coordinate_mode, mapped,
                                        type(exc).__name__, str(exc))
    if evaluation_callback is not None:
        evaluation_callback(result)
    return result


def evaluate_refinement_trial(
    state: RefinementState, trial_vector: np.ndarray, experimental: ExperimentalIntensityData, *,
    coordinate_mode="physical", phase_grating_method="fft", objective_mode="legacy",
    evaluation_callback=None, execution=None, before_evaluation=None,
    invalid_trial_policy="raise", bounds_policy=None, refinement_bounds=None,
) -> RefinementTrialResult | InvalidRefinementTrial:
    """Evaluate a mapped physical trial, with independently selected safety policies.

    bounds_policy=None inherits an explicit driver scope, otherwise means 'none'.
    Explicit 'none' disables bounds for this call. 'reject' checks inclusive
    physical limits before preparing the sample or running forward simulation.
    """
    from .refinement_bounds import _scope, refinement_bound_policy
    inherited_policy, inherited_bounds, numerical = _scope.get()
    policy = inherited_policy if bounds_policy is None else bounds_policy
    bounds = inherited_bounds if refinement_bounds is None else refinement_bounds
    with refinement_bound_policy(state, bounds_policy=policy, refinement_bounds=bounds):
        return _evaluate_with_invalid_policy(state, trial_vector, experimental,
            coordinate_mode=coordinate_mode, phase_grating_method=phase_grating_method,
            objective_mode=objective_mode, evaluation_callback=evaluation_callback,
            execution=execution, before_evaluation=before_evaluation,
            invalid_trial_policy=invalid_trial_policy if numerical is None else numerical)
