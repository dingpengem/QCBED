"""Evaluation/report data assembly, independent of plotting and optimizer logic."""

from dataclasses import dataclass, field

import numpy as np

from .experimental import ExperimentalIntensityData
from .multislice import _snapshot
from .objective import _point_contributions, compute_legacy_objective
from .output_io import SimulatedCbedData
from .refinement import pack_refinement_parameters
from .refinement_driver import run_legacy_refinement
from .refinement_trial import evaluate_refinement_trial
from .display_mask import valid_display_mask
from .masked_objective import compute_selected_objective


@dataclass(frozen=True, eq=False)
class CbedComparison:
    experimental: ExperimentalIntensityData
    initial_simulated: SimulatedCbedData
    final_simulated: SimulatedCbedData
    residual: np.ndarray = field(init=False)
    numerator_by_beam: np.ndarray = field(init=False)
    denominator_by_beam: np.ndarray = field(init=False)
    rms_by_beam: np.ndarray = field(init=False)
    rms_all_by_beam: np.ndarray = field(init=False)
    rms_masked_by_beam: np.ndarray = field(init=False)
    valid_display: np.ndarray = field(init=False)

    def __post_init__(self):
        e=self.experimental
        if self.initial_simulated.objective_mode != self.final_simulated.objective_mode:
            raise ValueError("comparison objective modes must match")
        for s in (self.initial_simulated,self.final_simulated):
            if (s.beam_indices_hk!=e.beam_indices_hk or
                (s.scan.scan_rows,s.scan.scan_columns)!=(e.scan_rows,e.scan_columns) or
                s.gradient_option!=e.gradient_option):
                raise ValueError("comparison requires identical beam order, scan dimensions and GRADIENT option")
        residual=self.final_simulated.scaled_intensity-e.intensities
        valid = valid_display_mask(e)
        def rms(selection):
            return np.array([np.sqrt(np.mean(residual[selection[:, b], b]**2))
                             if np.any(selection[:, b]) else np.nan
                             for b in range(residual.shape[1])], dtype=np.float64)
        if self.final_simulated.objective_mode == "legacy":
            errors,norms=_point_contributions(self.final_simulated.scaled_intensity,e)
        else:
            if np.any(e.variances[valid] <= 0):
                raise ValueError("masked contributions require positive variance at nonzero-weight points")
            errors,norms=np.zeros_like(residual),np.zeros_like(residual)
            w,v,x=e.weights[valid],e.variances[valid],e.intensities[valid]
            errors[valid]=w*residual[valid]*residual[valid]/v
            norms[valid]=w*x*x/v
        # Retain source tilt accumulation order within each beam.
        numerator=np.zeros(len(e.beam_indices_hk),dtype=np.float64)
        denominator=numerator.copy()
        for b in range(len(numerator)):
            for t in range(e.scan_rows*e.scan_columns):
                numerator[b]+=errors[t,b]
                denominator[b]+=norms[t,b]
        for name,value in (("residual",residual),("numerator_by_beam",numerator),
                           ("denominator_by_beam",denominator),
                           ("rms_by_beam",rms(valid)),
                           ("rms_masked_by_beam",rms(~valid)),
                           ("rms_all_by_beam",np.sqrt(np.mean(residual**2,axis=0))),
                           ("valid_display",valid)):
            object.__setattr__(self,name,_snapshot(value))


@dataclass(frozen=True)
class RefinementValidation:
    refinement: object
    comparison: CbedComparison


@dataclass(frozen=True)
class _EvaluationOnlyResult:
    """Diagnostic snapshot, not a successful legacy-driver execution."""
    initial_vector: np.ndarray
    initial_objective: float
    final_vector: np.ndarray
    final_objective: float
    final_state: object
    loop_count: int = 0
    stage_history: tuple = ()
    evaluations: int = 0


def _simulated(trial,experimental):
    return SimulatedCbedData(trial.physical_state.scan,experimental.beam_indices_hk,
        trial.calculated_intensities,trial.objective_intensities,trial.scaling,
        experimental.gradient_option,trial.objective,trial.objective_mode)


def validate_refinement(state, experimental, config, *, phase_grating_method="fft", objective_mode="legacy"):
    """Run/refine and independently evaluate states, or evaluate only at LOOPS<=0.

    Reevaluations are validation work, not driver history. If vectors are exactly
    equal the initial trial is reused for final maps, explicitly including the
    evaluation-only path. Baseline scan dimensions and physics are not reduced.
    LOOPS<=0 calls the trial evaluator directly: undefined objectives are retained
    for diagnostics, not fed into or accepted by the guarded optimizer driver.
    """
    # Evaluation-only diagnostics deliberately do not enter the optimizer or
    # bypass its sentinel guard. Undefined objective values stay undefined.
    if config.loops <= 0:
        vector=_snapshot(pack_refinement_parameters(state))
        trial=evaluate_refinement_trial(state,vector,experimental,
                                        phase_grating_method=phase_grating_method,objective_mode=objective_mode)
        check=compute_selected_objective(trial.objective_intensities,experimental,iref=state.iref,objective_mode=objective_mode)
        if not (np.isclose(trial.objective,check.objective,equal_nan=True) and
                np.isclose(trial.scaling,check.scaling,equal_nan=True)):
            raise ValueError("evaluation-only objective/YNORM inconsistent")
        simulated=_simulated(trial,experimental)
        result=_EvaluationOnlyResult(vector,trial.objective,vector,trial.objective,trial.physical_state)
        return RefinementValidation(result,CbedComparison(experimental,simulated,simulated))
    result=run_legacy_refinement(state,experimental,config,phase_grating_method=phase_grating_method,objective_mode=objective_mode)
    initial=evaluate_refinement_trial(state,pack_refinement_parameters(state),experimental,
                                      phase_grating_method=phase_grating_method,objective_mode=objective_mode)
    final=(initial if np.array_equal(result.initial_vector,result.final_vector) else
           evaluate_refinement_trial(state,result.final_vector,experimental,
                                      phase_grating_method=phase_grating_method,objective_mode=objective_mode))
    for label,trial,score in (("initial",initial,result.initial_objective),
                             ("final",final,result.final_objective)):
        check=compute_selected_objective(trial.objective_intensities,experimental,iref=state.iref,objective_mode=objective_mode)
        if not (np.isclose(trial.objective,score,rtol=1e-12,atol=1e-14) and
                trial.objective==check.objective and trial.scaling==check.scaling):
            raise ValueError(f"{label} objective/YNORM inconsistent with independent validation")
    return RefinementValidation(result,CbedComparison(experimental,_simulated(initial,experimental),
                                                      _simulated(final,experimental)))
