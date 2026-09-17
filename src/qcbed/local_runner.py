"""Local orchestration and observational history, using the existing driver."""

from dataclasses import dataclass
from numbers import Integral
from pathlib import Path
import json
import time

import numpy as np

from .execution import ForwardExecutionConfig
from .multislice import _snapshot
from .refinement import _entries, apply_refinement_parameters, build_refinement_parameter_layout
from .refinement_driver import run_legacy_refinement
from .refinement_trial import evaluate_refinement_trial
from .validation import CbedComparison, RefinementValidation, _simulated


@dataclass(frozen=True, eq=False)
class RefinementProgressRecord:
    evaluation_index: int
    stage_index: int
    stage_name: str
    loop_index: int
    kind: str  # evaluated trial or stage-retained state, never conflated
    objective: float
    ynorm: float
    parameter_vector: np.ndarray  # physical, no simplex coordinates

    def __post_init__(self):
        object.__setattr__(self,"parameter_vector",_snapshot(self.parameter_vector))


@dataclass(frozen=True)
class RefinementSnapshot:
    label: str
    progress: RefinementProgressRecord
    simulated: object


@dataclass(frozen=True)
class LocalOptimizationRunResult:
    validation: RefinementValidation
    progress: tuple[RefinementProgressRecord, ...]
    snapshots: tuple[RefinementSnapshot, ...]
    objective_mode: str
    elapsed_seconds: float
    execution: ForwardExecutionConfig = ForwardExecutionConfig()

    @property
    def final_physical_parameters(self):
        return self.validation.refinement.final_state


def _physical_values(state):
    return np.asarray([v for _,_,v in _entries(state)],dtype=np.float64)


def run_local_optimization(state, experimental, config, *, objective_mode="legacy",
                           snapshot_every=None, record_progress=True, output_directory=None,
                           reference_input_name="", phase_grating_method="fft"):
    """Run the existing driver; record trials and retained stage endpoints.

    Cadence is global 1-based *driver objective evaluations*, including pre-stage
    and initialization calls, not optimizer iterations. None disables intermediate
    snapshots. Initial/final are mandatory. A final diagnostic reevaluation occurs
    after the driver, regardless of recording settings, and is not an optimizer call.
    Output directories must be new to avoid overwriting prior runs or input files.
    """
    if objective_mode not in ("legacy","masked"):
        raise ValueError("objective_mode must be 'legacy' or 'masked'")
    if snapshot_every is not None and (isinstance(snapshot_every,bool) or
            not isinstance(snapshot_every,Integral) or snapshot_every < 1):
        raise ValueError("snapshot_every must be None or a positive integer")
    if not isinstance(record_progress,bool):
        raise TypeError("record_progress must be bool")
    destination=Path(output_directory).resolve() if output_directory is not None else None
    if destination is not None:
        root=Path(__file__).resolve().parents[2]
        if any(destination.is_relative_to(root/p) for p in ("legacy_fortran","reference_data")):
            raise ValueError("output_directory must not be inside protected source/reference directories")
        destination.mkdir(parents=True,exist_ok=False)
    started=time.perf_counter()
    baseline=repr(state)
    history=[]
    snapshots=[]
    stage_values=[]  # small scalar/vector metadata only, cleared at each stage
    count=0
    initial=None

    def record(record):
        if record_progress:
            history.append(record)
            if destination is not None:
                row={k:(v.tolist() if isinstance(v,np.ndarray) else v)
                     for k,v in vars(record).items()}
                with (destination/"progress.jsonl").open("a",encoding="utf-8") as stream:
                    stream.write(json.dumps(row,allow_nan=False)+"\n")

    def snapshot(label,record,simulated):
        snap=RefinementSnapshot(label,record,simulated)
        snapshots.append(snap)
        if destination is not None:
            np.savez_compressed(destination/f"{label}.npz",raw=simulated.raw_intensity,
                objective_intensity=simulated.objective_intensity,scaled=simulated.scaled_intensity,
                tilts=simulated.tilts,beams=np.asarray(simulated.beam_indices_hk),
                parameters=record.parameter_vector,objective=record.objective,ynorm=record.ynorm,
                objective_mode=objective_mode)
            # Immediate per-beam PNGs allow inspection while optimization runs.
            from .visualization import plot_cbed_comparison
            c=CbedComparison(experimental,initial,simulated)
            for b in range(len(simulated.beam_indices_hk)):
                fig=plot_cbed_comparison(c,beam_index=b)
                if label != "final":
                    fig.axes[2].set_title("Initial simulation" if label=="initial" else "Checkpoint trial")
                    fig.axes[3].set_title("Checkpoint − experimental")
                fig.savefig(destination/f"{label}_beam_{b+1:02d}.png",dpi=100)
                fig.clear()

    def observe(kind,stage_index,name,loop,payload):
        nonlocal count,initial
        if kind=="evaluation":
            count+=1
            if not np.isfinite(payload.objective) or payload.objective > 1e20:
                return  # let the existing driver/optimizer raise its domain error
            vector=_physical_values(payload.physical_state)
            r=RefinementProgressRecord(count,stage_index,name,loop,"evaluation",
                                       float(payload.objective),float(payload.scaling),vector)
            stage_values.append(r)
            record(r)
            if initial is None:
                initial=_simulated(payload,experimental)
                snapshot("initial",r,initial)
            elif snapshot_every is not None and count % snapshot_every == 0:
                snapshot(f"evaluation_{count:06d}",r,_simulated(payload,experimental))
        else:
            physical=_physical_values(apply_refinement_parameters(state,payload.retained_vector))
            matches=[r for r in stage_values if r.objective==payload.retained_objective and
                     np.array_equal(r.parameter_vector,physical)]
            if not matches:
                raise ValueError("retained stage state has no matching evaluated YNORM")
            r=RefinementProgressRecord(count,stage_index,name,loop,"retained",
                float(payload.retained_objective),matches[-1].ynorm,payload.retained_vector)
            record(r)
            stage_values.clear()

    result=run_legacy_refinement(state,experimental,config,objective_mode=objective_mode,
        phase_grating_method=phase_grating_method,progress_callback=observe)
    final_trial=evaluate_refinement_trial(state,result.final_vector,experimental,
        objective_mode=objective_mode,phase_grating_method=phase_grating_method,execution=config.execution)
    if not np.isclose(final_trial.objective,result.final_objective,rtol=1e-12,atol=1e-14):
        raise ValueError("final diagnostic objective disagrees with retained driver objective")
    final=_simulated(final_trial,experimental)
    if repr(state) != baseline:
        raise RuntimeError("baseline problem changed during optimization")
    stage=result.stage_history[-1] if result.stage_history else None
    final_record=RefinementProgressRecord(count,len(result.stage_history),
        stage.optimizer if stage else "EVALUATION",stage.loop_index if stage else 0,
        "final",float(final.objective),float(final.ynorm),result.final_vector)
    snapshot("final",final_record,final)
    validation=RefinementValidation(result,CbedComparison(experimental,initial,final))
    run=LocalOptimizationRunResult(validation,tuple(history),tuple(snapshots),objective_mode,
                                   time.perf_counter()-started,config.execution)
    if destination is not None:
        from .validation_report import generate_refinement_validation_report
        generate_refinement_validation_report(destination/"report.html",validation,
            layout=build_refinement_parameter_layout(state),config=config,
            reference_input_name=reference_input_name,local_run=run)
        from dataclasses import asdict
        summary=dict(objective_mode=objective_mode,execution_backend=config.execution.execution_backend,
            batch_size=config.execution.batch_size,process_workers=None,seed=config.seed,stages=list(config.stages),
            loops=config.loops,maxfsimp=config.maxfsimp,maxfpara=config.maxfpara,
            initial_objective=float(initial.objective),final_objective=float(final.objective),
            improvement=float(initial.objective-final.objective),initial_ynorm=float(initial.ynorm),
            final_ynorm=float(final.ynorm),initial_parameters=result.initial_vector.tolist(),
            final_parameters=result.final_vector.tolist(),
            parameter_changes=(result.final_vector-result.initial_vector).tolist(),
            driver_evaluations=result.evaluations,baseline_unchanged=True,
            elapsed_seconds=time.perf_counter()-started)
        (destination/"summary.json").write_text(json.dumps(summary,indent=2,allow_nan=False),encoding="utf-8")
        (destination/"final_physical_parameters.json").write_text(
            json.dumps(asdict(result.final_state),indent=2,default=lambda v:v.item() if isinstance(v,np.generic) else str(v)),
            encoding="utf-8")
    return run
