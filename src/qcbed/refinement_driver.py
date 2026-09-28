"""Source-backed oneitpara.f stage orchestration; no optimizer algorithms."""

from dataclasses import dataclass
from numbers import Integral

import numpy as np

from .execution import ForwardExecutionConfig, _execution
from .anneal import run_legacy_anneal
from .amoeba import build_legacy_simplex, run_legacy_amoeba
from .paramin import run_legacy_paramin
from .simpara import run_legacy_simpara
from .multislice import _snapshot
from .refinement import (build_refinement_parameter_layout, pack_refinement_parameters,
                         apply_refinement_parameters, _physical_vector, _writeback_mapped_parameters)
from .refinement_categories import LEGACY_OPTIMIZER_PARAMETER_LIMIT
from .refinement_io import _parse_refinement_input
from .refinement_trial import evaluate_refinement_trial, _invalid_policy
from .structure_factors import _finite


_STAGES = {"ANNE": ("ANNEAL",), "SIMP": ("AMOEBA",), "PARA": ("PARAMIN",),
           "SPAR": ("SIMPARA",), "PAAN": ("PARAMIN","ANNEAL"),
           "ANPA": ("ANNEAL","PARAMIN"), "SIPA": ("AMOEBA","PARAMIN"),
           "PASI": ("PARAMIN","AMOEBA"), "SIAN": ("AMOEBA","ANNEAL"),
           "ANSI": ("ANNEAL","AMOEBA")}


@dataclass(frozen=True)
class LegacyRefinementConfig:
    """Legacy control names; seed is explicit Python reproducibility policy."""

    method: str
    loops: int
    maxfsimp: int
    maxfpara: int
    ftol: float
    chimin: float
    maxtemp: float
    mintemp: float
    coolrate: float
    plamda: tuple[float, ...]
    seed: int
    paraopt: int = 0
    loud: int = 0
    execution: ForwardExecutionConfig = ForwardExecutionConfig()

    def __post_init__(self):
        object.__setattr__(self, "execution", _execution(self.execution))
        for name in ("loops","maxfsimp","maxfpara","seed","paraopt","loud"):
            value = getattr(self,name)
            if isinstance(value,(bool,np.bool_)) or not isinstance(value,Integral):
                raise ValueError(f"{name} must be an integer excluding bool")
        if min(self.maxfsimp,self.maxfpara,self.seed) < 0:
            raise ValueError("optimizer limits and seed must be nonnegative")
        if self.paraopt not in (0,1) or self.loud not in (0,1):
            raise ValueError("PARAOPT and LOUD must be 0 or 1")
        for name in ("ftol","chimin","maxtemp","mintemp","coolrate"):
            object.__setattr__(self,name,_finite(getattr(self,name),name))
        object.__setattr__(self,"plamda",tuple(_finite(x,"PLAMDA") for x in self.plamda))
        if self.ftol < 0 or self.maxtemp <= 0 or self.mintemp <= 0 or not 0 < self.coolrate < 1:
            raise ValueError("require FTOL>=0, positive temperatures and 0<COOLRATE<1")
        if not isinstance(self.method,str):
            raise TypeError("METHOD must be text")
        self.stages  # Validate dispatch without accepting mixed case.

    @property
    def stages(self):
        prefix = self.method[:4]
        if prefix in _STAGES:
            return _STAGES[prefix]
        if prefix in {x.lower() for x in _STAGES}:
            return _STAGES[prefix.upper()]
        if self.loops <= 0:
            return ()
        raise ValueError(f"unsupported legacy METHOD {self.method!r}")


@dataclass(frozen=True, eq=False)
class LegacyRefinementStage:
    optimizer: str
    loop_index: int
    input_vector: np.ndarray
    input_objective: np.float64
    retained_vector: np.ndarray
    retained_objective: np.float64
    optimizer_result: object
    seed: int | None

    def __post_init__(self):
        for name in ("input_vector","retained_vector"):
            object.__setattr__(self,name,_snapshot(getattr(self,name)))


@dataclass(frozen=True, eq=False)
class LegacyRefinementResult:
    initial_vector: np.ndarray
    initial_objective: np.float64
    final_vector: np.ndarray
    final_objective: np.float64
    final_state: object
    loop_count: int
    stage_history: tuple[LegacyRefinementStage, ...]
    evaluations: int

    def __post_init__(self):
        for name in ("initial_vector","final_vector"):
            object.__setattr__(self,name,_snapshot(getattr(self,name)))
        object.__setattr__(self,"stage_history",tuple(self.stage_history))


@dataclass(frozen=True, eq=False)
class LegacyRefinementRestart:
    """Stage-boundary physical vector; does not contain optimizer internals.

    Keep the original baseline state/config: simplex scaling and seeds depend
    on it. loop_index is one-based, stage_position zero-based.
    """
    loop_index: int
    stage_position: int
    vector: np.ndarray

    def __post_init__(self):
        for name, minimum in (("loop_index", 1), ("stage_position", 0)):
            v = getattr(self, name)
            if isinstance(v, (bool, np.bool_)) or not isinstance(v, Integral) or v < minimum:
                raise ValueError(f"{name} must be an integer >= {minimum}")
        v = self.vector
        if not isinstance(v, np.ndarray) or v.dtype != np.float64 or v.ndim != 1 or not np.all(np.isfinite(v)):
            raise ValueError("restart vector must be finite one-dimensional float64")
        object.__setattr__(self, "vector", _snapshot(v))


def parse_legacy_refinement(path, *, seed, experimental_header_path=None):
    """Return (RefinementState, LegacyRefinementConfig) from source records."""
    state, controls, records, record = _parse_refinement_input(
        path, experimental_header_path=experimental_header_path)
    try:
        config = LegacyRefinementConfig(**controls,seed=seed)
    except (TypeError,ValueError) as exc:
        raise records.error("refinement controls",repr(controls),str(exc),record) from exc
    return state,config


def _run_legacy_refinement(state, experimental, config: LegacyRefinementConfig, *,
                          phase_grating_method="fft", objective_mode="legacy",
                          progress_callback=None, restart=None, stage_start_callback=None,
                          before_evaluation=None, invalid_trial_policy="raise") -> LegacyRefinementResult:
    """Run complete selected sequences, carrying PARMIN through UPDPARS/packing.

    Baseline state (and hence AMOEBA SCAL) stays anchored across all stages.
    No extra final objective evaluation or all-time best selection is added.
    """
    _invalid_policy(invalid_trial_policy)
    if not isinstance(config,LegacyRefinementConfig):
        raise TypeError("config must be LegacyRefinementConfig")
    if objective_mode not in ("legacy", "masked"):
        raise ValueError("objective_mode must be 'legacy' or 'masked'")
    layout = build_refinement_parameter_layout(state)
    n = len(layout.parameters)
    if n > LEGACY_OPTIMIZER_PARAMETER_LIMIT or (config.loops > 0 and n == 0):
        raise ValueError(f"driver requires 1..{LEGACY_OPTIMIZER_PARAMETER_LIMIT} active parameters for optimization")
    if len(config.plamda) != n:
        raise ValueError("PLAMDA count must match active parameter count")
    initial = pack_refinement_parameters(state)
    current = initial.copy()
    start_loop, start_position = 1, 0
    if restart is not None:
        if not isinstance(restart, LegacyRefinementRestart):
            raise TypeError("restart must be LegacyRefinementRestart")
        if config.loops <= 0 or restart.loop_index > config.loops or restart.stage_position >= len(config.stages):
            raise ValueError("restart cursor is outside configured stages")
        if restart.vector.shape != initial.shape:
            raise ValueError("restart vector must match active parameter layout")
        current = restart.vector.copy()
        start_loop, start_position = restart.loop_index, restart.stage_position
    history = []
    evaluations = 0
    first_score = None
    selected = config.stages if config.loops > 0 else (None,)
    for loop in range(start_loop,max(config.loops,1)+1):
        for position,name in enumerate(selected):
            if loop == start_loop and position < start_position:
                continue
            stage_index=(loop-1)*len(selected)+position+1
            def observe(trial):
                if progress_callback is not None:
                    progress_callback("evaluation",stage_index,name or "EVALUATION",loop,trial)
            pre = evaluate_refinement_trial(state,current,experimental,
                coordinate_mode="physical",phase_grating_method=phase_grating_method,
                objective_mode=objective_mode,execution=config.execution,before_evaluation=before_evaluation,evaluation_callback=observe if progress_callback else None)
            evaluations += 1
            score = np.float64(pre.objective)
            if not np.isfinite(score) or score > 1e20:
                raise ValueError("driver objective must be finite and <= source CHISMIN sentinel 1e20")
            # Preserve ONEITPARA's supported PT writebacks, without repacking.
            retained = current.copy()
            _writeback_mapped_parameters(retained, pre.physical_state)
            if first_score is None:
                first_score = score
            start = retained.copy()
            if stage_start_callback is not None:
                stage_start_callback(LegacyRefinementRestart(loop, position, start), pre)
            incoming_score = score
            stage_seed = None
            common = dict(invalid_trial_policy=invalid_trial_policy,phase_grating_method=phase_grating_method,objective_mode=objective_mode,
                          execution=config.execution,before_evaluation=before_evaluation,
                          evaluation_callback=observe if progress_callback else None)
            if name in ("ANNEAL","SIMPARA"):
                stage_seed = int(np.random.SeedSequence([config.seed,loop,position]).generate_state(1,dtype=np.uint64)[0])
            if name == "ANNEAL":
                result = run_legacy_anneal(state,start,experimental,
                    max_iterations=config.maxfsimp,seed=stage_seed,
                    maximum_temperature=config.maxtemp,minimum_temperature=config.mintemp,
                    cooling_rate=config.coolrate,**common)
            elif name == "PARAMIN":
                stage_seed = int(np.random.SeedSequence([config.seed,loop,position]).generate_state(1,dtype=np.uint64)[0]) if config.paraopt else None
                result = run_legacy_paramin(state,start,experimental,max_cycles=config.maxfpara,
                    tolerance=config.ftol,random_order=bool(config.paraopt),seed=stage_seed,**common)
            elif name == "SIMPARA":
                result = run_legacy_simpara(state,start,experimental,
                    uncertainties=np.asarray(config.plamda,dtype=np.float64),
                    max_evaluations=config.maxfsimp,high_temperature=config.maxtemp,
                    low_temperature=config.mintemp,cooling=config.coolrate,
                    seed=stage_seed,random_order=bool(config.paraopt),**common)
            elif name == "AMOEBA":
                simplex = build_legacy_simplex(state,start,np.asarray(config.plamda,dtype=np.float64))
                result = run_legacy_amoeba(state,simplex,experimental,ftol=config.ftol,
                    max_evaluations=config.maxfsimp,fmin=config.chimin,**common)
                # AMOEBA does not reset callback CHISMIN after the pre-stage call.
                if result.best_evaluated_objective <= score:
                    retained = _physical_vector(layout,result.best_evaluated_vector,"legacy_simplex")
                    score = result.best_evaluated_objective
                    # ONEIT restores simplex PT, constructs physical PARMIN,
                    # then UPDPARS(PARMIN) writes selected projected entries.
                    # Depth/voltage fold the model only; do not repack/floor P.
                    mapped = apply_refinement_parameters(state, retained)
                    _writeback_mapped_parameters(retained, mapped)
            if name is not None:
                if name != "AMOEBA":
                    retained,score = result.best_vector.copy(),result.best_objective
                evaluations += result.evaluations
                history.append(LegacyRefinementStage(name,loop,start,incoming_score,
                    retained,score,result,stage_seed))
                if progress_callback is not None:
                    progress_callback("stage",stage_index,name,loop,history[-1])
            final_state = apply_refinement_parameters(state,retained)
            # Source P and PARMIN differ: only P is packed/floored for next stage.
            current = pack_refinement_parameters(final_state)
    return LegacyRefinementResult(initial,first_score,retained,score,final_state,
                                  max(config.loops,0),tuple(history),evaluations)


def run_legacy_refinement(state, experimental, config: LegacyRefinementConfig, *,
                          phase_grating_method="fft", objective_mode="legacy",
                          progress_callback=None, restart=None, stage_start_callback=None,
                          before_evaluation=None, invalid_trial_policy="raise",
                          bounds_policy="none", refinement_bounds=None):
    """Run the unchanged driver with optional mapped-physical bound rejection."""
    from .refinement_bounds import refinement_bound_policy
    with refinement_bound_policy(state, bounds_policy=bounds_policy,
            refinement_bounds=refinement_bounds, invalid_trial_policy=invalid_trial_policy):
        return _run_legacy_refinement(state, experimental, config,
            phase_grating_method=phase_grating_method, objective_mode=objective_mode,
            progress_callback=progress_callback, restart=restart,
            stage_start_callback=stage_start_callback, before_evaluation=before_evaluation,
            invalid_trial_policy="reject" if bounds_policy == "reject" else invalid_trial_policy)
