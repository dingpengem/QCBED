"""Durable orchestration of the existing driver; no new search algorithm.

Only driver-confirmed stage boundaries are restart points. Trial callbacks do
not expose optimizer acceptance, RNG state or in-progress retained state.
"""

from dataclasses import asdict, dataclass
from html import escape
import hashlib
import json
from numbers import Integral
import os
from pathlib import Path
from time import perf_counter
from types import SimpleNamespace

import numpy as np

from .refinement import (apply_refinement_parameters, build_refinement_parameter_layout,
                         pack_refinement_parameters, _entries)
from .refinement_driver import (LegacyRefinementConfig, LegacyRefinementRestart,
                                LegacyRefinementResult, LegacyRefinementStage,
                                run_legacy_refinement)
from .refinement_trial import evaluate_refinement_trial, InvalidRefinementTrial, BoundRejectedTrial, _invalid_policy
from .refinement_bounds import RefinementBounds, RefinementBoundRules, resolve_refinement_bounds, _validate_policy
from .validation import CbedComparison, RefinementValidation, _simulated
from .output_io import SimulatedCbedData, write_legacy_out_dat


POLICY = ("Masked objective is an explicit non-legacy numerical policy in which zero-weight "
          "points are excluded before variance division and before YNORM computation.")
RESTART_POLICY = ("Only initial/pre-stage/completed-stage retained states are durable. "
    "Interrupted stages restart from their last driver-confirmed boundary with the original "
    "stage seed. In-stage accepted improvements and optimizer/RNG internals are unavailable; "
    "this is not exact mid-optimizer continuation.")


@dataclass(frozen=True)
class FullRefinementConfig:
    refinement: LegacyRefinementConfig
    objective_mode: str = "masked"
    invalid_trial_policy: str = "reject"
    refinement_categories: tuple[str, ...] | None = None
    snapshot_every: int | None = 100
    checkpoint_every: int = 25
    max_total_evaluations: int | None = 1000
    max_wall_time: float | None = 3600.
    bounds_policy: str = "none"
    refinement_bounds: RefinementBounds = RefinementBounds()
    bound_rules: RefinementBoundRules | None = None

    def __post_init__(self):
        _invalid_policy(self.invalid_trial_policy)
        if self.bounds_policy not in ("none", "reject"):
            raise ValueError("bounds_policy must be 'none' or 'reject'")
        if not isinstance(self.refinement_bounds, RefinementBounds):
            raise TypeError("refinement_bounds must be RefinementBounds")
        if self.bound_rules is not None and not isinstance(self.bound_rules, RefinementBoundRules):
            raise TypeError("bound_rules must be RefinementBoundRules or None")
        if self.refinement_categories is not None:
            object.__setattr__(self, "refinement_categories", tuple(self.refinement_categories))
        if not isinstance(self.refinement, LegacyRefinementConfig):
            raise TypeError("refinement must be LegacyRefinementConfig")
        if self.objective_mode not in ("legacy", "masked"):
            raise ValueError("objective_mode must be 'legacy' or 'masked'")
        for name, minimum in (("snapshot_every", 1), ("checkpoint_every", 1), ("max_total_evaluations", 2)):
            value = getattr(self, name)
            if value is None and name != "checkpoint_every":
                continue
            if isinstance(value, (bool, np.bool_)) or not isinstance(value, Integral) or value < minimum:
                raise ValueError(f"{name} must be an integer >= {minimum}")
        if self.max_wall_time is not None and (isinstance(self.max_wall_time, (bool, np.bool_)) or
                not np.isfinite(self.max_wall_time) or self.max_wall_time <= 0):
            raise ValueError("max_wall_time must be finite and positive or None")


@dataclass(frozen=True)
class FullRefinementRun:
    output_directory: Path
    status: str
    summary: dict


class _BudgetStop(Exception):
    pass


def _json(value):
    return json.dumps(value, indent=2, allow_nan=False,
        default=lambda v: v.tolist() if isinstance(v, np.ndarray) else v.item() if isinstance(v, np.generic) else str(v))


def _atomic_json(path, value):
    temporary = path.with_name(path.name + ".tmp")
    with temporary.open("w", encoding="utf-8") as stream:
        stream.write(_json(value))
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, path)


def _save_arrays(path, **arrays):
    temporary = path.with_name(path.name + ".tmp")
    with temporary.open("wb") as stream:
        np.savez_compressed(stream, **arrays)
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, path)


def _digest(state, experimental):
    h = hashlib.sha256(_json(asdict(state)).encode())
    h.update(_json((experimental.scan_rows, experimental.scan_columns,
                   experimental.gradient_option, experimental.beam_indices_hk)).encode())
    for a in (experimental.intensities, experimental.weights, experimental.variances):
        h.update(a.tobytes())
    return h.hexdigest()


def load_refinement_checkpoint(path):
    """Load versioned JSON, never pickle; input/config identity is checked by runner."""
    path = Path(path)
    if path.is_dir():
        path = path / "checkpoint_latest.json"
    value = json.loads(path.read_text(encoding="utf-8"))
    if value.get("format") != "qcbed-retained-stage" or value.get("version") != 1:
        raise ValueError("unsupported refinement checkpoint format/version")
    vector = np.asarray(value["retained_vector"], dtype=np.float64)
    if vector.ndim != 1 or not np.all(np.isfinite(vector)):
        raise ValueError("checkpoint retained vector must be finite and one-dimensional")
    for name in ("evaluation_count", "completed_evaluations", "driver_evaluations"):
        if type(value[name]) is not int or value[name] < 0:
            raise ValueError(f"checkpoint {name} must be a nonnegative integer")
    for name in ("retained_objective", "retained_ynorm"):
        if value[name] is not None and not np.isfinite(value[name]):
            raise ValueError(f"checkpoint {name} must be finite")
    return value


def _preview(path, experimental, initial, simulated, metadata):
    from .visualization import plot_cbed_comparison
    path.mkdir(exist_ok=True)
    comparison = CbedComparison(experimental, initial, simulated)
    rows = ["<!doctype html><meta charset='utf-8'><title>QCBED retained preview</title>",
        "<style>body{font-family:sans-serif;max-width:1500px;margin:auto}img{width:100%}</style>",
        "<h1>Retained-state CBED comparison</h1><pre>" + escape(_json(metadata)) + "</pre>",
        "<p>Shared weight mask and experimental/simulated scale; residual centered on zero. "
        "No independent display normalization. Initial and current retained simulations "
        "occupy the existing four-panel comparison layout. No scientific convergence claim.</p>"]
    for b, beam in enumerate(experimental.beam_indices_hk):
        fig = plot_cbed_comparison(comparison, beam_index=b, mask_zero_weight=True)
        fig.axes[2].set_title("Current retained simulation")
        fig.axes[3].set_title("Retained minus experimental")
        name = f"beam_{b+1:02d}.png"
        fig.savefig(path / name, dpi=90)
        fig.clear()
        rows.append(f"<h2>Beam H={beam[0]}, K={beam[1]}</h2><img loading='lazy' src='{name}'>")
    (path / "index.html").write_text("\n".join(rows), encoding="utf-8")


def run_full_refinement(state, experimental, config: FullRefinementConfig, *,
                        output_directory, input_paths=(), resume=False, progress_hook=None,
                        checkpoint_identity=None, final_output_path=None, results_metadata=None):
    """Run with durable boundary checkpoints and optional safe stage restart.

    Limits count attempted evaluations, including initial, snapshot diagnostics
    and final verification. One slot is reserved for final verification. Wall
    limits are cooperative and may be exceeded by an in-flight call/finalization.
    progress_hook observes completed log records; it must not mutate the problem.
    """
    if not isinstance(config, FullRefinementConfig):
        raise TypeError("config must be FullRefinementConfig")
    destination = Path(output_directory).resolve()
    root = Path(__file__).resolve().parents[2]
    if any(destination.is_relative_to(root / p) for p in ("legacy_fortran", "reference_data", "tests/fixtures")):
        raise ValueError("output_directory is inside protected source/reference directories")
    export_path = None if final_output_path is None else Path(final_output_path).resolve()
    if export_path is not None:
        if any(export_path.is_relative_to(root / p) for p in ("legacy_fortran", "reference_data", "tests/fixtures")):
            raise ValueError("final_output_path is inside protected source/reference directories")
        if not resume and export_path.exists():
            raise FileExistsError(export_path)
    # Optional frontend identity is JSON-normalized; absent preserves legacy checkpoints.
    checkpoint_identity = json.loads(_json(checkpoint_identity))
    files = {str(Path(p).resolve()): hashlib.sha256(Path(p).read_bytes()).hexdigest() for p in input_paths}
    cfg = config.refinement
    layout = build_refinement_parameter_layout(state)
    _validate_policy(config.bounds_policy, config.refinement_bounds, state)
    if config.bound_rules is not None and resolve_refinement_bounds(state, config.bound_rules) != config.refinement_bounds:
        raise ValueError("resolved bounds must exactly match bound_rules and mapped starting layout")
    identity = _digest(state, experimental)
    settings = json.loads(_json(asdict(config)))
    labels = [f"{p.vector_index+1}: {p.category} {p.target}" for p in layout.parameters]
    if len(cfg.plamda) != len(labels) or (cfg.loops > 0 and not labels):
        raise ValueError("optimizer configuration does not match active parameter layout")
    if resume:
        cp = load_refinement_checkpoint(destination)
        if cp.get("checkpoint_identity") != checkpoint_identity:
            raise ValueError("checkpoint frontend/configuration identity mismatch")
        if cp["problem_digest"] != identity or cp["input_files"] != files:
            raise ValueError("checkpoint input/problem fingerprint mismatch")
        if cp["config"]["refinement"] != settings["refinement"] or cp["objective_mode"] != config.objective_mode:
            raise ValueError("resume requires the original optimizer/execution configuration and objective mode")
        if cp["layout"] != json.loads(_json(asdict(layout))):
            raise ValueError("checkpoint active parameter layout mismatch")
        for name, default in (("bounds_policy", "none"), ("refinement_bounds", {"entries": []}), ("bound_rules", None)):
            if cp["config"].get(name, default) != settings[name]:
                raise ValueError("resume requires the original hard-bound policy and bounds")
        # A hard exit can leave completed trial log rows newer than the latest
        # periodic checkpoint. Keep its safe boundary, but never reuse indices.
        progress_path = destination / "progress.jsonl"
        completed_indices, driver_count = set(), 0
        invalid_count, invalid_stages, invalid_reasons = 0, {}, {}
        bound_count, bound_stages = 0, {}
        if progress_path.exists():
            with progress_path.open(encoding="utf-8") as stream:
                for line_number, line in enumerate(stream, 1):
                    try:
                        row = json.loads(line)
                    except json.JSONDecodeError as exc:
                        raise ValueError(f"progress.jsonl line {line_number} is incomplete/corrupt; preserve it for explicit recovery") from exc
                    cp["evaluation_count"] = max(cp["evaluation_count"], row["global_evaluation_index"])
                    cp["elapsed_seconds"] = max(cp["elapsed_seconds"], row["elapsed_wall_time"] if "elapsed_wall_time" in row else row["elapsed_seconds"])
                    if row["classification"] in ("trial", "initial_retained", "verification", "snapshot_verification"):
                        completed_indices.add(row["global_evaluation_index"])
                    driver_count += row["classification"] in ("trial", "invalid_trial", "bound_rejection")
                    if row["classification"] == "bound_rejection":
                        bound_count += 1
                        bound_stages[row["stage"]] = bound_stages.get(row["stage"], 0)+1
                    if row["classification"] == "invalid_trial":
                        invalid_count += 1
                        invalid_stages[row["stage"]] = invalid_stages.get(row["stage"], 0)+1
                        reason = row["exception_type"]+": "+row["reason"]
                        invalid_reasons[reason] = invalid_reasons.get(reason, 0)+1
        cp.update(invalid_trial_count=invalid_count, invalid_trials_by_stage=invalid_stages, invalid_trial_reasons=invalid_reasons)
        cp.update(bound_rejection_count=bound_count, bound_rejections_by_stage=bound_stages)
        cp["completed_evaluations"] = max(cp["completed_evaluations"], len(completed_indices))
        cp["driver_evaluations"] = max(cp["driver_evaluations"], driver_count)
    else:
        destination.mkdir(parents=True, exist_ok=False)
        cp = dict(format="qcbed-retained-stage", version=1, problem_digest=identity,
            input_files=files, layout=asdict(layout), parameter_labels=labels,
            retained_vector=pack_refinement_parameters(state).tolist(), retained_objective=None,
            retained_ynorm=None, retained_at_evaluation=0, evaluation_count=0,
            completed_evaluations=0, driver_evaluations=0, elapsed_seconds=0.,
            loop_index=1, stage_position=0, stage_name="INITIAL", next_vector=pack_refinement_parameters(state).tolist(),
            config=settings, objective_mode=config.objective_mode, seed=cfg.seed,
            execution_backend=cfg.execution.execution_backend, batch_size=cfg.execution.batch_size,
            rng_state=None, restart_policy=RESTART_POLICY, stages=[], completed=False, restart_count=0)
        if checkpoint_identity is not None:
            cp["checkpoint_identity"] = checkpoint_identity
        _atomic_json(destination / "run_config.json", dict(**settings, input_files=files,
            parameter_labels=labels, layout=asdict(layout), masked_policy=POLICY, restart_policy=RESTART_POLICY,
            results_metadata=results_metadata))
    # Refuse concurrent writers. A hard process kill may leave this lock behind;
    # manual removal is allowed only after verifying the old process has exited.
    lock = destination / ".run.lock"
    with lock.open("x", encoding="ascii") as stream:
        stream.write(str(os.getpid()))
    (destination / "progress.jsonl").touch(exist_ok=True)
    started = perf_counter()
    previous_elapsed = cp["elapsed_seconds"]
    cp["config"] = settings
    cp["restart_count"] += int(resume)
    cp.setdefault("invalid_trial_count", 0)
    cp.setdefault("invalid_trials_by_stage", {})
    cp.setdefault("invalid_trial_reasons", {})
    cp.setdefault("bound_rejection_count", 0)
    cp.setdefault("bound_rejections_by_stage", {})
    stage_values = {}
    local_index = 0
    observed_stage = None
    current_stage = cp["stage_name"]
    pending_snapshot = False
    initial = None
    status, reason = "completed", "configured stages completed; convergence not inferred"

    def elapsed():
        return previous_elapsed + perf_counter() - started

    def checkpoint(label="running"):
        cp["status"] = label
        cp["elapsed_seconds"] = elapsed()
        # Publish immutable numbered file, then latest pointer content.
        name = f"checkpoint_{cp['evaluation_count']:06d}_{len(cp['stages']):04d}_{cp['restart_count']:03d}_{label}.json"
        _atomic_json(destination / name, cp)
        _atomic_json(destination / "checkpoint_latest.json", cp)

    latest_display = None

    def record_completed_evaluation(trial):
        # Output only: copy completed trial values, never request another trial.
        nonlocal latest_display
        sim = _simulated(trial, experimental)
        latest_display = dict(evaluation_index=cp["evaluation_count"],
            physical_active_vector=[v for _, _, v in _entries(trial.physical_state)])
        _save_arrays(destination / "latest_evaluation.npz",
            evaluation_index=cp["evaluation_count"], vector=latest_display["physical_active_vector"],
            objective=trial.objective, ynorm=trial.scaling, scaled=sim.scaled_intensity,
            experimental=experimental.intensities, beams=np.asarray(experimental.beam_indices_hk),
            scan_shape=[experimental.scan_rows, experimental.scan_columns],
            scan_json=_json(asdict(trial.physical_state.scan)))

    def log(kind, *, vector, objective, ynorm, stage=None, event_loop=None, **extra):
        if (kind in ("trial", "initial_retained", "retained_boundary", "snapshot_verification", "verification")
                and latest_display is not None
                and latest_display["evaluation_index"] == cp["evaluation_count"]):
            extra["physical_active_vector"] = latest_display["physical_active_vector"]
            extra["completed_evaluation"] = True
        row = dict(global_evaluation_index=cp["evaluation_count"],
            stage=current_stage if stage is None else stage,
            loop_index=cp["loop_index"] if event_loop is None else event_loop,
            stage_local_index=local_index, classification=kind, objective="inf" if np.isposinf(objective) else float(objective),
            ynorm=None if ynorm is None else float(ynorm), active_vector=[float(v) if np.isfinite(v) else str(v) for v in np.asarray(vector)], elapsed_seconds=elapsed(),
            accepted=extra.pop("accepted", None), restart_count=cp["restart_count"], **extra)
        with (destination / "progress.jsonl").open("a", encoding="utf-8") as stream:
            stream.write(json.dumps(row, allow_nan=False) + "\n")
            stream.flush()
            os.fsync(stream.fileno())
        if progress_hook is not None:
            progress_hook(row)

    def before():
        if config.max_total_evaluations is not None and cp["evaluation_count"] >= config.max_total_evaluations - 1:
            raise _BudgetStop("evaluation budget reached (one final-verification slot reserved)")
        if config.max_wall_time is not None and elapsed() >= config.max_wall_time:
            raise _BudgetStop("wall-time budget reached")
        cp["evaluation_count"] += 1

    def evaluate(vector, *, final=False):
        if final:
            if config.max_total_evaluations is not None and cp["evaluation_count"] >= config.max_total_evaluations:
                raise _BudgetStop("no evaluation budget remains for final verification")
            cp["evaluation_count"] += 1
        else:
            before()
        trial = evaluate_refinement_trial(state, np.asarray(vector, dtype=np.float64), experimental,
            objective_mode=config.objective_mode, execution=cfg.execution,
            bounds_policy=config.bounds_policy, refinement_bounds=config.refinement_bounds)
        if isinstance(trial, BoundRejectedTrial):
            raise ValueError(f"retained state violates configured hard bounds: {trial.violations}")
        if not np.isfinite(trial.objective) or trial.objective > 1e20:
            raise ValueError("retained objective must be finite and <= source CHISMIN sentinel")
        cp["completed_evaluations"] += 1
        record_completed_evaluation(trial)
        return trial

    def retain(vector, trial):
        cp.update(retained_vector=np.asarray(vector).tolist(), retained_objective=float(trial.objective),
                  retained_ynorm=float(trial.scaling), retained_at_evaluation=cp["evaluation_count"],
                  retained_loop_index=cp["loop_index"], retained_stage_name=current_stage)

    def snapshot(label, trial):
        sim = _simulated(trial, experimental)
        metadata = dict(evaluation_index=cp["evaluation_count"], retained_at_evaluation=cp["retained_at_evaluation"],
            stage=cp["retained_stage_name"], loop_index=cp["retained_loop_index"], objective=float(trial.objective),
            ynorm=float(trial.scaling), retained_active_vector=cp["retained_vector"], classification="retained")
        _save_arrays(destination / f"{label}_snapshot.npz",
            evaluation_index=cp["evaluation_count"],
            physical_vector=[v for _, _, v in _entries(trial.physical_state)],
            raw=sim.raw_intensity,
            objective_intensity=sim.objective_intensity, scaled=sim.scaled_intensity,
            vector=np.asarray(cp["retained_vector"]), objective=sim.objective, ynorm=sim.ynorm,
            experimental=experimental.intensities, beams=np.asarray(experimental.beam_indices_hk),
            scan_shape=np.asarray([experimental.scan_rows, experimental.scan_columns]),
            scan_json=_json(asdict(trial.physical_state.scan)))
        _atomic_json(destination / f"{label}_snapshot.json", metadata)
        _preview(destination / f"{label}_preview", experimental, initial or sim, sim, metadata)
        return sim

    def key(physical, objective):
        return (tuple(float(v) for _, _, v in _entries(physical)), float(objective))

    def stage_start(cursor, trial):
        nonlocal current_stage
        current_stage = cfg.stages[cursor.stage_position] if cfg.loops > 0 else "EVALUATION"
        cp.update(loop_index=cursor.loop_index, stage_position=cursor.stage_position,
            stage_name=current_stage, next_vector=cursor.vector.tolist())
        retain(cursor.vector, trial)
        checkpoint("stage_start")
        log("retained_boundary", vector=cursor.vector, objective=trial.objective, ynorm=trial.scaling)

    def observe(kind, stage_index, name, loop, payload):
        nonlocal local_index, pending_snapshot, current_stage, observed_stage
        current_stage = name
        if kind == "evaluation":
            if observed_stage != (loop, stage_index):
                local_index = 0
                observed_stage = (loop, stage_index)
            local_index += 1
            cp["driver_evaluations"] += 1
            if isinstance(payload, BoundRejectedTrial):
                cp["bound_rejection_count"] += 1
                counts = cp["bound_rejections_by_stage"]
                counts[name] = counts.get(name, 0)+1
                log("bound_rejection", vector=payload.trial_vector, objective=np.inf, ynorm=None,
                    event_loop=loop, driver_stage_index=stage_index, accepted=False,
                    coordinate_mode=payload.coordinate_mode, reason="bound_rejection",
                    raw_trial_vector=payload.trial_vector.tolist() if name in ("AMOEBA", "SIMPARA") else None,
                    retained_state_unchanged=True, violations=[asdict(v) for v in payload.violations],
                    mapped_physical_parameters=json.loads(_json(asdict(payload.physical_state))))
                print(f"Warning: rejected bound trial {cp['evaluation_count']} in {name}: {payload.violations}", flush=True)
                checkpoint("bound_rejection")
                return
            if isinstance(payload, InvalidRefinementTrial):
                cp["invalid_trial_count"] += 1
                counts = cp["invalid_trials_by_stage"]
                counts[name] = counts.get(name, 0)+1
                reasons = cp["invalid_trial_reasons"]
                reason = payload.exception_type+": "+payload.reason
                reasons[reason] = reasons.get(reason, 0)+1
                log("invalid_trial", vector=payload.trial_vector, objective=np.inf, ynorm=None,
                    event_loop=loop, driver_stage_index=stage_index, accepted=False,
                    coordinate_mode=payload.coordinate_mode, exception_type=payload.exception_type,
                    reason=payload.reason, retained_state_unchanged=True,
                    mapped_physical_parameters=None if payload.physical_state is None else json.loads(_json(asdict(payload.physical_state))))
                print(f"Warning: rejected invalid trial {cp['evaluation_count']} in {name}: {reason}", flush=True)
                checkpoint("invalid_trial")
                return
            cp["completed_evaluations"] += 1
            record_completed_evaluation(payload)
            stage_values[key(payload.physical_state, payload.objective)] = float(payload.scaling)
            log("trial", vector=[v for _, _, v in _entries(payload.physical_state)],
                objective=payload.objective, ynorm=payload.scaling, driver_stage_index=stage_index,
                event_loop=loop, vector_representation="mapped physical active values")
            if cp["evaluation_count"] % config.checkpoint_every == 0:
                checkpoint()
            if config.snapshot_every and cp["evaluation_count"] % config.snapshot_every == 0:
                pending_snapshot = True
        else:
            physical = apply_refinement_parameters(state, payload.retained_vector)
            ynorm = stage_values[key(physical, payload.retained_objective)]
            cp.update(retained_vector=payload.retained_vector.tolist(),
                retained_objective=float(payload.retained_objective), retained_ynorm=ynorm,
                retained_at_evaluation=cp["evaluation_count"], retained_loop_index=loop, retained_stage_name=name)
            cp["stages"].append(dict(optimizer=name, loop_index=loop,
                input_vector=payload.input_vector.tolist(), input_objective=float(payload.input_objective),
                retained_vector=payload.retained_vector.tolist(), retained_objective=float(payload.retained_objective),
                evaluations=payload.optimizer_result.evaluations, seed=payload.seed))
            position = (stage_index-1) % len(cfg.stages) + 1
            next_loop = loop + int(position == len(cfg.stages))
            cp.update(loop_index=next_loop, stage_position=position % len(cfg.stages),
                next_vector=pack_refinement_parameters(physical).tolist(), completed=next_loop > cfg.loops)
            checkpoint("stage_complete")
            log("stage_retained", vector=payload.retained_vector,
                objective=payload.retained_objective, ynorm=ynorm, event_loop=loop,
                selection="driver-confirmed stage result")
            stage_values.clear()
            if pending_snapshot:
                check = evaluate(cp["retained_vector"])
                verify(check)
                log("snapshot_verification", vector=cp["retained_vector"],
                    objective=check.objective, ynorm=check.scaling, event_loop=loop)
                snapshot(f"retained_{cp['evaluation_count']:06d}", check)
                pending_snapshot = False

    def verify(trial):
        if trial.objective != cp["retained_objective"] or trial.scaling != cp["retained_ynorm"]:
            raise ValueError("fresh retained evaluation disagrees with saved objective/YNORM")

    try:
        checkpoint("starting")
        try:
            if cp["retained_objective"] is None:
                vector = pack_refinement_parameters(state)
                first = evaluate(vector)
                initial = _simulated(first, experimental)
                np.save(destination / "initial_active_vector.npy", vector)
                np.save(destination / "initial_simulated_intensity.npy", first.calculated_intensities)
                _atomic_json(destination / "initial_parameters.json", asdict(first.physical_state))
                cp["initial_objective"], cp["initial_ynorm"] = float(first.objective), float(first.scaling)
                _save_arrays(destination / "initial_snapshot.npz", raw=initial.raw_intensity,
                    objective_intensity=initial.objective_intensity, scaled=initial.scaled_intensity,
                    vector=vector, objective=initial.objective, ynorm=initial.ynorm)
                retain(vector, first)
                checkpoint("initial")
                snapshot("initial", first)
                log("initial_retained", vector=vector, objective=first.objective, ynorm=first.scaling)
            else:
                with np.load(destination / "initial_snapshot.npz", allow_pickle=False) as saved:
                    initial_scan = apply_refinement_parameters(state, saved["vector"]).scan
                    initial = SimulatedCbedData(initial_scan, experimental.beam_indices_hk,
                        saved["raw"], saved["objective_intensity"], float(saved["ynorm"]),
                        experimental.gradient_option, float(saved["objective"]), config.objective_mode)
            if not cp["completed"]:
                restart = (LegacyRefinementRestart(cp["loop_index"], cp["stage_position"],
                    np.asarray(cp["next_vector"], dtype=np.float64)) if cfg.loops > 0 else None)
                run_legacy_refinement(state, experimental, cfg, objective_mode=config.objective_mode,
                    restart=restart, stage_start_callback=stage_start, progress_callback=observe,
                    before_evaluation=before, invalid_trial_policy=config.invalid_trial_policy,
                    bounds_policy=config.bounds_policy, refinement_bounds=config.refinement_bounds)
                cp["completed"] = True
        except KeyboardInterrupt:
            status, reason = "interrupted", "KeyboardInterrupt; latest confirmed retained boundary saved"
        except _BudgetStop as exc:
            status, reason = "limited", str(exc)
        checkpoint(status)
        # Even interrupted/budget-limited runs finalize the last safe boundary.
        if cp["retained_objective"] is None:
            summary = dict(status=status, reason=reason, output_directory=str(destination),
                           evaluations=cp["evaluation_count"], final_verification=False)
        else:
            final_trial = evaluate(cp["retained_vector"], final=True)
            verify(final_trial)
            log("verification", vector=cp["retained_vector"],
                objective=final_trial.objective, ynorm=final_trial.scaling,
                event_loop=cp["retained_loop_index"])
            final = snapshot("final", final_trial)
            vector = np.asarray(cp["retained_vector"], dtype=np.float64)
            _atomic_json(destination / "final_parameters.json", asdict(final_trial.physical_state))
            np.save(destination / "final_active_vector.npy", vector)
            np.save(destination / "final_simulated_intensity.npy", final.raw_intensity)
            np.save(destination / "final_scaled_intensity.npy", final.scaled_intensity)
            write_legacy_out_dat(destination / "final_out.dat", final)
            if export_path is not None and export_path != destination / "final_out.dat":
                export_path.parent.mkdir(parents=True, exist_ok=True)
                write_legacy_out_dat(export_path, final)
            history = tuple(LegacyRefinementStage(h["optimizer"], h["loop_index"],
                np.asarray(h["input_vector"]), h["input_objective"], np.asarray(h["retained_vector"]),
                h["retained_objective"], SimpleNamespace(evaluations=h["evaluations"]), h["seed"]) for h in cp["stages"])
            result = LegacyRefinementResult(pack_refinement_parameters(state), cp["initial_objective"],
                vector, final.objective, final_trial.physical_state, cfg.loops, history, cp["driver_evaluations"])
            from .visualization import plot_refinement_convergence, plot_parameter_trajectories
            figure = plot_refinement_convergence(result, objective_mode=config.objective_mode)
            figure.savefig(destination / "convergence.png", dpi=100)
            figure.clear()
            # Reuse the existing plot that explicitly distinguishes evaluated
            # trials (acceptance unknown) from confirmed retained boundaries.
            from .visualization import plot_local_progress
            plot_rows = []
            invalid_indices = []
            with (destination / "progress.jsonl").open(encoding="utf-8") as stream:
                for line in stream:
                    row = json.loads(line)
                    if row["classification"] in ("invalid_trial", "bound_rejection"):
                        invalid_indices.append(row["global_evaluation_index"])
                    if row["classification"] not in ("trial", "stage_retained"):
                        continue
                    stage_position = cfg.stages.index(row["stage"]) if row["stage"] in cfg.stages else 0
                    stage_number = (row["restart_count"]*max(cfg.loops, 1)+row["loop_index"]-1)*max(len(cfg.stages), 1)+stage_position+1
                    plot_rows.append(SimpleNamespace(kind="evaluation" if row["classification"] == "trial" else "retained",
                        evaluation_index=row["global_evaluation_index"], objective=row["objective"],
                        stage_index=stage_number, stage_name=row["stage"], loop_index=row["loop_index"]))
            figure = plot_local_progress(SimpleNamespace(progress=plot_rows, objective_mode=config.objective_mode))
            if invalid_indices:
                ax = figure.axes[0]
                ax.scatter(invalid_indices, [.96]*len(invalid_indices), marker="x", color="red",
                    transform=ax.get_xaxis_transform(), label="Rejected: numerical or bound (+inf; marker at plot top)")
                # Mixed axes/data transforms do not autoscale the x limits.
                left, right = ax.get_xlim()
                ax.set_xlim(min(left, min(invalid_indices)-.1),
                            max(right, max(invalid_indices)+.1))
                ax.legend()
            figure.savefig(destination / "trial_progress.png", dpi=100)
            figure.clear()
            if labels:
                figure = plot_parameter_trajectories(result, layout)
                stage_labels = ["Initial"] + [f"{h.optimizer}\nloop {h.loop_index}" for h in history]
                for ax in figure.axes:
                    ax.set_xticks(range(len(stage_labels)), stage_labels, rotation=30, ha="right")
                    for i in range(1, len(history)):
                        if history[i].loop_index != history[i-1].loop_index:
                            ax.axvline(i+.5, color="grey", linestyle=":")
                figure.savefig(destination / "parameter_trajectories.png", dpi=90)
                figure.clear()
            improvement = cp["initial_objective"] - float(final.objective)
            summary = dict(status=status, reason=reason, objective_mode=config.objective_mode,
                execution_backend=cfg.execution.execution_backend, batch_size=cfg.execution.batch_size,
                optimizer_mode=cfg.method, stages=cfg.stages, seed=cfg.seed,
                invalid_trial_policy=config.invalid_trial_policy, invalid_trial_count=cp["invalid_trial_count"],
                invalid_trials_by_stage=cp["invalid_trials_by_stage"], invalid_trial_reasons=cp["invalid_trial_reasons"],
                bounds_policy=config.bounds_policy, refinement_bounds=asdict(config.refinement_bounds),
                bound_rules=None if config.bound_rules is None else asdict(config.bound_rules),
                bound_rejection_count=cp["bound_rejection_count"], bound_rejections_by_stage=cp["bound_rejections_by_stage"],
                bound_policy_provenance="Hard bounds are a user-selected production refinement policy, not legacy Fortran behavior.",
                initial_objective=cp["initial_objective"], final_objective=float(final.objective),
                improvement=improvement, relative_improvement=improvement/cp["initial_objective"] if cp["initial_objective"] else None,
                initial_ynorm=cp["initial_ynorm"], final_ynorm=float(final.ynorm),
                evaluations=cp["evaluation_count"], completed_evaluations=cp["completed_evaluations"],
                driver_evaluations=cp["driver_evaluations"], runtime_seconds=elapsed(),
                output_directory=str(destination), report_path=str(destination / "report.html"),
                final_verification=True, restart_count=cp["restart_count"], restart_policy=RESTART_POLICY,
                parameter_changes=(vector-result.initial_vector).tolist())
            from .validation_report import generate_refinement_validation_report
            generate_refinement_validation_report(destination / "report.html",
                RefinementValidation(result, CbedComparison(experimental, initial, final)),
                layout=layout, config=cfg, reference_input_name=", ".join(files),
                limitations=(POLICY if config.objective_mode == "masked" else "Legacy objective unchanged.",
                    RESTART_POLICY, f"Run status: {status}. {reason}",
                    "Completed configured stages are not proof of scientific convergence."))
            report_path = destination / "report.html"
            text = report_path.read_text(encoding="utf-8")
            summary["runtime_seconds"] = elapsed()
            active_targets = {(p.category, p.target) for p in layout.parameters}
            bound_table = ("<h2>User-selected physical hard bounds</h2><p>Production policy, not legacy Fortran behavior. "
                "Endpoints are inclusive; inactive targets are metadata only.</p>"
                "<table><tr><th>Category</th><th>Target</th><th>Lower</th><th>Upper</th><th>Active</th></tr>" +
                "".join("<tr>" + "".join("<td>"+escape(str(v))+"</td>" for v in
                    (b.category, b.target, b.lower, b.upper, (b.category, b.target) in active_targets)) + "</tr>"
                    for b in config.refinement_bounds.entries) + "</table>")
            extra = ("<h2>Production run / restart summary</h2><pre>" + escape(_json(summary)) + "</pre>"
                "<p>Invalid numerical trials are rejected explicitly; acceptance of finite trials is unknown in callbacks. "
                "Only driver-confirmed boundaries appear on retained trajectories.</p>"
                + bound_table +
                "<h2>Evaluated trials versus retained stage states</h2><img src='trial_progress.png' "
                "alt='Evaluated trials, acceptance unknown, and retained stage markers'>"
                "<p><a href='initial_preview/index.html'>Initial all-beam gallery</a> | "
                "<a href='final_preview/index.html'>Final all-beam gallery</a> | "
                "<a href='progress.jsonl'>Append-only progress</a></p>"
                "<p>Parameter changes (same source-ordered labels as the parameter table):</p><pre>" +
                escape(_json(dict(zip(labels, summary["parameter_changes"])))) + "</pre>")
            report_path.write_text(text.replace("</html>", extra + "</html>"), encoding="utf-8")
        summary["runtime_seconds"] = elapsed()
        _atomic_json(destination / "summary.json", summary)
        checkpoint(status)
        return FullRefinementRun(destination, status, summary)
    except KeyboardInterrupt:
        checkpoint("interrupted_during_finalization")
        summary = dict(status="interrupted", reason="interrupted during finalization; resume to regenerate outputs",
                       evaluations=cp["evaluation_count"], output_directory=str(destination))
        _atomic_json(destination / "summary.json", summary)
        return FullRefinementRun(destination, "interrupted", summary)
    except _BudgetStop as exc:
        checkpoint("limited")
        summary = dict(status="limited", reason=str(exc), evaluations=cp["evaluation_count"],
                       output_directory=str(destination), final_verification=False)
        _atomic_json(destination / "summary.json", summary)
        return FullRefinementRun(destination, "limited", summary)
    except Exception as exc:
        checkpoint("failed")
        _atomic_json(destination / "summary.json", dict(status="failed", reason=f"{type(exc).__name__}: {exc}",
            evaluations=cp["evaluation_count"], output_directory=str(destination), final_verification=False))
        raise
    finally:
        lock.unlink()
