"""Read-only physics diagnosis of seed-2026 first ANNEAL proposal.

The isolation runs use the first 64 real tilts, not a shortened slice sequence.
Only scoped observation wrappers are installed; no kernel is replaced.
"""
from contextlib import ExitStack
from dataclasses import asdict, replace
import json
from pathlib import Path
import sys
from unittest.mock import patch

import numpy as np
from qcbed.anneal import _normal
from qcbed import ForwardExecutionConfig
from qcbed.refinement_driver import parse_legacy_refinement
from qcbed.refinement import pack_refinement_parameters, apply_refinement_parameters, build_refinement_parameter_layout
from qcbed.refinement_io import _CaptureRecords
from qcbed.legacy_io import _parse_for007_records
from qcbed.experimental import parse_experimental_intensities
import qcbed.refinement_trial as trial
import qcbed.simulation as simulation
import qcbed.tilt_series as tilts
from qcbed.full_runner import FullRefinementConfig, run_full_refinement


class ForwardDone(Exception):
    pass


def stats(a):
    a = np.asarray(a)
    finite = np.isfinite(a)
    with np.errstate(all="ignore"):
        magnitudes = np.abs(a[finite])
    def number(x):
        return float(x) if np.isfinite(x) else str(x)
    return dict(shape=list(a.shape), nonfinite=int(a.size-np.count_nonzero(finite)),
        abs_min=number(magnitudes.min()) if magnitudes.size else None,
        abs_max=number(magnitudes.max()) if magnitudes.size else None)


def diagnose(state, vector, experimental):
    info = dict(arrays=[])
    def record(name, a):
        s = stats(a)
        if s["nonfinite"] and "first_nonfinite" not in info:
            info["first_nonfinite"] = dict(name=name, **s)
        return s
    def wrap(name, fn):
        def measured(*a, **k):
            value = fn(*a, **k)
            info["arrays"].append(dict(name=name, **record(name, value)))
            return value
        return measured
    old_trace = sys.gettrace()
    last_finite = None
    last_identity = None
    def trace(frame, event, arg):
        nonlocal last_finite, last_identity
        if frame.f_code is not tilts._chunk.__code__:
            return None
        if event == "line" and "wave" in frame.f_locals:
            wave = frame.f_locals["wave"]
            identity = id(wave)
            if identity != last_identity:
                last_identity = identity
                s = stats(wave)
                occurrence = frame.f_locals.get("occurrence", 0)
                if s["nonfinite"] and "first_nonfinite" not in info:
                    info["first_nonfinite"] = dict(name="wave after grating multiplication", occurrence_zero_based=occurrence,
                        slice_type_zero_based=frame.f_locals.get("index"), previous_finite_wave=last_finite, **s)
                if not s["nonfinite"]:
                    last_finite = dict(occurrence_zero_based=occurrence, **s)
        return trace
    original_forward = trial.simulate_tilt_series_batched
    def forward(prepared, values, **kw):
        info["sequence_length"] = len(prepared.parameters.slice_sequence_zero_based)
        info["gratings"] = [stats(s.phase_grating) for s in prepared.slice_types]
        result = original_forward(prepared, values[:64], **kw)
        info["forward_amplitudes"] = stats(result.amplitudes)
        info["forward_intensities"] = stats(result.intensities)
        raise ForwardDone
    with ExitStack() as stack:
        for module in (simulation, trial):
            for name in ("elastic_structure_factors", "absorptive_structure_factors", "phase_grating_from_structure_factors"):
                stack.enter_context(patch.object(module, name, wrap(module.__name__+"."+name, getattr(module, name))))
        for name in ("batched_wave_real_to_reciprocal", "batched_wave_reciprocal_to_real"):
            original = getattr(tilts, name)
            def fft(*a, _fn=original, _name=name, **k):
                record(_name+" input", a[0])
                value = _fn(*a, **k)
                record(_name+" output", value)
                return value
            stack.enter_context(patch.object(tilts, name, fft))
        stack.enter_context(patch.object(trial, "simulate_tilt_series_batched", forward))
        sys.settrace(trace)
        try:
            trial.evaluate_refinement_trial(state, vector, experimental, objective_mode="masked",
                execution=ForwardExecutionConfig("batched", 64), invalid_trial_policy="raise")
        except ForwardDone:
            info["status"] = "finite first-64-tilt forward; full sequence retained; objective not evaluated"
        except (ValueError, FloatingPointError, OverflowError) as exc:
            info.update(status="numerical failure", exception_type=type(exc).__name__, reason=str(exc))
        finally:
            sys.settrace(old_trace)
    return info


def main():
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--continuation-only", action="store_true", help="reuse completed category diagnostics")
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    out = root / "tmp" / "invalid_trial_audit"
    out.mkdir(exist_ok=True)
    source = root / "tests/fixtures/legacy_reference" / "FOR007.DAT"
    state, original = parse_legacy_refinement(source, seed=2026)
    records = _CaptureRecords(source); _parse_for007_records(records)
    beam = records.values["output beam count and indices"]
    experimental = parse_experimental_intensities(source.parent / "in.dat",
        beam_indices_hk=tuple(zip(beam[1::2], beam[2::2])), scan=state.scan)
    initial = pack_refinement_parameters(state)
    stage_seed = int(np.random.SeedSequence([2026, 1, 0]).generate_state(1, dtype=np.uint64)[0])
    rng = np.random.Generator(np.random.PCG64(stage_seed))
    raw = np.array([_normal(rng, x, x*100.) for x in initial])
    projected = pack_refinement_parameters(apply_refinement_parameters(state, raw))
    layout = build_refinement_parameter_layout(state)
    results = dict(interpreter=sys.executable, stage_seed=stage_seed, initial=initial.tolist(),
        raw_proposal=raw.tolist(), projected_proposal=projected.tolist(), layout=asdict(layout), cases={})
    if args.continuation_only:
        results = json.loads((out / "diagnosis.json").read_text(encoding="utf-8"))
    for category in (() if args.continuation_only else (None, "absorption", "debye_waller", "reflection_multiplier", "depth", "scan_coordinate")):
        vector = projected.copy()
        if category:
            for p in layout.parameters:
                if p.category == category:
                    vector[p.vector_index] = initial[p.vector_index]
        name = category or "original_failed_proposal"
        print(f"Diagnosing: {name}", flush=True)
        results["cases"][name] = diagnose(state, vector, experimental)
        (out / "diagnosis.json").write_text(json.dumps(results, indent=2), encoding="utf-8")
        print(json.dumps(results["cases"][name]), flush=True)
    # Same source optimizer settings as failed production run; only a short
    # external evaluation budget stops this validation, not a new temperature.
    config = replace(original, method="ANSI", loops=1, maxfsimp=20, maxfpara=3,
                     execution=ForwardExecutionConfig("batched", 64))
    from datetime import datetime
    directory = out / ("continuation_"+datetime.now().strftime("%Y%m%d_%H%M%S"))
    result = run_full_refinement(state, experimental,
        FullRefinementConfig(config, invalid_trial_policy="reject", snapshot_every=None,
            checkpoint_every=1, max_total_evaluations=5, max_wall_time=600.),
        output_directory=directory, input_paths=(source, source.parent / "in.dat"))
    results["continuation"] = result.summary
    (out / "diagnosis.json").write_text(json.dumps(results, indent=2), encoding="utf-8")
    print(json.dumps(result.summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
