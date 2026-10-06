"""Short real-data reproduction of the historical 840 A, seed-2026 proposal.

The current input is read only. Restore only historical depth in memory; no long
fit, output/report timings, or physics substitutions enter measured evaluations.
"""
from dataclasses import replace
from datetime import datetime
import json
from pathlib import Path
from time import perf_counter
from unittest.mock import patch
import warnings
import sys

import numpy as np

from qcbed import ForwardExecutionConfig, RefinementBound, RefinementBounds, BoundRejectedTrial
from qcbed.anneal import run_legacy_anneal, _normal
from qcbed.refinement import pack_refinement_parameters, apply_refinement_parameters
from qcbed.refinement_driver import parse_legacy_refinement
from qcbed.refinement_bounds import refinement_bound_policy
from qcbed.refinement_trial import evaluate_refinement_trial
from qcbed.refinement_io import _CaptureRecords
from qcbed.legacy_io import _parse_for007_records
from qcbed.experimental import parse_experimental_intensities
import qcbed.refinement_trial as trial


def main():
    root = Path(__file__).resolve().parents[1]
    source = root/"tests/fixtures/legacy_reference/FOR007.DAT"
    state, _ = parse_legacy_refinement(source, seed=2026)
    current_depth = state.regions[0].depth_angstrom
    state = replace(state, regions=(replace(state.regions[0], depth_angstrom=840.),)+state.regions[1:])
    original_state = repr(state)
    records = _CaptureRecords(source); _parse_for007_records(records)
    b = records.values["output beam count and indices"]
    experimental = parse_experimental_intensities(source.parent/"in.dat",
        beam_indices_hk=tuple(zip(b[1::2], b[2::2])), scan=state.scan)
    initial = pack_refinement_parameters(state)
    assert len(initial) == 18
    stage_seed = int(np.random.SeedSequence([2026, 1, 0]).generate_state(1, dtype=np.uint64)[0])
    rng = np.random.Generator(np.random.PCG64(stage_seed))
    raw = np.array([_normal(rng, x, x*100.) for x in initial])
    projected = pack_refinement_parameters(apply_refinement_parameters(state, raw))
    assert projected[0] == 27925.729123683803
    bounds = RefinementBounds((RefinementBound("depth", (0,), 600., 1200.),))
    kwargs = dict(objective_mode="masked", execution=ForwardExecutionConfig("batched", 64),
                  invalid_trial_policy="reject")
    with patch.object(trial, "prepare_multislice_sample", side_effect=AssertionError("forward entered")) as spy:
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            start = perf_counter()
            bounded = evaluate_refinement_trial(state, raw, experimental, bounds_policy="reject",
                refinement_bounds=bounds, **kwargs)
            bounded_seconds = perf_counter()-start
        spy.assert_not_called()
        assert not caught and isinstance(bounded, BoundRejectedTrial)
    print(f"Bound rejection: {bounded_seconds:.6f} s; no preparation/forward calls", flush=True)
    start = perf_counter()
    numerical = evaluate_refinement_trial(state, raw, experimental, **kwargs)
    numerical_seconds = perf_counter()-start
    assert not isinstance(numerical, BoundRejectedTrial) and np.isposinf(numerical.objective)
    print(f"Propagation failure: {numerical_seconds:.6f} s; {numerical.reason}", flush=True)
    evaluations = []
    with refinement_bound_policy(state, bounds_policy="reject", refinement_bounds=bounds,
                                  invalid_trial_policy="reject"):
        with patch.object(trial, "prepare_multislice_sample", wraps=trial.prepare_multislice_sample) as forward:
            start = perf_counter()
            result = run_legacy_anneal(state, initial, experimental, max_iterations=1,
                seed=stage_seed, maximum_temperature=100., minimum_temperature=1., cooling_rate=.9,
                evaluation_callback=evaluations.append, **kwargs)
            continuation_seconds = perf_counter()-start
            forward_calls = forward.call_count
    assert len(evaluations) == 3 and isinstance(evaluations[1], BoundRejectedTrial)
    assert isinstance(evaluations[2], BoundRejectedTrial) and forward_calls == 1
    np.testing.assert_array_equal(evaluations[1].trial_vector, projected)
    np.testing.assert_array_equal(result.best_vector, initial)
    assert result.best_objective == evaluations[0].objective and repr(state) == original_state
    out = root/"tmp/refinement_bounds"/datetime.now().strftime("%Y%m%d_%H%M%S_%f")
    out.mkdir(parents=True)
    summary = dict(interpreter=sys.executable, current_input_depth=current_depth,
        historical_depth_restored_in_memory=840., seed=2026, stage_seed=stage_seed,
        raw_proposal=raw.tolist(), mapped_proposal=projected.tolist(),
        bounded_seconds=bounded_seconds, propagation_failure_seconds=numerical_seconds,
        seconds_saved=numerical_seconds-bounded_seconds, propagation_failure_reason=numerical.reason,
        bounded_forward_calls=0, bounded_warning_count=0, continuation_seconds=continuation_seconds,
        continuation_evaluations=result.evaluations, continuation_forward_calls=forward_calls,
        bound_rejections=sum(isinstance(r, BoundRejectedTrial) for r in evaluations),
        retained_unchanged=True, initial_objective=evaluations[0].objective, retained_objective=result.best_objective,
        scope="One historical proposal timing per policy; initialization plus two source ANNEAL proposals. No convergence claim.")
    (out/"summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2), flush=True)
    print(f"Output: {out}", flush=True)


if __name__ == "__main__":
    main()
