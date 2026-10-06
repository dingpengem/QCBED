from dataclasses import replace
import json
from pathlib import Path

import numpy as np
import pytest

from qcbed import ForwardExecutionConfig, InvalidRefinementTrial, select_refinement_categories
from qcbed.refinement import pack_refinement_parameters, apply_refinement_parameters, build_refinement_parameter_layout
from qcbed.refinement_trial import evaluate_refinement_trial
from qcbed.refinement_driver import run_legacy_refinement, parse_legacy_refinement
from qcbed.full_runner import FullRefinementConfig, run_full_refinement, load_refinement_checkpoint
from test_refinement_driver import problem, config
from test_refinement import data
from test_full_runner import cheap_render


@pytest.mark.parametrize("reason", ["phase grating must remain finite in float64", "real waves must be finite",
    "propagated reciprocal waves must be finite", "beam intensities must remain finite in float64"])
def test_numerical_failure_policy_and_immutability(monkeypatch, reason):
    import qcbed.refinement_trial as module
    s = problem(); vector = pack_refinement_parameters(s); baseline = repr(s)
    before = vector.copy()
    def broken(*args, **kwargs):
        raise ValueError(reason)
    monkeypatch.setattr(module, "simulate_tilt_series_batched", broken)
    execution = ForwardExecutionConfig("batched", 3)
    with pytest.raises(ValueError, match=reason):
        evaluate_refinement_trial(s, vector, data(), execution=execution, invalid_trial_policy="raise")
    seen = []
    rejected = evaluate_refinement_trial(s, vector, data(), execution=execution,
        invalid_trial_policy="reject", evaluation_callback=seen.append)
    assert isinstance(rejected, InvalidRefinementTrial) and np.isposinf(rejected.objective)
    assert rejected.exception_type == "ValueError" and rejected.reason == reason
    assert rejected.physical_state is not None and seen == [rejected]
    assert not hasattr(rejected, "amplitudes") and not hasattr(rejected, "scaling")
    np.testing.assert_array_equal(vector, before)
    assert repr(s) == baseline
    with pytest.raises(ValueError):
        rejected.trial_vector.setflags(write=True)


@pytest.mark.parametrize("failure", [FloatingPointError("overflow encountered in multiply"),
    ValueError("unsupported geometry"), RuntimeError("programming error")])
def test_only_numerical_errors_rejected(monkeypatch, failure):
    import qcbed.refinement_trial as module
    def broken(*a, **kw):
        raise failure
    monkeypatch.setattr(module, "_evaluate_refinement_trial", broken)
    if isinstance(failure, FloatingPointError):
        assert np.isposinf(evaluate_refinement_trial(problem(), np.array([2.]), data(), invalid_trial_policy="reject").objective)
    else:
        with pytest.raises(type(failure), match=str(failure)):
            evaluate_refinement_trial(problem(), np.array([2.]), data(), invalid_trial_policy="reject")


@pytest.mark.parametrize("field", ["objective", "scaling"])
def test_nonfinite_result_policy(monkeypatch, field):
    import qcbed.refinement_trial as module
    baseline = evaluate_refinement_trial(problem(), np.array([2.]), data())
    monkeypatch.setattr(module, "_evaluate_refinement_trial", lambda *a, **k: replace(baseline, **{field: np.nan}))
    assert np.isnan(getattr(evaluate_refinement_trial(problem(), np.array([2.]), data()), field))
    assert np.isposinf(evaluate_refinement_trial(problem(), np.array([2.]), data(), invalid_trial_policy="reject").objective)


@pytest.mark.parametrize("method", ["ANNE", "SIMP", "PARA", "SPAR", "ANSI", "PAAN"])
def test_optimizer_continues_and_keeps_valid_state(monkeypatch, method):
    import qcbed.refinement_trial as module
    original = module._evaluate_refinement_trial
    s = problem(); before = repr(s)
    cfg = config(method, maxfsimp=3, maxfpara=0, execution=ForwardExecutionConfig("batched", 3))
    def run_once():
        count = 0
        def faulty(*a, **k):
            nonlocal count
            count += 1
            if count == 3:
                raise ValueError("real waves must be finite")
            return original(*a, **k)
        monkeypatch.setattr(module, "_evaluate_refinement_trial", faulty)
        observed = []
        result = run_legacy_refinement(s, data(), cfg, invalid_trial_policy="reject",
            progress_callback=lambda *x: observed.append(x))
        assert count > 3 and np.isfinite(result.final_objective)
        invalid = [x for x in observed if isinstance(x[-1], InvalidRefinementTrial)]
        assert len(invalid) == 1
        assert all(np.isfinite(h.retained_objective) for h in result.stage_history)
        return result
    a, b = run_once(), run_once()
    assert a.final_objective == b.final_objective
    np.testing.assert_array_equal(a.final_vector, b.final_vector)
    assert repr(s) == before


def test_invalid_proposal_never_updates_anneal_retained():
    from qcbed.anneal import _run
    values = iter([10., np.inf, 11.])
    result = _run(np.array([2.]), lambda v: v.copy(), lambda v: next(values),
        max_iterations=1, maximum_temperature=1., minimum_temperature=.1, cooling_rate=.9,
        rng=np.random.Generator(np.random.PCG64(1)), allow_invalid=True)
    np.testing.assert_array_equal(result.best_vector, [2.])
    assert result.best_objective == 10. and result.evaluations == 3


@pytest.mark.parametrize("evaluation_cap", [None, 5])
def test_full_runner_invalid_json_and_summary(cheap_render, monkeypatch, tmp_path, capsys, evaluation_cap):
    import qcbed.refinement_trial as module
    original = module._evaluate_refinement_trial
    count = 0
    def faulty(*a, **kw):
        nonlocal count
        count += 1
        if count == 4:
            raise FloatingPointError("overflow encountered in multiply")
        return original(*a, **kw)
    monkeypatch.setattr(module, "_evaluate_refinement_trial", faulty)
    from matplotlib.figure import Figure
    savefig = Figure.savefig
    def check_marker(self, path, **kwargs):
        if Path(path).name == "trial_progress.png":
            ax = self.axes[0]
            assert ax.get_xlim()[0] < 4 < ax.get_xlim()[1]
        return savefig(self, path, **kwargs)
    monkeypatch.setattr(Figure, "savefig", check_marker)
    c = FullRefinementConfig(config("ANSI", maxfsimp=0, execution=ForwardExecutionConfig("batched", 3)),
        snapshot_every=None, checkpoint_every=1, max_total_evaluations=evaluation_cap)
    folder = tmp_path / "run"
    result = run_full_refinement(problem(), data(), c, output_directory=folder)
    assert result.status == ("completed" if evaluation_cap is None else "limited")
    assert result.summary["invalid_trial_count"] == 1
    assert result.summary["invalid_trials_by_stage"] == {"ANNEAL": 1}
    rows = [json.loads(line) for line in (folder / "progress.jsonl").read_text().splitlines()]
    row = next(r for r in rows if r["classification"] == "invalid_trial")
    assert row["objective"] == "inf" and row["ynorm"] is None and row["accepted"] is False
    assert row["retained_state_unchanged"] and row["mapped_physical_parameters"] is not None
    assert any(r["global_evaluation_index"] > row["global_evaluation_index"] for r in rows)
    assert "rejected invalid trial" in capsys.readouterr().out
    assert np.isfinite(load_refinement_checkpoint(folder)["retained_objective"])
    assert result.summary["completed_evaluations"] == result.summary["evaluations"]-1
    assert "invalid_trial_count" in (folder / "report.html").read_text()


def test_subset_preserves_excluded_reference_parameters():
    state, cfg = parse_legacy_refinement(Path(__file__).parents[1] / "tests/fixtures/legacy_reference" / "FOR007.DAT", seed=2026)
    original = repr(state)
    subset, controls = select_refinement_categories(state, cfg, ("depth", "scan_coordinate"))
    layout = build_refinement_parameter_layout(subset)
    assert len(layout.parameters) == 7 and len(controls.plamda) == 7
    assert {p.category for p in layout.parameters} == {"depth", "scan_coordinate"}
    v = pack_refinement_parameters(subset).copy(); v *= 1.01
    changed = apply_refinement_parameters(subset, v)
    assert changed.absorption == state.absorption
    assert changed.regions[0].debye_waller == state.regions[0].debye_waller
    assert [(r.elastic, r.absorptive) for r in changed.reflections[0]] == [(r.elastic, r.absorptive) for r in state.reflections[0]]
    assert changed.model == subset.model
    assert repr(state) == original


def test_callback_errors_not_rejected():
    def callback(result):
        raise ValueError("my callback must be finite")
    with pytest.raises(ValueError, match="my callback"):
        evaluate_refinement_trial(problem(), np.array([2.]), data(), invalid_trial_policy="reject", evaluation_callback=callback)


@pytest.mark.parametrize("optimizer", ["amoeba", "paramin", "simpara"])
def test_rejected_score_cannot_replace_retained_vector_or_score(optimizer):
    calls = []
    def objective(vector):
        calls.append(vector.copy())
        if len(calls) == 2:
            return np.inf
        return 10. if np.array_equal(vector, [2.]) else 11.
    if optimizer == "amoeba":
        from qcbed.amoeba import _run
        result = _run(np.array([[2.], [3.]]), objective, ftol=0.,
                      max_evaluations=1, fmin=0., allow_invalid=True)
        np.testing.assert_array_equal(result.best_evaluated_vector, [2.])
        assert result.best_evaluated_objective == 10.
    elif optimizer == "paramin":
        from qcbed.paramin import _run
        result = _run(np.array([2.]), np.array([.9]), np.array([1.1]),
                      lambda v: v.copy(), objective, max_cycles=0,
                      tolerance=0., allow_invalid=True)
    else:
        from qcbed.simpara import _run
        result = _run(np.array([2.]), np.array([.1]), objective,
                      max_evaluations=1, high_temperature=1., low_temperature=.1,
                      cooling=.9, random_order=False,
                      rng=np.random.Generator(np.random.PCG64(1)), allow_invalid=True)
    np.testing.assert_array_equal(result.best_vector, [2.])
    assert result.best_objective == 10. and len(calls) > 2


@pytest.mark.parametrize("vector", [np.array([1., 2.]), np.array([1], dtype=np.int64)])
def test_malformed_vector_is_not_an_invalid_physical_trial(vector):
    with pytest.raises((ValueError, TypeError)):
        evaluate_refinement_trial(problem(), vector, data(), invalid_trial_policy="reject")


@pytest.mark.parametrize("value", [np.nan, np.inf, -np.inf])
def test_nonfinite_proposal_is_rejected(value):
    result = evaluate_refinement_trial(problem(), np.array([value]), data(), invalid_trial_policy="reject")
    assert isinstance(result, InvalidRefinementTrial) and np.isposinf(result.objective)
    assert result.physical_state is None


def test_bad_experimental_variance_is_not_hidden(monkeypatch):
    import qcbed.refinement_trial as module
    def broken(*a, **kw):
        raise ValueError("masked objective: nonzero-weight point has invalid variance")
    monkeypatch.setattr(module, "_evaluate_refinement_trial", broken)
    with pytest.raises(ValueError, match="invalid variance"):
        evaluate_refinement_trial(problem(), np.array([2.]), data(), invalid_trial_policy="reject")
