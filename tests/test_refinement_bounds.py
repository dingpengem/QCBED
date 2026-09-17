from dataclasses import replace, FrozenInstanceError
import json
from pathlib import Path
from unittest.mock import Mock

import numpy as np
import pytest

from qcbed import RefinementBound, RefinementBounds, BoundRejectedTrial, scan_relative_bounds
from qcbed import RefinementBoundRules, resolve_refinement_bounds
from qcbed.refinement import RefinementParameterSpec, RefinementParameterLayout
from qcbed.refinement_bounds import refinement_bound_policy
from qcbed.refinement_trial import evaluate_refinement_trial
from qcbed.refinement_driver import run_legacy_refinement, parse_legacy_refinement
from qcbed.refinement import apply_refinement_parameters, pack_refinement_parameters
from qcbed.full_runner import FullRefinementConfig, run_full_refinement, load_refinement_checkpoint
from qcbed.refinement_subset import select_refinement_categories
from test_refinement_driver import problem, config
from test_refinement import data
from test_full_runner import cheap_render


def depth(lower=None, upper=None):
    return RefinementBounds((RefinementBound("depth", (0,), lower, upper),))


@pytest.mark.parametrize("value,lower,upper,reject", [
    (1., 2., None, True), (3., None, 2., True),
    (2., 2., 3., False), (3., 2., 3., False),
    (-2., 2., 3., False), (-4., 2., 3., True)])
def test_physical_bounds_and_forward_gate(monkeypatch, value, lower, upper, reject):
    import qcbed.refinement_trial as module
    original = module.prepare_multislice_sample
    forward = Mock(wraps=original)
    monkeypatch.setattr(module, "prepare_multislice_sample", forward)
    s = problem(); baseline = repr(s); vector = np.array([value])
    result = evaluate_refinement_trial(s, vector, data(), bounds_policy="reject",
                                      refinement_bounds=depth(lower, upper))
    assert isinstance(result, BoundRejectedTrial) == reject
    assert forward.call_count == (0 if reject else 1)
    assert repr(s) == baseline and vector[0] == value
    if reject:
        assert np.isposinf(result.objective) and result.reason == "bound_rejection"
        assert result.violations[0].value == abs(value)
        assert result.physical_state.regions[0].depth_angstrom == abs(value)
        assert result.exception_type == ""


def test_none_unchanged_and_unspecified_unbounded():
    s = problem(); v = np.array([3.])
    baseline = evaluate_refinement_trial(s, v, data())
    for kw in (dict(bounds_policy="none", refinement_bounds=depth(1., 2.)),
               dict(bounds_policy="reject", refinement_bounds=RefinementBounds())):
        out = evaluate_refinement_trial(s, v, data(), **kw)
        assert out.objective == baseline.objective
        np.testing.assert_array_equal(out.amplitudes, baseline.amplitudes)


def test_simplex_coordinates_are_converted_first():
    # Existing source conversion: SCAL a=-.5, b=3; q=-.5 maps to depth=4.
    r = evaluate_refinement_trial(problem(), np.array([-.5]), data(),
        coordinate_mode="legacy_simplex", bounds_policy="reject", refinement_bounds=depth(1., 3.))
    assert isinstance(r, BoundRejectedTrial) and r.violations[0].value == 4.


def test_multiple_parameters_one_violation():
    s = replace(problem(), scan_masks=("YY", "NN", "NN"), ranges=())
    bounds = RefinementBounds((RefinementBound("depth", (0,), 1., 3.),
                              RefinementBound("scan_coordinate", (0, 0), -.1, .1)))
    r = evaluate_refinement_trial(s, np.array([2., .2, 100.]), data(),
                                 bounds_policy="reject", refinement_bounds=bounds)
    assert len(r.violations) == 1 and r.violations[0].category == "scan_coordinate"
    assert r.physical_state.scan.corner_a == (.2, 100.)


@pytest.mark.parametrize("args", [("depth", (0,), 3., 2.), ("depth", (0,), np.nan, None),
    ("depth", (0,), None, np.inf), ("depth", (0,), True, None),
    ("other", (), 1., 2.), ("depth", (-1,), 1., 2.), ("depth", (), 1., 2.),
    ("depth", (0,), None, None)])
def test_malformed_bound(args):
    with pytest.raises((ValueError, TypeError)):
        RefinementBound(*args)


def test_configuration_validation_and_immutability():
    with pytest.raises(ValueError, match="unknown bound"):
        evaluate_refinement_trial(problem(), np.array([2.]), data(), bounds_policy="reject",
            refinement_bounds=RefinementBounds((RefinementBound("depth", (1,), 1., 3.),)))
    with pytest.raises(ValueError, match="duplicate"):
        RefinementBounds(depth(1., 3.).entries*2)
    with pytest.raises(ValueError, match="bounds_policy"):
        evaluate_refinement_trial(problem(), np.array([2.]), data(), bounds_policy="clip")
    b = depth(1., 3.)
    with pytest.raises(FrozenInstanceError):
        b.entries = ()


@pytest.mark.parametrize("mode", ["ANNE", "SIMP", "PARA", "SPAR", "ANSI", "PAAN"])
def test_optimizer_continuation_state_and_determinism(mode):
    s = problem(); before = repr(s)
    cfg = config(mode, maxfsimp=3, maxfpara=0)
    baseline = evaluate_refinement_trial(s, np.array([2.]), data()).objective
    def run():
        events = []
        result = run_legacy_refinement(s, data(), cfg, bounds_policy="reject",
            refinement_bounds=depth(2., 2.), progress_callback=lambda *a: events.append(a[-1]))
        rejected = [i for i, r in enumerate(events) if isinstance(r, BoundRejectedTrial)]
        assert rejected and rejected[0] < len(events)-1
        np.testing.assert_array_equal(result.final_vector, [2.])
        assert result.final_objective == baseline
        assert all(h.retained_objective == baseline for h in result.stage_history)
        return result, len(rejected)
    a, n = run(); b, m = run()
    assert n == m and a.final_objective == b.final_objective and repr(s) == before


def test_numeric_policy_remains_independent(monkeypatch):
    import qcbed.refinement_trial as module
    def broken(*args, **kwargs):
        raise FloatingPointError("overflow in simulation")
    monkeypatch.setattr(module, "prepare_multislice_sample", broken)
    with pytest.raises(FloatingPointError):
        run_legacy_refinement(problem(), data(), config(), bounds_policy="reject",
            refinement_bounds=depth(1., 3.), invalid_trial_policy="raise")
    # Scope restoration after exception: no bound policy leaks into later calls.
    monkeypatch.undo()
    assert not isinstance(evaluate_refinement_trial(problem(), np.array([4.]), data()), BoundRejectedTrial)


def test_scan_helper_negative_positive_zero_and_subset():
    s = replace(problem(), scan=replace(problem().scan, corner_a=(-.5, .5)),
                scan_masks=("YY", "YY", "YY"), ranges=())
    b = scan_relative_bounds(s, relative_fraction=.2, absolute_fallback=.01)
    assert [(v.lower, v.upper) for v in b.entries] == [(-.6, -.4), (.4, .6)]+[(-.01, .01)]*4
    sub, _ = select_refinement_categories(s, config(plamda=(.2,)*7), ("depth",))
    b.validate(sub)
    assert b.violations(sub) == ()  # Existing but inactive scan targets are metadata.


def test_runner_logs_reports_and_resume(cheap_render, tmp_path):
    s = problem()
    c = FullRefinementConfig(config("ANSI", maxfsimp=0), bounds_policy="reject",
        refinement_bounds=depth(2., 2.), snapshot_every=None, checkpoint_every=1,
        max_total_evaluations=5)
    folder = tmp_path/"run"
    first = run_full_refinement(s, data(), c, output_directory=folder)
    assert first.status == "limited"
    rows = [json.loads(x) for x in (folder/"progress.jsonl").read_text().splitlines()]
    rejection = next(r for r in rows if r["classification"] == "bound_rejection")
    assert rejection["objective"] == "inf" and rejection["retained_state_unchanged"]
    assert rejection["violations"][0]["category"] == "depth"
    assert rejection["mapped_physical_parameters"] and rejection["active_vector"]
    assert first.summary["bound_rejection_count"] == 1
    assert first.summary["invalid_trial_count"] == 0
    assert "bound_rejection_count" in (folder/"report.html").read_text()
    cp = load_refinement_checkpoint(folder)
    assert RefinementBounds.from_dict(cp["config"]["refinement_bounds"]) == c.refinement_bounds
    with pytest.raises(ValueError, match="original hard-bound"):
        run_full_refinement(s, data(), replace(c, refinement_bounds=depth(1., 3.)), output_directory=folder, resume=True)
    resumed = run_full_refinement(s, data(), replace(c, max_total_evaluations=100), output_directory=folder, resume=True)
    assert resumed.status == "completed"
    assert resumed.summary["bound_rejection_count"] >= first.summary["bound_rejection_count"]
    assert resumed.summary["final_objective"] == first.summary["initial_objective"]


def test_outside_start_fails_before_forward(cheap_render, monkeypatch, tmp_path):
    import qcbed.refinement_trial as module
    forward = Mock(side_effect=AssertionError("must not simulate"))
    monkeypatch.setattr(module, "prepare_multislice_sample", forward)
    with pytest.raises(ValueError, match="retained state violates"):
        run_full_refinement(problem(), data(), FullRefinementConfig(config(), bounds_policy="reject",
            refinement_bounds=depth(3., 4.)), output_directory=tmp_path/"run")
    forward.assert_not_called()


def layout_for(values):
    return RefinementParameterLayout(tuple(RefinementParameterSpec(i, c, t, v, None)
        for i, (c, t, v) in enumerate(values)))


def test_category_rule_expands_eight_reflections_and_specific_wins():
    layout = layout_for([("reflection_multiplier", (0, i, 0), 1.) for i in range(8)])
    rules = RefinementBoundRules(
        specific_bounds={("reflection_multiplier", (0, 2, 0)): (.8, 1.2)},
        category_bounds={"reflection_multiplier": (.5, 1.5)})
    resolved = resolve_refinement_bounds(layout, rules)
    assert len(resolved.entries) == 8
    for i, b in enumerate(resolved.entries):
        assert (b.lower, b.upper) == ((.8, 1.2) if i == 2 else (.5, 1.5))


def test_absolute_precedence_relative_signs_fallback_and_unbounded():
    layout = layout_for([("depth", (0,), 100.), ("scan_coordinate", (0, 0), -.5),
        ("scan_coordinate", (0, 1), .5), ("scan_coordinate", (1, 0), 0.),
        ("scan_coordinate", (1, 1), .01), ("voltage", (), 200.)])
    rules = RefinementBoundRules(category_bounds={"depth": (50., 200.)},
        category_relative_bounds={"depth": .1, "scan_coordinate": .2},
        relative_bound_absolute_fallback={"depth": 1., "scan_coordinate": .05})
    b = resolve_refinement_bounds(layout, rules)
    assert len(b.entries) == 5  # voltage is unbounded
    np.testing.assert_allclose([(x.lower, x.upper) for x in b.entries],
        [(50., 200.), (-.6, -.4), (.4, .6), (-.05, .05), (-.04, .06)])


def test_rule_resolution_uses_source_mapped_start():
    s = problem()
    s = replace(s, regions=(replace(s.regions[0], depth_angstrom=-2.),))
    b = resolve_refinement_bounds(s, RefinementBoundRules(category_relative_bounds={"depth": .2},
        relative_bound_absolute_fallback={"depth": .1}))
    assert (b.entries[0].lower, b.entries[0].upper) == (1.6, 2.4)


def test_all_reference_parameters_and_subset_without_per_target_configuration():
    state, cfg = parse_legacy_refinement(Path(__file__).parents[1]/"reference_data/FOR007.DAT", seed=2026)
    rules = RefinementBoundRules(category_bounds={"depth": (600., 1200.), "absorption": (.5, 1.5),
        "debye_waller": (.5, 1.5), "reflection_multiplier": (.5, 1.5)},
        category_relative_bounds={"scan_coordinate": .2},
        relative_bound_absolute_fallback={"scan_coordinate": .05})
    assert len(resolve_refinement_bounds(state, rules).entries) == 18
    sub, _ = select_refinement_categories(state, cfg, ("depth", "scan_coordinate"))
    resolved = resolve_refinement_bounds(sub, rules)
    assert len(resolved.entries) == 7
    assert {b.category for b in resolved.entries} == {"depth", "scan_coordinate"}


def test_rule_immutability_and_validation():
    source = {"depth": (1., 3.)}
    rules = RefinementBoundRules(category_bounds=source)
    source["depth"] = (0., 100.)
    assert rules.category_bounds == (("depth", (1., 3.)),)
    for kwargs in (dict(category_relative_bounds={"depth": .2}),
        dict(category_bounds={"unknown": (1., 2.)}),
        dict(category_relative_bounds={"depth": -.1}, relative_bound_absolute_fallback={"depth": .1})):
        with pytest.raises(ValueError):
            RefinementBoundRules(**kwargs)


def test_rules_and_resolved_bounds_persist_and_resume(cheap_render, tmp_path):
    rules = RefinementBoundRules(category_bounds={"depth": (2., 2.)})
    c = FullRefinementConfig(config("ANNE", maxfsimp=0), bounds_policy="reject", bound_rules=rules,
        refinement_bounds=resolve_refinement_bounds(problem(), rules), max_total_evaluations=5)
    folder = tmp_path/"rules"
    run_full_refinement(problem(), data(), c, output_directory=folder)
    saved = load_refinement_checkpoint(folder)["config"]
    assert RefinementBoundRules.from_dict(saved["bound_rules"]) == rules
    assert RefinementBounds.from_dict(saved["refinement_bounds"]) == c.refinement_bounds
    run_config = json.loads((folder/"run_config.json").read_text())
    assert run_config["bound_rules"] == saved["bound_rules"]
    assert "User-selected physical hard bounds" in (folder/"report.html").read_text()
    incompatible = RefinementBoundRules(specific_bounds={("depth", (0,)): (2., 2.)})
    with pytest.raises(ValueError, match="original hard-bound"):
        run_full_refinement(problem(), data(), replace(c, bound_rules=incompatible), output_directory=folder, resume=True)
    result = run_full_refinement(problem(), data(), replace(c, max_total_evaluations=100), output_directory=folder, resume=True)
    assert result.summary["refinement_bounds"] == {"entries": tuple(
        dict(category=b.category, target=b.target, lower=b.lower, upper=b.upper) for b in c.refinement_bounds.entries)}


def test_dry_run_displays_bounds_and_all_parameter_selection(monkeypatch, capsys):
    import importlib.util
    import sys
    path = Path(__file__).parents[1]/"examples/run_full_refinement.py"
    spec = importlib.util.spec_from_file_location("bound_example", path)
    example = importlib.util.module_from_spec(spec); spec.loader.exec_module(example)
    monkeypatch.setattr(sys, "argv", [str(path), "--dry-run"])
    monkeypatch.setattr(example, "REFINEMENT_CATEGORIES", None)
    example.main()
    output = capsys.readouterr().out
    assert "Active parameters: 18" in output
    assert "bounds=[600.0, 1200.0]" in output
    assert output.count("bounds=[0.5, 1.5]") == 11
    monkeypatch.setattr(example, "CATEGORY_BOUNDS", {})
    example.main()
    assert "bounds=unbounded" in capsys.readouterr().out
