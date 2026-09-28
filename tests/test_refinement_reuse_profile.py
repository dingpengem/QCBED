"""Benchmark instrumentation equivalence and exact dependency census."""
from dataclasses import replace
import importlib.util
from pathlib import Path
import sys

import numpy as np
import pytest

from qcbed import ForwardExecutionConfig
import qcbed.refinement_trial as trial
from qcbed.refinement_driver import run_legacy_refinement
from qcbed.refinement_bounds import RefinementBound, RefinementBounds
from test_refinement import data
from test_refinement_driver import problem, config
from test_refinement_execution import assert_exact


@pytest.fixture
def profiler():
    path = Path(__file__).resolve().parents[1]/'benchmarks/profile_refinement_reuse.py'
    spec = importlib.util.spec_from_file_location('reuse_profile', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def evaluate(state, depth, **kwargs):
    return trial.evaluate_refinement_trial(state, np.array([depth], dtype=np.float64), data(),
        execution=ForwardExecutionConfig('batched', 3), **kwargs)


def test_identity_is_exact_and_census_uses_full_keys(profiler):
    assert profiler.exact(1.) != profiler.exact(np.nextafter(1., 2.))
    assert profiler.exact(0.) != profiler.exact(-0.)
    assert profiler.exact(np.array([1.])) != profiler.exact(np.array([[1.]]))
    census = profiler.Census()
    census.add({'a': profiler.exact(1.)})
    census.add({'a': profiler.exact(np.nextafter(1., 2.))})
    census.add({'a': profiler.exact(1.)})
    assert census.report(3)['a'] == dict(eligible=3, unique=2, repeats=1,
        potential_hit_rate=1/3, repeats_per_all_trials=1/3)


def test_trial_equivalence_timing_and_hook_restoration(profiler):
    state = problem()
    baseline = evaluate(state, 5.)
    original = trial._evaluate_with_invalid_policy
    seen, before = [], []
    with profiler.TrialProfiler() as profile:
        measured = evaluate(state, 5., evaluation_callback=seen.append,
                            before_evaluation=lambda: before.append(1))
    assert_exact(baseline, measured)
    assert seen == [measured] and before == [1]
    assert trial._evaluate_with_invalid_policy is original
    assert sys.gettrace() is None
    row = profile.rows[0]
    for name in ('mapping', 'sequence', 'sample_preparation', 'tilt_generation',
                 'propagator_construction', 'propagation', 'beam_extraction', 'intensity', 'objective'):
        assert row['seconds'][name] > 0
    assert row['other_seconds'] >= 0
    assert row['total_seconds'] >= row['seconds']['forward_inclusive']


def test_mapped_sequence_scan_and_forward_dependencies(profiler):
    state = problem()
    with profiler.TrialProfiler() as profile:
        evaluate(state, 2.)
        evaluate(state, -2.)  # raw sign differs, mapped state is identical
        evaluate(state, 1.9)  # same exact sequence, different mapped depth
        evaluate(state, np.nextafter(2., 3.))  # cumulative boundary crossed
        shifted = replace(state, scan=replace(state.scan, corner_a=(.01, 0.)))
        evaluate(shifted, 2.)
        thermal = replace(state, regions=(replace(state.regions[0], debye_waller=(.1,)),))
        evaluate(thermal, 2.)
    counts = profile.report()['census']
    assert counts['mapped']['repeats'] == 1
    assert counts['sequence']['repeats'] == 4
    assert counts['scan']['repeats'] == 4
    assert counts['prepared']['repeats'] == 4
    assert counts['propagator']['repeats'] == 4
    assert counts['forward']['repeats'] == 2
    assert not profile.rows[4]['repeated']['forward']
    assert not profile.rows[5]['repeated']['prepared']
    assert profile.rows[5]['repeated']['propagator']


@pytest.mark.parametrize('method', ['ANSI', 'SPAR', 'PARA'])
def test_driver_decisions_counts_and_callbacks_unchanged(profiler, method):
    cfg = replace(config(method), execution=ForwardExecutionConfig('batched', 3))
    def run(events, checks):
        return run_legacy_refinement(problem(), data(), cfg,
            progress_callback=lambda *event: events.append(event),
            before_evaluation=lambda: checks.append(1))
    events, checks = [], []
    baseline = run(events, checks)
    actual_events, actual_checks = [], []
    with profiler.TrialProfiler() as profile:
        measured = run(actual_events, actual_checks)
    assert_exact(baseline, measured)
    assert_exact(events, actual_events)
    assert checks == actual_checks
    assert measured.evaluations == len(profile.rows) == len(checks)


def test_rejections_are_not_forward_hits_and_failure_restores_hooks(profiler):
    original = trial._evaluate_with_invalid_policy
    bounds = RefinementBounds((RefinementBound('depth', (0,), 1., 3.),))
    with profiler.TrialProfiler() as profile:
        result = evaluate(problem(), 5., bounds_policy='reject', refinement_bounds=bounds)
        assert isinstance(result, trial.BoundRejectedTrial)
    assert 'forward' not in profile.report()['census']
    assert profile.rows[0]['occurrences'] is None
    with pytest.raises(ValueError), profiler.TrialProfiler():
        evaluate(problem(), 0.)
    assert trial._evaluate_with_invalid_policy is original
    assert sys.gettrace() is None
