from dataclasses import replace, FrozenInstanceError

import numpy as np
import pytest

from qcbed.anneal import _normal, _run, run_legacy_anneal
from qcbed.refinement import RefinementReflection, pack_refinement_parameters
from qcbed.refinement_trial import evaluate_refinement_trial
from test_refinement import state, data


class Uniforms:
    def __init__(self, pair):
        self.pair, self.i = pair, 0

    def random(self):
        result = self.pair[self.i % len(self.pair)]
        self.i += 1
        return np.float64(result)


def engine(objective, *, initial=10., maximum=.1, minimum=.0001, cooling=.5, limit=0, pair=(.75, .5)):
    return _run(np.array([initial]), lambda v: v.copy(), objective, max_iterations=limit,
        maximum_temperature=maximum, minimum_temperature=minimum, cooling_rate=cooling, rng=Uniforms(pair))


def test_exact_first_normal_signed_width_and_rejection():
    rng = Uniforms((1., 1., .5, .5, .75, .5))  # outside disk, zero radius, valid
    assert _normal(rng, np.float64(10), np.float64(1)) == pytest.approx(11.665109222315396, abs=2e-15)
    assert rng.i == 6
    assert _normal(Uniforms((.75, .5)), np.float64(-10), np.float64(-1)) == pytest.approx(-11.665109222315396)


def test_improvement_tie_and_failure():
    improving = engine(lambda v: -v[0])
    expected = 11.665109222315396
    assert improving.best_vector[0] == pytest.approx(expected)
    assert improving.best_objective == improving.last_objective
    assert improving.amplitudes[0] == pytest.approx(1.48258*(expected-10)/10)
    tied = engine(lambda v: 1.)
    np.testing.assert_array_equal(tied.best_vector, improving.best_vector)
    failing = engine(lambda v: (v[0]-10)**2)
    assert failing.best_vector[0] == 10 and failing.best_objective == 0
    assert failing.last_vector[0] > 10 and failing.last_objective > 0
    assert failing.amplitudes[0] == .05
    assert failing.temperature == .1
    assert failing.iterations == 1 and failing.evaluations == 2  # MAXITS=0


def test_count_before_cooling_reset_and_temperature_stop():
    r = engine(lambda v: (v[0]-10)**2, maximum=.2, minimum=.1, limit=100)
    assert r.iterations == 3 and r.evaluations == 5
    assert r.temperature == .1 and r.amplitudes[0] == .1
    assert r.last_vector[0] == 10 and r.last_objective == 0


def test_periodic_and_temperature_resets_both_evaluate():
    r = engine(lambda v: 1., maximum=1., minimum=1e-8, cooling=.9, limit=9, pair=(.5, .75))
    assert r.iterations == 10 and r.evaluations == 22  # 1+10+10+1
    assert r.temperature == pytest.approx(.9**10)
    assert r.best_vector[0] == 10


def test_initial_temperature_threshold_and_empty_vector():
    r = engine(lambda v: 2., maximum=.1, minimum=.1, limit=10)
    assert r.iterations == 0 and r.evaluations == 1
    r = _run(np.empty(0), lambda v: v.copy(), lambda v: 3., max_iterations=2,
             maximum_temperature=1., minimum_temperature=.01, cooling_rate=.5, rng=Uniforms((.75, .5)))
    assert r.iterations == 3 and r.evaluations == 7
    assert r.best_vector.shape == r.amplitudes.shape == (0,)


def oracle(seed):
    # Independent flat scalar transcription, no production RNG/control helpers.
    rng = np.random.Generator(np.random.PCG64(seed))
    retained = -2.
    chi = (retained-1)**2
    threshold = chi
    amplitude, temperature = .4, .4
    evaluations, steps = 1, 0
    last, fit = retained, chi
    while steps <= 12 and temperature > .1:
        steps += 1
        while True:
            u, v = 2*rng.random()-1, 2*rng.random()-1
            radius = u*u+v*v
            if 0 < radius < 1:
                break
        last = (retained*amplitude)*np.sqrt(-2*np.log(radius)/radius)*u+retained
        candidate_width = 1.48258*abs((last-retained)/retained)
        fit = (last-1)**2
        evaluations += 1
        if fit <= chi:
            retained, chi = last, fit
        changed = fit <= threshold
        if changed:
            threshold, amplitude = fit, candidate_width
        small = amplitude < .1
        if not changed:
            amplitude *= .5
        if small:
            temperature *= .5
            amplitude = temperature
            last = retained
            fit = chi = threshold = (retained-1)**2
            evaluations += 1
        if steps % 10 == 0:
            last = retained
            fit = chi = threshold = (retained-1)**2
            evaluations += 1
    return retained, chi, last, fit, threshold, steps, evaluations, temperature, amplitude


@pytest.mark.parametrize("seed", [1, 8, 30])
def test_literal_source_oracle(seed):
    r = _run(np.array([-2.]), lambda v: v.copy(), lambda v: (v[0]-1)**2,
        max_iterations=12, maximum_temperature=.4, minimum_temperature=.1,
        cooling_rate=.5, rng=np.random.Generator(np.random.PCG64(seed)))
    actual = (r.best_vector[0], r.best_objective, r.last_vector[0], r.last_objective,
              r.adaptation_objective, r.iterations, r.evaluations, r.temperature, r.amplitudes[0])
    np.testing.assert_array_equal(actual, oracle(seed))


def problem():
    s = state()
    return replace(s, scan_masks=("YN", "NN", "NN"))


def run(s, vector, seed=4, **kwargs):
    settings = dict(max_iterations=2, seed=seed, maximum_temperature=.1,
                    minimum_temperature=.001, cooling_rate=.5)
    settings.update(kwargs)
    return run_legacy_anneal(s, vector, data(), **settings)


def test_real_trial_determinism_zero_input_and_immutability():
    s = problem()
    vector = np.array([0.])
    first = run(s, vector)
    second = run(s, vector)
    for name in ("best_vector", "last_vector", "amplitudes"):
        np.testing.assert_array_equal(getattr(first, name), getattr(second, name))
        assert getattr(first, name).dtype == np.float64
        with pytest.raises(ValueError):
            getattr(first, name).setflags(write=True)
    assert first.best_objective == evaluate_refinement_trial(s, first.best_vector, data()).objective
    assert first.best_vector[0] != 0  # ARRANGEPARS projects zero to +1e-10
    assert vector[0] == 0 and s.scan.corner_a == (0., 0.)
    with pytest.raises(FrozenInstanceError):
        first.iterations = 9


def test_seed_changes_path_and_no_global_rng_use():
    s = problem()
    first = run(s, np.array([-.2]), seed=1, max_iterations=0)
    np.random.seed(999)
    repeated = run(s, np.array([-.2]), seed=1, max_iterations=0)
    other = run(s, np.array([-.2]), seed=2, max_iterations=0)
    np.testing.assert_array_equal(first.last_vector, repeated.last_vector)
    assert not np.array_equal(first.last_vector, other.last_vector)


def test_capacity_and_mapping_failure_propagate():
    s = state()
    large = replace(s, reflections=(tuple(RefinementReflection(i, 0, 1., 1., "YY") for i in range(501)),))
    with pytest.raises(ValueError, match="1000"):
        run(large, pack_refinement_parameters(large))
    bad = replace(s, regions=(replace(s.regions[0], vacancy_values=(.1,)),))
    with pytest.raises(ValueError, match="hazards"):
        run(bad, np.empty(0))
    assert bad.regions[0].vacancy_values == (.1,)


@pytest.mark.parametrize("kwargs", [dict(max_iterations=-1), dict(max_iterations=True), dict(seed=-1),
    dict(seed=True), dict(cooling_rate=1), dict(cooling_rate=0), dict(minimum_temperature=0),
    dict(maximum_temperature=np.inf)])
def test_invalid_controls(kwargs):
    with pytest.raises((ValueError, TypeError)):
        run(problem(), np.array([1.]), **kwargs)


@pytest.mark.parametrize("value", [np.nan, np.inf, 1e21])
def test_undefined_objective_is_failure_not_penalty(value):
    with pytest.raises(ValueError, match="objective"):
        engine(lambda v: value)
