from dataclasses import replace, FrozenInstanceError

import numpy as np
import pytest

from qcbed.amoeba import build_legacy_simplex, run_legacy_amoeba, _run, _rank
from qcbed.refinement import pack_refinement_parameters, RefinementReflection
from qcbed.refinement_trial import evaluate_refinement_trial
from test_refinement import state, data


def problem(n=1):
    s = state()
    s = replace(s, regions=(replace(s.regions[0], thickness_mask=1),), ranges=((.5, 1.5),))
    if n == 2:
        s = replace(s, scan_masks=("YN", "NN", "NN"), ranges=((.5, 1.5), (.5, 1.5)))
    return s


def test_initialization_and_scal_anchor():
    s = problem()
    initial = np.array([3.])
    p = build_legacy_simplex(s, initial, np.array([-.25]))
    np.testing.assert_array_equal(p, [[1.], [.75]])  # SCAL=(.5,1), anchored to depth 2
    assert initial[0] == 3 and s.regions[0].depth_angstrom == 2
    s = problem(2)
    p = build_legacy_simplex(s, pack_refinement_parameters(s), np.array([.2, .3]))
    np.testing.assert_allclose(p, [[.5,.5],[.7,.5],[.5,.8]], atol=1e-15)
    # Packed zero scan coordinate is +1e-10; independent displacement stays finite.
    assert p.dtype == np.float64 and p.shape == (3,2)
    with pytest.raises(ValueError):
        p.setflags(write=True)
    np.testing.assert_array_equal(build_legacy_simplex(problem(), np.array([2.]), np.array([0.])), [[.5],[.5]])


@pytest.mark.parametrize("values,points,final", [
    ([0,2,4,1], [[0,0],[1,0],[0,1],[1,-1]], [[0,0],[1,0],[1,-1]]),
    ([0,2,4,-2,-1], [[0,0],[1,0],[0,1],[1,-1],[1.5,-2]], [[0,0],[1,0],[1.5,-2]]),
    ([0,2,4,-2,0], [[0,0],[1,0],[0,1],[1,-1],[1.5,-2]], [[0,0],[1,0],[1,-1]]),
    ([0,2,4,3,2.5], [[0,0],[1,0],[0,1],[1,-1],[.75,-.5]], [[0,0],[1,0],[.75,-.5]]),
    ([0,2,4,5,3], [[0,0],[1,0],[0,1],[1,-1],[.25,.5]], [[0,0],[1,0],[.25,.5]]),
    ([0,2,4,5,4,1,2], [[0,0],[1,0],[0,1],[1,-1],[.25,.5],[.5,0],[0,.5]], [[0,0],[.5,0],[0,.5]]),
    ([0,2,4,3,3,1,2], [[0,0],[1,0],[0,1],[1,-1],[.75,-.5],[.5,0],[.5,-.5]], [[0,0],[.5,0],[.5,-.5]]),
])
def test_source_branches_exact_trials(values, points, final):
    seen = []
    def objective(p):
        seen.append(p.copy())
        return values[len(seen)-1]
    initial = np.array([[0.,0.],[1.,0.],[0.,1.]])
    r = _run(initial, objective, ftol=0, max_evaluations=0, fmin=-100)
    np.testing.assert_array_equal(seen, points)
    np.testing.assert_array_equal(r.final_simplex, final)
    np.testing.assert_array_equal(initial, [[0,0],[1,0],[0,1]])
    assert r.evaluations == len(values) and r.search_evaluations == len(values)-3
    assert r.iterations == 1 and not r.converged
    assert r.best_evaluated_objective == min(values)
    if values == [0,2,4,-2,-1]:
        assert r.best_objective == -1 and r.best_evaluated_objective == -2
        np.testing.assert_array_equal(r.best_evaluated_vector, [1,-1])


def test_ranking_ties_and_no_best_swap():
    assert _rank([1,1,1]) == (0,1,0)
    assert _rank([3,0,0]) == (1,0,1)
    r = _run(np.array([[2.],[1.]]), lambda x:x[0]**2, ftol=0, max_evaluations=0, fmin=3)
    assert r.best_index == 1 and r.best_vector[0] == 1 and r.iterations == 0
    assert r.stop_reason == "objective_threshold" and r.evaluations == 2


def test_spread_strictness_zero_denominator_and_tied_expansion():
    p = np.array([[0.],[1.]])
    r = _run(p, lambda x: 2., ftol=.01, max_evaluations=0, fmin=-1)
    assert r.stop_reason == "objective_spread" and r.evaluations == 2
    r = _run(p, lambda x: 0., ftol=.01, max_evaluations=0, fmin=0)
    assert r.stop_reason == "objective_threshold"  # 0/0 not made into zero RTOL
    r = _run(p, lambda x: 2., ftol=0, max_evaluations=0, fmin=-1)
    assert r.evaluations == 4  # reflection tie triggers expansion, expansion tie rejected
    np.testing.assert_array_equal(r.final_simplex, [[0],[-1]])


def literal_one_dimensional_oracle(start, limit):
    # Independent scalar form of simpara.o; no production optimizer helpers.
    x = list(start)
    f = lambda t:(t-.3)**2+1
    y = list(map(f,x))
    calls = cycles = 0
    while True:
        low = 0 if y[0] <= y[1] else 1
        high = 1-low
        if 2*abs(y[high]-y[low])/(abs(y[high])+abs(y[low])) < 1e-10 or calls > limit:
            return x,y,calls,cycles
        reflected = 2*x[low]-x[high]
        fr = f(reflected)
        calls += 1
        if fr <= y[low]:
            expanded = 2*reflected-x[low]
            fe = f(expanded)
            calls += 1
            x[high],y[high] = (expanded,fe) if fe < y[low] else (reflected,fr)
        else:
            if fr < y[high]:
                x[high],y[high] = reflected,fr
            contracted = .5*x[high]+.5*x[low]
            fc = f(contracted)
            calls += 1
            if fc < y[high]:
                x[high],y[high] = contracted,fc
            else:
                x[high] = .5*(x[high]+x[low])
                y[high] = f(x[high])
                calls += 1
        cycles += 1


@pytest.mark.parametrize("start", [(4.,6.),(-3.,1.),(.3,.7)])
def test_literal_source_oracle(start):
    expected = literal_one_dimensional_oracle(start, 60)
    r = _run(np.array(start)[:,None], lambda x:(x[0]-.3)**2+1,
             ftol=1e-10, max_evaluations=60, fmin=0)
    np.testing.assert_array_equal(r.final_simplex[:,0], expected[0])
    np.testing.assert_array_equal(r.final_objectives, expected[1])
    assert (r.search_evaluations,r.iterations) == expected[2:]


def test_real_trial_coordinate_mode_reproducibility_and_immutability():
    s = problem()
    p = build_legacy_simplex(s, np.array([2.]), np.array([.1]))
    settings = dict(ftol=1e-8,max_evaluations=3,fmin=-1)
    r = run_legacy_amoeba(s,p,data(),**settings)
    again = run_legacy_amoeba(s,p,data(),**settings)
    np.testing.assert_array_equal(r.final_simplex,again.final_simplex)
    for row,y in zip(r.final_simplex,r.final_objectives):
        assert y == evaluate_refinement_trial(s,row.copy(),data(),coordinate_mode="legacy_simplex").objective
    assert s.regions[0].depth_angstrom == 2
    np.testing.assert_array_equal(p, [[.5],[.6]])
    with pytest.raises(FrozenInstanceError):
        r.iterations = 0
    with pytest.raises(ValueError):
        r.best_vector.setflags(write=True)


@pytest.mark.parametrize("p", [np.zeros((1,1)), np.zeros((2,2)),np.array([[0.],[np.inf]])])
def test_shape_validation(p):
    with pytest.raises(ValueError,match="shape"):
        run_legacy_amoeba(problem(),p,data(),ftol=.01,max_evaluations=0,fmin=0)


def test_capacity_and_unsupported_mapping():
    with pytest.raises(ValueError,match="1..1000"):
        build_legacy_simplex(state(),np.empty(0),np.empty(0))
    large = replace(state(), reflections=(tuple(RefinementReflection(i,0,1.,1.,"YY") for i in range(501)),))
    with pytest.raises(ValueError,match="1..1000"):
        build_legacy_simplex(large,np.zeros(1002),np.zeros(1002))
    s = problem()
    s = replace(s,regions=(replace(s.regions[0],vacancy_values=(.1,)),))
    with pytest.raises(ValueError,match="hazards"):
        run_legacy_amoeba(s,np.array([[.5],[.6]]),data(),ftol=.01,max_evaluations=0,fmin=0)


@pytest.mark.parametrize("value", [np.nan,np.inf,-np.inf])
def test_nonfinite_objective_rejected(value):
    with pytest.raises(ValueError,match="finite"):
        _run(np.array([[0.],[1.]]),lambda x:value,ftol=0,max_evaluations=0,fmin=0)
