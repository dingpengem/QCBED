from dataclasses import replace, FrozenInstanceError

import numpy as np
import pytest

from qcbed.paramin import _run, _fraction, _order, run_legacy_paramin
from qcbed.refinement import RefinementReflection
from qcbed.refinement_trial import evaluate_refinement_trial
from test_refinement import state, data


def engine(objective, initial=(10.,), bottom=.8, top=1.2, **kwargs):
    settings = dict(max_cycles=0, tolerance=0)
    settings.update(kwargs)
    n = len(initial)
    return _run(np.array(initial), np.full(n,bottom), np.full(n,top),
                lambda x:x.copy(), objective, **settings)


def test_exact_quadratic_samples_and_count():
    seen=[]
    def f(x):
        seen.append(x.copy())
        return (x[0]-11)**2
    r=engine(f)
    np.testing.assert_allclose(np.array(seen)[:,0],[10,8,12,11,11],atol=2e-14)
    assert r.cycles==1 and r.parameter_steps==1 and r.evaluations==5
    assert r.best_objective < 1e-26
    np.testing.assert_allclose(r.bottom_factors,[.9],atol=2e-15)
    assert not r.converged


def test_multiple_parameters_fixed_order():
    seen=[]
    def f(x):
        seen.append(x.copy())
        return (x[0]-11)**2+(x[1]-22)**2
    r=engine(f, initial=(10.,20.))
    expected=[[10,20],[8,20],[12,20],[11,20],[11,16],[11,24],[11,22],[11,22]]
    np.testing.assert_allclose(seen,expected,atol=1e-13)
    assert r.evaluations==8 and r.parameter_steps==2


@pytest.mark.parametrize("center,expected", [(30.,20.),(-10.,11.),(0.,11.)])
def test_minimum_cap_and_negative_or_zero_fallback(center,expected):
    seen=[]
    def f(x):
        seen.append(x[0])
        return (x[0]-center)**2
    engine(f)
    assert seen[3]==pytest.approx(expected,abs=1e-12)


def test_flat_ties_and_cycle_parity():
    r=engine(lambda x: 1.,max_cycles=1)
    # Flat interpolation uses alternating fallback, ties retain latest trial.
    np.testing.assert_allclose(r.best_vector,[10.45],atol=1e-13)
    assert r.evaluations==9 and r.cycles==2
    np.testing.assert_array_equal(r.uncertainties,[abs(1-.8)])


def test_degenerate_and_negative_curvature_formula():
    assert _fraction(np.float64(1),np.float64(1),np.float64(1),1.,1.,1,0)==1
    # Concave -(f-1.1)^2: AA=-1, BB=2.2 -> 2 + BB/(2AA)=.9.
    assert _fraction(np.float64(-.01),np.float64(-.09),np.float64(-.01),.8,1.2,1,0)==pytest.approx(.9)


def test_worse_parabolic_point_retains_sample_best():
    scores=iter([4.,1.,3.,5.,1.])
    r=engine(lambda x:next(scores))
    np.testing.assert_array_equal(r.best_vector,[8.])
    assert r.best_objective==1 and r.last_objective==1


def test_convergence_uses_old_range_and_initial_average():
    r=engine(lambda x:(x[0]-10)**2, tolerance=.15,max_cycles=20)
    assert r.cycles==2 and r.evaluations==9 and r.converged
    assert r.average_displacement==pytest.approx(.1)
    r=engine(lambda x:1.,tolerance=1,max_cycles=20)
    assert r.cycles==0 and r.evaluations==1 and r.converged


def test_reset_replaces_score_not_all_time_archive():
    scores=iter([4.,1.,3.,5.,7.])
    r=engine(lambda x:next(scores))
    assert r.best_objective==7 and r.last_objective==7
    np.testing.assert_array_equal(r.best_vector,[8.])


def test_literal_scalar_oracle_repeated_cycles():
    # Independent Lagrange divided differences for the parabola, scalar state
    # and explicit source evaluation/update sequence; no production helpers.
    x,b,t=10.,.8,1.2
    objective=lambda v:(v-11.3)**2+1
    best=objective(x)
    expected=[x]
    for cycle in range(1,5):
        start=x
        c0=best
        points=[b*start,t*(b*start)/b]
        c1,c2=map(objective,points)
        a=((c1-c0)/(b-1)-(c2-c0)/(t-1))/(b-t)
        bb=(c2-c0)/(t-1)-a*(t+1)
        frac=-bb/(2*a) if a>0 else 2+bb/(2*a)
        if not np.isfinite(frac) or frac==0 or frac<0:
            frac=1+(-1)**(cycle+1)*(1-b)/2
        frac=min(frac,2)
        points.append(frac*points[1]/t)
        for point in points:
            expected.append(point)
            score=objective(point)
            if score<=best:
                x,best=point,score
        change=abs(start-x)/abs(start)
        delta=min(change,.5)
        if change<abs(1-b)/2:
            delta=abs(1-b)/2
        b,t=1-delta,1+delta
        expected.append(x)
        best=objective(x)
    seen=[]
    def f(v):
        seen.append(v[0])
        return objective(v[0])
    r=engine(f,max_cycles=3)
    np.testing.assert_allclose(seen,expected,rtol=1e-12,atol=1e-12)
    assert r.evaluations==17 and r.cycles==4


def test_random_order_rejection_sampling_and_reproducibility():
    class Draws:
        def __init__(self): self.values=iter([.8,.8,.1,.5])
        def random(self): return next(self.values)
    assert _order(3,True,Draws())==[2,0,1]
    for seed in (0,7):
        a=_order(20,True,np.random.Generator(np.random.PCG64(seed)))
        b=_order(20,True,np.random.Generator(np.random.PCG64(seed)))
        assert a==b and sorted(a)==list(range(20))
    assert _order(20,True,np.random.default_rng(0))!=_order(20,True,np.random.default_rng(7))


def problem():
    s=state()
    return replace(s,regions=(replace(s.regions[0],thickness_mask=1),),ranges=((1.05,.95),))


def test_real_trial_and_immutability():
    s=problem()
    initial=np.array([2.])
    r=run_legacy_paramin(s,initial,data(),max_cycles=0,tolerance=0)
    assert r.best_objective==evaluate_refinement_trial(s,r.best_vector,data(),coordinate_mode="physical").objective
    assert initial[0]==2 and s.regions[0].depth_angstrom==2 and s.ranges==((1.05,.95),)
    assert r.evaluations==5
    with pytest.raises(FrozenInstanceError): r.cycles=0
    for a in (r.best_vector,r.last_vector,r.bottom_factors,r.top_factors,r.uncertainties):
        assert a.dtype==np.float64
        with pytest.raises(ValueError): a.setflags(write=True)


def test_capacity_missing_ranges_and_unsupported():
    large=replace(state(),reflections=(tuple(RefinementReflection(i,0,1.,1.,"YY") for i in range(501)),))
    for s,vector in [(large,np.zeros(1002)),(state(),np.empty(0))]:
        with pytest.raises(ValueError,match="1..1000"):
            run_legacy_paramin(s,vector,data(),max_cycles=0,tolerance=0)
    with pytest.raises(ValueError,match="range factors"):
        run_legacy_paramin(replace(problem(),ranges=()),np.array([2.]),data(),max_cycles=0,tolerance=0)
    s=problem()
    s=replace(s,regions=(replace(s.regions[0],vacancy_values=(.1,)),))
    with pytest.raises(ValueError,match="hazards"):
        run_legacy_paramin(s,np.array([2.]),data(),max_cycles=0,tolerance=0)


@pytest.mark.parametrize("kwargs", [dict(max_cycles=-1),dict(max_cycles=True),dict(tolerance=-1),
    dict(tolerance=np.nan),dict(random_order=1),dict(random_order=True),dict(seed=True)])
def test_invalid_controls(kwargs):
    settings=dict(max_cycles=0,tolerance=0)
    settings.update(kwargs)
    with pytest.raises((ValueError,TypeError)):
        run_legacy_paramin(problem(),np.array([2.]),data(),**settings)


@pytest.mark.parametrize("score", [np.nan,np.inf,1e21])
def test_invalid_objective(score):
    with pytest.raises(ValueError,match="objective"):
        engine(lambda x:score)


def test_zero_divisors_rejected():
    with pytest.raises(ValueError,match="zero"):
        engine(lambda x:1.,bottom=0.)
    with pytest.raises(ValueError,match="zero"):
        engine(lambda x:1.,initial=(0.,))


def test_projection_is_sequential_and_independent_adaptation_ifs():
    seen=[]
    def project(x):
        return np.maximum(x,9.)
    def evaluate(x):
        seen.append(x[0])
        return (x[0]-10)**2
    _run(np.array([10.]),np.array([.8]),np.array([1.2]),project,evaluate,
         max_cycles=0,tolerance=0)
    # First point 8 is projected to 9; next sample is 1.2*9/.8, not 12.
    np.testing.assert_allclose(seen[:3],[10,9,13.5],atol=1e-14)
    r=engine(lambda x:1.,bottom=-2.,top=4.)
    np.testing.assert_array_equal(r.bottom_factors,[-.5])
    np.testing.assert_array_equal(r.top_factors,[2.5])
