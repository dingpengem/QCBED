from dataclasses import replace, FrozenInstanceError

import numpy as np
import pytest

import qcbed.simpara as sim
from qcbed.refinement import RefinementReflection
from qcbed.refinement_trial import evaluate_refinement_trial
from test_refinement import state, data


def engine(objective, initial=(10.,20.), uncertainty=(1.,2.), **kwargs):
    settings=dict(max_evaluations=1,high_temperature=1.,low_temperature=.01,
                  cooling=.5,random_order=False,rng=np.random.default_rng(2))
    settings.update(kwargs)
    return sim._run(np.array(initial),np.array(uncertainty),objective,**settings)


class Uniforms:
    def __init__(self): self.i=0
    def random(self):
        value=(.75,.5)[self.i%2]
        self.i+=1
        return np.float64(value)


def test_exact_first_proposals_controlled_rng():
    seen=[]
    def f(x):
        seen.append(x.copy())
        return 1.
    rng=Uniforms()
    r=engine(f,rng=rng)
    np.testing.assert_allclose(seen,[[10,20],[8.334890777684604,16.669781555369208],
                                    [11.665109222315396,23.330218444630792]],atol=1e-14)
    assert rng.i==4 and r.evaluations==3 and r.cycles==1
    # Driver callback retains ties even though SIMPARA's strict update does not.
    np.testing.assert_array_equal(r.best_vector,seen[-1])
    assert r.improvements==0 and r.temperature==.5


def test_improvement_heating_and_retained_vs_last(monkeypatch):
    monkeypatch.setattr(sim,"_normal",lambda rng,center,width:width)
    scores=iter([10.,20.,8.,9.])
    r=engine(lambda x:next(scores))
    np.testing.assert_array_equal(r.best_vector,[11,22])
    np.testing.assert_allclose(r.last_vector,[10.75,21.5],atol=1e-14)
    assert r.best_objective==8 and r.last_objective==9
    assert r.improvements==1 and r.temperature==2
    np.testing.assert_array_equal(r.perturbation_scales,[2,4])
    assert r.evaluations==4


def test_worse_rejected_and_cooling(monkeypatch):
    monkeypatch.setattr(sim,"_normal",lambda rng,center,width:width)
    r=engine(lambda x:1. if np.array_equal(x,[10.,20.]) else 2.)
    np.testing.assert_array_equal(r.best_vector,[10,20])
    assert r.best_objective==1 and r.temperature==.5 and r.improvements==0


def test_single_parameter_reset_forces_cooling_and_rechecks_score(monkeypatch):
    monkeypatch.setattr(sim,"_normal",lambda rng,center,width:width)
    scores=iter([5.,8.,4.,3.,7.])
    r=engine(lambda x:next(scores),initial=(10.,),uncertainty=(1.,))
    assert r.best_objective==7 and r.last_objective==7
    assert r.improvements==1 and r.temperature==.5 and r.evaluations==5
    np.testing.assert_array_equal(r.direction,[0])
    np.testing.assert_array_equal(r.central_vector,r.best_vector)


def test_threshold_inclusive_and_initial_count(monkeypatch):
    monkeypatch.setattr(sim,"_normal",lambda rng,center,width:width)
    r=engine(lambda x:1.,high_temperature=.5,low_temperature=1.,max_evaluations=100)
    assert r.cycles==1 and r.stop_reason=="temperature"
    r=engine(lambda x:1.,max_evaluations=0)
    assert r.evaluations==1 and r.cycles==0 and r.stop_reason=="max_evaluations"
    r=engine(lambda x:1.,high_temperature=.4,low_temperature=1.)
    assert r.cycles==0 and r.evaluations==1 and r.stop_reason=="temperature"


def test_zero_negative_coordinates_and_zero_width(monkeypatch):
    monkeypatch.setattr(sim,"_normal",lambda rng,center,width:width)
    seen=[]
    def f(x):
        seen.append(x.copy())
        return 1.
    r=engine(f,initial=(0.,-2.),uncertainty=(0.,1.))
    np.testing.assert_array_equal(seen,[[0,-2],[0,-3],[0,-1]])
    assert r.evaluations==3


def test_parabola_zero_rejection_and_branch_specific_cap():
    c=np.array([0.,20000.,-20000.])
    # all share a symmetric bowl: exact minima central; zero does not qualify.
    new,qualified,counter=sim._parabola(c,c-1,c+1,np.float64(0),np.float64(1),np.float64(1),3)
    np.testing.assert_array_equal(new,[.5,10000,-20000])
    assert qualified==2 and counter==6
    new,qualified,_=sim._parabola(c,c-1,c+1,np.float64(1),np.float64(1),np.float64(1),3)
    assert qualified==0 and new[1]>10000  # cap is not a global bound


def test_literal_single_coordinate_control_oracle(monkeypatch):
    # Independent scalar, symmetric divided differences; n=1 resets every cycle.
    normals=iter([.5,-1.,.75])
    monkeypatch.setattr(sim,"_normal",lambda rng,center,width:width*next(normals))
    expected=[]
    out=4.
    temp=1.
    f=lambda x:(x-5.2)**2+1
    score=f(out)
    expected.append(out)
    for z in (.5,-1.,.75):
        center=out
        d=temp*z
        inverse,trial=center-d,center+d
        y1,y2=f(inverse),f(trial)
        expected.extend([inverse,trial])
        a=(y1+y2-2*score)/(2*d*d)
        if a>0:
            slope=(y2-y1)/(2*d)
            new=center-slope/(2*a)
            expected.append(new)
            y3=f(new)
        else:
            new=center+d/2
            y3=1e20
        # Both callback <= writes and routine strict comparisons.
        chi=score
        for x,y in ((inverse,y1),(trial,y2),(new,y3)):
            if y<=chi:
                out,chi=x,y
        for x,y in ((inverse,y1),(trial,y2),(new,y3)):
            if y<score:
                out,score=x,y
        expected.append(out)
        score=f(out)
        temp*=.5
    seen=[]
    def evaluate(x):
        seen.append(x[0])
        return f(x[0])
    r=engine(evaluate,initial=(4.,),uncertainty=(1.,),max_evaluations=9)
    np.testing.assert_allclose(seen,expected,atol=1e-12,rtol=1e-12)
    assert r.evaluations==13 and r.cycles==3 and r.temperature==.125


def test_direction_rejection_and_parallel_branch(monkeypatch):
    # N=5 avoids a reset. First improvement sets direction; cycle 2 fails,
    # cycle 3 rejects the all-positive sign vector (dot=5), then accepts dot=1;
    # cycle 4 follows the even no-improvement parallel branch.
    signs=iter([1]*15+[1,1,1,-1,-1]+[1]*5)
    draws=[]
    def normal(rng,center,width):
        value=next(signs)*width
        draws.append(value)
        return value
    monkeypatch.setattr(sim,"_normal",normal)
    scores=iter([10,14,8,9]+[9,10,11]*3)
    seen=[]
    def evaluate(x):
        seen.append(x.copy())
        return next(scores)
    r=engine(evaluate,initial=(10.,)*5,uncertainty=(1.,)*5,max_evaluations=10)
    assert r.cycles==4 and len(draws)==25 and r.evaluations==13
    np.testing.assert_array_equal(r.direction,np.ones(5))
    np.testing.assert_array_equal(seen[7],[10,10,10,12,12])


def problem():
    s=state()
    return replace(s,scan_masks=("YN","NN","NN"))


def run(s,initial,seed=4,**kwargs):
    settings=dict(uncertainties=np.array([.1]),max_evaluations=4,high_temperature=1.,
                  low_temperature=.01,cooling=.5,seed=seed)
    settings.update(kwargs)
    return sim.run_legacy_simpara(s,initial,data(),**settings)


def test_real_trial_reproducibility_and_immutability():
    s=problem()
    initial=np.array([0.])
    a=run(s,initial)
    b=run(s,initial)
    np.testing.assert_array_equal(a.best_vector,b.best_vector)
    np.testing.assert_array_equal(a.last_vector,b.last_vector)
    assert a.best_objective==evaluate_refinement_trial(s,a.best_vector,data()).objective
    assert initial[0]==0 and s.scan.corner_a==(0.,0.)
    for array in (a.best_vector,a.last_vector,a.central_vector,a.direction,a.perturbation_scales):
        assert array.dtype==np.float64
        with pytest.raises(ValueError): array.setflags(write=True)
    with pytest.raises(FrozenInstanceError): a.cycles=0
    # Flat objective isolates random trajectory from acceptance.
    first=engine(lambda x:1.,rng=np.random.default_rng(1))
    other=engine(lambda x:1.,rng=np.random.default_rng(2))
    assert not np.array_equal(first.last_vector,other.last_vector)


def test_capacity_and_mapping_failure():
    large=replace(state(),reflections=(tuple(RefinementReflection(i,0,1.,1.,"YY") for i in range(501)),))
    with pytest.raises(ValueError,match="1..1000"):
        run(large,np.zeros(1002))
    with pytest.raises(ValueError,match="1..1000"):
        run(state(),np.empty(0))
    s=problem()
    s=replace(s,regions=(replace(s.regions[0],vacancy_values=(.1,)),))
    with pytest.raises(ValueError,match="hazards"):
        run(s,np.array([0.]))


def test_random_order_assigns_draws_in_source_order(monkeypatch):
    class Draws:
        def __init__(self): self.values=iter([.8,.8,.1,.5])
        def random(self): return next(self.values)
    normals=iter([1.,2.,3.])
    monkeypatch.setattr(sim,"_normal",lambda rng,center,width:next(normals))
    seen=[]
    def evaluate(x):
        seen.append(x.copy())
        return 1.
    engine(evaluate,initial=(10.,20.,30.),uncertainty=(1.,1.,1.),
           random_order=True,rng=Draws())
    np.testing.assert_array_equal(seen,[[10,20,30],[8,17,29],[12,23,31]])


def test_real_callback_parameter_writeback_without_zero_floor():
    s=state()
    s=replace(s,regions=(replace(s.regions[0],debye_waller_masks=(1,)),))
    initial=np.array([-.3])
    r=run(s,initial,uncertainties=np.array([0.]),max_evaluations=0)
    np.testing.assert_array_equal(r.best_vector,[.3])
    assert initial[0]==-.3 and s.regions[0].debye_waller==(0.,)
    r=run(problem(),np.array([0.]),max_evaluations=0)
    assert r.best_vector[0]==0  # no ARRANGEPARS projection inside SIMPARA


@pytest.mark.parametrize("kwargs",[dict(seed=True),dict(seed=-1),dict(max_evaluations=-1),
    dict(cooling=1),dict(low_temperature=0),dict(high_temperature=np.inf),
    dict(uncertainties=np.array([-.1])),dict(random_order=1)])
def test_controls(kwargs):
    with pytest.raises((TypeError,ValueError)):
        run(problem(),np.array([0.]),**kwargs)


@pytest.mark.parametrize("score",[np.nan,np.inf,1e21])
def test_invalid_objectives(score):
    with pytest.raises(ValueError,match="objective"):
        engine(lambda x:score)
