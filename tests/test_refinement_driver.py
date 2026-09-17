from dataclasses import replace, FrozenInstanceError
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

import qcbed.refinement_driver as driver
from qcbed.refinement import pack_refinement_parameters, build_refinement_parameter_layout, apply_refinement_parameters
from qcbed.refinement_trial import evaluate_refinement_trial
from test_refinement import state, data


MODES={"ANNE":("ANNEAL",),"PARA":("PARAMIN",),"SPAR":("SIMPARA",),"SIMP":("AMOEBA",),
       "PAAN":("PARAMIN","ANNEAL"),"ANPA":("ANNEAL","PARAMIN"),
       "PASI":("PARAMIN","AMOEBA"),"SIPA":("AMOEBA","PARAMIN"),
       "SIAN":("AMOEBA","ANNEAL"),"ANSI":("ANNEAL","AMOEBA")}


def problem():
    s=state()
    return replace(s,regions=(replace(s.regions[0],thickness_mask=1),),ranges=((1.5,.5),))


def config(method="ANNE",loops=1,**kwargs):
    values=dict(method=method,loops=loops,maxfsimp=2,maxfpara=3,ftol=.01,chimin=.001,
                maxtemp=1.,mintemp=.1,coolrate=.5,plamda=(.2,),seed=42)
    values.update(kwargs)
    return driver.LegacyRefinementConfig(**values)


def install_stubs(monkeypatch):
    calls=[]
    pres=[]
    def evaluate(s,v,e,**kw):
        assert kw["coordinate_mode"]=="physical"
        pres.append(v.copy())
        return SimpleNamespace(objective=(20-v[0])**2,physical_state=apply_refinement_parameters(s,v))
    monkeypatch.setattr(driver,"evaluate_refinement_trial",evaluate)
    for name in ("ANNEAL","PARAMIN","SIMPARA","AMOEBA"):
        def optimize(s,v,e,_name=name,**kw):
            calls.append((_name,v.copy(),kw))
            if _name=="AMOEBA":
                # Fixed baseline SCAL: a=-.5, b=3. Improve physical x by one.
                retained=np.array([v[0,0]-.5])
                x=retained[0]/(-.5)+3
                return SimpleNamespace(best_evaluated_vector=retained,best_evaluated_objective=(20-x)**2,
                                       best_vector=np.array([99.]),evaluations=7)
            out=v+1
            return SimpleNamespace(best_vector=out,best_objective=(20-out[0])**2,
                                   last_vector=v+99,last_objective=999.,evaluations=5)
        monkeypatch.setattr(driver,"run_legacy_"+name.lower(),optimize)
    return calls,pres


@pytest.mark.parametrize("mode",list(MODES))
def test_all_source_sequences_and_handoffs(mode,monkeypatch):
    calls,pres=install_stubs(monkeypatch)
    s=problem()
    initial=pack_refinement_parameters(s)
    r=driver.run_legacy_refinement(s,data(),config(mode,loops=2))
    order=MODES[mode]*2
    assert tuple(c[0] for c in calls)==order
    assert tuple(h.optimizer for h in r.stage_history)==order
    assert [h.loop_index for h in r.stage_history]==[1]*len(MODES[mode])+[2]*len(MODES[mode])
    assert len(pres)==len(order)  # no extra initial/final reevaluation
    np.testing.assert_array_equal(np.array(pres)[:,0],np.arange(2,2+len(order)))
    for i,(name,v,kw) in enumerate(calls):
        if name=="AMOEBA":
            q=(2+i-3)*(-.5)
            np.testing.assert_allclose(v,[[q],[q+.2]],atol=1e-15)
            assert kw["ftol"]==.01 and kw["fmin"]==.001 and kw["max_evaluations"]==2
        elif name=="ANNEAL":
            assert kw["max_iterations"]==2 and kw["maximum_temperature"]==1
        elif name=="PARAMIN":
            assert kw["max_cycles"]==3 and kw["tolerance"]==.01
        else:
            np.testing.assert_array_equal(kw["uncertainties"],[.2])
            assert kw["max_evaluations"]==2 and kw["cooling"]==.5
        h=r.stage_history[i]
        assert h.input_objective==(20-(2+i))**2
        assert h.retained_objective==(20-(3+i))**2
    assert r.final_vector[0]==2+len(order) and r.final_objective==(18-len(order))**2
    assert r.evaluations==sum(1+(7 if x=="AMOEBA" else 5) for x in order)
    assert r.initial_objective==324 and s.regions[0].depth_angstrom==2
    np.testing.assert_array_equal(initial,[2])
    with pytest.raises(ValueError): r.final_vector.setflags(write=True)
    with pytest.raises(FrozenInstanceError): r.loop_count=10


@pytest.mark.parametrize("loops",[0,-2])
def test_nonpositive_loops_one_evaluation(loops,monkeypatch):
    calls,pres=install_stubs(monkeypatch)
    r=driver.run_legacy_refinement(problem(),data(),config("unknown",loops))
    assert not calls and len(pres)==1 and r.evaluations==1 and r.loop_count==0
    assert r.stage_history==() and r.final_vector[0]==2


def test_amoeba_preserves_pre_stage_callback_best(monkeypatch):
    install_stubs(monkeypatch)
    monkeypatch.setattr(driver,"run_legacy_amoeba",lambda *a,**k:SimpleNamespace(
        best_evaluated_vector=np.array([.6]),best_evaluated_objective=999.,evaluations=2))
    r=driver.run_legacy_refinement(problem(),data(),config("SIMP"))
    assert r.final_vector[0]==2 and r.final_objective==324 and r.evaluations==3


def test_raw_retained_point_vs_next_packing(monkeypatch):
    s=replace(state(),scan_masks=("YN","NN","NN"))
    seen=[]
    def evaluate(s,v,e,**kw):
        seen.append(v.copy())
        return SimpleNamespace(objective=1.,physical_state=apply_refinement_parameters(s,v))
    monkeypatch.setattr(driver,"evaluate_refinement_trial",evaluate)
    monkeypatch.setattr(driver,"run_legacy_anneal",lambda *a,**k:SimpleNamespace(
        best_vector=np.array([0.]),best_objective=1.,evaluations=1))
    r=driver.run_legacy_refinement(s,data(),config(loops=2))
    assert r.final_vector[0]==0  # PARMIN
    assert seen[1][0]==1e-10  # ARRANGEPARS puts floor into P only


def test_no_invented_all_time_archive(monkeypatch):
    install_stubs(monkeypatch)
    monkeypatch.setattr(driver,"run_legacy_anneal",lambda *a,**k:SimpleNamespace(
        best_vector=np.array([1.]),best_objective=361.,evaluations=1))
    r=driver.run_legacy_refinement(problem(),data(),config("PAAN"))
    assert r.stage_history[0].retained_objective==289
    assert r.final_vector[0]==1 and r.final_objective==361


def test_parsed_invalid_selector_has_context(tmp_path):
    source=Path(__file__).parents[1]/"reference_data"/"FOR007.DAT"
    text=source.read_text()
    text=text.replace("SPAR","SpAr").replace("spar","SpAr")
    path=tmp_path/"bad.dat"
    path.write_text(text)
    with pytest.raises(ValueError,match="bad.dat: section 'refinement controls', record"):
        driver.parse_legacy_refinement(path,seed=0,experimental_header_path=source.parent/"in.dat")


def test_real_hybrid_and_seed_determinism():
    s=problem()
    c=config("ANSI",maxfsimp=0,ftol=1e-4)
    a=driver.run_legacy_refinement(s,data(),c)
    b=driver.run_legacy_refinement(s,data(),c)
    np.testing.assert_array_equal(a.final_vector,b.final_vector)
    assert a.final_objective==evaluate_refinement_trial(s,a.final_vector,data()).objective
    for h in a.stage_history:
        assert h.retained_objective==evaluate_refinement_trial(s,h.retained_vector,data()).objective
    assert s.regions[0].depth_angstrom==2


@pytest.mark.parametrize("filename",["FOR007.DAT","STARTINGFOR007.DAT"])
def test_reference_controls_and_first_invocation(filename,monkeypatch):
    path=Path(__file__).parents[1]/"reference_data"/filename
    s,c=driver.parse_legacy_refinement(path,seed=17)
    assert len(build_refinement_parameter_layout(s).parameters)==18
    assert c.stages==("SIMPARA",)
    called=[]
    def evaluate(s,v,e,**kw):
        return SimpleNamespace(objective=1.,physical_state=apply_refinement_parameters(s,v))
    def optimizer(s,v,e,**kw):
        called.append((v.copy(),kw))
        return SimpleNamespace(best_vector=v.copy(),best_objective=1.,evaluations=1)
    monkeypatch.setattr(driver,"evaluate_refinement_trial",evaluate)
    monkeypatch.setattr(driver,"run_legacy_simpara",optimizer)
    r=driver.run_legacy_refinement(s,None,replace(c,loops=1))
    np.testing.assert_array_equal(r.initial_vector,pack_refinement_parameters(s))
    np.testing.assert_array_equal(called[0][0],r.initial_vector)
    kw=called[0][1]
    assert kw["max_evaluations"]==c.maxfsimp
    assert kw["high_temperature"]==c.maxtemp and kw["low_temperature"]==c.mintemp
    assert kw["cooling"]==c.coolrate and kw["random_order"]==bool(c.paraopt)
    np.testing.assert_array_equal(kw["uncertainties"],c.plamda)


@pytest.mark.parametrize("kwargs",[dict(method="AnNe"),dict(method="NONE"),dict(loops=True),
    dict(paraopt=2),dict(maxfsimp=-1),dict(coolrate=1),dict(ftol=np.nan)])
def test_config_rejection(kwargs):
    with pytest.raises((ValueError,TypeError)):
        config(**kwargs)


def test_lowercase_prefix_and_unsupported_trial():
    assert config("anne_suffix").stages==("ANNEAL",)
    s=problem()
    s=replace(s,regions=(replace(s.regions[0],vacancy_values=(.1,)),))
    with pytest.raises(ValueError,match="hazards"):
        driver.run_legacy_refinement(s,data(),config())
    with pytest.raises(ValueError,match="PLAMDA"):
        driver.run_legacy_refinement(problem(),data(),config(plamda=()))
