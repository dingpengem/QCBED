"""Literal restored GRADIENT and ONEIT callback/dataflow regressions."""
from dataclasses import replace
from types import SimpleNamespace
import numpy as np
import pytest

from qcbed.refinement_trial import _gradient, evaluate_refinement_trial
from qcbed.refinement import pack_refinement_parameters
from qcbed.refinement_driver import run_legacy_refinement
from qcbed.masked_objective import compute_selected_objective
from test_refinement import data
from test_refinement_driver import problem, config


def literal_gradient(values, rows, columns):
    # gradients.f:13-48. Explicit one-based x,y storage and per-beam loops;
    # no production reshape/padding/stencil helper participates in this oracle.
    out=np.empty_like(values)
    for beam in range(values.shape[1]):
        reflection={}
        pos=0
        for y in range(1,rows+1):
            for x in range(1,columns+1):
                reflection[x+1,y+1]=values[pos,beam]
                pos+=1
        pos=0
        for y in range(2,rows+2):
            for x in range(2,columns+2):
                total=reflection.get((x+1,y),0.)+reflection.get((x,y+1),0.)
                for dx,dy in ((-1,0),(0,-1),(1,1),(1,-1),(-1,1),(-1,-1)):
                    total=total+reflection.get((x+dx,y+dy),0.)
                out[pos,beam]=total/8-reflection[x,y]
                pos+=1
    return out


@pytest.mark.parametrize('rows,columns',[(1,1),(1,5),(4,1),(3,5),(4,4)])
def test_gradient_literal_boundaries_corners_and_beam_independence(rows,columns):
    values=np.arange(rows*columns*2,dtype=np.float64).reshape(-1,2)/7
    saved=values.copy()
    np.testing.assert_array_equal(_gradient(values,rows,columns),literal_gradient(values,rows,columns))
    changed=values.copy(); changed[:,1]*=13
    np.testing.assert_array_equal(_gradient(values,rows,columns)[:,0],_gradient(changed,rows,columns)[:,0])
    np.testing.assert_array_equal(values,saved)


def test_gradient_no_mask_or_variance_preprocessing_and_oneit_scal():
    s=replace(problem(),iref=-1.)
    e=replace(data(gradient=1),weights=np.array([[0.],[1.],[.5],[2.]]),
              variances=np.array([[0.],[1.],[2.],[3.]]))
    saved=[a.copy() for a in (e.intensities,e.variances,e.weights)]
    q=np.array([-1.]); before=q.copy()
    # SCAL a=-.5,b=3: physical PT=q/a+b=5, outside D1/D2 endpoints.
    scaled=evaluate_refinement_trial(s,q,e,coordinate_mode='legacy_simplex',objective_mode='masked')
    physical=evaluate_refinement_trial(s,np.array([5.]),e,objective_mode='masked')
    assert scaled.physical_state.regions[0].depth_angstrom==5.
    np.testing.assert_array_equal(scaled.amplitudes,physical.amplitudes)
    expected=literal_gradient(scaled.calculated_intensities,2,2)
    np.testing.assert_array_equal(scaled.objective_intensities,expected)
    assert expected[0,0]!=0  # zero-weight/zero-variance position was NOT skipped
    check=compute_selected_objective(expected,e,iref=-1.,objective_mode='masked')
    assert (scaled.scaling,scaled.objective)==(check.scaling,check.objective)
    legacy=compute_selected_objective(expected,e,iref=-1.,objective_mode='legacy')
    assert np.isfinite(legacy.scaling) and np.isnan(legacy.objective)
    for current,old in zip((e.intensities,e.variances,e.weights),saved):
        np.testing.assert_array_equal(current,old)
    np.testing.assert_array_equal(q,before)


@pytest.mark.parametrize('method',['ANNE','SIMP','PARA','SPAR','ANSI','PASI'])
def test_corrected_callback_gradient_and_retained_state(method):
    s=problem(); e=data(gradient=1); events=[]; original=repr(s)
    result=run_legacy_refinement(s,e,config(method,maxfsimp=0,maxfpara=0),
        objective_mode='masked',progress_callback=lambda *args:events.append(args))
    trials=[event[-1] for event in events if event[0]=='evaluation']
    assert len(trials)==result.evaluations
    for trial in trials:
        expected=literal_gradient(trial.calculated_intensities,2,2)
        np.testing.assert_array_equal(trial.objective_intensities,expected)
        objective=compute_selected_objective(expected,e,iref=trial.physical_state.iref,objective_mode='masked')
        assert trial.objective==objective.objective and trial.scaling==objective.scaling
    final=evaluate_refinement_trial(s,result.final_vector,e,objective_mode='masked')
    assert final.objective==result.final_objective
    assert final.physical_state==result.final_state
    assert repr(s)==original


def test_oneit_best_parmin_writeback_is_not_raw_simplex_conversion(monkeypatch):
    import qcbed.refinement_driver as driver
    from qcbed.refinement import apply_refinement_parameters
    s=problem()
    s=replace(s,regions=(replace(s.regions[0],thickness_mask=0,
        debye_waller=(2.,),debye_waller_masks=(1,)),))
    # a=-.5,b=3; q=2.5 -> raw physical -2. ONEIT:214-221 creates
    # PARMIN, then UPDPARS writes abs(-2)=2 into its DW entry (953-957).
    monkeypatch.setattr(driver,'evaluate_refinement_trial',lambda state,v,e,**kw:
        SimpleNamespace(objective=1.,physical_state=apply_refinement_parameters(state,v)))
    monkeypatch.setattr(driver,'run_legacy_amoeba',lambda *a,**kw:
        SimpleNamespace(best_evaluated_vector=np.array([2.5]),best_evaluated_objective=.5,evaluations=2))
    result=driver.run_legacy_refinement(s,data(),config('SIMP'))
    np.testing.assert_array_equal(result.final_vector,[2.])
    np.testing.assert_array_equal(result.stage_history[0].retained_vector,[2.])
    assert result.final_state.regions[0].debye_waller==(2.,)
