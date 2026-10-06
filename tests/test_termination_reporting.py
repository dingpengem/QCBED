"""Small scalar optimizer tests; no refinement campaign."""
from dataclasses import replace
import numpy as np
import pytest
from test_projected_elastic import config
from qcbed import compile_config_from_sources
from qcbed.anneal_local_refinement import _run_stage
from qcbed.termination import termination_summary


def run_stage(tmp_path,*,cap=None,budget=10000,constant=False,cancel=False,resume=False,method='AMOEBA'):
    tmp_path.mkdir(parents=True,exist_ok=True)
    c=config(tmp_path);c=replace(c,run=replace(c.run,max_stage_iterations=cap),refinement=replace(c.refinement,optimizer=replace(c.refinement.optimizer,maxfsimp=budget,ftol=1e-30,chimin=0.)))
    compiled=compile_config_from_sources(c)
    return _run_stage(compiled,lambda v:dict(objective=1. if constant else float(1+np.sum((v-1)**2)),ynorm=1.),tmp_path/'checkpoints'/('0000_'+method),
        optimizer=method,initial=compiled.packed_starting_vector,seed=2026,settings=c.refinement.optimizer,identity={'test':'counts'},
        dynamic=True,parameter_bounds=[(r.lower,r.upper) for r in compiled.bounds.entries],local_budget=budget,
        resume=resume,cancel_requested=(lambda:True) if cancel else None)


@pytest.mark.parametrize('cap',[25,30])
def test_explicit_iteration_limit_is_exact(tmp_path,cap):
    cp=run_stage(tmp_path,cap=cap)
    assert cp['optimizer_result']['iterations']==cap
    assert cp['optimizer_result']['stop_reason']=='max_iterations'
    assert cp['optimizer_result']['evaluations']>cap
    summary=termination_summary(tmp_path,{'plan':['AMOEBA'],'status':cp['status']})
    assert summary['iterations']==cap and summary['objective_evaluations']==sum(r['objective'] is not None for r in cp['history'])
    assert summary['stop_reason']=='Maximum iterations reached'
    again=run_stage(tmp_path,cap=cap,resume=True)
    assert again['history']==cp['history']


def test_evaluation_budget_convergence_and_cancel_separate(tmp_path):
    cp=run_stage(tmp_path/'budget',budget=0)
    assert cp['optimizer_result']['stop_reason']=='max_evaluations'
    cp=run_stage(tmp_path/'convergence',constant=True)
    assert cp['optimizer_result']['converged'] and cp['optimizer_result']['iterations']==0
    summary=termination_summary(tmp_path/'convergence',{'plan':['AMOEBA'],'status':'completed'})
    assert summary['stop_reason']=='Converged'
    cp=run_stage(tmp_path/'cancel',cancel=True)
    assert cp['status']=='cancelled'
    summary=termination_summary(tmp_path/'cancel',{'plan':['AMOEBA'],'status':'cancelled'})
    assert summary['objective_evaluations']==0 and summary['stop_reason']=='Cancelled by user'


def test_amoeba_legacy_budget_excludes_initialization():
    from qcbed.amoeba import _run
    a=_run(np.array([[0.,0.,0.],[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]]),lambda v:float(1+np.sum((v-3)**2)),ftol=0.,fmin=0.,max_evaluations=20)
    assert a.stop_reason=='max_evaluations'
    assert a.search_evaluations>20 and a.evaluations==4+a.search_evaluations


@pytest.mark.parametrize('method',['ANNEAL','AMOEBA','PARAMIN','SIMPARA'])
def test_exact_cap_available_for_all_modern_engines(tmp_path,method):
    cp=run_stage(tmp_path,cap=1,method=method)
    assert cp['optimizer_result']['iterations']==1
    assert cp['optimizer_result']['stop_reason']=='max_iterations'


def test_larger_evaluation_budget_continues_past_25(tmp_path):
    cp=run_stage(tmp_path,budget=40)
    assert len(cp['history'])>25 and cp['optimizer_result']['stop_reason']=='max_evaluations'
