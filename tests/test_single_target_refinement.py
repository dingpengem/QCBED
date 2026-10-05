"""Controlled B7b integration; no multi-target search or GUI execution."""
import json
from dataclasses import replace
import numpy as np
import pytest
from qcbed.single_target_refinement import (audit_binary_masks, _run_journal,
    load_single_target_checkpoint, run_single_target_refinement)
from examples.a9_real_variance_closure import prepare_reference_variance
from examples.b7a_refinement_layout import real_config
from qcbed import compile_config_from_sources, load_yaml_config


@pytest.fixture(scope='module')
def inputs():return compile_config_from_sources(real_config()),prepare_reference_variance()[0]


def test_binary_masks_and_gradient_threshold(inputs):
    from qcbed.data_preparation.extraction import threshold_support
    _,data=inputs
    audit=audit_binary_masks(data)
    assert audit['conventional']['active_count']==7974
    assert audit['gradient']['active_count']==7281
    assert all(v['active_unique_weights']==[1.] for v in audit.values())
    np.testing.assert_array_equal(threshold_support([.875,np.nextafter(.875,1),8/9,.9,1]),[0,1,1,1,1])
    w=data['gradient'].weights.copy();w[w>0]=8/9
    with pytest.raises(ValueError,match='binary'):
        audit_binary_masks({'gradient':replace(data['gradient'],weights=w)})


def cheap(b):return dict(ynorm=2.,objective=(b-.7)**2+.1)


def journal(path,**kw):return _run_journal(.86,.95,1.05,cheap,path,dict(test='identity'),**kw)


def test_journal_resume_and_cancel_and_best(tmp_path):
    direct=journal(tmp_path/'direct')
    calls=[]
    def evaluate(b):calls.append(b);return cheap(b)
    partial=_run_journal(.86,.95,1.05,evaluate,tmp_path/'resume',dict(test='identity'),
        cancel_requested=lambda:len(calls)>=2)
    assert partial['status']=='cancelled' and len(partial['history'])==2
    assert partial['last_evaluated']==partial['history'][-1]
    resumed=journal(tmp_path/'resume',resume=True)
    for key in ('history','best','last_evaluated','retained_final'):
        assert direct[key]==resumed[key]
    assert direct['best']['objective']==min(r['objective'] for r in direct['history'] if r['objective'] is not None)
    assert load_single_target_checkpoint(tmp_path/'resume')['version']==2
    with pytest.raises(FileExistsError):journal(tmp_path/'resume')
    with pytest.raises(ValueError,match='identity'):
        _run_journal(.86,.95,1.05,cheap,tmp_path/'resume',dict(test='wrong'),resume=True)


def test_bounds_are_not_evaluated(tmp_path):
    calls=[]
    def evaluate(b):
        assert .5<=b<=1.5
        calls.append(b);return cheap(b)
    result=_run_journal(.86,.1,3.,evaluate,tmp_path/'bounds',dict(test='bounds'))
    assert any(r['bound_status']=='rejected' for r in result['history'])
    assert len(calls)==sum(r['bound_status']=='accepted' for r in result['history'])
    assert all(r['objective'] is None for r in result['history'] if r['bound_status']=='rejected')


def test_checkpoint_corruption_failure_and_pre_evaluation_cancel(tmp_path):
    cp=journal(tmp_path/'cancel',cancel_requested=lambda:True)
    assert cp['status']=='cancelled' and cp['history']==[] and cp['best'] is None
    p=tmp_path/'cancel/checkpoint_latest.json';saved=json.loads(p.read_text());saved['status']='completed'
    p.write_text(json.dumps(saved))
    with pytest.raises(ValueError,match='checksum'):load_single_target_checkpoint(p)
    with pytest.raises(ValueError,match='Cannot read'):load_single_target_checkpoint(tmp_path/'missing.json')
    def fail(b):raise RuntimeError('worker-like numerical failure')
    with pytest.raises(RuntimeError):_run_journal(.86,.95,1.05,fail,tmp_path/'failed',dict(test=1))
    assert load_single_target_checkpoint(tmp_path/'failed')['status']=='failed'
    assert load_single_target_checkpoint(tmp_path/'failed')['history']==[]


def test_single_target_guard(inputs,tmp_path):
    c,data=inputs
    raw=c.source_configuration
    disabled=replace(raw,refinement=replace(raw.refinement,parameters=()))
    with pytest.raises(ValueError,match='exactly'):
        run_single_target_refinement(compile_config_from_sources(disabled),data,tmp_path/'bad')
    assert not (tmp_path/'bad').exists()


def test_real_controlled_resume_and_fixed_state(inputs,tmp_path):
    compiled,data=inputs
    initial=compiled.source_configuration
    direct=run_single_target_refinement(compiled,data,tmp_path/'direct')
    paused=run_single_target_refinement(compiled,data,tmp_path/'resume',stop_after=3)
    assert paused['status']=='interrupted' and len(paused['history'])==3
    resumed=run_single_target_refinement(compiled,data,tmp_path/'resume',resume=True)
    assert resumed['status']==direct['status']=='completed'
    for key in ('history','best','last_evaluated','retained_final'):
        assert direct[key]==resumed[key]
    d=direct['initial_diagnostics']
    assert d['baseline']['ynorm']==14870390.0181119
    assert d['baseline']['objective']==.939675432747781
    assert all(v==0 for v in d['repeat_errors'].values())
    assert [r['B'] for r in d['local_sensitivity']]==[.8598,.86,.8602]
    for row in direct['history']:
        if row['bound_status']=='accepted':
            assert row['active_count']==7974 and row['active_unique_weights']==[1.]
            assert row['zero_weight_excluded'] and row['unrelated_state_unchanged']
            assert row['original_sites_affected']==4 and row['simulation_atom_occurrences_affected']==830
            assert np.isfinite(row['ynorm']) and np.isfinite(row['objective'])
            assert sum(b['squared_objective_contribution'] for b in row['decomposition']['beams'])==pytest.approx(row['objective']**2,rel=1e-14)
    assert direct['physical_evaluation_count']<=9
    assert compiled.source_configuration==initial
    assert initial.crystallography.site_groups[0].debye_waller_b==.86
    output=load_yaml_config(tmp_path/'direct/refined_snapshot.yml')
    assert output.crystallography.site_groups[0].debye_waller_b==direct['best']['B']
    assert 'diagnostic refinement output' in output.title
