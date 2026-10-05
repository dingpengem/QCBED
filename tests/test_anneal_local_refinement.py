"""B7d workflow uses existing algorithms and exact completed-prefix replay."""
from dataclasses import replace
import json
import numpy as np
import pytest
from examples.b7a_refinement_layout import real_config
from qcbed import compile_config_from_sources
from qcbed.depth_refinement import enable_requested_depth,sample_extent
from qcbed.crystallography_refinement import apply_crystallographic_parameters,physical_target_layout
from qcbed.anneal_local_refinement import (_run_stage,case_settings,boundary_hit,
    select_global_best,load_stage_checkpoint)


@pytest.fixture(scope='module')
def compiled():return compile_config_from_sources(enable_requested_depth(real_config()))


@pytest.fixture
def problem(compiled):
    calls=[]
    def evaluate(d,b):
        calls.append((d,b))
        # Explicit smooth synthetic objective, not a claim about real thickness.
        from qcbed.depth_geometry import requested_depth_to_sample_geometry
        g=requested_depth_to_sample_geometry(d,(2.0245,2.0245))
        return dict(requested_depth=d,B=b,effective_depth=g.effective_depth,
            plateau_lower_exclusive=g.plateau_lower_exclusive,occurrence_count=g.occurrence_count,
            objective=.2+((d-1060)/100)**2+(b-.8)**2,ynorm=42.)
    def run(path,**kw):
        defaults=dict(optimizer='ANNEAL',initial=[1050.,.86],seed=2026,
            settings=case_settings(compiled.source_configuration),identity={'test':'synthetic'},local_budget=10)
        defaults.update(kw)
        return _run_stage(compiled,evaluate,path,**defaults)
    return run,calls


def test_layout_settings_and_group_propagation(compiled):
    assert compiled.packed_starting_vector.tolist()==[840.,.86]
    assert [x['target_id'] for x in physical_target_layout(compiled)]==[
        'sample.requested_depth','site_group.Al_fcc.debye_waller_b']
    settings=case_settings(compiled.source_configuration)
    assert (settings.maxtemp,settings.mintemp,settings.coolrate,settings.maxfsimp)==(.05,.0001,.9,59)
    changed=apply_crystallographic_parameters(compiled,np.array([1050.,.8]))
    assert all(s.debye_waller_b==.8 and s.physical_group_id=='Al_fcc'
        for slab in changed.crystallography.simulation_sample.depth_sample.slices for s in slab.sites)
    assert sample_extent(changed)['effective_depth']>=1050
    assert compiled.packed_starting_vector.tolist()==[840.,.86]


@pytest.mark.parametrize('optimizer',['ANNEAL','AMOEBA'])
def test_exact_resume_and_no_replayed_physics(problem,tmp_path,optimizer):
    run,calls=problem
    full=run(tmp_path/'full',optimizer=optimizer)
    split=run(tmp_path/'split',optimizer=optimizer,stop_after=4)
    assert split['status']=='interrupted'
    before=len(calls)
    resumed=run(tmp_path/'split',optimizer=optimizer,resume=True)
    assert len(calls)-before==sum(x['bound_status']=='accepted' for x in full['history'][4:])
    for k in ('history','best','last_evaluated','retained_final','optimizer_result'):
        assert full[k]==resumed[k]
    assert load_stage_checkpoint(tmp_path/'split')['status']=='completed'


@pytest.mark.parametrize('optimizer',['ANNEAL','AMOEBA'])
@pytest.mark.parametrize('count',[0,2])
def test_cancel_then_resume(problem,tmp_path,optimizer,count):
    run,calls=problem
    result=run(tmp_path/'run',optimizer=optimizer,cancel_requested=lambda:len(calls)>=count)
    assert result['status']=='cancelled'
    assert sum(x['bound_status']=='accepted' for x in result['history'])==count
    resumed=run(tmp_path/'run',optimizer=optimizer,resume=True)
    assert resumed['status']=='completed'
    assert resumed['history'][:len(result['history'])]==result['history']


def test_seed_reproducibility_and_global_handoff(problem,tmp_path):
    run,_=problem
    a=run(tmp_path/'a');b=run(tmp_path/'b');c=run(tmp_path/'c',seed=2027)
    assert a['history']==b['history'] and a['history']!=c['history']
    best=select_global_best([a,c])
    assert best['objective']==min(a['best']['objective'],c['best']['objective'])
    local=run(tmp_path/'local',optimizer='AMOEBA',initial=[best['requested_depth'],best['B']])
    first=local['history'][0]
    np.testing.assert_allclose([first['requested_depth'],first['B']],
        [best['requested_depth'],best['B']],rtol=1e-15,atol=1e-15)
    assert local['status']=='completed'


def test_bounds_reject_without_physics(problem,tmp_path,compiled):
    run,calls=problem
    settings=replace(case_settings(compiled.source_configuration),maxtemp=1.,maxfsimp=8)
    cp=run(tmp_path/'run',settings=settings)
    rejected=[r for r in cp['history'] if r['bound_status']=='rejected']
    assert rejected and all(r['objective'] is None and r['ynorm'] is None for r in rejected)
    assert all(1000<=d<=1100 and .5<=b<=1.5 for d,b in calls)
    assert cp['best']['objective'] is not None


@pytest.mark.parametrize('point',[[840.,.86],[1050.,.49],[1101.,.86],[np.nan,.86]])
def test_invalid_initial(problem,tmp_path,point):
    run,_=problem
    with pytest.raises(ValueError,match='initial'):run(tmp_path/'run',initial=point)
    assert not (tmp_path/'run').exists()


def test_checkpoint_protection(problem,tmp_path):
    run,_=problem
    run(tmp_path/'run',stop_after=2)
    with pytest.raises(FileExistsError):run(tmp_path/'run')
    with pytest.raises(ValueError,match='identity'):run(tmp_path/'run',resume=True,seed=9)
    p=tmp_path/'run'/'checkpoint_latest.json';value=json.loads(p.read_text())
    value['history'][0]['B']=.99;p.write_text(json.dumps(value))
    with pytest.raises(ValueError,match='checksum'):run(tmp_path/'run',resume=True)


def test_failure_not_cancel(compiled,tmp_path):
    def fail(d,b):raise RuntimeError('worker problem')
    with pytest.raises(RuntimeError):
        _run_stage(compiled,fail,tmp_path/'run',optimizer='ANNEAL',initial=[1050.,.86],seed=1,
            settings=case_settings(compiled.source_configuration),identity={})
    assert load_stage_checkpoint(tmp_path/'run')['status']=='failed'


@pytest.mark.parametrize('d,lo,hi,hit',[(1000,999,1001,True),(1001,999,1001,True),
    (1099,1098,1100.1,True),(1050,1049,1051,False)])
def test_boundary_geometry(d,lo,hi,hit):
    assert boundary_hit(dict(requested_depth=d,plateau_lower_exclusive=lo,effective_depth=hi)) is hit


@pytest.mark.parametrize('changes',[dict(maxtemp=0.),dict(mintemp=-1.),dict(coolrate=1.),dict(maxfsimp=-1)])
def test_invalid_settings_use_production_validation(problem,compiled,tmp_path,changes):
    run,_=problem
    with pytest.raises(ValueError):
        run(tmp_path/'invalid',settings=replace(case_settings(compiled.source_configuration),**changes))
    assert not (tmp_path/'invalid').exists()


def test_boundary_immediately_adjacent_plateau():
    assert boundary_hit(dict(requested_depth=1001.09,plateau_lower_exclusive=1000.103,effective_depth=1002.1275))
    assert not boundary_hit(dict(requested_depth=1004.831,plateau_lower_exclusive=1004.152,effective_depth=1006.1765))


def test_wider_case_window_identity_and_bounds(problem,tmp_path):
    run,calls=problem
    result=run(tmp_path/'wide',initial=[1000.,.86],window=(900.,1100.),stop_after=5)
    assert result['identity']['window']==[900.,1100.]
    assert all(900<=d<=1100 for d,b in calls)
    with pytest.raises(ValueError,match='identity'):
        run(tmp_path/'wide',initial=[1000.,.86],resume=True)
    resumed=run(tmp_path/'wide',initial=[1000.,.86],window=(900.,1100.),resume=True)
    assert resumed['status']=='completed'
    assert not boundary_hit(dict(requested_depth=1000.,plateau_lower_exclusive=999.,effective_depth=1001.),(900.,1100.))
    with pytest.raises(ValueError,match='legal'):
        run(tmp_path/'invalid_window',window=(500.,1100.))
