from dataclasses import replace
import numpy as np
import pytest
from examples.b7e_gradient_depth_b_refinement import gradient_config, local_config, physical_simplex, boundary_status
from examples.a9_real_variance_closure import prepare_reference_variance
from qcbed import compile_config_from_sources
from qcbed.depth_refinement import evaluate_depth_point
from qcbed.anneal_local_refinement import run_optimizer_stage, case_settings, select_global_best
from qcbed.crystallography_refinement import physical_target_layout
from qcbed.data_preparation import prepared_dataset_manifest


@pytest.fixture(scope='module')
def problem(tmp_path_factory):
    config=gradient_config(tmp_path_factory.mktemp('gradient')/'body.dat')
    c=compile_config_from_sources(config)
    datasets,_=prepare_reference_variance()
    return c,{'gradient':datasets['gradient']}


def test_real_gradient_baseline_simplex_and_invariants(problem):
    c,data=problem
    before=prepared_dataset_manifest(data['gradient'])
    assert [p['target_id'] for p in physical_target_layout(c)]==[
        'sample.requested_depth','site_group.Al_fcc.debye_waller_b']
    baseline=evaluate_depth_point(c,data,840.,.86,signal_mode='gradient')
    assert (baseline['ynorm'],baseline['objective'])==(13447206.106343817,1.1896911365294378)
    lc=compile_config_from_sources(local_config(c.source_configuration,dict(requested_depth=1000.,B=.86)))
    vertices=physical_simplex(lc,dict(requested_depth=1000.,B=.86))
    np.testing.assert_allclose(vertices,[[1000.,.86],[995.8,.86],[1000.,.8815]],rtol=0,atol=2e-13)
    rows=[evaluate_depth_point(lc,data,*v,signal_mode='gradient') for v in vertices]
    assert rows[0]['occurrence_count']!=rows[1]['occurrence_count']
    assert rows[0]['objective']!=rows[2]['objective']
    assert all(r['active_count']==7281 and r['active_unique_weights']==[1.] for r in rows)
    assert prepared_dataset_manifest(data['gradient'])==before
    assert c.packed_starting_vector.tolist()==[840.,.86]


@pytest.mark.parametrize('depth,b,expected',[(900.,.5,('near 900 Å','near 0.5 Å²')),
    (1100.,1.5,('near 1100 Å','near 1.5 Å²')),(1000.,.86,('interior','interior'))])
def test_boundary_guard_and_inward_simplex(problem,depth,b,expected):
    from qcbed.depth_geometry import requested_depth_to_sample_geometry
    g=requested_depth_to_sample_geometry(depth,(2.0245,2.0245))
    row=dict(requested_depth=depth,B=b,effective_depth=g.effective_depth,plateau_lower_exclusive=g.plateau_lower_exclusive)
    status=boundary_status(row)
    assert (status['depth'],status['B'])==expected
    c,_=problem;lc=compile_config_from_sources(local_config(c.source_configuration,row))
    vertices=physical_simplex(lc,row)
    # Source SCAL inverse may round an unchanged endpoint by one ulp.
    # Test the displaced coordinates; never clip production vertices.
    assert 900<=vertices[1,0]<=1100
    assert .5<=vertices[2,1]<=1.5
    np.testing.assert_allclose(vertices[0],[depth,b],rtol=1e-15,atol=0)


@pytest.mark.parametrize('optimizer',['ANNEAL','AMOEBA'])
@pytest.mark.parametrize('completed',[0,2])
def test_gradient_stage_cancel_resume_identity(problem,tmp_path,monkeypatch,optimizer,completed):
    c,data=problem
    import qcbed.anneal_local_refinement as module
    calls=[]
    def fake(compiled,datasets,depth,b,*,signal_mode):
        assert signal_mode=='gradient'
        calls.append((depth,b))
        return dict(objective=.2+((depth-1000)/100)**2+(b-.8)**2,ynorm=42.,
                    active_count=7281,active_unique_weights=[1.])
    monkeypatch.setattr(module,'evaluate_depth_point',fake)
    settings=replace(case_settings(c.source_configuration),maxfsimp=3)
    kw=dict(optimizer=optimizer,initial=[1000.,.86],seed=2026,settings=settings,
        window=(900.,1100.),signal_mode='gradient',local_budget=2)
    full=run_optimizer_stage(c,data,tmp_path/'full',**kw)
    calls.clear()
    partial=run_optimizer_stage(c,data,tmp_path/'split',**kw,cancel_requested=lambda:len(calls)>=completed)
    assert partial['status']=='cancelled' and len(partial['history'])==completed
    resumed=run_optimizer_stage(c,data,tmp_path/'split',**kw,resume=True)
    for k in ('history','best','last_evaluated','retained_final','optimizer_result'):
        assert resumed[k]==full[k]
    assert resumed['identity']['problem']['signal_mode']=='gradient'
    altered=compile_config_from_sources(local_config(c.source_configuration,dict(requested_depth=1000.,B=.86)))
    with pytest.raises(ValueError,match='identity mismatch'):
        run_optimizer_stage(altered,data,tmp_path/'split',**kw,resume=True)
    best=select_global_best([full,resumed])
    assert best['objective']==full['best']['objective']


def test_reject_fractional_gradient_weights(problem,tmp_path):
    c,data=problem;d=data['gradient'];w=d.weights.copy();w[w!=0]=.9
    bad=replace(d,weights=w)
    with pytest.raises(ValueError,match='binary'):
        run_optimizer_stage(c,{'gradient':bad},tmp_path/'bad',optimizer='ANNEAL',
            initial=[1000.,.86],seed=2026,window=(900.,1100.),signal_mode='gradient')
    assert not (tmp_path/'bad').exists()
