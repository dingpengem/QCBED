"""B7c depth/B mapping and deterministic discrete-profile control."""
from dataclasses import replace
import numpy as np
import pytest
from examples.b7a_refinement_layout import real_config
from qcbed import compile_config_from_sources
from qcbed.crystallography_refinement import apply_crystallographic_parameters, physical_target_layout
from qcbed.depth_refinement import enable_requested_depth, sample_extent, run_depth_profile, load_depth_checkpoint


@pytest.fixture(scope='module')
def compiled():
    return compile_config_from_sources(enable_requested_depth(real_config()))


def test_dynamic_source_order_and_independent_changes(compiled):
    assert compiled.parameter_count == 2
    assert compiled.packed_starting_vector.tolist() == [840., .86]
    assert [p['target_id'] for p in physical_target_layout(compiled)] == [
        'sample.requested_depth', 'site_group.Al_fcc.debye_waller_b']
    base = compiled.crystallography.simulation_sample
    depth = apply_crystallographic_parameters(compiled, np.array([838., .86]))
    thermal = apply_crystallographic_parameters(compiled, np.array([840., .88]))
    assert sample_extent(depth)['occurrence_count'] == 414
    assert sample_extent(thermal)['occurrence_count'] == 415
    for a, b in zip(base.prepared.slice_types, depth.crystallography.simulation_sample.prepared.slice_types):
        np.testing.assert_array_equal(a.phase_grating, b.phase_grating)
    assert thermal.crystallography.simulation_sample.parameters.slice_sequence_zero_based == base.parameters.slice_sequence_zero_based
    assert any(not np.array_equal(a.phase_grating, b.phase_grating) for a,b in zip(
        base.prepared.slice_types, thermal.crystallography.simulation_sample.prepared.slice_types))
    assert compiled.source_configuration.crystallography.depth_angstrom == 840.
    assert apply_crystallographic_parameters(compiled, np.array([840., .86])).source_configuration == compiled.source_configuration


@pytest.mark.parametrize('vector', [[599.,.86],[1201.,.86],[840.,.49],[840.,1.51],[np.nan,.86],[840.,np.inf]])
def test_both_bounds_reject(compiled, vector):
    with pytest.raises(ValueError):
        apply_crystallographic_parameters(compiled, np.array(vector))


@pytest.mark.parametrize('vector', [[600.,.5],[1200.,1.5]])
def test_bounds_inclusive(compiled, vector):
    assert apply_crystallographic_parameters(compiled,np.array(vector)).packed_starting_vector.tolist() == vector


@pytest.fixture
def profile(compiled, monkeypatch):
    import qcbed.depth_refinement as api
    extent = sample_extent(compiled)
    depths = [extent['plateau_lower_exclusive'],840.,float(np.nextafter(extent['effective_depth'],np.inf))]
    audit = dict(status='PASS',depth_gate='READY',candidate_depths=depths)
    calls=[]
    def evaluate(c, datasets, d, b):
        calls.append((d,b))
        e = sample_extent(apply_crystallographic_parameters(c,np.array([d,b])))
        objective=.939675432747781 + (b-.80)**2-(.86-.80)**2 + (840.-d)*.0001
        return dict(B=b,**e,objective=objective,ynorm=14870390.0181119,
                    active_count=7974,active_unique_weights=[1.],invariants='unchanged')
    monkeypatch.setattr(api,'evaluate_depth_point',evaluate)
    return audit,calls


def test_resume_exact_and_no_replayed_evaluations(compiled,profile,tmp_path):
    audit,calls=profile
    full=run_depth_profile(compiled,{},audit,tmp_path/'full')
    assert len(full['history'])==15 and len(full['sensitivity_grid'])==9
    partial=run_depth_profile(compiled,{},audit,tmp_path/'split',stop_after=7)
    assert partial['status']=='interrupted' and len(partial['history'])==7
    before=len(calls)
    resumed=run_depth_profile(compiled,{},audit,tmp_path/'split',resume=True)
    assert len(calls)-before==8
    for key in ('history','best','last_evaluated','retained_final'):
        assert resumed[key]==full[key]
    assert resumed['status']=='completed'
    assert load_depth_checkpoint(tmp_path/'split')['history']==full['history']


@pytest.mark.parametrize('after',[0,2,5])
def test_cancel_safe_boundary_and_resume(compiled,profile,tmp_path,after):
    audit,calls=profile
    # Baseline + grid are diagnostics, preceding the optimizer's safe boundaries.
    cancelled=run_depth_profile(compiled,{},audit,tmp_path/'run',cancel_requested=lambda:len(calls)>=10+after)
    assert cancelled['status']=='cancelled'
    assert len(cancelled['history'])==after
    resumed=run_depth_profile(compiled,{},audit,tmp_path/'run',resume=True)
    assert resumed['status']=='completed' and len(resumed['history'])==15
    assert resumed['history'][:after]==cancelled['history']


def test_identity_and_corruption_rejected(compiled,profile,tmp_path):
    audit,_=profile
    out=tmp_path/'run'
    run_depth_profile(compiled,{},audit,out,stop_after=1)
    changed=compile_config_from_sources(replace(compiled.source_configuration,title='different'))
    with pytest.raises(ValueError,match='identity'):
        run_depth_profile(changed,{},audit,out,resume=True)
    child=out/'depth_0'/'checkpoint_latest.json'
    child.unlink()
    with pytest.raises(ValueError,match='missing B7c child'):
        run_depth_profile(compiled,{},audit,out,resume=True)
    p=out/'checkpoint_latest.json'
    p.write_text(p.read_text().replace('interrupted','corrupted'))
    with pytest.raises(ValueError,match='checksum'):
        load_depth_checkpoint(out)
