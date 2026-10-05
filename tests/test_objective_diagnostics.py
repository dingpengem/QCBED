"""B6b fixed-input closure, not optimizer or historical unmatched regression."""
from dataclasses import replace
import json
import numpy as np
import pytest
from qcbed.objective_diagnostics import evaluate_preflight_objectives, residual_diagnostics
from qcbed.masked_objective import compute_selected_objective
from examples.b6a_real_input_preflight import synthetic_inputs
from examples.b6b_real_objective import evaluate_real, write_diagnostics


@pytest.fixture(scope='module')
def real():
    return evaluate_real()


def verify(report, arrays):
    assert report['status']=='PASS'
    assert report['active_parameter_count']==0
    for mode,r in report['modes'].items():
        assert r['primary_evaluations']==r['determinism_repeat_evaluations']==1
        assert all(v==0 for v in r['repeat_max_absolute_errors'].values())
        assert np.isfinite(r['ynorm']) and np.isfinite(r['objective']) and r['objective']>=0
        assert r['active_count']==(7974 if mode=='conventional' else 7281)
        a=arrays[mode]
        assert a['calculated_intensities'].shape==(1444,9)
        assert np.isfinite(a['complex_amplitudes']).all()
        assert np.all(a['calculated_intensities']>=0)
        assert np.isfinite(a['calculated_intensities']).all()
        assert np.all(a['weighted_squared_residual'][a['weights']==0]==0)
        d=r['decomposition']
        assert sum(b['active_count'] for b in d['beams'])==r['active_count']
        assert sum(b['squared_objective_contribution'] for b in d['beams'])==pytest.approx(r['objective']**2,rel=1e-14)
        assert sum(b['allocated_root_objective'] for b in d['beams'])==pytest.approx(r['objective'],rel=1e-14)
        assert [b['stable_id'] for b in d['beams']]==[b['stable_id'] for b in report['preflight']['beam_table']]
        assert d['beams'][0]['original_h']==d['beams'][0]['original_k']==0
        assert d['beams'][0]['active_count']>0
    np.testing.assert_array_equal(arrays['conventional']['complex_amplitudes'],arrays['gradient']['complex_amplitudes'])
    assert report['modes']['conventional']['gradient_option']==0
    assert report['modes']['gradient']['gradient_option']==1
    assert not np.array_equal(arrays['conventional']['experimental'],arrays['gradient']['experimental'])
    assert not np.array_equal(arrays['conventional']['variances'],arrays['gradient']['variances'])
    assert not np.array_equal(arrays['conventional']['weights'],arrays['gradient']['weights'])


def test_synthetic_determinism_identity_no_optimizer(monkeypatch):
    import qcbed.objective_diagnostics as diagnostic
    import qcbed.refinement_driver as driver
    def forbidden(*a,**kw):raise AssertionError('B6b must not run an optimizer')
    monkeypatch.setattr(driver,'run_legacy_refinement',forbidden)
    monkeypatch.setattr(driver,'_run_legacy_refinement',forbidden)
    compiled,data,baselines=synthetic_inputs()
    sample=compiled.crystallography.simulation_sample.prepared
    forward=diagnostic.simulate_tilt_series_batched
    objective=diagnostic.compute_selected_objective
    calls=[]; objectives=[]
    def checked(prepared,*a,**kw):
        assert prepared is sample
        calls.append(prepared)
        return forward(prepared,*a,**kw)
    def counted(*a,**kw):
        objectives.append(1)
        return objective(*a,**kw)
    monkeypatch.setattr(diagnostic,'simulate_tilt_series_batched',checked)
    monkeypatch.setattr(diagnostic,'compute_selected_objective',counted)
    before=compiled.packed_starting_vector.copy()
    report,arrays=evaluate_preflight_objectives(compiled,data,baseline_manifests=baselines,input_kind='synthetic')
    assert len(calls)==len(objectives)==4
    np.testing.assert_array_equal(before,compiled.packed_starting_vector)
    verify(report,arrays)


def test_real_smoke_both_modes(real):
    verify(*real)
    r,a=real
    assert r['preflight']['sample']['strict_b5_compatible']
    assert r['source_configuration']['crystallography']['cbed_basis']==dict(g1_hkl=(2,0,0),g2_hkl=(0,2,0))
    expected=[(0,0),(-2,0),(2,0),(0,-2),(0,2),(-2,-2),(2,2),(-2,2),(2,-2)]
    assert [(b['internal_H'],b['internal_K']) for b in r['preflight']['beam_table']]==expected


def test_zero_weight_placeholders_excluded():
    from qcbed.experimental import ExperimentalIntensityData
    c=np.array([[1.,2.],[3.,4.]])
    exp=ExperimentalIntensityData(2,1,0,((0,0),(2,0)),np.array([[2.,9.],[6.,8.]]),
        np.array([[1.,0.],[2.,3.]]),np.array([[1.,0.],[1.,1.]]))
    first=compute_selected_objective(c,exp,iref=-1,objective_mode='masked')
    changed=exp.intensities.copy();changed[0,1]=1e100
    altered=replace(exp,intensities=changed)
    second=compute_selected_objective(c,altered,iref=-1,objective_mode='masked')
    assert first==second
    d,a=residual_diagnostics(c,altered,second.scaling,second.objective)
    assert a['weighted_squared_residual'][0,1]==0
    assert d['zero_weight_contribution_exact']


def test_artifact_round_trip(real,tmp_path):
    report,arrays=real
    output=tmp_path/'report'
    write_diagnostics(report,arrays,output)
    saved=json.loads((output/'objective_summary.json').read_text())
    for mode in arrays:
        assert saved['modes'][mode]['objective']==report['modes'][mode]['objective']
        for name,a in arrays[mode].items():
            np.testing.assert_array_equal(np.load(output/mode/f'{name}.npy'),a)
        np.testing.assert_array_equal(np.load(output/mode/'lattice_0_0/experimental.npy'),arrays[mode]['experimental'][:,0].reshape(38,38))
        assert (output/mode/'comparison.png').is_file()
    with pytest.raises(FileExistsError):write_diagnostics(report,arrays,output)
