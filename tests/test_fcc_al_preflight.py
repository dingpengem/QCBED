"""Confirmed FCC indexing and compatibility-only legacy mapping; no objective."""
import json
import numpy as np
import pytest
from qcbed.crystallography import (
    legacy_for007_001_basis_to_conventional_hkl as forward,
    conventional_hkl_to_legacy_for007_001_basis as inverse,
    fcc_conventional_reflection_diagnostic as diagnostic,
)
from examples.b6a_fcc_al_preflight import real_inputs, main


@pytest.mark.parametrize('pq,hkl', [((1,1),(2,0,0)),((1,-1),(0,2,0)),
    ((1,0),(1,1,0)),((0,1),(1,-1,0)),((-2,3),(1,-5,0)),((0,0),(0,0,0))])
def test_legacy_exact_round_trip(pq,hkl):
    assert forward(*pq)==hkl
    assert inverse(hkl)==pq


@pytest.mark.parametrize('hkl', [(1,0,0),(2,0,1),(2.,0,0),(True,0,0),(0,0)])
def test_inverse_rejection(hkl):
    with pytest.raises(ValueError):inverse(hkl)


@pytest.mark.parametrize('value',[True,1.,'1',np.bool_(False)])
def test_forward_strict_integer(value):
    with pytest.raises(ValueError):forward(value,0)


@pytest.mark.parametrize('hkl,expected', [((0,0,0),'direct'),((2,-2,0),'allowed'),
    ((1,1,1),'allowed'),((-1,3,-5),'allowed'),((1,1,0),'forbidden'),((2,1,0),'forbidden')])
def test_fcc_diagnostic(hkl,expected):
    assert diagnostic(hkl)==expected


def test_real_preflight_without_propagation(monkeypatch):
    import qcbed.simulation as simulation
    import qcbed.multislice as multislice
    import qcbed.refinement_trial as trial
    def forbidden(*a,**kw):raise AssertionError('B6b forbidden')
    monkeypatch.setattr(simulation,'simulate_one_tilt',forbidden)
    monkeypatch.setattr(multislice,'simulate_slice_sequence',forbidden)
    monkeypatch.setattr(trial,'evaluate_refinement_trial',forbidden)
    monkeypatch.setattr(trial,'_evaluate_refinement_trial',forbidden)
    report,data=real_inputs()
    assert report['status']=='PASS'
    expected=[(0,0),(-2,0),(2,0),(0,-2),(0,2),(-2,-2),(2,2),(-2,2),(2,-2)]
    for r,(h,k) in zip(report['beam_table'],expected,strict=True):
        assert (r['original_h'],r['original_k'],r['original_l'])==(h,k,0)
        assert (r['internal_H'],r['internal_K'])==(h,k)
        assert r['zone_law_value']==0
        assert r['fcc_diagnostic']==('direct' if (h,k)==(0,0) else 'allowed')
    assert all(v==0 for v in report['max_errors'].values())
    assert report['sample']['atom_count']==4
    assert report['sample']['simulation_slice_count']==2
    assert not report['propagation_executed'] and not report['objective_executed']
    for mode,count in [('conventional',7974),('gradient',7281)]:
        assert report['modes'][mode]['active_count']==count
        assert report['modes'][mode]['signal_regression']=='exact'
        assert report['modes'][mode]['a10_regression']=='exact'


def test_output_body_round_trip(tmp_path):
    from qcbed import load_yaml_config,compile_config_from_sources
    from examples.b6a_fcc_al_preflight import ROOT
    from qcbed.experimental_config import load_experimental_data
    out=tmp_path/'report'
    assert main(['--output',str(out)])==0
    config=load_yaml_config(ROOT/'examples/fcc_al_real_preflight.yml')
    from dataclasses import replace
    config=replace(config,experiment=replace(config.experiment,path=str(out/'conventional_body.dat')))
    compiled=compile_config_from_sources(config)
    parsed=load_experimental_data(compiled, base_directory=ROOT)
    _,data=real_inputs()
    for name in ('intensities','variances','weights'):
        np.testing.assert_array_equal(getattr(parsed,name),getattr(data['conventional'],name))
    assert json.loads((out/'preflight.json').read_text())['status']=='PASS'
    with pytest.raises(FileExistsError):main(['--output',str(out)])
