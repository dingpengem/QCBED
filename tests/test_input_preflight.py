"""B6a synthetic gate tests; no real physical indexing is asserted."""
from dataclasses import replace
import copy
import json
from pathlib import Path
import numpy as np
import pytest
from examples.b6a_real_input_preflight import synthetic_inputs,write_report,main
from qcbed.input_preflight import preflight_inputs,PreflightError,require_physical_config
from qcbed import compile_config_from_sources,load_yaml_config
from qcbed.configuration import CbedReflectionBasisConfig
from qcbed.crystallography import ZoneAxisCell,ReciprocalLattice

ROOT=Path(__file__).resolve().parents[1]


@pytest.fixture
def inputs():
    return synthetic_inputs()


def check(inputs):
    compiled,data,baseline=inputs
    return preflight_inputs(compiled,data,baseline_manifests=baseline,input_kind='synthetic')


def test_full_preflight_and_no_execution(monkeypatch,tmp_path):
    import qcbed.simulation as simulation
    import qcbed.multislice as multislice
    import qcbed.refinement_trial as trial
    def forbidden(*a,**kw):raise AssertionError('B6a must not evaluate propagation/objective')
    monkeypatch.setattr(simulation,'simulate_one_tilt',forbidden)
    monkeypatch.setattr(multislice,'simulate_slice_sequence',forbidden)
    monkeypatch.setattr(trial,'evaluate_refinement_trial',forbidden)
    monkeypatch.setattr(trial,'_evaluate_refinement_trial',forbidden)
    inputs=synthetic_inputs()
    before={m:tuple(a.tobytes() for a in (d.intensities,d.variances,d.weights)) for m,d in inputs[1].items()}
    report=check(inputs)
    assert report['status']=='PASS' and report['input_kind']=='synthetic'
    assert len(report['beam_table'])==9
    assert all(v==0 for v in report['max_errors'].values())
    assert report['sample']['transform']==[[1,0,0],[0,1,0],[0,0,1]]
    assert report['sample']['simulation_slice_count']==2
    assert report['sample']['effective_depth_angstrom']==6
    for mode,count in [('conventional',886),('gradient',809)]:
        r=report['modes'][mode]
        assert r['signal_shape']==r['variance_shape']==r['weight_shape']==[1444,9]
        assert r['support_per_beam']==[count]*9
        assert r['active_count']==count*9
        assert r['active_variance']['nonpositive_count']==0
        assert before[mode]==tuple(a.tobytes() for a in (inputs[1][mode].intensities,inputs[1][mode].variances,inputs[1][mode].weights))
    assert report==check(inputs) # deterministic and read-only
    write_report(report,tmp_path/'pass')
    assert json.loads((tmp_path/'pass/preflight.json').read_text())['status']=='PASS'
    with pytest.raises(FileExistsError):write_report(report,tmp_path/'pass')


def test_missing_and_unconfirmed_physical_inputs(inputs,tmp_path):
    legacy=load_yaml_config(ROOT/'examples/qcbed_config.yml')
    with pytest.raises(PreflightError,match='Missing original'):require_physical_config(legacy)
    with pytest.raises(PreflightError,match='confirmation'):
        preflight_inputs(inputs[0],inputs[1],baseline_manifests=inputs[2],input_kind='confirmed_real')
    assert main(['--config',str(ROOT/'examples/qcbed_config.yml'),'--output',str(tmp_path/'blocked')])==1
    report=json.loads((tmp_path/'blocked/preflight.json').read_text())
    assert report['status']=='FAIL' and report['code']=='missing_physical_inputs'
    assert not report['propagation_executed'] and not report['objective_executed']


def test_wrong_zone_and_strict_b5_and_missing_thermal(inputs):
    config=inputs[0].source_configuration
    bad=replace(config,crystallography=replace(config.crystallography,
        cbed_basis=CbedReflectionBasisConfig((1,0,0),(0,0,1))))
    with pytest.raises(ValueError,match='outside'):compile_config_from_sources(bad)
    bad=replace(config,microscope=replace(config.microscope,zone_axis=(1,1,0)),
        crystallography=replace(config.crystallography,zone_axis=(1,1,0),
            cbed_basis=CbedReflectionBasisConfig((1,-1,0),(0,0,1))))
    with pytest.raises(ValueError,match='Unsupported depth-resolved'):compile_config_from_sources(bad)
    bad=replace(config,crystallography=replace(config.crystallography,debye_waller=()))
    with pytest.raises(ValueError,match='No complete thermal model'):compile_config_from_sources(bad)


def test_basis_swap_is_not_auto_corrected(inputs):
    compiled,data,baseline=inputs
    original=check(inputs)
    raw=compiled.source_configuration
    basis=raw.crystallography.cbed_basis
    swapped=replace(raw,crystallography=replace(raw.crystallography,
        cbed_basis=CbedReflectionBasisConfig(basis.g2_hkl,basis.g1_hkl)))
    changed=compile_config_from_sources(swapped)
    report=check((changed,data,baseline))
    assert original['beam_table'][1]['original_h']==-1
    assert report['beam_table'][1]['original_h']==0
    assert report['beam_table'][1]['original_k']==-1
    # A coherent user swap is a different valid assignment, not detectable indexing truth.
    mixed=replace(compiled,crystallography=replace(compiled.crystallography,reflection_basis=changed.crystallography.reflection_basis))
    with pytest.raises(PreflightError,match='ordered experimental basis'):check((mixed,data,baseline))


def test_equivalent_but_different_zone_cell_rejected(inputs):
    compiled,data,baseline=inputs;c=compiled.crystallography
    new=ZoneAxisCell(c.zone_cell.original_lattice,c.zone_axis)
    mixed=replace(compiled,crystallography=replace(c,zone_cell=new))
    with pytest.raises(PreflightError,match='exact same'):check((mixed,data,baseline))


def test_reciprocal_mismatch(inputs):
    compiled,data,baseline=inputs
    damaged=copy.deepcopy(compiled)
    cell=damaged.crystallography.zone_cell
    object.__setattr__(cell,'reciprocal_lattice',ReciprocalLattice(cell.reciprocal_lattice.vectors*1.001))
    with pytest.raises(PreflightError,match='reciprocal bases differ'):check((damaged,data,baseline))


def test_variance_invalid_and_array_regression(inputs):
    compiled,data,baseline=inputs
    d=data['conventional'];v=d.variances.copy();v[d.weights>0]=0
    bad=dict(data,conventional=replace(d,variances=v))
    with pytest.raises(PreflightError,match='positive'):check((compiled,bad,baseline))
    bad=dict(data,conventional=replace(d,intensities=d.intensities+1))
    with pytest.raises(PreflightError,match='intensities differs'):check((compiled,bad,baseline))


def test_beam_order_and_a10_regressions(inputs):
    compiled,data,baseline=inputs
    d=data['conventional']
    bad=dict(data,conventional=replace(d,beams=tuple(reversed(d.beams)),lattice_coordinates=tuple(reversed(d.lattice_coordinates))))
    with pytest.raises(PreflightError,match='beam IDs/order'):check((compiled,bad,baseline))
    bad=dict(data,conventional=replace(d,scan_geometry=replace(d.scan_geometry,corner_a=(0.,0.))))
    with pytest.raises(PreflightError,match='A10 geometry'):check((compiled,bad,baseline))


def test_runtime_sample_mismatch(inputs):
    compiled,data,baseline=inputs
    state=replace(compiled.state,regions=(replace(compiled.state.regions[0],depth_angstrom=7),))
    with pytest.raises(PreflightError,match='region depth'):check((replace(compiled,state=state),data,baseline))


def test_original_structure_cannot_be_replaced_without_sample(inputs):
    compiled,data,baseline=inputs
    c=compiled.crystallography
    original=c.original_structure
    site=replace(original.sites[0],fractional_position=(.25,.25,.25))
    modified=replace(original,sites=(site,)+original.sites[1:])
    corrupted=replace(compiled,crystallography=replace(c,original_structure=modified))
    with pytest.raises(PreflightError,match='Original structure positions'):
        check((corrupted,data,baseline))
