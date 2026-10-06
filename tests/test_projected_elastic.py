"""Physical E_g architecture and fixed-evaluation checks; no refinement campaign."""
from dataclasses import replace
from pathlib import Path
import numpy as np
import pytest
from test_modern_grouped_run import modern_config
from qcbed import compile_config,compile_config_from_sources,prepare_run_snapshot,run_snapshot,load_yaml_config,save_yaml_config
from qcbed.configuration import ProjectedElasticConfig
from qcbed.projected_elastic import catalogue,iam_grids,complex_value,coordinates,layout_rows,check_ready
from qcbed.crystallography_refinement import apply_crystallographic_parameters,physical_target_layout
from qcbed.grouped_run import evaluate_grouped,check_grouped,load_grouped_checkpoint
from qcbed.phase_grating import phase_grating_from_structure_factors
from qcbed.absorption import absorptive_structure_factors


def config(tmp_path,mode='complex',initial=None,reflections=((2,0,0),)):
    c=modern_config(tmp_path,0)
    reference=__import__('qcbed.projected_elastic',fromlist=['model_reference']).model_reference(compile_config_from_sources(c)) if initial is not None else None
    records=tuple(ProjectedElasticConfig(h,(0.,.5),mode,True,initial,reference,
        lower=(0. if mode=='magnitude' else -100.,-100.),upper=(100.,100.)) for h in reflections)
    return replace(c,crystallography=replace(c.crystallography,projected_elastic=records))


@pytest.mark.parametrize('mode',['real','complex','magnitude'])
def test_iam_baseline_exact_and_dynamic_count(tmp_path,mode):
    c=config(tmp_path,mode);s=prepare_run_snapshot(c,base_directory=tmp_path)
    bare=replace(c,crystallography=replace(c.crystallography,projected_elastic=()))
    baseline=prepare_run_snapshot(bare,base_directory=tmp_path)
    n=2 if mode=='complex' else 1
    assert s.compiled.parameter_count==baseline.compiled.parameter_count+n
    check_grouped(s.compiled,s.experimental)
    row,a=evaluate_grouped(s,s.compiled.packed_starting_vector)
    other,b=evaluate_grouped(baseline,baseline.compiled.packed_starting_vector)
    assert row['objective']==other['objective']
    for name in ('amplitudes','calculated','scaled'):np.testing.assert_array_equal(a[name],b[name])
    for x,y in zip(s.compiled.crystallography.simulation_sample.prepared.slice_types,baseline.compiled.crystallography.simulation_sample.prepared.slice_types):
        np.testing.assert_array_equal(x.phase_grating,y.phase_grating)


@pytest.mark.parametrize('mode',['real','complex','magnitude'])
def test_only_intended_slice_and_friedel_pair_change(tmp_path,mode):
    s=prepare_run_snapshot(config(tmp_path,mode),base_directory=tmp_path)
    c=s.compiled;v=c.packed_starting_vector.copy();index=next(i for i,p in enumerate(c.layout.parameters) if p.category=='projected_elastic');v[index]+=.25
    if mode=='complex':v[index+1]+=.17
    changed=apply_crystallographic_parameters(c,v)
    before=c.crystallography.simulation_sample;after=changed.crystallography.simulation_sample
    row=c.crystallography.projected_elastic[0];r=changed.crystallography.projected_elastic[0]
    if mode=='real': assert r['initial'].imag==row['initial'].imag
    if mode=='magnitude':assert np.angle(r['initial'])==pytest.approx(np.angle(row['iam']))
    # Reconstruct independently using exactly the specified complex coefficient and its conjugate.
    elastic=iam_grids(c)[0].copy();plus,minus=row['positions'];elastic[plus]=r['initial'];elastic[minus]=r['initial'].conjugate()
    stored=before.prepared.slice_types[0];p=before.parameters;a=p.slices[0]
    absorption=absorptive_structure_factors(before.model.species,before.model.atoms_by_slice[0],stored.reciprocal_grid,
        voltage_kv=p.applied_voltage_kv,applied_debye_waller_angstrom_squared=a.debye_waller_angstrom_squared,
        coefficients=p.absorption_coefficients,vacancy_fraction=0.,absorption_raw=before.model.absorption_raw)
    expected=phase_grating_from_structure_factors(elastic,absorption,stored.reciprocal_grid,voltage_kv=p.applied_voltage_kv,
        dsap_inverse_angstrom=p.dsap_inverse_angstrom,method=p.phase_grating_method)
    np.testing.assert_array_equal(after.prepared.slice_types[0].phase_grating,expected)
    np.testing.assert_array_equal(after.prepared.slice_types[1].phase_grating,before.prepared.slice_types[1].phase_grating)
    assert not np.array_equal(expected,stored.phase_grating)


def test_fcc_teaching_values_all_original_reflections(tmp_path):
    c=modern_config(tmp_path)
    root=Path(__file__).resolve().parents[1]
    c=replace(c,crystallography=replace(c.crystallography,structure=replace(c.crystallography.structure,path=str(root/'reference_data/Al_FCC.cif'),format='cif'),
        slice_boundaries_fractional=(0.,1.),site_groups=tuple(replace(g,absorption=None) for g in c.crystallography.site_groups)))
    compiled=compile_config_from_sources(c);rows=catalogue(compiled)
    r=next(r for r in rows if r['hkl']==(2,0,0));assert r['mode']=='real'
    assert r['iam']==pytest.approx(19.591754648483377-1.7592060366859781e-9j,abs=1e-13)
    from qcbed.crystallography.index_transform import internal_HK_to_original_hkl
    assert {b.original_hkl for b in compiled.crystallography.indexed_beams}=={(0,0,0),(-2,0,0),(2,0,0),(0,-2,0),(0,2,0),(-2,-2,0),(2,2,0),(-2,2,0),(2,-2,0)}
    for b in compiled.crystallography.indexed_beams:
        assert internal_HK_to_original_hkl(b.internal_hk,compiled.crystallography.zone_cell)==b.original_hkl
    assert max(abs(r['iam'].imag) for r in rows)<1e-8


def test_noncentrosymmetric_complex_and_zero_safe_override(tmp_path):
    c=config(tmp_path)
    compiled=compile_config_from_sources(c);structure=compiled.crystallography.original_structure
    sites=list(structure.sites);sites[0]=replace(sites[0],fractional_position=(.123,.21,.1))
    complex_compiled=compile_config(c,structure=replace(structure,sites=tuple(sites)))
    assert abs(complex_compiled.crystallography.projected_elastic[0]['iam'].imag)>.1
    # An empty physical slice yields exactly zero IAM; no occupancy-policy change.
    zero=replace(c,crystallography=replace(c.crystallography,slice_boundaries_fractional=(0.,.1,.2,.5,1.),
        projected_elastic=(replace(c.crystallography.projected_elastic[0],interval=(.1,.2)),)))
    z=compile_config_from_sources(zero);assert z.crystallography.projected_elastic[0]['iam']==0
    v=z.packed_starting_vector.copy();v[-2:]=[.5,.2]
    updated=apply_crystallographic_parameters(z,v)
    assert updated.crystallography.projected_elastic[0]['initial']==.5+.2j
    record=replace(zero.crystallography.projected_elastic[0],mode='magnitude')
    with pytest.raises(ValueError,match='phase'):
        compile_config_from_sources(replace(zero,crystallography=replace(zero.crystallography,projected_elastic=(record,))))


def test_yaml_reordering_bounds_stale_and_reset(tmp_path):
    c=config(tmp_path,reflections=((2,0,0),(0,2,0)))
    compiled=compile_config_from_sources(c);row=compiled.crystallography.projected_elastic[0]
    records=tuple(replace(r,initial=(1.,2.),reference=row['reference']) for r in c.crystallography.projected_elastic)
    c=replace(c,crystallography=replace(c.crystallography,projected_elastic=records))
    p=tmp_path/'project.yml';save_yaml_config(c,p,base_directory=tmp_path);loaded=load_yaml_config(p)
    assert loaded==c
    reversed_c=replace(c,crystallography=replace(c.crystallography,beams=tuple(reversed(c.crystallography.beams)),projected_elastic=tuple(reversed(records))))
    assert physical_target_layout(compile_config_from_sources(c))==physical_target_layout(compile_config_from_sources(reversed_c))
    changed=replace(c,crystallography=replace(c.crystallography,site_groups=tuple(replace(g,debye_waller_b=.9) for g in c.crystallography.site_groups)))
    changed_compiled=compile_config_from_sources(changed)
    assert all(r['initial']==1+2j and r['stale'] for r in changed_compiled.crystallography.projected_elastic)
    with pytest.raises(ValueError,match='stale IAM'):check_ready(changed_compiled)
    reset=replace(changed,crystallography=replace(changed.crystallography,projected_elastic=tuple(replace(r,initial=None,reference=None) for r in records)))
    check_ready(compile_config_from_sources(reset))
    v=compiled.packed_starting_vector.copy();v[-1]=101
    with pytest.raises(ValueError,match='outside bounds'):apply_crystallographic_parameters(compiled,v)


@pytest.mark.parametrize('method',['ANNE','SIMP','PARA','SPAR'])
def test_all_optimizers_zero_imaginary_resume_results(tmp_path,method):
    from qcbed.gui.results import load_results
    c=config(tmp_path,initial=(9.,0.))
    c=replace(c,refinement=replace(c.refinement,optimizer=replace(c.refinement.optimizer,method=method,maxfsimp=0,maxfpara=0)))
    s=prepare_run_snapshot(c,base_directory=tmp_path);events=[]
    check_grouped(s.compiled,s.experimental)
    run_snapshot(s,progress_hook=events.append,cancel_requested=lambda:len(events)>=1)
    run_snapshot(s,resume=True)
    result=load_results(c.run.output_directory);assert result.summary['status']=='completed'
    r=result.summary['projected_elastic'][0]
    assert r['initial']==9+0j and r['latest'] is not None and r['best'] is not None and r['retained'] is not None
    assert r['hkl']==(2,0,0) and r['interval']==(0.,.5)
    assert 'Projected elastic coefficient' in result.parameters[-1]['label']
    other=replace(c,run=replace(c.run,output_directory=str(tmp_path/'fresh'),output_filename=str(tmp_path/'fresh.dat')))
    run_snapshot(prepare_run_snapshot(other,base_directory=tmp_path))
    fresh=load_results(other.run.output_directory)
    assert [(h['physical_active_vector'],h['objective']) for h in result.history]==[(h['physical_active_vector'],h['objective']) for h in fresh.history]
    incompatible=replace(c,crystallography=replace(c.crystallography,projected_elastic=(replace(c.crystallography.projected_elastic[0],hkl=(0,2,0)),)))
    with pytest.raises(ValueError,match='identity mismatch'):load_grouped_checkpoint(c.run.output_directory,prepare_run_snapshot(incompatible,base_directory=tmp_path))


def test_central_duplicate_negative_owner_rejected(tmp_path):
    with pytest.raises(ValueError,match='read-only'):ProjectedElasticConfig((0,0,0),(0.,1.))
    with pytest.raises(ValueError,match='Friedel'):ProjectedElasticConfig((-2,0,0),(0.,1.))
    c=config(tmp_path);c=replace(c,crystallography=replace(c.crystallography,projected_elastic=c.crystallography.projected_elastic*2))
    with pytest.raises(ValueError,match='Duplicate'):compile_config_from_sources(c)


def test_gui_editor_iam_save_reload_and_results_labels(tmp_path):
    from PySide6.QtWidgets import QApplication
    from PySide6.QtCore import Qt
    from qcbed.gui.window import ConfigWindow
    from qcbed.gui.document import EditorDocument
    from qcbed.gui.results_view import ResultsView
    app=QApplication.instance() or QApplication([])
    w=ConfigWindow();c=modern_config(tmp_path,0)
    w.document=EditorDocument(c,tmp_path/'editor.yml');w.refresh()
    w.workflow.accept(0,c);assert w.validate_stage(1);assert w.validate_c_part('simulation')
    panel=w.elastic_panel
    i=next(i for i,r in enumerate(panel.rows) if r['hkl']==(2,0,0) and r['interval']==(0.,.5))
    r=panel.rows[i]
    assert complex(panel.table.item(i,3).text().replace(' ','').replace('i','j'))==r['iam']
    panel.table.item(i,0).setCheckState(Qt.CheckState.Checked)
    panel.table.cellWidget(i,4).setCurrentIndex(1)
    panel.table.item(i,5).setText('9 +0i');panel.table.item(i,6).setText('-100, -100');panel.table.item(i,7).setText('100, 100')
    panel.table.setCurrentCell(i,5)
    assert 'phase=' in panel.details_text.text() and '|E_g|=' in panel.details_text.text()
    assert 'Imaginary' == panel.bound_labels[1].text()
    panel.bound_fields[1,6].setText('-50')
    panel.bound_fields[1,6].textEdited.emit('-50')
    assert panel.table.item(i,6).text()=='-100, -50'
    assert panel.apply()
    assert w.document.config.crystallography.projected_elastic[0].lower==(-100.,-50.)
    assert w.document.config.crystallography.projected_elastic[0].initial==(9.,0.)
    assert w.save_document(tmp_path/'saved.yml')
    assert not w.pending, 'Saving committed coefficient values must clear draft state'
    loaded=load_yaml_config(tmp_path/'saved.yml')
    assert loaded.crystallography.projected_elastic==w.document.config.crystallography.projected_elastic
    loaded=replace(loaded,run=replace(loaded.run,output_directory=str(tmp_path/'recorded'),output_filename=str(tmp_path/'recorded/final_out.dat')),
        refinement=replace(loaded.refinement,optimizer=replace(loaded.refinement.optimizer,method='SIMP',maxfsimp=0)))
    run_snapshot(prepare_run_snapshot(loaded,base_directory=tmp_path))
    view=ResultsView();view.load(loaded.run.output_directory)
    assert view.elastic.rowCount()==1
    assert view.elastic.item(0,0).text()=='(2, 0, 0)'
    view.select_elastic(0)
    assert '(2, 0, 0)' in view.physical_target.currentText()
    assert '[2, 0, 0]' in view.beams.currentText()
    w.document.saved_config=w.document.config;w.close();view.close()


def test_equivalent_legacy_multiplier_math_including_amplitudes(tmp_path):
    from qcbed.refinement import RefinementReflection
    from qcbed.refinement_trial import evaluate_refinement_trial
    from qcbed.projected_elastic import model_reference
    c=modern_config(tmp_path,0)
    c=replace(c,crystallography=replace(c.crystallography,slice_boundaries_fractional=(0.,1.)))
    base=prepare_run_snapshot(c,base_directory=tmp_path)
    row=next(r for r in catalogue(base.compiled) if r['hkl']==(2,0,0))
    value=1.25*row['iam']
    record=ProjectedElasticConfig((2,0,0),(0.,1.),'complex',True,(value.real,value.imag),model_reference(base.compiled),(-100.,-100.),(100.,100.))
    c=replace(c,crystallography=replace(c.crystallography,projected_elastic=(record,)))
    s=prepare_run_snapshot(c,base_directory=tmp_path)
    _,new=evaluate_grouped(s,s.compiled.packed_starting_vector)
    state=replace(base.compiled.state,reflections=((RefinementReflection(2,0,1.25,1.,'NN'),RefinementReflection(-2,0,1.25,1.,'NN')),))
    old=evaluate_refinement_trial(state,base.compiled.packed_starting_vector,base.experimental,
        execution=base.compiled.execution,objective_mode='masked')
    # Scalar complex assignment vs NumPy array multiplication can differ by last-bit rounding.
    np.testing.assert_allclose(new['amplitudes'],old.amplitudes,rtol=1e-14,atol=1e-16)
    np.testing.assert_allclose(new['calculated'],old.calculated_intensities,rtol=1e-14,atol=1e-16)


def test_recast_identity_and_coefficient_only_layout(tmp_path):
    c=config(tmp_path,reflections=((0,2,0),))
    c=replace(c,microscope=replace(c.microscope,zone_axis=(1,0,0)),
        crystallography=replace(c.crystallography,zone_axis=(1,0,0),cbed_basis=replace(c.crystallography.cbed_basis,g1_hkl=(0,2,0),g2_hkl=(0,0,2))),
        refinement=replace(c.refinement,parameters=(),bounds=replace(c.refinement.bounds,specific=())))
    s=prepare_run_snapshot(c,base_directory=tmp_path)
    assert s.compiled.parameter_count==2
    layout=check_grouped(s.compiled,s.experimental)
    assert all(r['projected_elastic']['hkl']==[0,2,0] for r in layout)
    row=s.compiled.crystallography.projected_elastic[0]
    sample=s.compiled.crystallography.simulation_sample;grid=sample.prepared.slice_types[0].reciprocal_grid
    hk=tuple(int(a[row['positions'][0]]) for a in (grid.indices.h,grid.indices.k))
    assert hk!=(0,2)
    vector=s.compiled.packed_starting_vector.copy();vector[0]+=.1
    changed=apply_crystallographic_parameters(s.compiled,vector)
    assert changed.crystallography.projected_elastic[0]['hkl']==(0,2,0)
    assert changed.crystallography.projected_elastic[0]['initial'].real==vector[0]


def test_mixed_zero_scan_start_is_not_legacy_packing_floor(tmp_path):
    from qcbed.configuration import ParameterBinding,ParameterTarget,ParameterCategory,TargetBound
    c=config(tmp_path)
    target=ParameterTarget(ParameterCategory.SCAN_COORDINATE,('A','first'))
    c=replace(c,scan=replace(c.scan,corner_a=(0.,c.scan.corner_a[1])),refinement=replace(c.refinement,
        parameters=c.refinement.parameters+(ParameterBinding(target,(1.05,.95),.0005),),
        bounds=replace(c.refinement.bounds,specific=c.refinement.bounds.specific+(TargetBound(target,-1.,1.),))))
    s=prepare_run_snapshot(c,base_directory=tmp_path)
    assert len(s.compiled.range_factors)==s.compiled.parameter_count
    i=next(i for i,p in enumerate(s.compiled.layout.parameters) if p.category=='scan_coordinate')
    assert s.compiled.packed_starting_vector[i]==0.
    check_grouped(s.compiled,s.experimental)
    changed=apply_crystallographic_parameters(s.compiled,s.compiled.packed_starting_vector)
    assert changed.state.scan.corner_a[0]==0.


def test_zero_magnitude_start_retains_iam_phase_not_zero_phase(tmp_path):
    from qcbed.projected_elastic import results_coefficients
    c=config(tmp_path,'magnitude',initial=(0.,0.))
    s=prepare_run_snapshot(c,base_directory=tmp_path);check_grouped(s.compiled,s.experimental)
    row=s.compiled.crystallography.projected_elastic[0]
    initial=s.compiled.packed_starting_vector
    metadata=results_coefficients(physical_target_layout(s.compiled),initial)[0]
    assert metadata['initial']==0 and metadata['phase_radians'] is None
    v=initial.copy();v[-1]=2.
    changed=apply_crystallographic_parameters(s.compiled,v)
    value=changed.crystallography.projected_elastic[0]['initial']
    assert abs(value)==pytest.approx(2.)
    assert np.angle(value)==pytest.approx(np.angle(row['iam']),abs=1e-20)
    result=results_coefficients(physical_target_layout(s.compiled),initial,displayed=v)[0]
    assert result['displayed']==value
