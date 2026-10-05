"""Small physical propagation and dynamic orchestration regressions; no campaign."""
from dataclasses import replace
import numpy as np
import pytest
from test_modern_grouped_run import modern_config
from qcbed.configuration import ParameterTarget,ParameterCategory,ParameterBinding,TargetBound,BoundsConfig
from qcbed import compile_config_from_sources,prepare_run_snapshot,run_snapshot
from qcbed.crystallography_refinement import apply_crystallographic_parameters,physical_target_layout
from qcbed.config_targets import available_refinement_targets
from qcbed.grouped_run import check_grouped,evaluate_grouped
from qcbed.gui.results import load_results


def selection(config, specifications):
    return replace(config,refinement=replace(config.refinement,active_categories=None,
        parameters=tuple(ParameterBinding(t,(.9,1.1),.01) for t,_,_ in specifications),
        bounds=replace(config.refinement.bounds,specific=tuple(TargetBound(t,lo,hi) for t,lo,hi in specifications))))


def target(category,*ids):return ParameterTarget(ParameterCategory(category),ids)
DEPTH=target('depth','sample','requested_depth')
THERMAL=target('debye_waller','site_group','Al_fcc')
VOLTAGE=target('voltage')
SCAN=target('scan_coordinate','A','first')
ABS=target('absorption','site_group','Al_fcc','b')


@pytest.mark.parametrize('specs',[
    [(THERMAL,.5,1.5)],[(DEPTH,1.,20.)],[(DEPTH,1.,20.),(THERMAL,.5,1.5)],
    [(THERMAL,.5,1.5),(DEPTH,1.,20.),(VOLTAGE,100.,300.)],
    [(ABS,.5,1.5)],[(SCAN,-2.,2.)],
])
def test_dynamic_layout_and_physical_propagation(tmp_path,specs):
    c=modern_config(tmp_path,0)
    if any(t==VOLTAGE for t,_,_ in specs):
        c=replace(c,crystallography=replace(c.crystallography,site_groups=tuple(replace(g,absorption=None) for g in c.crystallography.site_groups)))
    c=selection(c,specs);s=prepare_run_snapshot(c,base_directory=tmp_path)
    assert len(check_grouped(s.compiled,s.experimental))==len(specs)
    baseline,arrays=evaluate_grouped(s,s.compiled.packed_starting_vector)
    vector=s.compiled.packed_starting_vector.copy()
    # Modify every active physical value without changing its binding or identity.
    vector*=1.01
    changed=apply_crystallographic_parameters(s.compiled,vector)
    np.testing.assert_allclose(changed.packed_starting_vector,vector,rtol=0,atol=0)
    assert [p['target_id'] for p in physical_target_layout(changed)]==[p['target_id'] for p in physical_target_layout(s.compiled)]
    row,other=evaluate_grouped(s,vector)
    if specs[0][0]!=DEPTH or len(specs)>1:
        assert not np.array_equal(arrays['amplitudes'],other['amplitudes'])
    assert row['physical_active_vector']==vector.tolist()
    assert baseline['physical_active_vector']==s.compiled.packed_starting_vector.tolist()


@pytest.mark.parametrize('method',['ANNE','SIMP','PARA','SPAR'])
def test_one_target_engine_resume_and_results(tmp_path,method):
    c=selection(modern_config(tmp_path,0),[(THERMAL,.5,1.5)])
    c=replace(c,refinement=replace(c.refinement,optimizer=replace(c.refinement.optimizer,method=method,maxfsimp=0,maxfpara=0)))
    s=prepare_run_snapshot(c,base_directory=tmp_path);events=[]
    run_snapshot(s,progress_hook=events.append,cancel_requested=lambda:len(events)>=1)
    run_snapshot(s,resume=True)
    result=load_results(c.run.output_directory)
    assert result.summary['status']=='completed'
    assert len(result.parameters)==1
    assert result.parameters[0]['label']=='Site group Al_fcc → Debye–Waller B'
    assert result.parameters[0]['best'] is not None
    assert all(len(h['physical_active_vector'])==1 for h in result.history)
    fresh=replace(c,run=replace(c.run,output_directory=str(tmp_path/'fresh'),output_filename=str(tmp_path/'fresh.dat')))
    run_snapshot(prepare_run_snapshot(fresh,base_directory=tmp_path))
    uninterrupted=load_results(fresh.run.output_directory)
    assert [(h['physical_active_vector'],h['objective']) for h in result.history]==[(h['physical_active_vector'],h['objective']) for h in uninterrupted.history]


def test_subset_order_yaml_and_constraints(tmp_path):
    from qcbed.config_serialization import save_yaml_config
    from qcbed import load_yaml_config
    c=selection(modern_config(tmp_path),[(SCAN,-2.,2.),(THERMAL,.5,1.5),(ABS,.5,1.5)])
    compiled=compile_config_from_sources(c)
    assert [p.category for p in compiled.layout.parameters]==['absorption','debye_waller','scan_coordinate']
    saved=save_yaml_config(c,tmp_path/'project.yml',base_directory=tmp_path)
    reopened=load_yaml_config(tmp_path/'project.yml')
    assert saved==reopened
    np.testing.assert_array_equal(compile_config_from_sources(reopened).packed_starting_vector,compiled.packed_starting_vector)
    v=compiled.packed_starting_vector.copy();v[0]=2.
    with pytest.raises(ValueError,match='outside bounds'):apply_crystallographic_parameters(compiled,v)
    infos=available_refinement_targets(compiled)
    assert not next(i for i in infos if i.category.value=='voltage').available
    assert all(not i.available for i in infos if i.category.value in ('atom_coordinate','dilation','constant_vacancy'))


def test_version_one_journal_keeps_legacy_replay_branch(tmp_path,monkeypatch):
    import qcbed.grouped_run as runner
    c=modern_config(tmp_path,0);s=prepare_run_snapshot(c,base_directory=tmp_path)
    original=runner._save
    def old_save(directory,cp):
        cp['version']=1
        original(directory,cp)
    monkeypatch.setattr(runner,'_save',old_save)
    events=[]
    run_snapshot(s,progress_hook=events.append,cancel_requested=lambda:len(events)>=1)
    monkeypatch.setattr(runner,'_save',original)
    assert runner.load_grouped_checkpoint(c.run.output_directory,s)['version']==1
    run_snapshot(s,resume=True)
    assert load_results(c.run.output_directory).summary['status']=='completed'


def test_zero_start_is_not_legacy_packing_floor(tmp_path):
    c=selection(modern_config(tmp_path,0),[(target('absorption','site_group','Al_fcc','a'),0.,1.)])
    c=replace(c,refinement=replace(c.refinement,optimizer=replace(c.refinement.optimizer,method='SPAR',maxfsimp=0)))
    snapshot=prepare_run_snapshot(c,base_directory=tmp_path)
    run_snapshot(snapshot,cancel_requested=lambda:True)
    from qcbed.grouped_run import load_grouped_checkpoint
    assert load_grouped_checkpoint(c.run.output_directory,snapshot)['initial']==[0.]


def test_final_output_uses_retained_scan_not_starting_scan(tmp_path,monkeypatch):
    import qcbed.grouped_run as runner
    import qcbed.output_io as output_io
    recorded=[]
    writer=output_io.write_legacy_out_dat
    def write(path,data):
        recorded.append(data)
        return writer(path,data)
    monkeypatch.setattr(output_io,'write_legacy_out_dat',write)
    c=selection(modern_config(tmp_path,0),[(SCAN,-2.,2.)])
    c=replace(c,refinement=replace(c.refinement,optimizer=replace(c.refinement.optimizer,method='SIMP',maxfsimp=0)))
    snapshot=prepare_run_snapshot(c,base_directory=tmp_path)
    # Deterministic output-layer fixture; no synthetic objective claimed as physics.
    def evaluate(snapshot,vector):
        return dict(objective=10.+float(vector[0]),ynorm=1.,physical_active_vector=list(vector)),dict(scaled=np.ones((4,9)))
    monkeypatch.setattr(runner,'evaluate_grouped',evaluate)
    run_snapshot(snapshot)
    cp=runner.load_grouped_checkpoint(c.run.output_directory,snapshot)
    retained=cp['retained_final']['physical_active_vector'][0]
    assert retained!=c.scan.corner_a[0]
    assert recorded[0].scan.corner_a[0]==retained
    output=output_io.parse_legacy_out_dat(c.run.output_filename,scan=recorded[0].scan,beam_indices_hk=recorded[0].beam_indices_hk)
    np.testing.assert_array_equal(output.scaled_intensity,np.ones((4,9)))
