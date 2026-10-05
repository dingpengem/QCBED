from dataclasses import replace
from pathlib import Path
import numpy as np
import pytest
from qcbed import prepare_run_snapshot,run_snapshot,compile_config_from_sources
from qcbed.configuration import TargetBound,MeshConfig
from examples.b7a_refinement_layout import real_config
from qcbed.depth_refinement import enable_requested_depth
from qcbed.grouped_run import evaluate_grouped,load_grouped_checkpoint,grouped_history
from qcbed.gui.results import load_results


def modern_config(tmp_path,mode=1):
    c=enable_requested_depth(real_config())
    body=tmp_path/'data.dat'
    values=np.linspace(1.,2.,36).reshape(4,9);weights=np.ones((4,9));weights[0,0]=0
    with body.open('w') as f:
        for a in (values,np.ones_like(values),weights):np.savetxt(f,a.T.ravel(),fmt='%.17g')
    return replace(c,crystallography=replace(c.crystallography,depth_angstrom=4.,mesh=MeshConfig(8,8)),
        scan=replace(c.scan,scan_columns=2,scan_rows=2),experiment=replace(c.experiment,path=str(body),gradient_option=mode),
        refinement=replace(c.refinement,optimizer=replace(c.refinement.optimizer,method='ANSI',loops=1,
            maxfsimp=2,maxtemp=.01,mintemp=.0001,coolrate=.9),
            bounds=replace(c.refinement.bounds,specific=tuple(replace(b,lower=1.,upper=20.) if b.target.category.value=='depth' else b for b in c.refinement.bounds.specific))),
        run=replace(c.run,output_directory=str(tmp_path/'run'),output_filename=str(tmp_path/'run'/'final_out.dat'),max_total_evaluations=100,max_wall_time_seconds=60.))


@pytest.mark.parametrize('mode',[0,1])
def test_fresh_cancel_resume_results_snapshot(tmp_path,mode):
    c=modern_config(tmp_path,mode);snapshot=prepare_run_snapshot(c,base_directory=tmp_path)
    events=[]
    result=run_snapshot(snapshot,progress_hook=events.append,cancel_requested=lambda:len(events)>=2)
    assert result.status=='cancelled'
    cp=load_grouped_checkpoint(tmp_path/'run',snapshot)
    assert len(grouped_history(tmp_path/'run',cp))==2
    partial=load_results(tmp_path/'run')
    assert partial.summary['status']=='cancelled' and partial.summary['partial']
    assert partial.summary['best'] and partial.summary['last_evaluated']
    assert partial.comparison['mask'].sum()==35
    assert np.isnan(partial.comparison['experimental'][0,0])
    assert partial.summary['signal_mode']==('gradient' if mode else 'conventional')
    resumed=run_snapshot(snapshot,resume=True)
    assert resumed.status=='completed'
    complete=load_results(tmp_path/'run')
    assert complete.summary['retained_final'] and complete.summary['optimizer_result']
    assert len(complete.parameters)==2 and complete.parameters[0]['label']=='Sample → requested depth'
    fresh=replace(c,run=replace(c.run,output_directory=str(tmp_path/'other'),output_filename=str(tmp_path/'other.dat')))
    uninterrupted=prepare_run_snapshot(fresh,base_directory=tmp_path)
    run_snapshot(uninterrupted)
    other=load_results(tmp_path/'other')
    assert [(h['physical_active_vector'],h['objective']) for h in complete.history]==[(h['physical_active_vector'],h['objective']) for h in other.history]
    assert c.crystallography.site_groups[0].debye_waller_b==.86
    changed=replace(snapshot,experimental=replace(snapshot.experimental,intensities=snapshot.experimental.intensities+1))
    with pytest.raises(ValueError,match='identity mismatch'):run_snapshot(changed,resume=True)


def test_real_gradient_canonical_execution_equivalence(tmp_path):
    from examples.a9_real_variance_closure import prepare_reference_variance
    from qcbed.canonical_run import RunSnapshot
    c=compile_config_from_sources(enable_requested_depth(real_config()))
    datasets,_=prepare_reference_variance()
    exp=c.crystallography.adapt_dataset(datasets['gradient'],objective_mode='masked').experimental
    row,arrays=evaluate_grouped(RunSnapshot(c,exp),[840.,.86])
    assert row['active_count']==7281
    assert (row['ynorm'],row['objective'])==(13447206.106343817,1.1896911365294378)
    np.testing.assert_array_equal(arrays['experimental'],exp.intensities)
    assert sum(b['active_count'] for b in row['per_beam'])==7281


def test_public_validation_uses_grouped_resume_identity(tmp_path):
    from qcbed import validate_config,ValidationPhase
    c=modern_config(tmp_path)
    snapshot=prepare_run_snapshot(c,base_directory=tmp_path)
    run_snapshot(snapshot,cancel_requested=lambda:True)
    phases=(ValidationPhase.FIELD,ValidationPhase.COMPILE,ValidationPhase.EXECUTION,ValidationPhase.FILES)
    assert not validate_config(c,base_directory=tmp_path,phases=phases,resume=True)
    changed=replace(c,crystallography=replace(c.crystallography,depth_angstrom=5.))
    issues=validate_config(changed,base_directory=tmp_path,phases=phases,resume=True)
    assert any(i.code=='resume_incompatible' and 'identity mismatch' in i.message for i in issues)
    result=load_results(tmp_path/'run')
    assert result.summary['partial'] and not result.history and result.comparison is None


def test_modern_invalid_membership_diagnostic_and_capability_gate(tmp_path):
    from qcbed.modern_frontend import validate_modern
    c=modern_config(tmp_path);group=c.crystallography.site_groups[0]
    invalid=replace(c,crystallography=replace(c.crystallography,site_groups=(replace(group,site_indices=(0,1,2)),)))
    issues=validate_modern(invalid,tmp_path)
    assert issues and issues[0].section=='crystallography'
    assert 'cover every original site' in issues[0].message
    one=replace(c,refinement=replace(c.refinement,parameters=c.refinement.parameters[:1]))
    issues=validate_modern(one,tmp_path)
    assert not issues  # A single supported target no longer requires depth.


def test_nonreference_group_identity_and_readonly_results(tmp_path,monkeypatch):
    from qcbed.crystallography.site_groups import group_target
    from qcbed.grouped_run import check_grouped
    from qcbed import available_refinement_targets
    c=modern_config(tmp_path);old=group_target('Al_fcc');new=group_target('sample_Al')
    c=replace(c,crystallography=replace(c.crystallography,site_groups=(replace(c.crystallography.site_groups[0],id='sample_Al'),)),
        refinement=replace(c.refinement,parameters=tuple(replace(b,target=new) if b.target==old else b for b in c.refinement.parameters),
            bounds=replace(c.refinement.bounds,specific=tuple(replace(b,target=new) if b.target==old else b for b in c.refinement.bounds.specific))))
    snapshot=prepare_run_snapshot(c,base_directory=tmp_path)
    assert check_grouped(snapshot.compiled)[1]['target_id']=='site_group.sample_Al.debye_waller_b'
    events=[];run_snapshot(snapshot,progress_hook=events.append,cancel_requested=lambda:bool(events))
    files={p:p.read_bytes() for p in (tmp_path/'run').rglob('*') if p.is_file()}
    monkeypatch.setattr('qcbed.grouped_run.evaluate_grouped',lambda *args:pytest.fail('Results must not evaluate physics'))
    result=load_results(tmp_path/'run')
    assert result.parameters[1]['label']=='Site group sample_Al → Debye–Waller B'
    assert all(p.read_bytes()==value for p,value in files.items())


def test_malformed_checkpoint_readiness(tmp_path):
    from qcbed.gui.readiness import readiness
    c=modern_config(tmp_path)
    directory=Path(c.run.output_directory);directory.mkdir()
    (directory/'checkpoint_latest.json').write_text('[]')
    issues=readiness(c,tmp_path,resume=True)
    assert any(i.code=='resume_incompatible' and 'must contain an object' in i.message for i in issues)
