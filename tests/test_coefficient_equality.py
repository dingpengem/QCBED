"""Optional physical equality constraints; tiny deterministic architecture tests."""
from dataclasses import replace
import numpy as np
import pytest
from test_projected_elastic import config
from qcbed import compile_config_from_sources,save_yaml_config,load_yaml_config,prepare_run_snapshot,run_snapshot
from qcbed.projected_elastic import results_coefficients,layout_rows,iam_grids
from qcbed.crystallography_refinement import apply_crystallographic_parameters,physical_target_layout
from qcbed.canonical_run import canonical_checkpoint_identity


def linked(tmp_path,mode='real'):
    c=config(tmp_path,mode,reflections=((2,0,0),(0,2,0)))
    return replace(c,crystallography=replace(c.crystallography,projected_elastic=tuple(replace(r,equality_group='Group 1') for r in c.crystallography.projected_elastic)))


@pytest.mark.parametrize('mode,n',[('real',1),('complex',2),('magnitude',1)])
def test_shared_layout_propagation_results_and_unlink(tmp_path,mode,n):
    c=linked(tmp_path,mode);a=compile_config_from_sources(c)
    ordinary=sum(p.category!='projected_elastic' for p in a.layout.parameters)
    assert a.parameter_count==ordinary+n
    v=a.packed_starting_vector.copy();v[-n]+=.2
    changed=apply_crystallographic_parameters(a,v)
    rows=changed.crystallography.projected_elastic
    assert rows[0]['initial']==pytest.approx(rows[1]['initial'])
    meta=results_coefficients(physical_target_layout(a),a.packed_starting_vector,latest=v)
    assert len(meta)==2 and all(m['equality_group']=='Group 1' for m in meta)
    assert meta[0]['parameter_indices']==meta[1]['parameter_indices']
    # Same propagation as explicit independent coefficients; Friedel is handled once per owner.
    independent=replace(changed.source_configuration,crystallography=replace(changed.source_configuration.crystallography,
        projected_elastic=tuple(replace(r,equality_group=None) for r in changed.source_configuration.crystallography.projected_elastic)))
    b=compile_config_from_sources(independent)
    for x,y in zip(changed.crystallography.simulation_sample.prepared.slice_types,b.crystallography.simulation_sample.prepared.slice_types):np.testing.assert_array_equal(x.phase_grating,y.phase_grating)
    assert b.parameter_count==ordinary+2*n
    assert len(rows)==2


@pytest.mark.parametrize('change,reason',[
    ({'interval':(.5,1.)},'slice'),({'mode':'complex'},'mode'),({'initial':(3.,2.)},'starting'),
    ({'lower':(80.,None),'upper':(90.,None)},'intersection'),({'refine':False},'activation')])
def test_incompatible_groups_rejected(tmp_path,change,reason):
    c=linked(tmp_path);r=list(c.crystallography.projected_elastic);r[1]=replace(r[1],**change)
    with pytest.raises(ValueError,match=reason):compile_config_from_sources(replace(c,crystallography=replace(c.crystallography,projected_elastic=tuple(r))))


def test_bound_intersection_and_stable_persistence(tmp_path):
    c=linked(tmp_path);r=list(c.crystallography.projected_elastic);r[0]=replace(r[0],lower=(-50.,None));r[1]=replace(r[1],upper=(50.,None))
    c=replace(c,crystallography=replace(c.crystallography,projected_elastic=tuple(r)))
    a=compile_config_from_sources(c);row=layout_rows(a)[0];assert (row['lower'],row['upper'])==(-50.,50.)
    p=tmp_path/'group.yml';save_yaml_config(c,p,base_directory=tmp_path);assert load_yaml_config(p)==c
    reordered=replace(c,crystallography=replace(c.crystallography,beams=tuple(reversed(c.crystallography.beams)),projected_elastic=tuple(reversed(r))))
    assert physical_target_layout(compile_config_from_sources(reordered))==physical_target_layout(a)
    independent=replace(c,crystallography=replace(c.crystallography,projected_elastic=tuple(replace(x,equality_group=None) for x in r)))
    assert canonical_checkpoint_identity(a)!=canonical_checkpoint_identity(compile_config_from_sources(independent))


def test_real_group_retains_own_imaginary_and_magnitude_own_phase(tmp_path):
    c=linked(tmp_path);r=list(c.crystallography.projected_elastic)
    r=[replace(x,initial=(10.,float(i+1))) for i,x in enumerate(r)]
    a=compile_config_from_sources(replace(c,crystallography=replace(c.crystallography,projected_elastic=tuple(r))))
    v=a.packed_starting_vector.copy();v[-1]=11
    b=apply_crystallographic_parameters(a,v)
    assert {z['initial'] for z in b.crystallography.projected_elastic}=={11+1j,11+2j}
    # Results use member phase, not the group's first member phase.
    layout=list(physical_target_layout(a));p=layout[-1];m=p['projected_elastic'];m['mode']='magnitude'
    for j,member in enumerate(m['members']):member.update(mode='magnitude',iam=[0.,2.] if j else [2.,0.])
    result=results_coefficients(layout,v)
    assert {z['initial'] for z in result}=={11+0j,11j}


def test_group_resume_and_gui_link_unlink(tmp_path):
    from qcbed.gui.results import load_results
    from qcbed.grouped_run import load_grouped_checkpoint
    c=linked(tmp_path);c=replace(c,refinement=replace(c.refinement,optimizer=replace(c.refinement.optimizer,method='SIMP',maxfsimp=0)))
    snapshot=prepare_run_snapshot(c,base_directory=tmp_path)
    run_snapshot(snapshot,cancel_requested=lambda:True);run_snapshot(snapshot,resume=True)
    result=load_results(c.run.output_directory);assert len(result.summary['projected_elastic'])==2
    altered=replace(c,crystallography=replace(c.crystallography,projected_elastic=tuple(replace(r,equality_group=None) for r in c.crystallography.projected_elastic)))
    with pytest.raises(ValueError,match='identity mismatch'):load_grouped_checkpoint(c.run.output_directory,prepare_run_snapshot(altered,base_directory=tmp_path))
    from PySide6.QtWidgets import QApplication
    from qcbed.gui.window import ConfigWindow
    from qcbed.gui.document import EditorDocument
    app=QApplication.instance() or QApplication([]);w=ConfigWindow();w.document=EditorDocument(altered,tmp_path/'gui.yml');w.refresh()
    w.workflow.accept(0,altered);assert w.validate_stage(1);assert w.validate_c_part('simulation')
    panel=w.elastic_panel;indices=[i for i,r in enumerate(panel.rows) if r['hkl'] in ((2,0,0),(0,2,0)) and r['interval']==(0.,.5)]
    panel.link_rows(indices,'User pair')
    assert {r.equality_group for r in w.document.config.crystallography.projected_elastic}=={'User pair'}
    panel=w.elastic_panel;panel.table.selectRow(indices[0]);panel.unlink_selected()
    assert all(r.equality_group is None for r in w.document.config.crystallography.projected_elastic)
    w.document.saved_config=w.document.config;w.close()


@pytest.mark.parametrize('mode,expected',[('real',2),('complex',3),('magnitude',2)])
def test_depth_plus_pair_dimension(tmp_path,mode,expected):
    from qcbed.refinement_categories import ParameterCategory
    c=linked(tmp_path,mode);c=replace(c,refinement=replace(c.refinement,active_categories=(ParameterCategory.DEPTH,ParameterCategory.PROJECTED_ELASTIC)))
    assert compile_config_from_sources(c).parameter_count==expected
    independent=replace(c,crystallography=replace(c.crystallography,projected_elastic=tuple(replace(r,equality_group=None) for r in c.crystallography.projected_elastic)))
    assert compile_config_from_sources(independent).parameter_count==1+2*(expected-1)
