"""Compiled coordinate counts, GUI validation, and physical group identities."""
from dataclasses import replace
import pytest
from test_projected_elastic import config
from qcbed import compile_config_from_sources,save_yaml_config,load_yaml_config,prepare_run_snapshot,run_snapshot
from qcbed.gui.resolved_layout import resolved_layout
from qcbed.refinement_categories import ParameterCategory


@pytest.mark.parametrize('mode,n',[('real',1),('complex',2),('magnitude',1)])
def test_exact_dimensions_and_shared_bounds(tmp_path,mode,n):
    c=config(tmp_path,mode,initial=(10.,0.),reflections=((2,0,0),(0,2,0)))
    c=replace(c,refinement=replace(c.refinement,active_categories=(ParameterCategory.PROJECTED_ELASTIC,)))
    assert compile_config_from_sources(c).parameter_count==2*n
    c=replace(c,crystallography=replace(c.crystallography,projected_elastic=tuple(replace(r,equality_group='Pair') for r in c.crystallography.projected_elastic)))
    compiled=compile_config_from_sources(c)
    counts,text=resolved_layout(compiled)
    assert counts['active_physical_targets']==2 and counts['optimizer_parameters']==n
    assert len(compiled.packed_starting_vector)==n and len(compiled.bounds.entries)==n
    assert counts['mandatory_friedel_relationships']==2
    c=replace(c,refinement=replace(c.refinement,active_categories=(ParameterCategory.DEPTH,ParameterCategory.PROJECTED_ELASTIC)))
    assert compile_config_from_sources(c).parameter_count==1+n
    c=replace(c,crystallography=replace(c.crystallography,projected_elastic=tuple(replace(r,equality_group=None) for r in c.crystallography.projected_elastic)))
    assert compile_config_from_sources(c).parameter_count==1+2*n


def test_two_groups_writeback_yaml_resume_results(tmp_path):
    from qcbed.crystallography_refinement import apply_crystallographic_parameters
    from qcbed.gui.results import load_results
    from qcbed.grouped_run import load_grouped_checkpoint
    c=config(tmp_path,'real',initial=(10.,0.),reflections=((2,0,0),(2,2,0),(0,2,0),(2,-2,0)))
    records=tuple(replace(r,equality_group='Axial' if r.hkl in ((2,0,0),(0,2,0)) else 'Diagonal') for r in c.crystallography.projected_elastic)
    c=replace(c,crystallography=replace(c.crystallography,projected_elastic=records),refinement=replace(c.refinement,active_categories=(ParameterCategory.DEPTH,ParameterCategory.PROJECTED_ELASTIC)))
    a=compile_config_from_sources(c);counts,text=resolved_layout(a)
    assert counts==dict(active_physical_targets=5,shared_parameter_groups=2,independent_targets=1,mandatory_friedel_relationships=4,optimizer_parameters=3)
    assert 'Shared Parameter Group: Axial' in text and 'Shared Parameter Group: Diagonal' in text
    v=a.packed_starting_vector.copy();v[-2:]=[11.,12.]
    altered=apply_crystallographic_parameters(a,v)
    values={r['record'].equality_group:r['initial'].real for r in altered.crystallography.projected_elastic}
    assert set(values.values())=={11.,12.}
    assert all(r['initial'].real==values[r['record'].equality_group] for r in altered.crystallography.projected_elastic)
    path=tmp_path/'groups.yml';save_yaml_config(c,path,base_directory=tmp_path)
    loaded=load_yaml_config(path);assert resolved_layout(compile_config_from_sources(loaded))==(counts,text)
    snap=prepare_run_snapshot(c,base_directory=tmp_path);run_snapshot(snap,cancel_requested=lambda:True)
    load_grouped_checkpoint(c.run.output_directory,prepare_run_snapshot(loaded,base_directory=tmp_path))
    result=load_results(c.run.output_directory)
    assert len(result.parameters)==3 and len(result.summary['projected_elastic'])==4
    assert {r['equality_group'] for r in result.summary['projected_elastic']}=={'Axial','Diagonal'}
    assert sum('Shared refinement parameter group' in p['label'] for p in result.parameters)==2


def test_three_member_group(tmp_path):
    c=config(tmp_path,'complex',initial=(10.,1.),reflections=((2,0,0),(0,2,0),(2,2,0)))
    c=replace(c,crystallography=replace(c.crystallography,projected_elastic=tuple(replace(r,equality_group='Triple') for r in c.crystallography.projected_elastic)),refinement=replace(c.refinement,active_categories=(ParameterCategory.PROJECTED_ELASTIC,)))
    counts,_=resolved_layout(compile_config_from_sources(c))
    assert counts['active_physical_targets']==3 and counts['optimizer_parameters']==2


def test_duplicate_membership_rejected(tmp_path):
    c=config(tmp_path,'real',initial=(10.,0.),reflections=((2,0,0),(0,2,0)))
    records=tuple(replace(r,equality_group='One') for r in c.crystallography.projected_elastic)
    c=replace(c,crystallography=replace(c.crystallography,projected_elastic=records+(replace(records[0],equality_group='Two'),)))
    with pytest.raises(ValueError,match='Duplicate projected elastic hkl/slice identity'):
        compile_config_from_sources(c)
