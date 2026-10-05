"""Schema/persistence equivalence only: fixed evaluations, no optimization."""
from dataclasses import asdict,replace
import json
from pathlib import Path
import numpy as np
import pytest
import yaml
from qcbed import load_yaml_config,save_yaml_config,prepare_run_snapshot,compile_config_from_sources
from qcbed.config_serialization import _plain
from qcbed.grouped_run import grouped_identity,evaluate_grouped
from qcbed.modern_yaml import target_name
from test_modern_grouped_run import modern_config


def save(c,p):
    save_yaml_config(c,p,base_directory=p.parent)
    return yaml.safe_load(p.read_text())


def test_abc_roundtrip_compiled_fixed_numerical_equivalence(tmp_path):
    c=modern_config(tmp_path);p=tmp_path/'user.yml';doc=save(c,p);loaded=load_yaml_config(p)
    assert loaded==c
    assert set(doc)=={'schema_version','data_prep','atomic_model','qcbed'}
    assert 'output_filename' not in doc['qcbed']['output']
    before=prepare_run_snapshot(c,base_directory=tmp_path)
    after=prepare_run_snapshot(loaded,base_directory=tmp_path)
    assert grouped_identity(before)==grouped_identity(after)
    assert before.compiled.state==after.compiled.state
    assert before.compiled.layout==after.compiled.layout
    assert before.compiled.bounds==after.compiled.bounds
    first,arrays=evaluate_grouped(before,[4.,.86])
    second,other=evaluate_grouped(after,[4.,.86])
    assert first==second
    for k,v in arrays.items():np.testing.assert_array_equal(v,other[k])
    for a,b in zip(before.compiled.crystallography.simulation_sample.prepared.slice_types,
                   after.compiled.crystallography.simulation_sample.prepared.slice_types):
        np.testing.assert_array_equal(a.phase_grating,b.phase_grating)
    assert not Path(c.run.output_directory).exists()
    for forbidden in ('internal_hk','effective_depth','occurrence_count','ynorm','checkpoint_identity','slice_count'):
        assert forbidden not in p.read_text().lower()


def test_verbose_compatibility_and_compact_default_output(tmp_path):
    c=modern_config(tmp_path)
    c=replace(c,run=replace(c.run,output_filename=str(Path(c.run.output_directory)/'final_out.dat')))
    p=tmp_path/'old.yml';p.write_text(yaml.safe_dump(_plain(asdict(c)),sort_keys=False))
    assert load_yaml_config(p)==c
    doc=save(c,tmp_path/'compact.yml')
    assert 'filename' not in doc['qcbed']['output']
    assert load_yaml_config(tmp_path/'compact.yml')==c
    # Versioned defaults are actual codec defaults, not case-specific inference.
    from qcbed.modern_yaml import OPTIMIZER_DEFAULTS
    c=replace(c,refinement=replace(c.refinement,optimizer=replace(c.refinement.optimizer,**OPTIMIZER_DEFAULTS)))
    doc=save(c,tmp_path/'minimal.yml')
    assert 'optimizer' not in doc['qcbed']
    assert load_yaml_config(tmp_path/'minimal.yml')==c


def test_mask_origin_provenance_save_as(tmp_path):
    from qcbed.gui.data_model import attach_overview,edit_mask,select_zone_axis,provenance,read_data
    from test_gui_data_model import overview_files
    c=modern_config(tmp_path);pattern,manifest=overview_files(c,tmp_path)
    c=attach_overview(c,tmp_path,pattern,manifest)
    c=select_zone_axis(c,tmp_path,beam=0)
    c=edit_mask(c,tmp_path,beam=1,box=(0,0,1,1))
    baseline=read_data(c,tmp_path)
    doc=save(c,tmp_path/'a.yml');descriptor=Path(doc['data_prep']['prepared_dataset'])
    stored=json.loads(descriptor.read_text());assert 'preparation_metadata' in stored
    sub=tmp_path/'elsewhere';sub.mkdir();save_yaml_config(c,sub/'copy.yml',base_directory=tmp_path)
    loaded=load_yaml_config(sub/'copy.yml');assert loaded==c
    assert provenance(loaded)==provenance(c)
    assert provenance(loaded)['zone_axis_selection']['zone_axis_origin_px']==[3.25,4.75]
    np.testing.assert_array_equal(read_data(loaded,sub).weights,baseline.weights)
    assert read_data(loaded,sub).weights[0,1]==0


def test_dynamic_group_identity_and_target_order(tmp_path):
    from qcbed.crystallography.site_groups import group_target
    c=modern_config(tmp_path);old=group_target('Al_fcc');new=group_target('user.group')
    c=replace(c,crystallography=replace(c.crystallography,site_groups=(replace(c.crystallography.site_groups[0],id='user.group'),)),
        refinement=replace(c.refinement,parameters=tuple(replace(b,target=new) if b.target==old else b for b in reversed(c.refinement.parameters)),
        bounds=replace(c.refinement.bounds,specific=tuple(replace(b,target=new) if b.target==old else b for b in c.refinement.bounds.specific))))
    doc=save(c,tmp_path/'dynamic.yml');assert 'site_group.user.group.debye_waller_b' in doc['qcbed']['refinement']['targets']
    assert load_yaml_config(tmp_path/'dynamic.yml')==c
    empty=replace(c,refinement=replace(c.refinement,parameters=()))
    save(empty,tmp_path/'empty.yml');assert load_yaml_config(tmp_path/'empty.yml')==empty


@pytest.mark.parametrize('mutate',[
    lambda d:d['qcbed'].update(effective_depth=5),
    lambda d:d['atomic_model'].update(internal_hk=[0,0]),
    lambda d:d['qcbed'].update(signal_mode='unknown'),
    lambda d:d['qcbed']['optimizer'].update(cycles=True),
    lambda d:d['qcbed'].update(accelerating_voltage_kv=float('nan')),
    lambda d:d['qcbed']['refinement']['targets'].update(unknown={'range_factors':[1,1],'plamda':.1}),
    lambda d:d['qcbed'].update(workflow=['AMOEBA','ANNEAL']),
])
def test_strict_compact_validation(tmp_path,mutate):
    p=tmp_path/'config.yml';d=save(modern_config(tmp_path),p);mutate(d)
    p.write_text(yaml.safe_dump(d))
    with pytest.raises((ValueError,TypeError)):load_yaml_config(p)


def test_descriptor_metadata_tampering_and_duplicate_keys(tmp_path):
    from qcbed.gui.data_model import attach_overview
    from test_gui_data_model import overview_files
    c=modern_config(tmp_path);pattern,manifest=overview_files(c,tmp_path);c=attach_overview(c,tmp_path,pattern,manifest)
    p=tmp_path/'config.yml';d=save(c,p)
    side=Path(c.experiment.path+'.prep.json');side.write_text(side.read_text()+' ')
    with pytest.raises(ValueError,match='checksum'):load_yaml_config(p)
    p.write_text('schema_version: 2\nschema_version: 2\n')
    with pytest.raises(ValueError,match='duplicate YAML'):load_yaml_config(p)


def test_shipped_raw_fcc_al_template_compiles_after_preparation(tmp_path):
    """The teaching YAML must decode after A supplies its real prepared descriptor.

    The small synthetic descriptor exercises configuration only, not a reference fit.
    The packaged acceptance uses the actual raw pattern for end-to-end preparation.
    """
    from qcbed.modern_yaml import decode_modern, publish_prepared_descriptor
    root = Path(__file__).resolve().parents[1]
    project = root / 'examples/fcc_al.yml'
    raw = yaml.safe_load(project.read_text())
    assert (project.parent / raw['data_prep']['pattern']).is_file()
    assert (project.parent / raw['atomic_model']['structure']['path']).is_file()
    raw['data_prep'] = {'prepared_dataset': str(publish_prepared_descriptor(modern_config(tmp_path), tmp_path))}
    config = decode_modern(raw, project.parent)
    assert config.crystallography.depth_angstrom == 1000
    assert config.crystallography.slice_boundaries_fractional == (0., 1.)
    assert config.crystallography.site_groups[0].debye_waller_b == .86
    assert config.crystallography.site_groups[0].occupancy == 1.
    assert config.execution.execution_backend == 'batched'
    assert config.execution.batch_size == 64
    assert config.run.max_total_evaluations == 2
    assert len(config.refinement.parameters) == 2
    saved = tmp_path / 'prepared.yml'
    save_yaml_config(config, saved, base_directory=project.parent)
    assert load_yaml_config(saved) == config
