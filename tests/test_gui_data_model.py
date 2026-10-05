from dataclasses import replace
import hashlib,json
from pathlib import Path
import numpy as np
import pytest
from qcbed.gui.data_model import read_data,edit_mask,provenance,attach_overview,restored_overview
from qcbed import save_yaml_config,load_yaml_config,prepare_run_snapshot,run_snapshot
from qcbed.gui.readiness import readiness
from test_modern_grouped_run import modern_config


def test_mask_roundtrip_restore_resume_and_exact_data(tmp_path):
    config=modern_config(tmp_path);old=read_data(config,tmp_path);source=Path(config.experiment.path).read_bytes()
    masked=edit_mask(config,tmp_path,beam=1,box=(0,0,2,1));data=read_data(masked,tmp_path)
    np.testing.assert_array_equal(data.intensities,old.intensities)
    np.testing.assert_array_equal(data.variances,old.variances)
    assert data.weights[:,1].tolist()==[0,0,1,1]
    assert Path(config.experiment.path).read_bytes()==source
    save_yaml_config(masked,tmp_path/'saved.yml',base_directory=tmp_path)
    loaded=load_yaml_config(tmp_path/'saved.yml');assert loaded==masked
    sub=tmp_path/'save-as';sub.mkdir();save_yaml_config(loaded,sub/'copy.yml',base_directory=tmp_path)
    assert load_yaml_config(sub/'copy.yml')==masked
    snapshot=prepare_run_snapshot(loaded,base_directory=tmp_path)
    events=[];run_snapshot(snapshot,progress_hook=events.append,cancel_requested=lambda:bool(events))
    assert not readiness(loaded,tmp_path,resume=True)
    assert any(i.code=='resume_incompatible' for i in readiness(config,tmp_path,resume=True))
    restored=edit_mask(loaded,tmp_path,reset=True)
    np.testing.assert_array_equal(read_data(restored,tmp_path).weights,old.weights)
    # Restoring a source-invalid point must never enable it.
    restored=edit_mask(loaded,tmp_path,beam=0,box=(0,0,1,1),restore=True)
    assert read_data(restored,tmp_path).weights[0,0]==0


def overview_files(config,tmp_path):
    data=read_data(config,tmp_path)
    pattern=tmp_path/'pattern.dat';np.savetxt(pattern,np.arange(200).reshape(10,20))
    manifest=dict(signal_mode='gradient',arrays=dict(intensities=dict(sha256=hashlib.sha256(data.intensities.tobytes()).hexdigest())),
        beam_mapping=[dict(refinement_column=i,stable_beam_id=f'b{i}',center_px=[i+1,4]) for i in range(9)],
        rois=[dict(beam_id=f'b{i}',center_px=[i+1,4],geometry=dict(radius_x_px=.4,radius_y_px=.8)) for i in range(9)],
        reciprocal_frame=dict(zone_axis_origin_px=[3.25,4.75]))
    path=tmp_path/'overview.json';path.write_text(json.dumps(manifest))
    return pattern,path


def test_overview_metadata_binding_and_tamper_rejection(tmp_path):
    c=modern_config(tmp_path);pattern,manifest=overview_files(c,tmp_path)
    c=attach_overview(c,tmp_path,pattern,manifest)
    image,m=restored_overview(c,tmp_path)
    assert image.shape==(10,20) and m['reciprocal_frame']['zone_axis_origin_px']==[3.25,4.75]
    c=edit_mask(c,tmp_path,beam=0,box=(0,0,1,1));assert restored_overview(c,tmp_path)[0].shape==(10,20)
    manifest.write_text('{}')
    with pytest.raises(ValueError,match='source changed'):restored_overview(c,tmp_path)


def test_invalid_mask_coordinates_do_not_write(tmp_path):
    c=modern_config(tmp_path)
    with pytest.raises(ValueError):edit_mask(c,tmp_path,beam=0,box=(0,0,3,1))
    assert not (tmp_path/'.qcbed-data').exists()


def test_legacy_duplicate_beams_and_header_are_preserved(tmp_path):
    from examples.reference_configuration import reference_configuration
    from qcbed.config_paths import resolve_config_paths
    c=reference_configuration(legacy_header=True)
    root=Path(__file__).resolve().parents[1]
    c=replace(c,experiment=replace(c.experiment,path=str(root/"reference_data/in.dat")),legacy_export=replace(c.legacy_export,experimental_filename_raw=None))
    c=resolve_config_paths(c,base_directory=root)
    c=replace(c,experiment=replace(c.experiment,beam_indices_hk=c.experiment.beam_indices_hk[:-1]+c.experiment.beam_indices_hk[:1]))
    data=read_data(c,root)
    assert len(set(data.beam_indices_hk))<len(data.beam_indices_hk)
    masked=edit_mask(c,tmp_path,beam=len(data.beam_indices_hk)-1,box=(0,0,1,1))
    changed=read_data(masked,tmp_path)
    assert changed.beam_indices_hk==data.beam_indices_hk
    np.testing.assert_array_equal(changed.intensities,data.intensities)
    np.testing.assert_array_equal(changed.variances,data.variances)
    assert masked.experiment.header=='absent'
