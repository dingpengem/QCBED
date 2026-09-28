"""Atomic YAML saving and identical canonical/YAML snapshot execution."""
from dataclasses import replace, FrozenInstanceError, asdict
import hashlib
import json
from pathlib import Path

import numpy as np
import pytest

from qcbed import (compile_config, save_yaml_config, load_yaml_config, resolve_config_paths,
    prepare_run_snapshot, run_snapshot, run_config, validate_config, ValidationPhase as Phase)
from qcbed.canonical_run import canonical_checkpoint_identity
from qcbed.run import run_yaml_refinement
from qcbed.full_runner import load_refinement_checkpoint
from test_yaml_runner import setup
from test_full_runner import cheap_render
from test_multiregion_configuration import configuration
from test_configuration import example


@pytest.mark.parametrize('reference',[False,True])
def test_exact_save_load_and_explicit_paths(tmp_path,monkeypatch,reference):
    raw = example('reference_configuration').reference_configuration() if reference else configuration(((2.,),(1.,3.),(4.,)),(4.,7.,8.),(0,2))
    raw = replace(raw,microscope=replace(raw.microscope,energy_kv=200.12345678901235))
    base = tmp_path/'document'; base.mkdir()
    # Path behavior cannot change with caller CWD.
    monkeypatch.chdir(tmp_path)
    resolved = resolve_config_paths(raw,base_directory=base)
    result = save_yaml_config(raw,'saved.yml',base_directory=base)
    assert result == resolved == load_yaml_config(base/'saved.yml')
    assert compile_config(result) == compile_config(resolved)
    assert result.experiment.path == str((base/raw.experiment.path).resolve())
    assert result.run.output_filename == str((base/raw.run.output_filename).resolve())
    assert not (tmp_path/'saved.yml').exists()
    # An already resolved canonical snapshot round-trips with complete equality.
    save_yaml_config(result,'again.yml',base_directory=base)
    assert load_yaml_config(base/'again.yml') == result
    with pytest.raises(ValueError,match='absolute'):
        save_yaml_config(raw,'bad.yml',base_directory='relative')


def test_save_as_and_atomic_replacement_failure(tmp_path,monkeypatch):
    raw = configuration(((2.,),),(4.,))
    path = tmp_path/'config.yml'
    original = save_yaml_config(raw,path,base_directory=tmp_path)
    content = path.read_bytes()
    changed = replace(original,title='changed title')
    with pytest.raises(FileExistsError): save_yaml_config(changed,path,base_directory=tmp_path)
    assert path.read_bytes() == content
    import qcbed.config_serialization as serializer
    original_replace = serializer.os.replace
    def fail(*args): raise OSError('publication failed')
    monkeypatch.setattr(serializer.os,'replace',fail)
    with pytest.raises(OSError,match='publication failed'):
        save_yaml_config(changed,path,base_directory=tmp_path,overwrite=True)
    assert path.read_bytes() == content
    assert not list(tmp_path.glob('.config.yml.*'))
    monkeypatch.setattr(serializer.os,'replace',original_replace)
    save_yaml_config(changed,path,base_directory=tmp_path,overwrite=True)
    assert load_yaml_config(path) == changed
    assert not list(tmp_path.glob('.config.yml.*'))


def test_save_never_overwrites_experimental_or_reference_files(tmp_path):
    raw,path = setup(tmp_path)
    data = Path(raw.experiment.path); before = data.read_bytes()
    with pytest.raises(ValueError,match='experimental'):
        save_yaml_config(raw,data,base_directory=tmp_path,overwrite=True)
    assert data.read_bytes() == before
    root = Path(__file__).resolve().parents[1]
    with pytest.raises(ValueError,match='protected'):
        save_yaml_config(raw,root/'reference_data/new.yml',base_directory=tmp_path)


def test_file_readiness_and_export_diagnostics(tmp_path):
    raw,path = setup(tmp_path)
    assert validate_config(raw,phases=(Phase.FILES,),base_directory=tmp_path) == ()
    assert validate_config(raw,phases=(Phase.FILES,))[0].code == 'base_directory_required'
    bad = replace(raw,experiment=replace(raw.experiment,path='missing.dat'))
    issues = validate_config(bad,phases=(Phase.FILES,),base_directory=tmp_path)
    assert [(i.code,i.section,i.field) for i in issues] == [('experimental_data_error','experiment','path')]
    (tmp_path/'run').mkdir()
    issues = validate_config(raw,phases=(Phase.FILES,),base_directory=tmp_path)
    assert any(i.code == 'output_exists' and i.field == 'output_directory' for i in issues)
    issues = validate_config(raw,phases=(Phase.FILES,),base_directory=tmp_path,resume=True)
    assert any(i.code == 'resume_incompatible' for i in issues)


def test_canonical_to_yaml_resume_and_old_identity_formula(tmp_path,cheap_render):
    raw,path = setup(tmp_path)
    snapshot = prepare_run_snapshot(raw,base_directory=tmp_path)
    compiled = snapshot.compiled
    old_identity = dict(format='qcbed-canonical-run',version=1,digest=hashlib.sha256(json.dumps(dict(
        state=asdict(compiled.state),runner=asdict(compiled.runner),layout=asdict(compiled.layout),
        experiment=asdict(compiled.experiment),run=asdict(compiled.configuration.run)),
        sort_keys=True,allow_nan=False).encode()).hexdigest())
    assert canonical_checkpoint_identity(compiled) == old_identity
    def stop(row):
        if row['classification'] == 'stage_retained': raise KeyboardInterrupt
    interrupted = run_snapshot(snapshot,progress_hook=stop)
    assert interrupted.status == 'interrupted'
    assert load_refinement_checkpoint(tmp_path/'run')['checkpoint_identity'] == old_identity
    assert not validate_config(raw,phases=(Phase.FILES,),base_directory=tmp_path,resume=True)
    assert run_yaml_refinement(path,resume=True).status == 'completed'
    assert run_config(raw,base_directory=tmp_path,resume=True).status == 'completed'
    direct = replace(raw,run=replace(raw.run,output_directory=str(tmp_path/'direct'),output_filename=str(tmp_path/'direct.dat')))
    assert run_config(direct,base_directory=tmp_path).status == 'completed'
    np.testing.assert_array_equal(np.load(tmp_path/'run/final_active_vector.npy'),np.load(tmp_path/'direct/final_active_vector.npy'))
    np.testing.assert_array_equal(np.load(tmp_path/'run/final_simulated_intensity.npy'),np.load(tmp_path/'direct/final_simulated_intensity.npy'))
    one = json.loads((tmp_path/'run/summary.json').read_text())
    two = json.loads((tmp_path/'direct/summary.json').read_text())
    assert one['final_objective'] == two['final_objective'] and one['final_ynorm'] == two['final_ynorm']


def test_snapshot_isolation_and_no_yaml_dependency(tmp_path,cheap_render,monkeypatch):
    raw,path = setup(tmp_path,1)
    snapshot = prepare_run_snapshot(raw,base_directory=tmp_path)
    before = snapshot.experimental.intensities.copy()
    edited = replace(raw,sample=replace(raw.sample,regions=(replace(raw.sample.regions[0],depth_angstrom=10.),)))
    Path(raw.experiment.path).write_text('2 1\n1 1\n'+'1 1\n1 1\n'*2)
    assert snapshot.compiled.state.regions[0].depth_angstrom == 4.
    np.testing.assert_array_equal(snapshot.experimental.intensities,before)
    with pytest.raises(ValueError): snapshot.experimental.intensities.setflags(write=True)
    with pytest.raises(FrozenInstanceError): snapshot.compiled = compile_config(edited)
    import sys
    monkeypatch.setitem(sys.modules,'yaml',None)
    assert run_snapshot(snapshot).status == 'completed'
    with pytest.raises(ValueError,match='fingerprint'):
        run_config(raw,base_directory=tmp_path,resume=True)
    with pytest.raises(ValueError,match='identity'):
        run_config(edited,base_directory=tmp_path,resume=True)
    assert any(i.code == 'resume_incompatible' for i in validate_config(raw,phases=(Phase.FILES,),base_directory=tmp_path,resume=True))
