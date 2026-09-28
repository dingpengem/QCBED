"""Production adapter tests use real refinement and durable output operations."""
from dataclasses import asdict, replace
import json
from pathlib import Path

import numpy as np
import pytest

from qcbed.run import run_yaml_refinement, main
from qcbed import compile_config
from qcbed.full_runner import run_full_refinement, load_refinement_checkpoint
from qcbed.configuration import BoundsConfig, TargetBound, ParameterTarget, ParameterCategory as C
from test_yaml_config import plain, write
from test_multiregion_configuration import configuration
from test_dynamic_configuration import experimental
from test_full_runner import cheap_render
from test_configuration import example, legacy_full
from qcbed.experimental import parse_experimental_intensities


def setup(tmp_path, count=2, active=None):
    if active is None: active = tuple(range(count))
    raw = configuration(tuple((1.,2.) for _ in range(count)),tuple(4.+i for i in range(count)),active)
    body = tmp_path/'body.dat'
    body.write_text('1 1\n1 1\n'*3)
    raw = replace(raw,experiment=replace(raw.experiment,path=str(body)),
        run=replace(raw.run,output_directory=str(tmp_path/'run'),output_filename=str(tmp_path/'output.dat'),
                    snapshot_every=1,checkpoint_every=1,max_total_evaluations=150),
        refinement=replace(raw.refinement,optimizer=replace(raw.refinement.optimizer,loops=2)))
    return raw,write(tmp_path,plain(asdict(raw)))


@pytest.mark.parametrize('count,active',[(1,(0,)),(2,(1,)),(3,(0,1,2)),(5,(0,2,4))])
def test_yaml_only_fresh_and_outputs(tmp_path,cheap_render,monkeypatch,count,active):
    raw,path = setup(tmp_path,count,active)
    # Fail immediately if any legacy FOR007 entry point is used.
    import qcbed.legacy_io as legacy
    def forbidden(*args,**kwargs): raise AssertionError('FOR007 dependency')
    monkeypatch.setattr(legacy,'parse_for007',forbidden)
    monkeypatch.setattr(legacy,'_parse_for007_records',forbidden)
    assert not list(tmp_path.glob('*FOR007*'))
    result = run_yaml_refinement(path)
    assert result.status == 'completed'
    folder = tmp_path/'run'
    cp = load_refinement_checkpoint(folder)
    assert cp['checkpoint_identity']['format'] == 'qcbed-canonical-run'
    assert cp['input_files'] == {}
    assert len(cp['retained_vector']) == len(active)
    for name in ('checkpoint_latest.json','initial_snapshot.npz','final_snapshot.npz',
                 'final_out.dat','final_parameters.json','summary.json','report.html','progress.jsonl'):
        assert (folder/name).exists()
    assert (tmp_path/'output.dat').read_bytes() == (folder/'final_out.dat').read_bytes()
    compiled = compile_config(raw)
    direct = run_full_refinement(compiled.state,experimental(),compiled.runner,output_directory=tmp_path/'direct')
    assert direct.summary['final_objective'] == result.summary['final_objective']
    np.testing.assert_array_equal(np.load(folder/'final_active_vector.npy'),np.load(tmp_path/'direct/final_active_vector.npy'))
    assert any(json.loads(line)['classification'] == 'snapshot_verification' for line in (folder/'progress.jsonl').read_text().splitlines())


def test_interrupted_resume_formatting_and_identity(tmp_path,cheap_render):
    raw,path = setup(tmp_path)
    def stop(row):
        if row['classification'] == 'stage_retained': raise KeyboardInterrupt
    first = run_yaml_refinement(path,progress_hook=stop)
    assert first.status == 'interrupted'
    # Mapping key order and comments do not affect identity.
    import yaml
    path.write_text('# reformatted\n'+yaml.safe_dump(plain(asdict(raw)),sort_keys=True))
    result = run_yaml_refinement(path,resume=True)
    assert result.status == 'completed'
    compiled = compile_config(raw)
    run_full_refinement(compiled.state,experimental(),compiled.runner,output_directory=tmp_path/'direct')
    np.testing.assert_array_equal(np.load(tmp_path/'run/final_active_vector.npy'),np.load(tmp_path/'direct/final_active_vector.npy'))
    assert json.loads((tmp_path/'run/final_parameters.json').read_text()) == json.loads((tmp_path/'direct/final_parameters.json').read_text())
    cp = load_refinement_checkpoint(tmp_path/'run')
    with pytest.raises(ValueError,match='identity mismatch'):
        run_full_refinement(compiled.state,experimental(),compiled.runner,output_directory=tmp_path/'run',resume=True)
    checkpoint = tmp_path/'run/checkpoint_latest.json'
    old = dict(cp)
    old.pop('checkpoint_identity')
    checkpoint.write_text(json.dumps(old))
    with pytest.raises(ValueError,match='identity mismatch'):
        run_yaml_refinement(path,resume=True)
    checkpoint.write_text(json.dumps(cp))
    changed = replace(raw,run=replace(raw.run,checkpoint_every=2))
    write(tmp_path,plain(asdict(changed)))
    with pytest.raises(ValueError,match='identity mismatch'):
        run_yaml_refinement(path,resume=True)
    write(tmp_path,plain(asdict(raw)))
    (tmp_path/'body.dat').write_text('2 1\n1 1\n'+'1 1\n1 1\n'*2)
    with pytest.raises(ValueError,match='fingerprint mismatch'):
        run_yaml_refinement(path,resume=True)


def test_yaml_bound_rejections(tmp_path,cheap_render):
    raw,path = setup(tmp_path,1)
    raw = replace(raw,refinement=replace(raw.refinement,bounds=BoundsConfig('reject',specific=(
        TargetBound(ParameterTarget(C.DEPTH,('region-0',)),4.,4.),))))
    write(tmp_path,plain(asdict(raw)))
    result = run_yaml_refinement(path)
    assert result.summary['bound_rejection_count'] > 0
    np.testing.assert_array_equal(np.load(tmp_path/'run/final_active_vector.npy'),[4.])


def test_legacy_reference_vs_yaml(tmp_path,cheap_render):
    root = Path(__file__).resolve().parents[1]
    state,runner = legacy_full()
    raw = example('reference_configuration').reference_configuration(legacy_header=True)
    # Explicit test-sized reference workload, identical on both entry paths.
    raw = replace(raw,sample=replace(raw.sample,regions=(replace(raw.sample.regions[0],depth_angstrom=4.),)),
        experiment=replace(raw.experiment,path=str(root/'reference_data/in.dat')),
        run=replace(raw.run,output_directory=str(tmp_path/'yaml'),output_filename=str(tmp_path/'output.dat')),
        legacy_export=replace(raw.legacy_export,experimental_filename_raw=None,output_filename_raw=None),
        refinement=replace(raw.refinement,bounds=BoundsConfig(),optimizer=replace(raw.refinement.optimizer,loops=0)))
    state = replace(state,regions=(replace(state.regions[0],depth_angstrom=4.),))
    c = compile_config(raw)
    # Only filename metadata differs from the parsed legacy state.
    assert replace(c.state,model=replace(c.state.model,
        output_filename=state.model.output_filename,output_filename_raw=state.model.output_filename_raw,
        experimental_intensity_filename=state.model.experimental_intensity_filename,
        experimental_intensity_filename_raw=state.model.experimental_intensity_filename_raw)) == state
    data = parse_experimental_intensities(root/'reference_data/in.dat',beam_indices_hk=raw.experiment.beam_indices_hk,scan=state.scan)
    legacy = run_full_refinement(state,data,c.runner,output_directory=tmp_path/'legacy')
    canonical = run_yaml_refinement(write(tmp_path,plain(asdict(raw))))
    assert canonical.summary['final_objective'] == legacy.summary['final_objective']
    assert canonical.summary['final_ynorm'] == legacy.summary['final_ynorm']
    np.testing.assert_array_equal(np.load(tmp_path/'yaml/final_simulated_intensity.npy'),np.load(tmp_path/'legacy/final_simulated_intensity.npy'))


def test_cli(tmp_path,cheap_render,capsys):
    _,path = setup(tmp_path,1)
    assert main(['--config',str(path)]) == 0
    assert 'completed' in capsys.readouterr().out
    assert main(['--config',str(path),'--resume']) == 0
