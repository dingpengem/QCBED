"""Real spawned canonical runs; no Qt required."""
from dataclasses import replace
import os
from pathlib import Path
import time

import numpy as np
import pytest
from qcbed import prepare_run_snapshot, run_config
from qcbed.gui.run_controller import RunController, readiness, _worker
from test_dynamic_configuration import small_config


def run_config_fixture(tmp_path):
    raw = small_config()
    body = tmp_path/'body.dat'
    body.write_text('1 1\n1 1\n'*3)
    return replace(raw, experiment=replace(raw.experiment, path=str(body)),
        run=replace(raw.run, output_directory=str(tmp_path/'run'), output_filename=str(tmp_path/'out.dat')))


def wait(controller):
    events = []
    deadline = time.monotonic()+90
    while controller.active:
        events.extend(controller.poll())
        assert time.monotonic() < deadline, 'worker timeout'
        time.sleep(.01)
    return events


def test_spawn_run_snapshot_progress_resume_and_equivalence(tmp_path):
    config = run_config_fixture(tmp_path)
    controller = RunController()
    controller.start(config, base_directory=tmp_path)
    assert controller.state == 'Running'
    changed = replace(config, title='future snapshot')
    assert controller.snapshot.compiled.configuration.title != changed.title
    body = Path(config.experiment.path)
    original_data = body.read_text()
    body.write_text('999 999\n999 999\n'*3)
    events = wait(controller)
    assert controller.state == 'Completed'
    body.write_text(original_data)
    assert next(p for k,p in events if k == 'started')['pid'] != os.getpid()
    assert any(k == 'progress' for k,p in events)
    assert controller.result['checkpoint']
    assert controller.result['final_state']
    assert any(i.code == 'output_exists' for i in readiness(config,tmp_path))
    with pytest.raises(ValueError, match='new output directory'):
        controller.start(config,base_directory=tmp_path)
    assert not [i for i in readiness(config,tmp_path,resume=True) if i.severity == 'error']
    controller.start(config,base_directory=tmp_path,resume=True)
    wait(controller)
    assert controller.state == 'Completed'
    incompatible = replace(config, run=replace(config.run, checkpoint_every=config.run.checkpoint_every+1))
    assert any(i.code == 'resume_incompatible' for i in readiness(incompatible,tmp_path,resume=True))
    with pytest.raises(ValueError): controller.start(incompatible,base_directory=tmp_path,resume=True)
    direct = replace(config, run=replace(config.run,output_directory=str(tmp_path/'direct'),output_filename=str(tmp_path/'direct.dat')))
    result = run_config(direct,base_directory=tmp_path)
    assert result.summary['final_objective'] == controller.result['summary']['final_objective']
    np.testing.assert_array_equal(np.load(tmp_path/'run/final_active_vector.npy'),np.load(tmp_path/'direct/final_active_vector.npy'))
    assert (tmp_path/'out.dat').read_bytes() == (tmp_path/'direct.dat').read_bytes()


def test_worker_exception_and_interrupted_delivery(tmp_path, monkeypatch):
    import qcbed.gui.run_controller as module
    snapshot = prepare_run_snapshot(run_config_fixture(tmp_path),base_directory=tmp_path)
    class Connection:
        def __init__(self): self.events = []
        def send(self,event): self.events.append(event)
        def close(self): pass
    connection = Connection()
    def fail(*args, **kwargs): raise RuntimeError('test worker failure')
    monkeypatch.setattr(module,'run_snapshot',fail)
    _worker(snapshot,False,connection)
    failure = next(p for k,p in connection.events if k == 'failure')
    assert 'Traceback' in failure and 'test worker failure' in failure


def test_process_failure_state(tmp_path):
    # Existing output introduced after launch exercises runner protection in worker.
    # Use a canonical snapshot and worker directly so the race is deterministic.
    import multiprocessing
    controller = RunController()
    snapshot = prepare_run_snapshot(run_config_fixture(tmp_path),base_directory=tmp_path)
    Path(snapshot.compiled.configuration.run.output_directory).mkdir()
    context = multiprocessing.get_context('spawn')
    receiver,sender = context.Pipe(False)
    process = context.Process(target=_worker,args=(snapshot,False,sender))
    process.start()
    sender.close()
    controller.process,controller.connection = process,receiver
    controller.state = 'Running'
    events = wait(controller)
    assert controller.state == 'Failed'
    assert any(k == 'failure' and 'Traceback' in p for k,p in events)


def test_budget_interrupted_state(tmp_path):
    from qcbed.configuration import ParameterBinding, ParameterTarget, ParameterCategory
    config = run_config_fixture(tmp_path)
    target = ParameterTarget(ParameterCategory.DEPTH,(config.sample.regions[0].id,))
    config = replace(config, refinement=replace(config.refinement,
        parameters=(ParameterBinding(target,(1.1,.9),.1),),
        optimizer=replace(config.refinement.optimizer,loops=1)),
        run=replace(config.run,max_total_evaluations=2))
    controller = RunController()
    controller.start(config,base_directory=tmp_path)
    wait(controller)
    assert controller.state == 'Interrupted'
    assert controller.result['status'] == 'limited'
    assert controller.result['checkpoint']
    assert not [i for i in readiness(config,tmp_path,resume=True) if i.severity == 'error']


@pytest.mark.parametrize('existing', [False, True])
def test_missing_checkpoint_does_not_attempt_resume_load(tmp_path, monkeypatch, existing):
    import qcbed.gui.readiness as module
    config = run_config_fixture(tmp_path)
    config = replace(config, run=replace(config.run, output_directory='relative-run'))
    if existing:
        (tmp_path / 'relative-run').mkdir()
    original = module.validate_config
    def validate(*args, **kwargs):
        assert not kwargs.get('resume', False), 'Missing checkpoint must not be loaded'
        return original(*args, **kwargs)
    monkeypatch.setattr(module, 'validate_config', validate)
    issues = readiness(config, tmp_path, resume=True)
    assert any(i.code == 'resume_checkpoint_missing' for i in issues)
    assert not (tmp_path / 'relative-run/checkpoint_latest.json').exists()
