"""Cancellation delegates to canonical retained-boundary interruption/resume."""
from dataclasses import replace
from pathlib import Path
from threading import Event
import numpy as np
import pytest
from qcbed import prepare_run_snapshot, run_config
from qcbed.full_runner import load_refinement_checkpoint
from qcbed.gui.run_controller import _worker, RunController, readiness
from test_gui_run import run_config_fixture, wait


class Connection:
    def __init__(self, event, boundary=None, directory=None):
        self.event, self.boundary, self.directory = event,boundary,directory
        self.events = []
        self.checkpoints = {}
    def send(self, message):
        self.events.append(message)
        kind,payload = message
        if kind == 'progress' and payload['classification'] == self.boundary and not self.event.is_set():
            self.checkpoints = {p:p.read_bytes() for p in self.directory.glob('checkpoint_*.json') if p.name != 'checkpoint_latest.json'}
            self.event.set()
    def close(self): pass


def test_cancel_before_execution_has_no_fabricated_checkpoint(tmp_path):
    config = run_config_fixture(tmp_path)
    event = Event(); event.set()
    connection = Connection(event)
    _worker(prepare_run_snapshot(config,base_directory=tmp_path),False,connection,event)
    result = next(p for k,p in connection.events if k == 'result')
    assert result['status'] == 'cancelled'
    assert result['summary']['evaluations'] == 0
    assert not result['checkpoint']
    assert not Path(config.run.output_directory).exists()
    assert not any(k == 'progress' for k,p in connection.events)


@pytest.mark.parametrize('boundary',['initial_retained','trial','retained_boundary','stage_retained'])
def test_cancel_safe_boundaries_preserves_checkpoint_and_resume(tmp_path,boundary):
    from qcbed.configuration import ParameterBinding,ParameterTarget,ParameterCategory
    config = run_config_fixture(tmp_path)
    target = ParameterTarget(ParameterCategory.DEPTH,(config.sample.regions[0].id,))
    config = replace(config, refinement=replace(config.refinement,
        parameters=(ParameterBinding(target,(1.1,.9),.1),),
        optimizer=replace(config.refinement.optimizer,loops=1)))
    event = Event()
    directory = Path(config.run.output_directory)
    connection = Connection(event,boundary,directory)
    snapshot = prepare_run_snapshot(config,base_directory=tmp_path)
    _worker(snapshot,False,connection,event)
    result = next(p for k,p in connection.events if k == 'result')
    assert result['status'] == 'cancelled'
    assert result['cancellation']['boundary'] == boundary
    assert result['checkpoint']
    assert len([k for k,p in connection.events if k == 'cancellation_effective']) == 1
    cp = load_refinement_checkpoint(directory)
    assert cp['retained_objective'] is not None
    assert not (directory/'.run.lock').exists()
    for path,content in connection.checkpoints.items(): assert path.read_bytes() == content
    assert not [i for i in readiness(config,tmp_path,resume=True) if i.severity == 'error']
    resumed = run_config(config,base_directory=tmp_path,resume=True)
    direct = replace(config,run=replace(config.run,output_directory=str(tmp_path/'direct'),output_filename=str(tmp_path/'direct.dat')))
    uninterrupted = run_config(direct,base_directory=tmp_path)
    assert resumed.status == 'completed'
    assert resumed.summary['final_objective'] == uninterrupted.summary['final_objective']
    np.testing.assert_array_equal(np.load(directory/'final_active_vector.npy'),np.load(tmp_path/'direct/final_active_vector.npy'))
    # Cancelling a resume before starting must not rewrite the existing checkpoint.
    before = (directory/'checkpoint_latest.json').read_bytes()
    early = Connection(event)
    _worker(snapshot,True,early,event)
    assert (directory/'checkpoint_latest.json').read_bytes() == before


def test_cancel_pending_failure_is_failed_not_cancelled(tmp_path,monkeypatch):
    import qcbed.gui.run_controller as module
    event = Event()
    connection = Connection(event)
    def fail(*args,**kwargs):
        event.set()
        raise RuntimeError('failure with cancellation pending')
    monkeypatch.setattr(module,'run_snapshot',fail)
    _worker(prepare_run_snapshot(run_config_fixture(tmp_path),base_directory=tmp_path),False,connection,event)
    assert any(k == 'failure' and 'failure with cancellation pending' in p for k,p in connection.events)
    assert not any(k == 'result' for k,p in connection.events)


def test_repeated_cancel_and_spawn_state(tmp_path):
    controller = RunController()
    assert not controller.request_cancel()
    controller.start(run_config_fixture(tmp_path),base_directory=tmp_path)
    assert controller.request_cancel()
    assert not controller.request_cancel()
    events = wait(controller)
    assert controller.state == 'Cancelled'
    assert any(k == 'cancellation_effective' for k,p in events)
    assert not controller.request_cancel()


def test_failure_after_effective_cancellation_still_reports_failure(tmp_path,monkeypatch):
    import qcbed.gui.run_controller as module
    event = Event()
    connection = Connection(event)
    def fail(snapshot, *, progress_hook, **kwargs):
        event.set()
        try:
            progress_hook(dict(classification='trial',global_evaluation_index=2,stage='ANNEAL'))
        except KeyboardInterrupt:
            raise RuntimeError('finalization failure')
    monkeypatch.setattr(module,'run_snapshot',fail)
    _worker(prepare_run_snapshot(run_config_fixture(tmp_path),base_directory=tmp_path),False,connection,event)
    assert any(k == 'cancellation_effective' for k,p in connection.events)
    assert any(k == 'failure' and 'finalization failure' in p for k,p in connection.events)
    assert not any(k == 'result' for k,p in connection.events)
