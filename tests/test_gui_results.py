"""Artifact-only results inspection; never evaluates numerical work on load."""
from dataclasses import replace
import json
from pathlib import Path

import numpy as np
import pytest

from qcbed import run_config
from qcbed.configuration import ParameterBinding, ParameterTarget, ParameterCategory
from qcbed.gui.results import load_results
from test_gui_run import run_config_fixture


@pytest.fixture
def recorded(tmp_path):
    config = run_config_fixture(tmp_path)
    target = ParameterTarget(ParameterCategory.DEPTH, (config.sample.regions[0].id,))
    config = replace(config, refinement=replace(config.refinement,
        parameters=(ParameterBinding(target, (1.1, .9), .1),)))
    run_config(config, base_directory=tmp_path)
    return config, Path(config.run.output_directory)


def test_completed_values_and_read_only(recorded, monkeypatch):
    config, directory = recorded
    before = {p: p.read_bytes() for p in directory.rglob('*') if p.is_file()}
    import qcbed.full_runner
    monkeypatch.setattr(qcbed.full_runner, 'evaluate_refinement_trial',
                        lambda *a, **k: pytest.fail('Viewer must not compute a trial'))
    result = load_results(directory)
    assert not result.summary['partial']
    assert result.parameters[0]['category'] == 'depth'
    assert config.sample.regions[0].id in result.parameters[0]['label']
    assert result.parameters[0]['final'] == np.load(directory / 'final_active_vector.npy')[0]
    assert result.parameters[0]['delta'] == result.parameters[0]['final'] - result.parameters[0]['initial']
    with np.load(directory / 'final_snapshot.npz') as arrays:
        np.testing.assert_array_equal(result.comparison['experimental'], arrays['experimental'])
        np.testing.assert_array_equal(result.comparison['scaled'], arrays['scaled'])
        assert result.summary['final_ynorm'] == arrays['ynorm']
        assert result.summary['final_objective'] == arrays['objective']
    assert not result.comparison['scaled'].flags.writeable
    assert result.history == tuple(json.loads(line) for line in (directory / 'progress.jsonl').read_text().splitlines())
    assert before == {p: p.read_bytes() for p in directory.rglob('*') if p.is_file()}


def test_resume_history_continuity(recorded, tmp_path):
    config, directory = recorded
    previous = load_results(directory).history
    run_config(config, base_directory=tmp_path, resume=True)
    result = load_results(directory)
    assert result.history[:len(previous)] == previous
    assert result.summary['restart_count'] == 1
    assert {row['restart_count'] for row in result.history} == {0, 1}


def test_partial_and_truncated_history(recorded):
    _, directory = recorded
    (directory / 'summary.json').unlink()
    with (directory / 'progress.jsonl').open('a') as stream:
        stream.write('{"incomplete":')
    result = load_results(directory)
    assert result.summary['partial']
    assert result.summary['status'] == 'incomplete'
    assert any('Unreadable progress' in w for w in result.warnings)
    assert result.parameters[0]['final'] is not None


def test_older_output_and_bounds_rejection(recorded):
    _, directory = recorded
    path = directory / 'run_config.json'
    config = json.loads(path.read_text())
    config.pop('results_metadata')
    (directory / 'latest_evaluation.npz').unlink()
    config['refinement_bounds'] = {'entries': [
        dict(category='depth', target=[0], lower=1., upper=10.)]}
    path.write_text(json.dumps(config))
    for path in directory.glob('*_snapshot.npz'):
        with np.load(path) as arrays:
            values = {k: arrays[k] for k in arrays.files if k not in ('experimental', 'beams', 'scan_shape')}
        np.savez(path, **values)
    with (directory / 'progress.jsonl').open('a') as stream:
        stream.write(json.dumps(dict(global_evaluation_index=99, classification='bound_rejection',
                                    objective='inf', restart_count=1)) + '\n')
    result = load_results(directory)
    assert result.parameters[0]['bounds']['lower'] == 1.
    assert result.history[-1]['classification'] == 'bound_rejection'
    assert result.comparison is None
    assert any('Stable target IDs' in w for w in result.warnings)


@pytest.mark.parametrize('filename,content', [
    ('checkpoint_latest.json', '{}'), ('run_config.json', '[]'),
    ('checkpoint_latest.json', '{broken'), ('summary.json', '[]')])
def test_malformed_artifacts(recorded, filename, content):
    _, directory = recorded
    (directory / filename).write_text(content)
    with pytest.raises(ValueError, match='Cannot inspect canonical run'):
        load_results(directory)


def test_unrecognized_directory(tmp_path):
    with pytest.raises(ValueError, match='Cannot inspect canonical run'):
        load_results(tmp_path)
    assert not list(tmp_path.iterdir())


def test_cancelled_retained_state(tmp_path):
    from qcbed.gui.run_controller import RunController
    from test_gui_run import wait
    config = run_config_fixture(tmp_path)
    controller = RunController()
    controller.start(config, base_directory=tmp_path)
    controller.request_cancel()
    wait(controller)
    directory = Path(config.run.output_directory)
    if (directory / 'checkpoint_latest.json').exists():
        result = load_results(directory)
        assert result.summary['partial']
    else:
        with pytest.raises(ValueError):
            load_results(directory)


def test_cooperative_cancel_retained_results_and_resume(tmp_path):
    from threading import Event
    from qcbed import prepare_run_snapshot
    from qcbed.gui.run_controller import _worker
    from test_gui_cancel import Connection
    config = run_config_fixture(tmp_path)
    target = ParameterTarget(ParameterCategory.DEPTH, (config.sample.regions[0].id,))
    config = replace(config, refinement=replace(config.refinement,
        parameters=(ParameterBinding(target, (1.1, .9), .1),),
        optimizer=replace(config.refinement.optimizer, loops=1)))
    directory = Path(config.run.output_directory)
    event = Event()
    connection = Connection(event, 'initial_retained', directory)
    _worker(prepare_run_snapshot(config, base_directory=tmp_path), False, connection, event)
    assert next(p for k, p in connection.events if k == 'result')['status'] == 'cancelled'
    partial = load_results(directory)
    assert partial.summary['partial']
    assert partial.parameters[0]['final'] is not None
    run_config(config, base_directory=tmp_path, resume=True)
    final = load_results(directory)
    assert not final.summary['partial']
    assert final.history[:len(partial.history)] == partial.history
