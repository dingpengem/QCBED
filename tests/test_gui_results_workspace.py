"""Completed-event provenance, not optimizer-retained or mutable worker state."""
from dataclasses import replace
import json
from pathlib import Path

import numpy as np
import pytest
from qcbed import run_config
from qcbed.configuration import ParameterBinding, ParameterTarget, ParameterCategory
from qcbed.gui.results import load_results, completed_events
from test_gui_run import run_config_fixture
from test_gui_results import recorded


def test_live_completed_evaluation_and_cancel(tmp_path):
    config = run_config_fixture(tmp_path)
    target = ParameterTarget(ParameterCategory.DEPTH, (config.sample.regions[0].id,))
    config = replace(config, refinement=replace(config.refinement,
        parameters=(ParameterBinding(target, (1.1, .9), .1),),
        optimizer=replace(config.refinement.optimizer, loops=1)))
    directory = Path(config.run.output_directory)
    seen = []
    def observe(event):
        if event['classification'] not in ('initial_retained', 'trial', 'verification'):
            return
        view = load_results(directory)
        assert view.summary['displayed_evaluation'] == event['global_evaluation_index']
        assert [p['final'] for p in view.parameters] == event['physical_active_vector']
        assert view.summary['retained_objective'] == event['objective']
        assert view.summary['retained_ynorm'] == event['ynorm']
        assert view.comparison is not None
        with np.load(directory / 'latest_evaluation.npz') as arrays:
            np.testing.assert_array_equal(view.comparison['scaled'], arrays['scaled'])
        seen.append(event)
        if event['classification'] == 'trial' and sum(e['classification'] == 'trial' for e in seen) == 1:
            raise KeyboardInterrupt('test safe interruption')
    run_config(config, base_directory=tmp_path, progress_hook=observe)
    assert any(r['classification'] == 'trial' for r in seen)
    result = load_results(directory)
    assert result.summary['partial']
    latest = completed_events(result.history)[-1]
    assert result.summary['displayed_evaluation'] == latest['global_evaluation_index']
    assert [p['final'] for p in result.parameters] == latest['physical_active_vector']
    run_config(config, base_directory=tmp_path, resume=True)
    resumed = load_results(directory)
    assert resumed.history[:len(result.history)] == result.history


def test_history_selection_and_incomplete_record(recorded):
    _, directory = recorded
    result = load_results(directory)
    first = completed_events(result.history)[0]
    historical = load_results(directory, evaluation=first['global_evaluation_index'])
    assert [p['final'] for p in historical.parameters] == first['physical_active_vector']
    assert historical.summary['retained_objective'] == first['objective']
    before = result.summary['displayed_evaluation']
    with (directory / 'progress.jsonl').open('a') as stream:
        stream.write(json.dumps(dict(global_evaluation_index=999, classification='trial',
                                    objective=-1, physical_active_vector=[999.])))
    after = load_results(directory)
    assert after.summary['displayed_evaluation'] == before
    assert any('Unreadable progress' in w for w in after.warnings)
    with pytest.raises(ValueError):
        load_results(directory, evaluation=999)


def test_missing_historical_pattern_is_not_recomputed(recorded):
    _, directory = recorded
    result = load_results(directory)
    first = completed_events(result.history)[0]
    for p in directory.glob('*_snapshot.npz'):
        p.unlink()
    # The atomic current file belongs to a different completed evaluation.
    with np.load(directory / 'latest_evaluation.npz') as arrays:
        assert int(arrays['evaluation_index']) != first['global_evaluation_index']
    result = load_results(directory, evaluation=first['global_evaluation_index'])
    assert result.comparison is None
    assert result.parameters[0]['final'] == first['physical_active_vector'][0]


def test_uncommitted_current_pattern_not_used(recorded):
    _, directory = recorded
    path = directory / 'latest_evaluation.npz'
    with np.load(path) as arrays:
        values = {k: arrays[k] for k in arrays.files}
    values['evaluation_index'] = 999
    np.savez(path, **values)
    for p in directory.glob('*_snapshot.npz'):
        p.unlink()
    assert load_results(directory).comparison is None
