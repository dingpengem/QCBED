import os
os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')
import pytest
pytest.importorskip('PySide6')
from PySide6.QtWidgets import QApplication, QAbstractItemView
from qcbed.gui.results_view import ResultsView
from test_gui_results import recorded


def test_results_tables_and_native_actions(recorded, monkeypatch):
    from PySide6.QtGui import QDesktopServices
    app = QApplication.instance() or QApplication([])
    config, directory = recorded
    view = ResultsView()
    view.load(directory)
    assert 'Completed' in view.summary.toPlainText()
    assert view.parameters.rowCount() == 1
    assert view.parameters.item(0, 0).text() == view.results.parameters[0]['label']
    assert view.history.rowCount() == len(view.results.history)
    assert view.comparison.rowCount() == view.results.comparison['scaled'].shape[0]
    assert view.parameters.editTriggers() == QAbstractItemView.EditTrigger.NoEditTriggers
    urls = []
    monkeypatch.setattr(QDesktopServices, 'openUrl', lambda url: urls.append(url.toLocalFile()) or True)
    assert view.open_artifact('')
    assert view.open_artifact('report.html')
    assert urls == [str(directory), str(directory / 'report.html')]
    assert not view.open_artifact('../other')
    old = view.results
    with pytest.raises(ValueError):
        view.load(directory / 'missing')
    assert view.results is old
    view.close()


def test_window_existing_results_does_not_edit_document(recorded):
    from qcbed.gui.window import ConfigWindow
    from qcbed.gui.document import EditorDocument
    app = QApplication.instance() or QApplication([])
    config, directory = recorded
    window = ConfigWindow()
    window.document = EditorDocument(config, directory.parent / 'config.yml')
    window.refresh()
    original = window.document.config
    assert window.open_results(directory)
    assert window.document.config is original
    assert not window.document.dirty
    window.close()


def test_bounds_and_rejection_presentation(recorded):
    import json
    app = QApplication.instance() or QApplication([])
    _, directory = recorded
    path = directory / 'run_config.json'
    content = json.loads(path.read_text())
    content['refinement_bounds'] = {'entries': [
        dict(category='depth', target=[0], lower=1., upper=10.)]}
    path.write_text(json.dumps(content))
    with (directory / 'progress.jsonl').open('a') as stream:
        stream.write(json.dumps(dict(global_evaluation_index=99, classification='bound_rejection',
                                    objective='inf', stage='ANNEAL', restart_count=1)) + '\n')
    view = ResultsView()
    view.load(directory)
    assert "'lower': 1.0" in view.parameters.item(0, 6).text()
    assert view.history.item(view.history.rowCount()-1, 4).text() == 'BOUND REJECTION'
    assert view.history.item(view.history.rowCount()-1, 5).text() == '1'
    view.close()


def test_persistent_menu_workspace_no_automatic_popup(recorded):
    from types import SimpleNamespace
    from qcbed import prepare_run_snapshot
    from qcbed.gui.window import ConfigWindow
    app = QApplication.instance() or QApplication([])
    config, directory = recorded
    window = ConfigWindow()
    assert list(window.results_actions) == ['Show Results', 'Open Existing Run…', 'Open HTML Report', 'Open Output Folder']
    from qcbed.gui.document import EditorDocument
    window.document = EditorDocument(config, directory.parent / 'config.yml')
    window.refresh()
    initial = window.workspace.currentWidget()
    window.controller.snapshot = prepare_run_snapshot(config, base_directory=directory.parent)
    window.controller.state = 'Completed'
    window.poll_run()
    assert window.workspace.currentWidget() is initial
    assert window.results_view.results.directory == directory
    window.results_actions['Show Results'].trigger()
    assert window.workspace.currentWidget() is window.results_view
    view = window.results_view
    window.refresh_results(force=True)
    assert window.results_view is view
    assert window.workspace.currentWidget() is view
    assert window.results_actions['Open HTML Report'].isEnabled()
    window.controller.state = 'Running'
    window.controller.process = object()
    window.controller.cancel_pending = True
    window.refresh_results(force=True)
    assert 'Cancel pending' in view.summary.toPlainText()
    window.controller.process = None
    window.controller.state = 'Cancelled'
    window.refresh_results(force=True)
    assert 'Cancelled' in view.summary.toPlainText()
    assert 'Latest valid completed evaluation:' in view.summary.toPlainText()
    assert view.history.rowCount() > 0
    window.close()


def test_pattern_plot_duplicate_beams_and_history(recorded):
    import numpy as np
    app = QApplication.instance() or QApplication([])
    _, directory = recorded
    path = directory / 'latest_evaluation.npz'
    with np.load(path) as arrays:
        values = {k: arrays[k] for k in arrays.files}
    values['beams'] = np.repeat(values['beams'], 2, axis=0)
    values['experimental'] = np.repeat(values['experimental'], 2, axis=1)
    values['scaled'] = np.repeat(values['scaled'], 2, axis=1)
    np.savez(path, **values)
    view = ResultsView()
    view.load(directory)
    assert [view.tabs.tabText(i) for i in range(view.tabs.count())] == [
        'Overview', 'Refined Parameters', 'χ² History', 'Physical / YNORM history', 'Beam details', 'Pattern Matching', 'Report / Files']
    assert view.beams.count() == 2
    assert view.beams.itemText(0) != view.beams.itemText(1)
    assert 'Beam (0, 0)' in view.beams.itemText(1)
    view.beams.setCurrentIndex(1)
    assert len(view.pattern_figure.axes) == 6
    np.testing.assert_array_equal(view.pattern_figure.axes[0].collections[0].get_array(),
                                  values['experimental'][:, 1].reshape(2, 2))
    view.common_scale.setChecked(False)
    assert len(view.pattern_figure.axes) == 6
    assert len(view.history_figure.axes[0].collections[0].get_offsets()) > 0
    view.evaluations.setCurrentIndex(1)
    assert view.results.summary['historical']
    assert view.selection == view.results.summary['displayed_evaluation']
    view.close()


def test_switch_existing_and_current_run(recorded, tmp_path):
    from dataclasses import replace
    from qcbed import run_config, prepare_run_snapshot
    from qcbed.gui.window import ConfigWindow
    app = QApplication.instance() or QApplication([])
    config, directory = recorded
    other = replace(config, run=replace(config.run,
        output_directory=str(tmp_path / 'other'), output_filename=str(tmp_path / 'other.dat')))
    run_config(other, base_directory=tmp_path)
    window = ConfigWindow()
    window.controller.snapshot = prepare_run_snapshot(config, base_directory=tmp_path)
    window.controller.state = 'Completed'
    window.follow_run_results()
    assert window.results_view.results.directory == directory
    assert window.open_results(tmp_path / 'other')
    window.refresh_results(force=True)
    assert window.results_view.results.directory == tmp_path / 'other'
    window.follow_run_results()
    assert window.results_view.results.directory == directory
    window.close()
