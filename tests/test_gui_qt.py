"""Optional offscreen Qt interaction tests; no pytest-qt dependency."""
import os
from pathlib import Path

os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')
import pytest
pytest.importorskip('PySide6')
from PySide6.QtCore import Qt
from PySide6.QtWidgets import QApplication
from qcbed import Id, load_yaml_config
from qcbed.gui.window import ConfigWindow

REFERENCE = Path(__file__).resolve().parent / 'data/reference_config.yml'


@pytest.fixture(scope='module')
def app():
    return QApplication.instance() or QApplication([])


@pytest.fixture
def window(app):
    window = ConfigWindow(REFERENCE)
    window.show()
    app.processEvents()
    yield window
    window.pending = False
    window.document.saved_config = window.document.config
    window.close()
    window.deleteLater()
    app.processEvents()


def test_open_and_unchanged_apply(window):
    original = load_yaml_config(REFERENCE)
    assert window.document.config == original
    assert 'Active parameters: 7' in window.summary.toPlainText()
    assert window.refinement_table.rowCount() > 18
    window.mark_pending()
    assert window.apply_edits()
    assert window.document.config == original


def test_edits_save_reload_preserve_unedited(window, tmp_path):
    old = window.document.config
    region = old.sample.regions[0]
    window.widgets[('sample','regions',Id(region.id),'depth_angstrom')].setText('912.1234567890123')
    window.widgets[('scan','scan_rows')].setText('3')
    window.widgets[('scan','scan_columns')].setText('4')
    window.widgets[('scan','corner_b')].setText('0.41234567890123456, 0.52')
    window.widgets[('experiment','path')].setText(str(tmp_path/'measurements.dat'))
    assert window.document.config is old
    assert window.apply_edits()
    current = window.document.config
    assert current.sample.regions[0].depth_angstrom == 912.1234567890123
    assert current.sample.regions[0].slices == region.slices
    assert current.refinement == old.refinement
    assert current.microscope == old.microscope
    assert 'Tilts: 12' in window.summary.toPlainText()
    destination = tmp_path/'saved.yml'
    assert window.save_document(destination)
    assert load_yaml_config(destination) == current
    window.open_document(destination)
    assert window.document.config == current


def find_row(window, target):
    return next(row for row in range(window.refinement_table.rowCount())
        if window.refinement_table.item(row, 0).data(Qt.ItemDataRole.UserRole) == target)


def test_sorted_target_binding_identity(window):
    old = window.document.config
    binding = old.refinement.parameters[0]
    table = window.refinement_table
    table.sortItems(0, Qt.SortOrder.DescendingOrder)
    row = find_row(window, binding.target)
    table.item(row, 4).setText('-0.125')
    table.item(row, 5).setText('0.375')
    table.item(row, 6).setText('0.0012345678901234567')
    assert window.apply_edits()
    changed = next(b for b in window.document.config.refinement.parameters if b.target == binding.target)
    assert changed.range_factors == (-0.125, 0.375)
    assert changed.plamda == 0.0012345678901234567
    assert window.document.config.sample == old.sample
    row = find_row(window, binding.target)
    window.refinement_table.item(row, 1).setCheckState(Qt.CheckState.Unchecked)
    assert window.apply_edits()
    row = find_row(window, binding.target)
    window.refinement_table.item(row, 1).setCheckState(Qt.CheckState.Checked)
    for col, value in [(4,'-0.1'),(5,'0.2'),(6,'0.01')]: window.refinement_table.item(row,col).setText(value)
    assert window.apply_edits()
    assert any(b.target == binding.target for b in window.document.config.refinement.parameters)


def test_invalid_draft_issue_navigation_and_revert(window, tmp_path):
    original = window.document.config
    widget = window.widgets[('scan','scan_rows')]
    widget.setText('2.5')
    assert not window.save_document(tmp_path/'invalid.yml')
    assert not (tmp_path/'invalid.yml').exists()
    assert window.document.config is original
    assert widget.text() == '2.5'
    assert window.issues.topLevelItemCount()
    item = window.issues.topLevelItem(0)
    window.navigate_issue(item)
    assert window.navigation.currentItem().text(0) == 'Scan'
    window.refresh()
    assert window.widgets[('scan','scan_rows')].text() == str(original.scan.scan_rows)


def test_file_readiness_issues(window, tmp_path):
    window.widgets[('experiment','path')].setText(str(tmp_path/'missing.dat'))
    window.validate_document()
    assert any(window.issues.topLevelItem(i).text(1) == 'External file/output readiness' for i in range(window.issues.topLevelItemCount()))
    assert 'External file/output readiness: Blocked' in window.summary.toPlainText()
    assert 'Configuration validity: Valid' in window.summary.toPlainText()


def test_common_widgets_and_revert_clear_draft_errors(window, tmp_path):
    window.widgets[('title',)].setText('Edited title')
    window.widgets[('microscope','energy_kv')].setText('210.12345678901235')
    window.widgets[('microscope','zone_axis')].setText('1, 0, 0')
    window.widgets[('scan','basis_1')].setText('1, 1, 0')
    window.widgets[('scan','basis_2')].setText('0, 1, 1')
    window.widgets[('scan','corner_a')].setText('-0.5, 0.3')
    window.widgets[('scan','corner_c')].setText('-0.4, -0.5')
    window.widgets[('execution','execution_backend')].setCurrentText('serial')
    window.widgets[('execution','batch_size')].setText('')
    window.widgets[('run','output_directory')].setText(str(tmp_path/'out'))
    window.widgets[('run','output_filename')].setText(str(tmp_path/'result.dat'))
    assert window.apply_edits()
    assert window.document.config.title == 'Edited title'
    assert window.document.config.microscope.energy_kv == 210.12345678901235
    assert window.document.config.scan.basis_2 == (0,1,1)
    assert window.document.config.execution.execution_backend == 'serial'
    window.widgets[('scan','scan_rows')].setText('invalid')
    assert not window.apply_edits()
    window.refresh()
    assert not any(i.code == 'invalid_editor_input' for i in window.document.issues)


def test_existing_output_readiness_presentation(window, tmp_path):
    window.widgets[('run','output_directory')].setText(str(tmp_path))
    window.validate_document()
    assert not any(i.code=='output_exists' for i in window.document.issues)
    assert not window.run_action.isEnabled()  # A/B/C confirmation is still required.
    assert 'Configuration validity: Valid' in window.summary.toPlainText()
    assert window.save_document(tmp_path/'valid.yml')


def test_run_readiness_controls(window, tmp_path):
    window.widgets[('run','output_directory')].setText(str(tmp_path))
    window.validate_document()
    assert not window.run_action.isEnabled()
    assert 'Validate the current configuration' in window.run_action.toolTip()  # Legacy fixture has no modern A/B/C receipts.
    assert not window.resume_action.isEnabled()
    assert window.resume_action.toolTip()


def test_gui_spawn_responsive_completion(window, app, tmp_path):
    import time
    from PySide6.QtCore import QTimer
    from test_gui_run import run_config_fixture
    from qcbed.gui.document import EditorDocument
    window.document = EditorDocument(run_config_fixture(tmp_path),tmp_path/'config.yml')
    window.refresh()
    for stage in range(3):window.workflow.accept(stage,window.document.config)
    window.validate_document()
    assert window.run_action.isEnabled()
    window.start_run()
    assert not window.run_action.isEnabled()
    assert not window.resume_action.isEnabled()
    assert not window.pages.isEnabled()
    ticks = []
    timer = QTimer()
    timer.setInterval(10)
    timer.timeout.connect(lambda: ticks.append(1))
    timer.start()
    deadline = time.monotonic()+90
    while window.controller.active or window.timer.isActive():
        app.processEvents()
        assert time.monotonic() < deadline
        time.sleep(.005)
    timer.stop()
    assert ticks
    assert window.controller.state == 'Completed'
    assert 'Completed' in window.run_status.text()
    assert 'classification' in window.run_log.toPlainText()
    assert window.pages.isEnabled()
    assert window.resume_action.isEnabled()


def test_gui_worker_failure_display(window, monkeypatch):
    window.controller.state = 'Failed'
    monkeypatch.setattr(window.controller,'poll',lambda: [('failure','Traceback: worker failure test')])
    window.poll_run()
    assert window.run_status.text().startswith('Failed')
    assert 'Traceback: worker failure test' in window.run_log.toPlainText()
    assert window.pages.isEnabled()


def test_cooperative_cancel_restores_gui(window,app,tmp_path):
    import time
    from test_gui_run import run_config_fixture
    from qcbed.gui.document import EditorDocument
    window.document = EditorDocument(run_config_fixture(tmp_path),tmp_path/'cancel.yml')
    window.refresh()
    for stage in range(3):window.workflow.accept(stage,window.document.config)
    assert not window.cancel_action.isEnabled()
    window.start_run()
    assert window.cancel_action.isEnabled()
    window.cancel_action.trigger()
    assert not window.cancel_action.isEnabled()
    assert window.controller.cancel_pending
    assert 'pending' in window.run_status.text()
    assert not window.pages.isEnabled()
    window.cancel_run()
    assert window.run_log.toPlainText().count('Cancellation requested;') == 1
    deadline = time.monotonic()+90
    ticks = 0
    while window.controller.active or window.timer.isActive():
        app.processEvents()
        ticks += 1
        assert time.monotonic() < deadline
        time.sleep(.005)
    # A pre-evaluation cancellation may finish in the first event-loop turn.
    assert ticks >= 1
    assert window.controller.state == 'Cancelled'
    assert window.run_status.text().startswith('Cancelled')
    assert window.pages.isEnabled()
    assert not window.cancel_action.isEnabled()
    assert 'boundary' in window.run_log.toPlainText()


@pytest.mark.parametrize('state', ['absent', 'compatible', 'missing', 'incompatible', 'corrupt'])
def test_output_directory_run_resume_states(window, tmp_path, state):
    from dataclasses import replace
    from pathlib import Path
    from qcbed import run_config
    from qcbed.gui.document import EditorDocument
    from test_gui_run import run_config_fixture
    config = run_config_fixture(tmp_path)
    directory = Path(config.run.output_directory)
    if state == 'missing':
        directory.mkdir()
    elif state in ('compatible', 'incompatible', 'corrupt'):
        run_config(config, base_directory=tmp_path)
        if state == 'incompatible':
            config = replace(config, run=replace(config.run,
                checkpoint_every=config.run.checkpoint_every + 1))
        elif state == 'corrupt':
            (directory / 'checkpoint_latest.json').write_text('{broken json')
    before = {p: p.read_bytes() for p in tmp_path.rglob('*') if p.is_file()}
    window.document = EditorDocument(config, tmp_path / 'config.yml')
    window.refresh()
    for stage in range(3):window.workflow.accept(stage,window.document.config)
    window.validate_document()
    assert window.run_action.isEnabled()
    assert window.resume_action.isEnabled() == (state == 'compatible')
    text = window.run_action.toolTip() + window.resume_action.toolTip()
    assert 'FileNotFoundError' not in text and '[Errno 2]' not in text
    if state == 'missing':
        assert 'Resume: Not available — no previous run found' in text
        assert 'or use Resume' not in text
    if state == 'incompatible':
        assert 'identity mismatch' in text
    if state == 'corrupt':
        assert 'Expecting property name' in text
    if state == 'absent':
        assert 'output directory does not exist' in text
    assert before == {p: p.read_bytes() for p in tmp_path.rglob('*') if p.is_file()}
