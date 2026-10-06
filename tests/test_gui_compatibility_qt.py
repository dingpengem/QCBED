import os
os.environ.setdefault('QT_QPA_PLATFORM','offscreen')
import pytest
pytest.importorskip('PySide6')
from PySide6.QtWidgets import QDialog, QDialogButtonBox
from qcbed.gui.compatibility_widgets import ImportDialog,ExportDialog,ReportDialog
from qcbed.gui.compatibility import compatibility_report,import_document
from test_gui_qt import app,window
from test_legacy_config import ROOT,policy


def fill_import(dialog,tmp_path):
    values = dict(source=str(ROOT/'tests/fixtures/legacy_reference/FOR007.DAT'),experimental=str(ROOT/'tests/fixtures/legacy_reference/in.dat'),
        seed='123',batch_size='8',execution='batched',objective='masked',invalid='reject',bounds='none',
        output_directory=str(tmp_path/'run'),output_filename=str(tmp_path/'out.dat'),
        snapshot_every='0',checkpoint_every='1',max_evaluations='50',max_seconds='60')
    for key,value in values.items():
        control = dialog.controls[key]
        if hasattr(control,'setText'): control.setText(value)
        else: control.setCurrentText(value)


def test_explicit_policy_import_and_save_as(window,tmp_path,monkeypatch):
    dialog = ImportDialog(window)
    dialog.submit()
    assert 'explicit' in dialog.error.text()
    assert dialog.result_values is None
    fill_import(dialog,tmp_path)
    dialog.controls['bounds'].setCurrentIndex(0)
    dialog.submit()
    assert dialog.result_values is None
    dialog.controls['bounds'].setCurrentText('none')
    dialog.controls['batch_size'].clear()
    dialog.submit()
    assert 'explicit batch' in dialog.error.text()
    dialog.controls['batch_size'].setText('8')
    dialog.submit()
    source,experimental,chosen = dialog.result_values
    assert chosen.objective_mode == 'masked'
    assert chosen.invalid_trial_policy == 'reject'
    window.import_legacy_document(source,experimental,chosen)
    assert window.document.dirty and window.document.needs_save_as
    assert 'Optimizer parameters: 18' in window.summary.toPlainText()
    requests = []
    monkeypatch.setattr(window,'choose_save_as',lambda: requests.append('save_as') or True)
    assert window.save_current()
    assert requests == ['save_as']
    assert window.save_document(tmp_path/'canonical.yml')
    assert not window.document.dirty
    assert not window.document.needs_save_as


def test_import_failure_preserves_current_document(window,tmp_path):
    original = window.document
    source = tmp_path/'bad.dat'
    source.write_text((ROOT/'tests/fixtures/legacy_reference/FOR007.DAT').read_text()+'\nAMBIGUOUS\n')
    with pytest.raises(ValueError,match='trailing'):
        window.import_legacy_document(source,ROOT/'tests/fixtures/legacy_reference/in.dat',policy(tmp_path))
    assert window.document is original


def test_export_actions_review_omissions_and_protect_history(window,tmp_path,monkeypatch):
    window.document = import_document(ROOT/'tests/fixtures/legacy_reference/FOR007.DAT',ROOT/'tests/fixtures/legacy_reference/in.dat',policy(tmp_path))
    window.refresh()
    config = window.document.config
    history = list(window.document.undo_stack)
    reports = []
    def options(dialog):
        dialog.destination.setText(str(tmp_path/('bundle' if dialog.bundle else 'FOR007.dat')))
        dialog.exp_width.setText('4096'); dialog.out_width.setText('4096')
        dialog.submit()
        return dialog.result()
    def review(dialog):
        reports.append(dialog.text.toPlainText())
        return QDialog.DialogCode.Accepted
    monkeypatch.setattr(ExportDialog,'exec',options)
    monkeypatch.setattr(ReportDialog,'exec',review)
    window.export_legacy(False)
    window.export_legacy(True)
    assert (tmp_path/'FOR007.dat').exists()
    assert (tmp_path/'bundle/FOR007.dat').exists()
    assert (tmp_path/'bundle/in.dat').exists()
    assert (tmp_path/'bundle/compatibility.json').exists()
    assert len(reports) == 2
    assert all('FAITHFULLY REPRESENTABLE' in t and 'OMITTED' in t and 'UNSUPPORTED' in t for t in reports)
    assert 'seed' in reports[0] and 'bounds' in reports[0]
    assert window.document.config is config
    assert window.document.undo_stack == history
    assert window.document.dirty
    protected = ExportDialog(tmp_path,bundle=False)
    protected.destination.setText(str(tmp_path/'FOR007.dat'))
    protected.submit()
    assert 'already exists' in protected.error.text()
    assert protected.result_values is None


def test_readonly_report_blocks_unsupported_export(app,tmp_path):
    from dataclasses import replace
    from qcbed.configuration import LegacyExportConfig
    document = import_document(ROOT/'tests/fixtures/legacy_reference/FOR007.DAT',ROOT/'tests/fixtures/legacy_reference/in.dat',policy(tmp_path))
    report = compatibility_report(replace(document.config,legacy_export=LegacyExportConfig()),base_directory=tmp_path)
    dialog = ReportDialog(report,exporting=True)
    assert dialog.text.isReadOnly()
    assert 'DSMAX' in dialog.text.toPlainText()
    assert not dialog.findChild(QDialogButtonBox).button(QDialogButtonBox.StandardButton.Ok).isEnabled()
