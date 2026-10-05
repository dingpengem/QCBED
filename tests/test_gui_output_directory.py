from dataclasses import replace
import yaml
from qcbed.gui.document import EditorDocument,editable_fields
from qcbed import load_yaml_config
from test_modern_grouped_run import modern_config


def test_modern_directory_only_output_and_old_yaml(tmp_path):
    config=modern_config(tmp_path)
    doc=EditorDocument(config,tmp_path/'project.yml')
    assert ('run','output_filename') not in {s.key for s in editable_fields(doc.config)}
    directory=tmp_path/'new-run'
    assert doc.apply({('run','output_directory'):str(directory)})
    assert doc.config.run.output_filename==str(directory/'final_out.dat')
    doc.save(tmp_path/'saved.yml')
    raw=yaml.safe_load((tmp_path/'saved.yml').read_text())
    assert 'filename' not in raw['qcbed']['output']
    assert 'output_file' not in raw['qcbed']['output']
    assert load_yaml_config(tmp_path/'saved.yml')==doc.config
    raw['qcbed']['output']['filename']=str(tmp_path/'old-external.dat')
    (tmp_path/'older.yml').write_text(yaml.safe_dump(raw))
    loaded=load_yaml_config(tmp_path/'older.yml')
    assert loaded.run.output_filename==str(directory/'final_out.dat')

import errno
from pathlib import Path
import pytest
from qcbed.gui.output_paths import validate_output_directory,modern_execution_config
from qcbed.gui.readiness import readiness,NO_CHECKPOINT


@pytest.mark.parametrize('contents',[(),('notes.txt',),('checkpoint_latest.json','final_out.dat')])
def test_existing_directory_valid_and_preserved(tmp_path,contents):
    for name in contents:(tmp_path/name).write_text('user data')
    assert validate_output_directory(str(tmp_path),tmp_path,create=True)==tmp_path
    assert {p.name:p.read_text() for p in tmp_path.iterdir()}=={n:'user data' for n in contents}


def test_missing_directory_created_and_file_rejected(tmp_path):
    new=tmp_path/'parent'/'results'
    assert validate_output_directory(str(new),tmp_path,create=True)==new
    assert new.is_dir()
    file=tmp_path/'notes';file.write_text('untouched')
    with pytest.raises(ValueError,match='file; select a directory'):
        validate_output_directory(str(file),tmp_path,create=True)
    with pytest.raises(ValueError,match='ancestor is a file'):
        validate_output_directory(str(file/'child'),tmp_path,create=True)
    with pytest.raises(ValueError,match='Select an output directory'):
        validate_output_directory(' ',tmp_path)


def test_nonwritable_directory_clear_error(tmp_path,monkeypatch):
    import tempfile
    def denied(*a,**kw):raise PermissionError(errno.EACCES,'Permission denied')
    monkeypatch.setattr(tempfile,'TemporaryFile',denied)
    with pytest.raises(ValueError,match='Cannot write/create output directory: Permission denied'):
        validate_output_directory(str(tmp_path),tmp_path,create=True)


def test_fresh_ready_independently_of_resume_and_old_outputs(tmp_path):
    c=modern_config(tmp_path);out=Path(c.run.output_directory)
    out.mkdir();(out/'notes.txt').write_text('keep')
    (out/'final_out.dat').write_text('old results')
    assert not [i for i in readiness(c,tmp_path) if i.severity=='error']
    assert [i.message for i in readiness(c,tmp_path,resume=True)]==[NO_CHECKPOINT]
    first=modern_execution_config(c,tmp_path);second=modern_execution_config(c,tmp_path)
    assert first.run.output_directory!=second.run.output_directory
    assert Path(first.run.output_directory).parent==out
    assert sorted(p.name for p in out.iterdir())==['final_out.dat','notes.txt']
    # A corrupt old checkpoint does not block a new Fresh Run.
    (out/'checkpoint_latest.json').write_text('{}')
    assert not [i for i in readiness(c,tmp_path) if i.severity=='error']
    assert modern_execution_config(c,tmp_path,resume=True).run.output_directory==str(out)
    child=Path(first.run.output_directory);child.mkdir();(child/'checkpoint_latest.json').write_text('{}')
    assert modern_execution_config(c,tmp_path,resume=True).run.output_directory==str(child)


def test_dashboard_output_is_single_editor_and_separate_step(tmp_path):
    pytest.importorskip('PySide6')
    from PySide6.QtWidgets import QApplication,QLineEdit
    from qcbed.gui.window import ConfigWindow
    app=QApplication.instance() or QApplication([])
    w=ConfigWindow();w.document=EditorDocument(modern_config(tmp_path),tmp_path/'project.yml');w.refresh()
    try:
        c=w.document.config
        w.workflow.accept(0,c);assert w.validate_stage(1)
        assert w.validate_c_part('simulation');assert w.validate_stage(2)
        assert not w.run_action.isEnabled()
        editor=w.widgets[('run','output_directory')]
        assert w.run_dashboard.findChildren(QLineEdit).count(editor)==1
        assert editor not in w.c_tabs.findChildren(QLineEdit)
        new=tmp_path/'new-output'
        editor.setText(str(new))
        assert not w.run_action.isEnabled() and w.workflow.ready
        assert w.validate_output_directory()
        assert new.is_dir() and w.workflow.ready
        w.check_run_readiness()
        assert w.run_action.isEnabled() and not w.resume_action.isEnabled()
        assert w.run_dashboard.status['Output directory'].text()=='Valid'
        assert w.run_dashboard.status['Run'].text()=='Ready'
        assert w.run_dashboard.status['Resume'].text()=='Not available — no previous run found'
        assert w.document.config.run.output_directory==str(new)
        editor.setText(str(tmp_path/'project-file'));(tmp_path/'project-file').write_text('keep')
        assert not w.validate_output_directory()
        assert 'file' in w.run_dashboard.status['Output directory'].text()
    finally:
        w.document.saved_config=w.document.config;w.pending=False;w.output_draft=False;w.close();w.deleteLater();app.processEvents()


def test_uncreatable_directory_clear_error(tmp_path,monkeypatch):
    def denied(*a,**kw):raise PermissionError(errno.EACCES,'Permission denied')
    monkeypatch.setattr(Path,'mkdir',denied)
    with pytest.raises(ValueError,match='Cannot write/create output directory: Permission denied'):
        validate_output_directory(str(tmp_path/'missing'),tmp_path,create=True)


def test_lightweight_validation_does_not_compile_or_inspect_metadata(tmp_path,monkeypatch):
    def forbidden(*a,**kw):pytest.fail('Directory validation must not compile or read checkpoints')
    monkeypatch.setattr('qcbed.modern_frontend.validate_modern',forbidden)
    monkeypatch.setattr('qcbed.grouped_run.load_grouped_checkpoint',forbidden)
    (tmp_path/'checkpoint_latest.json').write_text('not even valid JSON')
    assert validate_output_directory(str(tmp_path),tmp_path,create=True)==tmp_path


def test_c_valid_without_output_and_d_owns_readiness(tmp_path):
    pytest.importorskip('PySide6')
    from PySide6.QtWidgets import QApplication,QLineEdit,QComboBox
    from qcbed.gui.window import ConfigWindow
    app=QApplication.instance() or QApplication([])
    c=modern_config(tmp_path);c=replace(c,run=replace(c.run,output_directory=''))
    w=ConfigWindow();w.document=EditorDocument(c,tmp_path/'project.yml');w.refresh()
    try:
        w.workflow.accept(0,c);assert w.validate_stage(1)
        assert w.validate_c_part('simulation');assert w.validate_stage(2)
        assert w.workflow.ready and not w.document.files_checked
        assert not w.run_action.isEnabled()
        assert 'Output directory not selected' in w.run_dashboard.reason.text()
        assert w.widgets[('run','output_directory')].text()==''
        assert w.c_tabs.count()==3
        assert w.save_document(tmp_path/'without-output.yml')
        assert load_yaml_config(tmp_path/'without-output.yml').run.output_directory==''
        for key in (('run','output_directory'),('execution','execution_backend'),('execution','batch_size')):
            assert w.run_dashboard.isAncestorOf(w.widgets[key])
            assert not w.c_tabs.isAncestorOf(w.widgets[key])
        # C does not consume or validate an independent D output draft.
        file=tmp_path/'not-directory';file.write_text('keep')
        w.widgets[('run','output_directory')].setText(str(file))
        assert w.validate_stage(2) and w.workflow.ready
        assert not w.validate_output_directory() and w.workflow.ready
        w.widgets[('run','output_directory')].setText(str(tmp_path/'results'))
        assert w.validate_stage(2) and not w.run_action.isEnabled()
        assert w.validate_output_directory();w.check_run_readiness()
        assert w.workflow.ready and w.run_action.isEnabled()
        assert w.save_document(tmp_path/'saved.yml')
        assert load_yaml_config(tmp_path/'saved.yml').run.output_directory==str(tmp_path/'results')
    finally:
        w.document.saved_config=w.document.config;w.pending=False;w.output_draft=False;w.close();w.deleteLater();app.processEvents()


def test_batched_64_advanced_d_persists_without_changing_c(tmp_path):
    pytest.importorskip('PySide6')
    from PySide6.QtWidgets import QApplication
    from qcbed.gui.window import ConfigWindow
    app=QApplication.instance() or QApplication([])
    c=modern_config(tmp_path)
    c=replace(c,execution=replace(c.execution,batch_size=32),
        run=replace(c.run,checkpoint_every=7,snapshot_every=11))
    w=ConfigWindow();w.document=EditorDocument(c,tmp_path/'project.yml');w.refresh()
    try:
        original=c.execution
        assert not w.execution_advanced.toggle.isChecked()
        assert w.document.config.execution==original
        w.workflow.accept(0,c);assert w.validate_stage(1)
        assert w.validate_c_part('simulation');assert w.validate_stage(2)
        w.workflow_tabs.setCurrentIndex(3);w.show();w.execution_advanced.set_expanded(True);app.processEvents()
        mode=w.widgets[('execution','execution_backend')];batch=w.widgets[('execution','batch_size')]
        for _ in range(5):app.processEvents()
        for widget in (mode,batch,w.widgets[('run','max_total_evaluations')],w.widgets[('run','max_wall_time_seconds')]):
            w.run_dashboard.content_scroll.ensureWidgetVisible(widget)
            for _ in range(5):app.processEvents()
            assert not widget.visibleRegion().isEmpty()
        assert mode.isVisible() and batch.isVisible()
        w.controller.process=object();w.run_dashboard.update_state()
        assert not mode.isEnabled() and not batch.isEnabled()
        w.controller.process=None;w.run_dashboard.update_state()
        assert w.execution_advanced.isAncestorOf(mode) and not w.c_tabs.isAncestorOf(mode)
        mode.setCurrentText('batched');batch.setText('64')
        assert w.workflow.ready
        assert w.apply_edits();assert w.workflow.ready
        assert w.document.config.execution.execution_backend=='batched'
        assert w.document.config.execution.batch_size==64
        assert w.save_document(tmp_path/'batched.yml')
        loaded=load_yaml_config(tmp_path/'batched.yml')
        assert loaded.execution==w.document.config.execution
        assert loaded.crystallography==c.crystallography
        assert loaded.run.checkpoint_every==7 and loaded.run.snapshot_every==11
        w.open_document(tmp_path/'batched.yml')
        assert w.widgets[('execution','batch_size')].text()=='64'
        assert w.widgets[('execution','execution_backend')].currentText()=='batched'
        assert not w.execution_advanced.toggle.isChecked()
    finally:
        w.document.saved_config=w.document.config;w.pending=False;w.output_draft=False;w.close();w.deleteLater();app.processEvents()
