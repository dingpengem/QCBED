from pathlib import Path
import pytest
pytest.importorskip('PySide6')
from PySide6.QtWidgets import QApplication,QFileDialog
from qcbed.gui.window import ConfigWindow
from qcbed.gui.document import EditorDocument
from qcbed import load_yaml_config
from test_modern_grouped_run import modern_config


def test_one_project_yaml_lifecycle(tmp_path,monkeypatch):
    app=QApplication.instance() or QApplication([])
    w=ConfigWindow();w.show()
    assert w.document is None
    assert w.project_path_label.text()=='Project configuration: Unsaved'
    assert any(a.text()=='Load Configuration…' for a in w.edit_actions)
    w.document=EditorDocument(modern_config(tmp_path),tmp_path/'placeholder.yml')
    w.document.needs_save_as=True;w.refresh()
    assert not w.save_current()  # incomplete validation never publishes YAML
    for stage in range(3):w.workflow.accept(stage,w.document.config)
    dest=tmp_path/'project.yml';calls=[]
    def choose(*a,**kw):
        calls.append(1);return str(dest),''
    monkeypatch.setattr(QFileDialog,'getSaveFileName',choose)
    assert w.save_current()
    assert calls==[1] and w.document.path==dest
    assert w.project_path_label.text()==f'Project configuration: {dest}'
    for stage in range(3):w.workflow.accept(stage,w.document.config)
    assert w.save_current() and calls==[1]
    assert list(tmp_path.glob('*.yml'))==[dest]
    saved=load_yaml_config(dest)
    assert saved.crystallography.site_groups[0].debye_waller_b==.86
    assert saved.crystallography.slice_boundaries_fractional==w.document.config.crystallography.slice_boundaries_fractional
    assert saved.run.output_directory==w.document.config.run.output_directory
    w.site_group_pending=True
    assert not w.save_current()
    w.site_group_pending=False
    monkeypatch.setattr(QFileDialog,'getOpenFileName',lambda *a,**kw:(str(dest),''))
    w.choose_open()
    assert w.document.config==saved
    w.preparation.reset_preparation()
    assert w.preparation.associated_project_path==dest
    assert w.project_path_label.text()==f'Project configuration: {dest}'
    w.preparation.dirty=False;w.close()
