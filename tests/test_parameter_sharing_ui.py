"""Presentation-only sharing patch: coordinates and persistence remain canonical."""
import pytest
from test_projected_elastic import config
from qcbed import compile_config_from_sources,save_yaml_config,load_yaml_config


@pytest.mark.parametrize('mode,before,after',[('real',3,2),('complex',5,3),('magnitude',3,2)])
def test_sharing_controls_counts_confirmation_and_independence(tmp_path,mode,before,after,monkeypatch):
    from PySide6.QtWidgets import QApplication,QPushButton,QInputDialog,QMessageBox
    from PySide6.QtCore import QItemSelectionModel
    from qcbed.gui.window import ConfigWindow
    from qcbed.gui.document import EditorDocument
    app=QApplication.instance() or QApplication([])
    c=config(tmp_path,mode,reflections=((2,0,0),(0,2,0)))
    ordinary=sum(p.category!='projected_elastic' for p in compile_config_from_sources(c).layout.parameters)
    before+=ordinary-1;after+=ordinary-1
    w=ConfigWindow();w.document=EditorDocument(c,tmp_path/'project.yml');w.refresh()
    w.workflow.accept(0,c);assert w.validate_stage(1);assert w.validate_c_part('simulation'), w.document.issues
    panel=w.elastic_panel;w.show();app.processEvents()
    buttons={b.text():b for b in panel.findChildren(QPushButton)}
    assert not buttons['Share One Refinement Parameter'].isHidden()
    assert not buttons['Make Independent'].isHidden()
    assert 'Link Selected' not in buttons and 'Unlink' not in buttons
    assert panel.table.horizontalHeaderItem(10).text()=='Parameter sharing'
    assert panel.dimension_counts()==(ordinary+2,before,before)
    indices=[i for i,r in enumerate(panel.rows) if r['hkl'] in ((2,0,0),(0,2,0)) and r['interval']==(0.,.5)]
    for i in indices:panel.table.selectionModel().select(panel.table.model().index(i,0),QItemSelectionModel.SelectionFlag.Select|QItemSelectionModel.SelectionFlag.Rows)
    prompts=[]
    monkeypatch.setattr(QInputDialog,'getText',lambda *a,**k:('Group 1',True))
    original=w.document.config;pending=w.pending;dirty=panel.dirty
    monkeypatch.setattr(QMessageBox,'question',lambda *a:QMessageBox.StandardButton.Cancel)
    buttons['Share One Refinement Parameter'].click()
    assert w.document.config==original and w.pending==pending and panel.dirty==dirty
    assert panel.dimension_counts()==(ordinary+2,before,before)
    def confirm(*args):prompts.append(args[2]);return QMessageBox.StandardButton.Yes
    monkeypatch.setattr(QMessageBox,'question',confirm)
    buttons['Share One Refinement Parameter'].click()
    assert f'{before} → {after}' in prompts[-1]
    assert ('one real and one imaginary' in prompts[-1]) if mode=='complex' else ('one '+('real' if mode=='real' else 'magnitude')+' optimizer parameter' in prompts[-1])
    panel=w.elastic_panel
    assert panel.dimension_counts()==(ordinary+2,before,after)
    assert f'Optimizer parameters: {after}' in panel.dimension_summary.text()
    assert all(panel.table.item(i,10).text()=='Shared: Group 1' for i in indices)
    compiled=compile_config_from_sources(w.document.config);assert compiled.parameter_count==after
    # Persisted labels remain unchanged because they are part of historical checkpoint identity.
    from qcbed.crystallography_refinement import physical_target_layout
    assert physical_target_layout(compiled)[-1]['label'].startswith('Linked coefficient group ')
    path=tmp_path/'shared.yml';save_yaml_config(w.document.config,path,base_directory=tmp_path)
    assert load_yaml_config(path)==w.document.config
    assert 'equality_group: Group 1' in path.read_text()
    assert w.validate_c_part('refinement'), w.document.issues
    panel=w.elastic_panel
    assert f'Final optimizer-vector dimension: {after}' in panel.resolved.toPlainText()
    assert f'Optimizer parameters: {after}' in w.summary.toPlainText()
    assert w.run_dashboard.status['Optimizer parameters'].text()==str(after)+' committed'
    assert w.run_dashboard.status['Active physical targets'].text()==str(ordinary+2)+' committed'
    panel.table.selectRow(indices[0])
    next(b for b in panel.findChildren(QPushButton) if b.text()=='Make Independent').click()
    assert w.elastic_panel.dimension_counts()==(ordinary+2,before,before)
    assert compile_config_from_sources(w.document.config).parameter_count==before
    w.document.saved_config=w.document.config;w.pending=False;w.close()


def test_results_labels_present_shared_without_rewriting_checkpoint(tmp_path):
    from test_coefficient_equality import linked
    from qcbed import prepare_run_snapshot,run_snapshot
    from qcbed.gui.results import load_results
    from pathlib import Path
    import hashlib
    c=linked(tmp_path);s=prepare_run_snapshot(c,base_directory=tmp_path)
    run_snapshot(s,cancel_requested=lambda:True)
    checkpoint=Path(c.run.output_directory)/'checkpoint_latest.json';before=hashlib.sha256(checkpoint.read_bytes()).hexdigest()
    r=load_results(c.run.output_directory)
    assert any(p['label'].startswith('Shared refinement parameter group ') for p in r.parameters)
    assert not any('Linked coefficient group' in p['label'] for p in r.parameters)
    assert hashlib.sha256(checkpoint.read_bytes()).hexdigest()==before
