"""Opt-in zero-evaluation navigation acceptance with multiple shared groups."""
import json
from PySide6.QtCore import Qt,QTimer,QItemSelectionModel
from PySide6.QtWidgets import QPushButton,QInputDialog,QMessageBox,QScrollArea
from sharing_state_smoke import exercise


def run(w,app,out):
    assert w.c_tabs.currentIndex()==0
    w.c_next['simulation'].click();assert w.c_tabs.currentIndex()==1
    exercise(w,app,out,'real')
    assert w.c_tabs.currentIndex()==1
    p=w.elastic_panel
    identities=((2,0,0),(0,2,0),(2,2,0),(2,-2,0))
    for i,r in enumerate(p.rows):
        if r['hkl'] not in identities:continue
        p.table.item(i,0).setCheckState(Qt.CheckState.Checked)
        p.table.item(i,5).setText('10');p.table.item(i,6).setText('0');p.table.item(i,7).setText('40')
    w.validate_refinement_button.click();assert w.workflow.refinement_receipt is not None
    for name,members in (('Axial',identities[:2]),('Diagonal',identities[2:])):
        p=w.elastic_panel;p.table.clearSelection()
        for i,r in enumerate(p.rows):
            if r['hkl'] in members:p.table.selectionModel().select(p.table.model().index(i,0),QItemSelectionModel.SelectionFlag.Select|QItemSelectionModel.SelectionFlag.Rows)
        timer=QTimer()
        def accept():
            dialog=app.activeModalWidget()
            if isinstance(dialog,QInputDialog):dialog.setTextValue(name);dialog.accept()
            elif isinstance(dialog,QMessageBox):dialog.button(QMessageBox.StandardButton.Yes).click()
        timer.timeout.connect(accept);timer.start(10)
        next(b for b in p.findChildren(QPushButton) if b.text()=='Share One Refinement Parameter').click();timer.stop()
        assert w.c_tabs.currentIndex()==1
        assert name in w.elastic_panel.live_groups.text()
        assert w.workflow.refinement_receipt is None and not w.run_action.isEnabled()
    w.validate_refinement_button.click();assert w.c_tabs.currentIndex()==1
    assert w.workflow.refinement_receipt is not None
    assert 'Shared parameter groups: 2' in w.elastic_panel.dimension_summary.text()
    assert 'Optimizer parameters: 3' in w.elastic_panel.dimension_summary.text()
    for scroll in w.findChildren(QScrollArea):
        if scroll.isAncestorOf(w.elastic_panel):scroll.ensureWidgetVisible(w.elastic_panel.group_summary)
    app.processEvents();w.grab().save(str(out/'shared-groups.png'))
    w.c_next['refinement'].click();assert w.c_tabs.currentIndex()==2
    next(b for b in w.c_tabs.widget(2).findChildren(QPushButton) if b.text()=='Validate Optimizer').click();assert w.workflow.receipts[2] is not None,w.stage_status[2].text()
    assert w.c_tabs.currentIndex()==2
    w.c_next['optimizer'].click();assert w.workflow_tabs.currentIndex()==3
    w.widgets[('run','output_directory')].setText(str(out/'run'))
    assert w.validate_output_directory();w.check_run_readiness()
    assert w.run_action.isEnabled(),w.run_dashboard.reason.text()
    assert 'Run is ready' in w.guidance.next.text()
    app.processEvents();w.grab().save(str(out/'workflow-ready.png'))
    (out/'result.json').write_text(json.dumps(dict(status='PASS',frozen=True,optimizer_evaluations=0,
        checks=['B Structure → explicit Site Groups; g1/g2 preserved','C validation remains in place','Actual Share/Make Independent widget route','Two independent user groups, five physical targets, three coordinates','Explicit Next to Optimizer and D','Persistent workflow guidance','Writable output and Run readiness']),indent=2))
    return 0
