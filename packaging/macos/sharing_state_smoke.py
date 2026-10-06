"""Exercise public widget actions; record every committed state without optimization."""
import json
from PySide6.QtCore import Qt,QTimer
from PySide6.QtWidgets import QPushButton,QInputDialog,QMessageBox
from PySide6.QtTest import QTest
from qcbed import compile_config_from_sources
from qcbed.crystallography_refinement import physical_target_layout,apply_crystallographic_parameters


def exercise(w,app,out,mode='real'):
    states=[]
    def select():
        w.workflow_tabs.setCurrentIndex(2);w.c_tabs.setCurrentIndex(1);app.processEvents()
        panel=w.elastic_panel
        indices=[i for i,r in enumerate(panel.rows) if r['hkl'] in ((2,0,0),(0,2,0)) and r['interval']==(0.,1.)]
        panel.table.clearSelection()
        for n,i in enumerate(indices):
            item=panel.table.item(i,1);panel.table.scrollToItem(item);app.processEvents()
            QTest.mouseClick(panel.table.viewport(),Qt.MouseButton.LeftButton,Qt.KeyboardModifier.ControlModifier if n else Qt.KeyboardModifier.NoModifier,panel.table.visualItemRect(item).center())
        return panel,indices
    def snapshot(name):
        panel=w.elastic_panel;c=compile_config_from_sources(w.document.config);layout=physical_target_layout(c)
        vector=c.packed_starting_vector.copy()
        for i,p in enumerate(c.layout.parameters):
            if p.category=='projected_elastic':vector[i]+=.01
        try:
            trial=apply_crystallographic_parameters(c,vector);writeback_error=None
        except ValueError as exc:
            trial=c;writeback_error=str(exc)
        states.append(dict(name=name,mode=mode,selected=[panel.rows[i]['hkl'] for i in panel.selected_rows()],
            visible=[dict(hkl=r['hkl'],refine=panel.table.item(i,0).checkState()==Qt.CheckState.Checked,sharing=panel.table.item(i,10).text()) for i,r in enumerate(panel.rows) if r['hkl'] in ((2,0,0),(0,2,0))],
            committed=[dict(hkl=r.hkl,group=r.equality_group,mode=r.mode,refine=r.refine) for r in w.document.config.crystallography.projected_elastic],
            receipt=repr(w.workflow.refinement_receipt),simulation_receipt=repr(w.workflow.simulation_receipt),message=panel.message.text(),validation_status=w.stage_status[2].text(),layout=layout,vector_length=len(vector),bounds_length=len(c.bounds.entries),
            displayed=w.run_dashboard.status['Optimizer parameters'].text(),physical=w.run_dashboard.status['Active physical targets'].text(),run_enabled=w.run_action.isEnabled(),
            writeback_error=writeback_error,writeback=[dict(hkl=r['hkl'],group=r['record'].equality_group,value=[r['initial'].real,r['initial'].imag]) for r in trial.crystallography.projected_elastic]))
        (out/f'sharing-{mode}.json').write_text(json.dumps(states,indent=2,default=str))
    w.workflow_tabs.setCurrentIndex(2);w.c_tabs.setCurrentIndex(1);app.processEvents()
    t=w.refinement_table
    for i in range(t.rowCount()):
        depth=t.item(i,0).data(Qt.ItemDataRole.UserRole).category=='depth'
        t.item(i,1).setCheckState(Qt.CheckState.Checked if depth else Qt.CheckState.Unchecked)
        if depth:
            value=w.document.config.crystallography.depth_angstrom
            t.item(i,8).setText(str(value*.5));t.item(i,9).setText(str(value*1.5))
    panel,indices=select()
    for i in indices:
        panel.table.item(i,0).setCheckState(Qt.CheckState.Checked)
        panel.table.cellWidget(i,4).setCurrentIndex(panel.table.cellWidget(i,4).findData(mode))
        panel.table.item(i,5).setText('19.5 +0i' if mode=='complex' else '19.5')
        panel.table.item(i,6).setText('0, -5' if mode=='complex' else '0')
        panel.table.item(i,7).setText('40, 5' if mode=='complex' else '40')
    w.validate_refinement_button.click();app.processEvents();snapshot('independent_validated')
    panel,indices=select();next(b for b in panel.findChildren(QPushButton) if b.text()=='Make Independent').click()
    w.validate_refinement_button.click();app.processEvents();snapshot('made_independent_validated')
    panel,indices=select();snapshot('selected_before_share')
    timer=QTimer()
    def dialogs():
        d=app.activeModalWidget()
        if isinstance(d,QInputDialog):d.setTextValue('Audit pair');d.accept()
        elif isinstance(d,QMessageBox):d.button(QMessageBox.StandardButton.Yes).click()
    timer.timeout.connect(dialogs);timer.start(10)
    next(b for b in panel.findChildren(QPushButton) if b.text()=='Share One Refinement Parameter').click()
    timer.stop();app.processEvents();snapshot('shared_before_validation')
    w.validate_refinement_button.click();app.processEvents();snapshot('shared_validated')
    panel,indices=select();next(b for b in panel.findChildren(QPushButton) if b.text()=='Make Independent').click()
    app.processEvents();snapshot('independent_before_validation')
    w.validate_refinement_button.click();app.processEvents();snapshot('independent_revalidated')
    return states
