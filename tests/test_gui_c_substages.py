"""Real widget transitions; no optimizer/objective evaluation."""
from dataclasses import replace
import pytest
pytest.importorskip('PySide6')
from PySide6.QtWidgets import QApplication, QPushButton
from PySide6.QtCore import Qt
from qcbed.gui.window import ConfigWindow
from qcbed.gui.document import EditorDocument
from test_modern_grouped_run import modern_config


def test_c_simulation_refinement_transitions(tmp_path):
    app=QApplication.instance() or QApplication([])
    c=modern_config(tmp_path)
    w=ConfigWindow();w.document=EditorDocument(c,tmp_path/'project.yml')
    w.refresh();w.show();app.processEvents()
    try:
        assert [w.c_tabs.tabText(i) for i in range(w.c_tabs.count())]==['Simulation Setup','Refinement Targets','Optimizer']
        assert w.navigation.topLevelItem(2).childCount()==3
        assert w.refinement_table.rowCount()==0
        assert not w.validate_refinement_button.isEnabled()
        assert 'Validate Simulation Setup' in w.refinement_notice.text()
        w.workflow.accept(0,c) # A scientific validation has independent coverage.
        assert w.validate_stage(1)
        w.workflow_tabs.setCurrentIndex(2)
        w.widgets[('crystallography','depth_angstrom')].setText('10')
        w.findChild(QPushButton,'validateSimulationSetup').click();app.processEvents()
        assert w.workflow.simulation_receipt is not None,w.stage_status[2].text()
        assert w.c_tabs.isTabEnabled(1)
        table=w.refinement_table
        row=next(r for r in range(table.rowCount()) if table.item(r,0).data(Qt.ItemDataRole.UserRole).category.value=='depth')
        assert table.item(row,3).text()=='10.0'
        old=w.document.config.refinement.parameters
        table.item(row,4).setText('')
        assert 'Incomplete — advanced search settings required' in table.item(row,2).text()
        assert not w.validate_c_part('refinement')
        assert w.document.config.refinement.parameters==old
        table.item(row,4).setText('.9')
        assert w.validate_c_part('refinement'),w.stage_status[2].text()
        assert w.workflow.refinement_receipt is not None
        binding=w.document.config.refinement.parameters
        w.widgets[('crystallography','depth_angstrom')].setText('12')
        assert w.workflow.simulation_receipt is None
        assert w.workflow.refinement_receipt is None
        assert not w.validate_stage(2)
        assert w.validate_c_part('simulation')
        assert w.document.config.refinement.parameters==binding
        assert w.workflow.refinement_receipt is None
        assert w.validate_c_part('refinement')
        assert w.validate_stage(2),w.stage_status[2].text()
        assert not w.document.files_checked
        assert w.workflow.ready
        # Execution edits do not discard physical/refinement receipts.
        w.widgets[('refinement','optimizer','seed')].setText('7')
        assert w.workflow.simulation_receipt is not None
        assert w.workflow.refinement_receipt is not None
        assert not w.workflow.ready
        assert w.validate_stage(2)
        w.workflow.invalidate(1)
        assert w.workflow.simulation_receipt is None
        assert w.workflow.refinement_receipt is None
    finally:
        w.document.saved_config=w.document.config;w.pending=False;w.close()


def test_refinement_bounds_and_execution_are_separate(tmp_path):
    app=QApplication.instance() or QApplication([])
    c=modern_config(tmp_path)
    w=ConfigWindow();w.document=EditorDocument(c,tmp_path/'project.yml');w.refresh()
    try:
        w.workflow.accept(0,c);assert w.validate_stage(1)
        assert w.validate_c_part('simulation')
        table=w.refinement_table
        row=next(r for r in range(table.rowCount()) if table.item(r,0).data(Qt.ItemDataRole.UserRole).category.value=='depth')
        table.item(row,8).setText('')
        assert 'Incomplete — bounds required' in table.item(row,2).text()
        assert 'D1, D2' not in table.item(row,2).text()
        old=w.document.config
        assert not w.validate_c_part('refinement')
        assert w.document.config==old
        table.item(row,8).setText('1')
        table.item(row,6).setText('-.01')
        assert w.validate_c_part('refinement') # Signed simplex metadata is preserved.
        assert any(b.plamda==-.01 for b in w.document.config.refinement.parameters)
        from pathlib import Path
        Path(c.run.output_directory).mkdir()
        w.refresh()
        from qcbed.gui.execution_details import external_checks
        assert 'Output directory: PASS' in external_checks(w.document)
        assert not any('confirmation required' in line for line in external_checks(w.document))
        assert w.workflow.simulation_receipt is not None
        assert w.workflow.refinement_receipt is not None
    finally:
        w.document.saved_config=w.document.config;w.pending=False;w.close()
