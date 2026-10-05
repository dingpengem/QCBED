from dataclasses import replace
import pytest
pytest.importorskip('PySide6')
from PySide6.QtWidgets import QApplication
from PySide6.QtCore import Qt
from qcbed.gui.window import ConfigWindow
from qcbed.gui.document import EditorDocument
from qcbed import compile_config_from_sources
from test_modern_grouped_run import modern_config


def test_refinement_drafts_commit_and_readiness(tmp_path):
    app=QApplication.instance() or QApplication([])
    c=modern_config(tmp_path)
    c=replace(c,crystallography=replace(c.crystallography,depth_angstrom=100.),
        refinement=replace(c.refinement,parameters=(),bounds=replace(c.refinement.bounds,specific=())))
    w=ConfigWindow();w.document=EditorDocument(c,tmp_path/'project.yml')
    w.refresh();w.show();w.workflow_tabs.setCurrentIndex(2)
    w.navigation.setCurrentItem(w.section_items['Refinement']);app.processEvents()
    w.workflow.accept(0,c);assert w.validate_stage(1)
    w.widgets[('crystallography','depth_angstrom')].setText('1000')
    assert w.validate_c_part('simulation')
    table=w.refinement_table
    rows={table.item(r,0).data(Qt.ItemDataRole.UserRole).category.value:r for r in range(table.rowCount())}
    row=rows['depth']
    table.item(row,1).setCheckState(Qt.CheckState.Checked)
    assert 'Incomplete' in table.item(row,2).text()
    assert table.item(row,3).text()=='1000.0'
    assert not w.apply_edits()
    assert all(table.item(row,col).text() for col in (4,5,6))
    assert w.document.config.crystallography.depth_angstrom==1000.
    for category in ('depth','debye_waller'):
        r=rows[category]
        table.item(r,1).setCheckState(Qt.CheckState.Checked)
        for col,text in ((4,'.8'),(5,'1.2'),(6,'.1'),(8,'600' if category=='depth' else '.5'),(9,'1200' if category=='depth' else '1.5')):
            table.item(r,col).setText(text)
    assert w.validate_c_part('refinement')
    compiled=compile_config_from_sources(w.document.config)
    from qcbed.refinement import pack_refinement_parameters
    assert tuple(pack_refinement_parameters(compiled.state))==(.86,1000.) or tuple(pack_refinement_parameters(compiled.state))==(1000.,.86)
    table=w.refinement_table
    values={table.item(r,0).data(Qt.ItemDataRole.UserRole).category.value:table.item(r,3).text() for r in range(table.rowCount())}
    assert values['depth']=='1000.0' and values['debye_waller']=='0.86'
    w.workflow.accept(0,w.document.config)
    assert w.validate_stage(1)
    assert w.validate_c_part('simulation')
    assert w.validate_c_part('refinement')
    assert w.validate_stage(2),w.stage_status[2].text()
    assert w.workflow.ready and not w.run_action.isEnabled()
    assert w.validate_output_directory()
    w.check_run_readiness()
    assert w.run_action.isEnabled()
    assert not w.resume_action.isEnabled()
    w.document.saved_config=w.document.config;w.pending=False;w.close()


def test_zero_targets_are_not_a_file_failure_and_report_exact_states(tmp_path):
    from qcbed.modern_frontend import validate_modern
    from qcbed.config_validation import ValidationPhase
    from qcbed.gui.execution_details import describe_execution
    app=QApplication.instance() or QApplication([])
    c=modern_config(tmp_path)
    c=replace(c,refinement=replace(c.refinement,parameters=()))
    issues=validate_modern(c,tmp_path,check_files=True)
    assert any(i.code=='grouped_capability' for i in issues)
    assert not any(i.code=='experimental_data_error' for i in issues)
    w=ConfigWindow();w.document=EditorDocument(c,tmp_path/'project.yml');w.refresh()
    w.workflow.accept(0,c);assert w.validate_stage(1)
    assert w.validate_c_part('simulation')
    report=describe_execution(w.document,w.refinement_table)
    assert 'committed binding: False' in report
    assert 'not selected' in report
    assert 'Prepared experimental data '+c.experiment.path+': PASS' in report
    assert 'Output directory: PASS' in report
    broken=replace(w.document.config,experiment=replace(w.document.config.experiment,path=str(tmp_path/'missing.dat')))
    w.document.config=broken
    from qcbed.gui.execution_details import external_checks
    assert any('Prepared experimental data' in r and 'BLOCKED' in r for r in external_checks(w.document))
    w.document.saved_config=w.document.config;w.close()


def test_dynamic_absorption_and_scan_bindings_through_widgets(tmp_path):
    from qcbed.grouped_run import check_grouped
    app=QApplication.instance() or QApplication([])
    c=modern_config(tmp_path)
    c=replace(c,refinement=replace(c.refinement,parameters=(),active_categories=None,
        bounds=replace(c.refinement.bounds,specific=())))
    w=ConfigWindow();w.document=EditorDocument(c,tmp_path/'project.yml');w.refresh();w.show()
    w.workflow_tabs.setCurrentIndex(2);w.navigation.setCurrentItem(w.section_items['Refinement']);app.processEvents()
    w.workflow.accept(0,c);assert w.validate_stage(1)
    assert w.validate_c_part('simulation')
    table=w.refinement_table
    for row in range(table.rowCount()):
        t=table.item(row,0).data(Qt.ItemDataRole.UserRole)
        if (t.category.value=='absorption' and t.identifiers[-1]=='b') or (t.category.value=='scan_coordinate' and t.identifiers==('A','first')):
            assert table.item(row,1).flags() & Qt.ItemFlag.ItemIsEnabled
            table.item(row,1).setCheckState(Qt.CheckState.Checked)
            for col,value in ((4,'.9'),(5,'1.1'),(6,'.01'),(8,'.5' if t.category.value=='absorption' else '-2'),(9,'1.5' if t.category.value=='absorption' else '2')):
                table.item(row,col).setText(value)
    assert w.apply_edits()
    c=compile_config_from_sources(w.document.config)
    assert [p['category'] for p in check_grouped(c)]==['absorption','scan_coordinate']
    assert all(w.refinement_table.item(r,2).text()=='Active' for r in range(w.refinement_table.rowCount()) if w.refinement_table.item(r,1).checkState()==Qt.CheckState.Checked)
    w.document.saved_config=w.document.config;w.pending=False;w.close()
