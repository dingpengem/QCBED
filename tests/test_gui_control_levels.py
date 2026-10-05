from dataclasses import replace
import pytest
pytest.importorskip('PySide6')
from PySide6.QtWidgets import QApplication
from qcbed.gui.window import ConfigWindow
from qcbed.gui.document import EditorDocument
from test_modern_grouped_run import modern_config


def test_classified_sections_preserve_physical_inputs(tmp_path):
    app=QApplication.instance() or QApplication([])
    c=modern_config(tmp_path)
    w=ConfigWindow();w.document=EditorDocument(c,tmp_path/'project.yml');w.refresh();w.show();app.processEvents()
    try:
        assert not w.preparation.advanced_settings.toggle.isChecked()
        assert w.preparation.binning.value()==2
        assert w.data_panel.mask_layer.count()==4
        assert not w.atomic_expert.toggle.isChecked()
        assert w.site_group_table.isColumnHidden(2)
        assert not w.site_group_table.isColumnHidden(3)
        assert not w.numerical_advanced.toggle.isChecked()
        assert w.recommended_numerics.isChecked()
        assert not w.optimizer_advanced_settings.toggle.isChecked()
        assert not w.execution_advanced.toggle.isChecked()
        assert w.document.config==c
        w.workflow.accept(0,c);assert w.validate_stage(1)
        w.workflow_tabs.setCurrentIndex(2)
        assert w.validate_c_part('simulation')
        assert w.refinement_table.rowCount()>0
        assert w.document.config==c
    finally:
        w.document.saved_config=w.document.config;w.pending=False;w.close()


def test_missing_numerics_are_exposed_not_invented(tmp_path):
    app=QApplication.instance() or QApplication([])
    c=modern_config(tmp_path);c=replace(c,crystallography=replace(c.crystallography,slice_boundaries_fractional=()))
    w=ConfigWindow();w.document=EditorDocument(c,tmp_path/'project.yml');w.refresh()
    try:
        assert not w.numerical_advanced.toggle.isChecked()
        assert w.recommended_numerics.isEnabled()
        assert w.widgets[('crystallography','slice_boundaries_fractional')].text()==''
        assert w.document.config==c
    finally:
        w.document.saved_config=w.document.config;w.pending=False;w.close()


def test_endpoint_only_boundary_is_visible_editable_and_valid(tmp_path):
    app=QApplication.instance() or QApplication([])
    c=modern_config(tmp_path)
    c=replace(c,crystallography=replace(c.crystallography,slice_boundaries_fractional=(0.,1.)))
    w=ConfigWindow();w.document=EditorDocument(c,tmp_path/'project.yml');w.refresh()
    try:
        key=('crystallography','slice_boundaries_fractional')
        assert not w.numerical_advanced.toggle.isChecked()
        assert w.recommended_numerics.isEnabled()
        assert w.widgets[key].text()=='0.0, 1.0'
        assert not w.widgets[key].isReadOnly()
        assert 'one slice' in w.slice_boundary_help.text()
        w.workflow.accept(0,c);assert w.validate_stage(1)
        w.workflow_tabs.setCurrentIndex(2)
        assert w.validate_c_part('simulation')
        assert w.workflow.simulation_receipt is not None
        assert not w.run_action.isEnabled()
        assert w.document.config==c
        w.widgets[key].setText('0, 0.5, 1')
        assert w.validate_c_part('simulation'),w.stage_status[2].text()
        assert w.document.config.crystallography.slice_boundaries_fractional==(0.,.5,1.)
    finally:
        w.document.saved_config=w.document.config;w.pending=False;w.close()
