"""Physical selection and advanced search controls through real Qt widgets."""
from dataclasses import replace
import pytest
pytest.importorskip('PySide6')
from PySide6.QtWidgets import QApplication
from PySide6.QtCore import Qt
from qcbed.gui.window import ConfigWindow
from qcbed.gui.document import EditorDocument
from qcbed.optimizer_defaults import recommended_binding
from qcbed import compile_config_from_sources,load_yaml_config
from test_modern_grouped_run import modern_config


def window(tmp_path,empty=True):
    c=modern_config(tmp_path)
    if empty:c=replace(c,refinement=replace(c.refinement,parameters=(),bounds=replace(c.refinement.bounds,specific=())))
    w=ConfigWindow();w.document=EditorDocument(c,tmp_path/'project.yml');w.refresh();w.show()
    w.workflow.accept(0,c);assert w.validate_stage(1)
    w.workflow_tabs.setCurrentIndex(2);assert w.validate_c_part('simulation')
    return w


def close(w):
    w.document.saved_config=w.document.config;w.pending=False;w.close()


@pytest.mark.parametrize('categories',[('depth',),('debye_waller',),('depth','debye_waller')])
def test_normal_targets_need_only_bounds(tmp_path,categories):
    app=QApplication.instance() or QApplication([]);w=window(tmp_path)
    try:
        assert [w.c_tabs.tabText(i) for i in range(w.c_tabs.count())]==['Simulation Setup','Refinement Targets','Optimizer']
        assert w.recommended_optimizer.isChecked()
        assert not w.advanced_optimizer.isVisible()
        t=w.refinement_table
        assert all(t.isColumnHidden(c) for c in (4,5,6))
        for r in range(t.rowCount()):
            target=t.item(r,0).data(Qt.ItemDataRole.UserRole)
            if target.category.value not in categories:continue
            t.item(r,1).setCheckState(Qt.CheckState.Checked)
            t.item(r,8).setText('1' if target.category.value=='depth' else '.5')
            t.item(r,9).setText('20' if target.category.value=='depth' else '1.5')
        assert w.validate_c_part('refinement'),w.stage_status[2].text()
        c=w.document.config
        assert len(c.refinement.parameters)==len(categories)
        assert all(b==recommended_binding(b.target) for b in c.refinement.parameters)
        from qcbed.configuration import ParameterBinding
        # Explicit pre-redesign reference metadata, independent of defaults resolver.
        expected=replace(c,refinement=replace(c.refinement,parameters=tuple(
            ParameterBinding(b.target,(1.05,.95),.5 if b.target.category.value=='depth' else .0002)
            for b in c.refinement.parameters)))
        a,b=compile_config_from_sources(c),compile_config_from_sources(expected)
        import numpy as np
        np.testing.assert_array_equal(a.packed_starting_vector,b.packed_starting_vector)
        assert a.range_factors==b.range_factors and a.plamda==b.plamda
        assert all(w.refinement_table.item(r,2).text()=='Active' for r in range(w.refinement_table.rowCount()) if w.refinement_table.item(r,1).checkState()==Qt.CheckState.Checked)
        assert w.save_document(tmp_path/'saved.yml')
        assert load_yaml_config(tmp_path/'saved.yml')==w.document.config
        import yaml
        targets=yaml.safe_load((tmp_path/'saved.yml').read_text())['qcbed']['refinement']['targets']
        assert all(value=={} for value in targets.values())
    finally:close(w)


def test_advanced_override_and_optimizer_controls(tmp_path):
    app=QApplication.instance() or QApplication([]);w=window(tmp_path,False)
    try:
        w.c_tabs.setCurrentIndex(2)
        w.advanced_optimizer_toggle.click();w.recommended_optimizer.setChecked(False)
        app.processEvents();assert w.advanced_optimizer.isVisible()
        t=w.optimizer_metadata
        t.item(0,3).setText('-.125')
        assert w.validate_c_part('refinement')
        assert any(b.plamda==-.125 for b in w.document.config.refinement.parameters)
        assert w.save_document(tmp_path/'override.yml')
        loaded=load_yaml_config(tmp_path/'override.yml')
        assert loaded==w.document.config
        w.document=EditorDocument(loaded,tmp_path/'override.yml');w.refresh()
        assert not w.recommended_optimizer.isChecked()
        before=w.document.config.refinement.parameters
        combo=w.widgets[('refinement','optimizer','method')]
        combo.setCurrentText('ANNEAL')
        assert w.widgets[('refinement','optimizer','maxtemp')].isEnabled()
        assert not w.widgets[('refinement','optimizer','ftol')].isEnabled()
        assert w.optimizer_metadata.item(0,3).text()=='Not used by the current optimizer'
        combo.setCurrentText('AMOEBA')
        assert not w.widgets[('refinement','optimizer','maxtemp')].isEnabled()
        assert w.widgets[('refinement','optimizer','ftol')].isEnabled()
        assert w.document.config.refinement.parameters==before
        combo.setCurrentText('PARAMIN')
        assert w.widgets[('refinement','optimizer','maxfpara')].isEnabled()
        assert w.optimizer_metadata.item(0,3).text()=='Not used by the current optimizer'
    finally:close(w)
