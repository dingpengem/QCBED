from dataclasses import replace
import pytest
from qcbed.gui.document import EditorDocument,editable_fields,parse_text
from test_modern_grouped_run import modern_config


@pytest.mark.parametrize('text',['0','0,.5,.5,1','0,.5,.9','-.1,.5,1','0,nan,1',''])
def test_invalid_boundaries_rejected_before_compile(text):
    with pytest.raises(ValueError):parse_text(text,'boundaries')


def test_boundaries_editor_and_incomplete_draft(tmp_path):
    c=modern_config(tmp_path)
    key=('crystallography','slice_boundaries_fractional')
    spec=next(s for s in editable_fields(c) if s.key==key)
    assert spec.section=='Numerical / Simulation cell'
    doc=EditorDocument(replace(c,crystallography=replace(c.crystallography,
        slice_boundaries_fractional=(0.,1.))),tmp_path/'draft.yml')
    assert doc.summary()['compile'].startswith('Valid')
    assert doc.apply({key:'0,1'})
    assert doc.apply({key:'[0, 0.5, 1]'})
    assert doc.config.crystallography.slice_boundaries_fractional==(0.,.5,1.)


def test_boundary_field_in_c(tmp_path):
    pytest.importorskip('PySide6')
    from PySide6.QtWidgets import QApplication
    from qcbed.gui.window import ConfigWindow
    app=QApplication.instance() or QApplication([])
    w=ConfigWindow()
    w.document=EditorDocument(modern_config(tmp_path),tmp_path/'config.yml')
    w.refresh()
    section=w.section_items['Numerical / Simulation cell']
    assert section.parent().text(0)=='C. Main QCBED Parameters'
    field=w.widgets[('crystallography','slice_boundaries_fractional')]
    assert not field.isReadOnly()
    field.setText('0, 0.25, 0.5, 1')
    assert w.apply_edits()
    assert w.document.config.crystallography.slice_boundaries_fractional==(0.,.25,.5,1.)
    w.document.saved_config=w.document.config
    w.close()


def test_validate_c_invalid_draft_and_reference_roundtrip(tmp_path):
    pytest.importorskip('PySide6')
    from PySide6.QtWidgets import QApplication
    from qcbed.gui.window import ConfigWindow
    from qcbed import load_yaml_config
    from pathlib import Path
    app=QApplication.instance() or QApplication([])
    w=ConfigWindow()
    c=modern_config(tmp_path)
    w.document=EditorDocument(c,tmp_path/'config.yml')
    w.refresh()
    key=('crystallography','slice_boundaries_fractional')
    w.workflow.accept(0,c);assert w.validate_stage(1)
    w.workflow_tabs.setCurrentIndex(2)
    w.widgets[key].setText('')
    assert not w.validate_c_part('simulation')
    assert w.stage_status[2].text().startswith('Slice boundaries must start at 0, end at 1, and contain at least 2 strictly increasing values.')
    assert not w.run_action.isEnabled()
    # Replacing the CIF clears structure metadata, but retains C slicing.
    w.document.load_structure_file(Path(__file__).resolve().parents[1]/'reference_data/Al_FCC.cif')
    assert w.document.config.crystallography.slice_boundaries_fractional==c.crystallography.slice_boundaries_fractional
    w.document.config=replace(c,crystallography=replace(c.crystallography,slice_boundaries_fractional=(0.,.5,1.)))
    w.refresh()
    assert w.save_document(tmp_path/'saved.yml')
    assert load_yaml_config(tmp_path/'saved.yml').crystallography.slice_boundaries_fractional==(0.,.5,1.)
    w.document.saved_config=w.document.config;w.close()


def test_actual_visible_c_editor_after_tab_reparent(tmp_path):
    pytest.importorskip('PySide6')
    from PySide6.QtWidgets import QApplication
    from qcbed.gui.window import ConfigWindow
    app=QApplication.instance() or QApplication([])
    w=ConfigWindow();w.document=EditorDocument(modern_config(tmp_path),tmp_path/'config.yml')
    w.refresh();w.show();app.processEvents()
    w.workflow_tabs.setCurrentIndex(2)
    w.navigation.setCurrentItem(w.section_items['Numerical / Simulation cell'])
    app.processEvents()
    assert not w.numerical_advanced.toggle.isChecked();app.processEvents()
    widget=w.widgets[('crystallography','slice_boundaries_fractional')]
    assert widget.isVisible()
    assert widget.visibleRegion().boundingRect().width()>100
    widget.setText('0, 0.5, 1')
    w.workflow.accept(0,w.document.config);assert w.validate_stage(1)
    assert w.validate_c_part('simulation')
    assert w.document.config.crystallography.slice_boundaries_fractional==(0.,.5,1.)
    w.document.saved_config=w.document.config;w.pending=False;w.close()
