"""Actual widgets: geometry, regeneration, sharing, and information disclosure."""
import importlib.util
from pathlib import Path
from PySide6.QtWidgets import QPushButton,QApplication
from PySide6.QtCore import Qt
from test_workflow_navigation import window


def harness():
    path=Path(__file__).parents[1]/'packaging/macos/c_layout_smoke.py'
    spec=importlib.util.spec_from_file_location('layout_smoke',path)
    m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m


def test_c_responsive_regeneration_sharing(window,tmp_path):
    w=window;app=QApplication.instance()
    assert w.validate_stage(1);assert w.validate_c_part('simulation')
    w.workflow_tabs.setCurrentIndex(2);w.c_tabs.setCurrentIndex(1)
    for _ in range(2):
        w.refresh();assert w.c_tabs.count()==3
    p=w.elastic_panel
    assert 'Check Refine' in p.sharing_instructions.text() and 'Select two or more' in p.sharing_instructions.text()
    assert 'Make Independent' in p.sharing_semantics.text()
    assert not p.sharing_example.toggle.isChecked()
    indices=[i for i,r in enumerate(p.rows) if r['record']]
    p.link_rows(indices,'Pair')
    assert 'Shared Pair:' in w.elastic_panel.group_summary.text()
    w.validate_refinement_button.click();assert w.workflow.refinement_receipt is not None
    assert w.c_tabs.currentIndex()==1
    assert not w.elastic_panel.layout_details.toggle.isChecked()
    harness().audit(w,app,tmp_path)


def test_action_row_stacks_without_collision(window):
    from qcbed.gui.responsive_actions import ResponsiveActions
    from PySide6.QtWidgets import QApplication
    app=QApplication.instance()
    row=ResponsiveActions([QPushButton('Share One Refinement Parameter'),QPushButton('Make Independent')]);row.show()
    try:
        row.resize(260,80);harness().settle(app)
        assert row.items[0].geometry().bottom()<row.items[1].geometry().top()
        for b in row.items:assert row.rect().contains(b.geometry())
        row.resize(800,80);harness().settle(app)
        assert row.items[0].geometry().right()<row.items[1].geometry().left()
    finally:row.close()
